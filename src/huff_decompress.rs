//! This module provides huffman encoding and decoding routines

use std::io::{ Read};

use crate::bitstream::BitStreamReader;
pub use crate::constants::LIMIT;
use crate::utils::REVERSED_BITS;

pub fn build_tree(table: &mut [u16; 1 << LIMIT], code_lengths: &[u8; LIMIT + 1], symbols: &[u8])
{
    /*
     * Build a Huffman tree from the code and the code lengths.
     *
     * The one thing to understand it that the table we create is LSB based
     * and not MSB based
     *
     * The poop here is that we do everything in MSB fashion except the
     * writing to the lookup table(because it's easier that way)
     *
     * MSB first decoding looks like
     * [code]_{extra_bits}
     * While lsb is
     * {extra_bits}_code.
     *
     * Filling a MSB first table is easy , and has the best locality since
     * it walks the table consecutively. Each code is gets 2^(LIMIT-code_len)
     * times but for LSB tables, entries are separated in 2^(LIMIT-n) times.
     * which kinda sucks since we start doing strided access , but the table is small
     * and LSB saves us some instructions in the main loop.
     *
     * Now to do LSB table lookup setting, I simply do an MSB table lookup but when writing
     * i do a MSB->LSB bit reversal on the table.
     * E.g say we have a code 11011_0000, corresponding to a symbol 234.
     *
     * What we do is simply reverse the bits to 0000_11011 (the encoding part reversed the bit-stream too, see utils.rs/to_u32()
     *
     * Effectively converting the MSB to an LSB quite cheaply
     *
     */
    let mut code = 0;

    let mut p = 0;

    for i in 1..=LIMIT
    {
        for _ in 0..code_lengths[i]
        {
            let look_bits = code << (LIMIT - i);

            for k in 0..(1 << (LIMIT - i))
            {
                // the top 5-16 bits contain the reversed stream, so we need
                // to remove the bottom ones.
                let reversed_bits =
                    (REVERSED_BITS[(look_bits + k) as usize] >> (16 - LIMIT)) as usize;
                // Do the bit-reversal here.
                let entry = &mut table[reversed_bits];

                // format
                // symbol-> 8-16
                // bits_consumed -> 0-8 (makes it easier to retrieve on x86, just an al)
                *entry = u16::from(symbols[p]) << 8 | (i as u16);
            }
            p += 1;

            code += 1;
        }
        // do not go one code length higher
        assert!(code <= 1 << i);
        code *= 2;
    }
}

fn decompress_huff_inner(
    buf: &[u8], entries: &[u16; 1 << LIMIT], offsets: &[usize; 5], block_size: usize,
    dest: &mut [u8],
)
{
    /*
     * The difference with compress_huff are small but subtle hence deserve their
     * own explanation
     *
     * In, we compress 5 bytes per go because the CPU has more work to do
     * here in decompression, we enter a fighting state where we try to decompress
     * as much as possible from different streams in flight(main loop).
     *
     * Why is because of data dependencies.
     *
     * Serially, decoding a symbol requires 5-7 cycles per byte
     * [https://fgiesen.wordpress.com/2021/08/30/entropy-coding-in-oodle-data-huffman-coding/]
     *
     * So we can't let 5 bytes to decode serially before switching to another stream
     * (e.g how encode works) since we create a large serial dependency (the CPU doesn't
     * have an infinite out of order window). So that's why we decode per byte
     * serially hoping we reach a steady state where we are limited by instruction
     * throughput and not latency
     *
     */

    // initialize streams with pointers to compressed data.
    let mut stream1 = BitStreamReader::new(&buf[0..offsets[0]]);

    let mut stream2 = BitStreamReader::new(&buf[offsets[0]..offsets[1]]);

    let mut stream3 = BitStreamReader::new(&buf[offsets[1]..offsets[2]]);

    let mut stream4 = BitStreamReader::new(&buf[offsets[2]..offsets[3]]);

    let mut stream5 = BitStreamReader::new(&buf[offsets[3]..offsets[4]]);

    let stream_size = (block_size + 4) / 5;

    // Split the large vector into small ones
    let (dest1, rem) = dest.split_at_mut(stream_size);

    let (dest2, rem) = rem.split_at_mut(stream_size);

    let (dest3, rem) = rem.split_at_mut(stream_size);

    let (dest4, dest5) = rem.split_at_mut(stream_size);

    let mut min_pos = 0;

    unsafe {
        for ((((a, b), c), d), e) in dest1
            .chunks_exact_mut(10)
            .zip(dest2.chunks_exact_mut(10))
            .zip(dest3.chunks_exact_mut(10))
            .zip(dest4.chunks_exact_mut(10))
            .zip(dest5.chunks_exact_mut(10))
        {
            /*
             * There is no need to check if we're nearing the end of the loop.
             * The above zip does it for us.
             */
            // decode bytes, 50 symbols per loop.

            // Unrolling by 2x improves speed from 1635 Mb/s to 1700 Mb/s(mac-os M1 arm 2020)
            macro_rules! decode_single {
                ($index:tt) => {
                    stream1.decode_single(a.get_mut($index).unwrap(), entries);

                    stream2.decode_single(b.get_mut($index).unwrap(), entries);

                    stream3.decode_single(c.get_mut($index).unwrap(), entries);

                    stream4.decode_single(d.get_mut($index).unwrap(), entries);

                    stream5.decode_single(e.get_mut($index).unwrap(), entries);
                };
            }
            macro_rules! refill {
                () => {
                    stream1.refill_fast();

                    stream2.refill_fast();

                    stream3.refill_fast();

                    stream4.refill_fast();

                    stream5.refill_fast();
                };
            }
            refill!();

            decode_single!(0);

            decode_single!(1);

            decode_single!(2);

            decode_single!(3);

            decode_single!(4);

            refill!();

            decode_single!(5);

            decode_single!(6);

            decode_single!(7);

            decode_single!(8);

            decode_single!(9);

            min_pos += 10;
        }
    }
    /*
     *  Do the last bytes separately
     * note, we don't check if a stream overwrote to the next start of
     * the stream
     * its practically impossible (thank you Rust) because of zips and
     * the fact that slices track their sizes and panic if you write beyond
     */
    let mut last_pos = min_pos;

    //stream1-stream4 decode equal bits, (since they are equal,
    // we can use one loop for that
    for (((a, b), c), d) in dest1[min_pos..]
        .chunks_mut(5)
        .zip(dest2[min_pos..].chunks_mut(5))
        .zip(dest3[min_pos..].chunks_mut(5))
        .zip(dest4[min_pos..].chunks_mut(5))
    {
        // we know lengths are equal, so we can do this
        let len = a.len();

        min_pos += len;
        // refill.
        unsafe {
            stream1.refill_fast();

            stream2.refill_fast();

            stream3.refill_fast();

            stream4.refill_fast();
        }
        for i in 0..len
        {
            stream1.decode_single(&mut a[i], entries);

            stream2.decode_single(&mut b[i], entries);

            stream3.decode_single(&mut c[i], entries);

            stream4.decode_single(&mut d[i], entries);
        }
    }
    // the remainder is easy
    // fill stream 5 up to what's left
    for a in dest5[last_pos..].chunks_mut(5)
    {
        unsafe {
            stream5.refill_fast();
        }
        last_pos += a.len();
        for i in a.iter_mut()
        {
            //decode each one slowly by slowly
            stream5.decode_single(i, entries);
        }
    }
    // Do the assertions
    // Assertion 1
    // All bytes were decoded.
    // min_pos contains bytes for stream 1-4, last pos contains bytes
    // decoded for stream 5, they should be equal to the end.
    assert_eq!(
        (min_pos * 4) + last_pos,
        block_size,
        "Not all bytes were decoded, internal error"
    );

    // check that no stream read more data than needed
    assert!(
        stream1.check_final(),
        "Stream 1 was possibly corrupted, position:{},length:{}",
        stream1.get_position(),
        stream1.get_src_len()
    );

    assert!(
        stream2.check_final(),
        "Stream 2 was possibly corrupted, current position:{},expected position:{}",
        stream2.get_position(),
        stream2.get_src_len()
    );

    assert!(
        stream3.check_final(),
        "Stream 3 was possibly corrupted, current position:{},expected position:{}",
        stream3.get_position(),
        stream3.get_src_len()
    );

    assert!(
        stream4.check_final(),
        "Stream 4 was possibly corrupted, current position:{},expected position:{}",
        stream4.get_position(),
        stream4.get_src_len()
    );

    assert!(
        stream5.check_final(),
        "Stream 5 was possibly corrupted, current position:{},expected position:{}",
        stream5.get_position(),
        stream5.get_src_len()
    );
    // everything is good, we won Mr Stark.
}
/// Read Huffman Compressed data from src and write into dest
///
/// Caveats to understand for dest
/// we will write from `dest.len()`  going forward, after decompression
/// is done , dest.len() will contain uncompressed data
pub fn huff_decompress_4x<R: Read>(src: &mut R, dest: &mut Vec<u8>)
{
    let mut length = [0, 0, 0, 0];
    // read the length
    src.read_exact(&mut length[0..3]).unwrap();

    let mut block_length = u32::from_le_bytes(length);

    // not all bytes will be used
    let mut source = vec![0; block_length as usize];

    // we know that blocks have equal sizes except the last one
    // so if we read the first one we can determine how much space we will need

    let mut tbl = [0; 1 << LIMIT];

    let mut checksum = [0; 3];

    let mut jump_table = [0; 10];

    let mut code_lengths = [0; 12];

    let mut symbols = [0; 256];

    loop
    {
        block_length = u32::from_le_bytes(length);

        if dest.capacity() <= (block_length as usize + dest.len()) 
        {
            
            dest.reserve(block_length as usize);
        }
        // read checksum
        src.read_exact(&mut checksum).unwrap();
        // read jump table
        src.read_exact(&mut jump_table).unwrap();

        // two bytes per jump table, stored in little endian form.
        let tbl1 = u32::from(jump_table[0]) + (u32::from(jump_table[1]) << 8);

        let tbl2 = u32::from(jump_table[2]) + (u32::from(jump_table[3]) << 8) + tbl1;

        let tbl3 = u32::from(jump_table[4]) + (u32::from(jump_table[5]) << 8) + tbl2;

        let tbl4 = u32::from(jump_table[6]) + (u32::from(jump_table[7]) << 8) + tbl3;

        // end of this stream.
        let end = u32::from(jump_table[8]) + (u32::from(jump_table[9]) << 8) + tbl4;

        let offsets = [tbl1, tbl2, tbl3, tbl4, end].map(|x| x as usize);

        // read code lengths
        src.read_exact(&mut code_lengths[1..]).unwrap();

        let codes = code_lengths.iter().map(|x| *x as usize).sum::<usize>();
        // read symbols
        src.read_exact(&mut symbols[0..codes]).unwrap();

        // Build the Huffman tree
        build_tree(&mut tbl, &code_lengths, &symbols);

        let huff_source = &mut source[0..end as usize];

        src.read_exact(huff_source).unwrap();
        // new length
        let start = dest.len();

        let new_len = dest.len() + block_length as usize;
        // set length to be the capacity
        if new_len > dest.capacity()
        {
            let cap = dest.capacity();
            dest.reserve(new_len - cap);
        }
        // Don't continue if we don't have capacity to create a new write
        assert!(
            new_len <= dest.capacity(),
            "{},{}",
            new_len,
            dest.capacity()
        );

        // this is the laziest way to use uninitialized memory
        // 1. We know decompress_huff_inner will write to the region
        // between
        unsafe { dest.set_len(new_len) };

        decompress_huff_inner(
            huff_source,
            &tbl,
            &offsets,
            block_length as usize,
            // write to unintialized memory :)
            dest.get_mut(start..).unwrap(),
        );
        // read the length for the next iteration
        if src.read_exact(length.get_mut(0..3).unwrap()).is_err()
        {
            // if there is no more data, we can effectively say we are done
            // no more streams.
            //
            // Although we should be able to indicate a last block using some bit in the
            // spec
            break;
        }
    }
}

#[test]
fn huff_decompress()
{
    use std::fs::OpenOptions;
    use std::io::{BufReader, BufWriter};

    use crate::huff_compress_4x;
    {
        let fs = OpenOptions::new()
            .create(true)
            .write(true)
            .truncate(true)
            .open("/Users/calebe/CLionProjects/zcif/tests.zcif")
            .unwrap();
        let mut fs = BufWriter::with_capacity(1 << 24, fs);

        let fd = OpenOptions::new()
            .read(true)
            .open("/Users/calebe/git/FiniteStateEntropy/programs/enwiki.small")
            .unwrap();
        let mut fd = BufReader::new(fd);

        huff_compress_4x(&mut fd, &mut fs);
    }
    {
        let fs = OpenOptions::new()
            .create(false)
            .write(false)
            .read(true)
            .open("/Users/calebe/CLionProjects/zcif/tests.zcif")
            .unwrap();
        let mut fs = BufReader::with_capacity(1 << 24, fs);

        let fd = OpenOptions::new()
            .create(true)
            .write(true)
            .open("/Users/calebe/git/FiniteStateEntropy/programs/enwiki.xml")
            .unwrap();
        let mut fd = BufWriter::new(fd);

        //huff_decompress_4x(&mut fs, &mut fd);
    }
}
