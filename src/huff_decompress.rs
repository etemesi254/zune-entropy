//! This module provides huffman encoding and decoding routines

pub use crate::constants::LIMIT;
use crate::errors::EntropyErrors;
use crate::huff_bitstream::BitStreamReader;
use crate::unsafe_utils::extend;
use crate::utils::{read_rle, read_uncompressed, REVERSED_BITS};

pub fn build_tree(
    table: &mut [u16; 1 << LIMIT], code_lengths: &[u8; LIMIT + 1], symbols: &[u8],
) -> Result<(), EntropyErrors>
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
        if code > 1 << i
        {
            return Err(EntropyErrors::CorruptHeader(
                "Code length above expected code".to_string(),
            ));
        }
        assert!(code <= 1 << i);
        code *= 2;
    }
    Ok(())
}
fn decode_symbols(
    buf: &[u8], entries: &[u16; 1 << LIMIT], offsets: &[usize; 5], block_size: usize,
    dest: &mut [u8],
) -> Result<(), EntropyErrors>
{
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        unsafe {
            if is_x86_feature_detected!("bmi2")
            {
                return decode_symbols_bmi2(buf, entries, offsets, block_size, dest);
            }
        }
    }
    decode_symbols_fallback(buf, entries, offsets, block_size, dest)
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "bmi2")]
unsafe fn decode_symbols_bmi2(
    buf: &[u8], entries: &[u16; 1 << LIMIT], offsets: &[usize; 5], block_size: usize,
    dest: &mut [u8],
) -> Result<(), EntropyErrors>
{
    decode_symbols_fallback(buf, entries, offsets, block_size, dest)
}

#[inline(always)]
fn decode_symbols_fallback(
    buf: &[u8], entries: &[u16; 1 << LIMIT], offsets: &[usize; 5], block_size: usize,
    dest: &mut [u8],
) -> Result<(), EntropyErrors>
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
    if (min_pos * 4) + last_pos != block_size
    {
        return Err(EntropyErrors::CorruptStream(format!(
            "Not all bytes were decoded,expected:{},decoded:{}",
            block_size,
            (min_pos * 4) + last_pos
        )));
    }

    // check that no stream read more data than needed
    if !stream1.check_final()
        | !stream2.check_final()
        | !stream3.check_final()
        | !stream4.check_final()
        | !stream5.check_final()
    {
        return Err(EntropyErrors::CorruptStream(
            "A stream  was possibly corrupted, integrity of data is compromised".to_string(),
        ));
    }
    //everything is good, we won Mr Stark.
    Ok(())
}
/// Read Huffman Compressed data from src and write into dest
///
/// Caveats to understand for dest
/// we will write from `dest.len()`  going forward, after decompression
/// is done , dest.len() will contain uncompressed data
pub fn huff_decompress(src: &[u8], dest: &mut Vec<u8>) -> Result<(), EntropyErrors>
{
    let mut block_length;
    let mut symbols: &[u8];
    let mut jump_table: &[u8; 10];

    let mut tbl = [0; 1 << LIMIT];
    let mut code_lengths: [u8; 12] = [0; 12];

    // current position from where we are reading src_buf from.
    let mut src_position = 4;

    // read block information
    let mut block_info = src[0];

    // read the length
    let mut length = [0, 0, 0, 0];
    length[0..3].copy_from_slice(&src[1..4]);

    // we know that blocks have equal sizes except the last one
    // so if we read the first one we can determine how much space we will need
    loop
    {
        block_length = u32::from_le_bytes(length);

        if dest.capacity() <= (block_length as usize + dest.len())
        {
            dest.reserve(block_length as usize);
        }
        // 0b10 uncompressed
        if (block_info >> 6) == 0b10
        {
            read_uncompressed(&src[src_position..], block_length, dest);
            src_position += block_length as usize;
        }
        else if (block_info >> 6) == 0b01
        {
            // RLE block
            read_rle(&src[src_position..], block_length, dest);
            src_position += 1;
        }
        else if (block_info >> 6) == 0b11
        {
            // fse compressed block
            panic!("FSE block encountered where Huffman Was expected");
        }
        else
        {
            // compressed block.

            // read jump table
            jump_table = src[src_position..src_position + 10].try_into().unwrap();
            src_position += 10;
            //src.read_exact(&mut jump_table)?;

            // two bytes per jump table, stored in little endian form.
            let tbl1 = usize::from(u16::from_le_bytes(jump_table[0..2].try_into().unwrap()));

            let tbl2 = usize::from(u16::from_le_bytes(jump_table[2..4].try_into().unwrap())) + tbl1;

            let tbl3 = usize::from(u16::from_le_bytes(jump_table[4..6].try_into().unwrap())) + tbl2;

            let tbl4 = usize::from(u16::from_le_bytes(jump_table[6..8].try_into().unwrap())) + tbl3;

            // end of this stream.
            let end = usize::from(u16::from_le_bytes(jump_table[8..10].try_into().unwrap())) + tbl4;

            let offsets = [tbl1, tbl2, tbl3, tbl4, end];

            // read code lengths
            code_lengths[1..].copy_from_slice(&src[src_position..src_position + 11]);
            src_position += 11;
            //src.read_exact(&mut code_lengths[1..])?;

            let codes = code_lengths.iter().map(|x| *x as usize).sum::<usize>();
            if codes > 256
            {
                return Err(EntropyErrors::CorruptHeader(
                    "Code lengths are not less than 256".to_string(),
                ));
            }
            // read symbols
            symbols = &src[src_position..src_position + codes];
            src_position += codes;

            // Build the Huffman tree
            build_tree(&mut tbl, &code_lengths, symbols)?;

            let huff_source = &src[src_position..src_position + end];
            src_position += end;
            // new length
            let start = dest.len();

            let new_len = dest.len() + block_length as usize;
            // UNSAFE
            extend(dest,new_len);

            decode_symbols(
                huff_source,
                &tbl,
                &offsets,
                block_length as usize,
                // write to uninitialized memory :)
                dest.get_mut(start..).unwrap(),
            )?;
        }
        // top bit in block info indicates if block is the last
        // block.
        if (block_info >> 5 & 1) == 1
        {
            // end of block
            break;
        }
        // not last block, pull in more bytes.
        block_info = src[src_position];
        src_position += 1;

        // read the length for the next iteration
        length[0..3].copy_from_slice(&src[src_position..src_position + 3]);
        src_position += 3;
    }
    Ok(())
}
