//! This module provides huffman encoding and decoding routines

use std::io::{Read, Write};

use crate::{bitstream::BitStreamReader, utils::REVERSED_BITS};
pub const LIMIT: usize = 11;

#[derive(Copy, Clone, Debug, Default)]
pub struct SingleEntry
{
    pub symbol: u8,
    pub bits_consumed: u8,
}

struct HuffmanSingleDecompTable<'a>
{
    // this table is 8K bytes in memory
    pub table: &'a mut [u16;1<<LIMIT] ,
}

impl<'a> HuffmanSingleDecompTable<'a>
{
    /// Create a new Huffman Decompression instance
    fn new(
        code_lengths: &[u8; LIMIT + 1], symbols: &[u8], table:&'a mut [u16;1<<LIMIT],
    ) -> HuffmanSingleDecompTable<'a>
    {
        let mut tbl = HuffmanSingleDecompTable {
            table,
        };
        tbl.build_tree(code_lengths, symbols);
        tbl
    }

    pub fn build_tree(
        &mut self, code_lengths: &[u8; LIMIT + 1], symbols: &[u8],
    )
    {
        let mut code = 0;

        let mut p = 0;

        for i in 1..=LIMIT
        {
            for _ in 0..code_lengths[i]
            {
                let look_bits = code << (LIMIT - i);

                for k in 0..(1 << (LIMIT - i))
                {
                    let entry = &mut self.table
                        [(REVERSED_BITS[(look_bits + k) as usize] >> (16 - LIMIT)) as usize];


                    *entry |= u16::from(symbols[p])<<8 |(i as u16);
                }
                p += 1;

                code += 1;
            }
            code *= 2;
        }
    }
}

fn decompress_huff_inner(
    buf: &[u8], table: &HuffmanSingleDecompTable, offsets: &[usize; 5], block_size: usize,
    dest: &mut [u8],
)
{
    // initialize streams with pointers to compressed data.
    let mut stream1 = BitStreamReader::new(&buf[0..offsets[0]]);

    let mut stream2 = BitStreamReader::new(&buf[offsets[0]..offsets[1]]);

    let mut stream3 = BitStreamReader::new(&buf[offsets[1]..offsets[2]]);

    let mut stream4 = BitStreamReader::new(&buf[offsets[2]..offsets[3]]);

    let mut stream5 = BitStreamReader::new(&buf[offsets[3]..offsets[4]]);

    let stream_size = (block_size + 4) / 5;

    let (dest1, rem) = dest.split_at_mut(stream_size);
    let (dest2, rem) = rem.split_at_mut(stream_size);
    let (dest3, rem) = rem.split_at_mut(stream_size);
    let (dest4, dest5) = rem.split_at_mut(stream_size);

    unsafe {
        let entries = &table.table;

        for ((((a, b), c), d), e) in dest1
            .chunks_exact_mut(20)
            .zip(dest2.chunks_exact_mut(20))
            .zip(dest3.chunks_exact_mut(20))
            .zip(dest4.chunks_exact_mut(20))
            .zip(dest5.chunks_exact_mut(20))
        {
            // check that there is at least 11*20(220 ) bytes in the buffer
            if !stream1.check_first()
                || !stream2.check_first()
                || !stream3.check_first()
                || !stream4.check_first()
            {
                // use the safer decoder loops
                break;
            }
            // decode bytes, up to twenty symbols per loop.
            macro_rules! decode_single {
                ($a:tt,$b:tt) => {
                    stream1.decode_multi(a.get_mut($a..$b).unwrap().try_into().unwrap(), entries);
                    stream2.decode_multi(b.get_mut($a..$b).unwrap().try_into().unwrap(), entries);
                    stream3.decode_multi(c.get_mut($a..$b).unwrap().try_into().unwrap(), entries);
                    stream4.decode_multi(d.get_mut($a..$b).unwrap().try_into().unwrap(), entries);
                    stream5.decode_multi(e.get_mut($a..$b).unwrap().try_into().unwrap(), entries);
                };
            }

            decode_single!(0, 5);
            decode_single!(5, 10);
            decode_single!(10, 15);
            decode_single!(15, 20);
        }
    }
}

pub fn huff_decompress_4x<R: Read, W: Write>(src: &mut R, dest: &mut W)
{
    let mut length = [0, 0, 0, 0];
    // read the length
    src.read_exact(&mut length[0..3]).unwrap();

    let mut block_length = u32::from_le_bytes(length);

    // not all bytes will be used
    let mut source = vec![0; block_length as usize];

    // we know that blocks have equal sizes except the last one
    // so if we read the first one we can determine how much space we will need
    let mut dest_temp = vec![0; block_length as usize];
    
    let mut tbl = [0;1<<LIMIT];
    
    loop
    {
        block_length = u32::from_le_bytes(length);

        // read checksum
        let mut checksum = [0; 3];

        src.read_exact(&mut checksum).unwrap();

        // read jump table
        let mut jump_table = [0; 10];

        src.read_exact(&mut jump_table).unwrap();

        // two bytes per jump table, stored in little endian form.
        let tbl1 = u32::from(jump_table[0]) + (u32::from(jump_table[1]) << 8);

        let tbl2 = u32::from(jump_table[2]) + (u32::from(jump_table[3]) << 8) + tbl1;

        let tbl3 = u32::from(jump_table[4]) + (u32::from(jump_table[5]) << 8) + tbl2;

        let tbl4 = u32::from(jump_table[6]) + (u32::from(jump_table[7]) << 8) + tbl3;

        // end of this stream.
        let end = u32::from(jump_table[8]) + (u32::from(jump_table[9]) << 8) + tbl4;

        let offsets = [tbl1, tbl2, tbl3, tbl4, end].map(|x| x as usize);

        //read code lengths
        let mut code_lengths = [0; 12];

        src.read_exact(&mut code_lengths[1..]).unwrap();

        let codes = code_lengths.iter().map(|x| *x as usize).sum::<usize>();

        let mut symbols = vec![0; codes];

        src.read_exact(&mut symbols).unwrap();

        let huff_table = HuffmanSingleDecompTable::new(&code_lengths, &symbols, &mut tbl);

        let huff_source = &mut source[0..end as usize];

        src.read_exact(huff_source).unwrap();

        // now we have offsets, streams and everything we need to decode. Let's go
        decompress_huff_inner(
            huff_source,
            &huff_table,
            &offsets,
            block_length as usize,
            &mut dest_temp,
        );
        // copy dest_temp to writer
        dest.write_all(&dest_temp[0..block_length as usize])
            .unwrap();

        // read the length for the next iteration
        if src.read_exact(&mut length[0..3]).is_err()
        {
            break;
        }
    }
}

#[test]
fn huff_decompress()
{
    use std::fs::OpenOptions;
    use std::io::{BufReader, BufWriter};

    let fs = OpenOptions::new()
        .create(false)
        .write(false)
        .read(true)
        .open("/Users/calebe/CLionProjects/zcif/tests.zcif")
        .unwrap();
    let mut fs = BufReader::with_capacity(1 << 24, fs);

    let fd = OpenOptions::new()
        .read(true)
        .open("/Users/calebe/git/FiniteStateEntropy/programs/enwiki.small")
        .unwrap();
    let mut fd = BufWriter::new(fd);

    huff_decompress_4x(&mut fs, &mut fd);
}
