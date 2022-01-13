//! This module provides huffman encoding and decoding routines

use std::io::{Cursor, Read};

use crate::bitstream::{BitStreamReader, Flags};

//https://play.rust-lang.org/?version=stable&mode=debug&edition=2021&gist=8296911b0e1ac889c5255aa91f5c39fc

pub const LIMIT: usize = 11;
/// Controls how many multi_decode symbols we can have.
///
/// It SHOULD REMAIN @ 2 for the best asm output.
const MAX_MULTI_DECODE: usize = 2;

/// Controls the minimum amount of bits allowed in the lookup table
const MIN_BITS: usize = 0;

#[derive(Copy, Clone, Default, Debug)]
pub struct DoubleEntry
{
    /*
    * This is just another trick.
    * The entry is the same size as a u64 so it's transmuted
    * To that layout, hence getting bits consumed and symbols resolved
    * Is simply shifts and not ugly loads

    * Also there is some weird stuff if this is shortened to u32,
    * the compiler emits spills, and I don't know how to fix them
    */
    /// symbols represent actual Huffman symbols
    /// which can be between 1 and 4 resolved depending on the bit entry
    pub symbols: [u8; MAX_MULTI_DECODE],
    // how many bits the actual symbols consume.
    pub bits_consumed: u8,
    // how many symbols were resolved during decoding
    pub symbols_resolved: u8,
}

#[derive(Copy, Clone, Default)]
pub struct SingleEntry
{
    pub symbol: u8,
    pub bits_consumed: u8,
}
// https://godbolt.org/z/9TP7jsnxn

#[derive(Copy, Clone, Default)]
struct Lookup
{
    bits_consumed: u8,
    symbol: u8,
}

pub struct HuffmanSingleDecompTable
{
    // this table is 8K bytes in memory
    pub table: [SingleEntry; 1 << LIMIT],
}

impl HuffmanSingleDecompTable
{
    /// Create a new Huffman Decompression instance
    fn new(code_lengths: &[u8; LIMIT], symbols: &[u8]) -> HuffmanSingleDecompTable
    {
        let mut tbl = HuffmanSingleDecompTable {
            table: [SingleEntry::default(); 1 << LIMIT],
        };
        tbl.build_tree(code_lengths, symbols);
        tbl
    }

    pub fn build_tree(&mut self, code_lengths: &[u8; LIMIT], symbols: &[u8])
    {
        const LIMIT_POWERED: usize = (1 << LIMIT) - 1;
        // An example is at @ https://play.rust-lang.org/?version=stable&mode=debug&edition=2021&gist=98631850955fbe2b236d28105b0e37aa

        /*
         * Build the symbol tree from code lengths
         * This is a nifty technique from jpeg where codes
         * Are actually built from code-lengths so there
         * is no need to transmit the codes themselves.
         * But since we have a hard limit on codes 12 bits,
         * we can do further stuff and  create a multi symbol decode
         * because we're amazing.
         */

        let mut code = 0;

        let mut p = 0;

        // Build the lookup table itself
        for i in 0..=LIMIT
        {
            for _ in 0..code_lengths[i]
            {
                // i-> its current code length
                // code -> the current code.
                let look_bits = ((code) as usize) << (LIMIT - i);

                for k in 0..(1 << (LIMIT - i))
                {
                    let entry = &mut self.table[(look_bits + k) & LIMIT_POWERED];
                    entry.bits_consumed = i as u8;
                    entry.symbol = symbols[p];
                }

                code += 1;

                p += 1;
            }
            code *= 2;
        }
    }
    /// Read Huffman Headers and parse output
    pub(crate) fn read_headers(buf: &mut Cursor<&[u8]>)
        -> std::io::Result<HuffmanSingleDecompTable>
    {
        // code lengths
        let mut code_lengths = [0; LIMIT];
        buf.read_exact(&mut code_lengths)?;
        // sum code lengths to get code sum.
        let sum_symbols = code_lengths.iter().map(|x| *x as usize).sum::<usize>();
        if sum_symbols > 256
        {
            // corrupt instance
            // Todo: Change this to an error when we have an error instance.
            panic!("Corrupt Huffman Header");
        }
        // actual symbols
        let mut codes = vec![0; sum_symbols];

        buf.read_exact(&mut codes)?;

        return Ok(HuffmanSingleDecompTable::new(&code_lengths, &codes));
    }
}

pub struct HuffmanDoubleDecompTable
{
    // this table is 16K bytes in memory
    pub table: [DoubleEntry; 1 << LIMIT],
}

impl HuffmanDoubleDecompTable
{
    fn new(code_lengths: &[u8; LIMIT], symbols: &[u8]) -> HuffmanDoubleDecompTable
    {
        let mut tbl = HuffmanDoubleDecompTable {
            table: [DoubleEntry::default(); 1 << LIMIT],
        };
        tbl.build_tree(code_lengths, symbols);
        tbl
    }
    /// Read Huffman Headers and parse output
    pub(crate) fn read_headers(buf: &mut Cursor<&[u8]>)
        -> std::io::Result<HuffmanDoubleDecompTable>
    {
        // code lengths
        let mut code_lengths = [0; LIMIT];
        buf.read_exact(&mut code_lengths)?;
        // sum code lengths to get code sum.
        let sum_symbols = code_lengths.iter().map(|x| *x as usize).sum::<usize>();
        if sum_symbols > 256
        {
            // corrupt instance
            // Todo: Change this to an error when we have an error instance.
            panic!("Corrupt Huffman Header");
        }
        // actual symbols
        let mut codes = vec![0; sum_symbols];
        // read the symbols.
        buf.read_exact(&mut codes)?;

        return Ok(HuffmanDoubleDecompTable::new(&code_lengths, &codes));
    }
    fn build_tree(&mut self, code_lengths: &[u8; LIMIT], symbols: &[u8])
    {
        const LIMIT_POWERED: usize = (1 << LIMIT) - 1;
        // An example is at @ https://play.rust-lang.org/?version=stable&mode=debug&edition=2021&gist=98631850955fbe2b236d28105b0e37aa

        /*
         * Build the symbol tree from code lengths
         * This is a nifty technique from jpeg where codes
         * Are actually built from code-lengths so there
         * is no need to transmit the codes themselves.
         * But since we have a hard limit on codes 12 bits,
         * we can do further stuff and  create a multi symbol decode
         * because we're amazing.
         */

        let mut code = 1;

        let mut p = 0;

        let mut temp_table = [Lookup::default(); 1 << LIMIT];
        // Build the lookup table itself
        for i in 0..LIMIT
        {
            for _ in 0..code_lengths[i]
            {
                // i-> its current code length
                // code -> the current code.
                let look_bits = ((code) as usize) << (LIMIT - i);

                for k in 0..(1 << (LIMIT - i))
                {
                    let entry = &mut self.table[(look_bits + k) & LIMIT_POWERED];
                    entry.bits_consumed = i as u8;
                    entry.symbols[0] = symbols[p];

                    let mut tbl = &mut temp_table[(look_bits + k) & LIMIT_POWERED];
                    tbl.bits_consumed = i as u8;
                    tbl.symbol = symbols[p]
                }

                code += 1;

                p += 1;
            }
            code *= 2;
        }

        // do the actual table filling
        // The first pass filled in top bits, now here we add multi symbol
        // decode

        /*
         * Okay here it gets complicated
         *
         * For multi-symbol decode, we use a referral technique,  initially we created a temporary
         * table(temp_table) which contains entries for a single table decode.
         * But if an entry consumes 5 bits we have 7 bits wasted in the loop. so here we try to
         * fill 7 bits with another bit. Yann Collet has a wonderful explanation here
         * https://fastcompression.blogspot.com/2015/10/huffman-revisited-part-4-multi-bytes.html
         *
         *
         * Basically we check if we first have enough space in our first temp_table to store additional bits
         * e.g if we found out that a symbol consumes 4 bits, it means the additional  7 bits(11-4)
         * are unused, so from there we simply look at the bit sequence, remove the bits already consumed
         * and then check if the sum of the current partially vacant entry and this new entry are below the limit
         * if they are, we can add this new symbol.
         *
         * Let's walk through one to make sense
         *
         * stream      =  0000_0001
         * symbol      = 0
         * code_length = 2
         * hence the remaining unused bits are (MSB fashion)
         * padded_stream      =    00_0001 ( First two zeroes are gone)
         *
         * pad streams with zeroes by shifting up until its the same size as stream(8)
         * stream_padded  = 0000_10_00 (two last zeroes are new)
         *
         * check if stream_padded resolves to a symbol
         * new_symbol       = temp_table[stream_padded].symbol
         * new_code_length  = temp_table[stream_padded].code_length
         *
         * if the sum of the new code length and initial symbol is less than LIMIT, create a new entry
         * entry[stream] =
         * {
         *  symbols:[symbol1,symbol2]
         *  bits_consumed:code_length+new_code_length,
         *  bytes:2
         * }
         */

        for i in 0..(1 << LIMIT)
        {
            let copy = temp_table[i];

            let mut entry = &mut self.table[i];

            let remaining_bits = LIMIT - usize::from(copy.bits_consumed);

            if (remaining_bits) > MIN_BITS
            {
                // We have space for another symbol
                let bits = i & ((1 << remaining_bits) - 1);

                let value = temp_table[(bits << copy.bits_consumed) & LIMIT_POWERED];
                // Fill entry with another symbol.

                if usize::from(entry.bits_consumed + value.bits_consumed) < LIMIT
                    && value.bits_consumed != 0
                {
                    // Increment bits consumed in symbol
                    entry.bits_consumed += value.bits_consumed;
                    // put the new symbol
                    entry.symbols[1] = value.symbol;
                    // increment symbols reserved
                    entry.symbols_resolved += 1;
                }
            }
        }
    }
}

pub fn decompress_huff_double_symbols(buf: &mut Cursor<&[u8]>)
{
    //https://godbolt.org/z/xchM5cY5q

    // read headers
    // 1. Headers start with 0xFF 0xDA to indicate start of stream
    // 2. The next 3 bytes indicate how long this decompressed block is.
    // 3. The ne
    // After that we have Huffman symbols.
    // then the data
    let mut huff_header = [0; 14];

    buf.read_exact(&mut huff_header).unwrap();
    if &huff_header[0..2] != &[0xFF, 0xDA]
    {
        panic!("Invalid Huffman header");
    }

    // The total block size is stored in the next 3 bytes, hence a block can be from 0..2^24 bytes
    let mut block_size = u32::from_le_bytes(huff_header[2..6].try_into().unwrap()) as usize;

    // remove the last byte
    block_size >>= 8;

    let huff_table = HuffmanDoubleDecompTable::read_headers(buf).unwrap();

    let stride = (block_size + 3) / 4;
    let mut bsize = stride;

    // get the remaining bits
    let block = &buf.get_ref()[buf.position() as usize..];

    // initialize streams with pointers to compressed data.
    let mut stream1 = BitStreamReader::new(&block[0..stride]);

    let mut stream2 = BitStreamReader::new(&block[bsize..bsize + stride]);
    bsize += stride;
    let mut stream3 = BitStreamReader::new(&block[bsize..bsize + stride]);
    bsize += stride;
    let mut stream4 = BitStreamReader::new(&block[bsize..]);

    // Variable to tell us to end the main decoder loop
    let mut end_signal = true;
    // destination is block size
    let mut destination = vec![0; block_size];

    // we know how long the dest block is, so to get how long each stream is
    // we can divide by 4.

    let stream_size = block_size >> 2;

    let (dest1, rem) = destination.split_at_mut(stream_size);
    let (dest2, rem) = rem.split_at_mut(stream_size);
    let (dest3, dest4) = rem.split_at_mut(stream_size);

    let mut offset1 = 0;
    let mut offset2 = 0;
    let mut offset3 = 0;
    let mut offset4 = 0;

    let mut entry1;
    let mut entry2;
    let mut entry3;
    let mut entry4;

    let mut buf1 = 0;
    let mut bits1 = 0;

    let mut buf2 = 0;
    let mut bits2 = 0;

    let mut buf3 = 0;
    let mut bits3 = 0;

    let mut buf4 = 0;
    let mut bits4 = 0;
    unsafe {
        end_signal |= stream1.refill_fast(&mut buf1, &mut bits1) == Flags::Normal;

        end_signal |= stream2.refill_fast(&mut buf2, &mut bits2) == Flags::Normal;

        end_signal |= stream3.refill_fast(&mut buf3, &mut bits3) == Flags::Normal;

        end_signal |= stream4.refill_fast(&mut buf4, &mut bits4) == Flags::Normal;
    }
    macro_rules! decode_huff_fast_multi_double {
        ($entry:tt,$reg:tt,$bits_remaining:tt,$dest:tt,$dest_position:tt) => {
            // peek bits
            $entry = huff_table.table[($reg >> (64 - LIMIT)) as usize];
            let temp = ($entry).bits_consumed as u8;

            // drop bits
            $reg <<= temp;
            $bits_remaining -= temp;

            // write the symbols
            let symbols_resolved = ($entry).symbols_resolved as usize;
            // copy 2 symbols
            ($dest)
                .as_mut_ptr()
                .add(($dest_position))
                .copy_from(($entry).symbols.as_ptr(), MAX_MULTI_DECODE);
            // but update position to how many bits were consumed
            ($dest_position) += symbols_resolved;
        };
    }
    unsafe {
        while end_signal
        {
            // the inline version generates too many register spills
            // and reloads, which are expensive(not tested), so I prefer the macro version
            // here.
            //

            decode_huff_fast_multi_double!(entry1, buf1, bits1, dest1, offset1);
            decode_huff_fast_multi_double!(entry2, buf2, bits2, dest2, offset2);
            decode_huff_fast_multi_double!(entry3, buf3, bits3, dest3, offset3);
            decode_huff_fast_multi_double!(entry4, buf4, bits4, dest4, offset4);

            decode_huff_fast_multi_double!(entry1, buf1, bits1, dest1, offset1);
            decode_huff_fast_multi_double!(entry2, buf2, bits2, dest2, offset2);
            decode_huff_fast_multi_double!(entry3, buf3, bits3, dest3, offset3);
            decode_huff_fast_multi_double!(entry4, buf4, bits4, dest4, offset4);

            decode_huff_fast_multi_double!(entry1, buf1, bits1, dest1, offset1);
            decode_huff_fast_multi_double!(entry2, buf2, bits2, dest2, offset2);
            decode_huff_fast_multi_double!(entry3, buf3, bits3, dest3, offset3);
            decode_huff_fast_multi_double!(entry4, buf4, bits4, dest4, offset4);

            decode_huff_fast_multi_double!(entry1, buf1, bits1, dest1, offset1);
            decode_huff_fast_multi_double!(entry2, buf2, bits2, dest2, offset2);
            decode_huff_fast_multi_double!(entry3, buf3, bits3, dest3, offset3);
            decode_huff_fast_multi_double!(entry4, buf4, bits4, dest4, offset4);

            decode_huff_fast_multi_double!(entry1, buf1, bits1, dest1, offset1);
            decode_huff_fast_multi_double!(entry2, buf2, bits2, dest2, offset2);
            decode_huff_fast_multi_double!(entry3, buf3, bits3, dest3, offset3);
            decode_huff_fast_multi_double!(entry4, buf4, bits4, dest4, offset4);

            // Refill buffer.

            end_signal |= stream1.refill_fast(&mut buf1, &mut bits1) == Flags::Normal;

            end_signal |= stream2.refill_fast(&mut buf2, &mut bits2) == Flags::Normal;

            end_signal |= stream3.refill_fast(&mut buf3, &mut bits3) == Flags::Normal;

            end_signal |= stream4.refill_fast(&mut buf4, &mut bits4) == Flags::Normal;
        }
    }
}

pub fn decompress_huff_single_symbols(buf: &mut Cursor<&[u8]>)
{
    // read headers
    // headers start with 0xFF 0xDA to indicate start of stream
    // The next 4 bytes indicate how long this block is
    // The final 12 bytes tell us bitstream offsets.
    // After that we have Huffman code sizes and code lengths
    // then the data
    let mut huff_header = [0; 14];

    buf.read_exact(&mut huff_header).unwrap();
    if &huff_header[0..2] != &[0xFF, 0xDA]
    {
        panic!("Invalid Huffman header");
    }

    // The total block size is stored in the next 3 bytes, hence a block can be from 0..2^24 bytes
    let mut block_size = u32::from_le_bytes(huff_header[2..6].try_into().unwrap()) as usize;

    // remove the last byte
    block_size >>= 8;

    let huff_table = HuffmanSingleDecompTable::read_headers(buf).unwrap();

    let stride = (block_size + 3) / 4;
    let mut bsize = stride;

    // After Huffman tables what follows is the block of data.,
    let block = &buf.get_ref()[buf.position() as usize..];

    // initialize streams with pointers to compressed data.
    let mut stream1 = BitStreamReader::new(&block[0..stride]);

    let mut stream2 = BitStreamReader::new(&block[bsize..bsize + stride]);

    bsize += stride;
    let mut stream3 = BitStreamReader::new(&block[bsize..bsize + stride]);

    bsize += stride;
    let mut stream4 = BitStreamReader::new(&block[bsize..]);

    // Variable to tell us to end the main decoder loop
    let mut end_signal = true;
    // destination is block size
    let mut destination = vec![0; block_size];

    // we know how long the dest block is, so to get how long each stream is
    // we can divide by 4.

    let stream_size = (block_size + 3) >> 2;

    let (dest1, rem) = destination.split_at_mut(stream_size);
    let (dest2, rem) = rem.split_at_mut(stream_size);
    let (dest3, dest4) = rem.split_at_mut(stream_size);

    let mut entry1;
    let mut entry2;
    let mut entry3;
    let mut entry4;

    let mut buf1 = 0;
    let mut bits1 = 0;

    let mut buf2 = 0;
    let mut bits2 = 0;

    let mut buf3 = 0;
    let mut bits3 = 0;

    let mut buf4 = 0;
    let mut bits4 = 0;
    unsafe {
        end_signal |= stream1.refill_fast(&mut buf1, &mut bits1) == Flags::Normal;

        end_signal |= stream2.refill_fast(&mut buf2, &mut bits2) == Flags::Normal;

        end_signal |= stream3.refill_fast(&mut buf3, &mut bits3) == Flags::Normal;

        end_signal |= stream4.refill_fast(&mut buf4, &mut bits4) == Flags::Normal;
    }
    macro_rules! decode_huff_fast_multi_single {
        ($entry:tt,$reg:tt,$bits_remaining:tt,$dest:tt) => {
            // peek bits
            $entry = huff_table.table[($reg >> (64 - LIMIT)) as usize];
            let temp = ($entry).bits_consumed as u8;
            // drop bits
            $reg <<= temp;
            $bits_remaining -= temp;

            // copy 1 symbol
            ($dest) = $entry.symbol;
        };
    }
    unsafe {
        for (((a, b), c), d) in dest1
            .chunks_exact_mut(5)
            .zip(dest2.chunks_exact_mut(5))
            .zip(dest3.chunks_exact_mut(5))
            .zip(dest4.chunks_exact_mut(5))
        {
            decode_huff_fast_multi_single!(entry1, buf1, bits1, (a[0]));
            decode_huff_fast_multi_single!(entry2, buf2, bits2, (b[0]));
            decode_huff_fast_multi_single!(entry3, buf3, bits3, (c[0]));
            decode_huff_fast_multi_single!(entry4, buf4, bits4, (d[0]));

            decode_huff_fast_multi_single!(entry1, buf1, bits1, (a[1]));
            decode_huff_fast_multi_single!(entry2, buf2, bits2, (b[1]));
            decode_huff_fast_multi_single!(entry3, buf3, bits3, (c[1]));
            decode_huff_fast_multi_single!(entry4, buf4, bits4, (d[1]));

            decode_huff_fast_multi_single!(entry1, buf1, bits1, (a[2]));
            decode_huff_fast_multi_single!(entry2, buf2, bits2, (b[2]));
            decode_huff_fast_multi_single!(entry3, buf3, bits3, (c[2]));
            decode_huff_fast_multi_single!(entry4, buf4, bits4, (d[2]));

            decode_huff_fast_multi_single!(entry1, buf1, bits1, (a[3]));
            decode_huff_fast_multi_single!(entry2, buf2, bits2, (b[3]));
            decode_huff_fast_multi_single!(entry3, buf3, bits3, (c[3]));
            decode_huff_fast_multi_single!(entry4, buf4, bits4, (d[3]));

            decode_huff_fast_multi_single!(entry1, buf1, bits1, (a[4]));
            decode_huff_fast_multi_single!(entry2, buf2, bits2, (b[4]));
            decode_huff_fast_multi_single!(entry3, buf3, bits3, (c[4]));
            decode_huff_fast_multi_single!(entry4, buf4, bits4, (d[4]));

            // Refill buffer.

            end_signal |= stream1.refill_fast(&mut buf1, &mut bits1) == Flags::Normal;

            end_signal |= stream2.refill_fast(&mut buf2, &mut bits2) == Flags::Normal;

            end_signal |= stream3.refill_fast(&mut buf3, &mut bits3) == Flags::Normal;

            end_signal |= stream4.refill_fast(&mut buf4, &mut bits4) == Flags::Normal;

            if !end_signal
            {
                // something happened.
                break;
            }
        }
    }
}
