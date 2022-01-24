//! BitStreamReader API
//!
//! This module provides an interface to read and write bits (and bytes)

use crate::constants::LIMIT;

pub struct BitStreamReader<'src>
{
    // buffer from which we are pulling in bits from
    // used in decompression.
    src: &'src [u8],
    // position in our buffer,
    position: usize,

    bits_left: u8,
    buffer: u64,
}

impl<'src> BitStreamReader<'src>
{
    /// Create a new BitStreamReader instance
    ///
    /// # Expectations
    /// The buffer must be padded with fill bytes in the end,
    /// if not, this becomes UB in the refill phase.
    pub fn new(in_buffer: &'src [u8]) -> BitStreamReader<'src>
    {
        BitStreamReader {
            bits_left: 0,
            buffer: 0,
            src: in_buffer,
            position: 0,
        }
    }
    /// Refill the bitstream ensuring the buffer has bits between
    /// 56 and 63.
    ///
    #[inline(always)]
    pub unsafe fn refill_fast(&mut self)
    {
        /*
         * The refill always guarantees refills between 56-63
         * This version is a modification of Variant 4 (tho I had come up with
         * A similar but not so optimal method) found in Fabian's Giesen
         * Reading bits in far too many ways.
         * @ https://fgiesen.wordpress.com/2018/02/20/reading-bits-in-far-too-many-ways-part-2/
         *
         * Bits stored will never go above 63 and if bits are in the range 56-63 no refills occur.
         */

        let mut buf = [0; 8];

        // position points to the end initially, so subtracting means we
        // are reading the buffer from end to start.
        std::ptr::copy_nonoverlapping(self.src.as_ptr().add(self.position), buf.as_mut_ptr(), 8);

        // create a u64 from an array of u8's
        let new_buffer = u64::from_le_bytes(buf);
        // num indicates how many bytes we actually consumed.
        let num = 63 - self.bits_left;
        // offset position
        self.position += (num >> 3) as usize;
        // shift number of bits
        self.buffer |= new_buffer << self.bits_left;

        // update bits left
        // bits left are now between 56-63
        self.bits_left |= 56;
    }

    #[inline(always)]
    pub const fn peek_bits<const LOOKAHEAD: usize>(&self) -> usize
    {
        (self.buffer & ((1 << LOOKAHEAD) - 1)) as usize
    }
    /// Decode a single symbol
    #[inline(always)]
    pub fn decode_single(&mut self, dest: &mut u8, table: &[u16; (1 << LIMIT)])
    {
        let entry = table[(self.peek_bits::<LIMIT>())];

        let bits = (entry & 0xFF) as u8;

        // remove bits read.
        self.buffer >>= bits;

        self.bits_left -= bits;
        // write to position
        *dest = (entry >> 8) as u8;
    }

    // Check that we didn't read past our buffer
    pub fn check_final(&self) -> bool
    {
        /*
         * Here we take into account the position, and actual bits consumed
         * since we may overshoot so we need to account for the bytes
         */

        self.position - (usize::from(self.bits_left >> 3)) == self.src.len()
    }
    /// Get current position of the inner buffer
    pub fn get_position(&self) -> usize
    {
        // take into account the bytes not consumed in the current
        // bitstream.
        self.position - (usize::from(self.bits_left >> 3))
    }
    /// Get the length  of the inner buffer
    pub fn get_src_len(&self) -> usize
    {
        self.src.len()
    }
}

/// A compact bit writer for the Huffman encoding.
/// algorithm.
pub struct BitStreamWriter
{
    // Number of actual bits in the bit buffer.
    bits_in_buffer: u8,
    // should I use usize?
    buf: u64,
    // position to write this in the output buffer
    position: usize,
}

impl BitStreamWriter
{
    pub fn new() -> BitStreamWriter
    {
        BitStreamWriter {
            buf: 0,
            bits_in_buffer: 0,
            position: 0,
        }
    }
    /// Deal with initial bits which can't be handled by the fast path.
    pub fn write_bits_slow(&mut self, symbols: &[u8], entries: &[u32; 256], out_buf: &mut [u8])
    {
        let mut flush_bit = 0;

        for symbol in symbols
        {
            let entry = entries[usize::from(*symbol)];

            // add to the top bits
            self.buf |= u64::from(entry >> 8) << self.bits_in_buffer;

            self.bits_in_buffer += (entry & 0xFF) as u8;

            flush_bit += 1;

            if flush_bit == 5
            {
                flush_bit = 0;
                unsafe {
                    self.flush_fast(out_buf);
                }
            }
        }
        // flush any left
        unsafe {
            self.flush_fast(out_buf);
        }
    }
    #[inline(always)]
    pub unsafe fn write_bits_fast(
        &mut self, symbols: &[u8; 5], entry: &[u32; 256], out_buf: &mut [u8],
    )
    {
        /*
         *The limit is 11 bits per symbol, therefore we can go
         * to 5 symbols per encode (55) before flushing.
         *
         * The bits are added in a fifo manner with the bits representing the first
         * symbol
         *
         */
        macro_rules! encode_single {
            ($pos:tt) => {
                let entry = entry[symbols[$pos] as usize];

                // add to the top bits
                self.buf |= (u64::from(entry >> 8) << self.bits_in_buffer);

                self.bits_in_buffer += (entry & 0xFF) as u8;
            };
        }
        // TODO: Benchmarks report faster speeds when using
        // movzxd than when using the earlier shift and get, benchmark
        // on linux to see.
        encode_single!(0);

        encode_single!(1);

        encode_single!(2);

        encode_single!(3);

        encode_single!(4);

        // flush to output buffer
        self.flush_fast(out_buf);
    }

    /// Flush bits to the output buffer
    ///
    /// After calling this routine, the bit buffer is guaranteed
    /// to have less than 8 bits.
    #[inline(always)]
    unsafe fn flush_fast(&mut self, out_buf: &mut [u8])
    {
        /*
        * Most of this is from Eric Biggers libdeflate
        * @ https://github.com/ebiggers/libdeflate/blob/master/lib/deflate_compress.c
        *
        * Which has some nice properties , branchless writes high ILP etc etc.

        */

        let buf = self.buf.to_le_bytes();
        // write 8 bytes
        out_buf
            .as_mut_ptr()
            .add(self.position)
            .copy_from(buf.as_ptr(), 8);
        // but update position to point to the full number of symbols we read
        let bytes_written = self.bits_in_buffer & 56;
        // remove those bits we wrote.
        self.buf >>= bytes_written;
        // increment position
        self.position += (bytes_written >> 3) as usize;

        self.bits_in_buffer &= 7;
    }
    #[cold]
    pub unsafe fn flush_final(&mut self, out_buf: &mut [u8])
    {
        // shift buffer by 7, so that the remaining bits are in the range 8-15
        self.buf <<= 7;
        self.bits_in_buffer += 7;
        self.flush_fast(out_buf);
    }
    /// Set everything to zero.
    pub fn reset(&mut self)
    {
        self.bits_in_buffer = 0;
        self.buf = 0;
        self.position = 0;
    }
    pub fn get_position(&self) -> usize
    {
        self.position
    }
}
