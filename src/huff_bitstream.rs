//! `BitStreamReader` API
//!
//! This module provides an interface to read and write bits (and bytes) for
//! huffman

use crate::constants::{LIMIT, TABLE_SIZE};

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
    /// Create a new `BitStreamReader` instance
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
         *
         * Bits stored will never go above 63 and if bits are in the range 56-63 no refills occur.
         */

        let mut buf = [0; 8];

        std::ptr::copy_nonoverlapping(self.src.as_ptr().add(self.position), buf.as_mut_ptr(), 8);

        // create a u64 from an array of u8's
        let new_buffer = u64::from_le_bytes(buf);
        // num indicates how many bytes we actually consumed.
        let num = 63 ^ self.bits_left;
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
    #[cfg(not(all(target_arch = "x86_64", target_feature = "bmi2")))] // bmi support
    pub fn decode_single(&mut self, dest: &mut u8, table: &[u16; TABLE_SIZE])
    {
        let entry = table[(self.peek_bits::<LIMIT>())];

        let bits = (entry & 0xFF) as u8;

        // remove bits read.
        self.buffer >>= bits;
        //
        // Rust generates the worst codegen here , some weird instructions which
        // hurt performance,(888 Mb/s) , explicitly subtracting bumps up performance to a cool
        // 1010 Mb/s
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            // this forces the compiler to keep values in register
            unsafe {
                use std::arch::asm;
                asm!("sub {0},{1}", inout(reg_byte) self.bits_left, in(reg_byte) bits);
            }
        }
        #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
        {
            // other architectures, normal subtraction,
            self.bits_left -= bits;
        }

        // write to position
        *dest = (entry >> 8) as u8;
    }

    #[cfg(all(target_arch = "x86_64", target_feature = "bmi2"))]
    pub fn decode_single(&mut self, dest: &mut u8, table: &[u16; TABLE_SIZE])
    {
        /*
         * Generate better code for CPUs supporting BMI2,
         * this is a compile time constraint and should be enabled via
         * RUSTFLAGS ="-C target-features=+bmi2" cargo build --release to enable it.
         *
         * Performance difference is about 90 Mb/s(1100 Mb/s vs 1010 Mb/s) better compared to  the decode_single one
         * (AMD Ryzen 3600U)
         *
         * The main issue is that Rust spills values to memory even when not needed.
         * forcing use of asm tells Rust to keep them in registers and not memory,
         */
        let entry = table[(self.peek_bits::<LIMIT>())];

        unsafe {
            use std::arch::asm;

            // keep values in register Rust.
            asm!(
                // shrx instruction uses the lower 0-6 bits of the last(entry:r) register, that's why ther
                // is no explicit masking
                 "shrx {buf}, {buf}, {entry:r}", // self.buf >>= (entry & 0xFF);
                 "sub {bits_left},{entry:l}", // self.bits_left -= (entry & 0xFF);

                buf= inout(reg) self.buffer,
                entry = in(reg_abcd) entry,
                bits_left = inout(reg_byte) self.bits_left,

            );
        }
        *dest = (entry >> 8) as u8;
    }

    pub fn get_bits(&mut self, num_bits: u8) -> u64
    {
        let mask = (1_u64 << num_bits) - 1;

        let value = self.buffer & mask;

        self.buffer >>= num_bits;

        self.bits_left -= num_bits;

        value
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
}

/// A compact bit writer for the Huffman encoding.
/// algorithm.
pub struct BitStreamWriter<'dest>
{
    // Number of actual bits in the bit buffer.
    bits: u8,
    // should I use usize?
    buf: u64,
    // position to write this in the output buffer
    position: usize,
    // dest
    dest: &'dest mut [u8],
}

impl<'dest> BitStreamWriter<'dest>
{
    pub fn new(dest: &mut [u8]) -> BitStreamWriter
    {
        BitStreamWriter {
            buf: 0,
            bits: 0,
            position: 0,
            dest,
        }
    }
    /// Deal with initial bits which can't be handled by the fast path.
    pub fn write_bits_slow(&mut self, symbols: &[u8], entries: &[u32; 256])
    {
        let mut flush_bit = 0;

        for symbol in symbols
        {
            let entry = entries[usize::from(*symbol)];

            self.add_bits(u64::from(entry >> 8), (entry & 0xFF) as u8);

            flush_bit += 1;

            if flush_bit == 5
            {
                flush_bit = 0;
                unsafe {
                    self.flush_fast();
                }
            }
        }
        // flush any left
        unsafe {
            self.flush_fast();
        }
    }
    /// Add new bits into the bit buffer variable
    #[inline(always)]
    pub fn add_bits(&mut self, value: u64, nbits: u8)
    {
        self.buf |= value << self.bits;

        self.bits += nbits;
    }
    #[inline(always)]
    #[allow(clippy::trivially_copy_pass_by_ref)] // Passing by value is detrimental here
    pub unsafe fn encode_symbols(&mut self, symbols: &[u8; 5], entry: &[u32; 256])
    {
        /*
         * The limit is 11 bits per symbol, therefore we can go
         * to 5 symbols per encode (55) before flushing.
         *
         */
        macro_rules! encode_single {
            ($pos:tt) => {
                let entry = entry[symbols[$pos] as usize];

                // add to the top bits
                self.add_bits(u64::from(entry >> 8), (entry & 0xFF) as u8)
            };
        }

        encode_single!(0);

        encode_single!(1);

        encode_single!(2);

        encode_single!(3);

        encode_single!(4);

        // flush to output buffer
        self.flush_fast();
    }
    /// Flush bits to the output buffer
    ///
    /// After calling this routine, the bit buffer is guaranteed
    /// to have less than 8 bits.
    #[inline(always)]
    pub(crate) unsafe fn flush_fast(&mut self)
    {
        let buf = self.buf.to_le_bytes();
        // write 8 bytes
        self.dest
            .as_mut_ptr()
            .add(self.position)
            .copy_from(buf.as_ptr(), 8);
        // but update position to point to the full number of symbols we read
        let bytes_written = self.bits & 56;

        // remove those bits we wrote.
        self.buf >>= bytes_written;
        // increment position
        self.position += (bytes_written >> 3) as usize;

        self.bits &= 7;
    }

    #[cold]
    pub unsafe fn flush_final(&mut self)
    {
        let diff = (-i16::from(self.bits) & 7) as u8;

        self.bits += diff;

        self.flush_fast();
    }

    pub fn get_position(&self) -> usize
    {
        self.position
    }
    pub fn get_output(&self) -> &[u8]
    {
        &self.dest[0..self.position]
    }
}
