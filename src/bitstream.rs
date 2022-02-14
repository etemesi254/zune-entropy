//! BitStreamReader API
//!
//! This module provides an interface to read and write bits (and bytes)

use crate::constants::{LIMIT, TABLE_LOG, TABLE_SIZE};

use crate::utils::Symbols;

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
        /*
         * TODO: Is it better to use ctlz here?
         *
         * These are the benefits I can see:
         *
         * 1. We don't need to maintain bits_left variable(more available registers)
         * 2. We save a sub instruction in the main loop(more available registers)
         *
         * The problems on the other hand aren't good.
         *
         * 1. Increased complexity at this refill stage
         * 2. Requires me to switch to nightly to use unsafe lzcnt that is UB if 0(note can't be zero)
         * 3. Count leading zeroes have a latency of 3 hence this might be a tad slower(and its in the serial dependency).
         * 4. ILP in the inner loop favours sub instruction as it can execute in parallel.
         * 5. When we want Count leading zeroes, we want to be in 64-bit,but in 64 bit, we don't spill to stack space(
         *    I checked)
         *
         * -> 6 it's SLOW on Arm, reduces speed by 200Mb/s idk why
         *  TODO: Investigate on x86_64
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
    pub fn decode_single(&mut self, dest: &mut u8, table: &[u16; (1 << LIMIT)])
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

const MASK: [u16; 17] = [
    0, 1, 3, 7, 15, 31, 63, 127, 255, 511, 1023, 2047, 4095, 8191, 16383, 32767, 65535,
];

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
    pub unsafe fn encode_bits_huffman(
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

    #[inline(always)]
    pub unsafe fn encode_bits_fse(
        &mut self, symbol: u8, entries: &[Symbols; 256], next_states: &[u16; TABLE_SIZE],
        curr_state: &mut u16,
    )
    {
        let symbol = entries[usize::from(symbol)];

        // TODO: Make num bits be determined by a state threshold like
        // Cbloom FSE
        // TODO:See if I can shift down with a power of two, to determine
        // num_bits
        let num_bits = (symbol.x + (*curr_state as u32)) >> TABLE_LOG;

        // write bits
        self.buf |= u64::from(*curr_state & *MASK.get_unchecked(num_bits as usize))
            << self.bits_in_buffer;

        self.bits_in_buffer += num_bits as u8;

        // Determine next state
        let offset = (*curr_state >> num_bits) as u16;

        *curr_state = *next_states.get_unchecked((symbol.y + offset) as usize);
    }
    /// Flush bits to the output buffer
    ///
    /// After calling this routine, the bit buffer is guaranteed
    /// to have less than 8 bits.
    #[inline(always)]
    pub(crate) unsafe fn flush_fast(&mut self, out_buf: &mut [u8])
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
        // align the bit buffer to a byte

        // (why did it take me soo many hours to figure this out :-( )
        self.bits_in_buffer += (-i16::from(self.bits_in_buffer) & 7) as u8;

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
