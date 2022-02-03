//! BMI optimized bitstream reader.
//!
//! This is an optimized BMI2 variant of a BitStream reader with
//! super powers

use crate::constants::LIMIT;

pub struct BmiBitStreamReader<'src>
{
    // buffer from which we are pulling in bits from
    // used in decompression.
    src: &'src [u8],
    // position in our buffer,
    position: usize,
    // buffer containing current bits.
    buffer: u64,
}

impl <'src> BmiBitStreamReader<'src> {
    pub fn new(in_buffer: &'src [u8]) -> BmiBitStreamReader<'src>
    {
        BmiBitStreamReader {
            buffer: 0,
            src: in_buffer,
            position: 0,
        }
    }
    //#[cfg(target_feature = "lzcnt")]
    pub unsafe fn refill(&mut self){
        // count leading zeroes
        // leading zeroes tell us how many empty bits we have
        let mut bits_consumed = self.buffer.leading_zeros();

        let mut buf = [0; 8];

        std::ptr::copy_nonoverlapping(self.src.as_ptr().add(self.position), buf.as_mut_ptr(), 8);

        // create a u64 from an array of u8's
        let mut new_buffer = u64::from_le_bytes(buf);
        // the number of full bytes we can read
        self.position += (bits_consumed >> 3) as usize;
        // actual bits in the stream not consumed
        let actual_bits = bits_consumed ^ 63;
        // shift up so that we have new bits are added on top
        // of the old bits
        new_buffer<<=actual_bits;
        let marker = 1<<(bits_consumed+1);

        // or old bits and new bits
        self.buffer |= new_buffer;

        // add the marker
        self.buffer |= marker;

    }
    #[inline(always)]
    pub const fn peek_bits<const LOOKAHEAD: usize>(&self) -> usize
    {
        (self.buffer & ((1 << LOOKAHEAD) - 1)) as usize
    }


    #[inline(always)]
    #[cfg(all(target_arch = "x86_64", target_feature = "bmi2"))] // bmi support
    pub fn decode_single(&mut self, dest: &mut u8, table: &[u16; (1 << LIMIT)])
    {
        unsafe {
            use std::arch::asm;

            // keep values in register Rust.
            asm!(
            // shrx instruction uses the lower 0-6 bits of the last(entry:r) register, that's why ther
            // is no explicit masking
            "shrx {buf}, {buf}, {entry:r}", // self.buf >>= (entry & 0xFF);

            buf= inout(reg) self.buffer,
            entry = in(reg_abcd) entry,

            );
        }
    }

    // Check that we didn't read past our buffer
    pub fn check_final(&self) -> bool
    {
        /*
         * Here we take into account the position, and actual bits consumed
         * since we may overshoot so we need to account for the bytes
         */

        self.position - ((self.buffer.leading_zeros() >> 3) as usize) == self.src.len()
    }
    /// Get current position of the inner buffer
    pub fn get_position(&self) -> usize
    {
        // take into account the bytes not consumed in the current
        // bitstream.
        self.position - (usize::from((self.buffer.leading_zeros() >> 3) as usize))
    }
    /// Get the length  of the inner buffer
    pub fn get_src_len(&self) -> usize
    {
        self.src.len()
    }
}
