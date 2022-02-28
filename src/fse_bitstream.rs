//! FSE bit-readers and writers
//! 
//! This file features a standalone Bit/Io implementation for 
//! tANS/FSE bitreaders/bitwriters
//! 
//! It is Little Endian and entirely lifted from Eric Biggers xpack
//! (see https://github.com/ebiggers/xpack) and it's pretty fast.
//! (well it's the fastest I could come up with)
//! 
//! And as is with the Huffman , the BitStream takes care of encoding and decoding
//! inside itself. Encoding is simply passing a symbol to encode_symbol()
//! decoding happens when you call decode_symbol() (each with appropriate arguments)
//! 
//! Also all functions found in the inner encoder and decoder loops are inlined.
//! 
//! 
//! # Safety
//! The only unsafe function is the `flush_fast` and `refill_fast` since 
//! they both write/read 8 bytes to/from memory , if it happens that the memory region is
//! out of bounds, they will obviously be UB.
use crate::constants::{MAX_TABLE_LOG, TABLE_LOG, TABLE_SIZE};
use crate::utils::Symbols;
/// Compact FSE bit-stream writer.
pub struct FseStreamWriter<'dest>
{
    // Number of actual bits in the bit buffer.
    bits: u8,
    // Stores current unflushed bits
    buf: u64,
    // position to write this in the output buffer
    pub(crate) position: usize,

    dest: &'dest mut [u8],
}
impl<'dest> FseStreamWriter<'dest>
{
    /// Create a new stream writer for a FSE
    pub fn new(out_dest: &mut [u8]) -> FseStreamWriter
    {
        // start 8 bits from the end
        // if we are on a 32 bit arch probably change this.
        let position = out_dest.len() - (u64::BITS / u8::BITS) as usize;
        return FseStreamWriter {
            bits: 0,
            buf: 0,
            position,
            dest: out_dest,
        };
    }
    /// Add new bits to the buffer
    /// and update bits in the buffer.
    #[inline(always)]
    fn add_bits(&mut self, nbits: u8, value: u16)
    {
        // check that adding the value won't corrupt top bits
        // indicating value wasn't masked
        debug_assert!((nbits as u32) <= value.leading_zeros());

        // new bits are added to the lower bits of the bit buffer
        self.buf = (self.buf << nbits) | u64::from(value);

        self.bits += nbits as u8;
    }
    /// Encode symbols and update current state(curr_state) 
    /// to point to the next state.
    #[inline(always)]
    pub fn encode_symbol(
        &mut self, symbol: u8, entries: &[Symbols; 256], next_states: &[u16; TABLE_SIZE],
        curr_state: &mut u16,
    )
    {
        let symbol = entries[usize::from(symbol)];
        // How number of bits evolves is a bit tricky..
        let num_bits = ((symbol.x + (*curr_state as u32)) >> TABLE_LOG) as u8;

        let mask = (1 << num_bits) - 1;

        let low_bits = *curr_state & mask;

        self.add_bits(num_bits, low_bits);

        // Determine next state
        let offset = (*curr_state >> num_bits) as i16;

        *curr_state = next_states[((symbol.z + offset) as usize & (TABLE_SIZE - 1))];
    }

    /// Encode final values of states to the bitstream
    /// The decoder will read final states from the bitstream
    /// and work its way to the initial state
    pub fn encode_final_states(&mut self, c1: u16, c2: u16, c3: u16, c4: u16, c5: u16)
    {
        // this should be zero since we did a flush earlier
        assert_eq!(self.bits, 0);

        self.buf = 0;
        macro_rules! encode_single {
            ($state:tt,$bits:tt) => {
                self.add_bits($bits as u8, $state & ((1 << MAX_TABLE_LOG) - 1));
            };
        }
        encode_single!(c1, MAX_TABLE_LOG);
        encode_single!(c2, MAX_TABLE_LOG);
        encode_single!(c3, MAX_TABLE_LOG);
        encode_single!(c4, MAX_TABLE_LOG);
        encode_single!(c5, MAX_TABLE_LOG);

        // add a dummy bit so that we have 56 bits in the buffer
        self.add_bits(1, 0);

        unsafe {
            self.flush_fast();
        }

        // no bits left
        assert_eq!(self.bits, 0);
    }

    pub fn flush_final(&mut self)
    {
        assert!(self.bits <= 7);
        // terminate the encoder with a 1 bit so that the decoder
        // knows where to start

        self.buf <<= 8 - self.bits;
        self.buf |= 1 << (7 - self.bits);
        self.bits = 8;

        unsafe {
            self.flush_fast();
        }
    }

    pub fn get_output(&self) -> &[u8]
    {
        // the offset(8) is because position points 8 bytes from the bits written
        let offset = (u64::BITS / u8::BITS) as usize;

        &self.dest[self.position + offset..]
    }

    #[inline(always)]
    pub(crate) unsafe fn flush_fast(&mut self)
    {
        // align bits to the top of the  buffer
        let buf = (self.buf << ((64 - self.bits) & (u64::BITS - 1) as u8)).to_le_bytes();
        // write 8 bytes
        self.dest
            .as_mut_ptr()
            .add(self.position)
            .copy_from(buf.as_ptr(), 8);
        // but update position to point to the full number of symbols we read
        let bytes_written = self.bits & 56;
        // Decrement position
        self.position -= (bytes_written >> 3) as usize;
        self.bits &= 7;
    }

    pub fn get_position(&self) -> usize
    {
        // 8 is because position points 8 bytes from actual bits
        self.dest.len() - (self.position + (u64::BITS / u8::BITS) as usize)
    }
}

pub struct FSEStreamReader<'src>
{
    // buffer from which we are pulling in bits from
    // used in decompression.
    src: &'src [u8],
    // position in our buffer,
    position: usize,

    bits_left: u8,
    buffer: u64,
}

impl<'src> FSEStreamReader<'src>
{
    pub fn new(in_buffer: &'src [u8]) -> FSEStreamReader<'src>
    {
        FSEStreamReader {
            bits_left: 0,
            buffer: 0,
            src: in_buffer,
            position: 0,
        }
    }

    /// Read some bytes from the input buffer.
    #[inline(always)]
    pub(crate) unsafe fn refill_fast(&mut self)
    {
        /*
         * The refill always guarantees refills between 56-63
         * This version is a modification of Variant 4 found in Fabian's Giesen
         * Reading bits in far too many ways.
         * @ https://fgiesen.wordpress.com/2018/02/20/reading-bits-in-far-too-many-ways-part-2/
         *
         * Bits stored will never go above 63 and if bits are in the range 56-63 no refills occur
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
        self.bits_left |= 56;
    }

    /// Retrieve a variable amount of bits from tha bit-buffer
    ///
    /// # Arguments
    /// nbits: Number of bits to fetch from bit buffer
    /// # Returns
    /// The number of bits read
    #[inline(always)]
    fn get_bits(&mut self, nbits: u8) -> u16
    {
        debug_assert!(self.bits_left >= nbits);

        let mask = (1 << nbits) - 1;

        let bytes = (self.buffer & mask) as u16;

        self.bits_left -= nbits;

        self.buffer >>= nbits;

        return bytes;
    }
    /// Align the input bitstream to start where the
    /// the encoder left at
    pub fn align_decoder(&mut self)
    {
        unsafe {
            self.refill_fast();
        }
        // we added a marker in the encoder to tell us where the encoding ends
        // remove the marker
        let padding_bits = 1 + self.buffer.trailing_zeros();

        self.bits_left -= padding_bits as u8;

        self.buffer >>= padding_bits;
    }
    /// Decode a single symbol from a state and update next state
   
    #[inline(always)]
    pub fn decode_symbol(&mut self, state: &mut u16, dest: &mut u8, states: &[u32; TABLE_SIZE])
    {
        // It's plain cute how this is the inverse of encode_symbol, code by code.

        // format of states array
        // symbol     -> 0..8 bits
        // num_bits   -> 8..16 bits.
        // next_state -> 16..32 bits.

        // It's not faster to use get_unchecked(well on Mac OS), the ' & (TABLE_SIZE-1)'
        // kinda does the same thing
        let next_state = states[(*state & (TABLE_SIZE - 1) as u16) as usize];

        // extract symbol
        *dest = (next_state & 0xFF) as u8;

        // number of bits
        let num_bits = ((next_state >> 8) & 255) as u8;

        let low_bits = self.get_bits(num_bits);

        // Determine next state from decoder's perspective
        // (previous state from encoder's perspective)
        *state = ((next_state >> 16) as u16) + low_bits;
    }

    ///Initialize last states
    ///
    /// The encoder has to store final states and the
    /// decoder has to read final states and wound back into initial states
    pub fn init_states(&mut self) -> (u16, u16, u16, u16, u16)
    {
        unsafe {
            // refill exactly 56 bits
            self.refill_fast();

            // we added a dummy byte to ensure the flush works, so let's discard that
            self.get_bits(1);

            // each state is stored into 11 bits
            let c5 = self.get_bits(11);

            let c4 = self.get_bits(11);

            let c3 = self.get_bits(11);

            let c2 = self.get_bits(11);

            let c1 = self.get_bits(11);

            return (c1, c2, c3, c4, c5);
        }
    }
}
