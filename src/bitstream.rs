//! BitStreamReader API
//!
//! This module provides an interface to read and write bits (and bytes)


#[derive(Eq, PartialEq, Clone, Copy)]
pub enum Flags
{
    UnderFlow,
    Normal,
}
pub struct BitStreamReader<'src>
{
    // buffer from which we are pulling in bits from
    // used in decompression.
    src: &'src [u8],
    // position in our buffer,
    position: usize,
    // End of buffer flag
    flag: Flags,
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
        let stream = BitStreamReader {
            src: in_buffer,
            // buffer is read from end to start.
            position: in_buffer.len(),
            flag: Flags::Normal,
        };

        // return
        stream
    }
    /// Refill the bitstream ensuring the buffer has bits between
    /// 56 and 63.
    ///
    #[inline(always)]
    pub unsafe fn refill_fast(
        &mut self, // current buffer of our bits
        buffer: &mut u64,
        // how many bits are left
        bits_left: &mut u8,
    ) -> Flags
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

        // read 8 bytes to a temporary buffer
        if self.position <= 8
        {
            // switch to the careful reload
            self.flag = Flags::UnderFlow;
            return self.flag;
        }

        let mut buf = [0; 8];

        // position points to the end initially, so subtracting means we
        // are reading the buffer from end to start.
        std::ptr::copy_nonoverlapping(
            self.src.as_ptr().sub(self.position - 8),
            buf.as_mut_ptr(),
            8,
        );

        // create a u64 from an array of u8's
        let new_buffer = u64::from_le_bytes(buf);
        // num indicates how many bytes we actually consumed.
        // since bytes occupy 8 bits, we have to consume bits in multiples of 8.
        let num = (i32::from(63 - *bits_left) & (-8)) as u8;
        // Subtract what we read from pointer
        self.position -= (num >> 3) as usize;
        // shift number of bits
        *buffer |= new_buffer >> num;
        // update bits left
        // bits left are now between 56-63
        *bits_left |= 56;

        return self.flag;
    }
}

pub struct BitStreamWriter
{
    // kept below
    empty_bits: u8,
    // should I use usize?
    //
    buf: u64,
    position: usize,
}

impl BitStreamWriter
{
    pub fn new() -> BitStreamWriter
    {
        BitStreamWriter {
            buf: 0,
            empty_bits: 0,
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
            self.buf |= u64::from(entry >> 8) << self.empty_bits;
            
            self.empty_bits += (entry & 0xFF) as u8;

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
        * to 5 symbols per encode (55).
        *
        * The symbols are read in little endian order, symbol[0] is stored
        * at bits 24..32(for variable large_entry). symbol[4] is stored at bit position(0..8)
        *
        * Loads are optimized to one  large variable , (single load is slower).
        *
        *
        * The fifth symbol is the weird one since it doesn't fit into a u32 hence 
        * it is handles in the last bit
        */

        let large_entry = u32::from_le_bytes(symbols[0..4].try_into().unwrap());

        macro_rules! encode_single {
            ($pos:tt) => {
                let entry = entry[((large_entry >> ($pos)) & 255) as usize];

                // add to the top bits
                self.buf |= (u64::from(entry >> 8) << self.empty_bits);

                 self.empty_bits += (entry & 0xFF) as u8;

            };
        }

        encode_single!(0);

        encode_single!(8);

        encode_single!(16);

        encode_single!(24);
        // the black sheep (can't fit into u32)
        let entry = entry[usize::from(symbols[4])];

        // add to the top bits
        self.buf |= u64::from(entry >> 8) << self.empty_bits;

        self.empty_bits += (entry & 0xFF) as u8;

        // flush to output buffer
        self.flush_fast(out_buf);
    }

    pub unsafe fn flush_fast(&mut self, out_buf: &mut [u8])
    {
        // bits are in the top buffer arranged in the top buffer following each other
        // the first bit is in the MSB.
        // take the big endian representation
        let buf = self.buf.to_le_bytes();
        // write 8 bytes
        out_buf
            .as_mut_ptr()
            .add(self.position)
            .copy_from(buf.as_ptr(), 8);
        // but update position to point to the full number of symbols we read
        let bytes_written = self.empty_bits & 56;
        // remove those bits we read.
        self.buf >>= bytes_written;
        // increment position
        self.position += (bytes_written >> 3) as usize;

        self.empty_bits &= 7;
    }
    pub fn reset(&mut self)
    {
        self.empty_bits = 0;
        self.buf = 0;
        self.position = 0;
    }
    pub fn get_position(&self) -> usize
    {
        self.position
    }
}
