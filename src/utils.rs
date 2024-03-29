//! Small function utilities for compression and decompression

use std::io::Write;

use crate::huff_decompress::LIMIT;

pub const REVERSED_BITS: [u32; 1 << LIMIT] = reverse_bits();

#[derive(Copy, Clone, Default, Debug)]
pub struct Symbols
{
    /// Symbol in Huffman
    /// Next state in FSE
    pub z: i16,
    /// Can represent two things
    /// 1. Code lengths in huffman
    /// 2. State counts in  FSE
    pub y: u16,

    /// Can represent two things
    /// 1. Histogram count
    /// 2. Actual Huffman code
    pub x: u32,
}

impl Symbols
{
    /// Create a needed representation
    pub fn to_u32(self) -> u32
    {
        //
        // 0-8: code length
        // 8-19: code
        (REVERSED_BITS[self.x as usize] >> (16 - self.y)) << 8 | u32::from(self.y)
    }
    pub fn to_u64(self) -> u64
    {
        // 0-32 sym.z
        // 32-64 sym.x
        u64::from(self.z as u32) << 32 | u64::from(self.x)
    }
}

/// Calculate the occurrences of a byte in a distribution
pub fn histogram(data: &[u8]) -> [Symbols; 256]
{
    /*
     * The reason for splitting 4x is largely represented by
     * how lazy I am and the compiler, I found it to be the sweet spot
     * for loop unrolling in x86(arm just wont), and ILP. Most x86 CPU's
     * can add up to 4 numbers per cycle (can't find one for Arm tho).
     * and perform two loads per cycle so 4 seems good enough.
     * And the compiler removes the remainder loop all-together
     */

    // contains our count values
    let mut val = [Symbols::default(); 256];
    // allocate 4x size
    let mut counts = [0_u32; 256 * 4];
    // break into  4.
    let (start1, counts) = counts.split_at_mut(256);

    let (start2, counts) = counts.split_at_mut(256);

    let (start3, start4) = counts.split_at_mut(256);

    let chunks = data.chunks_exact(8);

    let remainder = chunks.remainder();

    for i in chunks
    {
        // count as fast as possible
        // This is the fastest platform independent histogram function I could find.
        //
        // Probably attributed to powturbo and Nathan Kurtz but it's also in
        // FSE/lib/hist.c

        let tmp1 = u64::from_le_bytes(i[0..8].try_into().unwrap());

        start1[((tmp1 >> 56) & 255) as usize] += 1;

        start2[((tmp1 >> 48) & 255) as usize] += 1;

        start3[((tmp1 >> 40) & 255) as usize] += 1;

        start4[((tmp1 >> 32) & 255) as usize] += 1;

        start1[((tmp1 >> 24) & 255) as usize] += 1;

        start2[((tmp1 >> 16) & 255) as usize] += 1;

        start3[((tmp1 >> 8) & 255) as usize] += 1;

        start4[(tmp1 & 255) as usize] += 1;
    }
    for i in remainder
    {
        start1[usize::from(*i)] += 1;
    }
    // add them together
    for (((((i, a), b), c), d), e) in val
        .iter_mut()
        .enumerate()
        .zip(start1.iter())
        .zip(start2.iter())
        .zip(start3.iter())
        .zip(start4.iter())
    {
        a.x += b + c + d + e;

        a.z = i as i16;
    }

    val
}

/// Reverse bits from MSB to LSB
const fn reverse_bits() -> [u32; 1 << LIMIT]
{
    /* this is an implementation from Eric Biggers libdeflate
     * see https://github.com/ebiggers/libdeflate
     * which is a variation of the famous stanford website bit hacks
     * @ https://graphics.stanford.edu/~seander/bithacks.html#BitReverseObvious
     * and a question asked in Stackoverflow
     * @ https://stackoverflow.com/questions/746171/efficient-algorithm-for-bit-reversal-from-msb-lsb-to-lsb-msb-in-c
     *
     * The only difference is that we do this during compilation and store the result in binary.
     */
    let mut results = [0_u32; 2048];
    let mut i = 0_u32;
    while i < 2048
    {
        let mut codeword = i as u32;

        /* Flip adjacent 1-bit fields. */
        codeword = ((codeword & 0x5555) << 1) | ((codeword & 0xAAAA) >> 1);

        /* Flip adjacent 2-bit fields. */
        codeword = ((codeword & 0x3333) << 2) | ((codeword & 0xCCCC) >> 2);

        /* Flip adjacent 4-bit fields. */
        codeword = ((codeword & 0x0F0F) << 4) | ((codeword & 0xF0F0) >> 4);

        /* Flip adjacent 8-bit fields. */
        codeword = ((codeword & 0x00FF) << 8) | ((codeword & 0xFF00) >> 8);

        results[i as usize] = codeword;

        i += 1;
    }
    results
}

pub fn write_uncompressed<W: Write>(buf: &[u8], dest: &mut W, is_last: bool)
{
    let mut info_bit = [0];
    // write info bit

    // indicate it's uncompressed
    info_bit[0] |= 1 << 7 | u8::from(is_last) << 5;

    // write info bit
    dest.write_all(&info_bit).unwrap();
    // total block size, in little endian
    dest.write_all(&buf.len().to_le_bytes()[0..3])
        .expect("Could not write block size");
    // write as uncompressed
    dest.write_all(buf).unwrap();
}

pub fn read_uncompressed(src: &[u8], block_length: u32, dest: &mut Vec<u8>)
{
    // block was uncompressed
    // assert that the above reserve actually worked
    // read to the dest buffer
    dest.extend_from_slice(&src[0..block_length as usize]);
}

pub fn read_rle(src: &[u8], block_length: u32, dest: &mut Vec<u8>)
{
    let rle = src[0];
    // read the byte

    let old_len = dest.len();
    let new_len = old_len + block_length as usize;
    // SAFETY:
    //  1. New len must be equal to or less than capacity -> confirmed by the assert
    //     statement
    //  2. Elements between old_len-new_len must be initialized -> Done straight after
    //     setting up the new length.
    assert!(dest.capacity() >= new_len, "Internal error, report to repo");

    unsafe {
        dest.set_len(new_len);
    }
    // fill with rle
    dest[old_len..new_len].fill(rle);
    // done
}

pub fn write_rle<W: Write>(src: &[u8], dest: &mut W, is_last: bool)
{
    // RLE block
    let mut info_bit = [0];
    // Set 6'th bit to indicate this is an RLE block
    info_bit[0] |= 1 << 6 | u8::from(is_last) << 5;
    // write info bit
    dest.write_all(&info_bit).unwrap();
    // total block size, in little endian
    dest.write_all(&src.len().to_le_bytes()[0..3])
        .expect("Could not write block size");
    // write a single byte.
    dest.write_all(&[src[0]]).unwrap();
}
