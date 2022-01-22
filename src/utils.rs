//! Small function utilities for compression and decompression

use crate::huff_decompress::LIMIT;

pub const REVERSED_BITS:[u32;1<<LIMIT] = reverse_bits();

#[derive(Copy, Clone, Default, Debug)]
pub struct Symbols
{
    pub symbol: u16,
    /// At most uses 11 bits
    pub code_length: u16,

    /// Can represent two things
    /// 1. Histogram count
    /// 2. Actual code length.
    pub x: u32,
}

impl Symbols
{
    /// Create a needed representation
    pub fn to_u32(self) -> u32
    {
        //
        //
        // 0-8: code length
        // 8-19: code
        (REVERSED_BITS[self.x as usize] >> (16 - self.code_length)) << 8 | u32::from(self.code_length)
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

        //*hasher = hash(*hasher);

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
    for (((((i,a), b), c), d), e) in val
        .iter_mut()
        .enumerate()
        .zip(start1.iter())
        .zip(start2.iter())
        .zip(start3.iter())
        .zip(start4.iter())
    {
        a.x += b + c + d + e;
        a.symbol = i as u16;
    }

    val
}

/// Reverse bits from MSB to LSB
const fn reverse_bits() -> [u32; 1 << LIMIT]
{
    let mut results = [0_u32; 2048];
    let mut i = 0_u32;
    while i<2048{
    
        // https://stackoverflow.com/questions/746171/efficient-algorithm-for-bit-reversal-from-msb-lsb-to-lsb-msb-in-c
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
        i+=1;
    }
    return results;
}
