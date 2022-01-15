//! Small function utilities for compression and decompression

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
    pub fn to_u32(&self) -> u32
    {
        // format
        // shift up code lengths by 4 bits
        (self.x) << 8 | u32::from(self.code_length)
    }
}

/// Calculate the occurrences of a byte in a distribution
pub fn histogram(data: &[u8]) -> (u32,[Symbols; 256])
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

        let tmp1 = u64::from_be_bytes(i[0..8].try_into().unwrap());

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
    let mut i = 0;
    let mut sum = 0;
    // add them together
    for ((((a, b), c), d), e) in val
        .iter_mut()
        .zip(start1.iter())
        .zip(start2.iter())
        .zip(start3.iter())
        .zip(start4.iter())
    {
        a.x += b + c + d + e;
        sum += a.x;
        a.symbol = i;
        i += 1;
    }

    return (sum,val);
}
