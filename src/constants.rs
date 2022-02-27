/// Upper bits to be used  for encoding
/// and decoding
pub const LIMIT: usize = 11;

pub const SMALL_CHUNK_SIZE: usize = 20;

/// Recommended table
pub const TABLE_LOG: usize = 11;
/// Maximum size of table.
/// If the table goes above this, some invaraints in encoding and decoding won't work
pub const MAX_TABLE_LOG: usize = 11;

pub const TABLE_SIZE: usize = 1 << TABLE_LOG;

pub fn state_generator(num_states: usize) -> usize
{
    // play with this

    /*  Apparently this state distribution works well for  11 state bits
     *  which is what we use, mainly
     *  But for smaller state counts, e.g if we decide to use 9 bits, switch
     *  to the more well distributed state of
     *  of (num_states >> 1) | (num_states>>3) | 3 (fse/xpack)
     *  or probably use a different distribution hueristic.
     */
    return (num_states >> 2) | (num_states >> 3) | 5;
}
