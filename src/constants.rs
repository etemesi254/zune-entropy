/// Upper bits to be used  for encoding
/// and decoding
pub const LIMIT: usize = 11;

pub const SMALL_CHUNK_SIZE: usize = 20;

/// Recommended table
pub const TABLE_LOG: usize = 10;
/// Maximum size of table.
/// If the table goes above this, some invaraints in encoding and decoding won't work
pub const MAX_TABLE_LOG: usize = 11;

pub const TABLE_SIZE: usize = 1 << TABLE_LOG;
