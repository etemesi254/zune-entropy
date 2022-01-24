#![allow(clippy::assertions_on_constants)]
#![warn(clippy::perf)]
mod bitstream;
mod constants;
mod huff_compress;
mod huff_decompress;
mod io;
mod utils;
pub use huff_compress::huff_compress_4x;
pub use huff_decompress::huff_decompress_4x;
