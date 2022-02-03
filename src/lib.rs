#![allow(clippy::assertions_on_constants)]
#![warn(clippy::perf, clippy::pedantic, clippy::perf)]
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::inline_always
)]

mod bitstream;
mod constants;
mod huff_compress;
mod huff_decompress;
mod huff_decompress_bmi;
mod io;
mod utils;
pub use huff_compress::huff_compress_4x;
pub use huff_decompress::huff_decompress_4x;
