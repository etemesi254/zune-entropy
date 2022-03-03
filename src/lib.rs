//! A bunch of really fast entropy coders
#![feature(stdsimd)]
#![allow(clippy::assertions_on_constants)]
#![warn(clippy::perf, clippy::pedantic, clippy::perf)]
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::inline_always
)]
extern crate core;

mod constants;

mod fse_bitstream;
mod huff_bitstream;

mod fse_compress;
mod huff_compress;

mod fse_decompress;
mod huff_decompress;

mod utils;

pub use fse_compress::fse_compress;
pub use fse_decompress::fse_decompress;
pub use huff_compress::huff_compress;
pub use huff_decompress::huff_decompress;
