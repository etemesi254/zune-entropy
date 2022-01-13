#![feature(portable_simd)]
mod bitstream;
mod huff_compress;
mod huff_decompress;
mod io;
mod utils;
pub use huff_compress::huff_compress_4x;
pub use huff_decompress::{decompress_huff_double_symbols, decompress_huff_single_symbols};
#[cfg(test)]
mod tests
{
    #[test]
    fn it_works()
    {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}
