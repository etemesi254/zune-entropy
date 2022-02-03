//! Huffman Compression algorithm
//!
//!  # Format
//! 2.  block information (specific to each block)
//!     1 byte- Block information
//!         Block information is in order, I.e you can't have an uncompressed
//!         and RLE block at the same time
//!
//!         7th bit-> last block ?
//!         6th bit -> Uncompressed?
//!         5th bit -> RLE?
//!
//!      0-24 bits(3 bytes) - Total Block size.
//!
//!      3 bytes - Checksum (should probably use xxhash)
//!
//!      10 bytes jump table
//!
//!        2 bytes per jump table size.
//!
//!        Cannot go above 65536
//!
//!     10 bytes - Code lengths for symbols
//!              - There can be no symbol with code length 0 or 1
//!                code_length[i] represents codes of length i+2,
//!                so code_length[9] == codes of length 11
//!
//!     n bytes - Sum of the code lengths give the total number of
//!     symbols in ascending order.
//!
//!     Actual bits for stream 1.
//!
use std::cmp::min;
use std::convert::TryInto;
use std::io::Write;

use crate::bitstream::BitStreamWriter;
use crate::constants::LIMIT;
use crate::utils::{histogram, Symbols};

const SMALL_CHUNK_SIZE: usize = 20;

fn fast_log2(x: f32) -> f32
{
    /*
     *
     * Fast log approximation from
     * https://github.com/romeric/fastapprox
     *
     * Some pretty good stuff.
     */
    let vx = x.to_bits();
    let mx = (vx & 0x007F_FFFF) | 0x3f00_0000;
    let mx_f = f32::from_bits(mx);

    let mut y = vx as f32;
    // 1/(1<<23)
    y *= 1.192_092_9e-7;

    y - 124.225_52 - 1.498_030_3 * mx_f - 1.725_88 / (0.352_088_72 + mx_f)
}
#[allow(clippy::mut_range_bound)]
fn limited_kraft(histogram: &mut [Symbols; 256], hist_sum: u32)
{
    /*
     * Length limiting things since 2002
     *
     * This is an implementation based on steven brummes blog
     * Length-limited prefix codes (https://create.stephan-brumme.com/length-limited-prefix-codes)
     * With some tweaks for performance.
     * I avoid the heap based one due to branch mis-predictions and me not trusting my gut
     * But I've found this to be quite good(faster
     *
     *
     */
    const INITIAL_THRESHOLD: f32 = 0.4375; // (28/64)

    const STEP_THRESHOLD: f32 = 0.015_625; // (1/64)

    const ONE: usize = 1 << LIMIT;

    let mut code_lengths = [0_u8; 256];

    let inv_hist_sum = 1.0 / (hist_sum as f32);

    let mut fast_log: [f32; 256] = [0.0; 256];

    let mut spent = 0_usize;

    for i in 0..256
    {
        if histogram[i].x == 0
        {
            continue;
        }
        // Compute theoretical number of bits
        let entropy = -fast_log2((histogram[i].x as f32) * inv_hist_sum);

        fast_log[i] = entropy;

        let mut rounded = (entropy + 0.5).round() as u8;

        rounded = rounded.clamp(1, LIMIT as u8);

        code_lengths[i] = rounded;

        histogram[i].code_length = u16::from(rounded);

        spent += ONE >> rounded;
    }

    // zero the offset after
    let mut zero_offset = false;

    let mut offset = 0;
    // Kraft sum exceeds one shorten it
    let mut threshold = INITIAL_THRESHOLD;

    let mut step_threshold = STEP_THRESHOLD;

    while spent > ONE
    {
        for i in offset..256
        {
            if zero_offset
            {
                // zero out offset for the next iteration
                offset = 0;
            }
            let code_length = &mut code_lengths[i];

            if *code_length == 0
            {
                // empty symbol, no code assigned to it.
                continue;
            }
            let entropy = fast_log[i];

            let diff = entropy - f32::from(*code_length);

            let dfs = diff - threshold;

            step_threshold = step_threshold.max(dfs);

            if diff - threshold > f32::EPSILON
            {
                // diff is greater than threshold (to machine epsilon)
                // increase that code length
                *code_length += 1;

                spent -= ONE >> (*code_length);

                histogram[i].code_length = u16::from(*code_length);

                if spent <= ONE
                {
                    // detect if we are already done.
                    break;
                }
            }
            else if diff > threshold - step_threshold && offset == 0
            {
                // `i` will point to the next subsequent
                // iteration where the code length is
                // less than the threshold for that iteration.
                offset = i;
            }
        }
        if offset != 0
        {
            // make zero_offset true so that on
            // the inner loop, we can reset
            // offset once it's use.
            zero_offset = true;
        }

        threshold -= step_threshold;
    }
    // kraft sum is below 1, try to bring it to be as close
    // to one as possible
    if spent < ONE
    {
        for i in 0..256
        {
            let have = ONE >> code_lengths[i];

            if ONE - spent >= have && code_lengths[i] != 1
            {
                code_lengths[i] -= 1;

                histogram[i].code_length = u16::from(code_lengths[i]);

                spent += have;

                if ONE == spent
                {
                    break;
                }
            }
        }
    }

    histogram.sort_unstable_by(|a, b| a.code_length.cmp(&b.code_length));
}
fn generate_codes(symbols: &mut [Symbols; 256], non_zero: usize) -> [u8; LIMIT-1]
{
    /*
     * Generate actual code lengths that will be used for encoding.
     * This uses JPEG style code sizes defined in C2 of the JPEG Annex
     * https://www.w3.org/Graphics/JPEG/itu-t81.pdf.
     *
     * Do note that symbols are sorted in ascending order, and the non_zero
     * argument points to the position of the first non-zero symbol
     */

    let mut code_lengths = [0_u8; LIMIT-1];

    let mut current_size = symbols[0].code_length;

    let mut code = 0;

    for sym in symbols[non_zero..].iter_mut()
    {
        let size = sym.code_length;
        // change the current size, but only
        // when it's different from the previous code size.

        let size_diff = size - current_size;

        code <<= size_diff as usize;

        current_size = size;

        code_lengths[usize::from(sym.code_length) - 2] += 1;

        sym.x = code;

        code += 1;
    }

    let mut sym_new = [Symbols::default(); 256];
    // put symbols in the right position.
    for sym in symbols.iter()
    {
        sym_new[(sym.symbol & 255) as usize] = *sym;
    }
    *symbols = sym_new;

    code_lengths
}

/// Compress a buffer `src` and store it into
/// buffer `dest`
pub fn huff_compress_4x<W: Write>(src: &[u8], dest: &mut W)
{
    /*
     * Main code for compression.
     *
     * The code here is quite dense but it does what you expect for Huffman
     *
     * 1. Histogramming
     * 2. length limit
     * 3. Generate actual codes
     * 4. Compress
     *
     * From histogramming, we can deduce compression ratio
     * (bits for symbol * occurence) from which we may
     * decide if the block is worth compressing. If the compression ratio is small,
     * we won't compress and we emit an uncompressed block.
     *
     * If the ratio is high, we continue compressing.
     *
     * We then go and do length limiting, with 11 symbols max.
     * The choice for 11 is simple, it's a large number, and we can go for 5 codes without
     * flushing to the output buffer,(11*5 = 55 bits) hence less instructions per symbol output.
     *
     * After that we generate actual codes from the code lengths and then we move into the
     * hot loops.
     *
     * For the encoding loops, we first encode the initial symbols until
     * the src buffer is divisible by 25. When we reach this we can enter a steady
     * state of 4 unrolls each consuming 5 bytes per iteration(unrolling it further doesn't help)
     * which are encoded by the bit-stream, and after this loop ends, we simply are done.
     * we write headers, and write our compressed buffer and call it a day.
     *
     *
     */

    // Config parameter, large numbers use more memory
    // compress faster BUT hurts compression.
    // Smaller ones favour compression but hurt speed.

    // summary 323 mb enwiki file (head  -n 4000000 ./enwiki > enwiki.small)
    //
    // Machine
    // MacBook-Pro.local 21.2.0 Darwin Kernel Version 21.2.0: Sun Nov 28 20:29:10 PST 2021; root:xnu-8019.61.5~1/RELEASE_ARM64_T8101 arm64

    // |Chunk size| Speed       | Ratio    |
    // |----------|-------------|----------|
    // |1 << 18   | 1.541 Gb/s  |  0.68197 |
    // |1 << 17   | 1.486 Gb/s  |  0.67952 |
    // |1 << 16   | 1.378 Gb/s  |  0.67715 |
    // |1 << 15   | 1.248 Gb/s  |  0.67531 |
    //
    // TODO: If this is changed, there is a hard error on test files
    // investigate why
    const CHUNK_SIZE: usize = 1 << 17;

    // safety, if it goes above it can't be stored in the block.
    assert!(CHUNK_SIZE < 1 << 23);

    const START: usize = (CHUNK_SIZE + 4) / 5;

    // start is our current pointer to where we start compressing from
    let mut start = 0;
    // end is our current pointer to where we stop compressing at
    let mut end = min(CHUNK_SIZE, src.len());

    let mut src_chunk = &src[start..end];
    // Initialize stream writers
    let mut stream1 = BitStreamWriter::new();
    let mut stream2 = BitStreamWriter::new();
    let mut stream3 = BitStreamWriter::new();
    let mut stream4 = BitStreamWriter::new();
    let mut stream5 = BitStreamWriter::new();

    let mut buf = vec![0; CHUNK_SIZE + 1000];

    // Initialize destination buffers.
    // out buffer should be in buffer / 5 + 200 bytes extra for padding
    // use one large table +200 bytes(each get 200 bytes). and split into four smaller ones

    let (buf1, remainder) = buf.split_at_mut(START + 200);
    let (buf2, remainder) = remainder.split_at_mut(START + 200);
    let (buf3, remainder) = remainder.split_at_mut(START + 200);
    let (buf4, buf5) = remainder.split_at_mut(START + 200);

    let mut is_last = false;

    while !is_last
    {
        {
            let mut info_bit = [0];

            let start = (src_chunk.len() + 3) / 5;
            // 1. Count items in the buffer for histogram statistics
            let mut freq_counts = histogram(src_chunk);

            // length limit
            limited_kraft(&mut freq_counts, src_chunk.len() as u32);

            //find first non-zero element
            let non_zero = freq_counts
                .iter()
                .position(|x| x.code_length != 0)
                .unwrap_or(0);
            let last_sym = freq_counts[non_zero];

            // iterate using code lengths times symbols to determine if it's useful to compress
            let mut compressed_ratio = 0;

            for code in &freq_counts[non_zero..]
            {
                compressed_ratio += u32::from(code.code_length) * code.x;
            }

            if end == src.len()
            {
                // we reached the end,
                is_last = true;
                // indicate block is the last block
                info_bit[0] |= 1 << 7;
            }

            if compressed_ratio - 4096 > (src_chunk.len() as u32 * (u8::BITS))
            {
                // encoding didn't work, (codes were assigned a longer distribution of lengths, probably
                // a uniformly distributed data(limited-kraft doesn't like it )
                // set bit as uncompressed

                // indicate it's uncompressed
                info_bit[0] |= 1 << 6;

                // write info bit
                dest.write_all(&info_bit).unwrap();
                // total block size, in little endian
                dest.write_all(&src_chunk.len().to_le_bytes()[0..3])
                    .expect("Could not write block size");

                // write as uncompressed
                dest.write_all(src).unwrap();
            }
            else if last_sym.x == src_chunk.len() as u32
            {
                // RLE block

                // Set 5'th bit to indicate this is an RLE block
                info_bit[0] |= 1 << 5;
                // write info bit
                dest.write_all(&info_bit).unwrap();
                // total block size, in little endian
                dest.write_all(&src_chunk.len().to_le_bytes()[0..3])
                    .expect("Could not write block size");
                // write a single byte.
                dest.write_all(&[src_chunk[0]]).unwrap();
            }
            else
            {
                // generate actual codes
                // codes run from 1..=11 inclusive, e.g if a code length has been given 10 bits, it will be at
                // position 9.
                let code_lengths = generate_codes(&mut freq_counts, non_zero);

                let entries = freq_counts.map(crate::utils::Symbols::to_u32);

                // Initialize read buffers
                let (src1, remainder) = src_chunk.split_at(start);

                let (src2, remainder) = remainder.split_at(start);

                let (src3, remainder) = remainder.split_at(start);

                let (src4, src5) = remainder.split_at(start);
                // deal with symbols until all are aligned.
                let start1 = src1.len() % SMALL_CHUNK_SIZE;

                let start2 = src2.len() % SMALL_CHUNK_SIZE;

                let start3 = src3.len() % SMALL_CHUNK_SIZE;

                let start4 = src4.len() % SMALL_CHUNK_SIZE;

                let start5 = src5.len() % SMALL_CHUNK_SIZE;

                // write until all chunks are aligned to a 25 character boundary
                stream1.write_bits_slow(&src1[0..start1], &entries, buf1);

                stream2.write_bits_slow(&src2[0..start2], &entries, buf2);

                stream3.write_bits_slow(&src3[0..start3], &entries, buf3);

                stream4.write_bits_slow(&src4[0..start4], &entries, buf4);

                stream5.write_bits_slow(&src5[0..start5], &entries, buf5);

                // now chunks are aligned to 25, no need to check for remainders because they won't be there
                for ((((chunk1, chunk2), chunk3), chunk4), chunk5) in src1[start1..]
                    .chunks_exact(SMALL_CHUNK_SIZE)
                    .zip(src2[start2..].chunks_exact(SMALL_CHUNK_SIZE))
                    .zip(src3[start3..].chunks_exact(SMALL_CHUNK_SIZE))
                    .zip(src4[start4..].chunks_exact(SMALL_CHUNK_SIZE))
                    .zip(src5[start5..].chunks_exact(SMALL_CHUNK_SIZE))
                {
                    unsafe {
                        macro_rules! write_bits {
                            ($start:tt,$end:tt) => {
                                stream1.write_bits_fast(
                                    chunk1[$start..$end].try_into().unwrap(),
                                    &entries,
                                    buf1,
                                );
                                stream2.write_bits_fast(
                                    chunk2[$start..$end].try_into().unwrap(),
                                    &entries,
                                    buf2,
                                );
                                stream3.write_bits_fast(
                                    chunk3[$start..$end].try_into().unwrap(),
                                    &entries,
                                    buf3,
                                );
                                stream4.write_bits_fast(
                                    chunk4[$start..$end].try_into().unwrap(),
                                    &entries,
                                    buf4,
                                );
                                stream5.write_bits_fast(
                                    chunk5[$start..$end].try_into().unwrap(),
                                    &entries,
                                    buf5,
                                );
                            };
                        }
                        /*
                         * This version where each writes 5 bytes ia noticeably
                         * faster than where we write one bit with the streams
                         * fighting for writes , probably due to some micro-architectural issue
                         * I am yet to find.
                         *
                         *  -  I think it's cache  or too many unrolls
                         */

                        write_bits!(0, 5);

                        write_bits!(5, 10);

                        write_bits!(10, 15);

                        write_bits!(15, 20);
                    }
                }
                unsafe {
                    stream1.flush_final(buf1);
                    stream2.flush_final(buf2);
                    stream3.flush_final(buf3);
                    stream4.flush_final(buf4);
                    stream5.flush_final(buf5);
                }
                // write headers
                {
                    dest.write_all(&info_bit).unwrap();
                    // total block size, in little endian
                    dest.write(&src_chunk.len().to_le_bytes()[0..3])
                        .expect("Could not write block size");
                    // Todo, add checksum
                    dest.write_all(&[0, 0, 0])
                        .expect("Could not write checksum ");

                    // add jump tables
                    // 10 Bytes
                    dest.write_all(&stream1.get_position().to_le_bytes()[0..2])
                        .expect("Could not write jump table info");

                    dest.write_all(&stream2.get_position().to_le_bytes()[0..2])
                        .expect("Could not write jump table info");

                    dest.write_all(&stream3.get_position().to_le_bytes()[0..2])
                        .expect("Could not write jump table info");

                    dest.write_all(&stream4.get_position().to_le_bytes()[0..2])
                        .expect("Could not write jump table info");

                    dest.write_all(&stream5.get_position().to_le_bytes()[0..2])
                        .expect("Could not write jump table info");
                    // code lengths
                    dest.write_all(&code_lengths)
                        .expect("Could not write code lengths to buffer");

                    // sort again  by codes this time.
                    freq_counts.sort_unstable_by(|a, b| a.x.cmp(&b.x));
                    let sum = code_lengths.iter().map(|x| *x as usize).sum::<usize>();

                    // insert the last symbol, it was assigned a code length of 0000 hence the sort messed it up
                    freq_counts[non_zero] = last_sym;

                    // so now we have a bunch of non-zeroes, we know how many they were since we read them
                    let mut symbols = [0; 256];

                    for (sym, pos) in freq_counts[non_zero..].iter().zip(symbols.iter_mut())
                    {
                        *pos = (sym.symbol & 255) as u8;
                    }
                    // write symbols
                    dest.write_all(&symbols[0..sum])
                        .expect("Could not write symbols to destination buffer");
                }

                // write data
                // These serve two purpose
                // 1. Write data
                // 2. Check for data that is overwritten
                //      -> The write_bits_fast may write to the next stream position(since we use one large vector
                //          and subdivide it into smaller chunks) stream 1 may overwrite to stream2 position which corrupts stream
                //          2 data.
                dest.write_all(&buf1[0..stream1.get_position()])
                    .expect("Failed to write to destination buffer");

                dest.write_all(&buf2[0..stream2.get_position()])
                    .expect("Failed to write to destination buffer");

                dest.write_all(&buf3[0..stream3.get_position()])
                    .expect("Failed to write to destination buffer");

                dest.write_all(&buf4[0..stream4.get_position()])
                    .expect("Failed to write to destination buffer");

                dest.write_all(&buf5[0..stream5.get_position()])
                    .expect("Failed to write to destination buffer");

                stream1.reset();
                stream2.reset();
                stream3.reset();
                stream4.reset();
                stream5.reset();
            }
        }
        // read the next block
        start = end;
        // end increases either by chunk size or points to end of buffer.
        end = min(end + CHUNK_SIZE, src.len());

        src_chunk = &src[start..end];
    }
}

#[test]
fn huff_compress()
{
    use std::fs::{read, OpenOptions};
    use std::io::BufWriter;

    let fs = OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true)
        .open("/Users/calebe/CLionProjects/zcif/tests.zcif")
        .unwrap();
    let mut fs = BufWriter::with_capacity(1 << 24, fs);

    let fd = vec![221; 10000];

    huff_compress_4x(&fd, &mut fs);
    println!(
        "{:?} {:?}",
        fs.get_ref().metadata().unwrap().len() as f64,
        fd.len() as f64
    );
}
