//! Huffman Compression algorithm

use std::convert::TryInto;
use std::io::Write;

use crate::bitstream::BitStreamWriter;
use crate::huff_decompress::LIMIT;
use crate::utils::{histogram, Symbols};

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
    let mx = (vx & 0x007FFFFF) | 0x3f000000;
    let mx_f = f32::from_bits(mx);

    let mut y = vx as f32;
    // 1/(1<<23)
    y *= 1.1920928955078125e-7;

    return y - 124.22551499 - 1.498030302 * mx_f - 1.72587999 / (0.3520887068 + mx_f);
}
fn limited_kraft(histogram: &mut [Symbols; 256],hist_sum:u32)
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
    const STEP_THRESHOLD: f32 = 0.015625; // (1/64)

    let mut code_lengths = [0_u8; 256];

    let inv_hist_sum = 1.0 / (hist_sum as f32);

    const ONE: usize = 1 << LIMIT;

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

        histogram[i].code_length = rounded as u16;

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

            let diff = entropy - (*code_length as f32);

            let dfs = diff - threshold;

            step_threshold = step_threshold.max(dfs);

            if diff - threshold > f32::EPSILON
            {
                // diff is greater than threshold (to machine epsilon)
                // increase that code length
                *code_length += 1;

                spent -= ONE >> (*code_length);
                histogram[i].code_length = (*code_length) as u16;

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
fn generate_codes(symbols: &mut [Symbols; 256], non_zero: usize) -> [u8; LIMIT]
{
    /*
     * Generate actual code lengths that will be used for encoding.
     * This uses JPEG style code sizes defined in C2 of the JPEG Annex
     * https://www.w3.org/Graphics/JPEG/itu-t81.pdf.
     *
     * Do note that symbols are sorted in ascending order, and the non_zero
     * argument points to the position of the first non-zero symbol
     */

    let mut code_lengths = [0_u8; LIMIT];

    let mut current_size = symbols[non_zero].code_length;

    let mut code = 1 << current_size;


    for sym in symbols[non_zero..].iter_mut()
    {
        // change the current size, but only
        // when it's different from the previous code size.
        let ls = 1 & u8::from(current_size != sym.code_length);

        current_size = sym.code_length;

        code_lengths[usize::from(sym.code_length) - 1] += 1;

        code <<= ls;

        sym.x = code;

        code += 1;
    }

    symbols.sort_unstable_by(|a, b| a.symbol.cmp(&b.symbol));
    return code_lengths;
}

/// Compress a buffer `src` and store it into
/// buffer `dest`
pub fn huff_compress_4x<W: Write>(src: &[u8], dest: &mut W)
{
    const CHUNK_SIZE: usize = 1 << 16;

    // Initialize stream writers
    let mut stream1 = BitStreamWriter::new();
    let mut stream2 = BitStreamWriter::new();
    let mut stream3 = BitStreamWriter::new();
    let mut stream4 = BitStreamWriter::new();

    let mut buf = vec![0; CHUNK_SIZE + 800];

    const START: usize = (CHUNK_SIZE + 3) / 4;

    // Initialize destination buffers.
    // out buffer should be in buffer / 4 + 200 bytes extra for padding
    // use one large table +200 bytes(each get 200 bytes). and split into four smaller ones

    let (buf1, remainder) = buf.split_at_mut(START + 200);
    let (buf2, remainder) = remainder.split_at_mut(START + 200);
    let (buf3, buf4) = remainder.split_at_mut(START + 200);
    // chunk into regular sizes
    for src in src.chunks(CHUNK_SIZE)
    {
        let start = (src.len() + 3) / 4;
        // 1. Count items in the buffer for histogram statistics
        let (hist_sum,mut freq_counts) = histogram(src);

        // length limit
        limited_kraft(&mut freq_counts,hist_sum);

        //find first non-zero element
        let non_zero = freq_counts
            .iter()
            .position(|x| x.code_length != 0)
            .unwrap_or(0);

        // iterate using code lengths times symbols to determine if it's useful to compress
        let mut end = 0;
        for code in &freq_counts[non_zero..]
        {
            end += u32::from(code.code_length) * u32::from(code.x);
        }
        if end - 4096 > (src.len() as u32 * (u8::BITS))
        {
            // encoding didn't work, (codes were assigned a longer distribution of lengths, probably
            // a uniformly distributed data(limited-kraft doesn't like it )
            // TODO: Print some stats
            // emit as it is uncompressed
        }
        else
        {
            // generate actual codes
            // codes run from 1..=11 inclusive, e.g if a code length has been given 10 bits, it will be at
            // position 9.
            let code_lengths = generate_codes(&mut freq_counts, non_zero);

            let freq_counts_stream = freq_counts.iter().map(|x| x.to_u32()).collect::<Vec<u32>>();

            let freq_count_stream: [u32; 256] = freq_counts_stream.try_into().unwrap();

            // Initialize read buffers
            let (src1, remainder) = src.split_at(start);
            let (src2, remainder) = remainder.split_at(start);
            let (src3, src4) = remainder.split_at(start);
            // deal with symbols until all are aligned

            let start1 = src1.len() % 20;
            let start2 = src2.len() % 20;
            let start3 = src3.len() % 20;
            let start4 = src4.len() % 20;
            // write until all chunks are aligned to a 25 character boundary
            stream1.write_bits_slow(&src1[0..start1], &freq_count_stream, buf1);

            stream2.write_bits_slow(&src2[0..start2], &freq_count_stream, buf2);

            stream3.write_bits_slow(&src3[0..start3], &freq_count_stream, buf3);

            stream4.write_bits_slow(&src4[0..start4], &freq_count_stream, buf4);

            // now chunks are aligned to 25, no need to check for remainders because they won't be there
            for (((chunk1, chunk2), chunk3), chunk4) in src1[start1..]
                .chunks_exact(20)
                .zip(src2[start2..].chunks_exact(20))
                .zip(src3[start3..].chunks_exact(20))
                .zip(src4[start4..].chunks_exact(20))
            {
                unsafe {
                    stream1.write_bits_fast(
                        chunk1[0..5].try_into().unwrap(),
                        &freq_count_stream,
                        buf1,
                    );
                    stream2.write_bits_fast(
                        chunk2[0..5].try_into().unwrap(),
                        &freq_count_stream,
                        buf2,
                    );
                    stream3.write_bits_fast(
                        chunk3[0..5].try_into().unwrap(),
                        &freq_count_stream,
                        buf3,
                    );
                    stream4.write_bits_fast(
                        chunk4[0..5].try_into().unwrap(),
                        &freq_count_stream,
                        buf4,
                    );

                    stream1.write_bits_fast(
                        chunk1[5..10].try_into().unwrap(),
                        &freq_count_stream,
                        buf1,
                    );
                    stream2.write_bits_fast(
                        chunk2[5..10].try_into().unwrap(),
                        &freq_count_stream,
                        buf2,
                    );
                    stream3.write_bits_fast(
                        chunk3[5..10].try_into().unwrap(),
                        &freq_count_stream,
                        buf3,
                    );
                    stream4.write_bits_fast(
                        chunk4[5..10].try_into().unwrap(),
                        &freq_count_stream,
                        buf4,
                    );

                    stream1.write_bits_fast(
                        chunk1[10..15].try_into().unwrap(),
                        &freq_count_stream,
                        buf1,
                    );
                    stream2.write_bits_fast(
                        chunk2[10..15].try_into().unwrap(),
                        &freq_count_stream,
                        buf2,
                    );
                    stream3.write_bits_fast(
                        chunk3[10..15].try_into().unwrap(),
                        &freq_count_stream,
                        buf3,
                    );
                    stream4.write_bits_fast(
                        chunk4[10..15].try_into().unwrap(),
                        &freq_count_stream,
                        buf4,
                    );

                    stream1.write_bits_fast(
                        chunk1[15..20].try_into().unwrap(),
                        &freq_count_stream,
                        buf1,
                    );
                    stream2.write_bits_fast(
                        chunk2[15..20].try_into().unwrap(),
                        &freq_count_stream,
                        buf2,
                    );
                    stream3.write_bits_fast(
                        chunk3[15..20].try_into().unwrap(),
                        &freq_count_stream,
                        buf3,
                    );
                    stream4.write_bits_fast(
                        chunk4[15..20].try_into().unwrap(),
                        &freq_count_stream,
                        buf4,
                    );
                }
            }
            // write headers
            {
                // code lengths
                dest.write(&code_lengths)
                    .expect("Could not write code lengths to buffer");
                // sort again  by codes this time.
                freq_counts.sort_by(|a, b| a.x.cmp(&b.x));
                // so now we have a bunch of non-zeroes, we know how many they were since we read them
                let mut symbols = [0; 256];
                for (sym, pos) in freq_counts[non_zero..].iter().zip(symbols.iter_mut())
                {
                    *pos = (sym.symbol & 255) as u8;
                }
                dest.write(&symbols[0..(256 - non_zero)])
                    .expect("Could not write symbols to destination buffer");
            }
            // write data
            dest.write(&buf1[0..stream1.get_position()])
                .expect("Failed to write to destination buffer");

            dest.write(&buf2[0..stream2.get_position()])
                .expect("Failed to write to destination buffer");

            dest.write(&buf3[0..stream3.get_position()])
                .expect("Failed to write to destination buffer");

            dest.write(&buf4[0..stream4.get_position()])
                .expect("Failed to write to destination buffer");

            stream1.reset();
            stream2.reset();
            stream3.reset();
            stream4.reset();
        }
    }
}

#[test]
fn huff_compress()
{
    let mut random = vec![0_u8; 1 << 20];
    use rand::{thread_rng, Rng};
    thread_rng().fill(&mut random[..]);
    random.sort_unstable();
    let mut dest = vec![];
    //random.sort_unstable();
    huff_compress_4x(&random[..], &mut dest);
    println!("{:?}",random.len() as f32/dest.len() as f32);

    println!("{:?}, {:?}",&dest[0..11],&dest[11..11+(&dest[0..11]).iter().map(|a| usize::from(*a)).sum::<usize>()])
}
