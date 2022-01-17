//! Huffman Compression algorithm

use std::convert::TryInto;
use std::io::{Read, Write, Seek, SeekFrom};
use crate::bitstream::BitStreamWriter;
use crate::huff_decompress::LIMIT;
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
    let mx = (vx & 0x007FFFFF) | 0x3f000000;
    let mx_f = f32::from_bits(mx);

    let mut y = vx as f32;
    // 1/(1<<23)
    y *= 1.1920928955078125e-7;

    return y - 124.22551499 - 1.498030302 * mx_f - 1.72587999 / (0.3520887068 + mx_f);
}
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
    

        code_lengths[usize::from(sym.code_length) - 1] += 1;

        sym.x = code;

        code += 1;

    }

    
    let mut sym_new = [Symbols::default();256];
    // put symbols in the right position.
    for sym in symbols.iter(){
        sym_new[(sym.symbol & 255) as usize] = *sym; 
    }
    *symbols = sym_new;

    return code_lengths;
}

/// Compress a buffer `src` and store it into
/// buffer `dest`
pub fn huff_compress_4x<R:Read+Seek, W: Write>(src: &mut R, dest: &mut W)
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

    // This depends on the CPU cache size.
    // For mac os, we have it being 132 kb , since L1 is 192 kb 
    // increasing that will cause L1 cache thrashing. 

    // TODO:Pragammatically get L1 cache size?
    const CHUNK_SIZE: usize = 1 << 17;

    let mut src_buf = vec![0; CHUNK_SIZE];
    // Some weird bug would arise where the reader would
    // point to the end of the file when reading from in memory buffers.
    // I'm too lazy to try and understand it , so here is my workaround
    // Just seek it to be zero.
    src.seek(SeekFrom::Start(0)).unwrap();

    // size is how many bytes were actually read.
    let mut size = src.read(&mut src_buf).unwrap();
    
    // Initialize stream writers
    let mut stream1 = BitStreamWriter::new();
    let mut stream2 = BitStreamWriter::new();
    let mut stream3 = BitStreamWriter::new();
    let mut stream4 = BitStreamWriter::new();
    let mut stream5 = BitStreamWriter::new();

    let mut buf = vec![0; CHUNK_SIZE + 1000];

    const START: usize = (CHUNK_SIZE + 4) / 5;

    // Initialize destination buffers.
    // out buffer should be in buffer / 4 + 200 bytes extra for padding
    // use one large table +200 bytes(each get 200 bytes). and split into four smaller ones

    let (buf1, remainder) = buf.split_at_mut(START + 200);
    let (buf2, remainder) = remainder.split_at_mut(START + 200);
    let (buf3, remainder) = remainder.split_at_mut(START + 200);
    let (buf4,buf5) = remainder.split_at_mut(START + 200);

    loop
    {
        // chunk depending on how much data we read
        for src_chunk in src_buf[0..size].chunks(size)
        {
            let start = (src_chunk.len() + 3) / 5;
            // 1. Count items in the buffer for histogram statistics
            let (hist_sum, mut freq_counts) = histogram(src_chunk);

            // length limit
            limited_kraft(&mut freq_counts, hist_sum);

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
            if end - 4096 > (src_chunk.len() as u32 * (u8::BITS))
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

                let freq_counts_stream =
                    freq_counts.iter().map(|x| x.to_u32()).collect::<Vec<u32>>();

                let freq_count_stream: [u32; 256] = freq_counts_stream.try_into().unwrap();

                // Initialize read buffers
                let (src1, remainder) = src_chunk.split_at(start);
                let (src2, remainder) = remainder.split_at(start);
                let (src3, remainder) = remainder.split_at(start);
                let (src4,src5) = remainder.split_at(start);
                // deal with symbols until all are aligned.
                let start1 = src1.len() % SMALL_CHUNK_SIZE;
                let start2 = src2.len() % SMALL_CHUNK_SIZE;
                let start3 = src3.len() % SMALL_CHUNK_SIZE;
                let start4 = src4.len() % SMALL_CHUNK_SIZE;
                let start5 = src5.len() % SMALL_CHUNK_SIZE;

                // write until all chunks are aligned to a 25 character boundary
                stream1.write_bits_slow(&src1[0..start1], &freq_count_stream, buf1);

                stream2.write_bits_slow(&src2[0..start2], &freq_count_stream, buf2);

                stream3.write_bits_slow(&src3[0..start3], &freq_count_stream, buf3);

                stream4.write_bits_slow(&src4[0..start4], &freq_count_stream, buf4);

                stream5.write_bits_slow(&src5[0..start5], &freq_count_stream, buf5);

                // now chunks are aligned to 25, no need to check for remainders because they won't be there
                for ((((chunk1, chunk2), chunk3), chunk4),chunk5) in src1[start1..]
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
                                    &freq_count_stream,
                                    buf1,
                                );
                                stream2.write_bits_fast(
                                    chunk2[$start..$end].try_into().unwrap(),
                                    &freq_count_stream,
                                    buf2,
                                );
                                stream3.write_bits_fast(
                                    chunk3[$start..$end].try_into().unwrap(),
                                    &freq_count_stream,
                                    buf3,
                                );
                                stream4.write_bits_fast(
                                    chunk4[$start..$end].try_into().unwrap(),
                                    &freq_count_stream,
                                    buf4,
                                );
                                stream5.write_bits_fast(
                                    chunk5[$start..$end].try_into().unwrap(),
                                    &freq_count_stream,
                                    buf5,
                                );
                            };
                        }

                        write_bits!(0, 5);
                        write_bits!(5, 10);
                        write_bits!(10, 15);
                        write_bits!(15, 20);
                        
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


                dest.write(&buf5[0..stream5.get_position()])
                    .expect("Failed to write to destination buffer");
                
                stream1.reset();
                stream2.reset();
                stream3.reset();
                stream4.reset();
                stream5.reset();
            }
        }

        size = src.read(&mut src_buf).unwrap();

        if size == 0
        {
            break;
        }
    }
}

#[test]
fn huff_compress()
{
    use std::fs::OpenOptions;
    use std::io::{BufReader,BufWriter};
    
    let fs = OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true)
        .open("/Users/calebe/CLionProjects/zcif/tests")
        .unwrap();
    let mut fs = BufWriter::with_capacity(1 << 24, fs);

    let fd = OpenOptions::new()
        .read(true)
        .open("/Users/calebe/git/FiniteStateEntropy/programs/enwiki")
        .unwrap();
    let mut fd = BufReader::new(fd);

    huff_compress_4x(&mut fd, &mut fs);
    println!(
        "{:?}",
        fs.get_ref().metadata().unwrap().len() as f64
            / fd.get_ref().metadata().unwrap().len() as f64
    );
}
