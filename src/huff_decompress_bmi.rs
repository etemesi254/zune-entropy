#![cfg(all(target_arch = "x86_64"))]

use crate::constants::LIMIT;

#[target_feature(enable = "lzcnt")]
#[target_feature(enable = "bmi2")]
pub unsafe fn decompress_huff_inner_bmi(
    buf: &[u8], entries: &[u16; 1 << LIMIT], offsets: &[usize; 5], block_size: usize,
    dest: &mut [u8],
)
{
    /*
     * The difference with compress_huff are small but subtle hence deserve their
     * own explanation
     *
     * In, we compress 5 bytes per go because the CPU has more work to do
     * here in decompression, we enter a fighting state where we try to decompress
     * as much as possible from different streams in flight(main loop).
     *
     * Why is because of data dependencies.
     *
     * Serially, decoding a symbol requires 5-7 cycles per byte
     * [https://fgiesen.wordpress.com/2021/08/30/entropy-coding-in-oodle-data-huffman-coding/]
     *
     * So we can't let 5 bytes to decode serially before switching to another stream
     * (e.g how encode works) since we create a large serial dependency (the CPU doesn't
     * have an infinite out of order window). So that's why we decode per byte
     * serially hoping we reach a steady state where we are limited by instruction
     * throughput and not latency
     *
     */

    use crate::bitstream::bmi::BmiBitStreamReader;

    // initialize streams with pointers to compressed data.
    let mut stream1 = BmiBitStreamReader::new(&buf[0..offsets[0]]);

    let mut stream2 = BmiBitStreamReader::new(&buf[offsets[0]..offsets[1]]);

    let mut stream3 = BmiBitStreamReader::new(&buf[offsets[1]..offsets[2]]);

    let mut stream4 = BmiBitStreamReader::new(&buf[offsets[2]..offsets[3]]);

    let mut stream5 = BmiBitStreamReader::new(&buf[offsets[3]..offsets[4]]);

    let stream_size = (block_size + 4) / 5;

    // Split the large vector into small ones
    let (dest1, rem) = dest.split_at_mut(stream_size);

    let (dest2, rem) = rem.split_at_mut(stream_size);

    let (dest3, rem) = rem.split_at_mut(stream_size);

    let (dest4, dest5) = rem.split_at_mut(stream_size);

    let mut min_pos = 0;


    // stream1.initial_refill();
    //
    // stream2.initial_refill();
    //
    // stream3.initial_refill();
    //
    // stream4.initial_refill();
    //
    // stream5.initial_refill();

    unsafe {
        for ((((a, b), c), d), e) in dest1
            .chunks_exact_mut(10)
            .zip(dest2.chunks_exact_mut(10))
            .zip(dest3.chunks_exact_mut(10))
            .zip(dest4.chunks_exact_mut(10))
            .zip(dest5.chunks_exact_mut(10))
        {
            /*
             * There is no need to check if we're nearing the end of the loop.
             * The above zip does it for us.
             */
            // decode bytes, 50 symbols per loop.

            // Unrolling by 2x improves speed from 1635 Mb/s to 1700 Mb/s(mac-os M1 arm 2020)
            macro_rules! decode_single {
                ($index:tt) => {
                    stream1.decode_single(a.get_mut($index).unwrap(), entries);

                    stream2.decode_single(b.get_mut($index).unwrap(), entries);

                    stream3.decode_single(c.get_mut($index).unwrap(), entries);

                   stream4.decode_single(d.get_mut($index).unwrap(), entries);

                     stream5.decode_single(e.get_mut($index).unwrap(), entries);
                };
            }
            macro_rules! refill {
                () => {
                    stream1.refill_fast();

                 stream2.refill_fast();

                   stream3.refill_fast();

                   stream4.refill_fast();

                    stream5.refill_fast();
                };
            }
            refill!();

            decode_single!(0);

            decode_single!(1);

            decode_single!(2);

            decode_single!(3);

            decode_single!(4);

            refill!();

            decode_single!(5);

            decode_single!(6);

            decode_single!(7);

            decode_single!(8);

            decode_single!(9);

            refill!();
            min_pos += 10;
        }
    }
    /*
     *  Do the last bytes separately
     * note, we don't check if a stream overwrote to the next start of
     * the stream
     * its practically impossible (thank you Rust) because of zips and
     * the fact that slices track their sizes and panic if you write beyond
     */
    let mut last_pos = min_pos;

    //stream1-stream4 decode equal bits, (since they are equal,
    // we can use one loop for that
    for (((a, b), c), d) in dest1[min_pos..]
        .chunks_mut(5)
        .zip(dest2[min_pos..].chunks_mut(5))
        .zip(dest3[min_pos..].chunks_mut(5))
        .zip(dest4[min_pos..].chunks_mut(5))
    {
        // we know lengths are equal, so we can do this
        let len = a.len();

        min_pos += len;
        // refill.
        unsafe {
            stream1.refill_fast();

            stream2.refill_fast();

            stream3.refill_fast();

            stream4.refill_fast();
        }
        for i in 0..len
        {
            stream1.decode_single(&mut a[i], entries);

            stream2.decode_single(&mut b[i], entries);

            stream3.decode_single(&mut c[i], entries);

            stream4.decode_single(&mut d[i], entries);
        }
    }
    // the remainder is easy
    // fill stream 5 up to what's left
    for a in dest5[last_pos..].chunks_mut(5)
    {
        unsafe {
            stream5.refill_fast();
        }
        last_pos += a.len();
        for i in a.iter_mut()
        {
            //decode each one slowly by slowly
            stream5.decode_single(i, entries);
        }
    }
    // Do the assertions
    // Assertion 1
    // All bytes were decoded.
    // min_pos contains bytes for stream 1-4, last pos contains bytes
    // decoded for stream 5, they should be equal to the end.
    assert_eq!(
        (min_pos * 4) + last_pos,
        block_size,
        "Not all bytes were decoded, internal error"
    );

    // check that no stream read more data than needed
    assert!(
        stream1.check_final(),
        "Stream 1 was possibly corrupted, position:{},length:{}",
        stream1.get_position(),
        stream1.get_src_len()
    );

    assert!(
        stream2.check_final(),
        "Stream 2 was possibly corrupted, current position:{},expected position:{}",
        stream2.get_position(),
        stream2.get_src_len()
    );

    assert!(
        stream3.check_final(),
        "Stream 3 was possibly corrupted, current position:{},expected position:{}",
        stream3.get_position(),
        stream3.get_src_len()
    );

    assert!(
        stream4.check_final(),
        "Stream 4 was possibly corrupted, current position:{},expected position:{}",
        stream4.get_position(),
        stream4.get_src_len()
    );

    assert!(
        stream5.check_final(),
        "Stream 5 was possibly corrupted, current position:{},expected position:{}",
        stream5.get_position(),
        stream5.get_src_len()
    );
    //everything is good, we won Mr Stark.
}
