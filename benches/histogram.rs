//! Small function utilities for compression and decompression

use criterion::{criterion_group, criterion_main, Criterion};
use rand::Rng;

fn histogram_4x(data: &[u8])
{
    // allocate 4x size
    let mut counts = [0_u32; 256 * 4];
    // break into  4.
    let (start1, remainder) = counts.split_at_mut(256);
    let (start2, remainder) = remainder.split_at_mut(256);
    let (start3, start4) = remainder.split_at_mut(256);
    let chunks = data.chunks_exact(8);
    let remainder = chunks.remainder();
    for i in chunks
    {
        let tmp1 = u32::from_be_bytes(i[0..4].try_into().unwrap());
        let tmp2 = u32::from_be_bytes(i[4..8].try_into().unwrap());

        start1[((tmp1 >> 24) & 255) as usize] += 1;
        start2[((tmp1 >> 16) & 255) as usize] += 1;
        start3[((tmp1 >> 8) & 255) as usize] += 1;
        start4[(tmp1 & 255) as usize] += 1;

        start1[((tmp2 >> 24) & 255) as usize] += 1;
        start2[((tmp2 >> 16) & 255) as usize] += 1;
        start3[((tmp2 >> 8) & 255) as usize] += 1;
        start4[((tmp2 & 255) as usize)] += 1;
    }
    for i in remainder
    {
        start1[usize::from(*i)] += 1;
    }
    // add them together
    for (((a, b), c), d) in start1
        .iter_mut()
        .zip(start2.iter())
        .zip(start3.iter())
        .zip(start4.iter())
    {
        *a += b + c + d;
    }
}

fn histogram_4x_u64(data: &[u8])
{
    // allocate 4x size
    let mut counts = [0_u32; 256 * 4];
    // break into  4.
    let (start1, remainder) = counts.split_at_mut(256);
    let (start2, remainder) = remainder.split_at_mut(256);
    let (start3, start4) = remainder.split_at_mut(256);
    let chunks = data.chunks_exact(8);
    let remainder = chunks.remainder();
    for i in chunks
    {
        let tmp1 = u64::from_be_bytes(i[0..8].try_into().unwrap());

        start1[((tmp1 >> 56) & 255) as usize] += 1;
        start2[((tmp1 >> 48) & 255) as usize] += 1;
        start3[((tmp1 >> 40) & 255) as usize] += 1;
        start4[((tmp1 >> 32) & 255) as usize] += 1;

        start1[((tmp1 >> 24) & 255) as usize] += 1;
        start2[((tmp1 >> 16) & 255) as usize] += 1;
        start3[((tmp1 >> 8) & 255) as usize] += 1;
        start4[((tmp1 & 255) as usize)] += 1;
    }
    for i in remainder
    {
        start1[usize::from(*i)] += 1;
    }
    // add them together
    for (((a, b), c), d) in start1
        .iter_mut()
        .zip(start2.iter())
        .zip(start3.iter())
        .zip(start4.iter())
    {
        *a += b + c + d;
    }
}
fn criterion_benchmark(c: &mut Criterion)
{
    let mut rng = rand::thread_rng();
    let mut t = [0_u8; 1_000_000];

    rng.fill(&mut t[..]);

    c.bench_function("histogram 4x", |b| b.iter(|| histogram_4x(&t)));

    c.bench_function("histogram 4x u64", |b| b.iter(|| histogram_4x_u64(&t)));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
