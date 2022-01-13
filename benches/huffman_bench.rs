use std::fs::OpenOptions;
use std::io::{BufWriter, Seek, SeekFrom};
use std::time::Duration;

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use zcif::huff_compress_4x;
fn criterion_benchmark(c: &mut Criterion)
{
    // 1 mb
    let mut random = vec![0_u8; 1 << 24];
    let fs = OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true)
        .open("/Users/calebe/CLionProjects/zcif/tests")
        .unwrap();
    let mut fs = BufWriter::with_capacity(1 << 24, fs);
    use rand::{thread_rng, Rng};
    thread_rng().fill(&mut random[..]);
    random.sort_unstable();

    c.bench_function("Huff Compress 4X", |b| {
        b.iter(|| black_box(huff_compress_4x(&random, &mut fs)));
        // set length to be zero on all iterations
        fs.get_mut().set_len(0);
        fs.seek(SeekFrom::Start(0));
    });
}

criterion_group!(name=benches;
    config={
      let c = Criterion::default();
        c.measurement_time(Duration::from_secs(20))
      };
    targets=criterion_benchmark);
criterion_main!(benches);
