use std::time::Duration;

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use zcif::huff_compress_4x;
fn criterion_benchmark(c: &mut Criterion)
{
    let mut writer = Vec::with_capacity(85 * (1 << 20));

    let mut reader = std::io::Cursor::new(
        std::fs::read("/Users/calebe/git/FiniteStateEntropy/programs/enwiki.small").unwrap(),
    );

    c.bench_function("Huff Compress 4X", |b| {
        b.iter(|| {
            huff_compress_4x(&mut reader, &mut writer);
            black_box(());
            writer.clear();
            reader.set_position(0);
        });
        // set length to be zero on all iterations
    });
}

criterion_group!(name=benches;
    config={
      let c = Criterion::default();
        c.measurement_time(Duration::from_secs(20))
      };
    targets=criterion_benchmark);
criterion_main!(benches);
