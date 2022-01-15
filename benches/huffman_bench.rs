use std::fs::{File, OpenOptions};
use std::io::{BufReader, BufWriter, Seek, SeekFrom};
use std::time::Duration;

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use zcif::huff_compress_4x;
fn criterion_benchmark(c: &mut Criterion)
{
    let fs = OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true)
        .open("/Users/calebe/CLionProjects/zcif/tests")
        .unwrap();
    let mut fs = BufWriter::with_capacity(1 << 24, fs);

    let mut fd = OpenOptions::new().read(true).open("/Users/calebe/git/FiniteStateEntropy/programs/enwiki").unwrap();
    let mut fd = BufReader::new(fd);


    c.bench_function("Huff Compress 4X", |b| {
        b.iter(|| black_box(huff_compress_4x(&mut fd, &mut fs)));
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
