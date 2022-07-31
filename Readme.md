# Fast Entropy coders
This repo contains Rust implementations of Entropy coders , namely a Huffman Entropy coder and 
a Table Based Asymmetric numerical system.

Both of these are compression algorithms that work by using less bits to represent frequently occurring 
characters and more bits to represent less frequently occurring characters.

The library is very much in beta  and some files may not work(like `siliesa/mozilla`)

### Huffman 

Do note that the Huffman encoder and decoder implemented here is not the best compression wise,
you can achieve better compression ratio with Huffman encoders.

This is just meant to be a simple(and a fast) illustration of Huffman encoding and decoding.

### TANS
Table Based Asymmetric Numerical systems stem from Jarek Duda's ANS and Yann Collet's work on Finite State Entropy.

It features better compression ratio(Closer to the Shannon Limit) as compared to Huffman, and it really shines with 
skewed distributions.

It does have an asymptomatically slower speeds than Huffman as it does more work per symbol than a Huffman Encoder/Decoder, 
this is visible in really aggressive out-of-order machines(see benchmarks) but is still a good replacement for it.



## Benchmarks.

For running benchmarks of this library there is another library [here](https://github.com/etemesi254/zcif_bin/tree/main)

I'm going to compare it to Yann Collet's FSE library since they both do the same things extremely well
and to also point out a few gotchas.

#### Machine Specs
- Model name:          AMD Ryzen 5 4500U with Radeon Graphics
- CPU family:          23
- Model:               96
- Thread(s) per core:  1
- Core(s) per socket:  6


- L1d:                   192 KiB (6 instances)
- L1i:                   192 KiB (6 instances)
- L2:                    3 MiB (6 instances)
- L3:                    8 MiB (2 instances)

#### Speeds 

| File        | Library   | Encoder      | Decoder      | Compression ratio |
|-------------|-----------|--------------|--------------|-------------------|
| **enwiki8** |
|             | Huffman   | 940.97 MB/s  | 1310.28 MB/s | 65.00%            |
|             | TANS      | 602.56 MB/s  | 1022.93 MB/s | 63.90%            |
||
|             | FSE/Huff0 | 924.6 MB/s   | 1049.5 MB/s  | 63.26%            |
|             | FSE/TANS  | 457.7 MB/s   | 572.8 MB/s   | 63.12%            |


Note that FSE/* beats this library by compression ratio on both Huffman and TANS, it employs better and more 
sophisticated algorithms while I'm using the simplest ones I could come up with.

The TANS implementation is faster since this library employs more interleaving which pays off well

# Acknowledgements

This couldn't have been made possible without the great works of 

- Jarek Duda - The inventor of ANS based algorithms, and collector, most of his
work and others can be found [here](https://encode.su/threads/2078-List-of-Asymmetric-Numeral-Systems-implementations)

- Eric Biggers - [xpack](https://github.com/ebiggers/xpack) (lz77+tANS).
- Yann Collet - [FiniteStateEntropy](https://github.com/Cyan4973/FiniteStateEntropy),
- Yann Collet - [FSE algorithm blog posts](http://fastcompression.blogspot.com/2013/12/finite-state-entropy-new-breed-of.html)
- Charles Bloom - [Understanding TANS](http://cbloomrants.blogspot.com/2014/01/1-30-14-understanding-ans-1.html) blog series,
- Fabian Geisen - [Reading bits in far too many ways](https://fgiesen.wordpress.com/2018/02/19/reading-bits-in-far-too-many-ways-part-1/)
- Fabian Geisen -  [Interleaving rANS](https://arxiv.org/abs/1402.3392)
