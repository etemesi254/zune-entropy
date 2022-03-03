# A bunch of really fast entropy coders

This repo contains two types of static entropy coders, Huffman and FSE/tANS

# Summary 

File `enwiki.small`

Created from the command

```shell
head -n 400000 enwiki >> enwiki.small
```

The file size is **338929960 bytes**(338 Mb)

The results were obtained from an in memory benchmark( the whole file was read to memory).

Machine
```text
Darwin MacBook-Pro.local 21.2.0 Darwin Kernel Version 21.2.0: Sun Nov 28 20:29:10 PST 2021; root:xnu-8019.61.5~1/RELEASE_ARM64_T8101 arm64
```
(M1 Macbook 2020)

|Encoder/Decoder|Runs|Time    |Speed     |Ratio      |
|---------------|----|--------|----------|-----------|
| Huffman Encode| 20 |4409ms  | 1466 Mb/s|1.4717 to 1|
| Huffman Decode| 20 |3775ms  | 1721 Mb/s|           |
| tANS Encode   | 20 |8829ms  | 732  Mb/s|1.5000 to 1|
| tANS Decode   | 20 |6415ms  | 1016 Mb/s|           |

#Acknowledgements
This couldn't have been made possible without the great works of 

- Jarek Duda - The inventor of ANS based algorithms, and collector, most of his
work and others can be found [here](https://encode.su/threads/2078-List-of-Asymmetric-Numeral-Systems-implementations)

- Eric Biggers - [xpack](https://github.com/ebiggers/xpack) (lz77+tANS) ,

- Yann Collet - [FiniteStateEntropy](https://github.com/Cyan4973/FiniteStateEntropy),

- Yann Collet - [FSE algorithm blog posts](http://fastcompression.blogspot.com/2013/12/finite-state-entropy-new-breed-of.html)

- Charles Bloom - [Understanding TANS](http://cbloomrants.blogspot.com/2014/01/1-30-14-understanding-ans-1.html) blog series, 


- Fabian Geisen - [Reading bits in far too many ways](https://fgiesen.wordpress.com/2018/02/19/reading-bits-in-far-too-many-ways-part-1/)

- Fabian Geisen -  [Interleaving rANS](https://arxiv.org/abs/1402.3392)