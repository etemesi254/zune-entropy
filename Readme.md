### MSB and LSB

MSB is how we think about numbers.

LSB is how little-endian cpus think about numbers.

###### MSB

The number 256 , broken to bits -> `[1,0,0,0,0,0,0,0]`

###### LSB

The number 256 broken to bits. ->  `[0,0,0,0,0,0,0,1]`

## Huffman Internals

The one thing to understand is how this Huffman implementation works it is LSB based when reading from literals, and MSB
based on everything else

Because sorting is MSB and not LSB

## Why do it LSB instead of MSB?

Well to save instructions.

MSB has a weird refill strategy in the bitstream state.

It has a weird dependency when it comes to refilling.

But we need a very fast refill strategy because its in the hot path

|MSB|LSB|
|---|---|
|Bad refill Strategy| Amazing refill strategy|
|Easy to build tables(and fast too)| Strided accesses suck|

We should talk about MSB table access

```rust
const LIMIT: usize = 11;

fn build_msb(table: &mut [u16; 1 << LIMIT], code_lengths: &[u8; 11], symbols: &[u8]) {
    let mut code = 0;

    let mut p = 0;

    for i in 1..=LIMIT
    {
        for _ in 0..code_lengths[i]
        {
            let look_bits = code << (LIMIT - i);

            for k in 0..(1 << (LIMIT - i))
            {
                let entry = &mut table[look_bits + k];
                // store symbols in top 8 bits and length in bottom 8 bits
                *entry = u16::from(symbols[p]) << 8 | (i as u16);
            }
            p += 1;

            code += 1;
        }
        assert!(code <= 1 << i);

        code *= 2;
    }
}
```

The important thing here is to notice how cache/memory and SIMD friendly this operation is ideally the inner loop is a
convoluted `memset` (you can force the compiler to emit this using `fill` method on the array).

Now switching to LSB table initialization kinda sucks

```rust

const LIMIT: usize = 11;

fn build(table: &mut [u16; 1 << LIMIT], code_lengths: &[u8; 11], symbols: &[u8]) {
    let mut code = 0;

    let mut p = 0;

    for i in 1..=LIMIT
    {
        for _ in 0..code_lengths[i]
        {
            let look_bits = code << (LIMIT - i);

            for k in 0..(1 << (LIMIT - i))
            {
                // the top 5-16 bits contain the reversed stream, so we need
                // to remove the bottom ones.
                let reversed_bits =
                    (REVERSED_BITS[(look_bits + k) as usize] >> (16 - LIMIT)) as usize;

                let entry = &mut table[reversed_bits];

                *entry = u16::from(symbols[p]) << 8 | (i as u16);
            }
            p += 1;

            code += 1;
        }
        // do not go one code length higher
        assert!(code <= 1 << i);
        
        code *= 2;
    }
}
```

Well the LSB form has a nasty inner loop. By nasty I mean how tied up those operations are

First we have a weird REVERSED_BITS array lookup., then we find the entry and then we populate it.
