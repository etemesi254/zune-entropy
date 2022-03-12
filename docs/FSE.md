## FSE Format Documentation
Note, everything is Little L
### Header
1. **1 byte**- FSE block information.

| Bits   | Value                            |
|--------|----------------------------------|
| 5-7[1] | Identify the type of stream      |
| 4      | Identify if this is a last block |
| 0-3    | Currently nothing                |

The following types of streams are supported 

| Bit          | Type |
|--------------|------|
| Huffman      | 00   |
| RLE          | 01   |
| Uncompressed | 10   |
| FSE          | 11   |

This information is shared by the huffman section.

A  decoder must ensure that the first bits are 11 to ensure this is a valid FSE block.

A decoder should be able to handle different blocks inside it,
e.g a FSE decoder can handle decoding of RLE and uncompressed block but not Huffman block.

2. **3 bytes** - Block size

This indicates the total length of the symbols in this block, the set is 132kb(2^17) but can be anything between 
1 and 2^24(16 Mb).

The decoder can use this information to create an appropriate buffer to hold the decoded symbols.

After decoding, the decoder can ensure the decoded symbols match the decoded length.

3. **1 byte** - Information bit


| Bit  | Information                          |
|------|--------------------------------------|
| 4-8  | Maximum table log size               |
| 0-4  | Log2 size of the maximum state used  |


Maximum table log size is used by the decoder to set up the state table, e.g if it is 11, it means our state table
contains entries for 2048 (2^11) states, and if it's 5, means our state table contains entries for 32 (2^5) states.

Ranges between 5-11

The log2 size of the maximum state is used by the decoder while reading state bits from the stream, it tells the decoder
how many bits are used to encode normalized frequencies.

E.g. if the value retrieved is 10, all normalized frequencies will be encoded using 10 bits, for how they are laid out, see the 
header section.

5. **2 bytes** - tANS header Size

Indicates the header size for the tANS for this block, this size indicates how much data following it is to be
used to construct symbols with their respective slots as assigned by the encoder.

6. **1 byte**- Number of symbols stored in the headers

7. **n bytes**-Headers.

| Value                | Bits                             |
|----------------------|----------------------------------|
| Symbol               | 8 bits                           |
| Normalized frequency | lower 4 bits of information bits |

tANS headers contain symbol + normalized frequency of symbol.

The symbols are encoded using 8 bits, while the states are encoded using the value found in the lower 4 bits of the information bit(2 bytes above).

They are written in little endian from , from forwards.

After a decoder knows the size of the tANS headers, and the number of symbols, the decoding step becomes
read 8 bits, to get the symbol and read some more bits to get the normalized frequency of the symbol.

In the end, the headers return an array sorted in symbol order(0..255) containing symbols+normalized counts of the occuring
symbols.


8. **3 bytes** - Compressed size of the bitstream


9. **n bytes**-Compressed stream

The stream is encoded with 5 states, and the encoder ensures the buffer is padded up to 5 states, 
the decoder should discard extra padding added by the encoder.

Symbols are buffered up into CHUNK_SIZE(default is 132 kb) and encoded from backwards moving forwards, 
the stream is also written from backwards moving forwards.

Since ANS is a LIFO stack, the last symbol to encode is the first symbol to decode, meaning that the decoder
will write forwards.


