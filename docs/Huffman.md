## Huffman Format Documentation

1. **1 byte**- Huffman block information.

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

This information is shared by the FSE section.

A  decoder must ensure that the first bits are `00` to ensure this is a valid huffman block.

A decoder can also get a mixture of encoded blocks and RLE/Uncompressed blocks, such functions
should be handled with their respective functions.

A Huffman decoder should not receive a tANS block.


2. *3 bytes*- Total block size(Little Endian)

3. *10 bytes*- Jump table (Little Endian)

    2 bytes per table jump size. Cannot go above 2^16(65535)
4. **10 bytes**- Code lengths for symbols
   - There can be no symbol with code length 0
     `code_length[i]` represents codes of length `i+2`,
   so `code_length[9]` == codes of length 11
5. **n bytes** - Sum of the code lengths give the total number of symbols in ascending order
6. **n bytes** - Actual data for stream 1


The bits are stored  in little endian, with the codes assigned to a symbol bit reversed(like deflate).

Maximum code length is 11 bits, minimum is 1 bit.

