# Other Blocks possible

## 1. RLE

Format

| Bytes   | Use                                                           |
|---------|---------------------------------------------------------------|
| 1 byte  | Information Bit,identifies an RLE block, format `0b0100_0000` |
| 3 bytes | Little Endian, Identifies block size                          |
| 1 byte  | RLE byte                                                      |

The decoder will replicate the RLE byte block_size times (`memset`)

## 2. Uncompressed

| Bytes   | Use                                                           |
|---------|---------------------------------------------------------------|
| 1 byte  | Information Bit,identifies an RLE block, format `0b1000_0000` |
| 3 bytes | Little Endian, Identifies block size                          |
| n bytes | Uncompressed bytes in this block.                             |

The decoder will copy `block_size` bytes from the source to the resulting destination buffer