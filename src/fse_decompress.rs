use std::io::Read;

use crate::constants::{state_generator, TABLE_LOG, TABLE_SIZE};
use crate::fse_bitstream::FSEStreamReader;
use crate::huff_bitstream::BitStreamReader;
use crate::utils::Symbols;

fn spread_symbols(freq_counts: &[Symbols; 256]) -> [u32; TABLE_SIZE]
{
    // the start and the end of each state
    const INITIAL_STATE: usize = 0;

    let state_gen = state_generator(TABLE_SIZE);

    // if state is even,slots distribution won't work
    assert_eq!(state_gen & 1, 1, "State cannot be an even number");

    let mut state = INITIAL_STATE;

    let mut state_array = [Symbols::default(); TABLE_SIZE];

    // This table contains all previous cumulative state counts
    // sorted in ascending order from 0..non_zero
    let mut cumulative_state_counts = [0; 256];

    let mut slots = [0; 256];

    let mut c_count = 0;

    for sym in freq_counts.iter()
    {
        if sym.y == 0
        {
            continue;
        }

        /*
         * Allocate state to symbols
         */

        let symbol = sym.z;

        slots[symbol as usize] = sym.y;

        cumulative_state_counts[(symbol) as usize] = c_count;

        let mut count = sym.y;

        c_count += count;

        while count > 0
        {
            state_array[state].z = symbol;
            state_array[state].y = sym.y;

            state = (state + state_gen) & (TABLE_SIZE - 1);

            count -= 1;
        }
    }
    // Cumulative count should be equal to table size
    assert_eq!(
        c_count as usize, TABLE_SIZE,
        "Cumulative count is not equal to table size, internal error"
    );

    // state must be equal to the initial state, since its a cyclic state, otherwise we messed up
    assert_eq!(state, INITIAL_STATE, "Internal error, state is not zero");

    /*
     * Okay one more thing is now assigning bits to state
     * everyone agrees that this is the format
     *      Some threshold
     *           |
     *           v
     *|----------|-----------|
     *|min bits  | min_bits+1|
     *|----------|-----------|
     *
     * our state distribution.
     * During encoding, we assigned states based on normalized frequencies
     * and these are the same states that we used to get the number of bits for a symbol
     *
     * min_bits = TABLE_LOG-ceil(log2(assigned_counts));
     *
     * now with FSE, we sometimes encode as min_bits, sometimes as min_bits+1
     *
     * Each state now is associated with one of these choices.
     * Each state also contains destination_range of length 2^num_bits;
     * This is the state the encoder was prior to encoding this symbol and entering this state
     *
     * Okay to assign bits
     * the numerically first n slots are assigned `min_bits+1`, and all `slots-n` sockets are assigned
     * min_bits and are mapped to a destination range with states 0( this ensures after hitting a high state),
     * we hit a low state
     *
     *
     * Lets solve for n
     *
     *
     *	(2**(min_bits+1))n + (2**min_bits)(count - n) = num_states
     *	(2**min_bits)(2n + count - n) = num_states
     *	(2**min_bits)(n + count) = num_states
     *	n + count = num_states / (2**min_bits)
     *	n = num_states / (2**min_bits) - count
     *
     *
     *
     */
    for state in 0..TABLE_SIZE
    {
        let mut sym = &mut state_array[state];

        let counter = slots[sym.z as usize];

        slots[sym.z as usize] += 1;

        let num_bits = (TABLE_LOG as u32) - (15 - (counter | 1).leading_zeros());

        let destination_range_start = (counter << num_bits) - TABLE_SIZE as u16;

        // y stores start of next range
        sym.y = destination_range_start;

        // x stores number of bits
        sym.x = num_bits;
    }

    // pack into u32
    // symbol     -> 0..8 bits    [sym.z]
    // num_bits   -> 8..16 bits.  [sym.x]
    // next_state -> 16..32 bits. [sym.y]

    let next_state = state_array.map(|x| u32::from(x.y) << 16 | x.x << 8 | ((x.z) as u32));

    return next_state;
}

fn decode_symbols(src: &[u8], states: &[u32; TABLE_SIZE], dest: &mut [u8], block_size: usize)
{
    const SIZE: usize = 25;
    let mut stream = FSEStreamReader::new(src);

    // initialize states
    let (mut c1, mut c2, mut c3, mut c4, mut c5) = stream.init_states();

    let mut initial = [0_u8; 5];

    stream.align_decoder();

    macro_rules! decode_five {
        ($to:tt,$start:tt) => {
            stream.refill_fast();

            stream.decode_symbol(&mut c5, &mut $to[$start + 4], states);

            stream.decode_symbol(&mut c4, &mut $to[$start + 3], states);

            stream.decode_symbol(&mut c3, &mut $to[$start + 2], states);

            stream.decode_symbol(&mut c2, &mut $to[$start + 1], states);

            stream.decode_symbol(&mut c1, &mut $to[$start + 0], states);

        };
    }
    // so in our bad scheme, we added some bits that were unneeded.
    // so we need to know how many unneeded bits were written.
    // the unneeded bits are block size to the lower multiple of 5
    let rounded_down = ((block_size) / 5) * 5;

    unsafe {
        decode_five!(initial, 0);
    }

    let start = block_size-rounded_down;
    // now dest is aligned to a 5 byte boundary
    // let's goooooo
    let chunks = dest.get_mut(start..).unwrap().chunks_exact_mut(SIZE);

    unsafe {
        for chunk in chunks
        {
            decode_five!(chunk, 0);

            decode_five!(chunk, 5);

            decode_five!(chunk, 10);

            decode_five!(chunk, 15);

            decode_five!(chunk, 20);
        }
    }
    let remainder = dest
        .get_mut(start..)
        .unwrap()
        .chunks_exact_mut(SIZE)
        .into_remainder();

    for chunk in remainder.chunks_exact_mut(5)
    {
        unsafe {
            decode_five!(chunk, 0);
        }
    }

    dest[0..start].copy_from_slice(&initial[5-start..]);
}
fn read_headers(buf: &[u8], symbol_count: u8) -> [Symbols; 256]
{
    let mut symbols = [Symbols::default(); 256];

    let mut stream = BitStreamReader::new(buf);

    for _ in 0..(symbol_count / 2)
    {
        unsafe {
            stream.refill_fast();
        }
        let symbol = stream.get_bits(8) as usize;
        let state = stream.get_bits(11) as u16;
        symbols[symbol] = Symbols {
            z: symbol as i16,
            y: state,
            x: 0,
        };

        let symbol = stream.get_bits(8) as usize;
        let state = stream.get_bits(11) as u16;

        symbols[symbol] = Symbols {
            z: symbol as i16,
            y: state,
            x: 0,
        };
    }

    if (symbol_count & 1) != 0
    {
        // Do the last odd value
        unsafe {
            stream.refill_fast();
        }
        let symbol = stream.get_bits(8) as usize;
        let state = stream.get_bits(11) as u16;

        symbols[symbol] = Symbols {
            z: symbol as i16,
            y: state,
            x: 0,
        };
    }
    return symbols;
}
pub fn fse_decompress<R: Read>(src: &mut R, dest: &mut Vec<u8>)
{
    // read block information
    let mut block_info = [0];

    src.read_exact(&mut block_info).unwrap();

    let mut length = [0, 0, 0, 0];
    // read the length
    src.read_exact(&mut length[0..3]).unwrap();

    let mut block_length = u32::from_le_bytes(length);

    // not all bytes will be used
    let mut source = vec![0; (block_length + 200) as usize];

    // header cannot go beyond 600 bytes
    let mut header = vec![0; 600];

    let mut header_t = [0; 2];

    let mut compressed_size = [0; 4];

    let mut symbols_count = [0];

    loop
    {
        block_length = u32::from_le_bytes(length);

        src.read_exact(&mut header_t).unwrap();

        if dest.capacity() <= (block_length as usize + dest.len())
        {
            // no don't grow using armotized technique, this may become really big
            dest.reserve(block_length as usize);
        }

        let header_size = u16::from_le_bytes(header_t);

        src.read_exact(&mut symbols_count).unwrap();

        src.read_exact(&mut header[0..usize::from(header_size)])
            .unwrap();

        let symbols = read_headers(&header[0..usize::from(header_size)], symbols_count[0]);
        // reconstruct next_state
        let next_state = spread_symbols(&symbols);

        src.read_exact(&mut compressed_size[0..3]).unwrap();

        let compressed_length = (u32::from_le_bytes(compressed_size)) as usize;

        src.read_exact(&mut source[0..compressed_length]).unwrap();

        let start = dest.len();

        let new_len = dest.len() + block_length as usize + 5;

        unsafe { dest.set_len(dest.capacity()) };
        // set length to be the capacity
        if new_len > dest.capacity()
        {
            let cap = dest.capacity();
            dest.reserve_exact(new_len - cap);
        }
        unsafe { dest.set_len(new_len) };

        // Don't continue if we don't have capacity to create a new write
        assert!(
            new_len <= dest.capacity(),
            "{},{}",
            new_len,
            dest.capacity()
        );

        decode_symbols(
            &source[0..compressed_length],
            &next_state,
            dest.get_mut(start..).unwrap(),
            block_length as usize,
        );

        // we had to add a +5 length to allow for unused symbols to
        // be added to the state,
        // so do not forget to remove them since they contain dummy values.
        unsafe { dest.set_len(new_len - 5) };

        // we reached the last block.
        if (block_info[0] >> 6) & 1 == 1
        {
            break;
        }

        // not last block, pull in more bytes.
        src.read_exact(&mut block_info).unwrap();
        // read the length for the next iteration

        src.read_exact(length.get_mut(0..3).unwrap()).unwrap();
    }
}

#[test]
fn fse_decompress_test()
{
    use std::fs::read;
    use std::io::Cursor;

    use crate::fse_compress;

    let mut fs = Vec::with_capacity(1 << 24);

    let fd = read("/Users/calebe/Projects/Data/en").unwrap();

    fse_compress(&fd, &mut fs);

    let mut cursor = Cursor::new(fs);

    let mut vecr = vec![];

    fse_decompress(&mut cursor, &mut vecr)
}
