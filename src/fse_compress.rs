//! FSE compression.
//!
//! Format
//! See FSE.md docs for format
//!
use std::cmp::{max, min};
use std::io::Write;

use log::Level::Trace;
use log::{log_enabled, trace};

use crate::constants::{
    state_generator, CHUNK_SIZE, EXTRA, MAX_TABLE_LOG, MIN_TABLE_LOG, TABLE_LOG, TABLE_SIZE,
};
use crate::errors::EntropyErrors;
use crate::fse_bitstream::FseStreamWriter;
use crate::huff_bitstream::BitStreamWriter;
use crate::utils::{histogram, write_rle, write_uncompressed, Symbols};

const SIZE: usize = 25;

const THRESH0LD: u32 = 1024;

///
/// Scale up or down frequency occurrences from histogramming to fit in to.
/// fit into `2^TABLE_LOG(TABLE_SIZE`). while ensuring every non-zero frequency gets a frequency
/// >=1, even for low probability symbols
///
/// ## Arguments
///- sym:  symbol results from histogramming.
///   After this function completes, each non zero symbol will have its allocated slots at
///   Symbols.y, hence this function modifies all y values of sym
fn normalize_frequencies_fast(
    sym: &mut [Symbols; 256], tbl_size: usize, non_zero: usize, total: usize,
)
{
    // If two values sym[x] and sym[x-1] are separated by THRESHOLD
    // reduce sym[x] by (sym[x]-sym[x-1])/THRESHOLD.
    const THRESH0LD: u16 = 16;
    // An upper limit on how many corrections we can go through.
    // if it goes above this, shorten the most common data values
    const RECURSION_LIMIT: usize = 30;

    let ratio = tbl_size as f32 / total as f32;

    let slice = sym.get_mut((255 & non_zero)..).unwrap();

    let mut summed_bits = 0;

    for sl in slice
    {
        if sl.x == 0
        {
            continue;
        }

        // the number of slots is given by symbol occurrence (sl.x)
        // multiplied by the ratio(TABLE_SIZE/source length),rounded to the nearest number,
        // every one rounded to zero should be rounded up to 1.
        sl.y = max((ratio * (sl.x as f32)).floor() as u16, 1);

        summed_bits += sl.y;
    }
    // correction part
    if (summed_bits as usize) < tbl_size
    {
        // some more slots are available

        // give the largest symbol
        sym[255].y += (tbl_size) as u16 - summed_bits;

        summed_bits += (tbl_size) as u16 - summed_bits;
    } else if (summed_bits as usize) > tbl_size
    {
        /* do some hacky heuristic's
         * It isn't the best heuristic but it produces better output than a naive heuristic.
         *
         * The ideal goal is slot distribution.  We want to reduce slots for largely occurring symbols
         * to make space for those symbols rounded up to 1.
         *
         * We want a good slot distribution algorithm, we don't want to just massively take all slots
         * from the symbol with maximum slots, it will hurt compression
         *
         * So what heuristic did I choose?
         * -> Decrease a slot by (slot[i].available_slots - slot[i-1].available_slots)/threshold
         * which isn't actually a bad heuristic since it ensures a slot[i] will always have larger than slot[i-1]
         * symbols.
         *
         * -> Someone may notice that this may have a bad behaviour where it can recur forever if
         *    slots[i] and slots[i-1] aren't largely spaced
         *   The solution is to have a threshold, and if we cross that threshold during slot correction, just nuke the
         *   top elements to ensure the slots == TABLE_SIZE and call it a day. Not optimal but works
         */
        let mut excess_states = summed_bits - (tbl_size as u16);

        let mut recursion_depth = 0;

        let mut index = 255;

        let mut second = sym[index - 1];

        let mut first = sym.get_mut(index).unwrap();

        trace!("Excess states:{}", excess_states);

        loop
        {
            if first.y >=1
            {
                // Don't nuke a symbol with 1 slot, go start at the top.

                // reset index, and restart
                index = 255;

                second = sym[254];

                first = sym.get_mut(255).unwrap();
            }

            // ensure that this operation at least reduces a slot(the +1)
            let v = min((first.y.saturating_sub( second.y)) / THRESH0LD, excess_states);

            excess_states -= v;

            summed_bits -= v;

            first.y -= v;

            if excess_states == 0
            {
                // Done :)
                break;
            }
            index -= 1;

            second = sym[index - 1];

            first = sym.get_mut(index).unwrap();

            recursion_depth += 1;

            if recursion_depth == RECURSION_LIMIT
            {
                // okay we tried

                // nuke the best, second and third best

                let div = excess_states / 3;

                let rem = excess_states % 3;

                sym.get_mut(255).unwrap().y -= div + rem;

                sym.get_mut(254).unwrap().y -= div;

                sym.get_mut(253).unwrap().y -= div;

                summed_bits -= (div * 3) + rem;

                break;
            }
        }
    }
    assert_eq!(summed_bits as usize, tbl_size);
}

/// Allocate bits for each non zero frequency in the state
///
/// # Modifies
/// This modifies the value `x` of every non zero frequency counts
/// storing a value that allows us to branchlessly determine if we are using
/// `max_bits` or max_bits-1
///
/// # Returns
/// A number indicating how compressible this block is
/// This number  divided by the total block size will give the compression
/// ratio of the current block(in float). But since divisions are expensive,
/// you can simply subtract from a certain threshold to see if this block can be compressed,
fn generate_state_bits(freq_counts: &mut [Symbols; 256], non_zero: usize, table_log: usize) -> u32
{
    /*
     * Generate max_bits for every non-zero frequency in freq_counts
     * Okay this was done in the notebooks but i'll still do it
     * Ps-> Probability of symbol 's' occurring
     *
     * entropy = log2(1/Ps)
     * entorpy = log2(1/(actual_slots/TABLE_SIZE)
     *
     * Inverse the values in the brackets
     *
     * entropy = log2(TABLE_SIZE/actual_slots)
     *
     * apply log laws
     *
     * entropy = log2(TABLE_SIZE) -  log2(actual_slots)
     *
     * entropy = LOG_TABLE_SIZE - log2(actual_slots)
     *
     * But integer log2 can be implemented by a leading zeros count
     * (builtin_clz in gcc/clang) see https://stackoverflow.com/questions/11376288/fast-computing-of-log2-for-64-bit-integers
     *
     *
     *
     * Do note we do not store actual state symbols but we store a value that guarantees
     * us to branchlessly determine if we are using max_bits or max_bits-1
     */

    let mut compressibility = 0;

    for sym in &mut freq_counts[non_zero..].iter_mut()
    {
        let hist = sym.x;
        // y contains actual slots
        // x -> histogram
        // replace histogram(x) with actual bits

        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            #[cfg(not(target_feature = "lzcnt"))]
            {
                // the |1 is because bsr instruction has undefined behaviour when src is zero.
                // so we have to ensure it's never zero for Rust to generate a better code.
                // the or doesn't affect anything since all items are greater than or equal to 1,
                // and we only care about the top bit.
                sym.x = (table_log as u32) - (u16::BITS - (sym.y | 1).leading_zeros()) + 1;
            }
            #[cfg(target_feature = "lzcnt")]
            {
                // lzcnt instruction does not have undefined behaviour when it's input is zero.
                sym.x = (table_log as u32) - (u16::BITS - (sym.y).leading_zeros()) + 1;
            }
        }
        #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
        {
            // everyone else defines leading zeroes correctly
            sym.x = (table_log as u32) - (u16::BITS - (sym.y).leading_zeros()) + 1;
        }

        // store a value that allows us to determine if we are using max_bits or max_bits-1 bits

        // Remember
        // sym.y <- slots allocated to this symbol
        // sym.x <- max_bits

        compressibility += sym.x * hist;
        sym.x = (sym.x << MAX_TABLE_LOG) - (u32::from(sym.y) << (sym.x));
    }
    // compressibility contains the theoretical bits that will be used to represent the symbols
    compressibility
}

/// Generate a pseudo state which we will use to spread symbols
///
///
/// State generator must emit a cyclic property , in that a state can go from 0 to DEST_TABLE_SIZE-1
/// without repeating any state before it reaches a state that isn't repeated.
///
/// Okay a cyclic state is (mod n), e.g for a table of size 6, we can have state values
/// ```text
/// 0,1,2,3,4,5 | 0,1,2,3,4,5
/// [---------] | [ --------]
/// [    n    ] | [     n   ]
/// ```
///
/// But in TANS, a state like that is bound to lead to bad compression.
/// e.g AAAABBBBCCC is a bad symbol spread, since state transitions would be expensive,
/// you want ABCABCABABC kinda state,  so we need a good strategy for symbol spread,
/// there are many techniques, sorting, precise symbol spread and a lot of other techniques
/// but this lazy state generator is actually good.
fn spread_symbols(
    freq_counts: &mut [Symbols; 256], table_size: usize,
) -> Result<([u16; TABLE_SIZE], [i32; 256]), EntropyErrors>
{
    /*
     * Ps: table size is (1<<table_log)
     *
     * Ps you can store cumulative counts and build a decoding table
     * from there
     * Should I do that?
     */
    // the start and the end of each state
    const INITIAL_STATE: usize = 0;

    let state_gen = state_generator(table_size);

    // if state is even,slots distribution won't work
    assert_eq!(state_gen & 1, 1, "State cannot be an even number");

    let mut state = INITIAL_STATE;

    let mut state_array = [0; TABLE_SIZE];

    // This table contains all previous cumulative state counts
    // sorted in ascending order from 0..non_zero
    let mut cumulative_state_counts = [0; 256];

    let mut c_count = 0;

    let mut next_state_offset = [0; 256];

    for sym in freq_counts.iter_mut()
    {
        if sym.y == 0
        {
            continue;
        }
        // TODO-> Give symbols with a single probability closer to state

        /*
         * In the spread symbols, we allocated states to symbols, now we want to do the opposite,
         * allocate symbols to states
         */

        let symbol = sym.z;

        cumulative_state_counts[((symbol) as usize)] = c_count;

        let mut count = sym.y;

        next_state_offset[(sym.z) as usize] += i32::from(c_count) - i32::from(count);

        c_count += count;

        while count > 0
        {
            state_array[state] = symbol;

            state = (state + state_gen) & (table_size - 1);

            count -= 1;
        }
    }


    // Cumulative count should be equal to table size
    if c_count as usize != table_size
    {
        return Err(EntropyErrors::CorruptHeader(format!(
            "Cumulative count is not equal to table size, error. Cumulative count:{},table_size:{}",
            c_count, table_size
        )));
    }
    if state != INITIAL_STATE
    {
        // state must be equal to the initial state, since its a cyclic state, otherwise we messed up

        return Err(EntropyErrors::CorruptHeader(format!(
            "Internal error, cyclic state is not back to initial state state:{},expected:{}",
            state, INITIAL_STATE
        )));
    }

    /*
     * Build next_state[].  This array maps symbol occurrences in the
     * state table, ordered primarily by increasing symbol value and
     * secondarily by increasing state, to their states, adjusted upwards by
     * num_states.
     */
    let mut next_state = [0; TABLE_SIZE];

    for (i, symbol) in state_array.iter().take(table_size).enumerate()
    {
        let symbol = (*symbol) as usize;

        let pos = cumulative_state_counts[symbol];

        next_state[pos as usize] = (table_size + i) as u16;

        cumulative_state_counts[symbol] += 1;
    }
    Ok((next_state, next_state_offset))
}

/// Write headers, packing symbols and their state counts
/// into a bitstream format
///
/// The format is
/// |symbol| -> 8 bits
///
/// |State counts| -> `max_state_bits(see` below for explanation)
fn write_headers<W: Write>(
    symbols: &[Symbols; 256], non_zero: usize, block_size: usize, table_log: usize,
    last_block: bool, dest: &mut W,
) -> Result<(), EntropyErrors>
{
    /*
     * Okay we have a header problem
     * what do we need to re-construct next_state?
     *256
     * Well we need to know how many slots each symbol was given
     * and what symbol it was
     * So what to do?
     * okay let's try to give it the following construction
     *
     *
     * | symbol | slot |
     * But we know symbol is between 0..255
     * and slot is between 0..TABLE_SIZE
     * so we can use a bit-packing convention
     * slot occupies 0..max_state_bits bits.
     *
     * Okay max_states_bits is a bit complex, so we have a sorted
     * array of symbols, with the maximum state being in the 255th bit
     * So to know how many bits to use to encode states, simply take the
     * number of bits for the maximum symbol.
     * E.g if maximum symbol was given 255 slots, we know all other symbols
     * were given less than 255 slots, so can be represented with log2(255)
     * bits=8 bits, now all states will be encoded in 8 bits
     */

    let info_bit = 1 << 7 | 1 << 6 | u8::from(last_block) << 5;

    let chunks = symbols[non_zero..].chunks_exact(2);

    let remainder = chunks.remainder();
    // we can determine the size of the output buffer.
    // each symbol + state occupies up to 18 bits.
    // so output == (symbol * 2 bytes)+2 bits for each symbol.
    // So overallocate
    let mut output = vec![0_u8; (260 - non_zero) * 4];

    let mut stream = BitStreamWriter::new(&mut output);
    let state_sym = symbols[non_zero..]
        .iter()
        .max_by(|a, b| a.y.cmp(&b.y))
        .unwrap();

    // this won't go Past 10 bits
    let state_bits = (u16::BITS - state_sym.y.leading_zeros()) as u8;

    assert!(state_bits <= 11);

    macro_rules! encode_single {
        ($symbol:tt,$state:tt) => {
            stream.add_bits($symbol, u8::BITS as u8);
            stream.add_bits($state, state_bits as u8);
        };
    }

    let (mut symbol, mut state);
    for chunk in chunks
    {
        symbol = chunk[0].z as u64;
        state = u64::from(chunk[0].y);

        encode_single!(symbol, state);

        symbol = chunk[1].z as u64;
        state = u64::from(chunk[1].y);

        encode_single!(symbol, state);

        unsafe {
            stream.flush_fast();
        }
    }
    for chunk in remainder
    {
        symbol = chunk.z as u64;
        state = u64::from(chunk.y);

        encode_single!(symbol, state);

        unsafe {
            stream.flush_fast();
        }
    }
    unsafe {
        stream.flush_final();
    }
    // two bytes Little endian, indicate this header size
    let header_size = stream.get_position();

    trace!("tANS Header size: {} bytes", header_size);
    // write info bit
    dest.write_all(&[info_bit])?;
    // write block size
    dest.write_all(&block_size.to_le_bytes()[0..3])?;

    // write number used to do state bits
    let compact_log = (table_log << 4) as u8;
    // + number for table log.
    dest.write_all(&[compact_log | state_bits])?;

    // header size
    dest.write_all(&header_size.to_le_bytes()[0..2])?;

    // write number of symbols
    dest.write_all(&[(255 - non_zero) as u8])?;

    dest.write_all(stream.get_output())?;

    Ok(())
}

/// Calculate the maximum state log bits
/// to use for this block
fn max_log(src_size: usize) -> usize
{
    let high_bits = usize::BITS - (src_size - 1).leading_zeros();

    max(min(TABLE_LOG, high_bits as usize), MIN_TABLE_LOG)
}

#[inline(always)]
fn encode_symbols_fallback<W: Write>(
    src: &[u8], common_symbol: i16, symbols: &mut [Symbols; 256], table_size: usize,
    next_states: &[u16; TABLE_SIZE], next_states_offset: &[i32; 256], buf: &mut [u8], dest: &mut W,
) -> Result<(), EntropyErrors>
{
    /*
     * The one thing to keep in mind is that how different this is from Huffman encoding
     * in tANS/FSE interleaving is done by adding a new state variable
     * while in Huffman, we use different streams with chains in
     *
     * So another good question is why is MAX_TABLE_LOG 11?
     *
     * It's a performance choice.
     *
     * That's how we get an edge on FSE.
     * With 11 states and 5 interleaved streams, we have a maximum
     * of 55 bits to be written and our bit-buffer can hold 63 bits max.
     * A flush operation can leave at most 7 bits so we cannot write past (63-7)
     * = 55 bits branchlessly.
     * And 11*5 streams gives the minimum instruction/symbol.
     *
     *
     * What about a buffer not aligned to multiple of 5?
     * Fill it with the most common symbol, we can do with additional symbols
     * and it allows for a clean decoder loop.
     */
    // place next_state inside symbol definition
    let mut pos = 0;
    let symbols = &symbols.map(|mut x| {
        x.z = next_states_offset[pos] as i16;
        pos += 1;
        x.to_u64()
    });

    let c = table_size as u16;

    // states for each interleaved streams
    let (mut c1, mut c2, mut c3, mut c4, mut c5) = (c, c, c, c, c);

    let mut stream = FseStreamWriter::new(buf);
    macro_rules! encode_slice {
        ($start:tt,$chunk:tt) => {
            stream.encode_symbol($chunk[$start + 0], symbols, next_states, &mut c1);

            stream.encode_symbol($chunk[$start + 1], symbols, next_states, &mut c2);

            stream.encode_symbol($chunk[$start + 2], symbols, next_states, &mut c3);

            stream.encode_symbol($chunk[$start + 3], symbols, next_states, &mut c4);

            stream.encode_symbol($chunk[$start + 4], symbols, next_states, &mut c5);

            stream.flush_fast();
        };
    }

    for chunk in src.rchunks_exact(SIZE)
    {
        // do some unrolling
        unsafe {
            encode_slice!(20, chunk);
            encode_slice!(15, chunk);
            encode_slice!(10, chunk);
            encode_slice!(5, chunk);
            encode_slice!(0, chunk);
        }
    }

    {
        // deal with symbols that were left over
        // (not divisible by SIZE)

        let rem_chunks = src.rchunks_exact(SIZE).remainder();

        let mut start = 5 - (rem_chunks.len() % 5);

        if rem_chunks.len() % 5 == 0
        {
            // if chunk is divisible by 10, don't add 5 dummy zeros
            // if it is divisible by 5, is it okay to add?
            start = 0;
        }
        // duplicate  the common symbol
        let mut new_loc = vec![common_symbol as u8; rem_chunks.len() + start];

        let end = new_loc.len() - rem_chunks.len();

        new_loc[end..].copy_from_slice(rem_chunks);

        // new_loc looks like
        // a,a,a,a,[b,c,s,d] (`a`,is the most common symbol)
        // duplicate values are found in the start since encoding  moves from
        // the end to the beginning.
        // new loc is divisible by 5 since we have 5 states.

        // do the final encoding.
        for chunk in new_loc.rchunks_exact(5)
        {
            unsafe {
                encode_slice!(0, chunk);
            };
        }
    }

    stream.flush_final();

    // you may think it's over, but that is where you go wrong
    // we need to write the final state values

    stream.encode_final_states(c1, c2, c3, c4, c5, table_size);
    // write size of this block
    dest.write_all(&stream.get_position().to_le_bytes()[0..3])?;

    dest.write_all(stream.get_output())?;

    // Logging information
    if log_enabled!(Trace)
    {
        trace!("Size of original block: {}", src.len());
        trace!("Size of compressed block: {}", stream.get_position());
        trace!(
            "Ratio :{:.6}",
            (stream.get_position() as f32) / (src.len() as f32)
        );
        println!();
    }
    Ok(())
}

#[inline(always)]
#[rustfmt::skip]
fn encode_symbols<W: Write>(
    src: &[u8], common_symbol: i16, symbols: &mut [Symbols; 256], table_size: usize,
    next_states: &[u16; TABLE_SIZE], next_states_offset: &[i32; 256], buf: &mut [u8], dest: &mut W,
) -> Result<(), EntropyErrors>
{
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    {
        if is_x86_feature_detected!("bmi2")
        {
            return unsafe {
                encode_symbols_bmi(
                    src, common_symbol,
                    symbols, table_size, next_states,
                    next_states_offset,
                    buf, dest,
                )
            };
        }
    }
    encode_symbols_fallback(
        src, common_symbol,
        symbols, table_size,
        next_states, next_states_offset,
        buf, dest)
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[target_feature(enable = "bmi2")]
#[rustfmt::skip]
unsafe fn encode_symbols_bmi<W: Write>(
    src: &[u8], common_symbol: i16, symbols: &mut [Symbols; 256], table_size: usize,
    next_states: &[u16; TABLE_SIZE], next_states_offset: &[i32; 256], buf: &mut [u8], dest: &mut W,
) -> Result<(), EntropyErrors>
{
    encode_symbols_fallback(
        src,
        common_symbol,
        symbols,
        table_size,
        next_states,
        next_states_offset,
        buf,
        dest,
    )
}

pub fn fse_compress(src: &[u8], dest: &mut Vec<u8>) -> Result<(), EntropyErrors>
{
    /*
     * FSE compression, as usual, the steps
     *
     * 1.  Histogram
     * 2.  Normalize counts
     * 3.  Create Compression table
     * 4.  Compress using table.
     *
     * -----------------------------------------------------------
     * 1. Histogramming.
     *  Simply count  number of occurrences of a byte in a distribution.
     *
     * -----------------------------------------------------------
     * 2. Normalizing counts
     *
     * Maintain the constraint that the sum of the frequencies for histogramming
     * equal 2^TABLE_LOG, and ensure all symbols with a probability of >0 get a frequency >=1.
     * Giving some slots of high probable symbols to the low probable ones.
     *
     * ------------------------------------------------------------
     * 3. The next step is actually creating the compression table, where we first
     * distribute slots to symbols, with a crazy symbol distribution algorithm.
     * allocate bits to each symbol and finally, finally encode
     *
     */
    // asserts-> all are static

    assert!(TABLE_LOG <= MAX_TABLE_LOG);
    assert!(!src.is_empty(), "Oo oo");

    // A lot of invariants won't work if this isn't obeyed
    assert!(TABLE_SIZE.is_power_of_two());

    let mut is_last = false;

    // start is our current pointer to where we start compressing from
    let mut start = 0;
    // end is our current pointer to where we stop compressing at
    let mut end = min(CHUNK_SIZE, src.len());

    let mut src_chunk = &src[start..end];

    let mut buf = vec![0; CHUNK_SIZE + EXTRA];

    while !is_last
    {
        /*
         * Loop until data is exhausted,
         */

        let tbl_log = max_log(src_chunk.len());

        let tbl_size = 1_usize << tbl_log;

        if end == src.len()
        {
            is_last = true;
        }
        let mut freq_counts = histogram(src_chunk);
        // sort by occurrences
        freq_counts.sort_unstable_by(|a, b| a.x.cmp(&b.x));

        if freq_counts[255].x == src_chunk.len() as u32
        {
            // rle block
            trace!("Encountered RLE block, emitting as RLE");

            write_rle(src_chunk, dest, is_last);

            // read the next block
            // Todo, replace with goto when it becomes stable
            start = end;
            // end increases either by chunk size or points to end of buffer.
            end = min(end + CHUNK_SIZE, src.len());

            src_chunk = &src[start..end];
            // next block
            continue;
        }
        //find first non-zero element
        let non_zero = freq_counts.iter().position(|x| x.x != 0).unwrap_or(0);

        // normalize frequencies.

        normalize_frequencies_fast(&mut freq_counts, tbl_size, non_zero, src_chunk.len());

        // Generate actual bit codes for the states for symbols.
        // Do also the equivalent of detecting if this block is compressible.
        // returning how many bytes will be shaved off when we compress
        let compressibility = generate_state_bits(&mut freq_counts, non_zero, tbl_log);

        if compressibility + THRESH0LD > (src_chunk.len() as u32 * u8::BITS)
        {
            trace!("Block size:{}", src_chunk.len());
            trace!(
                "Theoretical compression ratio:{}",
                compressibility as f32 / src_chunk.len() as f32
            );
            trace!("Emitting block as uncompressed");

            write_uncompressed(src_chunk, dest, is_last);
        } else {
            let mut new_counts = [Symbols::default(); 256];

            // undo the sorting the earlier one did
            for i in freq_counts
            {
                new_counts[i.z as usize] = i;
            }
            // write headers expects sorted symbols
            write_headers(
                &freq_counts,
                non_zero,
                src_chunk.len(),
                tbl_log,
                is_last,
                dest,
            )?;

            let common_symbol = freq_counts[255].z;
            // dbg!(common_symbol);
            // spread symbols, generating next_state also
            // spread symbols requires symbols arranged in alphabetical order
            // to generate cumulative frequency
            //
            let (next_states, next_states_offset) = spread_symbols(&mut new_counts, tbl_size)?;

            encode_symbols(
                src_chunk,
                common_symbol,
                &mut new_counts,
                tbl_size,
                &next_states,
                &next_states_offset,
                &mut buf,
                dest,
            )?;
        }
        // read the next block
        start = end;
        // end increases either by chunk size or points to end of buffer.
        end = min(end + CHUNK_SIZE, src.len());

        src_chunk = &src[start..end];
    }
    Ok(())
}

