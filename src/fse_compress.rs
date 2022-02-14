use std::cmp::{max, min};
use std::io::Write;

use crate::bitstream::BitStreamWriter;
use crate::constants::{MAX_TABLE_LOG, TABLE_LOG, TABLE_SIZE};
use crate::utils::{histogram, Symbols};



const CHUNK_SIZE: usize = 1 << 17;

/// Scale up or down frequency occurrences from histogramming to fit in to.
/// fit into 2^TABLE_LOG(TABLE_SIZE). while ensuring every non-zero frequency gets a frequency
/// >=1, even for low probability symbols
///
/// ## Arguments
///- sym:  symbol results from histogramming.
///   After this function completes, each non_zero symbol will have its allocated slots at
///   Symbols.y, hence this function modifies all y values of sym
fn normalize_frequencies(sym: &mut [Symbols; 256], non_zero: usize, total: usize)
{
    const THRESH0LD: u16 = 4;
    // An upper limit on how many corrections we can go through.
    // if it goes above this, shorten the most common data values
    const RECURSION_LIMIT: usize = 30;

    let ratio = TABLE_SIZE as f32 / total as f32;

    let slice = sym.get_mut((255 & non_zero)..).unwrap();

    let mut summed_bits = 0;

    for sl in slice
    {
        if sl.x == 0
        {
            continue;
        }

        sl.y = max((ratio * (sl.x as f32)).round() as u16, 1);

        summed_bits += sl.y;
    }
    // correction part
    if (summed_bits as usize) < TABLE_SIZE
    {
        // some more slots are available

        // give the largest symbol
        sym[255].y += (TABLE_SIZE) as u16 - summed_bits;
        summed_bits += (TABLE_SIZE) as u16 - summed_bits;
    }
    else if (summed_bits as usize) > TABLE_SIZE
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
        let mut excess_states = summed_bits - (TABLE_SIZE as u16);

        let mut recursion_depth = 0;

        let mut index = 255;

        let mut second = sym[index - 1];

        let mut first = sym.get_mut(index).unwrap();

        loop
        {
            if first.y == 1
            {
                // Don't nuke a symbol with 1 slot, go start at the top.

                // reset index, and restart
                index = 255;

                second = sym[index - 1];

                first = sym.get_mut(index).unwrap();
            }

            // ensure that this operation at least reduces a slot(the +1)
            // The division is turned to a shift (strength reduction), don't try optimizing it
            let v = min(((first.y - second.y) / THRESH0LD) + 1, excess_states);

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
    assert_eq!(summed_bits as usize, TABLE_SIZE)
}
/// Allocate bits for each non zero frequency in the state
///
/// ## Modifies
/// This modifies the value `x` of every non zero frequency counts
/// storing a value that allows us to branchlessly determine if we are using
/// max_bits or max_bits-1
///
fn generate_state_bits(freq_counts: &mut [Symbols; 256], non_zero: usize)
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

    // TODO: Add next state stuff. before I forget.
    // Now I'm just cleaning the junk.

    for sym in &mut freq_counts[non_zero..].iter_mut()
    {
        // y contains actual slots
        // x -> histogram
        // replace histogram(x) with actual bits

        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            #[cfg(not(target_feature = "lzcnt"))]
            {
                // the |1 is because bsr instruction has undefined behaviour when src is zero.
                // so we have to ensure it's never zero for Rust to generate a better code length
                // the or doesn't affect anything since all items are greater than 1, and we only care about
                // the top bit.
                sym.x = (TABLE_LOG as u32) - (u16::BITS - (sym.y | 1).leading_zeros()) + 1;
            }
            #[cfg(target_feature = "lzcnt")]
            {
                // lzcnt instruction does not have undefined behaviour when it's input is zero.
                sym.x = (TABLE_LOG as u32) - (u16::BITS - (sym.y).leading_zeros()) + 1;
            }
        }
        #[cfg(not(any(target_arch = "x86", target_arch = "x86_64")))]
        {
            // everyone else defines leading zeroes correctly
            sym.x = (TABLE_LOG as u32) - (u16::BITS - (sym.y).leading_zeros()) + 1;
        }
        // store a value that allows us to determine if we are using max_bits or max_bits-1 bits
        // format
        // Remember
        // sym.y <- slots allocated to this symbol
        // sym.x <- max_bits

        //BUG: If slots exceeds 1<<(TABLE_LOG-1), this panics, fix that.
        sym.x = (sym.x << TABLE_LOG) - (u32::from(sym.y) << (sym.x));
    }
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
fn state_generator(num_states: usize) -> usize
{
    // play with this
    // Changing this , make sure compression ratio changes
    return (num_states >> 1) | (num_states >> 3) | 3;
}
fn spread_symbols(freq_counts: &mut [Symbols; 256]) -> [u16; 1 << TABLE_LOG]
{
    let state_gen = state_generator(TABLE_SIZE);

    assert_eq!(state_gen & 1, 1, "State cannot be an even number");

    let mut state = 0;
    let mut state_array = [0; TABLE_SIZE];

    // This table contains all previous cumulative state counts
    // sorted in ascending order from 0..non_zero
    let mut cumulative_state_counts = [0; 256];

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

        let symbol = sym.symbol;

        cumulative_state_counts[usize::from(symbol)] = c_count;

        let mut count = sym.y;

        c_count += count;

        while count > 0
        {
            state_array[state] = symbol;

            state = (state + state_gen) & (TABLE_SIZE - 1);

            count -= 1;
        }
    }
    // Cumulative count should be equal to table size
    assert_eq!(
        c_count as usize, TABLE_SIZE,
        "Cumulative count is not equal to table size, internal error"
    );

    // state must be zero at this point, since its a cyclic state, otherwise we messed up
    assert_eq!(state, 0, "Internal error, state is not zero");

    /*
     * Build next_state[].  This array maps symbol occurrences in the
     * state table, ordered primarily by increasing symbol value and
     * secondarily by increasing state, to their states, adjusted upwards by
     * num_states.
     */
    let mut next_state = [0; 1 << TABLE_LOG];

    for i in 0..TABLE_SIZE
    {
        let symbol = usize::from(state_array[i]);

        let pos = cumulative_state_counts[symbol];

        next_state[pos as usize] = (TABLE_SIZE + i) as u16;

        cumulative_state_counts[symbol] += 1;
    }

    return next_state;
}

fn encode_symbols<W: Write>(
    src: &[u8], symbols: &[Symbols; 256], next_states: &[u16; TABLE_SIZE], buf: &mut [u8],
    dest: &mut W,
)
{

    let mut stream = BitStreamWriter::new();

    let c = TABLE_SIZE as u16;

    const SIZE:usize = 25;

    // states for each interleaved streams
    let (mut c1, mut c2, mut c3, mut c4, mut c5) = (c, c, c, c, c);

    // chunk into
    for chunk in src.chunks_exact(SIZE)
    {
        macro_rules! encode_slice {
            ($start:tt) => {
                stream.encode_bits_fse(
                    chunk[$start].try_into().unwrap(),
                    symbols,
                    next_states,
                    &mut c1,
                );

                stream.encode_bits_fse(
                    chunk[$start+1].try_into().unwrap(),
                    symbols,
                    next_states,
                    &mut c2,
                );

                stream.encode_bits_fse(
                    chunk[$start+2].try_into().unwrap(),
                    symbols,
                    next_states,
                    &mut c3,
                );

                stream.encode_bits_fse(
                    chunk[$start+3].try_into().unwrap(),
                    symbols,
                    next_states,
                    &mut c4,
                );

                stream.encode_bits_fse(
                    chunk[$start+4].try_into().unwrap(),
                    symbols,
                    next_states,
                    &mut c5,
                );

            stream.flush_fast(buf);
            };
        }
        unsafe {
            encode_slice!(0);
            encode_slice!(5);
            encode_slice!(10);
            encode_slice!(15);
            encode_slice!(20);

        }
    }

    dest.write_all(&buf[0..stream.get_position()])
        .expect("Failed to write to destination buffer");
}
pub fn fse_compress<W: Write>(src: &[u8], dest: &mut W)
{
    /*
     * FSE compression, as usual, the steps
     *
     * 1.  Histogram
     * 2.  Normalize counts
     * 3.  Create Compression table
     * 4. Compress using  table.
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

    let mut is_last = false;

    // start is our current pointer to where we start compressing from
    let mut start = 0;
    // end is our current pointer to where we stop compressing at
    let mut end = min(CHUNK_SIZE, src.len());

    let mut src_chunk = &src[start..end];

    let mut buf = vec![0; CHUNK_SIZE + 10];
    while !is_last
    {
        /*
         * Loop until data is exhausted,
         */
        if end == src.len()
        {
            is_last = true;
        }
        let mut freq_counts = histogram(src_chunk);
        // sort by occurrences
        freq_counts.sort_unstable_by(|a, b| a.x.cmp(&b.x));
        //find first non-zero element
        let non_zero = freq_counts.iter().position(|x| x.x != 0).unwrap_or(0);

        // normalize frequencies.
        normalize_frequencies(&mut freq_counts, non_zero, src_chunk.len());

        // Generate actual bit codes for the states for symbols.
        generate_state_bits(&mut freq_counts, non_zero);
        // spread symbols

        let mut new_counts = [Symbols::default(); 256];

        for i in freq_counts
        {
            new_counts[i.symbol as usize] = i;
        }
        let next_states = spread_symbols(&mut new_counts);

        encode_symbols(src_chunk, &new_counts, &next_states, &mut buf, dest);

        // read the next block
        start = end;
        // end increases either by chunk size or points to end of buffer.
        end = min(end + CHUNK_SIZE, src.len());

        src_chunk = &src[start..end];
    }
}

#[test]
fn fse_compress_test()
{
    use std::fs::{read, OpenOptions};
    use std::io::BufWriter;

    let fs = OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true)
        .open("/Users/calebe/CLionProjects/zcif/tests.zcf")
        .unwrap();
    let mut fs = BufWriter::with_capacity(1 << 24, fs);

    let fd = read("/Users/calebe/Projects/Data/enwiki.smaller").unwrap();

    fse_compress(&fd, &mut fs);
}
