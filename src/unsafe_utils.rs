/// especially when dealing with fast decompressor.
pub fn extend(vec: &mut Vec<u8>, new_len: usize)
{
    // add space for the new decompressed block
    vec.reserve_exact(new_len);
    // okay begin cheating.
    assert!(new_len <= vec.capacity());
    // use uninitialized memory
    // SAFETY
    // 1. Enough space is there since we reserved above(and the assert)
    // 2. Initialized memory is upheld later
    unsafe { vec.set_len(new_len) }
}
