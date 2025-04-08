import numpy as np, lattice_symmetries as ls, pytest, hypothesis, hypothesis.strategies as st, hypothesis.extra.numpy

phases = (hypothesis.Phase.explicit, hypothesis.Phase.reuse, hypothesis.Phase.generate)
rng = np.random.default_rng(seed=123)

def test_offset_ranges():
    from lattice_symmetries.compiler import _offset_ranges

    representatives = np.array([0, 1, 2, 4, 8, 9, 16, 17], dtype=np.uint64)
    number_bits, shift = 3, 61
    offsets, range_size = _offset_ranges(representatives, number_bits, shift)
    np.testing.assert_equal(offsets, [0, 0, 0, 0, 0, 0, 0, 0, 8])
    assert range_size == 8
    
    number_bits, shift = 3, 2
    offsets, range_size = _offset_ranges(representatives, number_bits, shift)
    np.testing.assert_equal(offsets, [0, 3, 4, 5, 5, 5, 5, 5, 8])
    assert range_size == 3

    representatives = np.array([0, 1, 2, 3], dtype=np.uint64)
    number_bits = 0
    for shift in range(64):
        offsets, range_size = _offset_ranges(representatives, number_bits, shift)
        assert len(offsets) == 2
        assert offsets[0] == 0
        assert offsets[1] == len(representatives)
        assert range_size == len(representatives)

    representatives = np.array([42], dtype=np.uint64)
    number_bits, shift = 2, 4
    offsets, range_size = _offset_ranges(representatives, number_bits, shift)
    np.testing.assert_equal(offsets, [0, 0, 0, 0, 1])
    assert range_size == 1

    # Test with numbers that have large gaps between them
    representatives = np.array([0, (1 << 32), (2 << 32), (3 << 32)], dtype=np.uint64)
    number_bits, shift = 2, 32
    offsets, range_size = _offset_ranges(representatives, number_bits, shift)
    assert len(offsets) == 5  # 2^2 + 1
    assert offsets[0] == 0
    assert offsets[1] == 1
    assert offsets[2] == 2
    assert offsets[3] == 3
    assert offsets[4] == 4
    assert range_size == 1

def test_explicit():
    # Test with small numbers first
    reps = np.array([0, 1, 2, 4, 8, 9, 16, 17], dtype=np.int64)
    norm = np.ones_like(reps, dtype=np.uint16)
    ctx = ls.compiler.search_ctx_t(ls.BasisInfo(5), reps, norm, prefix_bits=1)
    alpha = np.array([0, 1, 4, 5], dtype=np.uint64)
    out = ls.compiler.state_to_index(alpha, ctx)
    np.testing.assert_equal(out, [0, 1, 3, -1])

    # Test with larger numbers and more complex patterns
    reps = np.array([0, 3, 7, 15, 31, 63, 127, 255, 511, 1023], dtype=np.int64)
    norm = np.ones_like(reps, dtype=np.uint16)
    ctx = ls.compiler.search_ctx_t(ls.BasisInfo(10), reps, norm, prefix_bits=4)
    alpha = np.array([31, 0, 255, 7, 1024, 15, 63, 512, 3, 127, 8, 511, 32, 16, 256, 64, 1023], dtype=np.uint64)
    out = ls.compiler.state_to_index(alpha, ctx)
    np.testing.assert_equal(out, [4, 0, 7, 2, -1, 3, 5, -1, 1, 6, -1, 8, -1, -1, -1, -1, 9])

@hypothesis.given(
    st.lists(st.integers(min_value=0, max_value=2**20 - 1), min_size=0, max_size=10000, unique=True),
    st.integers(min_value=0, max_value=16),
    st.integers(min_value=1, max_value=10)
)
@hypothesis.example([], 5, 1)
@hypothesis.example([0], 2, 1)
@hypothesis.example([617904, 77273, 561870, 833131, 632964, 27154, 104622, 580730], 0, 3)
@hypothesis.example([48346, 107, 65536, 1], 1, 1)
@hypothesis.settings(max_examples=50, deadline=None, phases=phases)
def test_search(reps, prefix_bits, seed):
    rng = np.random.default_rng(seed=seed)
    reps = np.sort(np.asarray(reps, dtype=np.uint64))
    norm = np.ones(reps.size, dtype=np.uint16)
    if len(reps) == 0 or np.max(reps) == 0: bits = 1
    else: bits = int(1 + np.ceil(np.log2(np.max(reps))))
    search_ctx = ls.compiler.search_ctx_t(ls.BasisInfo(bits), reps, norm, prefix_bits)
    # print(search_ctx.keep_alive)

    def search_ref(arr, needle):
        if len(arr) == 0: return np.full(len(needle), fill_value=-1, dtype=np.int64)
        k = np.clip(np.searchsorted(arr, needle), 0, len(arr) - 1)
        return np.where(arr[k] == needle, k, -1)

    # Search for 0 elements
    out = ls.compiler.state_to_index([], search_ctx)
    assert out.tolist() == []

    # Search for various numbers of elements
    for c in [1, 10, 24, 943]:
        needle = rng.choice(reps, size=c) if len(reps) > 0 else np.full(c, fill_value=123, dtype=np.uint64)
        out = ls.compiler.state_to_index(needle, search_ctx)
        assert out.tolist() == search_ref(reps, needle).tolist()
        needle = rng.integers(low=-5, high=np.max(reps, initial=0) + 5, size=1)
        out = ls.compiler.state_to_index(needle, search_ctx)
        assert out.tolist() == search_ref(reps, needle).tolist()
