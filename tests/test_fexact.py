import numpy as np
import pytest

from fexact import fexact


def test_fexact_small():
    k = np.array([[5, 1, 5], [1, 1, 5]], dtype=np.int64)
    assert np.abs(0.4004524886877815 - fexact(k)) < 10e-5


def test_fexact_small_ws():
    k = np.array([[5, 1, 5], [1, 1, 5]], dtype=np.int64)
    assert np.abs(0.4004524886877815 - fexact(k, workspace=10000000)) < 10e-5


def test_fexact_hybrid():
    k = np.array([[10, 6, 7], [5, 7, 9]], dtype=np.int64)
    assert np.abs(0.3962 - fexact(k, workspace=10000000, hybrid=True)) < 10e-3
    assert np.abs(0.4158 - fexact(k, workspace=10000000, hybrid=False)) < 10e-3


def test_errors_neg():
    k = np.array([[-1, 6, 7] * 10, [5, 7, 9] * 10], dtype=np.int64)
    with pytest.raises(ValueError, match="nonnegative"):
        fexact(k, workspace=10000000, hybrid=False)


def test_errors_ws():
    k = np.array([[100000, 6, 7], [100000, 7, 9]], dtype=np.int64)
    with pytest.raises(ValueError, match="Out of workspace"):
        fexact(k, hybrid=False)


def test_wide_margins():
    # The C code's 32-bit node keys overflow here (FEXACT error 501).
    k = np.array([[100000, 6, 7], [100000, 7, 9]], dtype=np.int64)
    assert abs(fexact(k, workspace=10000000) - 0.917722) < 1e-6


def test_many_columns():
    # The C code's 32-bit path counts overflow here (it gives 0.329479).
    k = np.array([[10, 6, 7] * 10, [5, 7, 9] * 10], dtype=np.int64)
    p = fexact(k, workspace=20_000_000)
    assert abs(p - 0.928552) < 1e-6, p


def test_errors_ws2():
    k = np.array([[10, 6, 7] * 10, [5, 7, 9] * 10], dtype=np.int64)
    with pytest.raises(ValueError, match="Out of workspace"):
        fexact(k, workspace=2**3, hybrid=False)


def test_errors_ldstep():
    k = np.array(
        [
            [40, 0, 0, 0, 0],
            [0, 40, 10, 0, 10],
            [0, 0, 10, 0, 0],
            [10, 0, 0, 30, 0],
            [0, 0, 0, 20, 30],
        ],
        dtype=np.int64,
    )
    with pytest.raises(ValueError):
        fexact(k, workspace=int(3e7))


def test_input_validation():
    with pytest.raises(TypeError):
        fexact([[1, 2], [3, 4]])
    with pytest.raises(ValueError):
        fexact(np.array([[1, 2], [3, 4]], dtype=np.int32))
    with pytest.raises(ValueError):
        fexact(np.array([1, 2, 3], dtype=np.int64))
    with pytest.raises(ValueError, match="zero"):
        fexact(np.zeros((2, 2), dtype=np.int64))
    with pytest.raises(ValueError, match="NCOL"):
        fexact(np.array([[1], [2]], dtype=np.int64))


def test_strided_input():
    k = np.array([[5, 1, 5], [1, 1, 5]], dtype=np.int64)
    expected = fexact(k)
    assert fexact(np.asfortranarray(k)) == expected
    assert fexact(k.T.copy().T) == expected
    wide = np.zeros((2, 6), dtype=np.int64)
    wide[:, ::2] = k
    assert fexact(wide[:, ::2]) == expected


def test_many_matches_single():
    from fexact import fexact_many

    rng = np.random.default_rng(0)
    tables = [rng.integers(0, 6, size=(rng.integers(2, 4), rng.integers(2, 5))).astype(np.int64)
              for _ in range(60)]
    tables = [t for t in tables if t.sum() > 0 and t.shape[1] > 1]
    expected = [fexact(t, workspace=10**6, hybrid=True) for t in tables]
    assert fexact_many(tables, workspace=10**6, hybrid=True) == expected
    assert fexact_many(tables, workspace=10**6, hybrid=True, threads=1) == expected
    assert fexact_many([], threads=2) == []


def test_many_errors_name_the_table():
    from fexact import fexact_many

    good = np.array([[1, 2], [3, 4]], dtype=np.int64)
    bad = np.array([[-1, 2], [3, 4]], dtype=np.int64)
    with pytest.raises(ValueError, match="table 1: All elements of TABLE must be nonnegative"):
        fexact_many([good, bad, good])


def test_simulate_close_to_exact():
    from fexact import fexact_simulate

    k = np.array([[10, 6, 7], [5, 7, 9]], dtype=np.int64)
    p = fexact(k)
    sim = fexact_simulate(k, samples=200_000, seed=1)
    assert abs(sim - p) < 5 * (p * (1 - p) / 200_000) ** 0.5


def test_simulate_seed_and_threads():
    from fexact import fexact_simulate

    k = np.array([[10, 6, 7], [5, 7, 9], [3, 8, 2]], dtype=np.int64)
    a = fexact_simulate(k, samples=20_000, seed=42, threads=1)
    assert a == fexact_simulate(k, samples=20_000, seed=42, threads=7)
    assert a == fexact_simulate(k, samples=20_000, seed=42)
    assert 0 < fexact_simulate(k, samples=100) <= 1


def test_simulate_errors():
    from fexact import fexact_simulate

    with pytest.raises(ValueError, match="nonnegative"):
        fexact_simulate(np.array([[-1, 2], [3, 4]], dtype=np.int64))
    with pytest.raises(ValueError):
        fexact_simulate(np.array([1, 2], dtype=np.int64))


def test_values_beyond_int32():
    from fexact import fexact_simulate

    # The C extension silently wrapped these to 32 bits (2**32 + 5 -> 5).
    k = np.array([[2**32 + 5, 1], [1, 5]], dtype=np.int64)
    with pytest.raises(ValueError, match="Out of workspace"):
        fexact(k)
    assert fexact_simulate(k, samples=2000, seed=1) < 0.01
    k = np.array([[2**31, 0], [1, 1]], dtype=np.int64)
    with pytest.raises(ValueError, match="Out of workspace"):
        fexact(k)
    # P(x22 = 1) = 2 / (2**31 + 2): no random table is as extreme.
    assert fexact_simulate(k, samples=2000, seed=1) == 1 / 2001
