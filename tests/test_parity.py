"""Rust vs the legacy C extension: bit-identical p-values and identical errors
(on tables small enough that the C code's 32-bit counts do not overflow).

Skipped unless the C extension has been built (`python legacy/setup.py
build_ext --inplace`).
"""
import numpy as np
import pytest

import c_reference
from fexact import fexact as fexact_rs

_c = c_reference.load()
if _c is None:
    pytest.skip("C reference not built", allow_module_level=True)
fexact_c = _c.fexact


def _call(f, table, **kw):
    try:
        return ("ok", f(table, **kw))
    except ValueError as e:
        return ("err", str(e))


@pytest.mark.parametrize("seed", range(8))
def test_random_tables(seed):
    rng = np.random.default_rng(seed)
    for _ in range(150):
        r, c = rng.integers(1, 7, size=2)
        cap = max(1, int(rng.integers(1, 80)) // (r * c) + 1)
        t = rng.integers(0, cap + 1, size=(r, c)).astype(np.int64)
        if rng.random() < 0.2:
            t[rng.random(t.shape) < 0.3] = 0
        kw = dict(
            hybrid=bool(rng.random() < 0.4),
            workspace=int(rng.choice([200_000, 20_000, 3000, 700, 0])),
        )
        assert _call(fexact_c, t, **kw) == _call(fexact_rs, t, **kw), (t.tolist(), kw)


def test_f3xact_negative_fact_index():
    # f3xact reads fact[-1] for this table (out of bounds in C).
    t = np.array(
        [[0, 2, 1, 2, 1, 0], [2, 1, 1, 0, 0, 1], [1, 0, 0, 2, 1, 1],
         [2, 0, 0, 0, 2, 2], [1, 0, 0, 0, 0, 0], [0, 1, 2, 1, 2, 1]],
        dtype=np.int64,
    )
    for hybrid in (False, True):
        assert fexact_rs(t, hybrid=hybrid) == fexact_c(t, hybrid=hybrid)
