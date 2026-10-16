"""Regenerate tests/data/golden.txt from the legacy C extension.

Usage: python tests/gen_golden.py  (after `python legacy/setup.py build_ext --inplace`)

Line format: `nrow ncol hybrid workspace|cells (row-major)|ok <p>` or
`...|err <message>` with newlines in messages escaped as `\\n`. p-values are
Python reprs, which round-trip exactly.
"""
import multiprocessing as mp
import pathlib
import time

import numpy as np

import c_reference

fexact = c_reference.load().fexact

MAX_SECONDS = 0.004  # keep the corpus fast to replay
WORKSPACES = [200_000, 2_000_000, 20_000, 3000, 700, 0, 1, 330, 3600]
WEIGHTS = [0.35, 0.25, 0.15, 0.08, 0.05, 0.03, 0.03, 0.03, 0.03]


def gen(seed, n, big):
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        if big:
            r, c = rng.integers(2, 9, size=2)
            cap = max(1, int(rng.integers(1, 90)) // (r * c) + 1)
        else:
            r, c = rng.integers(1, 7, size=2)
            cap = max(1, int(rng.integers(1, 120)) // (r * c) + 1)
        t = rng.integers(0, cap + 1, size=(r, c)).astype(np.int64)
        if rng.random() < 0.2:
            t[rng.random(t.shape) < 0.3] = 0
        hybrid = bool(rng.random() < 0.4)
        ws = int(rng.choice(WORKSPACES, p=WEIGHTS))
        t0 = time.perf_counter()
        try:
            res = f"ok {fexact(t, workspace=ws, hybrid=hybrid)!r}"
        except ValueError as e:
            res = "err " + str(e).replace("\n", "\\n")
        if time.perf_counter() - t0 > MAX_SECONDS:
            continue
        cells = " ".join(map(str, t.ravel()))
        out.append(f"{r} {c} {int(hybrid)} {ws}|{cells}|{res}")
    return out


if __name__ == "__main__":
    jobs = [(s, 1500, False) for s in range(100, 108)] + [(s, 600, True) for s in range(200, 208)]
    with mp.Pool() as pool:
        lines = [l for chunk in pool.starmap(gen, jobs) for l in chunk]
    path = pathlib.Path(__file__).parent / "data" / "golden.txt"
    path.write_text("\n".join(lines) + "\n")
    print(f"{len(lines)} cases -> {path}")
