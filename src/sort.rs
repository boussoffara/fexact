//! Sorting with the C code's failure mode (`isort`).

use crate::error::Error;

/// Sorts `values` ascending.
///
/// The C code's quicksort keeps pending segments on a fixed 10-entry stack
/// and fails when that overflows. Since it always continues with the
/// smaller part, the segment split at depth `m` has at most
/// `len / 2^(m-1)` elements, so overflow needs `len >= 2048`; below that any
/// sort gives the same result. Above it, the original algorithm is run to
/// reproduce the error exactly.
pub(crate) fn sort_ascending(values: &mut [i32]) -> Result<(), Error> {
    if values.len() < 2048 {
        values.sort_unstable();
        Ok(())
    } else {
        quicksort_bounded(values)
    }
}

/// The C code's quicksort (pivot = first element, smaller part first).
fn quicksort_bounded(ix: &mut [i32]) -> Result<(), Error> {
    const DEPTH: usize = 10;
    let mut pending: Vec<(isize, isize)> = Vec::with_capacity(DEPTH);
    // Current segment ix[i..=j].
    let (mut i, mut j): (isize, isize) = (0, ix.len() as isize - 1);
    loop {
        if i >= j {
            match pending.pop() {
                Some(segment) => (i, j) = segment,
                None => return Ok(()),
            }
            continue;
        }

        // Partition around the pivot ix[lo].
        let (lo, hi) = (i, j);
        let pivot = ix[lo as usize];
        j += 1;
        loop {
            loop {
                i += 1;
                if !(i < j && pivot > ix[i as usize]) {
                    break;
                }
            }
            loop {
                j -= 1;
                if ix[j as usize] <= pivot {
                    break;
                }
            }
            if i >= j {
                break;
            }
            ix.swap(i as usize, j as usize);
        }
        ix.swap(lo as usize, j as usize);

        // Continue with the smaller part, saving the larger one.
        if pending.len() == DEPTH {
            return Err(Error::SortStackOverflow);
        }
        if j - lo < hi - j {
            pending.push((j + 1, hi));
            (i, j) = (lo, j - 1);
        } else {
            pending.push((lo, j - 1));
            (i, j) = (j + 1, hi);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pseudo_random(n: usize, seed: u64) -> Vec<i32> {
        let mut state = seed;
        (0..n)
            .map(|_| {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                ((state >> 33) % 50) as i32
            })
            .collect()
    }

    #[test]
    fn quicksort_matches_std_sort() {
        for n in (0..300).chain([1000, 2047]) {
            let mut v = pseudo_random(n, n as u64);
            let mut expected = v.clone();
            expected.sort_unstable();
            quicksort_bounded(&mut v).unwrap();
            assert_eq!(v, expected);
        }
    }

    #[test]
    fn quicksort_stack_overflow() {
        // Equal keys split every segment in half, so 100 000 of them need
        // about 16 nested segments.
        let mut v = vec![7; 100_000];
        assert_eq!(quicksort_bounded(&mut v), Err(Error::SortStackOverflow));
    }
}
