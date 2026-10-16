//! Replays tests/data/golden.txt, recorded from the legacy C extension
//! (see tests/gen_golden.py): p-values must match bit for bit and errors
//! must carry the same message.

use fexact::{fisher_exact, Hybrid, Options, Table};

#[test]
fn matches_c_extension() {
    let lines: Vec<&str> = include_str!("data/golden.txt")
        .lines()
        .filter(|l| !l.is_empty())
        .collect();
    let cases = lines.len();
    let threads = std::thread::available_parallelism().map_or(4, |n| n.get());
    let failures: Vec<String> = std::thread::scope(|s| {
        let handles: Vec<_> = lines
            .chunks(cases.div_ceil(threads))
            .map(|chunk| s.spawn(move || chunk.iter().filter_map(|l| check(l)).collect::<Vec<_>>()))
            .collect();
        handles
            .into_iter()
            .flat_map(|h| h.join().unwrap())
            .collect()
    });
    assert!(cases > 10_000, "corpus too small: {cases}");
    assert!(
        failures.is_empty(),
        "{} of {cases} cases differ from C, e.g.:\n{}",
        failures.len(),
        failures[..failures.len().min(5)].join("\n")
    );
}

/// Runs one corpus line; returns a description if it differs from C.
fn check(line: &str) -> Option<String> {
    let mut parts = line.splitn(3, '|');
    let (header, cells, expected) = (parts.next()?, parts.next()?, parts.next()?);
    let h: Vec<i32> = header.split(' ').map(|x| x.parse().unwrap()).collect();
    let (nrow, ncol, hybrid, workspace) = (h[0], h[1], h[2] == 1, h[3]);
    let cells: Vec<i64> = cells.split(' ').map(|x| x.parse().unwrap()).collect();
    let table = Table::from_row_major(&cells, nrow as usize, ncol as usize);
    let options = Options {
        workspace: workspace as u32,
        hybrid: hybrid.then_some(Hybrid::COCHRAN),
        ..Options::default()
    };
    let got = match fisher_exact(&table, &options) {
        Ok(r) => format!("ok {:?}", r.p_value),
        Err(e) => format!("err {}", e.to_string().replace('\n', "\\n")),
    };
    let same = match (expected.strip_prefix("ok "), got.strip_prefix("ok ")) {
        (Some(a), Some(b)) => {
            a.parse::<f64>().unwrap().to_bits() == b.parse::<f64>().unwrap().to_bits()
        }
        _ => expected == got,
    };
    (!same).then(|| format!("{line}\n    got: {got}"))
}
