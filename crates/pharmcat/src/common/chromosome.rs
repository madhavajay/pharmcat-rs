//! Chromosome ordering helpers ported from `pgkb-common`.

use std::cmp::Ordering;

/// Compares chromosome names like `pgkb-common`'s `ChromosomeNameComparator`.
pub fn compare_names(name1: Option<&str>, name2: Option<&str>) -> Ordering {
    match (name1, name2) {
        (None, None) => Ordering::Equal,
        (None, Some(_)) => Ordering::Less,
        (Some(_), None) => Ordering::Greater,
        (Some(name1), Some(name2)) => compare_non_null_names(name1, name2),
    }
}

fn compare_non_null_names(name1: &str, name2: &str) -> Ordering {
    let name1 = name1.strip_prefix("chr").unwrap_or(name1);
    let name2 = name2.strip_prefix("chr").unwrap_or(name2);

    match (parse_java_numeric(name1), parse_java_numeric(name2)) {
        (Some(n1), Some(n2)) => n1.cmp(&n2),
        (Some(_), None) => Ordering::Less,
        (None, Some(_)) => Ordering::Greater,
        (None, None) => name1.cmp(name2),
    }
}

fn parse_java_numeric(s: &str) -> Option<i32> {
    if s.is_empty() || !s.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }

    s.parse().ok()
}

/// Compares chromosome positions like `pgkb-common`'s
/// `ChromosomePositionComparator`.
pub fn compare_positions(
    pos1: Option<&str>,
    pos2: Option<&str>,
) -> Result<Ordering, PositionError> {
    match (pos1, pos2) {
        (None, None) => Ok(Ordering::Equal),
        (None, Some(_)) => Ok(Ordering::Less),
        (Some(_), None) => Ok(Ordering::Greater),
        (Some(pos1), Some(pos2)) => {
            let pos1 = parse_position(pos1)?;
            let pos2 = parse_position(pos2)?;

            let chr_order = compare_names(Some(pos1.chromosome), Some(pos2.chromosome));
            if chr_order == Ordering::Equal {
                Ok(pos1.position.cmp(&pos2.position))
            } else {
                Ok(chr_order)
            }
        }
    }
}

#[derive(Debug, Eq, PartialEq)]
struct ChromosomePosition<'a> {
    chromosome: &'a str,
    position: i32,
}

fn parse_position(s: &str) -> Result<ChromosomePosition<'_>, PositionError> {
    let (chromosome, position) = s
        .split_once(':')
        .ok_or_else(|| PositionError::InvalidFormat(s.to_owned()))?;

    let chromosome = chromosome.strip_prefix("chr").unwrap_or(chromosome);
    if chromosome.is_empty()
        || chromosome.len() > 2
        || !chromosome
            .bytes()
            .all(|b| b == b'_' || b.is_ascii_alphanumeric())
    {
        return Err(PositionError::InvalidFormat(s.to_owned()));
    }

    if position.is_empty() || !position.bytes().all(|b| b.is_ascii_digit()) {
        return Err(PositionError::InvalidFormat(s.to_owned()));
    }

    let position = position
        .parse()
        .map_err(|_| PositionError::InvalidFormat(s.to_owned()))?;

    Ok(ChromosomePosition {
        chromosome,
        position,
    })
}

/// Invalid chromosome-position input.
#[derive(Debug, Eq, PartialEq)]
pub enum PositionError {
    /// Input is not in the expected `chrX:1234` shape.
    InvalidFormat(String),
}

impl std::fmt::Display for PositionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidFormat(value) => {
                write!(
                    f,
                    "'{value}' is not in the expected chromosomal position format"
                )
            }
        }
    }
}

impl std::error::Error for PositionError {}

#[cfg(test)]
mod tests {
    use super::{PositionError, compare_names, compare_positions};
    use std::cmp::Ordering;

    fn sign(ordering: Ordering) -> i8 {
        match ordering {
            Ordering::Less => -1,
            Ordering::Equal => 0,
            Ordering::Greater => 1,
        }
    }

    #[test]
    fn chromosome_name_comparator_matches_pgkb_common_tests() {
        assert_eq!(sign(compare_names(None, None)), 0);
        assert_eq!(sign(compare_names(None, Some("chr1"))), -1);
        assert_eq!(sign(compare_names(Some("chr1"), None)), 1);

        assert_eq!(sign(compare_names(Some("chr1"), Some("chr1"))), 0);
        assert_eq!(sign(compare_names(Some("chr1"), Some("chr3"))), -1);
        assert_eq!(sign(compare_names(Some("chr13"), Some("chr1"))), 1);
        assert_eq!(sign(compare_names(Some("chr13"), Some("chr21"))), -1);
        assert_eq!(sign(compare_names(Some("chr13"), Some("chrX"))), -1);
        assert_eq!(sign(compare_names(Some("chr20"), Some("chrX"))), -1);
        assert_eq!(sign(compare_names(Some("chrX"), Some("chrY"))), -1);
        assert_eq!(sign(compare_names(Some("chrX"), Some("chrX"))), 0);

        assert_eq!(sign(compare_names(Some("1"), Some("1"))), 0);
        assert_eq!(sign(compare_names(Some("1"), Some("3"))), -1);
        assert_eq!(sign(compare_names(Some("13"), Some("1"))), 1);
        assert_eq!(sign(compare_names(Some("13"), Some("21"))), -1);
        assert_eq!(sign(compare_names(Some("13"), Some("X"))), -1);
        assert_eq!(sign(compare_names(Some("20"), Some("X"))), -1);
        assert_eq!(sign(compare_names(Some("X"), Some("20"))), 1);
        assert_eq!(sign(compare_names(Some("X"), Some("Y"))), -1);
        assert_eq!(sign(compare_names(Some("X"), Some("X"))), 0);
    }

    #[test]
    fn chromosome_position_comparator_matches_pgkb_common_tests() {
        assert_eq!(sign(compare_positions(None, None).unwrap()), 0);
        assert_eq!(sign(compare_positions(None, Some("chr1:1")).unwrap()), -1);
        assert_eq!(sign(compare_positions(Some("chr1:1"), None).unwrap()), 1);

        assert_eq!(
            sign(compare_positions(Some("chr1:4"), Some("chr1:100")).unwrap()),
            -1
        );
        assert_eq!(
            sign(compare_positions(Some("chr1:4"), Some("chr1:4")).unwrap()),
            0
        );
        assert_eq!(
            sign(compare_positions(Some("chr3:4"), Some("chr1:4")).unwrap()),
            1
        );
        assert_eq!(
            sign(compare_positions(Some("chr4:100"), Some("chr1:400")).unwrap()),
            1
        );
    }

    #[test]
    fn chromosome_position_comparator_rejects_bad_input_like_pgkb_common() {
        assert_eq!(
            compare_positions(Some("chr1"), Some("chr1:100")),
            Err(PositionError::InvalidFormat("chr1".to_owned()))
        );
        assert_eq!(
            compare_positions(Some("chr1:1"), Some(":100")),
            Err(PositionError::InvalidFormat(":100".to_owned()))
        );
    }
}
