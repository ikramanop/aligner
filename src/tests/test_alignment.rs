use crate::align::aligner_core::*;
use crate::align::enums::Direction::{Beginning, Diagonal, Left, Top};
use crate::align::enums::Protein;
use crate::files::*;
use ndarray::array;
use seq_io::fasta::Reader;

#[test]
fn test_global_alignment() {
    let matrix = load_blosum50();

    let expected_result = Box::new(GlobalAlignmentResult {
        alignment_matrix: array![
            [0, -8, -16, -24, -32, -40, -48, -56, -64, -72, -88],
            [-8, -2, -9, -17, -25, -33, -41, -49, -57, -65, -73],
            [-16, -10, -3, -4, -12, -20, -28, -36, -44, -52, -60],
            [-24, -18, -11, -6, -7, -15, -5, -13, -21, -29, -37],
            [-32, -14, -18, -13, -8, -9, -13, -7, -3, -11, -19],
            [-40, -22, -8, -16, -16, -9, -12, -15, -7, 3, -5],
            [-48, -30, -16, -3, -11, -11, -12, -12, -15, -5, 2],
            [-64, -38, -24, -11, -6, -12, -14, -15, -12, -9, 1]
        ],
        direction_matrix: array![
            [Beginning, Left, Left, Left, Left, Left, Left, Left, Left, Left, Left],
            [Top, Diagonal, Diagonal, Left, Left, Left, Left, Left, Left, Left, Left],
            [Top, Top, Diagonal, Diagonal, Left, Left, Left, Left, Left, Left, Left],
            [Top, Top, Top, Diagonal, Diagonal, Left, Diagonal, Left, Left, Left, Left],
            [
                Top, Diagonal, Diagonal, Diagonal, Diagonal, Diagonal, Top, Diagonal, Diagonal,
                Left, Left
            ],
            [Top, Top, Diagonal, Left, Top, Diagonal, Diagonal, Top, Diagonal, Diagonal, Left],
            [Top, Top, Top, Diagonal, Left, Diagonal, Diagonal, Diagonal, Top, Top, Diagonal],
            [
                Top, Top, Top, Top, Diagonal, Diagonal, Diagonal, Diagonal, Diagonal, Diagonal,
                Diagonal
            ]
        ],
        optimal_alignment: (
            vec![
                Protein::H,
                Protein::E,
                Protein::A,
                Protein::G,
                Protein::A,
                Protein::W,
                Protein::G,
                Protein::H,
                Protein::E,
                Protein::Blank,
                Protein::E,
            ],
            vec![
                Protein::Blank,
                Protein::P,
                Protein::A,
                Protein::Blank,
                Protein::Blank,
                Protein::W,
                Protein::Blank,
                Protein::H,
                Protein::E,
                Protein::A,
                Protein::E,
            ],
        ),
    });

    let contents = load_file_contents("./examples/book_example_1.fasta");
    let bytes = contents.as_bytes();

    let mut reader = Reader::new(&bytes[..]);

    let seqs: Vec<_> = match reader.records().collect() {
        Ok(seqs) => seqs,
        Err(_) => panic!("Test fails"),
    };

    let alignment =
        SimpleAligner::from_seqs(&seqs[0].seq, &seqs[1].seq).global_alignment(&8, &matrix);

    assert_eq!(
        alignment.get_alignment_matrix(),
        expected_result.alignment_matrix
    );
    assert_eq!(
        alignment.get_direction_matrix(),
        expected_result.direction_matrix
    );
    assert_eq!(
        alignment.get_optimal_alignment().0[..],
        expected_result.optimal_alignment.0[..]
    );
    assert_eq!(
        alignment.get_optimal_alignment().1[..],
        expected_result.optimal_alignment.1[..]
    );
}

#[test]
fn test_local_alignment() {
    let matrix = load_blosum50();

    let expected_result = Box::new(GlobalAlignmentResult {
        alignment_matrix: array![
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 5, 0, 5, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 2, 0, 20, 12, 4, 0, 0],
            [0, 10, 2, 0, 0, 0, 12, 18, 22, 14, 6],
            [0, 2, 16, 8, 0, 0, 4, 10, 18, 28, 20],
            [0, 0, 8, 21, 13, 5, 0, 4, 10, 20, 27],
            [0, 0, 6, 13, 18, 12, 4, 0, 4, 16, 26]
        ],
        direction_matrix: array![
            [
                Beginning, Beginning, Beginning, Beginning, Beginning, Beginning, Beginning,
                Beginning, Beginning, Beginning, Beginning
            ],
            [
                Beginning, Beginning, Beginning, Beginning, Beginning, Beginning, Beginning,
                Beginning, Beginning, Beginning, Beginning
            ],
            [
                Beginning, Beginning, Beginning, Diagonal, Beginning, Diagonal, Beginning,
                Beginning, Beginning, Beginning, Beginning
            ],
            [
                Beginning, Beginning, Beginning, Beginning, Diagonal, Beginning, Diagonal, Left,
                Left, Beginning, Beginning
            ],
            [
                Beginning, Diagonal, Left, Beginning, Beginning, Beginning, Top, Diagonal,
                Diagonal, Left, Left
            ],
            [
                Beginning, Top, Diagonal, Left, Beginning, Beginning, Top, Top, Diagonal, Diagonal,
                Left
            ],
            [
                Beginning, Beginning, Top, Diagonal, Left, Left, Beginning, Diagonal, Top, Top,
                Diagonal
            ],
            [
                Beginning, Beginning, Diagonal, Top, Diagonal, Diagonal, Left, Beginning, Diagonal,
                Diagonal, Diagonal
            ]
        ],
        optimal_alignment: (
            vec![Protein::A, Protein::W, Protein::G, Protein::H, Protein::E],
            vec![
                Protein::A,
                Protein::W,
                Protein::Blank,
                Protein::H,
                Protein::E,
            ],
        ),
    });

    let contents = load_file_contents("./examples/book_example_1.fasta");
    let bytes = contents.as_bytes();

    let mut reader = Reader::new(&bytes[..]);

    let seqs: Vec<_> = match reader.records().collect() {
        Ok(seqs) => seqs,
        Err(_) => panic!("Test fails"),
    };

    let alignment =
        SimpleAligner::from_seqs(&seqs[0].seq, &seqs[1].seq).local_alignment(&8, &matrix);

    assert_eq!(
        alignment.get_alignment_matrix(),
        expected_result.alignment_matrix
    );
    assert_eq!(
        alignment.get_direction_matrix(),
        expected_result.direction_matrix
    );
    assert_eq!(
        alignment.get_optimal_alignment().0[..],
        expected_result.optimal_alignment.0[..]
    );
    assert_eq!(
        alignment.get_optimal_alignment().1[..],
        expected_result.optimal_alignment.1[..]
    );
}
