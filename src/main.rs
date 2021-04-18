pub mod align;
pub mod files;

use crate::align::aligner::*;
use crate::files::*;
use clap::{load_yaml, App};
use seq_io::fasta::Reader;
use std::str::FromStr;

fn main() {
    let config = load_yaml!("../config/opts_config.yml");
    let matches = App::from(config).get_matches();

    let matrix: ndarray::Array2<i32>;
    let deletions: i32;

    if let Some(matrix_option) = matches.value_of("matrix") {
        match matrix_option {
            "PAM250" => matrix = load_pam250(),
            "BLOSUM50" => matrix = load_blosum50(),
            _ => panic!("No such matrix {}", matrix_option),
        }
    } else {
        panic!()
    }

    if let Some(deletions_option) = matches.value_of("deletions") {
        deletions = FromStr::from_str(deletions_option).unwrap()
    } else {
        deletions = 8
    }

    if let Some(input) = matches.value_of("INPUT") {
        let contents = load_file_contents(input);
        let bytes = contents.as_bytes();

        let mut reader = Reader::new(&bytes[..]);

        let seqs: Vec<_> = match reader.records().collect() {
            Ok(seqs) => seqs,
            Err(err) => panic!("Error with collecting files: {}", err),
        };

        if seqs.len() != 2 {
            panic!("There's should be 2 sequences, not {}", seqs.len())
        }

        let mut _aligner: SimpleAligner = SimpleAligner::from_seqs(&seqs[0].seq, &seqs[1].seq);

        let result = _aligner.local_alignment(&deletions, &matrix);

        println!("{:?}", result);
    }
}
