use aligner::align::aligner_core::{Aligner, AlignmentResult};
use aligner::align::enums::Protein;
use aligner::align::simple_aligner::*;
use aligner::files::*;
use clap::{load_yaml, App};
use seq_io::fasta::Reader;
use std::env;
use std::path::Path;
use std::str::{from_utf8, FromStr};

fn main() {
    let config = load_yaml!("../../config/opts_config.yml");
    let matches = App::from(config).get_matches();

    let matrix: ndarray::Array2<f64>;
    let deletions: f64;

    let db = match get_db("database/matrices") {
        Ok(db) => db,
        Err(err) => panic!("{}", err),
    };

    if let Some(matrix_option) = matches.value_of("matrix") {
        match matrix_option {
            "PAM250" => {
                matrix = match db.get("PAM250") {
                    Ok(matrix) => convert_csv_to_matrix(&matrix.unwrap().to_vec(), (25, 25)),
                    _ => {
                        panic!()
                    }
                }
            }
            "BLOSUM50" => {
                matrix = match db.get("BLOSUM50") {
                    Ok(matrix) => convert_csv_to_matrix(&matrix.unwrap().to_vec(), (25, 25)),
                    _ => {
                        panic!()
                    }
                }
            }
            _ => panic!("No such matrix {}", matrix_option),
        }
    } else {
        panic!()
    }

    match db.flush() {
        Ok(_) => {}
        Err(err) => {
            panic!("{}", err)
        }
    };

    if let Some(deletions_option) = matches.value_of("deletions") {
        deletions = FromStr::from_str(deletions_option).unwrap()
    } else {
        deletions = 8f64
    }

    if let Some(input) = matches.value_of("INPUT") {
        let base_dir = match env::var("BASE_SEQS_PATH") {
            Ok(var) => var,
            Err(_) => String::from_str("").unwrap(),
        };

        let path = Path::new(&base_dir).join(input);

        let contents = load_file_contents(path.as_path());
        let bytes = contents.as_bytes();

        let mut reader = Reader::new(bytes);

        let seqs: Vec<_> = match reader.records().collect() {
            Ok(seqs) => seqs,
            Err(err) => panic!("Error with collecting files: {}", err),
        };

        if seqs.len() != 2 {
            panic!("There's should be 2 sequences, not {}", seqs.len())
        }

        let mut aligner: SimpleAligner = SimpleAligner::from_seqs(&seqs[0].seq, &seqs[1].seq);

        let alignment: Box<dyn AlignmentResult>;

        if matches.is_present("global") {
            alignment = aligner.global_alignment(&deletions, &matrix);
        } else {
            alignment = aligner.local_alignment(&deletions, &matrix);
        }

        alignment.represent();

        if let Some(output_path) = matches.value_of("output") {
            let (alignment_1_u8, alignment_2_u8) = (
                Protein::protein_vec_to_u8_vec(alignment.get_optimal_alignment().0),
                Protein::protein_vec_to_u8_vec(alignment.get_optimal_alignment().1),
            );

            let (alignment_1, alignment_2) = (
                from_utf8(&alignment_1_u8).unwrap(),
                from_utf8(&alignment_2_u8).unwrap(),
            );

            let output = format!("{}\n{}", alignment_1, alignment_2);

            write_to_file(Path::new(output_path), &output);
        }
    }
}
