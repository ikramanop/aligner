use aligner_core::enums::Protein;
use aligner_core::simple::*;
use aligner_core::{get_blosum62, AlignerTrait, AlignmentTrait};
use aligner_helpers::files::*;
use clap::Parser;
use seq_io::fasta::Reader;
use std::env;
use std::path::Path;

pub mod args;

fn main() {
    let args = args::Args::parse();

    let path = format!(
        "{}/{}",
        env::current_dir().unwrap().to_str().unwrap(),
        args.input
    );

    let contents = load_file_contents(Path::new(&path));
    let bytes = contents.as_bytes();

    let mut reader = Reader::new(bytes);

    let seqs: Vec<_> = match reader.records().collect() {
        Ok(seqs) => seqs,
        Err(err) => panic!("Error with collecting files: {}", err),
    };

    if seqs.len() != 2 {
        panic!("There's should be 2 sequences, not {}", seqs.len())
    }

    let blosum62 = &get_blosum62();
    let query = std::str::from_utf8(&seqs[0].seq).unwrap().to_string();
    let target = std::str::from_utf8(&seqs[1].seq).unwrap().to_string();

    let result;

    if args.global {
        result = SimpleGlobalAligner::<Protein>::from_str_seqs(&query, &target)
            .unwrap()
            .perform_alignment(args.deletions, args.extension, blosum62, None)
            .unwrap();
    } else {
        result = SimpleLocalAligner::<Protein>::from_str_seqs(&query, &target)
            .unwrap()
            .perform_alignment(args.deletions, args.extension, blosum62, None)
            .unwrap();
    }

    println!("{:?}", result.alignment.get_alignment(blosum62))
}
