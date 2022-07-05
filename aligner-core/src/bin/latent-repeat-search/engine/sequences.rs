use crate::Error;
use aligner_helpers::csv::Record;
use aligner_helpers::files::load_file_contents;
use seq_io::fasta::Reader;
use std::path::Path;

const N: u8 = 78;

pub(crate) fn import_fasta_file(path: &Path) -> Result<Vec<seq_io::fasta::OwnedRecord>, Error> {
    let contents = load_file_contents(path);
    let bytes = contents.as_bytes();

    let mut reader = Reader::new(bytes);

    let seqs: Vec<_> = match reader.records().collect() {
        Ok(seqs) => seqs,
        Err(err) => {
            return Err(Error {
                msg: err.to_string(),
            })
        }
    };

    if seqs.is_empty() {
        return Err(Error {
            msg: "empty fasta file".to_string(),
        });
    }

    Ok(seqs)
}

pub(crate) fn prepare_sequence(raw_seq: &[u8], data: &[Record]) -> Vec<u8> {
    let mut result = Vec::from(raw_seq);

    for record in data.iter() {
        for elem in result[record.left_coord..record.right_coord].iter_mut() {
            *elem = N
        }
    }

    result
}
