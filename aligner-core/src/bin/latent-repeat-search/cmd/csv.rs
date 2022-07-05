use crate::cmd::{CMDOptions, CMDResult, DIRECT, INVERSE};
use crate::engine::calc::perform_calculation_per_sequence;
use crate::engine::sequences::{import_fasta_file, prepare_sequence};
use crate::error::Error;
use aligner_helpers::csv::read_csv;
use std::collections::HashMap;

pub(crate) fn run_csv_cmd(opts: &CMDOptions) -> CMDResult {
    info!("Entering exploring mode with csv support!!");

    let seqs = match import_fasta_file(opts.fasta_path.as_ref().unwrap().as_path()) {
        Ok(seqs) => seqs,
        Err(err) => return Err(err),
    };

    let sequence_data = match read_csv(opts.csv_path.as_ref().unwrap().as_path()) {
        Ok(data) => data,
        Err(err) => {
            return Err(Error {
                msg: err.to_string(),
            })
        }
    };

    let mut result = HashMap::new();

    for raw_seq in seqs.iter() {
        let head = std::str::from_utf8(&raw_seq.head).unwrap();

        let prepared_seq = match sequence_data.get(head) {
            Some(data) => prepare_sequence(&raw_seq.seq, data),
            None => raw_seq.seq.clone(),
        };

        let mut sequence_result = perform_calculation_per_sequence(opts, &prepared_seq, head);

        match sequence_result.remove(DIRECT) {
            Some(data) => result.insert(head.to_string(), data),
            None => None,
        };

        match sequence_result.remove(INVERSE) {
            Some(data) => result.insert(format!("{}-reversed", head).to_string(), data),
            None => None,
        };
    }

    Ok(result)
}
