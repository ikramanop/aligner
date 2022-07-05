use crate::cmd::{CMDOptions, CMDResult, DIRECT, INVERSE};
use crate::engine::calc::perform_calculation_per_sequence;
use crate::engine::sequences::import_fasta_file;
use std::collections::HashMap;

pub(crate) fn run_exploring_cmd(opts: &CMDOptions) -> CMDResult {
    info!("Entering exploring mode!!");

    let seqs = match import_fasta_file(opts.fasta_path.as_ref().unwrap().as_path()) {
        Ok(seqs) => seqs,
        Err(err) => return Err(err),
    };

    let mut result = HashMap::new();

    for raw_seq in seqs.iter() {
        let head = std::str::from_utf8(&raw_seq.head).unwrap();

        let mut sequence_result = perform_calculation_per_sequence(opts, &raw_seq.seq, head);

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
