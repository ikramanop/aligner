use crate::cmd::{CMDOptions, CMDResult};
use crate::engine::calc::{calculate_cycle, calculate_starting_values};
use crate::engine::generate_descendants;
use crate::engine::MutationPercent::Quarter;
use aligner_core::enums::{BioData, Index, DNA};
use aligner_core::heuristic::HeuristicPWMAligner;
use aligner_core::{get_random_pwm, AlignerTrait, Heuristics};
use std::collections::HashMap;

const TEST_SEQUENCE_LENGTH: usize = 100000;
const DESCENDANTS_AMOUNT: usize = 10;

pub(crate) fn run_testing_cmd(opts: &CMDOptions) -> CMDResult {
    info!("Entering testing mode!!");

    // testing chromosome
    let sequence_raw = DNA::random_seq(TEST_SEQUENCE_LENGTH).unwrap();

    // sample sequence to find
    let (query, freqs) =
        DNA::random_seq_with_freqs(opts.repeat_length + opts.query_offset).unwrap();

    // find matrix for randomly generated matrix
    debug!("Finding matrix for sample sequence");

    let mut matrix = get_random_pwm(opts.repeat_length);

    let mut aligner = HeuristicPWMAligner::from_seqs(&query, &[]).unwrap();
    let result = aligner
        .perform_alignment(
            opts.deletions,
            opts.extension,
            &matrix,
            Option::Some(Heuristics {
                r_squared: opts.rsquared,
                kd: opts.kd,
                frequencies: freqs,
            }),
        )
        .unwrap();

    matrix = result.matrix.unwrap();

    // get descendants of sample sequence
    debug!("Getting descendants of sample sequence");

    let mut descendants = generate_descendants(&query, DESCENDANTS_AMOUNT, Quarter).unwrap();

    let offset = sequence_raw.len() / (descendants.len() + 1);

    let mut sequence = sequence_raw[..offset].to_vec();
    for (i, descendant) in descendants.iter_mut().enumerate() {
        sequence.append(descendant);
        sequence.append(&mut sequence_raw[offset * i..offset * (i + 1)].to_vec())
    }

    // try to find descendants in random sequence
    debug!("Searching for descendants");

    let (mean, std) = calculate_starting_values(&sequence, &matrix, opts);

    debug!("Calculating exactly one cycle");

    let mut result = HashMap::new();

    result.insert(
        String::from("test"),
        (
            calculate_cycle(&sequence, &matrix, &Vec::<Index>::new(), mean, std, opts),
            matrix,
        ),
    );

    Ok(result)
}
