use crate::cmd::CMDOptions;
use crate::engine::task::Task;
use crate::engine::{index_coord, rotate_indices};
use crate::filter;
use aligner_core::enums::{BioData, Index, DNA};
use aligner_core::pwm::PWMAligner;
use aligner_core::AlignmentTrait;
use aligner_core::{get_random_pwm, AlignerTrait};
use aligner_helpers::matrices::transform_matrix;
use ndarray::Array2;
use rand::prelude::SliceRandom;
use rand::thread_rng;
use std::collections::HashMap;
use std::sync::mpsc::channel;
use std::thread;

const Z: f64 = 3f64;

pub(crate) fn calculate_starting_values(
    query: &[DNA],
    matrix: &Array2<f64>,
    opts: &CMDOptions,
) -> (f64, f64) {
    let length = query.len();
    let query_offset = opts.query_offset;
    let threads = opts.threads;
    let repeat_length = opts.repeat_length;
    let deletions = opts.deletions;
    let extension = opts.extension;

    let mut shuffled_query = Vec::from(query);
    let mut thread_rng = thread_rng();
    shuffled_query.shuffle(&mut thread_rng);

    let mut fs = vec![];

    let step = if opts.simple_init {
        query.len() / 1000
    } else {
        opts.query_offset
    };

    info!("Calculating starting values");

    let (tx, rx) = channel();

    for i in 0..opts.threads {
        let tx = tx.clone();
        let shuffled_query = shuffled_query.clone();
        let matrix = matrix.clone();

        thread::spawn(move || {
            for j in (i * query_offset..length).step_by(step * threads) {
                let border = if j + repeat_length + query_offset >= length {
                    length
                } else {
                    j + repeat_length + query_offset
                };

                let result = PWMAligner::<DNA>::from_seqs(&shuffled_query[j..border], &[])
                    .unwrap()
                    .perform_alignment(deletions, extension, &matrix, None)
                    .unwrap();

                tx.send(result.alignment.f).unwrap();
            }
            drop(tx);
        });
    }
    drop(tx);

    while let Ok(f) = rx.recv() {
        fs.push(f);
    }

    let mean = fs.iter().sum::<f64>() / fs.len() as f64;

    (
        mean,
        (fs.iter()
            .map(|value| (value - mean) * (value - mean))
            .sum::<f64>()
            / fs.len() as f64)
            .sqrt(),
    )
}

pub(crate) fn calculate_cycle(
    query: &[DNA],
    matrix: &Array2<f64>,
    indices: &[Index],
    mean: f64,
    std: f64,
    opts: &CMDOptions,
) -> Vec<Task> {
    let length = query.len();
    let query_offset = opts.query_offset;
    let threads = opts.threads;
    let repeat_length = opts.repeat_length;
    let deletions = opts.deletions;
    let extension = opts.extension;

    let mut tasks = vec![];

    let (tx, rx) = channel();

    for i in 0..opts.threads {
        let tx = tx.clone();
        let query = Vec::from(query);
        let matrix = matrix.clone();
        let indices = Vec::from(indices);

        thread::spawn(move || {
            for j in (i * query_offset..length).step_by(query_offset * threads) {
                let border = if j + repeat_length + query_offset >= length {
                    length
                } else {
                    j + repeat_length + query_offset
                };

                let result = PWMAligner::<DNA>::from_seqs(&query[j..border], &[])
                    .unwrap()
                    .perform_alignment(deletions, extension, &matrix, None)
                    .unwrap();

                tx.send(Task {
                    alignment: result.alignment.clone(),
                    left_coord: index_coord(j, &indices),
                    right_coord: index_coord(border, &indices),
                    z: (result.alignment.f - mean) / std,
                })
                .unwrap();
            }
            drop(tx);
        });
    }
    drop(tx);

    while let Ok(task) = rx.recv() {
        if task.z >= Z {
            debug!("{:?}", task);
            tasks.push(task);
        }
    }

    tasks
}

pub(crate) fn perform_calculation_per_sequence(
    opts: &CMDOptions,
    raw_seq: &[u8],
    head: &str,
) -> HashMap<String, (Vec<Task>, Array2<f64>)> {
    let (mut query, frequencies, indices) =
        DNA::from_u8_vec_with_freqs_and_indices(raw_seq).unwrap();
    let mut matrix = get_random_pwm(opts.repeat_length);

    matrix = transform_matrix(
        &matrix,
        0f64,
        opts.deletions * opts.extension as f64,
        &frequencies,
    )
    .unwrap();

    info!(
        "Successfully obtained a chromosome \"{}\" with length={}",
        head,
        query.len()
    );

    info!("Calculating direct of {}", head);

    let (mut mean, mut std) = calculate_starting_values(&query, &matrix, opts);

    info!("Calculated starting mean={} and std={}", mean, std);

    let mut result = HashMap::new();

    let mut tasks = vec![];

    for i in 0..opts.repeats {
        info!("Calculating cycle {}", i + 1);

        let new_tasks = calculate_cycle(&query, &matrix, &indices, mean, std, opts);

        if new_tasks.is_empty() {
            break;
        }
        tasks = filter(&new_tasks).unwrap();

        if i < opts.repeats - 1 {
            mean = tasks.iter().map(|value| value.alignment.f).sum::<f64>() / tasks.len() as f64;
            std = tasks
                .iter()
                .map(|value| (value.alignment.f - mean) * (value.alignment.f - mean))
                .sum::<f64>()
                / tasks.len() as f64;

            info!("mean={} and std={} for this cycle", mean, std);

            matrix = Array2::<f64>::zeros((matrix.shape()[0], matrix.shape()[1]));
            for task in tasks.iter() {
                matrix = matrix + task.alignment.get_frequency_matrix();
            }

            matrix = transform_matrix(
                &matrix,
                0f64,
                opts.deletions * opts.extension as f64,
                &frequencies,
            )
            .unwrap();
        }

        info!("Calculated cycle {}", i + 1);
    }

    result.insert("direct".to_string(), (tasks, matrix.clone()));

    if opts.reverse {
        info!("Calculating inverse of {}", head);

        query.reverse();

        let rotated_indices = rotate_indices(&indices, query.len());

        let mut tasks_inverted =
            calculate_cycle(&query, &matrix, &rotated_indices, mean, std, opts);

        tasks_inverted = filter(&tasks_inverted).unwrap();

        result.insert("inverse".to_string(), (tasks_inverted, matrix));
    }

    result
}
