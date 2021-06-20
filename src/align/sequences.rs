use crate::align::enums::Protein;
use ndarray::{arr1, Array1};
use ndarray_rand::rand::{thread_rng, Rng};
use ndarray_rand::rand_distr::Uniform;
use std::collections::HashMap;

const UTF8_SYMBOLS: &'static [u8] = &[
    65, 82, 78, 68, 67, 81, 69, 71, 72, 73, 76, 75, 77, 70, 80, 83, 84, 87, 89, 86,
];

pub fn generate_test_sequences_pair(
    avg_length: &i32,
    error: &i32,
) -> ((Vec<Protein>, Vec<Protein>), Array1<f64>) {
    let (mut sequence_1, mut sequence_2) = (Vec::<Protein>::new(), Vec::<Protein>::new());
    let mut counters = HashMap::<Protein, i32>::new();
    let mut rng = thread_rng();
    let u8_distr = Uniform::new(0, UTF8_SYMBOLS.len());
    let error_distr = Uniform::new_inclusive(-error, error);

    let seq_1_length = *avg_length + rng.sample(error_distr);
    let seq_2_length = *avg_length + rng.sample(error_distr);
    let total_length = seq_1_length + seq_2_length;

    for _ in 0..seq_1_length {
        let symbol = Protein::match_with_u8(UTF8_SYMBOLS[rng.sample(u8_distr)]);
        sequence_1.push(symbol);
        let counter = counters.entry(symbol).or_insert(0);
        *counter += 1;
    }

    for _ in 0..seq_2_length {
        let symbol = Protein::match_with_u8(UTF8_SYMBOLS[rng.sample(u8_distr)]);
        sequence_2.push(symbol);
        let counter = counters.entry(symbol).or_insert(0);
        *counter += 1;
    }

    let mut freqs = Vec::<f64>::new();

    for (_, counter) in counters.iter() {
        freqs.push(*counter as f64 / total_length as f64);
    }

    ((sequence_1, sequence_2), arr1(&freqs[..]))
}

pub fn generate_sequence(length: &i32) -> (Vec<Protein>, Array1<f64>) {
    let mut sequence = Vec::<Protein>::new();

    let mut rng = thread_rng();
    let distr = Uniform::new(0, UTF8_SYMBOLS.len());

    let mut counters = HashMap::<Protein, i32>::new();

    for _ in 0..*length {
        let symbol = Protein::match_with_u8(UTF8_SYMBOLS[rng.sample(distr)]);

        sequence.push(symbol);

        let counter = counters.entry(symbol).or_insert(0);
        *counter += 1;
    }

    let mut freqs = Vec::<f64>::new();

    for (_, counter) in counters.iter() {
        freqs.push(*counter as f64 / *length as f64);
    }

    (sequence, arr1(&freqs[..]))
}

pub fn generate_test_sequences(
    amount: &i32,
    avg_length: &i32,
    error: &i32,
) -> (Vec<Vec<Protein>>, Array1<f64>) {
    let mut sequences = Vec::<Vec<Protein>>::new();
    let mut counters = HashMap::<Protein, i32>::new();
    let mut total_length = 0;

    let mut rng = thread_rng();
    let u8_distr = Uniform::new(0, UTF8_SYMBOLS.len());
    let error_distr = Uniform::new_inclusive(-error, error);

    for _ in 0..*amount {
        let mut sequence = Vec::<Protein>::new();
        let seq_length = *avg_length + rng.sample(error_distr);
        total_length += seq_length;

        for _ in 0..seq_length {
            let symbol = Protein::match_with_u8(UTF8_SYMBOLS[rng.sample(u8_distr)]);
            sequence.push(symbol);
            let counter = counters.entry(symbol).or_insert(0);
            *counter += 1;
        }

        sequences.push(sequence);
    }

    let mut freqs = Vec::<f64>::new();

    for (_, counter) in counters.iter() {
        freqs.push(*counter as f64 / total_length as f64);
    }

    (sequences, arr1(&freqs[..]))
}
