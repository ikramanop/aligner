use crate::simple::SimpleLocalAligner;
use crate::{AlignerTrait, Error, Protein, Result};
use ndarray::{arr1, Array1, Array2, Zip};
use ndarray_stats::SummaryStatisticsExt;
use rand::seq::SliceRandom;
use rand::{thread_rng, Rng};
use std::thread;

const MAXITER: i32 = 10000;
const THREADS: i32 = 10;
const SEQUENCES: i32 = 5000;
const THRESHOLD_GLOBAL: f64 = 1e-6;
const THRESHOLD_LOCAL: f64 = 1e-4;

#[derive(Debug, Clone)]
pub struct DistributionParams {
    pub k: f64,
    pub lambda: f64,
    pub h: f64,
}

impl DistributionParams {
    fn get_p_value(&self, query_length: usize, target_length: usize, score: f64) -> f64 {
        let l = (self.k * query_length as f64 * target_length as f64).ln() / self.h;
        let nn = (query_length as f64 - l) * (target_length as f64 - l);

        debug!(
            "score={}, k={}, lambda={}, h={}",
            score, self.k, self.lambda, self.h
        );

        1f64 - (-self.k * nn * (-self.lambda * score).exp()).exp()
    }
}

pub fn calculate_distribution_params(
    query_length: usize,
    target_lengths: &Array1<usize>,
    scores: &Array1<f64>,
) -> Result<DistributionParams> {
    if scores.len() != target_lengths.len() {
        return Err(Error::ValidationError);
    }

    let sd = match scores.central_moment(2) {
        Ok(sd) => sd,
        Err(_) => return Err(Error::ValidationError),
    };

    let lambda = 1f64 / sd;
    let mut h = 1f64;

    let n = target_lengths.len() as f64;

    let mut nn_array = target_lengths.mapv(|t| (query_length * t) as f64);

    let k = n / (&nn_array * scores.mapv(|score| (-lambda * score).exp())).sum();

    let mut log_likelihood = n * (lambda * k).ln()
        + Zip::from(&nn_array)
            .and(scores)
            .map_collect(|nn, score| nn.ln() - lambda * score - k * nn * (-lambda * score).exp())
            .sum();

    let mut active_target_lengths = target_lengths.clone();
    let mut active_scores = scores.clone();

    for _ in 0..=MAXITER {
        let (k, lambda) = estimate_k_and_lambda_by_parameters(
            query_length,
            &active_target_lengths,
            &active_scores,
            k,
            lambda,
            h,
        );

        h = estimate_h_by_parameters(
            query_length,
            &active_target_lengths,
            &active_scores,
            k,
            lambda,
            h,
        );

        nn_array = target_lengths.mapv(|t| {
            let l = (k * query_length as f64 * t as f64).ln() / h;

            (query_length as f64 - l) * (t as f64 - l)
        });

        let log_likelihood_new = n * (lambda * k).log10()
            + (nn_array.mapv(|nn| nn.log10())
                - lambda * scores
                - k * &nn_array * scores.mapv(|score| (-lambda * score).exp()))
            .sum();

        if (log_likelihood_new - log_likelihood).abs() / log_likelihood < THRESHOLD_GLOBAL {
            return Ok(DistributionParams { k, lambda, h });
        }

        log_likelihood = log_likelihood_new;

        let mut target_lengths_buffer: Vec<usize> = vec![];
        let mut scores_buffer: Vec<f64> = vec![];

        Zip::from(scores)
            .and(target_lengths)
            .and(&nn_array)
            .for_each(|score, t, nn| {
                if n * (1f64 - (-k * nn * (-lambda * score).exp()).exp()) >= 1f64 {
                    target_lengths_buffer.push(*t);
                    scores_buffer.push(*score);
                }
            });

        active_target_lengths = arr1(&target_lengths_buffer);
        active_scores = arr1(&scores_buffer);
    }

    Ok(DistributionParams { k, lambda, h })
}

fn estimate_k_and_lambda_by_parameters(
    query_length: usize,
    target_lengths: &Array1<usize>,
    scores: &Array1<f64>,
    old_k: f64,
    old_lambda: f64,
    h: f64,
) -> (f64, f64) {
    let mut k = old_k;
    let mut lambda = old_lambda;

    let n = target_lengths.len() as f64;

    let mut nn_array = target_lengths.mapv(|t| {
        let l = (k * query_length as f64 * t as f64).ln() / h;

        (query_length as f64 - l) * (t as f64 - l)
    });

    let mut exponential_scores = scores.mapv(|score| (-lambda * score).exp());
    let mut sum = (&nn_array * &exponential_scores).sum();
    let mut weighted_sum = (&nn_array * scores * &exponential_scores).sum();

    for i in 0..=MAXITER {
        let lambda_f = 1f64 / lambda - scores.sum() / n + weighted_sum / sum;
        let lambda_fd = -lambda.powi(-2)
            - (&nn_array * &scores.mapv(|score| score * score) * &exponential_scores).sum() / sum
            + (weighted_sum / sum).powi(2);

        if !f64::is_finite(lambda_f) || !f64::is_finite(lambda_fd) {
            return (k, lambda);
        }

        let new_lambda = lambda - lambda_f / lambda_fd;

        exponential_scores = scores.mapv(|score| (-lambda * score).exp());
        sum = (&nn_array * &exponential_scores).sum();
        weighted_sum = (&nn_array * scores * &exponential_scores).sum();

        let new_k = n / sum;

        if !f64::is_finite(new_k) || new_k <= 0f64 {
            return (k, lambda);
        }

        k = new_k;
        lambda = new_lambda;

        trace!("Iteration {}, f={}, fd={}", i + 1, lambda_f, lambda_fd);

        trace!("Iteration {}, lambda={}, k={}", i + 1, lambda, k);

        if lambda_f.abs() < THRESHOLD_LOCAL {
            return (k, lambda);
        }

        nn_array = target_lengths.mapv(|t| {
            let l = (k * query_length as f64 * t as f64).ln() / h;

            (query_length as f64 - l) * (t as f64 - l)
        });
    }

    (k, lambda)
}

fn estimate_h_by_parameters(
    query_length: usize,
    target_lengths: &Array1<usize>,
    scores: &Array1<f64>,
    k: f64,
    lambda: f64,
    old_h: f64,
) -> f64 {
    let mut h = old_h;

    for i in 0..=MAXITER {
        let l_array = target_lengths.mapv(|t| (k * query_length as f64 * t as f64).ln() / h);

        let nn_array = Zip::from(target_lengths)
            .and(&l_array)
            .map_collect(|t, l| (query_length as f64 - l) * (*t as f64 - l));

        let a_array = 2f64 * &l_array - query_length as f64 - target_lengths.mapv(|t| t as f64);
        let b_array = 1f64 / &nn_array - k * scores.mapv(|score| (-lambda * score).exp());
        let c_array = -&l_array / h;

        let h_g = (&a_array * &b_array * &c_array).sum();
        let h_gd = (2f64 * &b_array * &c_array.mapv(|c| c * c)
            - (&a_array * &c_array / &nn_array).mapv(|u| u * u)
            - 2f64 * &a_array * &b_array * &c_array / h)
            .sum();

        if h_g.abs() < THRESHOLD_LOCAL {
            return h;
        }

        if h_gd > 0f64 {
            if h_g > 0f64 {
                h *= 2f64;
            } else {
                h /= 2f64;
            }
        } else if h_g <= 0f64 {
            h /= 2f64;
        } else {
            h -= h_g / h_gd;
        }

        trace!("Iteration {}: h={}", i + 1, h);
    }

    h
}

pub fn calculate_p_value(
    query: Vec<Protein>,
    target: Vec<Protein>,
    initial_score: f64,
    del: f64,
    ins: f64,
    matrix: &Array2<f64>,
) -> Result<f64> {
    debug!("Started shuffling and calculating input sequences!!");

    let mut scores = vec![initial_score];
    let mut lengths = vec![target.len()];

    let mut threads = vec![];

    for i in 0..THREADS {
        threads.push(thread::spawn({
            let query_clone = query.clone();
            let target_clone = target.clone();
            let matrix_clone = matrix.clone();
            let mut scores_scoped = Vec::<f64>::new();
            let mut lengths_scoped = Vec::<usize>::new();
            move || {
                let mut limit = SEQUENCES / THREADS;
                if i == 5 {
                    limit = SEQUENCES - (SEQUENCES / THREADS * (THREADS - 1)) - 1;
                }
                for _ in 0..limit {
                    let new_seq = shuffle_and_randomize_sequence(&target_clone);

                    let mut aligner =
                        SimpleLocalAligner::<Protein>::from_seqs(&query_clone, &new_seq).unwrap();

                    scores_scoped.push(
                        aligner
                            .perform_alignment(del, ins, &matrix_clone, None)
                            .unwrap()
                            .alignment
                            .f,
                    );
                    lengths_scoped.push(new_seq.len());
                }

                (scores_scoped, lengths_scoped)
            }
        }))
    }

    for t in threads {
        let mut items = t.join().unwrap();
        scores.append(&mut items.0);
        lengths.append(&mut items.1);
    }

    debug!("Calculating distribution params");

    let p_value = match calculate_distribution_params(
        query.len(),
        &Array1::from_vec(lengths),
        &Array1::from_vec(scores),
    ) {
        Ok(params) => params,
        Err(err) => return Err(err),
    }
    .get_p_value(query.len(), target.len(), initial_score);

    Ok(p_value)
}

fn shuffle_and_randomize_sequence(sequence: &[Protein]) -> Vec<Protein> {
    let mut thread_rng = thread_rng();

    let lock = thread_rng.gen_range(0..7);

    let target_bytes = &sequence[..sequence.len() - lock];
    let mut target_vec = Vec::<Protein>::from(target_bytes);

    target_vec.shuffle(&mut thread_rng);

    target_vec
}
