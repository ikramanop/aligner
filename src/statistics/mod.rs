use ndarray::{Array1, Zip};
use ndarray_stats::SummaryStatisticsExt;
use std::error::Error;
use std::fmt;
use std::result::Result;

const MAXITER: i32 = 10000;
const THRESHOLD_GLOBAL: f64 = 1e-6;
const THRESHOLD_LOCAL: f64 = 1e-4;

#[derive(Debug)]
pub struct DistributionParams {
    pub k: f64,
    pub lambda: f64,
    pub h: f64,
}

#[derive(Debug, Clone)]
struct CalculationError;

impl fmt::Display for CalculationError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "error calculating params")
    }
}

impl Error for CalculationError {}

#[derive(Debug, Clone)]
struct ValidationError;

impl fmt::Display for ValidationError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "error with input data")
    }
}

impl Error for ValidationError {}

pub fn calculate_distribution_params(
    query_length: &usize,
    target_lengths: &Array1<usize>,
    scores: &Array1<f64>,
) -> Result<DistributionParams, impl Error> {
    if scores.len() != target_lengths.len() {
        return Err(ValidationError);
    }

    let sd = match scores.central_moment(2) {
        Ok(sd) => sd,
        Err(_) => return Err(ValidationError),
    };

    let lambda = 1f64 / sd;
    let mut h = 1f64;

    let n = target_lengths.len() as f64;

    let mut nn_array = target_lengths.mapv(|t| (query_length * t) as f64);

    let k = n / (&nn_array * scores.mapv(|score| score.exp())).sum();

    let mut log_likelihood = n * (lambda * k).ln()
        + Zip::from(&nn_array)
            .and(scores)
            .map_collect(|nn, score| nn.ln() - lambda * score - k * nn * score.exp())
            .sum();

    let mut active_target_lengths = target_lengths.clone();
    let mut active_scores = scores.clone();

    for _ in 0..=MAXITER {
        let (k, lambda) = estimate_k_and_lambda_by_parameters(
            query_length,
            &active_target_lengths,
            &active_scores,
            &k,
            &lambda,
            &h,
        );

        h = estimate_h_by_parameters(
            query_length,
            &active_target_lengths,
            &active_scores,
            &k,
            &lambda,
            &h,
        );

        nn_array = target_lengths.mapv(|t| {
            let l = (k * *query_length as f64 * t as f64).ln() / h;

            (*query_length as f64 - l) * (t as f64 - l)
        });

        let log_likelihood_new = n * (lambda * k).log10()
            + (nn_array.mapv(|nn| nn.log10())
                - lambda * &active_scores
                - k * &nn_array * active_scores.mapv(|score| -lambda * score.exp()))
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

        active_target_lengths = Array1::from_vec(target_lengths_buffer);
        active_scores = Array1::from_vec(scores_buffer);
    }

    Ok(DistributionParams { k, lambda, h })
}

fn estimate_k_and_lambda_by_parameters(
    query_length: &usize,
    target_lengths: &Array1<usize>,
    scores: &Array1<f64>,
    old_k: &f64,
    old_lambda: &f64,
    h: &f64,
) -> (f64, f64) {
    let mut k = *old_k;
    let mut lambda = *old_lambda;

    let n = target_lengths.len() as f64;

    for i in 0..=MAXITER {
        let nn_array = target_lengths.mapv(|t| {
            let l = (k * *query_length as f64 * t as f64).ln() / h;

            (*query_length as f64 - l) * (t as f64 - l)
        });

        let exponential_scores = scores.mapv(|score| (-lambda * score).exp());
        let sum = (&nn_array * &exponential_scores).sum();
        let weighted_sum = (&nn_array * scores * &exponential_scores).sum();

        k = n / sum;

        let lambda_f = 1f64 / lambda - scores.sum() / n + weighted_sum / sum;
        let lambda_fd = -1f64 / lambda.powi(2)
            - (nn_array * &scores.mapv(|score| score * score) * &exponential_scores).sum() / sum
            + (weighted_sum / sum).powi(2);

        debug!("Iteration {}: f={}, fd={}", i + 1, lambda_f, lambda_fd);

        if lambda_f.abs() < THRESHOLD_LOCAL {
            return (k, lambda);
        }

        lambda -= lambda_f / lambda_fd;

        debug!("Iteration {}, lambda={}, k={}", i + 1, lambda, k);
    }

    (k, lambda)
}

fn estimate_h_by_parameters(
    query_length: &usize,
    target_lengths: &Array1<usize>,
    scores: &Array1<f64>,
    k: &f64,
    lambda: &f64,
    old_h: &f64,
) -> f64 {
    let mut h = *old_h;

    for i in 0..=MAXITER {
        let l_array = target_lengths.mapv(|t| (k * *query_length as f64 * t as f64).ln() / h);

        let nn_array = Zip::from(target_lengths)
            .and(&l_array)
            .map_collect(|t, l| (*query_length as f64 - l) * (*t as f64 - l));

        let a_array = 2f64 * &l_array - *query_length as f64 - target_lengths.mapv(|t| t as f64);
        let b_array = 1f64 / &nn_array - *k * scores.mapv(|score| (-lambda * score).exp());
        let c_array = -&l_array / h;

        let h_g = (&a_array * &b_array * &c_array).sum();
        let h_gd = (2f64 * &b_array * &c_array.mapv(|c| c.powi(2))
            - (&a_array * &c_array / &nn_array).mapv(|u| u.powi(2))
            - 2f64 * &a_array * &b_array * &c_array / h)
            .sum();

        debug!("Iteration {}: g={}, gd={}", i + 1, h_g, h_gd);

        if h_g.abs() < THRESHOLD_LOCAL {
            return h;
        }

        if h_gd > 0f64 {
            if h_g > 0f64 {
                h *= 2f64;
            } else {
                h /= 2f64;
            }
        } else {
            let hh = h;
            h -= h_g / h_gd;

            if h <= 0f64 {
                h = hh / 2f64;
            }
        }

        debug!("Iteration {}: h={}", i + 1, h);
    }

    h
}
