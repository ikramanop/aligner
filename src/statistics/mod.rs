use ndarray::Array1;
use ndarray_stats::SummaryStatisticsExt;
use std::error::Error;
use std::fmt;
use std::result::Result;

const MAXITER: i32 = 10000;
const THRESHOLD: f64 = 1e-6;

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

pub fn calculate_distribution_params<E>(
    scores: &Array1<f64>,
    t_lengths: Array1<i32>,
    q_length: usize,
) -> Result<DistributionParams, impl Error> {
    if scores.len() != t_lengths.len() {
        return Err(ValidationError);
    }

    let sd = match scores.central_moment(2) {
        Ok(sd) => sd,
        Err(_) => return Err(ValidationError),
    };
    let mut lambda = 1f64 / sd;
    let mut h = 1f64;

    let n = scores.len() as f64;

    let mut Ns = t_lengths.mapv(|t| q_length as f64 * t as f64);

    let mut k = n / (&Ns * scores.mapv(|score| score.exp())).sum();

    let log_likelihood = n * (lambda * k).log10()
        + (Ns.mapv(|N| N.log10())
            - lambda * scores
            - k * Ns * scores.mapv(|score| -lambda * score.exp()))
        .sum();

    let mut active_scores = scores.clone();

    for _ in 0..=MAXITER {
        Ns = t_lengths.mapv(move |t| {
            let l = (k * q_length as f64 * t as f64).log10() / h;

            (q_length as f64 - l) * (t as f64 - l)
        });

        let sum = (&Ns * &active_scores.mapv(|score| -lambda * score.exp())).sum();
        let weighted_sum =
            (&Ns * &active_scores * &active_scores.mapv(|score| -lambda * score.exp())).sum();

        k = n / sum;

        let lambda_f = 1f64 / lambda - &active_scores.sum() / n + weighted_sum / sum;
        let lambda_fd = -1f64 / lambda.powi(2)
            - (&Ns
                * &active_scores.mapv(|score| score * score)
                * active_scores.mapv(|score| -lambda * score.exp()))
            .sum()
                / sum
            + (weighted_sum / sum).powi(2);

        lambda = lambda - lambda_f / lambda_fd;

        h = estimate_h();

        let log_likelihood_new = n * (lambda * k).log10()
            + (Ns.mapv(|N| N.log10())
                - lambda * &active_scores
                - k * &Ns * active_scores.mapv(|score| -lambda * score.exp()))
            .sum();

        if (log_likelihood_new - log_likelihood).abs() / log_likelihood < THRESHOLD {
            return Ok(DistributionParams { k, lambda, h });
        }

        active_scores = active_scores.mapv(|_| 0f64);
        for (i, score) in scores.indexed_iter() {
            if n * (1f64 - (-k * &Ns[i] * (-lambda * score).exp()).exp()) >= 1f64 {
                active_scores[i] = *score;
            }
        }
    }

    Ok(DistributionParams { k, lambda, h })
}

fn estimate_h() -> f64 {
    0f64
}
