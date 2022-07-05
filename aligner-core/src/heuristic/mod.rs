use crate::alignment::{Alignment, PWMAlignment};
use crate::alignment_result::AlignmentResult;
use crate::pwm::PWMAligner;
use crate::simple::SimpleLocalAligner;
use crate::{AlignerTrait, AlignmentTrait, BioData, Error, Heuristics, Result};
use aligner_helpers::matrices::transform_matrix;
use ndarray::Array2;
use ndarray::Axis;

pub struct HeuristicAligner<T: BioData + Into<usize> + Copy + Eq> {
    pub query: Vec<T>,
    pub target: Vec<T>,
}

impl<T: BioData + Into<usize> + Copy + Eq> AlignerTrait<T, Alignment<T>> for HeuristicAligner<T> {
    fn from_str_seqs(query: &str, target: &str) -> Result<HeuristicAligner<T>> {
        Ok(HeuristicAligner {
            query: match T::str_to_vec(query) {
                Ok(query) => query,
                Err(err) => return Err(err),
            },
            target: match T::str_to_vec(target) {
                Ok(target) => target,
                Err(err) => return Err(err),
            },
        })
    }

    fn from_seqs(query: &[T], target: &[T]) -> Result<HeuristicAligner<T>> {
        Ok(HeuristicAligner {
            query: Vec::from(query),
            target: Vec::from(target),
        })
    }

    fn perform_alignment(
        &mut self,
        del: f64,
        ext: f64,
        matrix: &Array2<f64>,
        heuristics: Option<Heuristics>,
    ) -> Result<AlignmentResult<T, Alignment<T>>> {
        let mut params = match heuristics {
            Some(h) => h,
            None => return Err(Error::MissingArgument),
        };

        if (params.r_squared - 0 as f64).abs() < f64::EPSILON {
            params.r_squared = (matrix.len_of(Axis(0)) * matrix.len_of(Axis(1))) as f64;
        }

        let mut transformed_matrix =
            transform_matrix(matrix, params.kd, params.r_squared, &params.frequencies).unwrap();

        let mut max_f = 0f64;
        let mut result: AlignmentResult<T, Alignment<T>>;

        loop {
            let mut aligner = SimpleLocalAligner::from_seqs(&self.query, &self.target).unwrap();
            result = aligner
                .perform_alignment(del, ext, &transformed_matrix, None)
                .unwrap();

            if result.alignment.f > max_f {
                max_f = result.alignment.f;
                transformed_matrix = transform_matrix(
                    &result.alignment.get_frequency_matrix(),
                    params.kd,
                    params.r_squared,
                    &params.frequencies,
                )
                .unwrap();
            } else {
                result.matrix = Option::Some(transformed_matrix);
                return Ok(result);
            }
        }
    }
}

pub struct HeuristicPWMAligner<T: BioData + Into<usize> + Copy + Eq> {
    pub query: Vec<T>,
}

impl<T: BioData + Into<usize> + Copy + Eq> AlignerTrait<T, PWMAlignment<T>>
    for HeuristicPWMAligner<T>
{
    fn from_str_seqs(query: &str, _target: &str) -> Result<HeuristicPWMAligner<T>> {
        Ok(HeuristicPWMAligner {
            query: match T::str_to_vec(query) {
                Ok(query) => query,
                Err(err) => return Err(err),
            },
        })
    }

    fn from_seqs(query: &[T], _target: &[T]) -> Result<HeuristicPWMAligner<T>> {
        Ok(HeuristicPWMAligner {
            query: Vec::from(query),
        })
    }

    fn perform_alignment(
        &mut self,
        del: f64,
        ext: f64,
        matrix: &Array2<f64>,
        heuristics: Option<Heuristics>,
    ) -> Result<AlignmentResult<T, PWMAlignment<T>>> {
        let params = match heuristics {
            Some(h) => h,
            None => return Err(Error::MissingArgument),
        };

        let mut transformed_matrix =
            transform_matrix(matrix, params.kd, params.r_squared, &params.frequencies).unwrap();

        let mut max_f = 0f64;
        let mut result: AlignmentResult<T, PWMAlignment<T>>;

        loop {
            let mut aligner = PWMAligner::from_seqs(&self.query, &[]).unwrap();
            result = aligner
                .perform_alignment(del, ext, &transformed_matrix, None)
                .unwrap();

            if result.alignment.f > max_f {
                max_f = result.alignment.f;
                transformed_matrix = transform_matrix(
                    &result.alignment.get_frequency_matrix(),
                    params.kd,
                    params.r_squared,
                    &params.frequencies,
                )
                .unwrap();
            } else {
                result.matrix = Option::Some(transformed_matrix);
                return Ok(result);
            }
        }
    }
}
