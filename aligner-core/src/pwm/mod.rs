use crate::alignment::PWMAlignment;
use crate::alignment_result::AlignmentResult;
use crate::enums::Direction;
use crate::{AlignerTrait, BioData, Error, Heuristics, Result};
use ndarray::Array2;
use ndarray_stats::QuantileExt;
use std::marker::PhantomData;

pub struct PWMAligner<T: BioData + Into<usize> + Copy + Eq> {
    pub query: Vec<T>,
}

impl<T: BioData + Into<usize> + Copy + Eq> AlignerTrait<T, PWMAlignment<T>> for PWMAligner<T> {
    fn from_str_seqs(query: &str, _target: &str) -> Result<PWMAligner<T>> {
        Ok(PWMAligner {
            query: match T::str_to_vec(query) {
                Ok(query) => query,
                Err(err) => return Err(err),
            },
        })
    }

    fn from_seqs(query: &[T], _target: &[T]) -> Result<PWMAligner<T>> {
        Ok(PWMAligner {
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
        if heuristics.is_some() {
            return Err(Error::UnnecessaryArgument);
        };

        if matrix.shape()[0] != 4 {
            return Err(Error::MatrixShapeError);
        }

        let numbered_sequence: Vec<usize> = (1..=matrix.shape()[1]).collect();

        let dim = (self.query.len() + 1, numbered_sequence.len() + 1);

        let mut alignment_matrix = Array2::<f64>::zeros(dim);
        let mut direction_matrix =
            Array2::<Direction>::from_shape_fn(dim, |_| Direction::Beginning);

        let mut penalty = del;

        for elem_1 in numbered_sequence.iter() {
            for (y, elem_2) in self.query.iter().enumerate() {
                let y_real = y + 1;
                let seq_2_pos = (*elem_2).into();

                let assignment = Direction::get_direction_with_beginning(
                    alignment_matrix[[y_real - 1, *elem_1]] - penalty,
                    alignment_matrix[[y_real, elem_1 - 1]] - penalty,
                    alignment_matrix[[y_real - 1, elem_1 - 1]] + matrix[[seq_2_pos, elem_1 - 1]],
                );

                if assignment.1 != Direction::Beginning {
                    penalty = ext
                } else {
                    penalty = del
                }

                alignment_matrix[[y_real, *elem_1]] = assignment.0;
                direction_matrix[[y_real, *elem_1]] = assignment.1;
            }
        }

        let max_coords = alignment_matrix.argmax().unwrap();
        let (mut numbered_alignment, mut query_alignment) = (vec![], vec![]);
        let mut current_x = max_coords.1;
        let mut current_y = max_coords.0;

        loop {
            match direction_matrix[[current_y, current_x]] {
                Direction::Beginning => {
                    break;
                }
                Direction::Top => {
                    numbered_alignment.push(0);
                    query_alignment.push(self.query[current_y - 1]);
                    current_y -= 1;
                }
                Direction::Left => {
                    numbered_alignment.push(numbered_sequence[current_x - 1]);
                    query_alignment.push(T::blank());
                    current_x -= 1;
                }
                Direction::Diagonal => {
                    numbered_alignment.push(numbered_sequence[current_x - 1]);
                    query_alignment.push(self.query[current_y - 1]);
                    current_x -= 1;
                    current_y -= 1;
                }
            }
        }

        numbered_alignment.reverse();
        query_alignment.reverse();

        let f = *alignment_matrix.max().unwrap();

        Ok(AlignmentResult {
            alignment_matrix,
            direction_matrix,
            alignment: PWMAlignment {
                numbered: numbered_alignment,
                query: query_alignment,
                dim: matrix.shape()[1],
                coords: (
                    (current_x + 1, max_coords.1 + 1),
                    (current_y + 1, max_coords.0 + 1),
                ),
                f,
            },
            matrix: None,
            phantom: PhantomData,
        })
    }
}
