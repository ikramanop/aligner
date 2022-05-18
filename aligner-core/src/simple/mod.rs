use crate::alignment::Alignment;
use crate::alignment_result::AlignmentResult;
use crate::enums::Direction;
use crate::{AlignerTrait, BioData, Error, Heuristics, Result};
use ndarray::Array2;
use ndarray_stats::QuantileExt;
use std::marker::PhantomData;

pub struct SimpleGlobalAligner<T: BioData + Into<usize> + Copy + Eq> {
    pub query: Vec<T>,
    pub target: Vec<T>,
}

pub struct SimpleLocalAligner<T: BioData + Into<usize> + Copy + Eq> {
    pub query: Vec<T>,
    pub target: Vec<T>,
}

impl<T: BioData + Into<usize> + Copy + Eq> AlignerTrait<T, Alignment<T>>
    for SimpleGlobalAligner<T>
{
    fn from_str_seqs(query: &str, target: &str) -> Result<SimpleGlobalAligner<T>> {
        Ok(SimpleGlobalAligner {
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

    fn from_seqs(query: &[T], target: &[T]) -> Result<SimpleGlobalAligner<T>> {
        Ok(SimpleGlobalAligner {
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
        if heuristics.is_some() {
            return Err(Error::UnnecessaryArgument);
        };

        let dim = (self.target.len() + 1, self.query.len() + 1);

        let mut alignment_matrix = Array2::<f64>::zeros(dim);
        let mut direction_matrix =
            Array2::<Direction>::from_shape_fn(dim, |_| Direction::Beginning);

        for x in 1..self.query.len() + 1 {
            alignment_matrix[[0, x]] = -(x as f64) * del;
            direction_matrix[[0, x]] = Direction::Left
        }

        for y in 1..self.target.len() + 1 {
            alignment_matrix[[y, 0]] = -(y as f64) * del;
            direction_matrix[[y, 0]] = Direction::Top;
        }

        alignment_matrix[[0, self.query.len()]] = -(self.query.len() as f64 + 1f64) * del;
        alignment_matrix[[self.target.len(), 0]] = -(self.target.len() as f64 + 1f64) * del;

        let mut penalty = del;

        for (x, elem_1) in self.query.iter().enumerate() {
            for (y, elem_2) in self.target.iter().enumerate() {
                let x_real = x + 1;
                let y_real = y + 1;

                let seq_1_pos: usize = (*elem_1).into();
                let seq_2_pos: usize = (*elem_2).into();

                let assignment = Direction::get_direction(
                    alignment_matrix[[y_real - 1, x_real]] - penalty,
                    alignment_matrix[[y_real, x_real - 1]] - penalty,
                    alignment_matrix[[y_real - 1, x_real - 1]] + matrix[[seq_2_pos, seq_1_pos]],
                );

                if assignment.1 != Direction::Beginning {
                    penalty = ext
                } else {
                    penalty = del
                }

                alignment_matrix[[y_real, x_real]] = assignment.0;
                direction_matrix[[y_real, x_real]] = assignment.1;
            }
        }

        let mut current_x = self.query.len();
        let mut current_y = self.target.len();

        let (mut query_aligned, mut target_aligned) = (
            vec![*self.query.last().unwrap()],
            vec![*self.target.last().unwrap()],
        );

        loop {
            match direction_matrix[[current_y, current_x]] {
                Direction::Beginning => break,
                Direction::Top => {
                    query_aligned.push(T::blank());
                    target_aligned.push(self.target[current_y - 1]);
                    current_y -= 1;
                }
                Direction::Left => {
                    query_aligned.push(self.query[current_x - 1]);
                    target_aligned.push(T::blank());
                    current_x -= 1;
                }
                Direction::Diagonal => {
                    query_aligned.push(self.query[current_x - 1]);
                    target_aligned.push(self.target[current_y - 1]);
                    current_x -= 1;
                    current_y -= 1;
                }
            }
        }

        query_aligned.reverse();
        target_aligned.reverse();

        Ok(AlignmentResult {
            alignment_matrix,
            direction_matrix,
            alignment: Alignment {
                query: query_aligned,
                target: target_aligned,
                coords: ((1, self.query.len()), (1, self.target.len())),
                f: 0f64,
            },
            phantom: PhantomData,
            matrix: None,
        })
    }
}

impl<T: BioData + Into<usize> + Copy + Eq> AlignerTrait<T, Alignment<T>> for SimpleLocalAligner<T> {
    fn from_str_seqs(query: &str, target: &str) -> Result<SimpleLocalAligner<T>> {
        Ok(SimpleLocalAligner {
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

    fn from_seqs(query: &[T], target: &[T]) -> Result<SimpleLocalAligner<T>> {
        Ok(SimpleLocalAligner {
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
        if heuristics.is_some() {
            return Err(Error::UnnecessaryArgument);
        };

        let dim = (self.target.len() + 1, self.query.len() + 1);

        let mut alignment_matrix = Array2::<f64>::zeros(dim);
        let mut direction_matrix =
            Array2::<Direction>::from_shape_fn(dim, |_| Direction::Beginning);

        let mut penalty = del;

        for (x, elem_1) in self.query.iter().enumerate() {
            for (y, elem_2) in self.target.iter().enumerate() {
                let x_real = x + 1;
                let y_real = y + 1;

                let seq_1_pos = (*elem_1).into();
                let seq_2_pos = (*elem_2).into();

                let assignment = Direction::get_direction_with_beginning(
                    alignment_matrix[[y_real - 1, x_real]] - penalty,
                    alignment_matrix[[y_real, x_real - 1]] - penalty,
                    alignment_matrix[[y_real - 1, x_real - 1]] + matrix[[seq_2_pos, seq_1_pos]],
                );

                if assignment.1 != Direction::Beginning {
                    penalty = ext
                } else {
                    penalty = del
                }

                alignment_matrix[[y_real, x_real]] = assignment.0;
                direction_matrix[[y_real, x_real]] = assignment.1;
            }
        }

        let max_coords = alignment_matrix.argmax().unwrap();
        let (mut query_alignment, mut target_alignment) = (
            vec![self.query[max_coords.1 - 1]],
            vec![self.target[max_coords.0 - 1]],
        );
        let mut current_x = max_coords.1;
        let mut current_y = max_coords.0;

        loop {
            match direction_matrix[[current_y, current_x]] {
                Direction::Beginning => {
                    break;
                }
                Direction::Top => {
                    query_alignment.push(T::blank());
                    target_alignment.push(self.target[current_y - 1]);
                    current_y -= 1;
                }
                Direction::Left => {
                    query_alignment.push(self.query[current_x - 1]);
                    target_alignment.push(T::blank());
                    current_x -= 1;
                }
                Direction::Diagonal => {
                    query_alignment.push(self.query[current_x - 1]);
                    target_alignment.push(self.target[current_y - 1]);
                    current_x -= 1;
                    current_y -= 1;
                }
            }
        }

        query_alignment.reverse();
        target_alignment.reverse();

        let f = *alignment_matrix.max().unwrap();

        Ok(AlignmentResult {
            alignment_matrix,
            direction_matrix,
            alignment: Alignment {
                query: query_alignment,
                target: target_alignment,
                coords: (
                    (current_x + 1, max_coords.1 + 1),
                    (current_y + 1, max_coords.0 + 1),
                ),
                f,
            },
            phantom: PhantomData,
            matrix: None,
        })
    }
}
