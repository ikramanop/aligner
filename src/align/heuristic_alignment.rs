use crate::align::enums::Direction;
use crate::align::enums::Protein;
use crate::matrices::transform_matrix;
use ndarray::Array1;
use ndarray::Array2;
use ndarray::Axis;

pub struct HeuristicPairwiseAlignmentTool {
    sequences_pair: (Vec<Protein>, Vec<Protein>),
}

pub struct HeuristicPairwiseAlignmentResult {
    pub alignment_matrix: Array2<f64>,
    pub direction_matrix: Array2<Direction>,
    pub max_f: f64,
    optimal_alignment: (Vec<Protein>, Vec<Protein>),
}

impl HeuristicPairwiseAlignmentResult {
    pub fn get_frequency_matrix(&self) -> Array2<f64> {
        let mut frequency_matrix = Array2::<f64>::zeros((20, 20));

        for (x, y) in self
            .optimal_alignment
            .0
            .iter()
            .zip(self.optimal_alignment.1.iter())
        {
            if *x != Protein::Blank && *y != Protein::Blank {
                frequency_matrix[[*y as usize, *x as usize]] += 1f64;
            }
        }

        frequency_matrix
    }
}

impl HeuristicPairwiseAlignmentTool {
    pub fn from_sequences_pair(
        sequences: (Vec<Protein>, Vec<Protein>),
    ) -> HeuristicPairwiseAlignmentTool {
        HeuristicPairwiseAlignmentTool {
            sequences_pair: sequences,
        }
    }

    pub fn local_alignment(
        &mut self,
        del: f64,
        kd_value: f64,
        mut r_squared_value: f64,
        matrix: &Array2<f64>,
        frequences: &Array1<f64>,
    ) -> (f64, Array2<f64>, (Vec<Protein>, Vec<Protein>)) {
        if (r_squared_value - 0 as f64).abs() < f64::EPSILON {
            r_squared_value = (matrix.len_of(Axis(0)) * matrix.len_of(Axis(1))) as f64;
        }

        let mut transformed_matrix =
            transform_matrix(matrix, &kd_value, &r_squared_value, frequences).unwrap();

        let mut max_f = 0f64;
        let mut frequency_matrix: Array2<f64>;
        let optimal_alignment: (Vec<Protein>, Vec<Protein>);

        loop {
            let result = HeuristicPairwiseAlignmentTool::local_alignment_step(
                &self.sequences_pair.0,
                &self.sequences_pair.1,
                &del,
                &transformed_matrix,
            );

            frequency_matrix = result.get_frequency_matrix();

            println!("Max F: {}", result.max_f);

            if result.max_f > max_f {
                max_f = result.max_f;
                transformed_matrix =
                    transform_matrix(&frequency_matrix, &kd_value, &r_squared_value, frequences)
                        .unwrap();
            } else {
                optimal_alignment = result.optimal_alignment;
                break;
            }
        }

        (
            max_f,
            transform_matrix(&frequency_matrix, &kd_value, &r_squared_value, frequences).unwrap(),
            optimal_alignment,
        )
    }

    fn local_alignment_step(
        sequence_1: &[Protein],
        sequence_2: &[Protein],
        del: &f64,
        matrix: &Array2<f64>,
    ) -> HeuristicPairwiseAlignmentResult {
        let mut alignment_matrix =
            Array2::<f64>::zeros((sequence_2.len() + 1, sequence_1.len() + 1));

        let mut direction_matrix = Array2::<Direction>::from_shape_fn(
            (sequence_2.len() + 1, sequence_1.len() + 1),
            |_| Direction::Beginning,
        );

        let mut max_f: f64 = 0f64;
        let mut max_x: usize = 0;
        let mut max_y: usize = 0;

        for (x, elem_1) in sequence_1.iter().enumerate() {
            for (y, elem_2) in sequence_2.iter().enumerate() {
                let x_real = x + 1;
                let y_real = y + 1;

                let seq_1_pos = *elem_1 as usize;
                let seq_2_pos = *elem_2 as usize;

                let top = alignment_matrix[[y_real - 1, x_real]] - del;
                let left = alignment_matrix[[y_real, x_real - 1]] - del;
                let diagonal =
                    alignment_matrix[[y_real - 1, x_real - 1]] + matrix[[seq_2_pos, seq_1_pos]];

                let max = f64::max(f64::max(f64::max(top, left), diagonal), 0f64);

                alignment_matrix[[y_real, x_real]] = max;

                if max == 0f64 {
                    direction_matrix[[y_real, x_real]] = Direction::Beginning
                } else if (max - top).abs() < f64::EPSILON {
                    direction_matrix[[y_real, x_real]] = Direction::Top
                } else if (max - left).abs() < f64::EPSILON {
                    direction_matrix[[y_real, x_real]] = Direction::Left
                } else if (max - diagonal).abs() < f64::EPSILON {
                    direction_matrix[[y_real, x_real]] = Direction::Diagonal
                }

                if max >= max_f {
                    max_f = max;
                    max_x = x;
                    max_y = y;
                }
            }
        }

        let mut current_x = max_x;
        let mut current_y = max_y;
        let (mut optimal_alignment_1, mut optimal_alignment_2) =
            (vec![sequence_1[max_x]], vec![sequence_2[max_y]]);
        loop {
            match direction_matrix[[current_y, current_x]] {
                Direction::Beginning => break,
                Direction::Top => {
                    optimal_alignment_1.push(Protein::Blank);
                    optimal_alignment_2.push(sequence_2[current_y - 1]);
                    current_y -= 1;
                }
                Direction::Left => {
                    optimal_alignment_1.push(sequence_1[current_x - 1]);
                    optimal_alignment_2.push(Protein::Blank);
                    current_x -= 1;
                }
                Direction::Diagonal => {
                    optimal_alignment_1.push(sequence_1[current_x - 1]);
                    optimal_alignment_2.push(sequence_2[current_y - 1]);
                    current_x -= 1;
                    current_y -= 1;
                }
            }
        }

        optimal_alignment_1.reverse();
        optimal_alignment_2.reverse();

        HeuristicPairwiseAlignmentResult {
            alignment_matrix,
            direction_matrix,
            max_f,
            optimal_alignment: (optimal_alignment_1, optimal_alignment_2),
        }
    }
}
