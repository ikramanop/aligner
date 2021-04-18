use crate::align::enums::{Direction, Protein};
use ndarray::Array2;
use std::cmp::max;
use std::str::from_utf8;

#[derive(Debug, Clone)]
pub struct AlignResult {
    alignment_matrix: Array2<i32>,
    direction_matrix: Array2<Direction>,
    max_f: i32,
    optimal_alignment: String,
}

#[derive(Debug, Clone)]
pub struct SimpleAligner {
    u8_sequence_1: Vec<u8>,
    u8_sequence_2: Vec<u8>,
    pub protein_sequence_1: Vec<Protein>,
    pub protein_sequence_2: Vec<Protein>,
}

pub trait Aligner {
    // fn global_alignment(&mut self, del: &i32, matrix: &Array2<i32>) -> AlignResult;
    fn from_seqs(seq_1: &[u8], seq_2: &[u8]) -> Self;
    fn local_alignment(&mut self, del: &i32, matrix: &Array2<i32>) -> AlignResult;
    fn get_symbolic(&mut self) -> (String, String);
}

impl Aligner for SimpleAligner {
    fn from_seqs(seq_1: &[u8], seq_2: &[u8]) -> SimpleAligner {
        SimpleAligner {
            u8_sequence_1: seq_1.to_vec(),
            u8_sequence_2: seq_2.to_vec(),
            protein_sequence_1: Protein::u8_vec_to_protein_vec(seq_1),
            protein_sequence_2: Protein::u8_vec_to_protein_vec(seq_2),
        }
    }

    fn local_alignment(&mut self, del: &i32, matrix: &Array2<i32>) -> AlignResult {
        let mut alignment_matrix =
            Array2::<i32>::zeros((self.u8_sequence_2.len() + 1, self.u8_sequence_1.len() + 1));
        let mut direction_matrix = Array2::<Direction>::from_shape_fn(
            (self.u8_sequence_2.len() + 1, self.u8_sequence_1.len() + 1),
            |_| Direction::Beginning,
        );

        for (x, elem_1) in self.protein_sequence_1.iter().enumerate() {
            for (y, elem_2) in self.protein_sequence_2.iter().enumerate() {
                let x_real = x + 1;
                let y_real = y + 1;

                let seq_1_pos = *elem_1 as usize;
                let seq_2_pos = *elem_2 as usize;

                let up = alignment_matrix[[y_real - 1, x_real]] - del;
                let left = alignment_matrix[[y_real, x_real - 1]] - del;
                let diagonal =
                    alignment_matrix[[y_real - 1, x_real - 1]] + matrix[[seq_2_pos, seq_1_pos]];

                let max = max(max(max(up, left), diagonal), 0);
                alignment_matrix[[y_real, x_real]] = max;

                if max == 0 {
                    direction_matrix[[y_real, x_real]] = Direction::Beginning
                } else if max == up {
                    direction_matrix[[y_real, x_real]] = Direction::Up
                } else if max == left {
                    direction_matrix[[y_real, x_real]] = Direction::Left
                } else if max == diagonal {
                    direction_matrix[[y_real, x_real]] = Direction::Diagonal
                }
            }
        }

        AlignResult {
            alignment_matrix: alignment_matrix,
            direction_matrix: direction_matrix,
            max_f: 10,
            optimal_alignment: String::from("ABABA"),
        }
    }

    fn get_symbolic(&mut self) -> (String, String) {
        let symbolic_sequence_1 = match from_utf8(&self.u8_sequence_1) {
            Ok(sequence) => String::from(sequence),
            Err(err) => panic!("Error converting sequence to str: {}", err),
        };
        let symbolic_sequence_2 = match from_utf8(&self.u8_sequence_2) {
            Ok(sequence) => String::from(sequence),
            Err(err) => panic!("Error converting sequence to str: {}", err),
        };

        (symbolic_sequence_1, symbolic_sequence_2)
    }
}
