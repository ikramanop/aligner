use crate::align::enums::{Direction, Protein};
use ndarray::Array2;
use std::cmp::max;
use std::str::from_utf8;

pub trait AlignmentResult {
    fn represent(&self);
    fn get_alignment_matrix(&self) -> &Array2<i32>;
    fn get_direction_matrix(&self) -> &Array2<Direction>;
    fn get_optimal_alignment(&self) -> (&Vec<Protein>, &Vec<Protein>);
}

pub struct GlobalAlignmentResult {
    pub alignment_matrix: Array2<i32>,
    pub direction_matrix: Array2<Direction>,
    pub optimal_alignment: (Vec<Protein>, Vec<Protein>),
}

impl AlignmentResult for GlobalAlignmentResult {
    fn represent(&self) {
        println!(
            "Optimal Global Alignment:\n{}\n{}",
            from_utf8(&Protein::protein_vec_to_u8_vec(&self.optimal_alignment.0)).unwrap(),
            from_utf8(&Protein::protein_vec_to_u8_vec(&self.optimal_alignment.1)).unwrap()
        );
    }

    fn get_alignment_matrix(&self) -> &Array2<i32> {
        &self.alignment_matrix
    }

    fn get_direction_matrix(&self) -> &Array2<Direction> {
        &self.direction_matrix
    }

    fn get_optimal_alignment(&self) -> (&Vec<Protein>, &Vec<Protein>) {
        (&self.optimal_alignment.0, &self.optimal_alignment.1)
    }
}

#[derive(Debug)]
pub struct LocalAlignmentResult {
    pub alignment_matrix: Array2<i32>,
    pub direction_matrix: Array2<Direction>,
    #[allow(dead_code)]
    pub max_f: i32,
    pub optimal_alignment: (Vec<Protein>, Vec<Protein>),
}

impl AlignmentResult for LocalAlignmentResult {
    fn represent(&self) {
        println!(
            "Optimal Local Alignment:\n{}\n{}",
            from_utf8(&Protein::protein_vec_to_u8_vec(&self.optimal_alignment.0)).unwrap(),
            from_utf8(&Protein::protein_vec_to_u8_vec(&self.optimal_alignment.1)).unwrap()
        );
    }

    fn get_alignment_matrix(&self) -> &Array2<i32> {
        &self.alignment_matrix
    }

    fn get_direction_matrix(&self) -> &Array2<Direction> {
        &self.direction_matrix
    }

    fn get_optimal_alignment(&self) -> (&Vec<Protein>, &Vec<Protein>) {
        (&self.optimal_alignment.0, &self.optimal_alignment.1)
    }
}

pub struct SimpleAligner {
    u8_sequence_1: Vec<u8>,
    u8_sequence_2: Vec<u8>,
    pub protein_sequence_1: Vec<Protein>,
    pub protein_sequence_2: Vec<Protein>,
}

pub trait Aligner {
    fn from_seqs(seq_1: &[u8], seq_2: &[u8]) -> Self;
    fn local_alignment(&mut self, del: &i32, matrix: &Array2<i32>) -> Box<dyn AlignmentResult>;
    fn global_alignment(&mut self, del: &i32, matrix: &Array2<i32>) -> Box<dyn AlignmentResult>;
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

    fn global_alignment(&mut self, del: &i32, matrix: &Array2<i32>) -> Box<dyn AlignmentResult> {
        let mut alignment_matrix =
            Array2::<i32>::zeros((self.u8_sequence_2.len() + 1, self.u8_sequence_1.len() + 1));
        let mut direction_matrix = Array2::<Direction>::from_shape_fn(
            (self.u8_sequence_2.len() + 1, self.u8_sequence_1.len() + 1),
            |_| Direction::Beginning,
        );

        for x in 1..self.u8_sequence_1.len() + 1 {
            alignment_matrix[[0, x]] = -(x as i32) * del;
            direction_matrix[[0, x]] = Direction::Left
        }

        for y in 1..self.u8_sequence_2.len() + 1 {
            alignment_matrix[[y, 0]] = -(y as i32) * del;
            direction_matrix[[y, 0]] = Direction::Top;
        }

        alignment_matrix[[self.u8_sequence_2.len(), 0]] =
            -(self.u8_sequence_2.len() as i32 + 1) * del;
        alignment_matrix[[0, self.u8_sequence_1.len()]] =
            -(self.u8_sequence_1.len() as i32 + 1) * del;

        for (x, elem_1) in self.protein_sequence_1.iter().enumerate() {
            for (y, elem_2) in self.protein_sequence_2.iter().enumerate() {
                let x_real = x + 1;
                let y_real = y + 1;

                let seq_1_pos = *elem_1 as usize;
                let seq_2_pos = *elem_2 as usize;

                let top = alignment_matrix[[y_real - 1, x_real]] - del;
                let left = alignment_matrix[[y_real, x_real - 1]] - del;
                let diagonal =
                    alignment_matrix[[y_real - 1, x_real - 1]] + matrix[[seq_2_pos, seq_1_pos]];

                let max = max(max(top, left), diagonal);

                alignment_matrix[[y_real, x_real]] = max;

                if max == top {
                    direction_matrix[[y_real, x_real]] = Direction::Top
                } else if max == left {
                    direction_matrix[[y_real, x_real]] = Direction::Left
                } else if max == diagonal {
                    direction_matrix[[y_real, x_real]] = Direction::Diagonal
                }
            }
        }

        let mut current_x = self.protein_sequence_1.len() - 1;
        let mut current_y = self.protein_sequence_2.len() - 1;
        let (mut optimal_alignment_1, mut optimal_alignment_2) = (
            vec![*self.protein_sequence_1.last().unwrap()],
            vec![*self.protein_sequence_2.last().unwrap()],
        );

        loop {
            match direction_matrix[[current_y, current_x]] {
                Direction::Beginning => break,
                Direction::Top => {
                    optimal_alignment_1.push(Protein::Blank);
                    optimal_alignment_2.push(self.protein_sequence_2[current_y - 1]);
                    current_y -= 1;
                }
                Direction::Left => {
                    optimal_alignment_1.push(self.protein_sequence_1[current_x - 1]);
                    optimal_alignment_2.push(Protein::Blank);
                    current_x -= 1;
                }
                Direction::Diagonal => {
                    optimal_alignment_1.push(self.protein_sequence_1[current_x - 1]);
                    optimal_alignment_2.push(self.protein_sequence_2[current_y - 1]);
                    current_x -= 1;
                    current_y -= 1;
                }
            }
        }

        optimal_alignment_1.reverse();
        optimal_alignment_2.reverse();

        Box::new(GlobalAlignmentResult {
            alignment_matrix: alignment_matrix,
            direction_matrix: direction_matrix,
            optimal_alignment: (optimal_alignment_1, optimal_alignment_2),
        })
    }

    fn local_alignment(&mut self, del: &i32, matrix: &Array2<i32>) -> Box<dyn AlignmentResult> {
        let mut alignment_matrix =
            Array2::<i32>::zeros((self.u8_sequence_2.len() + 1, self.u8_sequence_1.len() + 1));
        let mut direction_matrix = Array2::<Direction>::from_shape_fn(
            (self.u8_sequence_2.len() + 1, self.u8_sequence_1.len() + 1),
            |_| Direction::Beginning,
        );

        let mut max_f: i32 = 0;
        let mut max_x: usize = 0;
        let mut max_y: usize = 0;

        for (x, elem_1) in self.protein_sequence_1.iter().enumerate() {
            for (y, elem_2) in self.protein_sequence_2.iter().enumerate() {
                let x_real = x + 1;
                let y_real = y + 1;

                let seq_1_pos = *elem_1 as usize;
                let seq_2_pos = *elem_2 as usize;

                let top = alignment_matrix[[y_real - 1, x_real]] - del;
                let left = alignment_matrix[[y_real, x_real - 1]] - del;
                let diagonal =
                    alignment_matrix[[y_real - 1, x_real - 1]] + matrix[[seq_2_pos, seq_1_pos]];

                let max = max(max(max(top, left), diagonal), 0);

                alignment_matrix[[y_real, x_real]] = max;

                if max == 0 {
                    direction_matrix[[y_real, x_real]] = Direction::Beginning
                } else if max == top {
                    direction_matrix[[y_real, x_real]] = Direction::Top
                } else if max == left {
                    direction_matrix[[y_real, x_real]] = Direction::Left
                } else if max == diagonal {
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
        let (mut optimal_alignment_1, mut optimal_alignment_2) = (
            vec![self.protein_sequence_1[max_x]],
            vec![self.protein_sequence_2[max_y]],
        );
        loop {
            match direction_matrix[[current_y, current_x]] {
                Direction::Beginning => break,
                Direction::Top => {
                    optimal_alignment_1.push(Protein::Blank);
                    optimal_alignment_2.push(self.protein_sequence_2[current_y - 1]);
                    current_y -= 1;
                }
                Direction::Left => {
                    optimal_alignment_1.push(self.protein_sequence_1[current_x - 1]);
                    optimal_alignment_2.push(Protein::Blank);
                    current_x -= 1;
                }
                Direction::Diagonal => {
                    optimal_alignment_1.push(self.protein_sequence_1[current_x - 1]);
                    optimal_alignment_2.push(self.protein_sequence_2[current_y - 1]);
                    current_x -= 1;
                    current_y -= 1;
                }
            }
        }

        optimal_alignment_1.reverse();
        optimal_alignment_2.reverse();

        Box::new(LocalAlignmentResult {
            alignment_matrix: alignment_matrix,
            direction_matrix: direction_matrix,
            max_f: max_f,
            optimal_alignment: (optimal_alignment_1, optimal_alignment_2),
        })
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
