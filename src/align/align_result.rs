use crate::align::aligner_core::AlignmentResult;
use crate::align::enums::{Direction, Protein};
use ndarray::Array2;
use std::str::from_utf8;

pub struct GlobalAlignmentResult {
    pub alignment_matrix: Array2<f64>,
    pub direction_matrix: Array2<Direction>,
    pub max_f: f64,
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

    fn get_alignment_matrix(&self) -> &Array2<f64> {
        &self.alignment_matrix
    }

    fn get_direction_matrix(&self) -> &Array2<Direction> {
        &self.direction_matrix
    }

    fn get_optimal_alignment(&self) -> (&Vec<Protein>, &Vec<Protein>) {
        (&self.optimal_alignment.0, &self.optimal_alignment.1)
    }

    fn get_frequency_matrix(&self) -> Array2<f64> {
        let mut frequency_matrix = Array2::<f64>::zeros((25, 25));

        for (x, y) in self
            .optimal_alignment
            .0
            .iter()
            .zip(self.optimal_alignment.1.iter())
        {
            if (*x != Protein::Del || *x != Protein::Ins) && (*y != Protein::Del || *y != Protein::Ins) {
                frequency_matrix[[*y as usize, *x as usize]] += 1f64;
            }
        }

        frequency_matrix
    }

    fn get_score(&self) -> f64 {
        self.max_f
    }
}

#[derive(Debug)]
pub struct LocalAlignmentResult {
    pub alignment_matrix: Array2<f64>,
    pub direction_matrix: Array2<Direction>,
    #[allow(dead_code)]
    pub max_f: f64,
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

    fn get_alignment_matrix(&self) -> &Array2<f64> {
        &self.alignment_matrix
    }

    fn get_direction_matrix(&self) -> &Array2<Direction> {
        &self.direction_matrix
    }

    fn get_optimal_alignment(&self) -> (&Vec<Protein>, &Vec<Protein>) {
        (&self.optimal_alignment.0, &self.optimal_alignment.1)
    }

    fn get_frequency_matrix(&self) -> Array2<f64> {
        let mut frequency_matrix = Array2::<f64>::zeros((25, 25));

        for (x, y) in self
            .optimal_alignment
            .0
            .iter()
            .zip(self.optimal_alignment.1.iter())
        {
            if (*x != Protein::Del || *x != Protein::Ins) && (*y != Protein::Del || *y != Protein::Ins) {
                frequency_matrix[[*y as usize, *x as usize]] += 1f64;
            }
        }

        frequency_matrix
    }

    fn get_score(&self) -> f64 {
        self.max_f
    }
}
