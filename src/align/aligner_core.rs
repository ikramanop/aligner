use crate::align::enums::{Direction, Protein};
use ndarray::Array2;

pub trait AlignmentResult {
    fn represent(&self);
    fn get_alignment_matrix(&self) -> &Array2<f64>;
    fn get_direction_matrix(&self) -> &Array2<Direction>;
    fn get_optimal_alignment(&self) -> (&Vec<Protein>, &Vec<Protein>);
    fn get_frequency_matrix(&self) -> Array2<f64>;
    fn get_score(&self) -> f64;
}

pub trait Aligner {
    fn local_alignment(&mut self, del: f64, ins: f64, matrix: &Array2<f64>) -> Box<dyn AlignmentResult>;
    fn global_alignment(&mut self, del: f64, ins: f64, matrix: &Array2<f64>) -> Box<dyn AlignmentResult>;
    fn get_symbolic(&mut self) -> (String, String);
}

