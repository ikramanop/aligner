use crate::enums::Direction;
use crate::{AlignmentTrait, BioData};
use ndarray::Array2;
use std::marker::PhantomData;

#[derive(Debug, Clone)]
pub struct AlignmentResult<T: BioData + Into<usize> + Copy + Eq, A: AlignmentTrait<T>> {
    pub alignment_matrix: Array2<f64>,
    pub direction_matrix: Array2<Direction>,
    pub alignment: A,
    pub matrix: Option<Array2<f64>>,
    pub(crate) phantom: PhantomData<T>,
}
