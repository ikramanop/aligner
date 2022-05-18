use aligner_core::alignment::PWMAlignment;
use aligner_core::enums::DNA;

#[derive(Debug, Clone)]
pub(crate) struct Task {
    pub(crate) alignment: PWMAlignment<DNA>,
    pub(crate) left_coord: usize,
    pub(crate) right_coord: usize,
    pub(crate) z: f64,
}