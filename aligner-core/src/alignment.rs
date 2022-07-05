use crate::{AlignmentTrait, BioData};
use ndarray::Array2;

#[derive(Debug, Clone)]
pub struct Alignment<T: BioData> {
    pub(crate) query: Vec<T>,
    pub(crate) target: Vec<T>,
    pub coords: ((usize, usize), (usize, usize)),
    pub f: f64,
}

impl<T: BioData + Into<usize> + Eq + Copy> AlignmentTrait<T> for Alignment<T> {
    fn get_frequency_matrix(&self) -> Array2<f64> {
        let mut frequency_matrix = Array2::<f64>::zeros((T::volume(), T::volume()));

        for (x, y) in self.query.iter().zip(self.target.iter()) {
            if *x != T::blank() && *y != T::blank() {
                frequency_matrix[[(*y).into(), (*x).into()]] += 1f64;
            }
        }

        frequency_matrix
    }

    fn get_alignment(&self, matrix: &Array2<f64>) -> Vec<T> {
        let mut alignment = Vec::<T>::new();

        for (x, y) in self.query.iter().zip(self.target.iter()) {
            if *x == *y {
                alignment.push(*x);
            } else if (*x != T::blank())
                && (*y != T::blank())
                && matrix[[(*y).into(), (*x).into()]] >= 0f64
            {
                alignment.push(T::pos())
            } else {
                alignment.push(T::blank())
            }
        }

        alignment
    }
}

#[derive(Debug, Clone)]
pub struct PWMAlignment<T: BioData> {
    pub(crate) numbered: Vec<usize>,
    pub(crate) query: Vec<T>,
    pub(crate) dim: usize,
    pub coords: ((usize, usize), (usize, usize)),
    pub f: f64,
}

impl<T: BioData + Into<usize> + Eq + Copy> AlignmentTrait<T> for PWMAlignment<T> {
    fn get_frequency_matrix(&self) -> Array2<f64> {
        let mut frequency_matrix = Array2::<f64>::zeros((T::volume(), self.dim));

        for (x, y) in self.numbered.iter().zip(self.query.iter()) {
            if *x != 0 && *y != T::blank() {
                frequency_matrix[[(*y).into(), *x - 1]] += 1f64;
            }
        }

        frequency_matrix
    }

    fn get_alignment(&self, _matrix: &Array2<f64>) -> Vec<T> {
        let mut alignment = Vec::<T>::new();

        for (x, y) in self.numbered.iter().zip(self.query.iter()) {
            if *x != 0 {
                alignment.push(*y)
            } else {
                alignment.push(T::blank())
            }
        }

        alignment
    }
}

impl<T: BioData + Into<usize> + Eq + Copy> PWMAlignment<T> {
    pub fn empty() -> Self {
        PWMAlignment {
            numbered: vec![],
            query: vec![],
            dim: 0,
            coords: ((0, 0), (0, 0)),
            f: 0f64,
        }
    }
}
