use crate::Error::WrongMatrixSpecified;
use crate::Result;
use ndarray::Axis;
use ndarray::{Array1, Array2};
use ndarray_stats::DeviationExt;
use roots::{find_roots_quadratic, Roots};

pub fn get_threshold(dim_1: usize) -> f64 {
    match dim_1 {
        20 => 22.6,
        21 => 23.1,
        22 => 23.6,
        23 => 24.1,
        24 => 24.6,
        _ => 0 as f64,
    }
}

pub fn transform_matrix(
    matrix: &Array2<f64>,
    k_d: f64,
    r_squared: f64,
    frequencies: &Array1<f64>,
) -> Result<Array2<f64>> {
    let f = Array1::<f64>::from_shape_fn(matrix.len_of(Axis(1)), |_| {
        1f64 / matrix.len_of(Axis(1)) as f64
    });

    let p = cross_product(frequencies, &f);

    let p_squared = p.mapv(|a| a * a).sum();

    let k_0 = (&p * matrix).sum();

    let a = (k_d - k_0) / p_squared;
    let b = k_d / p_squared;
    let difference = a - b;

    let denominator = (matrix + &p * difference).mapv(|a| a * a).sum();

    let a_coeff = (2f64 * b * (&p * (matrix + &p * difference)).sum()) / denominator;
    let b_coeff = (b * b * p_squared - r_squared) / denominator;

    match find_roots_quadratic(1f64, a_coeff, b_coeff) {
        Roots::No(_) => Err(WrongMatrixSpecified),
        Roots::One(roots) => Ok(&p * b + roots[0] * (matrix + &p * (a - b))),
        Roots::Two(roots) => {
            if roots[0] > 0f64 && roots[1] < 0f64 {
                Ok(&p * b + roots[0] * (matrix + &p * (a - b)))
            } else if roots[0] < 0f64 && roots[1] > 0f64 {
                Ok(&p * b + roots[1] * (matrix + &p * (a - b)))
            } else {
                let new_matrix_1 = &p * b + roots[0] * (matrix + &p * (a - b));
                let new_matrix_2 = &p * b + roots[1] * (matrix + &p * (a - b));

                let distance_1 = matrix.l2_dist(&new_matrix_1).unwrap();
                let distance_2 = matrix.l2_dist(&new_matrix_2).unwrap();

                if distance_1 < distance_2 {
                    Ok(new_matrix_1)
                } else {
                    Ok(new_matrix_2)
                }
            }
        }
        _ => panic!("This should not happen"),
    }
}

fn cross_product(array1: &Array1<f64>, array2: &Array1<f64>) -> Array2<f64> {
    let mut result = Array2::<f64>::zeros((array1.len(), array2.len()));

    for (y, elem_1) in array1.iter().enumerate() {
        for (x, elem_2) in array2.iter().enumerate() {
            result[[y, x]] = elem_1 * elem_2;
        }
    }

    result
}
