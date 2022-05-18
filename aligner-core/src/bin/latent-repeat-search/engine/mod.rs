use crate::engine::task::Task;
use aligner_core::enums::{BioData, Index};
use aligner_core::Result;
use rand::{thread_rng, Rng};

pub(crate) mod calc;
pub(crate) mod sequences;
pub(crate) mod task;

#[allow(dead_code)]
pub(crate) enum MutationPercent {
    Quarter = 4,
    Half = 2,
}

pub(crate) fn generate_descendants<T: BioData + From<usize> + Clone>(
    sequence: &[T],
    amount: usize,
    percent: MutationPercent,
) -> Result<Vec<Vec<T>>> {
    let mut result = vec![];

    let offset = percent as usize;

    for i in 0..amount {
        result.push(mutate(sequence, offset, i).unwrap())
    }

    Ok(result)
}

fn mutate<T: BioData + From<usize> + Clone>(
    sequence: &[T],
    offset: usize,
    start: usize,
) -> Result<Vec<T>> {
    let mut result = Vec::from(sequence);

    let mut thread_rng = thread_rng();

    for i in (start..sequence.len()).step_by(offset) {
        result[i] = thread_rng.gen_range(0..T::volume()).into()
    }

    Ok(result)
}

pub(crate) fn filter(tasks: &[Task]) -> Result<Vec<Task>> {
    if tasks.is_empty() {
        return Ok(vec![]);
    }
    if tasks.len() == 1 {
        return Ok(tasks.to_vec());
    }

    let mut result = vec![];

    let mut buffer = Vec::from(tasks);

    while !buffer.is_empty() {
        if buffer.len() == 1 {
            result.push(buffer.first().unwrap().clone());
            break;
        }
        let current = buffer.first().unwrap();
        let mut batch = vec![current];

        let mut new = vec![];

        for task in buffer[1..].iter() {
            if check_intersection(
                (current.left_coord, current.right_coord),
                (task.left_coord, task.right_coord),
            ) {
                batch.push(task);
            } else {
                new.push(task.clone());
            }
        }

        result.push(
            (*batch
                .iter()
                .max_by(|u, v| u.z.partial_cmp(&v.z).unwrap())
                .unwrap())
            .clone(),
        );

        buffer = new.clone();
    }

    Ok(result)
}

fn check_intersection(coords_1: (usize, usize), coords_2: (usize, usize)) -> bool {
    if coords_1.0 >= coords_2.0 && coords_1.0 <= coords_2.1 {
        return true;
    }
    if coords_1.1 >= coords_2.0 && coords_1.1 <= coords_2.1 {
        return true;
    }
    if coords_2.0 >= coords_1.0 && coords_2.1 <= coords_1.1 {
        return true;
    }
    if coords_1.0 >= coords_2.0 && coords_1.0 <= coords_2.1 {
        return true;
    }

    false
}

pub(crate) fn index_coord(target: usize, indices: &[Index]) -> usize {
    for index in indices.iter() {
        if target >= index.coord {
            return target + index.offset;
        }
    }

    target
}

pub(crate) fn rotate_indices(indices: &[Index], query_length: usize) -> Vec<Index> {
    let mut result = vec![];

    let _ref = match indices.first() {
        Some(_ref) => _ref,
        None => return result,
    };
    let full_length = query_length + _ref.offset;
    let mut offset = 0;

    for index in indices.iter() {
        offset += index.local_offset;
        result.push(Index {
            coord: full_length - index.coord - _ref.offset,
            offset,
            local_offset: index.local_offset,
        });
    }

    result.reverse();
    result
}
