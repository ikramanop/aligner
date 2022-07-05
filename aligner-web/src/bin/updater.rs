use aligner::repository::get_connection;

extern crate pretty_env_logger;
#[macro_use]
extern crate log;

fn main() {
    pretty_env_logger::init();

    let mut conn = get_connection(false);

    let hashes = conn.get_all_hashes().unwrap();

    for hash in hashes.iter() {
        let percentage = conn.get_percentage_by_hash(hash.clone()).unwrap();

        if (percentage - 100f64).abs() < f64::EPSILON {
            debug!("Dropping subtasks for {}", hash);

            let subtask = conn
                .get_subtask_with_max_f_value_by_hash(hash.clone())
                .unwrap();
            conn.insert_result_matrix_by_hash(subtask, hash.clone())
                .unwrap();
            conn.delete_subtasks_by_hash(hash.clone()).unwrap();
        } else {
            debug!("Calculated {}% for task with hash {}.", percentage, hash);
        }
    }
}
