use aligner::{database::get_connection, statistics::calculate_p_value};
use rdkafka::message::ToBytes;

extern crate pretty_env_logger;
#[macro_use]
extern crate log;

fn main() {
    pretty_env_logger::init();

    let mut conn = get_connection(false);

    let ids = conn.get_ids_with_null_p_value().unwrap();

    for id in ids.iter() {
        debug!("Calculating p-value for task with id {}", id);

        let subtask = conn.get_result_matrix_by_task_id(*id).unwrap();

        let p_value = calculate_p_value(
            subtask.query_sequence.to_bytes(),
            subtask.target_sequence.to_bytes(),
            subtask.f_value,
            subtask.del_value,
            subtask.del_value,
            &subtask.matrix,
        )
        .unwrap();

        debug!("Calculated p-value {} for task with id {}", p_value, id);

        conn.add_p_value_by_id(p_value, *id).unwrap();
    }
}
