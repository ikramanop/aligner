use crate::align::enums::Protein;
use crate::web::models::AlignJobRequest;
use core::panic;
use mysql::{error::MySqlError, prelude::*};
use mysql::{Conn, Opts, Result};
use ndarray::Array2;
use std::env;
use std::str::FromStr;

mod models;
mod queries;

fn init_database(conn: &mut Conn) -> Result<()> {
    for query in queries::INIT_QUERIES.iter() {
        match conn.query_drop(query) {
            Ok(_) => continue,
            Err(err) => return Err(err),
        }
    }

    Ok(())
}

pub struct Connector {
    conn: Conn,
}

impl Connector {
    pub fn get_base_matrices_with_limit(
        &mut self,
        dim: usize,
        limit: usize,
    ) -> Result<Vec<String>> {
        self.conn.exec_map(
            queries::GET_BASE_MATRICES_WITH_LIMIT,
            (dim, limit),
            |matrix| matrix,
        )
    }

    pub fn insert_base_matrix(&mut self, dim: usize, matrix: &Array2<f64>) -> Result<()> {
        match self.conn.exec_drop(
            queries::INSERT_BASE_MATRIX,
            (dim, serde_json::to_string(matrix).unwrap()),
        ) {
            Ok(_) => Ok(()),
            Err(err) => Err(err),
        }
    }

    pub fn insert_align_task(
        &mut self,
        job: &AlignJobRequest,
        query_sequence_id: String,
        query_sequence: String,
        target_sequence_id: String,
        target_sequence: String,
    ) -> Result<String> {
        let hash = job
            .to_hashing_struct(query_sequence.clone(), target_sequence.clone())
            .calculate_hash();

        match self.conn.exec_drop(
            queries::INSERT_ALIGN_TASK,
            (
                &hash,
                &query_sequence_id,
                &query_sequence,
                &target_sequence_id,
                &target_sequence,
                job.kd_value,
                job.r_squared_value,
                job.del_value,
                job.dim_value,
                job.matrices_volume_value,
                "active",
            ),
        ) {
            Ok(_) => Ok(hash),
            Err(err) => Err(err),
        }
    }

    pub fn get_align_task_id_by_hash(&mut self, hash: String) -> Result<i32> {
        match self
            .conn
            .exec_first(queries::GET_ALIGN_TASK_ID_BY_HASH, (hash,))
        {
            Ok(id) => Ok(id.unwrap()),
            Err(err) => Err(err),
        }
    }

    pub fn insert_align_subtask(
        &mut self,
        hash: String,
        f_value: f64,
        matrix: &Array2<f64>,
        sequences: (Vec<Protein>, Vec<Protein>),
    ) -> Result<()> {
        let id = match self.get_align_task_id_by_hash(hash) {
            Ok(id) => id,
            Err(err) => return Err(err),
        };

        match self.conn.exec_drop(
            queries::INSERT_ALIGN_SUBTASKS,
            (
                id,
                f_value,
                serde_json::to_string(matrix).unwrap(),
                String::from_utf8(Protein::protein_vec_to_u8_vec(&sequences.0)).unwrap(),
                String::from_utf8(Protein::protein_vec_to_u8_vec(&sequences.1)).unwrap(),
            ),
        ) {
            Ok(_) => Ok(()),
            Err(err) => Err(err),
        }
    }

    pub fn get_percentage_by_hash(&mut self, hash: String) -> Result<f64> {
        match self
            .conn
            .exec_first(queries::GET_PERCENTAGE_BY_HASH, (hash,))
        {
            Ok(percentage) => Ok(percentage.unwrap_or(0f64)),
            Err(err) => Err(err),
        }
    }

    pub fn get_subtask_with_max_f_value_by_hash(
        &mut self,
        hash: String,
    ) -> Result<models::AlignSubtask> {
        let row: (f64, String, String, String) = match self
            .conn
            .exec_first(queries::GET_SUBTASK_WITH_MAX_F_VALUE_BY_HASH, (hash,))
        {
            Ok(row) => match row {
                Some(row_tuple) => row_tuple,
                None => {
                    return Err(mysql::Error::MySqlError(MySqlError {
                        state: String::from_str("error").unwrap(),
                        message: String::from_str("empty row").unwrap(),
                        code: 228,
                    }))
                }
            },
            Err(err) => return Err(err),
        };

        Ok(models::AlignSubtask {
            f_value: row.0,
            matrix_json: row.1,
            result_query_sequence: row.2,
            result_target_sequence: row.3,
        })
    }

    pub fn insert_result_matrix_by_hash(
        &mut self,
        subtask: models::AlignSubtask,
        hash: String,
    ) -> Result<()> {
        let id = match self.get_align_task_id_by_hash(hash) {
            Ok(id) => id,
            Err(err) => return Err(err),
        };

        match self.conn.exec_drop(
            queries::INSERT_RESULT_MATRIX,
            (
                id,
                subtask.f_value,
                subtask.matrix_json,
                subtask.result_query_sequence,
                subtask.result_target_sequence,
            ),
        ) {
            Ok(_) => Ok(()),
            Err(err) => Err(err),
        }
    }

    pub fn delete_subtasks_by_hash(&mut self, hash: String) -> Result<()> {
        match self
            .conn
            .exec_drop(queries::DELETE_SUBTASKS_BY_HASH, (hash,))
        {
            Ok(_) => Ok(()),
            Err(err) => Err(err),
        }
    }

    pub fn get_ids_with_null_p_value(&mut self) -> Result<Vec<i32>> {
        match self.conn.query(queries::GET_IDS_WITH_NULL_P_VALUE) {
            Ok(result) => Ok(result),
            Err(err) => Err(err),
        }
    }

    pub fn get_result_matrix_by_task_id(&mut self, id: i32) -> Result<models::AlignTaskWithMatrix> {
        let row: (String, String, f64, f64, String) = match self
            .conn
            .exec_first(queries::GET_RESULT_MATRIX_BY_TASK_ID, (id,))
        {
            Ok(row) => match row {
                Some(row_tuple) => row_tuple,
                None => {
                    return Err(mysql::Error::MySqlError(MySqlError {
                        state: String::from_str("error").unwrap(),
                        message: String::from_str("empty row").unwrap(),
                        code: 228,
                    }))
                }
            },
            Err(err) => return Err(err),
        };

        Ok(models::AlignTaskWithMatrix {
            query_sequence: row.0,
            target_sequence: row.1,
            f_value: row.2,
            del_value: row.3,
            matrix: serde_json::from_str(&row.4).unwrap(),
        })
    }

    pub fn add_p_value_by_id(&mut self, p_value: f64, id: i32) -> Result<()> {
        match self
            .conn
            .exec_drop(queries::ADD_P_VALUE_BY_ID, (p_value, id))
        {
            Ok(_) => Ok(()),
            Err(err) => Err(err),
        }
    }

    pub fn get_all_hashes(&mut self) -> Result<Vec<String>> {
        self.conn.query_map(
            queries::GET_ALL_HASHES,
            |hash| hash,
        )
    }
}

pub fn get_connection(init: bool) -> Connector {
    let host = match env::var("MARIADB_HOST") {
        Ok(host) => host,
        Err(err) => panic!(
            "Failed with {}. Env value MARIADB_HOST might not be set.",
            err
        ),
    };
    let port = match env::var("MARIADB_PORT") {
        Ok(port) => port,
        Err(err) => panic!(
            "Failed with {}. Env value MARIADB_PORT might not be set.",
            err
        ),
    };
    let database = match env::var("MARIADB_DATABASE") {
        Ok(database) => database,
        Err(err) => panic!(
            "Failed with {}. Env value MARIADB_DATABASE might not be set.",
            err
        ),
    };
    let user = match env::var("MARIADB_USER") {
        Ok(user) => user,
        Err(err) => panic!(
            "Failed with {}. Env value MARIADB_USER might not be set.",
            err
        ),
    };
    let password = match env::var("MARIADB_PASSWORD") {
        Ok(password) => password,
        Err(err) => panic!(
            "Failed with {}. Env value MARIADB_PASSWORD might not be set.",
            err
        ),
    };

    let opts = Opts::from_url(
        format!(
            "mysql://{}:{}@{}:{}/{}",
            user, password, host, port, database
        )
        .as_str(),
    )
    .unwrap();

    let mut conn = match Conn::new(opts) {
        Ok(conn) => conn,
        Err(err) => panic!("Failed with {}. Unable to connect to mariadb", err),
    };

    if init {
        init_database(&mut conn).unwrap()
    }

    Connector { conn }
}
