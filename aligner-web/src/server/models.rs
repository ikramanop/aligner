use ndarray::{Array1, Array2};
use serde::{Deserialize, Serialize};
use std::collections::hash_map::DefaultHasher;
use std::collections::HashMap;
use std::hash::{Hash, Hasher};

#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct HealthCheck {
    pub nodes: Vec<HealthCheckUnit>,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct HealthCheckUnit {
    pub consumer_name: String,
    pub status: bool,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct AlignJobRequest {
    pub sequences: String,
    pub kd_value: f64,
    pub r_squared_value: f64,
    pub del_value: f64,
    pub dim_value: i32,
    pub matrices_volume_value: i32,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct AlignJob {
    pub sequence_1: Vec<u8>,
    pub sequence_2: Vec<u8>,
    pub matrix: Option<Array2<f64>>,
    pub frequences: Array1<f64>,
    pub kd_value: f64,
    pub r_squared_value: f64,
    pub del_value: f64,
    pub matrices_volume_value: i32,
    pub hash: String,
}

impl AlignJobRequest {
    pub fn to_hashing_struct(
        &self,
        query_sequence: String,
        target_sequence: String,
    ) -> HashingStruct {
        HashingStruct {
            query_sequence,
            target_sequence,
            kd_value: format!("{:.5}", self.kd_value),
            r_squared_value: format!("{:.5}", self.r_squared_value),
            del_value: format!("{:.5}", self.del_value),
            dim_value: self.dim_value,
            matrices_volume_value: self.matrices_volume_value,
        }
    }
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct AlignJobResult {
    pub matrix: Array2<f64>,
    pub max_f: f64,
    pub matrices_volume_value: i32,
    pub hash: String,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct ProgressEventResponse {
    pub progress: HashMap<String, f64>,
    pub message: String,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct AlignmentResultEventResponse {
    pub progress: f64,
    pub matrix: Array2<f64>,
    pub max_f: f64,
}

#[derive(Debug, Clone)]
pub struct AlignTask {
    pub query_sequence: String,
    pub target_sequence: String,
    pub kd_value: f64,
    pub r_squared_value: f64,
    pub del_value: f64,
    pub dim_value: i32,
    pub matrices_volume_value: i32,
    pub status: String,
    pub result_matrix: Option<Array2<f64>>,
    pub f_value: Option<f64>,
    pub result_query_sequence: Option<String>,
    pub result_target_sequence: Option<String>,
    pub p_value: Option<String>,
}

#[derive(Debug, Clone, Hash)]
pub struct HashingStruct {
    pub query_sequence: String,
    pub target_sequence: String,
    pub kd_value: String,
    pub r_squared_value: String,
    pub del_value: String,
    pub dim_value: i32,
    pub matrices_volume_value: i32,
}

impl HashingStruct {
    pub fn calculate_hash(&self) -> String {
        let mut s = DefaultHasher::new();
        self.hash(&mut s);
        s.finish().to_string()
    }
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct EmptySuccessfulResponse {}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct EmptySuccessfulResponseWithHashes {
    pub hashes: Vec<String>,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct ErroneousResponse {
    pub message: String,
}
