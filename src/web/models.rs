use ndarray::{Array1, Array2};
use serde::{Deserialize, Serialize};
use uuid::Uuid;

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
    pub uuid: Uuid,
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
    pub uuid: Uuid,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct AlignJobResult {
    pub matrix: Array2<f64>,
    pub max_f: f64,
    pub matrices_volume_value: i32,
    pub uuid: Uuid,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct ProgressEventResponse {
    pub progress: f64,
    pub message: String,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct AlignmentResultEventResponse {
    pub progress: f64,
    pub matrix: Array2<f64>,
    pub max_f: f64,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct EmptySuccessfulResponse {}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct ErroneousResponse {
    pub message: String,
}
