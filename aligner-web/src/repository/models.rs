use ndarray::Array2;

#[derive(Debug, Clone)]
pub struct AlignSubtask {
    pub f_value: f64,
    pub matrix_json: String,
    pub result_query_sequence: String,
    pub result_target_sequence: String,
}

#[derive(Debug, Clone)]
pub struct AlignTaskWithMatrix {
    pub query_sequence: String,
    pub target_sequence: String,
    pub f_value: f64,
    pub del_value: f64,
    pub matrix: Array2<f64>,
}
