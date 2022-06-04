use crate::args::Args;
use crate::cmd::{run_root_cmd, CMDOptions};
use crate::engine::filter;
use crate::error::Error;
use aligner_helpers::csv::{CsvInput, Record};
use clap::Parser;
use ndarray::Array2;
use std::collections::HashMap;
use std::path::PathBuf;

extern crate pretty_env_logger;
#[macro_use]
extern crate log;

mod args;
mod cmd;
mod engine;
mod error;

fn main() -> Result<(), Error> {
    pretty_env_logger::init();

    let args = Args::parse();

    let output_path = match args.output.clone() {
        Some(output) => PathBuf::from(output),
        None => std::env::current_dir().unwrap().join("output.csv"),
    };

    let matrices_output_path = match args.output.clone() {
        Some(output) => PathBuf::from(&format!("{}.matrices.json", output)),
        None => std::env::current_dir().unwrap().join("matrices.json"),
    };

    let opts = CMDOptions::from_args(&args);

    let result = match run_root_cmd(&opts) {
        Ok(result) => result,
        Err(err) => return Err(err),
    };

    let mut wtr = CsvInput::new(output_path.as_path()).unwrap();

    let mut matrices = HashMap::<String, Array2<f64>>::new();

    for (key, value) in result.iter() {
        for task in value.0.iter() {
            wtr.write(&Record {
                name: key.to_owned(),
                z_value: task.z,
                left_coord: task.left_coord,
                right_coord: task.right_coord,
            })
            .unwrap();
        }

        matrices.insert(key.clone(), value.1.clone());
    }

    std::fs::write(
        matrices_output_path.as_path(),
        serde_json::to_string(&matrices).unwrap(),
    )
    .unwrap();

    println!(
        "\nOutput written to:\n 1. Result: {}\n 2. Matrices: {}",
        output_path.into_os_string().into_string().unwrap(),
        matrices_output_path.into_os_string().into_string().unwrap()
    );

    Ok(())
}
