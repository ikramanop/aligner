use crate::cmd::csv::run_csv_cmd;
use crate::cmd::exploring::run_exploring_cmd;
use crate::cmd::testing::run_testing_cmd;
use crate::engine::task::Task;
use crate::error::Error;
use crate::Args;
use ndarray::Array2;
use std::collections::HashMap;
use std::path::PathBuf;

pub(crate) mod csv;
pub(crate) mod exploring;
pub(crate) mod testing;

const DIRECT: &str = "direct";
const INVERSE: &str = "inverse";

type CMDResult = Result<HashMap<String, (Vec<Task>, Array2<f64>)>, Error>;

pub(crate) struct CMDOptions {
    pub(crate) repeat_length: usize,

    pub(crate) query_offset: usize,

    pub(crate) deletions: f64,

    pub(crate) extension: f64,

    pub(crate) rsquared: f64,

    pub(crate) kd: f64,

    pub(crate) threads: usize,

    pub(crate) repeats: usize,

    pub(crate) simple_init: bool,

    pub(crate) reverse: bool,

    pub(crate) testing: bool,

    pub(crate) csv: bool,

    pub(crate) fasta_path: Option<PathBuf>,

    pub(crate) csv_path: Option<PathBuf>,
}

impl CMDOptions {
    pub(crate) fn from_args(args: &Args) -> Self {
        let mut testing = false;
        let mut csv = false;

        let fasta_path = match args.input.clone() {
            Some(fasta_path) => Some(PathBuf::from(fasta_path)),
            None => {
                testing = true;
                None
            }
        };

        let csv_path = match args.csv.clone() {
            Some(csv_path) => {
                csv = true;
                Some(PathBuf::from(csv_path))
            }
            None => None,
        };

        CMDOptions {
            repeat_length: args.repeat_length,
            query_offset: args.query_offset,
            deletions: args.deletions,
            extension: args.extension,
            rsquared: args.rsquared,
            kd: args.kd,
            threads: args.threads,
            repeats: args.repeats,
            simple_init: args.simple_init,
            reverse: args.reverse,
            testing,
            csv,
            fasta_path,
            csv_path,
        }
    }
}

pub(crate) fn run_root_cmd(opts: &CMDOptions) -> CMDResult {
    if opts.testing {
        return run_testing_cmd(opts);
    }
    if opts.csv {
        return run_csv_cmd(opts);
    }
    run_exploring_cmd(opts)
}
