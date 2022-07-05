extern crate log;
extern crate ndarray;
extern crate pretty_env_logger;

use std::result;

pub mod csv;
pub mod files;
pub mod matrices;

#[derive(Debug)]
pub enum Error {
    WrongMatrixSpecified,
}

pub type Result<T> = result::Result<T, Error>;
