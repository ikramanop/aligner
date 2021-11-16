extern crate pretty_env_logger;
#[macro_use] extern crate log;

pub mod align;
pub mod files;
pub mod matrices;
pub mod statistics;
pub mod web;

#[cfg(test)]
mod tests;
