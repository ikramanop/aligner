mod filters;
mod handlers;
use filters::get_api;
use std::env;
use warp::Filter;

extern crate pretty_env_logger;
#[macro_use]
extern crate log;

fn load_config(path: &'_ str) {
    debug!("Loading config from {:?}.", path);
    dotenv::from_filename(path).unwrap();
    debug!("Config successfully loaded.")
}

#[tokio::main]
async fn main() {
    pretty_env_logger::init();

    match env::var("CONFIG_PATH") {
        Ok(value) => load_config(&value),
        Err(err) => panic!("Unable to load config: {}", err),
    }

    let routes = get_api().with(warp::log("aligner"));

    warp::serve(routes).run(([127, 0, 0, 1], 3030)).await;
}
