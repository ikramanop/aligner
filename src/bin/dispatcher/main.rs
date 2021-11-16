mod filters;
mod handlers;

use std::env;

use filters::get_api;
use warp::Filter;

fn load_config(path: &'_ str) {
    println!("{:?}", path);
    dotenv::from_filename(path).unwrap();
}

#[tokio::main]
async fn main() {
    match env::var("CONFIG_PATH") {
        Ok(value) => load_config(&value),
        Err(err) => panic!("Unable to load config: {}", err),
    }

    println!("{:?}", env::var("SERVER_ENV").unwrap());

    pretty_env_logger::init();

    let routes = get_api().with(warp::log("aligner"));

    warp::serve(routes).run(([127, 0, 0, 1], 3030)).await;
}
