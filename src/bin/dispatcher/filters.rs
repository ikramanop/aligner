use crate::handlers;
use warp::{http::Method, Filter};

// API combined filter
pub fn get_api() -> impl Filter<Extract = impl warp::Reply, Error = warp::Rejection> + Clone {
    check_health().or(validate()).or(progress())
}

// GET health/check with no params
pub fn check_health() -> impl Filter<Extract = impl warp::Reply, Error = warp::Rejection> + Clone {
    let cors_health_check = warp::cors()
        .allow_any_origin()
        .allow_header("content-type")
        .allow_methods(&[Method::GET]);

    warp::get()
        .and(warp::path("health"))
        .and(warp::path("check"))
        .and_then(handlers::check_health)
        .with(cors_health_check)
}

// POST /validate with json body
pub fn validate() -> impl Filter<Extract = impl warp::Reply, Error = warp::Rejection> + Clone {
    let cors_validate = warp::cors()
        .allow_any_origin()
        .allow_header("content-type")
        .allow_methods(&[Method::POST]);

    warp::post()
        .and(warp::path("validate"))
        .and(warp::body::json())
        .and_then(handlers::validate)
        .with(cors_validate)
}

// GET {uuid}/progress for sse update process events
pub fn progress() -> impl Filter<Extract = impl warp::Reply, Error = warp::Rejection> + Clone {
    let cors_progress = warp::cors()
        .allow_any_origin()
        .allow_header("content-type")
        .allow_methods(&[Method::GET]);

    warp::get()
        .and(warp::path::param())
        .and(warp::path("progress"))
        .and_then(handlers::progress)
        .with(cors_progress)
}
