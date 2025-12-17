use axum::Router;

pub mod routes;
pub mod auth;
pub mod types;

pub fn build_router() -> Router {
    Router::new()
        .merge(routes::routes())
        .layer(auth::AuthLayer)
}
