/*
 * HK-core
 * Copyright (C) 2025 RAPIDENN
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 */

mod engine;
mod api;

use tracing_subscriber::EnvFilter;

#[tokio::main]
async fn main() {
    tracing_subscriber::fmt().with_env_filter(EnvFilter::from_default_env()).init();
    let app = api::build_router();
    let listener = tokio::net::TcpListener::bind(("0.0.0.0", 8080)).await.unwrap();
    tracing::info!("listening on 0.0.0.0:8080");
    axum::serve(listener, app).await.unwrap();
}
