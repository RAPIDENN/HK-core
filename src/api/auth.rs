use axum::body::Body;
use axum::{
    http::{Request, StatusCode},
    response::Response,
};
use futures::future::BoxFuture;
use tower_layer::Layer;
use tower_service::Service;

#[derive(Clone)]
pub struct AuthLayer;

#[derive(Clone)]
pub struct AuthMiddleware<S> {
    inner: S,
    token: String,
}

impl<S> Layer<S> for AuthLayer {
    type Service = AuthMiddleware<S>;

    fn layer(&self, inner: S) -> Self::Service {
        AuthMiddleware {
            inner,
            token: std::env::var("AUTH_TOKEN").unwrap_or_default(),
        }
    }
}

impl<S, ReqBody> Service<Request<ReqBody>> for AuthMiddleware<S>
where
    S: Service<Request<ReqBody>, Response = Response> + Clone + Send + 'static,
    S::Future: Send + 'static,
    ReqBody: Send + 'static,
{
    type Response = Response;
    type Error = S::Error;
    type Future = BoxFuture<'static, Result<Response, Self::Error>>;

    fn poll_ready(
        &mut self,
        cx: &mut std::task::Context<'_>,
    ) -> std::task::Poll<Result<(), Self::Error>> {
        self.inner.poll_ready(cx)
    }

    fn call(&mut self, req: Request<ReqBody>) -> Self::Future {
        let token = self.token.clone();
        let mut svc = self.inner.clone();
        Box::pin(async move {
            if authorized(&req, &token) {
                svc.call(req).await
            } else {
                let mut resp = Response::new(Body::from("Unauthorized"));
                *resp.status_mut() = StatusCode::UNAUTHORIZED;
                Ok(resp)
            }
        })
    }
}

fn authorized<B>(req: &Request<B>, token: &str) -> bool {
    if token.is_empty() {
        return false;
    }
    let Some(header) = req.headers().get(axum::http::header::AUTHORIZATION) else {
        return false;
    };
    match header.to_str() {
        Ok(value) => value == format!("Bearer {}", token),
        Err(_) => false,
    }
}
