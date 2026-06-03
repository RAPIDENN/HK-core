use crate::engine::mill::MillRefineOutput;
use crate::engine::mill::MillRunOutput;
use serde::Serialize;

#[derive(Serialize)]
pub struct ErrorResponse {
    pub error: String,
}

#[derive(Serialize)]
pub struct MillRunResponse {
    pub result: MillRunOutput,
}

#[derive(Serialize)]
pub struct MillRefineResponse {
    pub result: MillRefineOutput,
}
