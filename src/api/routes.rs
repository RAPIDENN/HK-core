use super::types::{ErrorResponse, MillRefineResponse, MillRunResponse};
use axum::{
    extract::State,
    http::StatusCode,
    routing::{get, post},
    Json, Router,
};
use serde::Deserialize;

#[derive(Clone)]
pub struct AppState;

pub fn routes() -> Router {
    Router::new()
        .route("/health", get(health))
        .route("/mill/run", post(mill_run))
        .route("/mill/refine", post(mill_refine))
        .with_state(AppState)
}

async fn health() -> Json<serde_json::Value> {
    Json(serde_json::json!({"status": "ok"}))
}

#[derive(Deserialize)]
struct MillRunRequest {
    l: usize,
    beta: f64,
    n_thermal_sweeps: usize,
    n_sweeps: usize,
    measure_every: usize,
    step_size: f64,
    seed: u64,
}

async fn mill_run(
    State(_): State<AppState>,
    axum::Json(payload): axum::Json<MillRunRequest>,
) -> Result<Json<MillRunResponse>, (StatusCode, Json<ErrorResponse>)> {
    if payload.l < 2 {
        return Err((
            StatusCode::BAD_REQUEST,
            Json(ErrorResponse {
                error: "l must be >= 2".into(),
            }),
        ));
    }
    if !(payload.beta.is_finite() && payload.beta > 0.0) {
        return Err((
            StatusCode::BAD_REQUEST,
            Json(ErrorResponse {
                error: "beta must be finite and > 0".into(),
            }),
        ));
    }
    if !(payload.step_size.is_finite() && payload.step_size > 0.0) {
        return Err((
            StatusCode::BAD_REQUEST,
            Json(ErrorResponse {
                error: "step_size must be finite and > 0".into(),
            }),
        ));
    }
    if payload.measure_every == 0 {
        return Err((
            StatusCode::BAD_REQUEST,
            Json(ErrorResponse {
                error: "measure_every must be >= 1".into(),
            }),
        ));
    }
    if payload.n_sweeps == 0 {
        return Err((
            StatusCode::BAD_REQUEST,
            Json(ErrorResponse {
                error: "n_sweeps must be >= 1".into(),
            }),
        ));
    }
    if payload.n_sweeps < payload.measure_every {
        return Err((
            StatusCode::BAD_REQUEST,
            Json(ErrorResponse {
                error: "n_sweeps must be >= measure_every".into(),
            }),
        ));
    }

    let cfg = crate::engine::mill::MillRunConfig {
        l: payload.l,
        beta: payload.beta,
        n_thermal_sweeps: payload.n_thermal_sweeps,
        n_sweeps: payload.n_sweeps,
        measure_every: payload.measure_every,
        step_size: payload.step_size,
        seed: payload.seed,
    };

    let out = crate::engine::mill::run_mill(cfg);
    Ok(Json(MillRunResponse { result: out }))
}

#[derive(Deserialize)]
struct MillRefineRequest {
    ls: Vec<usize>,
    beta: f64,
    n_thermal_sweeps: usize,
    n_sweeps: usize,
    measure_every: usize,
    step_size: f64,
    seed: u64,
    #[serde(default)]
    verdict_mode: Option<String>,
    #[serde(default)]
    gap_w_min: Option<usize>,
    #[serde(default)]
    gap_k_sigma: Option<f64>,
    #[serde(default)]
    plateau_rel_thresh: Option<f64>,
    #[serde(default)]
    plateau_mode: Option<String>,
    #[serde(default)]
    plateau_k: Option<f64>,
    #[serde(default)]
    plateau_chi2_max: Option<f64>,
    #[serde(default)]
    smeared_nonmax_fallback: Option<bool>,
}

async fn mill_refine(
    State(_): State<AppState>,
    axum::Json(payload): axum::Json<MillRefineRequest>,
) -> Result<Json<MillRefineResponse>, (StatusCode, Json<ErrorResponse>)> {
    if payload.ls.is_empty() {
        return Err((
            StatusCode::BAD_REQUEST,
            Json(ErrorResponse {
                error: "ls must be non-empty".into(),
            }),
        ));
    }
    if payload.ls.len() < 2 {
        return Err((
            StatusCode::BAD_REQUEST,
            Json(ErrorResponse {
                error: "ls must have len >= 2".into(),
            }),
        ));
    }
    if payload.ls.iter().any(|&l| l < 2) {
        return Err((
            StatusCode::BAD_REQUEST,
            Json(ErrorResponse {
                error: "all l must be >= 2".into(),
            }),
        ));
    }
    if !(payload.beta.is_finite() && payload.beta > 0.0) {
        return Err((
            StatusCode::BAD_REQUEST,
            Json(ErrorResponse {
                error: "beta must be finite and > 0".into(),
            }),
        ));
    }
    if !(payload.step_size.is_finite() && payload.step_size > 0.0) {
        return Err((
            StatusCode::BAD_REQUEST,
            Json(ErrorResponse {
                error: "step_size must be finite and > 0".into(),
            }),
        ));
    }
    if payload.measure_every == 0 {
        return Err((
            StatusCode::BAD_REQUEST,
            Json(ErrorResponse {
                error: "measure_every must be >= 1".into(),
            }),
        ));
    }
    if payload.n_sweeps == 0 {
        return Err((
            StatusCode::BAD_REQUEST,
            Json(ErrorResponse {
                error: "n_sweeps must be >= 1".into(),
            }),
        ));
    }
    if payload.n_sweeps < payload.measure_every {
        return Err((
            StatusCode::BAD_REQUEST,
            Json(ErrorResponse {
                error: "n_sweeps must be >= measure_every".into(),
            }),
        ));
    }

    if let Some(mode) = payload.verdict_mode.as_deref() {
        if !matches!(mode, "unanimous" | "ir_lmax") {
            return Err((
                StatusCode::BAD_REQUEST,
                Json(ErrorResponse {
                    error: "verdict_mode must be \"unanimous\" or \"ir_lmax\"".into(),
                }),
            ));
        }
    }

    if let Some(w_min) = payload.gap_w_min {
        if w_min < 1 {
            return Err((
                StatusCode::BAD_REQUEST,
                Json(ErrorResponse {
                    error: "gap_w_min must be >= 1".into(),
                }),
            ));
        }
    }
    if let Some(k_sigma) = payload.gap_k_sigma {
        if !(k_sigma.is_finite() && k_sigma > 0.0) {
            return Err((
                StatusCode::BAD_REQUEST,
                Json(ErrorResponse {
                    error: "gap_k_sigma must be finite and > 0".into(),
                }),
            ));
        }
    }
    if let Some(rel) = payload.plateau_rel_thresh {
        if !(rel.is_finite() && rel > 0.0 && rel <= 1.0) {
            return Err((
                StatusCode::BAD_REQUEST,
                Json(ErrorResponse {
                    error: "plateau_rel_thresh must be finite, > 0, and <= 1.0".into(),
                }),
            ));
        }
    }
    if let Some(mode) = payload.plateau_mode.as_deref() {
        if !matches!(mode, "legacy" | "stat") {
            return Err((
                StatusCode::BAD_REQUEST,
                Json(ErrorResponse {
                    error: "plateau_mode must be \"legacy\" or \"stat\"".into(),
                }),
            ));
        }
    }
    if let Some(k) = payload.plateau_k {
        if !(k.is_finite() && k > 0.0) {
            return Err((
                StatusCode::BAD_REQUEST,
                Json(ErrorResponse {
                    error: "plateau_k must be finite and > 0".into(),
                }),
            ));
        }
    }
    if let Some(x) = payload.plateau_chi2_max {
        if !(x.is_finite() && x > 0.0) {
            return Err((
                StatusCode::BAD_REQUEST,
                Json(ErrorResponse {
                    error: "plateau_chi2_max must be finite and > 0".into(),
                }),
            ));
        }
    }

    let cfg = crate::engine::mill::MillRefineConfig {
        ls: payload.ls,
        beta: payload.beta,
        n_thermal_sweeps: payload.n_thermal_sweeps,
        n_sweeps: payload.n_sweeps,
        measure_every: payload.measure_every,
        step_size: payload.step_size,
        seed: payload.seed,
        verdict_mode: payload.verdict_mode,
        gap_w_min: payload.gap_w_min,
        gap_k_sigma: payload.gap_k_sigma,
        plateau_rel_thresh: payload.plateau_rel_thresh,
        plateau_mode: payload.plateau_mode,
        plateau_k: payload.plateau_k,
        plateau_chi2_max: payload.plateau_chi2_max,
        smeared_nonmax_fallback: payload.smeared_nonmax_fallback,
    };

    let out = crate::engine::mill::run_mill_refine(cfg);
    Ok(Json(MillRefineResponse { result: out }))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::api;
    use axum::{body::Body, http::Request};
    use http_body_util::BodyExt;
    use serde_json::json;
    use tower::ServiceExt;

    fn router() -> Router {
        std::env::set_var("AUTH_TOKEN", "testtoken");
        api::build_router()
    }

    fn auth_header() -> (&'static str, &'static str) {
        (
            axum::http::header::AUTHORIZATION.as_str(),
            "Bearer testtoken",
        )
    }

    #[tokio::test]
    async fn health_ok() {
        let app = router();
        let response = app
            .oneshot(
                Request::builder()
                    .uri("/health")
                    .header(auth_header().0, auth_header().1)
                    .body(Body::empty())
                    .unwrap(),
            )
            .await
            .unwrap();
        assert_eq!(response.status(), StatusCode::OK);
    }

    #[tokio::test]
    async fn mill_refine_ok() {
        let app = router();
        let ls = vec![4usize, 8usize];
        let payload = json!({
            "ls": ls,
            "beta": 2.0,
            "n_thermal_sweeps": 2,
            "n_sweeps": 40,
            "measure_every": 2,
            "step_size": 0.3,
            "seed": 123
        });
        let response = app
            .clone()
            .oneshot(
                Request::builder()
                    .uri("/mill/refine")
                    .method("POST")
                    .header(auth_header().0, auth_header().1)
                    .header(axum::http::header::CONTENT_TYPE, "application/json")
                    .body(Body::from(payload.to_string()))
                    .unwrap(),
            )
            .await
            .unwrap();
        assert_eq!(response.status(), StatusCode::OK);

        let body_bytes = response.into_body().collect().await.unwrap().to_bytes();
        let v: serde_json::Value = serde_json::from_slice(&body_bytes).unwrap();

        let runs = v
            .get("result")
            .and_then(|r| r.get("runs"))
            .and_then(|r| r.as_array())
            .unwrap();
        assert_eq!(runs.len(), 2);

        let expected_ls = [4u64, 8u64];
        let expected_ls_str = ["4", "8"];
        let l_rp = 8u64;
        for (i, run) in runs.iter().enumerate() {
            let l = run.get("l").and_then(|x| x.as_u64()).unwrap();
            assert_eq!(l, expected_ls[i]);

            let pm = run.get("plaquette_mean").and_then(|x| x.as_f64()).unwrap();
            let ps = run.get("plaquette_std").and_then(|x| x.as_f64()).unwrap();
            assert!(pm.is_finite());
            assert!(ps.is_finite());

            let tests = run.get("tests").unwrap();
            let ti = tests
                .get("translation_invariance_max_row_dev")
                .and_then(|x| x.as_f64())
                .unwrap();
            let rp = tests
                .get("reflection_positivity_estimate")
                .and_then(|x| x.as_f64())
                .unwrap();
            assert!(ti.is_finite());
            assert!(rp.is_finite());
        }

        let deltas = v
            .get("result")
            .and_then(|r| r.get("convergence"))
            .and_then(|c| c.get("plaquette_mean_deltas"))
            .and_then(|d| d.as_array())
            .unwrap();
        assert_eq!(deltas.len(), 1);
        let deltas_f: Vec<f64> = deltas.iter().map(|x| x.as_f64().unwrap()).collect();
        assert!(deltas_f.iter().all(|x| x.is_finite()));
        let expected_max = deltas_f.iter().copied().fold(0.0, f64::max);

        let max_abs_delta = v
            .get("result")
            .and_then(|r| r.get("convergence"))
            .and_then(|c| c.get("max_abs_delta"))
            .and_then(|x| x.as_f64())
            .unwrap();
        assert!(max_abs_delta.is_finite());
        assert_eq!(max_abs_delta, expected_max);

        let inv = v
            .get("result")
            .and_then(|r| r.get("invariance_scaling"))
            .unwrap();
        let violations = inv.get("violations").and_then(|x| x.as_array()).unwrap();
        assert_eq!(violations.len(), runs.len());
        for (i, pt) in violations.iter().enumerate() {
            let l = pt.get("l").and_then(|x| x.as_u64()).unwrap();
            assert_eq!(l, expected_ls[i]);
            let val = pt.get("value").and_then(|x| x.as_f64()).unwrap();
            assert!(val.is_finite());
        }

        let trend = inv.get("trend").and_then(|x| x.as_str()).unwrap();
        assert!(matches!(trend, "decreasing" | "flat" | "increasing"));
        let slope = inv.get("slope_estimate").and_then(|x| x.as_f64()).unwrap();
        assert!(slope.is_finite());

        let rp = v
            .get("result")
            .and_then(|r| r.get("reflection_positivity"))
            .unwrap();
        let obs = rp.get("observables").and_then(|x| x.as_array()).unwrap();
        assert_eq!(obs.len(), 3);
        assert_eq!(obs[0].as_str().unwrap(), "F1");
        assert_eq!(obs[1].as_str().unwrap(), "F2");
        assert_eq!(obs[2].as_str().unwrap(), "F3");

        let matrix = rp.get("matrix").and_then(|x| x.as_array()).unwrap();
        assert_eq!(matrix.len(), 3);
        for row in matrix {
            let rowa = row.as_array().unwrap();
            assert_eq!(rowa.len(), 3);
            for x in rowa {
                let v = x.as_f64().unwrap();
                assert!(v.is_finite());
            }
        }

        let e = rp.get("eigenvalues").and_then(|x| x.as_array()).unwrap();
        assert_eq!(e.len(), 3);
        for x in e {
            let v = x.as_f64().unwrap();
            assert!(v.is_finite());
        }
        let min_ev = rp.get("min_eigenvalue").and_then(|x| x.as_f64()).unwrap();
        assert!(min_ev.is_finite());

        let classif = rp.get("classification").and_then(|x| x.as_str()).unwrap();
        assert!(matches!(classif, "positive" | "indefinite" | "negative"));
        let eps = 1e-8;
        let max_ev = e
            .iter()
            .map(|x| x.as_f64().unwrap())
            .fold(f64::NEG_INFINITY, f64::max);
        let expected = if min_ev >= -eps {
            "positive"
        } else if max_ev <= eps {
            "negative"
        } else {
            "indefinite"
        };
        assert_eq!(classif, expected);

        let me = v
            .get("result")
            .and_then(|r| r.get("mass_effective"))
            .unwrap();
        let l_me = me.get("l").and_then(|x| x.as_u64()).unwrap();
        assert_eq!(l_me, l_rp);
        let r_max = me.get("r_max").and_then(|x| x.as_u64()).unwrap();
        assert!(r_max >= 1);

        let corr = me.get("correlator").and_then(|x| x.as_array()).unwrap();
        assert_eq!(corr.len() as u64, r_max);
        for x in corr {
            let v = x.as_f64().unwrap();
            assert!(v.is_finite());
        }

        let m_eff = me.get("m_eff").and_then(|x| x.as_array()).unwrap();
        if r_max >= 2 {
            assert_eq!(m_eff.len() as u64, r_max - 1);
        } else {
            assert_eq!(m_eff.len(), 0);
        }
        for x in m_eff {
            let v = x.as_f64().unwrap();
            assert!(v.is_finite());
        }

        let pl = me.get("plateau").unwrap();
        let rs = pl.get("r_start").and_then(|x| x.as_u64()).unwrap();
        let re = pl.get("r_end").and_then(|x| x.as_u64()).unwrap();
        assert!(rs <= re);
        let pm = pl.get("m_eff_mean").and_then(|x| x.as_f64()).unwrap();
        let ps = pl.get("m_eff_std").and_then(|x| x.as_f64()).unwrap();
        assert!(pm.is_finite());
        assert!(ps.is_finite());

        let mes = v
            .get("result")
            .and_then(|r| r.get("mass_effective_scaling"))
            .unwrap();
        let series = mes.get("series").and_then(|x| x.as_array()).unwrap();
        assert_eq!(series.len(), runs.len());
        for (i, row) in series.iter().enumerate() {
            let l = row.get("l").and_then(|x| x.as_u64()).unwrap();
            assert_eq!(l, expected_ls[i]);
            let m = row.get("m_eff_mean").and_then(|x| x.as_f64()).unwrap();
            let s = row.get("m_eff_std").and_then(|x| x.as_f64()).unwrap();
            let w = row.get("plateau_width").and_then(|x| x.as_u64()).unwrap();
            assert!(m.is_finite());
            assert!(s.is_finite());
            let _ = w;
        }

        let delta_means = mes.get("delta_means").and_then(|x| x.as_array()).unwrap();
        assert_eq!(delta_means.len(), runs.len() - 1);
        for x in delta_means {
            let v = x.as_f64().unwrap();
            assert!(v.is_finite());
        }
        let madm = mes
            .get("max_abs_delta_mean")
            .and_then(|x| x.as_f64())
            .unwrap();
        assert!(madm.is_finite());

        let t = mes.get("trend_means").and_then(|x| x.as_str()).unwrap();
        assert!(matches!(t, "decreasing" | "flat" | "increasing"));
        let q = mes.get("plateau_quality").and_then(|x| x.as_str()).unwrap();
        assert!(matches!(q, "stable" | "marginal" | "unstable"));

        let gc = v
            .get("result")
            .and_then(|r| r.get("gap_compatibility"))
            .unwrap();
        let gc_op = gc.get("operator").and_then(|x| x.as_str()).unwrap();
        assert_eq!(gc_op, "raw");
        let tested = gc.get("tested_m0").and_then(|x| x.as_array()).unwrap();
        assert_eq!(tested.len(), 3);
        assert_eq!(tested[0].as_f64().unwrap(), 0.1);
        assert_eq!(tested[1].as_f64().unwrap(), 0.2);
        assert_eq!(tested[2].as_f64().unwrap(), 0.3);

        let per_l = gc.get("per_l").and_then(|x| x.as_object()).unwrap();
        assert_eq!(per_l.len(), runs.len());
        for l in expected_ls_str {
            let arr = per_l.get(l).and_then(|x| x.as_array()).unwrap();
            assert_eq!(arr.len(), tested.len());
            for s in arr {
                let st = s.as_str().unwrap();
                assert!(matches!(st, "compatible" | "incompatible" | "inconclusive"));
            }
        }

        let global = gc.get("global").and_then(|x| x.as_object()).unwrap();
        assert_eq!(global.len(), tested.len());
        for k in ["m0=0.1", "m0=0.2", "m0=0.3"] {
            let st = global.get(k).and_then(|x| x.as_str()).unwrap();
            assert!(matches!(st, "compatible" | "incompatible" | "inconclusive"));
        }

        let gcs = v
            .get("result")
            .and_then(|r| r.get("gap_compatibility_smeared"))
            .unwrap();
        let gcs_op = gcs.get("operator").and_then(|x| x.as_str()).unwrap();
        assert_eq!(gcs_op, "ape");
        let gcs_steps = gcs.get("steps").and_then(|x| x.as_u64()).unwrap();
        assert!(gcs_steps <= 3);
        let tested2 = gcs.get("tested_m0").and_then(|x| x.as_array()).unwrap();
        assert_eq!(tested2, tested);

        let per_l2 = gcs.get("per_l").and_then(|x| x.as_object()).unwrap();
        assert_eq!(per_l2.len(), runs.len());
        for l in expected_ls_str {
            let arr = per_l2.get(l).and_then(|x| x.as_array()).unwrap();
            assert_eq!(arr.len(), tested2.len());
            for s in arr {
                let st = s.as_str().unwrap();
                assert!(matches!(st, "compatible" | "incompatible" | "inconclusive"));
            }
        }

        let global2 = gcs.get("global").and_then(|x| x.as_object()).unwrap();
        assert_eq!(global2.len(), tested2.len());
        for k in ["m0=0.1", "m0=0.2", "m0=0.3"] {
            let st = global2.get(k).and_then(|x| x.as_str()).unwrap();
            assert!(matches!(st, "compatible" | "incompatible" | "inconclusive"));
        }

        let oc = v
            .get("result")
            .and_then(|r| r.get("operator_consistency"))
            .unwrap();
        let pair = oc.get("raw_vs_smeared").unwrap();
        let d = pair.get("delta_m_eff").and_then(|x| x.as_f64()).unwrap();
        let su = pair.get("sigma_units").and_then(|x| x.as_f64()).unwrap();
        let c2 = pair
            .get("consistent_2sigma")
            .and_then(|x| x.as_bool())
            .unwrap();
        assert!(d.is_finite());
        assert!(su.is_finite());
        let _ = c2;

        let mut v_mut = v.clone();
        v_mut["result"]["operator_smearing"]["ape"]["best"]["steps"] =
            serde_json::Value::from(3u64);
        assert_eq!(
            v["result"]["gap_compatibility"],
            v_mut["result"]["gap_compatibility"]
        );

        let os = v
            .get("result")
            .and_then(|r| r.get("operator_smearing"))
            .unwrap();
        let ape = os.get("ape").unwrap();
        let alpha = ape.get("alpha").and_then(|x| x.as_f64()).unwrap();
        assert_eq!(alpha, 0.5);
        let steps = ape.get("steps").and_then(|x| x.as_array()).unwrap();
        assert_eq!(steps.len(), 4);
        assert_eq!(steps[0].as_u64().unwrap(), 0);
        assert_eq!(steps[1].as_u64().unwrap(), 1);
        assert_eq!(steps[2].as_u64().unwrap(), 2);
        assert_eq!(steps[3].as_u64().unwrap(), 3);

        let results = ape.get("results").and_then(|x| x.as_object()).unwrap();
        assert_eq!(results.len(), 4);
        for k in ["0", "1", "2", "3"] {
            let row = results.get(k).unwrap();
            let w = row.get("plateau_width").and_then(|x| x.as_u64()).unwrap();
            let m = row.get("m_eff_mean").and_then(|x| x.as_f64()).unwrap();
            let s = row.get("m_eff_std").and_then(|x| x.as_f64()).unwrap();
            assert!(w <= l_rp);
            assert!(m.is_finite());
            assert!(s.is_finite());
        }

        let best = ape.get("best").unwrap();
        let bs = best.get("steps").and_then(|x| x.as_u64()).unwrap();
        assert!(bs <= 3);
        let crit = best.get("criterion").and_then(|x| x.as_str()).unwrap();
        assert_eq!(crit, "max_plateau_width");

        let fv = v
            .get("result")
            .and_then(|r| r.get("final_verdict"))
            .unwrap();
        let fv_status = fv.get("status").and_then(|x| x.as_str()).unwrap();
        assert!(matches!(
            fv_status,
            "compatible" | "incompatible" | "inconclusive"
        ));

        let basis = fv.get("basis").unwrap();
        let b_raw = basis.get("raw").and_then(|x| x.as_str()).unwrap();
        let b_sm = basis.get("smeared").and_then(|x| x.as_str()).unwrap();
        let b_ok = basis
            .get("consistency_ok")
            .and_then(|x| x.as_bool())
            .unwrap();

        let scalar_raw = {
            let g = gc.get("global").and_then(|x| x.as_object()).unwrap();
            let mut it = g.values();
            let first = it.next().unwrap().as_str().unwrap();
            if it.all(|v| v.as_str().unwrap() == first) {
                first
            } else {
                "inconclusive"
            }
        };
        let scalar_sm = {
            let g = gcs.get("global").and_then(|x| x.as_object()).unwrap();
            let mut it = g.values();
            let first = it.next().unwrap().as_str().unwrap();
            if it.all(|v| v.as_str().unwrap() == first) {
                first
            } else {
                "inconclusive"
            }
        };
        assert_eq!(b_raw, scalar_raw);
        assert_eq!(b_sm, scalar_sm);
        assert_eq!(b_ok, c2);

        let rule = fv.get("rule_applied").and_then(|x| x.as_str()).unwrap();
        assert!(matches!(rule, "R1" | "R2" | "R3"));
        let expl = fv.get("explanation").and_then(|x| x.as_str()).unwrap();
        assert!(!expl.is_empty());

        let expect_rule = if scalar_raw == scalar_sm
            && b_ok
            && (scalar_raw == "compatible" || scalar_raw == "incompatible")
        {
            "R1"
        } else if scalar_raw == scalar_sm && scalar_raw == "inconclusive" {
            "R3"
        } else {
            "R2"
        };
        assert_eq!(rule, expect_rule);
        if rule == "R1" {
            assert_eq!(fv_status, scalar_raw);
        } else {
            assert_eq!(fv_status, "inconclusive");
        }

        let response2 = app
            .oneshot(
                Request::builder()
                    .uri("/mill/refine")
                    .method("POST")
                    .header(auth_header().0, auth_header().1)
                    .header(axum::http::header::CONTENT_TYPE, "application/json")
                    .body(Body::from(payload.to_string()))
                    .unwrap(),
            )
            .await
            .unwrap();
        assert_eq!(response2.status(), StatusCode::OK);
        let body_bytes2 = response2.into_body().collect().await.unwrap().to_bytes();
        let v2: serde_json::Value = serde_json::from_slice(&body_bytes2).unwrap();
        assert_eq!(v, v2);
    }

    #[tokio::test]
    async fn mill_refine_rejects_bad_knobs() {
        let app = router();
        let payload = json!({
            "ls": [8, 16, 32, 64],
            "beta": 2.0,
            "n_thermal_sweeps": 1,
            "n_sweeps": 1,
            "measure_every": 1,
            "step_size": 0.3,
            "seed": 123,
            "gap_w_min": 0
        });
        let response = app
            .clone()
            .oneshot(
                Request::builder()
                    .uri("/mill/refine")
                    .method("POST")
                    .header(auth_header().0, auth_header().1)
                    .header(axum::http::header::CONTENT_TYPE, "application/json")
                    .body(Body::from(payload.to_string()))
                    .unwrap(),
            )
            .await
            .unwrap();
        assert_eq!(response.status(), StatusCode::BAD_REQUEST);

        let payload = json!({
            "ls": [8, 16, 32, 64],
            "beta": 2.0,
            "n_thermal_sweeps": 1,
            "n_sweeps": 1,
            "measure_every": 1,
            "step_size": 0.3,
            "seed": 123,
            "gap_k_sigma": 0.0
        });
        let response = app
            .clone()
            .oneshot(
                Request::builder()
                    .uri("/mill/refine")
                    .method("POST")
                    .header(auth_header().0, auth_header().1)
                    .header(axum::http::header::CONTENT_TYPE, "application/json")
                    .body(Body::from(payload.to_string()))
                    .unwrap(),
            )
            .await
            .unwrap();
        assert_eq!(response.status(), StatusCode::BAD_REQUEST);

        let payload = json!({
            "ls": [8, 16, 32, 64],
            "beta": 2.0,
            "n_thermal_sweeps": 1,
            "n_sweeps": 1,
            "measure_every": 1,
            "step_size": 0.3,
            "seed": 123,
            "plateau_rel_thresh": 2.0
        });
        let response = app
            .clone()
            .oneshot(
                Request::builder()
                    .uri("/mill/refine")
                    .method("POST")
                    .header(auth_header().0, auth_header().1)
                    .header(axum::http::header::CONTENT_TYPE, "application/json")
                    .body(Body::from(payload.to_string()))
                    .unwrap(),
            )
            .await
            .unwrap();
        assert_eq!(response.status(), StatusCode::BAD_REQUEST);

        let payload = json!({
            "ls": [8, 16, 32, 64],
            "beta": 2.0,
            "n_thermal_sweeps": 1,
            "n_sweeps": 1,
            "measure_every": 1,
            "step_size": 0.3,
            "seed": 123,
            "plateau_mode": "nope"
        });
        let response = app
            .clone()
            .oneshot(
                Request::builder()
                    .uri("/mill/refine")
                    .method("POST")
                    .header(auth_header().0, auth_header().1)
                    .header(axum::http::header::CONTENT_TYPE, "application/json")
                    .body(Body::from(payload.to_string()))
                    .unwrap(),
            )
            .await
            .unwrap();
        assert_eq!(response.status(), StatusCode::BAD_REQUEST);

        let payload = json!({
            "ls": [8, 16, 32, 64],
            "beta": 2.0,
            "n_thermal_sweeps": 1,
            "n_sweeps": 1,
            "measure_every": 1,
            "step_size": 0.3,
            "seed": 123,
            "verdict_mode": "nope"
        });
        let response = app
            .clone()
            .oneshot(
                Request::builder()
                    .uri("/mill/refine")
                    .method("POST")
                    .header(auth_header().0, auth_header().1)
                    .header(axum::http::header::CONTENT_TYPE, "application/json")
                    .body(Body::from(payload.to_string()))
                    .unwrap(),
            )
            .await
            .unwrap();
        assert_eq!(response.status(), StatusCode::BAD_REQUEST);
    }

    #[tokio::test]
    async fn mill_refine_ir_lmax_includes_report() {
        let app = router();
        let payload = json!({
            "ls": [8, 16],
            "beta": 2.0,
            "n_thermal_sweeps": 1,
            "n_sweeps": 2,
            "measure_every": 1,
            "step_size": 0.3,
            "seed": 123,
            "verdict_mode": "ir_lmax"
        });
        let response = app
            .oneshot(
                Request::builder()
                    .uri("/mill/refine")
                    .method("POST")
                    .header(auth_header().0, auth_header().1)
                    .header(axum::http::header::CONTENT_TYPE, "application/json")
                    .body(Body::from(payload.to_string()))
                    .unwrap(),
            )
            .await
            .unwrap();
        assert_eq!(response.status(), StatusCode::OK);
        let body_bytes = response.into_body().collect().await.unwrap().to_bytes();
        let v: serde_json::Value = serde_json::from_slice(&body_bytes).unwrap();
        let fv = v
            .get("result")
            .and_then(|r| r.get("final_verdict"))
            .unwrap();
        let rule = fv.get("rule_applied").and_then(|x| x.as_str()).unwrap();
        assert_eq!(rule, "ORDEN_015_IR_LMAX");
        let rep = fv.get("ir_lmax").unwrap();
        let l = rep.get("l").and_then(|x| x.as_u64()).unwrap();
        assert_eq!(l, 16);
        let channel = rep.get("channel").and_then(|x| x.as_str()).unwrap();
        assert!(matches!(channel, "raw" | "smeared_best"));
        let mean = rep.get("m_eff_mean").and_then(|x| x.as_f64()).unwrap();
        let std = rep.get("m_eff_std").and_then(|x| x.as_f64()).unwrap();
        assert!(mean.is_finite());
        assert!(std.is_finite());
        let w = rep.get("plateau_width").and_then(|x| x.as_u64()).unwrap();
        let _ = w;
        let per_m0 = rep.get("per_m0").and_then(|x| x.as_object()).unwrap();
        for k in ["m0=0.1", "m0=0.2", "m0=0.3"] {
            let st = per_m0.get(k).and_then(|x| x.as_str()).unwrap();
            assert!(matches!(st, "compatible" | "incompatible" | "inconclusive"));
        }
    }
}
