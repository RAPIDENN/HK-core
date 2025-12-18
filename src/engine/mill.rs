use rand::{rngs::StdRng, Rng, SeedableRng};
use serde::Serialize;
use std::collections::BTreeMap;
use std::f64::consts::{PI, TAU};

const JK_TARGET_N_BLOCKS: usize = 20;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum PlateauMode {
    Legacy,
    Stat,
}

#[derive(Clone, Copy, Debug)]
struct PlateauCfg {
    mode: PlateauMode,
    rel_thresh: f64,
    k: f64,
    chi2_max: f64,
}

fn parse_plateau_mode(mode: Option<&str>) -> PlateauMode {
    match mode {
        Some("stat") => PlateauMode::Stat,
        _ => PlateauMode::Legacy,
    }
}

#[derive(Clone, Debug)]
pub struct MillRunConfig {
    pub l: usize,
    pub beta: f64,
    pub n_thermal_sweeps: usize,
    pub n_sweeps: usize,
    pub measure_every: usize,
    pub step_size: f64,
    pub seed: u64,
}

fn plateau_width_from_r_range(r_start: usize, r_end: usize) -> usize {
    r_end.saturating_sub(r_start)
}

fn mass_scaling_plateau_from_acc_with_cfg(
    acc: &MassAccumulator,
    cfg: PlateauCfg,
) -> Option<MassPlateau> {
    acc.plateau_point(cfg)
}

fn analyze_mass_effective_scaling(series: &[MassEffectiveScalingPoint]) -> MassEffectiveScaling {
    let mut delta_means: Vec<f64> = Vec::new();
    for i in 1..series.len() {
        delta_means.push((series[i].m_eff_mean - series[i - 1].m_eff_mean).abs());
    }
    let max_abs_delta_mean = delta_means.iter().copied().fold(0.0, f64::max);

    let slope_estimate = linear_slope_l_vs_me(series);
    let trend_means = match slope_estimate {
        Some(s) if s < -1e-8 => "decreasing",
        Some(s) if s > 1e-8 => "increasing",
        _ => "flat",
    }
    .to_string();

    let w_min = 3usize;
    let delta_max = 0.2f64;

    let all_wide = series.iter().all(|p| p.plateau_width >= w_min);
    let any_plateau = series.iter().any(|p| p.plateau_width >= 1);
    let plateau_quality = if all_wide && max_abs_delta_mean <= delta_max {
        "stable"
    } else if any_plateau {
        "marginal"
    } else {
        "unstable"
    }
    .to_string();

    MassEffectiveScaling {
        series: series.to_vec(),
        delta_means,
        max_abs_delta_mean,
        trend_means,
        plateau_quality,
    }
}

fn linear_slope_l_vs_me(series: &[MassEffectiveScalingPoint]) -> Option<f64> {
    if series.len() < 2 {
        return None;
    }
    let n = series.len() as f64;
    let xs: Vec<f64> = series.iter().map(|p| p.l as f64).collect();
    let ys: Vec<f64> = series.iter().map(|p| p.m_eff_mean).collect();

    let x_mean = xs.iter().sum::<f64>() / n;
    let y_mean = ys.iter().sum::<f64>() / n;

    let mut num = 0.0;
    let mut den = 0.0;
    for i in 0..series.len() {
        let dx = xs[i] - x_mean;
        let dy = ys[i] - y_mean;
        num += dx * dy;
        den += dx * dx;
    }
    if den == 0.0 {
        return None;
    }
    let slope = num / den;
    if slope.is_finite() {
        Some(slope)
    } else {
        None
    }
}

#[derive(Serialize, Clone, Debug)]
pub struct MillRunOutput {
    pub trace_id: String,
    pub lattice: LatticeSummary,
    pub observables: ObservablesSummary,
    pub tests: TestsSummary,
}

#[derive(Serialize, Clone, Debug)]
pub struct LatticeSummary {
    pub dim: u8,
    pub l: usize,
    pub beta: f64,
    pub n_links: usize,
    pub n_plaquettes: usize,
    pub step_size: f64,
}

#[derive(Serialize, Clone, Debug)]
pub struct ObservablesSummary {
    pub n_measurements: usize,
    pub plaquette_mean: f64,
    pub plaquette_std: f64,
}

#[derive(Serialize, Clone, Debug)]
pub struct TestsSummary {
    pub translation_invariance_max_row_dev: f64,
    pub reflection_positivity_estimate: f64,
}

#[derive(Clone, Debug)]
pub struct MillRefineConfig {
    pub ls: Vec<usize>,
    pub beta: f64,
    pub n_thermal_sweeps: usize,
    pub n_sweeps: usize,
    pub measure_every: usize,
    pub step_size: f64,
    pub seed: u64,
    pub verdict_mode: Option<String>,
    pub gap_w_min: Option<usize>,
    pub gap_k_sigma: Option<f64>,
    pub plateau_rel_thresh: Option<f64>,
    pub plateau_mode: Option<String>,
    pub plateau_k: Option<f64>,
    pub plateau_chi2_max: Option<f64>,
    pub smeared_nonmax_fallback: Option<bool>,
}

#[derive(Serialize, Clone, Debug)]
pub struct MillRefineOutput {
    pub trace_id: String,
    pub runs: Vec<MillRefineRunRow>,
    pub convergence: ConvergenceSummary,
    pub invariance_scaling: InvarianceScaling,
    pub reflection_positivity: ReflectionPositivityReport,
    pub mass_effective: MassEffectiveReport,
    pub mass_effective_scaling: MassEffectiveScaling,
    pub gap_compatibility: GapCompatibility,
    pub gap_compatibility_smeared: GapCompatibilitySmeared,
    pub operator_consistency: OperatorConsistencyReport,
    pub operator_smearing: OperatorSmearingReport,
    pub final_verdict: FinalVerdict,
}

#[derive(Serialize, Clone, Debug)]
pub struct ConvergenceSummary {
    pub plaquette_mean_deltas: Vec<f64>,
    pub max_abs_delta: f64,
}

#[derive(Serialize, Clone, Debug)]
pub struct MillRefineRunRow {
    pub l: usize,
    pub plaquette_mean: f64,
    pub plaquette_std: f64,
    pub tests: TestsSummary,
}

#[derive(Serialize, Clone, Debug)]
pub struct InvarianceScaling {
    pub violations: Vec<InvarianceViolationPoint>,
    pub trend: String,
    pub slope_estimate: Option<f64>,
}

#[derive(Serialize, Clone, Debug)]
pub struct InvarianceViolationPoint {
    pub l: usize,
    pub value: f64,
}

#[derive(Serialize, Clone, Debug)]
pub struct ReflectionPositivityReport {
    pub l: usize,
    pub observables: Vec<String>,
    pub matrix: Vec<Vec<f64>>,
    pub eigenvalues: Vec<f64>,
    pub min_eigenvalue: f64,
    pub classification: String,
}

#[derive(Serialize, Clone, Debug)]
pub struct MassEffectiveReport {
    pub l: usize,
    pub r_max: usize,
    pub correlator: Vec<f64>,
    pub m_eff: Vec<f64>,
    pub plateau: MassPlateau,
}

#[derive(Serialize, Clone, Debug)]
pub struct MassPlateau {
    pub r_start: usize,
    pub r_end: usize,
    pub m_eff_mean: f64,
    pub m_eff_std: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub method: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub chi2_dof: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub n_points: Option<usize>,
}

#[derive(Serialize, Clone, Debug)]
pub struct MassEffectiveScaling {
    pub series: Vec<MassEffectiveScalingPoint>,
    pub delta_means: Vec<f64>,
    pub max_abs_delta_mean: f64,
    pub trend_means: String,
    pub plateau_quality: String,
}

#[derive(Serialize, Clone, Debug)]
pub struct MassEffectiveScalingPoint {
    pub l: usize,
    pub m_eff_mean: f64,
    pub m_eff_std: f64,
    pub plateau_width: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub plateau_r_start: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub plateau_r_end: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub plateau_method: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub plateau_chi2_dof: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub plateau_n_points: Option<usize>,
}

#[derive(Serialize, Clone, Debug)]
pub struct GapCompatibility {
    pub operator: String,
    pub tested_m0: Vec<f64>,
    pub per_l: BTreeMap<String, Vec<String>>,
    pub global: BTreeMap<String, String>,
}

#[derive(Serialize, Clone, Debug)]
pub struct GapCompatibilitySmeared {
    pub operator: String,
    pub steps: usize,
    pub tested_m0: Vec<f64>,
    pub per_l_stats: BTreeMap<String, OperatorSmearingResult>,
    pub per_l: BTreeMap<String, Vec<String>>,
    pub global: BTreeMap<String, String>,
}

#[derive(Serialize, Clone, Debug)]
pub struct OperatorConsistencyReport {
    pub raw_vs_smeared: OperatorConsistencyPair,
}

#[derive(Serialize, Clone, Debug)]
pub struct OperatorConsistencyPair {
    pub delta_m_eff: f64,
    pub sigma_units: f64,
    pub consistent_2sigma: bool,
}

#[derive(Serialize, Clone, Debug)]
pub struct FinalVerdict {
    pub status: String,
    pub basis: FinalVerdictBasis,
    pub rule_applied: String,
    pub explanation: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ir_lmax: Option<IrLmaxVerdictReport>,
}

#[derive(Serialize, Clone, Debug)]
pub struct FinalVerdictBasis {
    pub raw: String,
    pub smeared: String,
    pub consistency_ok: bool,
}

#[derive(Serialize, Clone, Debug)]
pub struct IrLmaxVerdictReport {
    pub l: usize,
    pub channel: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub smeared_steps: Option<usize>,
    pub width_min: usize,
    pub plateau_width: usize,
    pub k_sigma: f64,
    pub tested_m0: Vec<f64>,
    pub per_m0: BTreeMap<String, String>,
    pub m_eff_mean: f64,
    pub m_eff_std: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub plateau_method: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub chi2_dof: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub chi2_max: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub chi2_ok: Option<bool>,
}

#[derive(Serialize, Clone, Debug)]
pub struct OperatorSmearingReport {
    pub ape: ApeSmearingReport,
}

#[derive(Serialize, Clone, Debug)]
pub struct ApeSmearingReport {
    pub alpha: f64,
    pub steps: Vec<usize>,
    pub results: BTreeMap<String, OperatorSmearingResult>,
    pub best: OperatorSmearingBest,
}

#[derive(Serialize, Clone, Debug)]
pub struct OperatorSmearingResult {
    pub plateau_width: usize,
    pub m_eff_mean: f64,
    pub m_eff_std: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub plateau_r_start: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub plateau_r_end: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub plateau_method: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub plateau_chi2_dof: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub plateau_n_points: Option<usize>,
}

#[derive(Serialize, Clone, Debug)]
pub struct OperatorSmearingBest {
    pub steps: usize,
    pub criterion: String,
}

#[derive(Clone, Copy, Debug)]
struct Su2 {
    a0: f64,
    a1: f64,
    a2: f64,
    a3: f64,
}

impl Su2 {
    fn identity() -> Self {
        Self {
            a0: 1.0,
            a1: 0.0,
            a2: 0.0,
            a3: 0.0,
        }
    }

    fn dagger(self) -> Self {
        Self {
            a0: self.a0,
            a1: -self.a1,
            a2: -self.a2,
            a3: -self.a3,
        }
    }

    fn mul(self, b: Self) -> Self {
        // Quaternion multiplication for SU(2): U = a0 I + i a·σ.
        let a0 = self.a0;
        let a1 = self.a1;
        let a2 = self.a2;
        let a3 = self.a3;
        let b0 = b.a0;
        let b1 = b.a1;
        let b2 = b.a2;
        let b3 = b.a3;
        Self {
            a0: a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3,
            a1: a0 * b1 + a1 * b0 + a2 * b3 - a3 * b2,
            a2: a0 * b2 - a1 * b3 + a2 * b0 + a3 * b1,
            a3: a0 * b3 + a1 * b2 - a2 * b1 + a3 * b0,
        }
    }

    fn add(self, b: Self) -> Self {
        Self {
            a0: self.a0 + b.a0,
            a1: self.a1 + b.a1,
            a2: self.a2 + b.a2,
            a3: self.a3 + b.a3,
        }
    }

    fn scale(self, s: f64) -> Self {
        Self {
            a0: self.a0 * s,
            a1: self.a1 * s,
            a2: self.a2 * s,
            a3: self.a3 * s,
        }
    }

    fn norm2(self) -> f64 {
        self.a0 * self.a0 + self.a1 * self.a1 + self.a2 * self.a2 + self.a3 * self.a3
    }

    fn normalized_or(self, fallback: Self) -> Self {
        let n2 = self.norm2();
        if !(n2.is_finite() && n2 > 0.0) {
            return fallback;
        }
        let inv = 1.0 / n2.sqrt();
        Self {
            a0: self.a0 * inv,
            a1: self.a1 * inv,
            a2: self.a2 * inv,
            a3: self.a3 * inv,
        }
    }

    fn plaquette_value(self) -> f64 {
        // For SU(2), (1/2) Re Tr(U) == a0 in this parameterization.
        self.a0
    }

    fn projected(self) -> Self {
        self.normalized_or(Su2::identity())
    }
}

#[derive(Clone, Debug)]
struct Su2Gauge3D {
    l: usize,
    ux: Vec<Su2>,
    uy: Vec<Su2>,
    uz: Vec<Su2>,
}

impl Su2Gauge3D {
    fn new(l: usize) -> Self {
        let n = l * l * l;
        Self {
            l,
            ux: vec![Su2::identity(); n],
            uy: vec![Su2::identity(); n],
            uz: vec![Su2::identity(); n],
        }
    }

    #[inline]
    fn idx(&self, x: usize, y: usize, z: usize) -> usize {
        x + self.l * (y + self.l * z)
    }

    #[inline]
    fn wrap(&self, i: isize) -> usize {
        let l = self.l as isize;
        (((i % l) + l) % l) as usize
    }

    fn link_x(&self, x: usize, y: usize, z: usize) -> Su2 {
        self.ux[self.idx(x, y, z)]
    }

    fn link_y(&self, x: usize, y: usize, z: usize) -> Su2 {
        self.uy[self.idx(x, y, z)]
    }

    fn link_z(&self, x: usize, y: usize, z: usize) -> Su2 {
        self.uz[self.idx(x, y, z)]
    }

    fn set_link_x(&mut self, x: usize, y: usize, z: usize, v: Su2) {
        let idx = self.idx(x, y, z);
        self.ux[idx] = v;
    }

    fn set_link_y(&mut self, x: usize, y: usize, z: usize, v: Su2) {
        let idx = self.idx(x, y, z);
        self.uy[idx] = v;
    }

    fn set_link_z(&mut self, x: usize, y: usize, z: usize, v: Su2) {
        let idx = self.idx(x, y, z);
        self.uz[idx] = v;
    }

    fn plaquette_xy(&self, x: usize, y: usize, z: usize) -> Su2 {
        let xp = self.wrap(x as isize + 1);
        let yp = self.wrap(y as isize + 1);
        let u1 = self.link_x(x, y, z);
        let u2 = self.link_y(xp, y, z);
        let u3 = self.link_x(x, yp, z).dagger();
        let u4 = self.link_y(x, y, z).dagger();
        u1.mul(u2).mul(u3).mul(u4)
    }

    fn plaquette_xz(&self, x: usize, y: usize, z: usize) -> Su2 {
        let xp = self.wrap(x as isize + 1);
        let zp = self.wrap(z as isize + 1);
        let u1 = self.link_x(x, y, z);
        let u2 = self.link_z(xp, y, z);
        let u3 = self.link_x(x, y, zp).dagger();
        let u4 = self.link_z(x, y, z).dagger();
        u1.mul(u2).mul(u3).mul(u4)
    }

    fn plaquette_yz(&self, x: usize, y: usize, z: usize) -> Su2 {
        let yp = self.wrap(y as isize + 1);
        let zp = self.wrap(z as isize + 1);
        let u1 = self.link_y(x, y, z);
        let u2 = self.link_z(x, yp, z);
        let u3 = self.link_y(x, y, zp).dagger();
        let u4 = self.link_z(x, y, z).dagger();
        u1.mul(u2).mul(u3).mul(u4)
    }

    fn plaquette_cos_xy(&self, x: usize, y: usize, z: usize) -> f64 {
        self.plaquette_xy(x, y, z).plaquette_value()
    }

    fn plaquette_cos_xz(&self, x: usize, y: usize, z: usize) -> f64 {
        self.plaquette_xz(x, y, z).plaquette_value()
    }

    fn plaquette_cos_yz(&self, x: usize, y: usize, z: usize) -> f64 {
        self.plaquette_yz(x, y, z).plaquette_value()
    }

    fn plaquette_cos(&self, x: usize, y: usize, z: usize) -> f64 {
        let pxy = self.plaquette_cos_xy(x, y, z);
        let pxz = self.plaquette_cos_xz(x, y, z);
        let pyz = self.plaquette_cos_yz(x, y, z);
        (pxy + pxz + pyz) / 3.0
    }

    fn plaquette_mean(&self) -> f64 {
        let mut s = 0.0;
        let l = self.l;
        for z in 0..l {
            for y in 0..l {
                for x in 0..l {
                    s += self.plaquette_cos(x, y, z);
                }
            }
        }
        s / ((l * l * l) as f64)
    }

    fn plaquette_mean_by_row(&self) -> Vec<f64> {
        let l = self.l;
        let mut row = vec![0.0; l];
        for y in 0..l {
            let mut s = 0.0;
            for z in 0..l {
                for x in 0..l {
                    s += self.plaquette_cos(x, y, z);
                }
            }
            row[y] = s / ((l * l) as f64);
        }
        row
    }
}

pub fn run_mill_refine(cfg: MillRefineConfig) -> MillRefineOutput {
    let w_min = cfg.gap_w_min.unwrap_or(3);
    let k_sigma = cfg.gap_k_sigma.unwrap_or(2.0);
    let plateau_rel_thresh = cfg.plateau_rel_thresh.unwrap_or(0.05);
    let plateau_mode = parse_plateau_mode(cfg.plateau_mode.as_deref());
    let plateau_k = cfg.plateau_k.unwrap_or(2.0);
    let plateau_chi2_max = cfg.plateau_chi2_max.unwrap_or(2.0);
    let plateau_cfg = PlateauCfg {
        mode: plateau_mode,
        rel_thresh: plateau_rel_thresh,
        k: plateau_k,
        chi2_max: plateau_chi2_max,
    };
    let smeared_nonmax_fallback = cfg.smeared_nonmax_fallback.unwrap_or(false);
    let verdict_mode = cfg.verdict_mode.as_deref().unwrap_or("unanimous");

    let l_rp = cfg.ls.iter().copied().fold(0usize, usize::max);
    let expected_n_measurements = if cfg.measure_every > 0 {
        cfg.n_sweeps / cfg.measure_every
    } else {
        0
    };
    let mut rp_acc: Option<RpAccumulator> = None;
    let mut mass_acc: Option<MassAccumulator> = None;
    let mut mass_series: Vec<MassEffectiveScalingPoint> = Vec::new();
    let mut smear_by_l: BTreeMap<usize, OperatorSmearingAccumulator> = BTreeMap::new();

    let mut full_runs: Vec<MillRunOutput> = Vec::new();
    for &l in cfg.ls.iter() {
        if l == l_rp {
            let mut acc = RpAccumulator::new(l);
            let mut macc = MassAccumulator::new(l, expected_n_measurements);
            let mut sacc = OperatorSmearingAccumulator::new(l, expected_n_measurements);
            let out = run_mill_with_rp(
                &MillRunConfig {
                    l,
                    beta: cfg.beta,
                    n_thermal_sweeps: cfg.n_thermal_sweeps,
                    n_sweeps: cfg.n_sweeps,
                    measure_every: cfg.measure_every,
                    step_size: cfg.step_size,
                    seed: cfg.seed,
                },
                Some(&mut acc),
                Some(&mut macc),
                Some(&mut sacc),
            );
            rp_acc = Some(acc);
            if let Some(pt) = mass_scaling_plateau_from_acc_with_cfg(&macc, plateau_cfg) {
                let width = plateau_width_from_r_range(pt.r_start, pt.r_end);
                mass_series.push(MassEffectiveScalingPoint {
                    l,
                    m_eff_mean: pt.m_eff_mean,
                    m_eff_std: pt.m_eff_std,
                    plateau_width: width,
                    plateau_r_start: Some(pt.r_start),
                    plateau_r_end: Some(pt.r_end),
                    plateau_method: pt.method,
                    plateau_chi2_dof: pt.chi2_dof,
                    plateau_n_points: pt.n_points,
                });
            }
            mass_acc = Some(macc);
            smear_by_l.insert(l, sacc);
            full_runs.push(out);
        } else {
            let mut macc = MassAccumulator::new(l, expected_n_measurements);
            let mut sacc = OperatorSmearingAccumulator::new(l, expected_n_measurements);
            let out = run_mill_with_rp(
                &MillRunConfig {
                    l,
                    beta: cfg.beta,
                    n_thermal_sweeps: cfg.n_thermal_sweeps,
                    n_sweeps: cfg.n_sweeps,
                    measure_every: cfg.measure_every,
                    step_size: cfg.step_size,
                    seed: cfg.seed,
                },
                None,
                Some(&mut macc),
                Some(&mut sacc),
            );
            if let Some(pt) = mass_scaling_plateau_from_acc_with_cfg(&macc, plateau_cfg) {
                let width = plateau_width_from_r_range(pt.r_start, pt.r_end);
                mass_series.push(MassEffectiveScalingPoint {
                    l,
                    m_eff_mean: pt.m_eff_mean,
                    m_eff_std: pt.m_eff_std,
                    plateau_width: width,
                    plateau_r_start: Some(pt.r_start),
                    plateau_r_end: Some(pt.r_end),
                    plateau_method: pt.method,
                    plateau_chi2_dof: pt.chi2_dof,
                    plateau_n_points: pt.n_points,
                });
            }
            smear_by_l.insert(l, sacc);
            full_runs.push(out);
        }
    }

    let mut deltas: Vec<f64> = Vec::new();
    for i in 1..full_runs.len() {
        deltas.push(
            (full_runs[i].observables.plaquette_mean - full_runs[i - 1].observables.plaquette_mean)
                .abs(),
        );
    }
    let max_abs_delta = deltas.iter().copied().fold(0.0, f64::max);

    let runs: Vec<MillRefineRunRow> = full_runs
        .iter()
        .map(|r| MillRefineRunRow {
            l: r.lattice.l,
            plaquette_mean: r.observables.plaquette_mean,
            plaquette_std: r.observables.plaquette_std,
            tests: r.tests.clone(),
        })
        .collect();

    let invariance_scaling = analyze_invariance_scaling(&full_runs);

    let reflection_positivity = build_reflection_positivity_report(l_rp, rp_acc);

    let mass_effective = build_mass_effective_report(l_rp, mass_acc, plateau_cfg);

    let mass_effective_scaling = analyze_mass_effective_scaling(&mass_series);

    let gap_compatibility =
        analyze_gap_compatibility(&mass_effective_scaling.series, w_min, k_sigma);

    let operator_smearing = build_operator_smearing_report(smear_by_l.remove(&l_rp), plateau_cfg);

    let best_steps = operator_smearing.ape.best.steps;
    let mut per_l_smeared_stats: BTreeMap<String, OperatorSmearingResult> = BTreeMap::new();
    for p in mass_effective_scaling.series.iter() {
        let l = p.l;
        let key = l.to_string();
        if l == l_rp {
            let v = operator_smearing
                .ape
                .results
                .get(&best_steps.to_string())
                .cloned()
                .unwrap_or(OperatorSmearingResult {
                    plateau_width: 0,
                    m_eff_mean: 0.0,
                    m_eff_std: 0.0,
                    plateau_r_start: None,
                    plateau_r_end: None,
                    plateau_method: None,
                    plateau_chi2_dof: None,
                    plateau_n_points: None,
                });
            per_l_smeared_stats.insert(key, v);
        } else if let Some(acc) = smear_by_l.get(&l) {
            per_l_smeared_stats
                .insert(key, smearing_result_for_steps(acc, best_steps, plateau_cfg));
        } else {
            per_l_smeared_stats.insert(
                key,
                OperatorSmearingResult {
                    plateau_width: 0,
                    m_eff_mean: 0.0,
                    m_eff_std: 0.0,
                    plateau_r_start: None,
                    plateau_r_end: None,
                    plateau_method: None,
                    plateau_chi2_dof: None,
                    plateau_n_points: None,
                },
            );
        }
    }

    let gap_compatibility_smeared = analyze_gap_compatibility_smeared(
        &mass_effective_scaling.series,
        &operator_smearing,
        &per_l_smeared_stats,
        l_rp,
        w_min,
        k_sigma,
        smeared_nonmax_fallback,
    );

    let operator_consistency =
        build_operator_consistency(&mass_effective_scaling.series, &operator_smearing, l_rp);

    let final_verdict = synthesize_final_verdict_with_mode(
        &gap_compatibility,
        &gap_compatibility_smeared,
        &operator_consistency,
        &mass_effective_scaling.series,
        &operator_smearing,
        plateau_cfg,
        l_rp,
        k_sigma,
        verdict_mode,
    );

    MillRefineOutput {
        trace_id: format!("MILL_REFINE_SU2_3D_b{}_seed{}", cfg.beta, cfg.seed),
        runs,
        convergence: ConvergenceSummary {
            plaquette_mean_deltas: deltas,
            max_abs_delta,
        },
        invariance_scaling,
        reflection_positivity,
        mass_effective,
        mass_effective_scaling,
        gap_compatibility,
        gap_compatibility_smeared,
        operator_consistency,
        operator_smearing,
        final_verdict,
    }
}

fn scalar_global_verdict(global: &BTreeMap<String, String>) -> String {
    let mut it = global.values();
    let Some(first) = it.next() else {
        return "inconclusive".to_string();
    };
    if it.all(|v| v == first) {
        first.clone()
    } else {
        "inconclusive".to_string()
    }
}

fn synthesize_final_verdict(
    raw: &GapCompatibility,
    smeared: &GapCompatibilitySmeared,
    consistency: &OperatorConsistencyReport,
) -> FinalVerdict {
    let v_raw = scalar_global_verdict(&raw.global);
    let v_sm = scalar_global_verdict(&smeared.global);
    let consistency_ok = consistency.raw_vs_smeared.consistent_2sigma;

    let basis = FinalVerdictBasis {
        raw: v_raw.clone(),
        smeared: v_sm.clone(),
        consistency_ok,
    };

    let (status, rule_applied, explanation) = if v_raw == v_sm
        && consistency_ok
        && (v_raw == "compatible" || v_raw == "incompatible")
    {
        (
            v_raw.clone(),
            "R1".to_string(),
            "Raw and smeared operators agree on a decisive verdict and are consistent within 2σ."
                .to_string(),
        )
    } else if v_raw == v_sm && v_raw == "inconclusive" {
        (
            "inconclusive".to_string(),
            "R3".to_string(),
            "Both raw and smeared operators are inconclusive; current evidence does not support a gap decision.".to_string(),
        )
    } else {
        (
            "inconclusive".to_string(),
            "R2".to_string(),
            "Raw and smeared operators disagree or are inconsistent; no gap decision allowed."
                .to_string(),
        )
    };

    FinalVerdict {
        status,
        basis,
        rule_applied,
        explanation,
        ir_lmax: None,
    }
}

fn synthesize_final_verdict_with_mode(
    raw: &GapCompatibility,
    smeared: &GapCompatibilitySmeared,
    consistency: &OperatorConsistencyReport,
    raw_series: &[MassEffectiveScalingPoint],
    operator_smearing: &OperatorSmearingReport,
    plateau_cfg: PlateauCfg,
    l_max: usize,
    k_sigma: f64,
    verdict_mode: &str,
) -> FinalVerdict {
    if verdict_mode != "ir_lmax" {
        return synthesize_final_verdict(raw, smeared, consistency);
    }

    let v_raw = scalar_global_verdict(&raw.global);
    let v_sm = scalar_global_verdict(&smeared.global);
    let consistency_ok = consistency.raw_vs_smeared.consistent_2sigma;

    let basis = FinalVerdictBasis {
        raw: v_raw,
        smeared: v_sm,
        consistency_ok,
    };

    let (status, report) = verdict_ir_lmax(
        raw_series,
        operator_smearing,
        plateau_cfg,
        l_max,
        k_sigma,
        consistency_ok,
    );

    let explanation = format!(
        "ORDEN 015 (IR Lmax): channel={} Lmax={} width={} (min={}) k_sigma={}",
        report.channel, report.l, report.plateau_width, report.width_min, report.k_sigma
    );

    FinalVerdict {
        status,
        basis,
        rule_applied: "ORDEN_015_IR_LMAX".to_string(),
        explanation,
        ir_lmax: Some(report),
    }
}

fn verdict_ir_lmax(
    raw_series: &[MassEffectiveScalingPoint],
    operator_smearing: &OperatorSmearingReport,
    plateau_cfg: PlateauCfg,
    l_max: usize,
    k_sigma: f64,
    consistency_ok: bool,
) -> (String, IrLmaxVerdictReport) {
    let tested_m0 = vec![0.1, 0.2, 0.3];
    let width_min = 6usize;

    let (channel, smeared_steps, plateau_width, mean, std, method, chi2_dof) = if consistency_ok {
        let steps = operator_smearing.ape.best.steps;
        let r = operator_smearing
            .ape
            .results
            .get(&steps.to_string())
            .cloned()
            .unwrap_or(OperatorSmearingResult {
                plateau_width: 0,
                m_eff_mean: 0.0,
                m_eff_std: 0.0,
                plateau_r_start: None,
                plateau_r_end: None,
                plateau_method: None,
                plateau_chi2_dof: None,
                plateau_n_points: None,
            });
        (
            "smeared_best".to_string(),
            Some(steps),
            r.plateau_width,
            r.m_eff_mean,
            r.m_eff_std,
            r.plateau_method,
            r.plateau_chi2_dof,
        )
    } else {
        let p = raw_series.iter().find(|p| p.l == l_max).cloned();
        let (width, mean, std, method, chi2_dof) = if let Some(p) = p {
            (
                p.plateau_width,
                p.m_eff_mean,
                p.m_eff_std,
                p.plateau_method,
                p.plateau_chi2_dof,
            )
        } else {
            (0usize, 0.0f64, 0.0f64, None, None)
        };
        ("raw".to_string(), None, width, mean, std, method, chi2_dof)
    };

    let chi2_ok = if plateau_cfg.mode == PlateauMode::Stat {
        chi2_dof.map(|x| x.is_finite() && x <= plateau_cfg.chi2_max)
    } else {
        None
    };

    let mut per_m0: BTreeMap<String, String> = BTreeMap::new();
    for &m0 in tested_m0.iter() {
        let st = if plateau_width < width_min {
            "inconclusive"
        } else if !(mean.is_finite() && std.is_finite()) {
            "inconclusive"
        } else if std <= 0.0 {
            "inconclusive"
        } else if matches!(chi2_ok, Some(false)) {
            "inconclusive"
        } else if mean + k_sigma * std < m0 {
            "incompatible"
        } else if mean - k_sigma * std > m0 {
            "compatible"
        } else {
            "inconclusive"
        };
        per_m0.insert(format!("m0={}", m0), st.to_string());
    }

    let status = scalar_global_verdict(&per_m0);

    let report = IrLmaxVerdictReport {
        l: l_max,
        channel,
        smeared_steps,
        width_min,
        plateau_width,
        k_sigma,
        tested_m0,
        per_m0,
        m_eff_mean: mean,
        m_eff_std: std,
        plateau_method: method,
        chi2_dof,
        chi2_max: if plateau_cfg.mode == PlateauMode::Stat {
            Some(plateau_cfg.chi2_max)
        } else {
            None
        },
        chi2_ok,
    };

    (status, report)
}

fn analyze_gap_compatibility(
    series: &[MassEffectiveScalingPoint],
    w_min: usize,
    k_sigma: f64,
) -> GapCompatibility {
    let tested_m0 = vec![0.1, 0.2, 0.3];

    let mut per_l: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for p in series {
        let mut statuses: Vec<String> = Vec::new();
        for &m0 in tested_m0.iter() {
            statuses.push(
                classify_gap_compatibility(
                    p.m_eff_mean,
                    p.m_eff_std,
                    p.plateau_width,
                    m0,
                    w_min,
                    k_sigma,
                )
                .to_string(),
            );
        }
        per_l.insert(p.l.to_string(), statuses);
    }

    let mut global: BTreeMap<String, String> = BTreeMap::new();
    for (j, &m0) in tested_m0.iter().enumerate() {
        let mut any_incompatible = false;
        let mut all_compatible = true;
        let mut any_decisive = false;

        for p in series {
            let s = classify_gap_compatibility(
                p.m_eff_mean,
                p.m_eff_std,
                p.plateau_width,
                m0,
                w_min,
                k_sigma,
            );
            match s {
                "incompatible" => {
                    any_incompatible = true;
                    any_decisive = true;
                }
                "compatible" => {
                    any_decisive = true;
                }
                _ => {
                    all_compatible = false;
                }
            }
        }

        let g = if any_incompatible {
            "incompatible"
        } else if all_compatible && any_decisive {
            "compatible"
        } else {
            "inconclusive"
        };
        global.insert(format!("m0={}", m0), g.to_string());
        let _ = j;
    }

    GapCompatibility {
        operator: "raw".to_string(),
        tested_m0,
        per_l,
        global,
    }
}

fn analyze_gap_compatibility_smeared(
    series: &[MassEffectiveScalingPoint],
    operator_smearing: &OperatorSmearingReport,
    per_l_stats: &BTreeMap<String, OperatorSmearingResult>,
    l_max: usize,
    w_min: usize,
    k_sigma: f64,
    smeared_nonmax_fallback: bool,
) -> GapCompatibilitySmeared {
    let tested_m0 = vec![0.1, 0.2, 0.3];

    let best_steps = operator_smearing.ape.best.steps;
    let smeared = operator_smearing
        .ape
        .results
        .get(&best_steps.to_string())
        .cloned()
        .unwrap_or(OperatorSmearingResult {
            plateau_width: 0,
            m_eff_mean: 0.0,
            m_eff_std: 0.0,
            plateau_r_start: None,
            plateau_r_end: None,
            plateau_method: None,
            plateau_chi2_dof: None,
            plateau_n_points: None,
        });

    let mut per_l: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for p in series {
        let mut statuses: Vec<String> = Vec::new();
        for &m0 in tested_m0.iter() {
            let stats =
                per_l_stats
                    .get(&p.l.to_string())
                    .cloned()
                    .unwrap_or(OperatorSmearingResult {
                        plateau_width: 0,
                        m_eff_mean: 0.0,
                        m_eff_std: 0.0,
                        plateau_r_start: None,
                        plateau_r_end: None,
                        plateau_method: None,
                        plateau_chi2_dof: None,
                        plateau_n_points: None,
                    });
            let (mean, std, width) = if p.l == l_max {
                (smeared.m_eff_mean, smeared.m_eff_std, smeared.plateau_width)
            } else if smeared_nonmax_fallback {
                (stats.m_eff_mean, stats.m_eff_std, stats.plateau_width)
            } else {
                (0.0, 1.0, 0usize)
            };
            statuses
                .push(classify_gap_compatibility(mean, std, width, m0, w_min, k_sigma).to_string());
        }
        per_l.insert(p.l.to_string(), statuses);
    }

    let mut global: BTreeMap<String, String> = BTreeMap::new();
    for &m0 in tested_m0.iter() {
        let mut any_incompatible = false;
        let mut all_compatible = true;
        let mut any_decisive = false;

        for p in series {
            let stats =
                per_l_stats
                    .get(&p.l.to_string())
                    .cloned()
                    .unwrap_or(OperatorSmearingResult {
                        plateau_width: 0,
                        m_eff_mean: 0.0,
                        m_eff_std: 0.0,
                        plateau_r_start: None,
                        plateau_r_end: None,
                        plateau_method: None,
                        plateau_chi2_dof: None,
                        plateau_n_points: None,
                    });
            let (mean, std, width) = if p.l == l_max {
                (smeared.m_eff_mean, smeared.m_eff_std, smeared.plateau_width)
            } else if smeared_nonmax_fallback {
                (stats.m_eff_mean, stats.m_eff_std, stats.plateau_width)
            } else {
                (0.0, 1.0, 0usize)
            };
            let s = classify_gap_compatibility(mean, std, width, m0, w_min, k_sigma);
            match s {
                "incompatible" => {
                    any_incompatible = true;
                    any_decisive = true;
                }
                "compatible" => {
                    any_decisive = true;
                }
                _ => {
                    all_compatible = false;
                }
            }
        }

        let g = if any_incompatible {
            "incompatible"
        } else if all_compatible && any_decisive {
            "compatible"
        } else {
            "inconclusive"
        };
        global.insert(format!("m0={}", m0), g.to_string());
    }

    GapCompatibilitySmeared {
        operator: "ape".to_string(),
        steps: best_steps,
        tested_m0,
        per_l_stats: per_l_stats.clone(),
        per_l,
        global,
    }
}

fn build_operator_consistency(
    series: &[MassEffectiveScalingPoint],
    operator_smearing: &OperatorSmearingReport,
    l_max: usize,
) -> OperatorConsistencyReport {
    let raw = series.iter().find(|p| p.l == l_max).cloned();
    let best_steps = operator_smearing.ape.best.steps;
    let smeared = operator_smearing
        .ape
        .results
        .get(&best_steps.to_string())
        .cloned();

    let (m_raw, s_raw) = raw
        .map(|p| (p.m_eff_mean, p.m_eff_std))
        .unwrap_or((0.0, 0.0));
    let (m_sm, s_sm) = smeared
        .map(|r| (r.m_eff_mean, r.m_eff_std))
        .unwrap_or((0.0, 0.0));

    let delta = (m_raw - m_sm).abs();
    let denom = (s_raw * s_raw + s_sm * s_sm).sqrt().max(1e-12);
    let sigma = delta / denom;
    let consistent_2sigma = sigma <= 2.0;

    OperatorConsistencyReport {
        raw_vs_smeared: OperatorConsistencyPair {
            delta_m_eff: delta,
            sigma_units: sigma,
            consistent_2sigma,
        },
    }
}

fn classify_gap_compatibility(
    mean: f64,
    std: f64,
    plateau_width: usize,
    m0: f64,
    w_min: usize,
    k_sigma: f64,
) -> &'static str {
    if plateau_width < w_min {
        return "inconclusive";
    }
    if !(mean.is_finite() && std.is_finite()) {
        return "inconclusive";
    }
    if std <= 0.0 {
        return "inconclusive";
    }
    if mean + k_sigma * std < m0 {
        "incompatible"
    } else if mean - k_sigma * std >= m0 {
        "compatible"
    } else {
        "inconclusive"
    }
}

fn run_mill_with_rp(
    cfg: &MillRunConfig,
    mut rp: Option<&mut RpAccumulator>,
    mut mass: Option<&mut MassAccumulator>,
    mut smear: Option<&mut OperatorSmearingAccumulator>,
) -> MillRunOutput {
    let mut rng = StdRng::seed_from_u64(cfg.seed);
    let mut field = init_su2_field(cfg.l, &mut rng);

    let thermal = cfg.n_thermal_sweeps;
    for _ in 0..thermal {
        su2_sweep_update(&mut field, cfg.beta, cfg.step_size, &mut rng);
    }

    let measure_every = cfg.measure_every;
    let mut measurements: Vec<f64> = Vec::new();

    for sweep in 0..cfg.n_sweeps {
        su2_sweep_update(&mut field, cfg.beta, cfg.step_size, &mut rng);
        if (sweep + 1) % measure_every == 0 {
            let p = field.plaquette_mean();
            measurements.push(p);
            let _ = p;
            if let Some(acc) = rp.as_deref_mut() {
                acc.observe(&field);
            }
            if let Some(acc) = mass.as_deref_mut() {
                acc.observe(&field);
            }
            if let Some(acc) = smear.as_deref_mut() {
                acc.observe(&field);
            }
        }
    }

    let (plaq_mean, plaq_std) = mean_std(&measurements);

    MillRunOutput {
        trace_id: format!("MILL_SU2_3D_L{}_b{}_seed{}", cfg.l, cfg.beta, cfg.seed),
        lattice: LatticeSummary {
            dim: 3,
            l: cfg.l,
            beta: cfg.beta,
            n_links: 3 * cfg.l * cfg.l * cfg.l,
            n_plaquettes: 3 * cfg.l * cfg.l * cfg.l,
            step_size: cfg.step_size,
        },
        observables: ObservablesSummary {
            n_measurements: measurements.len(),
            plaquette_mean: plaq_mean,
            plaquette_std: plaq_std,
        },
        tests: TestsSummary {
            translation_invariance_max_row_dev: translation_invariance_row_dev(&field),
            reflection_positivity_estimate: reflection_positivity_estimate(&field),
        },
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum IcMode {
    Cold,
    Hot,
    Smooth,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum UpdateMode {
    Metropolis,
    Heatbath,
    HeatbathOverrelax,
}

fn parse_ic_mode_from_env() -> IcMode {
    match std::env::var("MILL_IC")
        .ok()
        .as_deref()
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
    {
        Some("hot") => IcMode::Hot,
        Some("smooth") => IcMode::Smooth,
        _ => IcMode::Cold,
    }
}

fn parse_update_mode_from_env() -> UpdateMode {
    match std::env::var("MILL_UPDATE")
        .ok()
        .as_deref()
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
    {
        Some("heatbath") | Some("hb") => UpdateMode::Heatbath,
        Some("hb_or") | Some("heatbath_overrelax") => UpdateMode::HeatbathOverrelax,
        _ => UpdateMode::Metropolis,
    }
}

fn init_su2_field(l: usize, rng: &mut StdRng) -> Su2Gauge3D {
    let mode = parse_ic_mode_from_env();
    let mut field = Su2Gauge3D::new(l);
    if mode == IcMode::Cold {
        return field;
    }

    for z in 0..l {
        for y in 0..l {
            for x in 0..l {
                field.set_link_x(x, y, z, su2_random_haar(rng));
                field.set_link_y(x, y, z, su2_random_haar(rng));
                field.set_link_z(x, y, z, su2_random_haar(rng));
            }
        }
    }

    if mode == IcMode::Smooth {
        // "Directed" IC: take a hot start and smooth it before thermalization.
        field = ape_smear_su2(&field, 0.5, 10);
    }
    field
}

pub fn analyze_invariance_scaling(runs: &[MillRunOutput]) -> InvarianceScaling {
    let violations: Vec<InvarianceViolationPoint> = runs
        .iter()
        .map(|r| InvarianceViolationPoint {
            l: r.lattice.l,
            value: r.tests.translation_invariance_max_row_dev,
        })
        .collect();

    let slope_estimate = linear_slope_l_vs_value(&violations);
    let trend = match slope_estimate {
        Some(s) if s < -1e-8 => "decreasing",
        Some(s) if s > 1e-8 => "increasing",
        _ => "flat",
    }
    .to_string();

    InvarianceScaling {
        violations,
        trend,
        slope_estimate,
    }
}

fn linear_slope_l_vs_value(points: &[InvarianceViolationPoint]) -> Option<f64> {
    if points.len() < 2 {
        return None;
    }
    let n = points.len() as f64;
    let xs: Vec<f64> = points.iter().map(|p| p.l as f64).collect();
    let ys: Vec<f64> = points.iter().map(|p| p.value).collect();

    let x_mean = xs.iter().sum::<f64>() / n;
    let y_mean = ys.iter().sum::<f64>() / n;

    let mut num = 0.0;
    let mut den = 0.0;
    for i in 0..points.len() {
        let dx = xs[i] - x_mean;
        let dy = ys[i] - y_mean;
        num += dx * dy;
        den += dx * dx;
    }
    if den == 0.0 {
        return None;
    }
    let slope = num / den;
    if slope.is_finite() {
        Some(slope)
    } else {
        None
    }
}

fn local_cos_sum_for_link_x(cfg: &Su2Gauge3D, x: usize, y: usize, z: usize) -> f64 {
    let ym = cfg.wrap(y as isize - 1);
    let zm = cfg.wrap(z as isize - 1);
    cfg.plaquette_cos_xy(x, y, z)
        + cfg.plaquette_cos_xy(x, ym, z)
        + cfg.plaquette_cos_xz(x, y, z)
        + cfg.plaquette_cos_xz(x, y, zm)
}

fn local_cos_sum_for_link_y(cfg: &Su2Gauge3D, x: usize, y: usize, z: usize) -> f64 {
    let xm = cfg.wrap(x as isize - 1);
    let zm = cfg.wrap(z as isize - 1);
    cfg.plaquette_cos_xy(x, y, z)
        + cfg.plaquette_cos_xy(xm, y, z)
        + cfg.plaquette_cos_yz(x, y, z)
        + cfg.plaquette_cos_yz(x, y, zm)
}

fn local_cos_sum_for_link_z(cfg: &Su2Gauge3D, x: usize, y: usize, z: usize) -> f64 {
    let xm = cfg.wrap(x as isize - 1);
    let ym = cfg.wrap(y as isize - 1);
    cfg.plaquette_cos_xz(x, y, z)
        + cfg.plaquette_cos_xz(xm, y, z)
        + cfg.plaquette_cos_yz(x, y, z)
        + cfg.plaquette_cos_yz(x, ym, z)
}

fn su2_sweep_update(cfg: &mut Su2Gauge3D, beta: f64, step_size: f64, rng: &mut StdRng) {
    match parse_update_mode_from_env() {
        UpdateMode::Metropolis => metropolis_sweep(cfg, beta, step_size, rng),
        UpdateMode::Heatbath => heatbath_sweep(cfg, beta, rng),
        UpdateMode::HeatbathOverrelax => heatbath_overrelax_sweep(cfg, beta, rng),
    }
}

fn metropolis_sweep(cfg: &mut Su2Gauge3D, beta: f64, step_size: f64, rng: &mut StdRng) {
    let l = cfg.l;
    let n_links = 3 * l * l * l;
    for _ in 0..n_links {
        let dir = rng.gen_range(0..3);
        let x = rng.gen_range(0..l);
        let y = rng.gen_range(0..l);
        let z = rng.gen_range(0..l);
        let delta = if step_size.is_finite() && step_size > 0.0 {
            rng.gen_range(-step_size..step_size)
        } else {
            0.0
        };
        let r = su2_random_near_identity(delta, rng);

        if dir == 0 {
            let old = cfg.link_x(x, y, z);
            let old_sum = local_cos_sum_for_link_x(cfg, x, y, z);
            cfg.set_link_x(x, y, z, r.mul(old).projected());
            let new_sum = local_cos_sum_for_link_x(cfg, x, y, z);
            let delta_s = -beta * (new_sum - old_sum);
            let accept = delta_s <= 0.0 || rng.gen::<f64>() < (-delta_s).exp();
            if !accept {
                cfg.set_link_x(x, y, z, old);
            }
        } else if dir == 1 {
            let old = cfg.link_y(x, y, z);
            let old_sum = local_cos_sum_for_link_y(cfg, x, y, z);
            cfg.set_link_y(x, y, z, r.mul(old).projected());
            let new_sum = local_cos_sum_for_link_y(cfg, x, y, z);
            let delta_s = -beta * (new_sum - old_sum);
            let accept = delta_s <= 0.0 || rng.gen::<f64>() < (-delta_s).exp();
            if !accept {
                cfg.set_link_y(x, y, z, old);
            }
        } else {
            let old = cfg.link_z(x, y, z);
            let old_sum = local_cos_sum_for_link_z(cfg, x, y, z);
            cfg.set_link_z(x, y, z, r.mul(old).projected());
            let new_sum = local_cos_sum_for_link_z(cfg, x, y, z);
            let delta_s = -beta * (new_sum - old_sum);
            let accept = delta_s <= 0.0 || rng.gen::<f64>() < (-delta_s).exp();
            if !accept {
                cfg.set_link_z(x, y, z, old);
            }
        }
    }
}

fn staple_sum_for_link_x(cfg: &Su2Gauge3D, x: usize, y: usize, z: usize) -> Su2 {
    let xp = cfg.wrap(x as isize + 1);
    let yp = cfg.wrap(y as isize + 1);
    let ym = cfg.wrap(y as isize - 1);
    let zp = cfg.wrap(z as isize + 1);
    let zm = cfg.wrap(z as isize - 1);

    let s_xy_f = cfg
        .link_y(xp, y, z)
        .mul(cfg.link_x(x, yp, z).dagger())
        .mul(cfg.link_y(x, y, z).dagger());
    let s_xy_b = cfg
        .link_y(xp, ym, z)
        .dagger()
        .mul(cfg.link_x(x, ym, z).dagger())
        .mul(cfg.link_y(x, ym, z));

    let s_xz_f = cfg
        .link_z(xp, y, z)
        .mul(cfg.link_x(x, y, zp).dagger())
        .mul(cfg.link_z(x, y, z).dagger());
    let s_xz_b = cfg
        .link_z(xp, y, zm)
        .dagger()
        .mul(cfg.link_x(x, y, zm).dagger())
        .mul(cfg.link_z(x, y, zm));

    s_xy_f.add(s_xy_b).add(s_xz_f).add(s_xz_b)
}

fn staple_sum_for_link_y(cfg: &Su2Gauge3D, x: usize, y: usize, z: usize) -> Su2 {
    let xp = cfg.wrap(x as isize + 1);
    let xm = cfg.wrap(x as isize - 1);
    let yp = cfg.wrap(y as isize + 1);
    let zp = cfg.wrap(z as isize + 1);
    let zm = cfg.wrap(z as isize - 1);

    let s_xy_f = cfg
        .link_x(x, yp, z)
        .mul(cfg.link_y(xp, y, z).dagger())
        .mul(cfg.link_x(x, y, z).dagger());
    let s_xy_b = cfg
        .link_x(xm, yp, z)
        .dagger()
        .mul(cfg.link_y(xm, y, z).dagger())
        .mul(cfg.link_x(xm, y, z));

    let s_yz_f = cfg
        .link_z(x, yp, z)
        .mul(cfg.link_y(x, y, zp).dagger())
        .mul(cfg.link_z(x, y, z).dagger());
    let s_yz_b = cfg
        .link_z(x, yp, zm)
        .dagger()
        .mul(cfg.link_y(x, y, zm).dagger())
        .mul(cfg.link_z(x, y, zm));

    s_xy_f.add(s_xy_b).add(s_yz_f).add(s_yz_b)
}

fn staple_sum_for_link_z(cfg: &Su2Gauge3D, x: usize, y: usize, z: usize) -> Su2 {
    let xp = cfg.wrap(x as isize + 1);
    let xm = cfg.wrap(x as isize - 1);
    let yp = cfg.wrap(y as isize + 1);
    let ym = cfg.wrap(y as isize - 1);
    let zp = cfg.wrap(z as isize + 1);

    let s_xz_f = cfg
        .link_x(x, y, zp)
        .mul(cfg.link_z(xp, y, z).dagger())
        .mul(cfg.link_x(x, y, z).dagger());
    let s_xz_b = cfg
        .link_x(xm, y, zp)
        .dagger()
        .mul(cfg.link_z(xm, y, z).dagger())
        .mul(cfg.link_x(xm, y, z));

    let s_yz_f = cfg
        .link_y(x, y, zp)
        .mul(cfg.link_z(x, yp, z).dagger())
        .mul(cfg.link_y(x, y, z).dagger());
    let s_yz_b = cfg
        .link_y(x, ym, zp)
        .dagger()
        .mul(cfg.link_z(x, ym, z).dagger())
        .mul(cfg.link_y(x, ym, z));

    s_xz_f.add(s_xz_b).add(s_yz_f).add(s_yz_b)
}

fn su2_modified_normal(param_exp: f64, rng: &mut StdRng) -> f64 {
    let r0 = rng.gen::<f64>().max(1e-12);
    let r1 = rng.gen::<f64>();
    let r2 = rng.gen::<f64>().max(1e-12);
    let c = (TAU * r1).cos();
    let v = -((r0.ln() + c * c * r2.ln()) / (2.0 * param_exp));
    v.max(0.0).sqrt()
}

fn su2_heatbath_sample_x0(param_exp: f64, rng: &mut StdRng) -> f64 {
    if !(param_exp.is_finite() && param_exp > 0.0) {
        return 1.0;
    }
    loop {
        let r = rng.gen::<f64>();
        let lambda = su2_modified_normal(param_exp, rng);
        if r * r <= 1.0 - lambda * lambda {
            return (1.0 - 2.0 * lambda * lambda).clamp(-1.0, 1.0);
        }
    }
}

fn su2_heatbath_sample(param_exp: f64, rng: &mut StdRng) -> Su2 {
    let x0 = su2_heatbath_sample_x0(param_exp, rng);
    let s = (1.0 - x0 * x0).max(0.0).sqrt();
    let mut x1 = rand_std_normal(rng);
    let mut x2 = rand_std_normal(rng);
    let mut x3 = rand_std_normal(rng);
    let mut n2 = x1 * x1 + x2 * x2 + x3 * x3;
    while !(n2.is_finite() && n2 > 1e-12) {
        x1 = rand_std_normal(rng);
        x2 = rand_std_normal(rng);
        x3 = rand_std_normal(rng);
        n2 = x1 * x1 + x2 * x2 + x3 * x3;
    }
    let inv = 1.0 / n2.sqrt();
    Su2 {
        a0: x0,
        a1: s * x1 * inv,
        a2: s * x2 * inv,
        a3: s * x3 * inv,
    }
    .projected()
}

fn su2_heatbath_update(staple_sum: Su2, beta: f64, rng: &mut StdRng) -> Su2 {
    let k2 = staple_sum.norm2();
    if !(k2.is_finite() && k2.is_normal() && k2 > 0.0) {
        return su2_random_haar(rng);
    }
    let k = k2.sqrt();
    let v_dag = staple_sum.dagger().scale(1.0 / k);
    let r = su2_heatbath_sample(beta * k, rng);
    r.mul(v_dag).projected()
}

fn su2_overrelax_update(u: Su2, staple_sum: Su2) -> Su2 {
    let v = staple_sum.projected();
    let v_dag = v.dagger();
    v_dag.mul(u.dagger()).mul(v_dag).projected()
}

fn heatbath_sweep(cfg: &mut Su2Gauge3D, beta: f64, rng: &mut StdRng) {
    let l = cfg.l;
    let n_links = 3 * l * l * l;
    for _ in 0..n_links {
        let dir = rng.gen_range(0..3);
        let x = rng.gen_range(0..l);
        let y = rng.gen_range(0..l);
        let z = rng.gen_range(0..l);
        if dir == 0 {
            let staple = staple_sum_for_link_x(cfg, x, y, z);
            cfg.set_link_x(x, y, z, su2_heatbath_update(staple, beta, rng));
        } else if dir == 1 {
            let staple = staple_sum_for_link_y(cfg, x, y, z);
            cfg.set_link_y(x, y, z, su2_heatbath_update(staple, beta, rng));
        } else {
            let staple = staple_sum_for_link_z(cfg, x, y, z);
            cfg.set_link_z(x, y, z, su2_heatbath_update(staple, beta, rng));
        }
    }
}

fn overrelax_sweep(cfg: &mut Su2Gauge3D, rng: &mut StdRng) {
    let l = cfg.l;
    let n_links = 3 * l * l * l;
    for _ in 0..n_links {
        let dir = rng.gen_range(0..3);
        let x = rng.gen_range(0..l);
        let y = rng.gen_range(0..l);
        let z = rng.gen_range(0..l);
        if dir == 0 {
            let old = cfg.link_x(x, y, z);
            let staple = staple_sum_for_link_x(cfg, x, y, z);
            cfg.set_link_x(x, y, z, su2_overrelax_update(old, staple));
        } else if dir == 1 {
            let old = cfg.link_y(x, y, z);
            let staple = staple_sum_for_link_y(cfg, x, y, z);
            cfg.set_link_y(x, y, z, su2_overrelax_update(old, staple));
        } else {
            let old = cfg.link_z(x, y, z);
            let staple = staple_sum_for_link_z(cfg, x, y, z);
            cfg.set_link_z(x, y, z, su2_overrelax_update(old, staple));
        }
    }
}

fn heatbath_overrelax_sweep(cfg: &mut Su2Gauge3D, beta: f64, rng: &mut StdRng) {
    heatbath_sweep(cfg, beta, rng);
    overrelax_sweep(cfg, rng);
}

fn su2_random_near_identity(step_size: f64, rng: &mut StdRng) -> Su2 {
    // Isotropic (Haar) SU(2) proposal, truncated near identity:
    // draw q ~ Haar(SU(2)) via uniform S^3, then take q^(t) with t=step_size/pi.
    if !(step_size.is_finite() && step_size > 0.0) {
        return Su2::identity();
    }
    let q = su2_random_haar(rng);
    let t = (step_size / PI).min(1.0);
    su2_pow(q, t)
}

fn su2_random_haar(rng: &mut StdRng) -> Su2 {
    // Uniform on S^3 via normalized 4D Gaussian -> Haar on SU(2).
    let x0 = rand_std_normal(rng);
    let x1 = rand_std_normal(rng);
    let x2 = rand_std_normal(rng);
    let x3 = rand_std_normal(rng);
    Su2 {
        a0: x0,
        a1: x1,
        a2: x2,
        a3: x3,
    }
    .projected()
}

fn rand_std_normal(rng: &mut StdRng) -> f64 {
    // Box-Muller transform (one sample).
    let u1: f64 = rng.gen::<f64>().max(1e-12);
    let u2: f64 = rng.gen::<f64>();
    (-2.0 * u1.ln()).sqrt() * (TAU * u2).cos()
}

fn su2_pow(u: Su2, t: f64) -> Su2 {
    if !(t.is_finite() && t >= 0.0) {
        return Su2::identity();
    }
    let u = u.projected();
    let a0 = u.a0.clamp(-1.0, 1.0);
    let phi = a0.acos();
    let sin_phi = (1.0 - a0 * a0).max(0.0).sqrt();
    if sin_phi < 1e-12 {
        return Su2::identity();
    }
    let nx = u.a1 / sin_phi;
    let ny = u.a2 / sin_phi;
    let nz = u.a3 / sin_phi;
    let phi2 = (t * phi).min(PI);
    let c = phi2.cos();
    let s = phi2.sin();
    Su2 {
        a0: c,
        a1: nx * s,
        a2: ny * s,
        a3: nz * s,
    }
    .projected()
}

fn mean_std(values: &[f64]) -> (f64, f64) {
    if values.is_empty() {
        return (f64::NAN, f64::NAN);
    }
    let n = values.len() as f64;
    let mean = values.iter().sum::<f64>() / n;
    let var = values.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / n;
    (mean, var.sqrt())
}

fn translation_invariance_row_dev(cfg: &Su2Gauge3D) -> f64 {
    let global = cfg.plaquette_mean();
    let rows = cfg.plaquette_mean_by_row();
    rows.into_iter()
        .map(|m| (m - global).abs())
        .fold(0.0, f64::max)
}

fn reflection_positivity_estimate(cfg: &Su2Gauge3D) -> f64 {
    let l = cfg.l;
    if l < 2 {
        return f64::NAN;
    }
    let mean_p = cfg.plaquette_mean();

    let y0 = 0usize;
    let y1 = l - 1 - y0;

    let mut f0 = 0.0;
    let mut f1 = 0.0;
    for z in 0..l {
        for x in 0..l {
            f0 += cfg.plaquette_cos(x, y0, z) - mean_p;
            f1 += cfg.plaquette_cos(x, y1, z) - mean_p;
        }
    }
    f0 * f1
}

pub fn run_mill(cfg: MillRunConfig) -> MillRunOutput {
    run_mill_with_rp(&cfg, None, None, None)
}

struct OperatorSmearingAccumulator {
    alpha: f64,
    steps: Vec<usize>,
    accs: Vec<MassAccumulator>,
}

impl OperatorSmearingAccumulator {
    fn new(l: usize, expected_n_measurements: usize) -> Self {
        let alpha = 0.5;
        let steps = vec![0usize, 1, 2, 3];
        let mut accs: Vec<MassAccumulator> = Vec::new();
        for _ in steps.iter() {
            accs.push(MassAccumulator::new(l, expected_n_measurements));
        }
        Self { alpha, steps, accs }
    }

    fn observe(&mut self, field: &Su2Gauge3D) {
        for (i, &k) in self.steps.iter().enumerate() {
            if k == 0 {
                self.accs[i].observe(field);
            } else {
                let smeared = ape_smear_su2(field, self.alpha, k);
                self.accs[i].observe(&smeared);
            }
        }
    }
}

fn build_operator_smearing_report(
    acc: Option<OperatorSmearingAccumulator>,
    plateau_cfg: PlateauCfg,
) -> OperatorSmearingReport {
    let alpha = 0.5;
    let steps = vec![0usize, 1, 2, 3];
    let mut results: BTreeMap<String, OperatorSmearingResult> = BTreeMap::new();
    let mut best_steps = 0usize;
    let mut best_width = 0usize;
    let mut best_std = f64::INFINITY;

    if let Some(a) = acc {
        for (i, &k) in a.steps.iter().enumerate() {
            let plateau = mass_scaling_plateau_from_acc_with_cfg(&a.accs[i], plateau_cfg);
            let (width, std) = if let Some(ref p) = plateau {
                (plateau_width_from_r_range(p.r_start, p.r_end), p.m_eff_std)
            } else {
                (0usize, 0.0f64)
            };
            let r = operator_smearing_result_from_plateau(plateau);
            results.insert(k.to_string(), r);
            if width > best_width || (width == best_width && std < best_std) {
                best_width = width;
                best_std = std;
                best_steps = k;
            }
        }
    } else {
        for &k in steps.iter() {
            results.insert(
                k.to_string(),
                OperatorSmearingResult {
                    plateau_width: 0,
                    m_eff_mean: 0.0,
                    m_eff_std: 0.0,
                    plateau_r_start: None,
                    plateau_r_end: None,
                    plateau_method: None,
                    plateau_chi2_dof: None,
                    plateau_n_points: None,
                },
            );
        }
    }

    OperatorSmearingReport {
        ape: ApeSmearingReport {
            alpha,
            steps,
            results,
            best: OperatorSmearingBest {
                steps: best_steps,
                criterion: "max_plateau_width".to_string(),
            },
        },
    }
}

fn smearing_result_for_steps(
    acc: &OperatorSmearingAccumulator,
    steps: usize,
    plateau_cfg: PlateauCfg,
) -> OperatorSmearingResult {
    let Some(i) = acc.steps.iter().position(|&k| k == steps) else {
        return OperatorSmearingResult {
            plateau_width: 0,
            m_eff_mean: 0.0,
            m_eff_std: 0.0,
            plateau_r_start: None,
            plateau_r_end: None,
            plateau_method: None,
            plateau_chi2_dof: None,
            plateau_n_points: None,
        };
    };
    let plateau = mass_scaling_plateau_from_acc_with_cfg(&acc.accs[i], plateau_cfg);
    operator_smearing_result_from_plateau(plateau)
}

fn operator_smearing_result_from_plateau(plateau: Option<MassPlateau>) -> OperatorSmearingResult {
    let Some(p) = plateau else {
        return OperatorSmearingResult {
            plateau_width: 0,
            m_eff_mean: 0.0,
            m_eff_std: 0.0,
            plateau_r_start: None,
            plateau_r_end: None,
            plateau_method: None,
            plateau_chi2_dof: None,
            plateau_n_points: None,
        };
    };
    OperatorSmearingResult {
        plateau_width: plateau_width_from_r_range(p.r_start, p.r_end),
        m_eff_mean: p.m_eff_mean,
        m_eff_std: p.m_eff_std,
        plateau_r_start: Some(p.r_start),
        plateau_r_end: Some(p.r_end),
        plateau_method: p.method,
        plateau_chi2_dof: p.chi2_dof,
        plateau_n_points: p.n_points,
    }
}

fn ape_smear_su2(field: &Su2Gauge3D, alpha: f64, steps: usize) -> Su2Gauge3D {
    if steps == 0 {
        return field.clone();
    }
    let mut cur = field.clone();
    let l = field.l;
    for _ in 0..steps {
        let mut next = cur.clone();
        for z in 0..l {
            for y in 0..l {
                for x in 0..l {
                    let xp = cur.wrap(x as isize + 1);
                    let xm = cur.wrap(x as isize - 1);
                    let yp = cur.wrap(y as isize + 1);
                    let ym = cur.wrap(y as isize - 1);
                    let zp = cur.wrap(z as isize + 1);
                    let zm = cur.wrap(z as isize - 1);

                    let ux = cur.link_x(x, y, z);
                    let staple_y_f = cur
                        .link_y(x, y, z)
                        .mul(cur.link_x(x, yp, z))
                        .mul(cur.link_y(xp, y, z).dagger());
                    let staple_y_b = cur
                        .link_y(x, ym, z)
                        .dagger()
                        .mul(cur.link_x(x, ym, z))
                        .mul(cur.link_y(xp, ym, z));
                    let staple_z_f = cur
                        .link_z(x, y, z)
                        .mul(cur.link_x(x, y, zp))
                        .mul(cur.link_z(xp, y, z).dagger());
                    let staple_z_b = cur
                        .link_z(x, y, zm)
                        .dagger()
                        .mul(cur.link_x(x, y, zm))
                        .mul(cur.link_z(xp, y, zm));
                    let staples = staple_y_f
                        .add(staple_y_b)
                        .add(staple_z_f)
                        .add(staple_z_b);
                    let ux_new = ux
                        .scale(1.0 - alpha)
                        .add(staples.scale(alpha * 0.25))
                        .projected();
                    next.set_link_x(x, y, z, ux_new);

                    let uy = cur.link_y(x, y, z);
                    let staple_x_f = cur
                        .link_x(x, y, z)
                        .mul(cur.link_y(xp, y, z))
                        .mul(cur.link_x(x, yp, z).dagger());
                    let staple_x_b = cur
                        .link_x(xm, y, z)
                        .dagger()
                        .mul(cur.link_y(xm, y, z))
                        .mul(cur.link_x(xm, yp, z));
                    let staple_z_f = cur
                        .link_z(x, y, z)
                        .mul(cur.link_y(x, y, zp))
                        .mul(cur.link_z(x, yp, z).dagger());
                    let staple_z_b = cur
                        .link_z(x, y, zm)
                        .dagger()
                        .mul(cur.link_y(x, y, zm))
                        .mul(cur.link_z(x, yp, zm));
                    let staples = staple_x_f
                        .add(staple_x_b)
                        .add(staple_z_f)
                        .add(staple_z_b);
                    let uy_new = uy
                        .scale(1.0 - alpha)
                        .add(staples.scale(alpha * 0.25))
                        .projected();
                    next.set_link_y(x, y, z, uy_new);

                    let uz = cur.link_z(x, y, z);
                    let staple_x_f = cur
                        .link_x(x, y, z)
                        .mul(cur.link_z(xp, y, z))
                        .mul(cur.link_x(x, y, zp).dagger());
                    let staple_x_b = cur
                        .link_x(xm, y, z)
                        .dagger()
                        .mul(cur.link_z(xm, y, z))
                        .mul(cur.link_x(xm, y, zp));
                    let staple_y_f = cur
                        .link_y(x, y, z)
                        .mul(cur.link_z(x, yp, z))
                        .mul(cur.link_y(x, y, zp).dagger());
                    let staple_y_b = cur
                        .link_y(x, ym, z)
                        .dagger()
                        .mul(cur.link_z(x, ym, z))
                        .mul(cur.link_y(x, ym, zp));
                    let staples = staple_x_f
                        .add(staple_x_b)
                        .add(staple_y_f)
                        .add(staple_y_b);
                    let uz_new = uz
                        .scale(1.0 - alpha)
                        .add(staples.scale(alpha * 0.25))
                        .projected();
                    next.set_link_z(x, y, z, uz_new);
                }
            }
        }
        cur = next;
    }
    cur
}

#[derive(Clone, Debug)]
struct MassBlock {
    n: usize,
    sum_p: f64,
    sum_pp: Vec<f64>,
}

struct MassAccumulator {
    l: usize,
    r_max: usize,
    count: usize,
    sum_p: f64,
    sum_pp: Vec<f64>,
    block_size: usize,
    cur_n: usize,
    cur_sum_p: f64,
    cur_sum_pp: Vec<f64>,
    blocks: Vec<MassBlock>,
}

impl MassAccumulator {
    fn new(l: usize, expected_n_measurements: usize) -> Self {
        let r_max = (l / 4).max(1);
        let mut block_size = if JK_TARGET_N_BLOCKS > 0 {
            (expected_n_measurements / JK_TARGET_N_BLOCKS).max(1)
        } else {
            expected_n_measurements.max(1)
        };
        if expected_n_measurements > 0 && expected_n_measurements / block_size < 2 {
            block_size = expected_n_measurements.max(1);
        }

        Self {
            l,
            r_max,
            count: 0,
            sum_p: 0.0,
            sum_pp: vec![0.0; r_max],
            block_size,
            cur_n: 0,
            cur_sum_p: 0.0,
            cur_sum_pp: vec![0.0; r_max],
            blocks: Vec::new(),
        }
    }

    fn observe(&mut self, field: &Su2Gauge3D) {
        let l = self.l;
        let mut pgrid: Vec<f64> = vec![0.0; l * l * l];
        let mut s = 0.0;
        for z in 0..l {
            for y in 0..l {
                for x in 0..l {
                    let v = field.plaquette_cos(x, y, z);
                    pgrid[x + l * (y + l * z)] = v;
                    s += v;
                }
            }
        }
        let mean_p = s / ((l * l * l) as f64);

        let mut mean_pp: Vec<f64> = vec![0.0; self.r_max];
        for r in 1..=self.r_max {
            let mut spp = 0.0;
            for z in 0..l {
                for y in 0..l {
                    let y2 = (y + r) % l;
                    for x in 0..l {
                        let a = pgrid[x + l * (y + l * z)];
                        let b = pgrid[x + l * (y2 + l * z)];
                        spp += a * b;
                    }
                }
            }
            mean_pp[r - 1] = spp / ((l * l * l) as f64);
        }

        self.observe_stats(mean_p, &mean_pp);
    }

    fn observe_stats(&mut self, mean_p: f64, mean_pp: &[f64]) {
        debug_assert_eq!(mean_pp.len(), self.r_max);

        self.sum_p += mean_p;
        for i in 0..self.r_max {
            self.sum_pp[i] += mean_pp[i];
        }
        self.count += 1;

        self.cur_sum_p += mean_p;
        for i in 0..self.r_max {
            self.cur_sum_pp[i] += mean_pp[i];
        }
        self.cur_n += 1;
        if self.cur_n >= self.block_size {
            let sum_pp = std::mem::take(&mut self.cur_sum_pp);
            let block = MassBlock {
                n: self.cur_n,
                sum_p: self.cur_sum_p,
                sum_pp,
            };
            self.blocks.push(block);
            self.cur_n = 0;
            self.cur_sum_p = 0.0;
            self.cur_sum_pp = vec![0.0; self.r_max];
        }
    }

    fn correlator_and_m_eff_from_sums(
        &self,
        count: usize,
        sum_p: f64,
        sum_pp: &[f64],
    ) -> (Vec<f64>, Vec<f64>) {
        if count == 0 {
            return (vec![0.0; self.r_max], Vec::new());
        }
        let denom = count as f64;
        let mean_p = sum_p / denom;
        let mut correlator = vec![0.0; self.r_max];
        for i in 0..self.r_max {
            correlator[i] = sum_pp[i] / denom - mean_p * mean_p;
        }

        let eps = 1e-12;
        let mut m_eff: Vec<f64> = Vec::new();
        if correlator.len() >= 2 {
            for i in 0..(correlator.len() - 1) {
                let c0 = correlator[i].abs().max(eps);
                let c1 = correlator[i + 1].abs().max(eps);
                m_eff.push((c0 / c1).ln());
            }
        }
        (correlator, m_eff)
    }

    fn plateau_point(&self, cfg: PlateauCfg) -> Option<MassPlateau> {
        if self.count == 0 {
            return None;
        }
        let (_correlator, m_eff) =
            self.correlator_and_m_eff_from_sums(self.count, self.sum_p, &self.sum_pp);
        if m_eff.len() < 2 {
            return Some(MassPlateau {
                r_start: 1,
                r_end: 2,
                m_eff_mean: 0.0,
                m_eff_std: 0.0,
                method: None,
                chi2_dof: None,
                n_points: None,
            });
        }

        let rel = if cfg.rel_thresh.is_finite() && cfg.rel_thresh > 0.0 {
            cfg.rel_thresh
        } else {
            0.05
        };
        let plateau_geom = find_mass_plateau(&m_eff, rel);

        if cfg.mode != PlateauMode::Stat || self.blocks.len() < 2 {
            return Some(self.plateau_with_jk_std(plateau_geom, &m_eff, None));
        }

        let k = if cfg.k.is_finite() && cfg.k > 0.0 {
            cfg.k
        } else {
            2.0
        };
        let chi2_max = if cfg.chi2_max.is_finite() && cfg.chi2_max > 0.0 {
            cfg.chi2_max
        } else {
            2.0
        };

        let sigma = self.jackknife_sigma_m_eff_points()?;
        if sigma.len() != m_eff.len() {
            return Some(self.plateau_with_jk_std(plateau_geom, &m_eff, None));
        }

        let Some(stat) = find_mass_plateau_stat(&m_eff, &sigma, k, chi2_max) else {
            let mut p = self.plateau_with_jk_std(plateau_geom, &m_eff, None);
            p.method = Some("stat_fallback_legacy".to_string());
            return Some(p);
        };

        let weights = stat.weights;
        let mut plateau = MassPlateau {
            r_start: stat.start + 1,
            r_end: stat.end + 1,
            m_eff_mean: stat.mu,
            m_eff_std: 0.0,
            method: Some("stat".to_string()),
            chi2_dof: Some(stat.chi2_dof),
            n_points: Some(stat.n_points),
        };
        plateau = self.plateau_with_jk_std(plateau, &m_eff, Some(&weights));
        Some(plateau)
    }

    fn plateau_with_jk_std(
        &self,
        mut plateau: MassPlateau,
        m_eff_full: &[f64],
        weights: Option<&[f64]>,
    ) -> MassPlateau {
        let start_idx = plateau.r_start.saturating_sub(1);
        let end_idx = plateau.r_end.saturating_sub(1);
        if !(start_idx <= end_idx && end_idx < m_eff_full.len()) {
            plateau.m_eff_std = 0.0;
            return plateau;
        }

        if self.blocks.len() < 2 {
            return plateau;
        }

        let mut thetas: Vec<f64> = Vec::new();
        for b in self.blocks.iter() {
            let Some(m_eff_loo) = self.m_eff_leave_one_block(b) else {
                continue;
            };
            if !(start_idx <= end_idx && end_idx < m_eff_loo.len()) {
                continue;
            }
            let window = &m_eff_loo[start_idx..=end_idx];
            let theta = match weights {
                Some(w) => weighted_mean_fixed_weights(window, w),
                None => {
                    let (m, _) = mean_std(window);
                    Some(m)
                }
            };
            if let Some(t) = theta {
                if t.is_finite() {
                    thetas.push(t);
                }
            }
        }

        if thetas.len() >= 2 {
            plateau.m_eff_std = jackknife_std_from_leave_one_out(&thetas);
        }
        plateau
    }

    fn m_eff_leave_one_block(&self, b: &MassBlock) -> Option<Vec<f64>> {
        let n_total = self.count;
        if b.n >= n_total {
            return None;
        }
        let n_loo = n_total - b.n;
        if n_loo == 0 {
            return None;
        }
        let mean_p = (self.sum_p - b.sum_p) / (n_loo as f64);
        let mut mean_pp = vec![0.0; self.r_max];
        for i in 0..self.r_max {
            mean_pp[i] = (self.sum_pp[i] - b.sum_pp[i]) / (n_loo as f64);
        }

        let mut corr_loo = vec![0.0; self.r_max];
        for i in 0..self.r_max {
            corr_loo[i] = mean_pp[i] - mean_p * mean_p;
        }
        let eps = 1e-12;
        let mut m_eff_loo: Vec<f64> = Vec::new();
        if corr_loo.len() >= 2 {
            for i in 0..(corr_loo.len() - 1) {
                let c0 = corr_loo[i].abs().max(eps);
                let c1 = corr_loo[i + 1].abs().max(eps);
                m_eff_loo.push((c0 / c1).ln());
            }
        }
        Some(m_eff_loo)
    }

    fn jackknife_sigma_m_eff_points(&self) -> Option<Vec<f64>> {
        if self.blocks.len() < 2 {
            return None;
        }
        let (_corr, m_eff) =
            self.correlator_and_m_eff_from_sums(self.count, self.sum_p, &self.sum_pp);
        if m_eff.is_empty() {
            return None;
        }

        let mut per_point: Vec<Vec<f64>> = vec![Vec::new(); m_eff.len()];
        for b in self.blocks.iter() {
            let Some(m_eff_loo) = self.m_eff_leave_one_block(b) else {
                continue;
            };
            if m_eff_loo.len() != m_eff.len() {
                continue;
            }
            for i in 0..m_eff.len() {
                let v = m_eff_loo[i];
                if v.is_finite() {
                    per_point[i].push(v);
                }
            }
        }

        let mut sigma = vec![0.0; m_eff.len()];
        for i in 0..m_eff.len() {
            if per_point[i].len() < 2 {
                sigma[i] = 0.0;
            } else {
                sigma[i] = jackknife_std_from_leave_one_out(&per_point[i]);
            }
        }
        Some(sigma)
    }
}

fn build_mass_effective_report(
    l: usize,
    acc: Option<MassAccumulator>,
    plateau_cfg: PlateauCfg,
) -> MassEffectiveReport {
    let mut r_max = (l / 4).max(1);
    let mut correlator = vec![0.0; r_max];
    let mut m_eff: Vec<f64> = Vec::new();
    if let Some(a) = acc {
        r_max = a.r_max;
        let (corr, me) = a.correlator_and_m_eff_from_sums(a.count, a.sum_p, &a.sum_pp);
        correlator = corr;
        m_eff = me;
        let plateau = a
            .plateau_point(plateau_cfg)
            .unwrap_or_else(|| find_mass_plateau(&m_eff, plateau_cfg.rel_thresh));

        return MassEffectiveReport {
            l,
            r_max,
            correlator,
            m_eff,
            plateau,
        };
    }

    let plateau = find_mass_plateau(&m_eff, plateau_cfg.rel_thresh);

    MassEffectiveReport {
        l,
        r_max,
        correlator,
        m_eff,
        plateau,
    }
}

fn find_mass_plateau(m_eff: &[f64], rel_thresh: f64) -> MassPlateau {
    if m_eff.len() < 2 {
        return MassPlateau {
            r_start: 1,
            r_end: 2,
            m_eff_mean: 0.0,
            m_eff_std: 0.0,
            method: None,
            chi2_dof: None,
            n_points: None,
        };
    }

    let mut best_start = 0usize;
    let mut best_end = 1usize;
    let mut cur_start = 0usize;

    for i in 0..(m_eff.len() - 1) {
        let a = m_eff[i];
        let b = m_eff[i + 1];
        let denom = a.abs().max(1e-12);
        let rel = (b - a).abs() / denom;
        if rel <= rel_thresh {
            let cur_end = i + 1;
            if (cur_end - cur_start) > (best_end - best_start) {
                best_start = cur_start;
                best_end = cur_end;
            }
        } else {
            cur_start = i + 1;
        }
    }

    if best_end <= best_start {
        best_start = 0;
        best_end = 1;
    }

    let slice = &m_eff[best_start..=best_end];
    let (mean, std) = mean_std(slice);

    MassPlateau {
        r_start: best_start + 1,
        r_end: best_end + 1,
        m_eff_mean: mean,
        m_eff_std: std,
        method: None,
        chi2_dof: None,
        n_points: None,
    }
}

#[derive(Clone, Debug)]
struct StatPlateau {
    start: usize,
    end: usize,
    mu: f64,
    chi2_dof: f64,
    n_points: usize,
    weights: Vec<f64>,
}

fn weighted_mean_fixed_weights(values: &[f64], weights: &[f64]) -> Option<f64> {
    if values.len() != weights.len() || values.is_empty() {
        return None;
    }
    let mut num = 0.0;
    let mut den = 0.0;
    for (v, w) in values.iter().copied().zip(weights.iter().copied()) {
        if !(v.is_finite() && w.is_finite() && w > 0.0) {
            return None;
        }
        num += w * v;
        den += w;
    }
    if den <= 0.0 || !den.is_finite() {
        return None;
    }
    let mu = num / den;
    if mu.is_finite() {
        Some(mu)
    } else {
        None
    }
}

fn find_mass_plateau_stat(
    m_eff: &[f64],
    sigma: &[f64],
    k: f64,
    chi2_max: f64,
) -> Option<StatPlateau> {
    if m_eff.len() < 2 || m_eff.len() != sigma.len() {
        return None;
    }
    if !(k.is_finite() && k > 0.0) {
        return None;
    }
    if !(chi2_max.is_finite() && chi2_max > 0.0) {
        return None;
    }

    let n = m_eff.len();
    let mut best: Option<StatPlateau> = None;

    for start in 0..(n - 1) {
        if !(m_eff[start].is_finite() && sigma[start].is_finite() && sigma[start] > 0.0) {
            continue;
        }
        for end in (start + 1)..n {
            if !(m_eff[end].is_finite() && sigma[end].is_finite() && sigma[end] > 0.0) {
                break;
            }

            let mut ok = true;
            for i in start..end {
                if !(m_eff[i].is_finite()
                    && m_eff[i + 1].is_finite()
                    && sigma[i].is_finite()
                    && sigma[i + 1].is_finite()
                    && sigma[i] > 0.0
                    && sigma[i + 1] > 0.0)
                {
                    ok = false;
                    break;
                }
                let d = (m_eff[i + 1] - m_eff[i]).abs();
                let tol = k * (sigma[i] * sigma[i] + sigma[i + 1] * sigma[i + 1]).sqrt();
                if !(tol.is_finite() && d <= tol) {
                    ok = false;
                    break;
                }
            }
            if !ok {
                break;
            }

            let window = &m_eff[start..=end];
            let sigw = &sigma[start..=end];
            let mut weights = Vec::with_capacity(sigw.len());
            for &s in sigw {
                if !(s.is_finite() && s > 0.0) {
                    ok = false;
                    break;
                }
                weights.push(1.0 / (s * s));
            }
            if !ok {
                break;
            }
            let Some(mu) = weighted_mean_fixed_weights(window, &weights) else {
                break;
            };

            let mut chi2 = 0.0;
            for (&v, &w) in window.iter().zip(weights.iter()) {
                let r = v - mu;
                chi2 += w * r * r;
            }
            let dof = (window.len() - 1) as f64;
            if dof <= 0.0 {
                continue;
            }
            let chi2_dof = chi2 / dof;
            if !(chi2_dof.is_finite() && chi2_dof <= chi2_max) {
                continue;
            }

            let cand = StatPlateau {
                start,
                end,
                mu,
                chi2_dof,
                n_points: window.len(),
                weights,
            };
            best = match best {
                None => Some(cand),
                Some(prev) => {
                    if cand.n_points > prev.n_points
                        || (cand.n_points == prev.n_points && cand.chi2_dof < prev.chi2_dof)
                    {
                        Some(cand)
                    } else {
                        Some(prev)
                    }
                }
            };
        }
    }

    best
}

fn jackknife_std_from_leave_one_out(theta_loo: &[f64]) -> f64 {
    let g = theta_loo.len();
    if g < 2 {
        return 0.0;
    }
    let g_f = g as f64;

    let theta_bar = theta_loo.iter().sum::<f64>() / g_f;
    let mut s = 0.0;
    for &t in theta_loo {
        let d = t - theta_bar;
        s += d * d;
    }
    let var = (g_f - 1.0) / g_f * s;
    var.max(0.0).sqrt()
}

struct RpAccumulator {
    l: usize,
    count: usize,
    sum: [[f64; 3]; 3],
}

impl RpAccumulator {
    fn new(l: usize) -> Self {
        Self {
            l,
            count: 0,
            sum: [[0.0; 3]; 3],
        }
    }

    fn observe(&mut self, field: &Su2Gauge3D) {
        let f = f_vec(field, self.l);
        let theta = theta_f_vec(field, self.l);
        for i in 0..3 {
            for j in 0..3 {
                self.sum[i][j] += theta[i] * f[j];
            }
        }
        self.count += 1;
    }
}

fn build_reflection_positivity_report(
    l: usize,
    acc: Option<RpAccumulator>,
) -> ReflectionPositivityReport {
    let observables = vec!["F1".to_string(), "F2".to_string(), "F3".to_string()];
    let mut matrix = vec![vec![f64::NAN; 3]; 3];
    if let Some(a) = acc {
        let denom = a.count as f64;
        let mut m = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                m[i][j] = a.sum[i][j] / denom;
            }
        }
        for i in 0..3 {
            for j in 0..3 {
                let sym = 0.5 * (m[i][j] + m[j][i]);
                matrix[i][j] = sym;
            }
        }
    }

    let mut a3 = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            a3[i][j] = matrix[i][j];
        }
    }
    let eigen = jacobi_eigenvalues_3x3(a3);
    let eigenvalues = vec![eigen[0], eigen[1], eigen[2]];
    let min_eigenvalue = eigenvalues.iter().copied().fold(f64::INFINITY, f64::min);
    let max_eigenvalue = eigenvalues
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let eps = 1e-8;
    let classification = if min_eigenvalue >= -eps {
        "positive"
    } else if max_eigenvalue <= eps {
        "negative"
    } else {
        "indefinite"
    }
    .to_string();

    ReflectionPositivityReport {
        l,
        observables,
        matrix,
        eigenvalues,
        min_eigenvalue,
        classification,
    }
}

fn region_plaquette_mean(field: &Su2Gauge3D, x0: usize, x1: usize, y0: usize, y1: usize) -> f64 {
    let mut s = 0.0;
    let mut n = 0usize;
    for z in 0..field.l {
        for y in y0..y1 {
            for x in x0..x1 {
                s += field.plaquette_cos(x, y, z);
                n += 1;
            }
        }
    }
    s / (n as f64)
}

fn region_plaquette_mean_reflect_y(
    field: &Su2Gauge3D,
    l: usize,
    x0: usize,
    x1: usize,
    y0: usize,
    y1: usize,
) -> f64 {
    let mut s = 0.0;
    let mut n = 0usize;
    for z in 0..field.l {
        for y in y0..y1 {
            let yr = (l - 1 - y) % l;
            for x in x0..x1 {
                s += field.plaquette_cos(x, yr, z);
                n += 1;
            }
        }
    }
    s / (n as f64)
}

fn region_neighbor_corr_x(field: &Su2Gauge3D, y0: usize, y1: usize) -> f64 {
    let l = field.l;
    let mut s = 0.0;
    let mut n = 0usize;
    for z in 0..l {
        for y in y0..y1 {
            for x in 0..l {
                let xp = (x + 1) % l;
                s += field.plaquette_cos(x, y, z) * field.plaquette_cos(xp, y, z);
                n += 1;
            }
        }
    }
    s / (n as f64)
}

fn region_neighbor_corr_x_reflect_y(field: &Su2Gauge3D, l: usize, y0: usize, y1: usize) -> f64 {
    let ll = field.l;
    let mut s = 0.0;
    let mut n = 0usize;
    for z in 0..ll {
        for y in y0..y1 {
            let yr = (l - 1 - y) % l;
            for x in 0..ll {
                let xp = (x + 1) % ll;
                s += field.plaquette_cos(x, yr, z) * field.plaquette_cos(xp, yr, z);
                n += 1;
            }
        }
    }
    s / (n as f64)
}

fn f_vec(field: &Su2Gauge3D, l: usize) -> [f64; 3] {
    let half = l / 2;
    let f1 = region_plaquette_mean(field, 0, half, 0, half);
    let f2 = region_plaquette_mean(field, half, l, 0, half);
    let f3 = region_neighbor_corr_x(field, 0, half);
    [f1, f2, f3]
}

fn theta_f_vec(field: &Su2Gauge3D, l: usize) -> [f64; 3] {
    let half = l / 2;
    let f1 = region_plaquette_mean_reflect_y(field, l, 0, half, 0, half);
    let f2 = region_plaquette_mean_reflect_y(field, l, half, l, 0, half);
    let f3 = region_neighbor_corr_x_reflect_y(field, l, 0, half);
    [f1, f2, f3]
}

fn jacobi_eigenvalues_3x3(mut a: [[f64; 3]; 3]) -> [f64; 3] {
    for _ in 0..64 {
        let mut p = 0usize;
        let mut q = 1usize;
        let mut max = a[0][1].abs();
        let v02 = a[0][2].abs();
        if v02 > max {
            max = v02;
            p = 0;
            q = 2;
        }
        let v12 = a[1][2].abs();
        if v12 > max {
            max = v12;
            p = 1;
            q = 2;
        }
        if !max.is_finite() || max < 1e-10 {
            break;
        }

        let app = a[p][p];
        let aqq = a[q][q];
        let apq = a[p][q];
        let phi = 0.5 * (2.0 * apq).atan2(aqq - app);
        let c = phi.cos();
        let s = phi.sin();

        let app_new = c * c * app - 2.0 * s * c * apq + s * s * aqq;
        let aqq_new = s * s * app + 2.0 * s * c * apq + c * c * aqq;
        a[p][p] = app_new;
        a[q][q] = aqq_new;
        a[p][q] = 0.0;
        a[q][p] = 0.0;

        for r in 0..3 {
            if r != p && r != q {
                let arp = a[r][p];
                let arq = a[r][q];
                a[r][p] = c * arp - s * arq;
                a[p][r] = a[r][p];
                a[r][q] = s * arp + c * arq;
                a[q][r] = a[r][q];
            }
        }
    }
    [a[0][0], a[1][1], a[2][2]]
}

#[cfg(test)]
mod mill_analysis_stats {
    use super::*;

    #[test]
    fn mill_analysis_stats_jackknife_detects_block_level_correlations() {
        let l = 24usize;
        let expected_n_measurements = 80usize;
        let mut acc = MassAccumulator::new(l, expected_n_measurements);

        let m = 0.2f64;
        let rel = 0.05f64;
        let r_max = acc.r_max;

        let mut base_corr: Vec<f64> = Vec::with_capacity(r_max);
        for r in 1..=r_max {
            base_corr.push((-m * (r as f64)).exp());
        }

        let mut a_r: Vec<f64> = Vec::with_capacity(r_max);
        for r in 1..=r_max {
            a_r.push(0.03 * (r as f64) / (r_max as f64));
        }

        let factors: Vec<Vec<f64>> = vec![
            a_r.iter().map(|a| 1.0 + a).collect(),
            a_r.iter().map(|a| 1.0 - a).collect(),
            a_r.iter().map(|a| 1.0 + 2.0 * a).collect(),
            a_r.iter().map(|a| 1.0 - 2.0 * a).collect(),
        ];

        for f in factors.iter() {
            for &x in f.iter() {
                assert!(x > 0.0);
            }
        }

        for block in 0..4 {
            for _ in 0..4 {
                let mut mean_pp = vec![0.0; r_max];
                for i in 0..r_max {
                    mean_pp[i] = base_corr[i] * factors[block][i];
                }
                acc.observe_stats(0.0, &mean_pp);
            }
        }

        let (_corr, m_eff) = acc.correlator_and_m_eff_from_sums(acc.count, acc.sum_p, &acc.sum_pp);
        let plateau_naive = find_mass_plateau(&m_eff, rel);
        let width_naive = plateau_naive.r_end.saturating_sub(plateau_naive.r_start);

        assert!(plateau_naive.m_eff_mean.is_finite());
        assert!(plateau_naive.m_eff_std.is_finite());
        assert!(plateau_naive.m_eff_std < 1e-12);

        let p = mass_scaling_plateau_from_acc_with_cfg(
            &acc,
            PlateauCfg {
                mode: PlateauMode::Legacy,
                rel_thresh: rel,
                k: 2.0,
                chi2_max: 2.0,
            },
        )
        .unwrap();
        let width = plateau_width_from_r_range(p.r_start, p.r_end);
        assert_eq!(width, width_naive);
        assert!(p.m_eff_mean.is_finite());
        assert!(p.m_eff_std.is_finite());
        assert!(p.m_eff_std > 1e-8);
        assert!((p.m_eff_mean - plateau_naive.m_eff_mean).abs() < 1e-3);
    }

    #[test]
    fn mill_plateau_stat_recovers_when_rel_thresh_fails() {
        let m_eff = vec![1.0, 1.2, 0.9, 1.1];
        let sigma = vec![0.3, 0.3, 0.3, 0.3];
        let stat = find_mass_plateau_stat(&m_eff, &sigma, 2.0, 2.0).unwrap();
        assert_eq!(stat.start, 0);
        assert_eq!(stat.end, 3);
        assert_eq!(stat.n_points, 4);
        assert!(stat.mu.is_finite());
        assert!(stat.chi2_dof.is_finite());
        assert!(stat.chi2_dof <= 2.0);
        assert_eq!(stat.weights.len(), 4);

        let legacy = find_mass_plateau(&m_eff, 0.05);
        let legacy_width = legacy.r_end.saturating_sub(legacy.r_start);
        assert_eq!(legacy_width, 1);
    }

    #[test]
    fn mill_plateau_stat_rejects_when_chi2_is_bad_even_if_local_passes() {
        let m_eff = vec![0.0, 1e6, -1e6, 2e6];
        let sigma = vec![1000.0, 1000.0, 1000.0, 1000.0];
        let stat = find_mass_plateau_stat(&m_eff, &sigma, 1e6, 2.0);
        assert!(stat.is_none());
    }
}

#[cfg(test)]
mod su2_kernel_tests {
    use super::*;

    fn approx_eq(a: f64, b: f64, eps: f64) {
        assert!(
            (a - b).abs() <= eps,
            "expected |{} - {}| <= {}",
            a,
            b,
            eps
        );
    }

    #[test]
    fn su2_plaquettes_are_orientation_consistent() {
        let l = 6usize;
        let mut rng = StdRng::seed_from_u64(123);
        let mut field = Su2Gauge3D::new(l);
        for z in 0..l {
            for y in 0..l {
                for x in 0..l {
                    field.set_link_x(x, y, z, su2_random_haar(&mut rng));
                    field.set_link_y(x, y, z, su2_random_haar(&mut rng));
                    field.set_link_z(x, y, z, su2_random_haar(&mut rng));
                }
            }
        }

        for z in 0..l {
            for y in 0..l {
                for x in 0..l {
                    let xp = field.wrap(x as isize + 1);
                    let yp = field.wrap(y as isize + 1);
                    let zp = field.wrap(z as isize + 1);

                    // Reverse orientation should give dagger; plaquette_value (a0) must match.
                    let p_xy = field.plaquette_xy(x, y, z).plaquette_value();
                    let p_yx = field
                        .link_y(x, y, z)
                        .mul(field.link_x(x, yp, z))
                        .mul(field.link_y(xp, y, z).dagger())
                        .mul(field.link_x(x, y, z).dagger())
                        .plaquette_value();
                    approx_eq(p_xy, p_yx, 1e-10);

                    let p_xz = field.plaquette_xz(x, y, z).plaquette_value();
                    let p_zx = field
                        .link_z(x, y, z)
                        .mul(field.link_x(x, y, zp))
                        .mul(field.link_z(xp, y, z).dagger())
                        .mul(field.link_x(x, y, z).dagger())
                        .plaquette_value();
                    approx_eq(p_xz, p_zx, 1e-10);

                    let p_yz = field.plaquette_yz(x, y, z).plaquette_value();
                    let p_zy = field
                        .link_z(x, y, z)
                        .mul(field.link_y(x, y, zp))
                        .mul(field.link_z(x, yp, z).dagger())
                        .mul(field.link_y(x, y, z).dagger())
                        .plaquette_value();
                    approx_eq(p_yz, p_zy, 1e-10);
                }
            }
        }
    }

    #[test]
    fn su2_ape_smearing_projects_back_to_su2() {
        let l = 6usize;
        let mut rng = StdRng::seed_from_u64(777);
        let mut field = Su2Gauge3D::new(l);
        for z in 0..l {
            for y in 0..l {
                for x in 0..l {
                    field.set_link_x(x, y, z, su2_random_haar(&mut rng));
                    field.set_link_y(x, y, z, su2_random_haar(&mut rng));
                    field.set_link_z(x, y, z, su2_random_haar(&mut rng));
                }
            }
        }

        let smeared = ape_smear_su2(&field, 0.5, 2);
        let mut max_dev = 0.0f64;
        for z in 0..l {
            for y in 0..l {
                for x in 0..l {
                    for u in [smeared.link_x(x, y, z), smeared.link_y(x, y, z), smeared.link_z(x, y, z)] {
                        max_dev = max_dev.max((u.norm2() - 1.0).abs());
                    }
                }
            }
        }
        assert!(max_dev < 1e-10, "max |norm2-1| too large: {}", max_dev);
    }

    #[test]
    fn su2_heatbath_sampler_stays_in_su2() {
        let mut rng = StdRng::seed_from_u64(42);
        let u = su2_heatbath_sample(3.0, &mut rng);
        approx_eq(u.norm2(), 1.0, 1e-12);
    }

    #[test]
    fn su2_overrelax_preserves_trace_against_unit_staple() {
        let mut rng = StdRng::seed_from_u64(1234);
        for _ in 0..200 {
            let u = su2_random_haar(&mut rng);
            let v = su2_random_haar(&mut rng);
            let t0 = u.mul(v).plaquette_value();
            let u2 = v.dagger().mul(u.dagger()).mul(v.dagger()).projected();
            let t1 = u2.mul(v).plaquette_value();
            approx_eq(t0, t1, 1e-12);
        }
    }

    #[test]
    fn su2_staple_sum_matches_local_cos_sum() {
        let l = 6usize;
        let mut rng = StdRng::seed_from_u64(999);
        let mut field = Su2Gauge3D::new(l);
        for z in 0..l {
            for y in 0..l {
                for x in 0..l {
                    field.set_link_x(x, y, z, su2_random_haar(&mut rng));
                    field.set_link_y(x, y, z, su2_random_haar(&mut rng));
                    field.set_link_z(x, y, z, su2_random_haar(&mut rng));
                }
            }
        }

        for _ in 0..200 {
            let x = rng.gen_range(0..l);
            let y = rng.gen_range(0..l);
            let z = rng.gen_range(0..l);

            let u = field.link_x(x, y, z);
            let local = local_cos_sum_for_link_x(&field, x, y, z);
            let staple = staple_sum_for_link_x(&field, x, y, z);
            let traced = u.mul(staple).plaquette_value();
            approx_eq(local, traced, 1e-10);

            let u = field.link_y(x, y, z);
            let local = local_cos_sum_for_link_y(&field, x, y, z);
            let staple = staple_sum_for_link_y(&field, x, y, z);
            let traced = u.mul(staple).plaquette_value();
            approx_eq(local, traced, 1e-10);

            let u = field.link_z(x, y, z);
            let local = local_cos_sum_for_link_z(&field, x, y, z);
            let staple = staple_sum_for_link_z(&field, x, y, z);
            let traced = u.mul(staple).plaquette_value();
            approx_eq(local, traced, 1e-10);
        }
    }
}
