//! Preconditioners for LSMR solver.
//!
//! This module defines the `RightPreconditioner` trait and provides
//! implementations for various preconditioning strategies.
//!
//! # Right Preconditioning
//!
//! For right preconditioning, we transform the problem:
//! ```text
//! min ||A x - b||  →  min ||A M^{-1} z - b||,  x = M^{-1} z
//! ```
//!
//! This requires applying M^{-1} (forward) and M^{-T} (transpose).
//! For symmetric preconditioners, M^{-T} = M^{-1}.

use crate::demean::types::DemeanContext;

/// Right preconditioner for LSMR.
///
/// Implementations must be:
/// - **Linear**: M^{-1}(ax + by) = a M^{-1}(x) + b M^{-1}(y)
/// - **Deterministic**: Same input always produces same output
/// - **Allocation-free** in `apply`/`apply_transpose`
///
/// For preconditioners with internal iteration (like P3 streaming),
/// the number of inner iterations must be fixed.
pub trait RightPreconditioner {
    /// Apply the preconditioner: z = M^{-1} * x.
    ///
    /// # Arguments
    /// * `x` - Input vector (length: n_coef)
    /// * `z` - Output vector (length: n_coef), overwritten
    fn apply(&self, x: &[f64], z: &mut [f64]);

    /// Apply the transpose: z = M^{-T} * x.
    ///
    /// Default implementation assumes symmetric preconditioner (M^{-T} = M^{-1}).
    /// Override for non-symmetric preconditioners.
    #[allow(dead_code)]
    fn apply_transpose(&self, x: &[f64], z: &mut [f64]) {
        self.apply(x, z);
    }
}

// =============================================================================
// P0: Identity Preconditioner (No preconditioning)
// =============================================================================

/// Identity preconditioner (no preconditioning).
///
/// Used as baseline for benchmarking and when no preconditioning is desired.
/// Simply copies input to output.
#[derive(Clone, Copy, Default)]
pub struct IdentityPreconditioner;

impl RightPreconditioner for IdentityPreconditioner {
    #[inline]
    fn apply(&self, x: &[f64], z: &mut [f64]) {
        z.copy_from_slice(x);
    }
}

// =============================================================================
// P1: Diagonal Preconditioner (Count Scaling)
// =============================================================================

/// Diagonal preconditioner using inverse square root of group counts.
///
/// For each coefficient, applies scaling by `sqrt(1 / group_count)`.
/// This normalizes coefficients by their statistical weight, which
/// improves conditioning when group sizes vary significantly.
///
/// # Cost
/// - Construction: O(n_coef) to compute sqrt(inv_group_weights)
/// - Apply: O(n_coef)
pub struct DiagonalPreconditioner {
    /// Precomputed sqrt(inv_group_weights) for each coefficient
    inv_sqrt_weights: Vec<f64>,
}

impl DiagonalPreconditioner {
    /// Create a diagonal preconditioner from a `DemeanContext`.
    ///
    /// Uses the precomputed `inv_group_weights` from each FE, taking
    /// the square root for symmetric preconditioning.
    pub fn new(ctx: &DemeanContext) -> Self {
        let mut inv_sqrt_weights = Vec::with_capacity(ctx.dims.n_coef);

        for fe in &ctx.fe_infos {
            for &w in &fe.inv_group_weights {
                inv_sqrt_weights.push(w.sqrt());
            }
        }

        debug_assert_eq!(inv_sqrt_weights.len(), ctx.dims.n_coef);

        Self { inv_sqrt_weights }
    }
}

impl RightPreconditioner for DiagonalPreconditioner {
    #[inline]
    fn apply(&self, x: &[f64], z: &mut [f64]) {
        debug_assert_eq!(x.len(), self.inv_sqrt_weights.len());
        debug_assert_eq!(z.len(), self.inv_sqrt_weights.len());

        for ((z_i, &x_i), &w) in z.iter_mut().zip(x.iter()).zip(self.inv_sqrt_weights.iter()) {
            *z_i = x_i * w;
        }
    }
}

// =============================================================================
// P2: Nullspace Deflation Preconditioner
// =============================================================================

/// Nullspace deflation preconditioner.
///
/// Projects out the constant (null space) component from each FE block.
/// For each FE, the coefficients are centered by subtracting their weighted mean.
///
/// This handles rank-deficiency cleanly: FE coefficients are only identified
/// up to a constant within each FE, and deflation removes this ambiguity.
///
/// # Cost
/// - Construction: O(n_fe) to store block boundaries
/// - Apply: O(n_coef) to compute and subtract means
pub struct NullspaceDeflator {
    /// Start and end indices for each FE block, plus their total weights
    /// Each tuple is (start, end, sum_of_weights)
    fe_blocks: Vec<(usize, usize, f64)>,
    /// Weights for computing weighted mean (inv_group_weights from context)
    weights: Vec<f64>,
}

impl NullspaceDeflator {
    /// Create a nullspace deflator from a `DemeanContext`.
    pub fn new(ctx: &DemeanContext) -> Self {
        let mut fe_blocks = Vec::with_capacity(ctx.dims.n_fe);
        let mut weights = Vec::with_capacity(ctx.dims.n_coef);

        for fe in &ctx.fe_infos {
            let start = fe.coef_start;
            let end = start + fe.n_groups;
            // Sum of weights (inverse of inv_group_weights = group counts)
            let sum_weights: f64 = fe.inv_group_weights.iter().map(|&w| 1.0 / w).sum();
            fe_blocks.push((start, end, sum_weights));

            // Store 1/inv_group_weights = group_count for weighted mean
            for &w in &fe.inv_group_weights {
                weights.push(1.0 / w); // group_count
            }
        }

        Self { fe_blocks, weights }
    }
}

impl RightPreconditioner for NullspaceDeflator {
    fn apply(&self, x: &[f64], z: &mut [f64]) {
        // First copy x to z
        z.copy_from_slice(x);

        // For each FE block, subtract the weighted mean
        // Formula: z_i = x_i - (sum_j x_j * w_j) / W
        for &(start, end, sum_weights) in &self.fe_blocks {
            // Compute weighted sum
            let weighted_sum: f64 = x[start..end]
                .iter()
                .zip(self.weights[start..end].iter())
                .map(|(&xi, &wi)| xi * wi)
                .sum();

            // Weighted mean
            let mean = weighted_sum / sum_weights;

            // Subtract mean from each coefficient in this block
            for z_i in z[start..end].iter_mut() {
                *z_i -= mean;
            }
        }
    }

    fn apply_transpose(&self, x: &[f64], z: &mut [f64]) {
        // The deflator matrix is: M = I - (1/W) * 1 * w^T
        // Its transpose is:       M^T = I - (1/W) * w * 1^T
        // So: apply_transpose(x)_i = x_i - (w_i / W) * sum(x)

        // First copy x to z
        z.copy_from_slice(x);

        // For each FE block, apply the transpose formula
        for &(start, end, sum_weights) in &self.fe_blocks {
            // Compute unweighted sum
            let total_sum: f64 = x[start..end].iter().sum();

            // Subtract (w_i / W) * sum for each coefficient
            for (z_i, &w_i) in z[start..end].iter_mut().zip(self.weights[start..end].iter()) {
                *z_i -= (w_i / sum_weights) * total_sum;
            }
        }
    }
}

// =============================================================================
// Composed Preconditioner
// =============================================================================

/// Composed preconditioner: applies P1 then P2.
///
/// For right preconditioning: M^{-1} = P2^{-1} * P1^{-1}
/// So apply(x) computes: P2(P1(x))
pub struct ComposedPreconditioner<P1, P2> {
    p1: P1,
    p2: P2,
    /// Scratch buffer for intermediate result (uses Mutex for thread safety)
    scratch: std::sync::Mutex<Vec<f64>>,
}

impl<P1: RightPreconditioner, P2: RightPreconditioner> ComposedPreconditioner<P1, P2> {
    /// Create a composed preconditioner.
    pub fn new(p1: P1, p2: P2, n_coef: usize) -> Self {
        Self {
            p1,
            p2,
            scratch: std::sync::Mutex::new(vec![0.0; n_coef]),
        }
    }
}

impl<P1: RightPreconditioner, P2: RightPreconditioner> RightPreconditioner
    for ComposedPreconditioner<P1, P2>
{
    fn apply(&self, x: &[f64], z: &mut [f64]) {
        // z = P2(P1(x))
        let mut scratch = self.scratch.lock().unwrap();
        self.p1.apply(x, &mut scratch);
        self.p2.apply(&scratch, z);
    }

    fn apply_transpose(&self, x: &[f64], z: &mut [f64]) {
        // For transpose: (P2 * P1)^T = P1^T * P2^T
        let mut scratch = self.scratch.lock().unwrap();
        self.p2.apply_transpose(x, &mut scratch);
        self.p1.apply_transpose(&scratch, z);
    }
}

// =============================================================================
// P3: Two-Block Streaming Preconditioner
// =============================================================================

/// Two-block streaming Gauss-Seidel preconditioner.
///
/// Performs fixed inner Gauss-Seidel iterations on a 2-FE subsystem.
/// The two FEs can be auto-selected based on cross-correlation (density of
/// the D_p^T * D_q matrix) or manually specified.
///
/// This preconditioner is most effective for 3+ FE problems where two FEs
/// have strong coupling.
///
/// # Cost
/// - Construction: O(1)
/// - Apply: O(n_obs × inner_iters) per application
pub struct TwoBlockStreamingPreconditioner {
    /// Index of first FE in the 2-block subsystem
    #[allow(dead_code)]
    fe_p: usize,
    /// Index of second FE in the 2-block subsystem
    #[allow(dead_code)]
    fe_q: usize,
    /// Number of inner GS iterations
    inner_iters: usize,
    /// Coefficient start offset for FE p
    coef_start_p: usize,
    /// Coefficient start offset for FE q
    coef_start_q: usize,
    /// Number of groups in FE p
    n_groups_p: usize,
    /// Number of groups in FE q
    n_groups_q: usize,
    /// Group IDs for FE p (for each observation)
    group_ids_p: Vec<usize>,
    /// Group IDs for FE q (for each observation)
    group_ids_q: Vec<usize>,
    /// Inverse group weights for FE p (1/count)
    inv_weights_p: Vec<f64>,
    /// Inverse group weights for FE q (1/count)
    inv_weights_q: Vec<f64>,
    /// Preallocated scratch buffers (use Mutex for interior mutability)
    scratch_p: std::sync::Mutex<Vec<f64>>,
    scratch_q: std::sync::Mutex<Vec<f64>>,
}

impl TwoBlockStreamingPreconditioner {
    /// Create a two-block streaming preconditioner.
    ///
    /// # Arguments
    /// * `ctx` - DemeanContext with FE information
    /// * `fe_p` - First FE index (or None for auto-select)
    /// * `fe_q` - Second FE index (or None for auto-select)
    /// * `inner_iters` - Number of inner GS iterations (typically 2-5)
    pub fn new(
        ctx: &DemeanContext,
        fe_p: Option<usize>,
        fe_q: Option<usize>,
        inner_iters: usize,
    ) -> Self {
        // Auto-select FE pair if not specified
        let (p, q) = match (fe_p, fe_q) {
            (Some(p), Some(q)) => (p, q),
            _ => select_best_fe_pair(ctx),
        };

        let fe_p_info = &ctx.fe_infos[p];
        let fe_q_info = &ctx.fe_infos[q];

        let n_groups_p = fe_p_info.n_groups;
        let n_groups_q = fe_q_info.n_groups;

        Self {
            fe_p: p,
            fe_q: q,
            inner_iters,
            coef_start_p: fe_p_info.coef_start,
            coef_start_q: fe_q_info.coef_start,
            n_groups_p,
            n_groups_q,
            group_ids_p: fe_p_info.group_ids.clone(),
            group_ids_q: fe_q_info.group_ids.clone(),
            inv_weights_p: fe_p_info.inv_group_weights.clone(),
            inv_weights_q: fe_q_info.inv_group_weights.clone(),
            scratch_p: std::sync::Mutex::new(vec![0.0; n_groups_p]),
            scratch_q: std::sync::Mutex::new(vec![0.0; n_groups_q]),
        }
    }

    /// Perform one GS sweep: update FE p coefficients given FE q coefficients.
    ///
    /// Solves: z_p[g] = (x_p[g] - sum_{i in g} z_q[group_q[i]]) / count_g
    fn sweep_p(&self, x: &[f64], z: &mut [f64], rhs: &mut [f64]) {
        // Initialize RHS from input
        rhs.iter_mut()
            .enumerate()
            .for_each(|(g, r)| *r = x[self.coef_start_p + g]);

        // Streaming matvec: rhs[gp] -= sum over obs with group_p=gp of z_q[group_q]
        for (&gp, &gq) in self.group_ids_p.iter().zip(self.group_ids_q.iter()) {
            rhs[gp] -= z[self.coef_start_q + gq];
        }

        // Divide by diagonal
        for (g, (&r, &inv_w)) in rhs.iter().zip(self.inv_weights_p.iter()).enumerate() {
            z[self.coef_start_p + g] = r * inv_w;
        }
    }

    /// Perform one GS sweep: update FE q coefficients given FE p coefficients.
    ///
    /// Solves: z_q[h] = (x_q[h] - sum_{i in h} z_p[group_p[i]]) / count_h
    fn sweep_q(&self, x: &[f64], z: &mut [f64], rhs: &mut [f64]) {
        // Initialize RHS from input
        rhs.iter_mut()
            .enumerate()
            .for_each(|(h, r)| *r = x[self.coef_start_q + h]);

        // Streaming matvec: rhs[gq] -= sum over obs with group_q=gq of z_p[group_p]
        for (&gp, &gq) in self.group_ids_p.iter().zip(self.group_ids_q.iter()) {
            rhs[gq] -= z[self.coef_start_p + gp];
        }

        // Divide by diagonal
        for (h, (&r, &inv_w)) in rhs.iter().zip(self.inv_weights_q.iter()).enumerate() {
            z[self.coef_start_q + h] = r * inv_w;
        }
    }
}

impl RightPreconditioner for TwoBlockStreamingPreconditioner {
    fn apply(&self, x: &[f64], z: &mut [f64]) {
        // Initialize z = 0 for GS iteration (solving M*z = x)
        z.iter_mut().for_each(|v| *v = 0.0);

        // Get scratch buffers
        let mut scratch_p = self.scratch_p.lock().unwrap();
        let mut scratch_q = self.scratch_q.lock().unwrap();

        // Perform inner_iters Gauss-Seidel iterations on the 2-FE subsystem
        // Each sweep uses the original input x as the RHS
        // Forward order: p then q
        for _ in 0..self.inner_iters {
            self.sweep_p(x, z, &mut scratch_p);
            self.sweep_q(x, z, &mut scratch_q);
        }

        // Copy non-preconditioned coefficients (other FEs) from input
        // The sweeps only modify FEs p and q; other FE coefficients pass through
        for (i, (&xi, zi)) in x.iter().zip(z.iter_mut()).enumerate() {
            let in_p = i >= self.coef_start_p && i < self.coef_start_p + self.n_groups_p;
            let in_q = i >= self.coef_start_q && i < self.coef_start_q + self.n_groups_q;
            if !in_p && !in_q {
                *zi = xi;
            }
        }
    }

    fn apply_transpose(&self, x: &[f64], z: &mut [f64]) {
        // The transpose of GS iteration requires reversing the sweep order.
        // If forward is: sweep_p then sweep_q (inner loop)
        // Then transpose is: sweep_q then sweep_p
        //
        // This is because for Gauss-Seidel solving M z = x where M = L + D + U:
        // Forward GS uses (L + D)^{-1}
        // Backward GS (transpose) uses (D + U)^{-1} = (L^T + D)^{-1}
        // which corresponds to reversing the update order.

        z.iter_mut().for_each(|v| *v = 0.0);

        // Get scratch buffers
        let mut scratch_p = self.scratch_p.lock().unwrap();
        let mut scratch_q = self.scratch_q.lock().unwrap();

        // Reverse order: q then p
        for _ in 0..self.inner_iters {
            self.sweep_q(x, z, &mut scratch_q);
            self.sweep_p(x, z, &mut scratch_p);
        }

        // Copy non-preconditioned coefficients from input
        for (i, (&xi, zi)) in x.iter().zip(z.iter_mut()).enumerate() {
            let in_p = i >= self.coef_start_p && i < self.coef_start_p + self.n_groups_p;
            let in_q = i >= self.coef_start_q && i < self.coef_start_q + self.n_groups_q;
            if !in_p && !in_q {
                *zi = xi;
            }
        }
    }
}

/// Select the best FE pair for two-block preconditioning.
///
/// Uses cross-correlation (density of D_p^T * D_q) to identify the most
/// coupled FE pair. Higher density means more observations share both FE
/// groups, indicating stronger coupling.
fn select_best_fe_pair(ctx: &DemeanContext) -> (usize, usize) {
    let n_fe = ctx.dims.n_fe;

    // For 2 FEs, the choice is obvious
    if n_fe == 2 {
        return (0, 1);
    }

    // For 1 FE, return (0, 0) - caller should handle this edge case
    if n_fe < 2 {
        return (0, 0);
    }

    let mut best_pair = (0, 1);
    let mut best_density = 0.0;

    for p in 0..n_fe {
        for q in (p + 1)..n_fe {
            let density = compute_pair_density(
                &ctx.fe_infos[p].group_ids,
                &ctx.fe_infos[q].group_ids,
                ctx.fe_infos[p].n_groups,
                ctx.fe_infos[q].n_groups,
            );
            if density > best_density {
                best_density = density;
                best_pair = (p, q);
            }
        }
    }

    best_pair
}

/// Compute the density of the cross-correlation matrix D_p^T * D_q.
///
/// Density = (number of distinct (group_p, group_q) pairs) / (n_groups_p * n_groups_q)
/// Higher density indicates stronger coupling between the two FEs.
fn compute_pair_density(
    group_ids_p: &[usize],
    group_ids_q: &[usize],
    n_groups_p: usize,
    n_groups_q: usize,
) -> f64 {
    use std::collections::HashSet;

    let mut pairs: HashSet<(usize, usize)> = HashSet::new();
    for (&gp, &gq) in group_ids_p.iter().zip(group_ids_q.iter()) {
        pairs.insert((gp, gq));
    }

    pairs.len() as f64 / (n_groups_p * n_groups_q) as f64
}

// =============================================================================
// P4: Sparse Gram Preconditioner
// =============================================================================

/// Sparse Gram preconditioner using Gauss-Seidel on sparse cross-correlation.
///
/// Similar to TwoBlockStreamingPreconditioner, but uses sparse matrix-vector
/// products instead of streaming over observations. This makes each GS sweep
/// O(nnz) where nnz is the number of non-zero entries in M_pq = D_p^T D_q,
/// rather than O(n_obs) for streaming.
///
/// For problems where nnz << n_obs (many groups have no shared observations),
/// this can be significantly faster.
///
/// # Cost
/// - Construction: O(n_obs) to build sparse cross-correlation matrix
/// - Apply: O((n_groups_p + n_groups_q + nnz) × inner_iters)
pub struct SparseGramPreconditioner {
    /// Coefficient start offset for FE p
    coef_start_p: usize,
    /// Coefficient start offset for FE q
    coef_start_q: usize,
    /// Number of groups in FE p
    n_groups_p: usize,
    /// Number of groups in FE q
    n_groups_q: usize,
    /// Inverse group weights for FE p (1/count)
    inv_diag_p: Vec<f64>,
    /// Inverse group weights for FE q (1/count)
    inv_diag_q: Vec<f64>,
    /// Sparse CSR matrix M_pq = D_p^T D_q
    mpq_row_ptr: Vec<usize>,
    mpq_col_idx: Vec<usize>,
    mpq_values: Vec<f64>,
    /// Sparse CSR matrix M_qp = M_pq^T
    mqp_row_ptr: Vec<usize>,
    mqp_col_idx: Vec<usize>,
    mqp_values: Vec<f64>,
    /// Number of inner GS iterations
    inner_iters: usize,
    /// Scratch buffers (protected by Mutex for interior mutability)
    scratch: std::sync::Mutex<SparseGSScratch>,
}

/// Scratch buffers for sparse GS iterations.
struct SparseGSScratch {
    z_p: Vec<f64>,
    z_q: Vec<f64>,
    temp_p: Vec<f64>,  // temp storage for M_pq @ z_q
    temp_q: Vec<f64>,  // temp storage for M_qp @ z_p
}

impl SparseGramPreconditioner {
    /// Create a sparse Gram preconditioner.
    pub fn new(
        ctx: &DemeanContext,
        fe_p: Option<usize>,
        fe_q: Option<usize>,
        inner_iters: usize,
    ) -> Self {
        // Auto-select FE pair if not specified
        let (p, q) = match (fe_p, fe_q) {
            (Some(p), Some(q)) => (p, q),
            _ => select_best_fe_pair(ctx),
        };

        let fe_p_info = &ctx.fe_infos[p];
        let fe_q_info = &ctx.fe_infos[q];

        let n_groups_p = fe_p_info.n_groups;
        let n_groups_q = fe_q_info.n_groups;

        // Store inverse diagonal (1/count) for GS sweeps
        let inv_diag_p: Vec<f64> = fe_p_info.inv_group_weights.to_vec();
        let inv_diag_q: Vec<f64> = fe_q_info.inv_group_weights.to_vec();

        // Build sparse CSR representation of M_pq = D_p^T D_q
        use std::collections::HashMap;
        let mut pq_counts: HashMap<(usize, usize), f64> = HashMap::new();

        for (&gp, &gq) in fe_p_info.group_ids.iter().zip(fe_q_info.group_ids.iter()) {
            *pq_counts.entry((gp, gq)).or_insert(0.0) += 1.0;
        }

        // Build CSR for M_pq
        let mut mpq_row_ptr = vec![0usize; n_groups_p + 1];
        let mut mpq_entries: Vec<(usize, usize, f64)> = pq_counts
            .iter()
            .map(|(&(gp, gq), &count)| (gp, gq, count))
            .collect();
        mpq_entries.sort_by_key(|&(gp, gq, _)| (gp, gq));

        for &(gp, _, _) in &mpq_entries {
            mpq_row_ptr[gp + 1] += 1;
        }
        for i in 1..=n_groups_p {
            mpq_row_ptr[i] += mpq_row_ptr[i - 1];
        }

        let mpq_col_idx: Vec<usize> = mpq_entries.iter().map(|&(_, gq, _)| gq).collect();
        let mpq_values: Vec<f64> = mpq_entries.iter().map(|&(_, _, v)| v).collect();

        // Build CSR for M_qp (transpose)
        let mut mqp_row_ptr = vec![0usize; n_groups_q + 1];
        let mut mqp_entries: Vec<(usize, usize, f64)> = pq_counts
            .iter()
            .map(|(&(gp, gq), &count)| (gq, gp, count))
            .collect();
        mqp_entries.sort_by_key(|&(gq, gp, _)| (gq, gp));

        for &(gq, _, _) in &mqp_entries {
            mqp_row_ptr[gq + 1] += 1;
        }
        for i in 1..=n_groups_q {
            mqp_row_ptr[i] += mqp_row_ptr[i - 1];
        }

        let mqp_col_idx: Vec<usize> = mqp_entries.iter().map(|&(_, gp, _)| gp).collect();
        let mqp_values: Vec<f64> = mqp_entries.iter().map(|&(_, _, v)| v).collect();

        // Initialize scratch buffers
        let scratch = SparseGSScratch {
            z_p: vec![0.0; n_groups_p],
            z_q: vec![0.0; n_groups_q],
            temp_p: vec![0.0; n_groups_p],
            temp_q: vec![0.0; n_groups_q],
        };

        Self {
            coef_start_p: fe_p_info.coef_start,
            coef_start_q: fe_q_info.coef_start,
            n_groups_p,
            n_groups_q,
            inv_diag_p,
            inv_diag_q,
            mpq_row_ptr,
            mpq_col_idx,
            mpq_values,
            mqp_row_ptr,
            mqp_col_idx,
            mqp_values,
            inner_iters,
            scratch: std::sync::Mutex::new(scratch),
        }
    }

    /// Sparse matrix-vector product: y = M_pq @ x_q
    #[inline]
    fn mpq_matvec(&self, x_q: &[f64], y_p: &mut [f64]) {
        for (i, y_i) in y_p.iter_mut().enumerate() {
            let start = self.mpq_row_ptr[i];
            let end = self.mpq_row_ptr[i + 1];
            let mut sum = 0.0;
            for idx in start..end {
                let j = self.mpq_col_idx[idx];
                sum += self.mpq_values[idx] * x_q[j];
            }
            *y_i = sum;
        }
    }

    /// Sparse matrix-vector product: y = M_qp @ x_p
    #[inline]
    fn mqp_matvec(&self, x_p: &[f64], y_q: &mut [f64]) {
        for (i, y_i) in y_q.iter_mut().enumerate() {
            let start = self.mqp_row_ptr[i];
            let end = self.mqp_row_ptr[i + 1];
            let mut sum = 0.0;
            for idx in start..end {
                let j = self.mqp_col_idx[idx];
                sum += self.mqp_values[idx] * x_p[j];
            }
            *y_i = sum;
        }
    }
}

impl RightPreconditioner for SparseGramPreconditioner {
    fn apply(&self, x: &[f64], z: &mut [f64]) {
        let mut scratch = self.scratch.lock().unwrap();
        let SparseGSScratch {
            z_p,
            z_q,
            temp_p,
            temp_q,
        } = &mut *scratch;

        // Extract x_p and x_q from the full coefficient vector
        let x_p = &x[self.coef_start_p..self.coef_start_p + self.n_groups_p];
        let x_q = &x[self.coef_start_q..self.coef_start_q + self.n_groups_q];

        // Initialize z_p, z_q = 0
        z_p.iter_mut().for_each(|v| *v = 0.0);
        z_q.iter_mut().for_each(|v| *v = 0.0);

        // Perform inner_iters Gauss-Seidel iterations using sparse matvec
        // Solving: diag(count_p) @ z_p + M_pq @ z_q = x_p
        //          M_qp @ z_p + diag(count_q) @ z_q = x_q
        for _ in 0..self.inner_iters {
            // sweep_p: z_p = (x_p - M_pq @ z_q) * inv_diag_p
            self.mpq_matvec(z_q, temp_p);
            for ((z_i, &t), (&x_i, &inv_d)) in z_p
                .iter_mut()
                .zip(temp_p.iter())
                .zip(x_p.iter().zip(self.inv_diag_p.iter()))
            {
                *z_i = (x_i - t) * inv_d;
            }

            // sweep_q: z_q = (x_q - M_qp @ z_p) * inv_diag_q
            self.mqp_matvec(z_p, temp_q);
            for ((z_i, &t), (&x_i, &inv_d)) in z_q
                .iter_mut()
                .zip(temp_q.iter())
                .zip(x_q.iter().zip(self.inv_diag_q.iter()))
            {
                *z_i = (x_i - t) * inv_d;
            }
        }

        // Copy result back to z, pass through other coefficients
        z.copy_from_slice(x);
        z[self.coef_start_p..self.coef_start_p + self.n_groups_p].copy_from_slice(z_p);
        z[self.coef_start_q..self.coef_start_q + self.n_groups_q].copy_from_slice(z_q);
    }

    fn apply_transpose(&self, x: &[f64], z: &mut [f64]) {
        // Transpose of GS requires reversing sweep order (q then p)
        let mut scratch = self.scratch.lock().unwrap();
        let SparseGSScratch {
            z_p,
            z_q,
            temp_p,
            temp_q,
        } = &mut *scratch;

        let x_p = &x[self.coef_start_p..self.coef_start_p + self.n_groups_p];
        let x_q = &x[self.coef_start_q..self.coef_start_q + self.n_groups_q];

        z_p.iter_mut().for_each(|v| *v = 0.0);
        z_q.iter_mut().for_each(|v| *v = 0.0);

        // Reverse order: q then p
        for _ in 0..self.inner_iters {
            // sweep_q first
            self.mqp_matvec(z_p, temp_q);
            for ((z_i, &t), (&x_i, &inv_d)) in z_q
                .iter_mut()
                .zip(temp_q.iter())
                .zip(x_q.iter().zip(self.inv_diag_q.iter()))
            {
                *z_i = (x_i - t) * inv_d;
            }

            // sweep_p second
            self.mpq_matvec(z_q, temp_p);
            for ((z_i, &t), (&x_i, &inv_d)) in z_p
                .iter_mut()
                .zip(temp_p.iter())
                .zip(x_p.iter().zip(self.inv_diag_p.iter()))
            {
                *z_i = (x_i - t) * inv_d;
            }
        }

        z.copy_from_slice(x);
        z[self.coef_start_p..self.coef_start_p + self.n_groups_p].copy_from_slice(z_p);
        z[self.coef_start_q..self.coef_start_q + self.n_groups_q].copy_from_slice(z_q);
    }
}

// =============================================================================
// Preconditioner Configuration
// =============================================================================

/// Available preconditioner types for LSMR.
#[derive(Clone, Copy, Debug, Default)]
#[allow(dead_code)]
pub enum PreconditionerKind {
    /// No preconditioning (identity).
    #[default]
    None,

    /// Diagonal scaling by inverse square root of group counts (P1).
    Diagonal,

    /// Nullspace deflation (DEPRECATED - produces incorrect results).
    ///
    /// This preconditioner was intended to improve conditioning by projecting
    /// out the nullspace (constant mode) from each FE block. However, as a
    /// right preconditioner in LSMR, it constrains the solution to have zero
    /// weighted mean within each FE block, which changes the optimization
    /// problem and produces different (incorrect) demeaned values.
    ///
    /// Use `Diagonal` or `TwoBlockStreaming` instead.
    DiagonalPlusDeflation,

    /// Two-block streaming Gauss-Seidel (P3).
    /// Performs fixed inner GS iterations on a 2-FE subsystem.
    TwoBlockStreaming {
        /// First FE index (or None for auto-select)
        fe_p: Option<usize>,
        /// Second FE index (or None for auto-select)
        fe_q: Option<usize>,
        /// Number of inner Gauss-Seidel iterations
        inner_iters: usize,
    },

    /// Sparse Gram preconditioner using CG iterations (P4).
    /// Builds sparse 2x2 block Gram matrix and uses PCG to solve.
    SparseGram {
        /// First FE index (or None for auto-select)
        fe_p: Option<usize>,
        /// Second FE index (or None for auto-select)
        fe_q: Option<usize>,
        /// Number of inner CG iterations
        inner_iters: usize,
    },
}

/// Boxed preconditioner for runtime polymorphism.
pub type BoxedPreconditioner = Box<dyn RightPreconditioner + Send + Sync>;

/// Create a preconditioner from configuration.
pub fn create_preconditioner(ctx: &DemeanContext, kind: PreconditionerKind) -> BoxedPreconditioner {
    match kind {
        PreconditionerKind::None => Box::new(IdentityPreconditioner),
        PreconditionerKind::Diagonal => Box::new(DiagonalPreconditioner::new(ctx)),
        PreconditionerKind::DiagonalPlusDeflation => {
            // Use deflation alone - the composition with diagonal scaling was
            // mathematically incorrect because deflation needs to operate in the
            // original coefficient space to project out the correct nullspace.
            Box::new(NullspaceDeflator::new(ctx))
        }
        PreconditionerKind::TwoBlockStreaming {
            fe_p,
            fe_q,
            inner_iters,
        } => {
            // Use streaming GS alone - it already approximates M^{-1} directly
            // Composing with diagonal scaling changes the normal equations structure
            Box::new(TwoBlockStreamingPreconditioner::new(ctx, fe_p, fe_q, inner_iters))
        }
        PreconditionerKind::SparseGram {
            fe_p,
            fe_q,
            inner_iters,
        } => {
            // Use sparse Gram preconditioner with PCG inner solver
            Box::new(SparseGramPreconditioner::new(ctx, fe_p, fe_q, inner_iters))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_identity_preconditioner() {
        let precond = IdentityPreconditioner;
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let mut z = vec![0.0; 5];

        precond.apply(&x, &mut z);
        assert_eq!(z, x);
    }

    #[test]
    fn test_identity_is_symmetric() {
        let precond = IdentityPreconditioner;
        let x = vec![1.0, 2.0, 3.0];
        let mut z1 = vec![0.0; 3];
        let mut z2 = vec![0.0; 3];

        precond.apply(&x, &mut z1);
        precond.apply_transpose(&x, &mut z2);
        assert_eq!(z1, z2);
    }
}
