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

use crate::demean::sweep::TwoFESweeper;
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
    fn apply(&mut self, x: &[f64], z: &mut [f64]);

    /// Apply the transpose: z = M^{-T} * x.
    ///
    /// Default implementation assumes symmetric preconditioner (M^{-T} = M^{-1}).
    /// Override for non-symmetric preconditioners.
    #[allow(dead_code)]
    fn apply_transpose(&mut self, x: &[f64], z: &mut [f64]) {
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
    fn apply(&mut self, x: &[f64], z: &mut [f64]) {
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
    fn apply(&mut self, x: &[f64], z: &mut [f64]) {
        debug_assert_eq!(x.len(), self.inv_sqrt_weights.len());
        debug_assert_eq!(z.len(), self.inv_sqrt_weights.len());

        for ((z_i, &x_i), &w) in z.iter_mut().zip(x.iter()).zip(self.inv_sqrt_weights.iter()) {
            *z_i = x_i * w;
        }
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
pub struct TwoBlockStreamingPreconditioner<'a> {
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
    /// Sweeper for FE p (computes p coefficients from q)
    sweeper_p: TwoFESweeper<'a>,
    /// Sweeper for FE q (computes q coefficients from p)
    sweeper_q: TwoFESweeper<'a>,
    /// Preallocated scratch buffers
    scratch_p: Vec<f64>,
    scratch_q: Vec<f64>,
}

impl<'a> TwoBlockStreamingPreconditioner<'a> {
    /// Create a two-block streaming preconditioner.
    ///
    /// # Arguments
    /// * `ctx` - DemeanContext with FE information
    /// * `fe_p` - First FE index (or None for auto-select)
    /// * `fe_q` - Second FE index (or None for auto-select)
    /// * `inner_iters` - Number of inner GS iterations (typically 2-5)
    pub fn new(
        ctx: &'a DemeanContext,
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
        let weights_ptr = ctx.weights.as_ref().map(|w| w.as_ptr());

        Self {
            inner_iters,
            coef_start_p: fe_p_info.coef_start,
            coef_start_q: fe_q_info.coef_start,
            n_groups_p,
            n_groups_q,
            sweeper_p: TwoFESweeper::new(ctx.dims.n_obs, weights_ptr, fe_p_info, fe_q_info),
            sweeper_q: TwoFESweeper::new(ctx.dims.n_obs, weights_ptr, fe_q_info, fe_p_info),
            scratch_p: vec![0.0; n_groups_p],
            scratch_q: vec![0.0; n_groups_q],
        }
    }

    /// Perform one GS sweep: update FE p coefficients given FE q coefficients.
    fn sweep_p(&mut self, x: &[f64], z: &mut [f64]) {
        let rhs = &x[self.coef_start_p..self.coef_start_p + self.n_groups_p];
        let other = &z[self.coef_start_q..self.coef_start_q + self.n_groups_q];

        self.sweeper_p.sweep(rhs, other, &mut self.scratch_p);

        z[self.coef_start_p..self.coef_start_p + self.n_groups_p].copy_from_slice(&self.scratch_p);
    }

    /// Perform one GS sweep: update FE q coefficients given FE p coefficients.
    fn sweep_q(&mut self, x: &[f64], z: &mut [f64]) {
        let rhs = &x[self.coef_start_q..self.coef_start_q + self.n_groups_q];
        let other = &z[self.coef_start_p..self.coef_start_p + self.n_groups_p];

        self.sweeper_q.sweep(rhs, other, &mut self.scratch_q);

        z[self.coef_start_q..self.coef_start_q + self.n_groups_q].copy_from_slice(&self.scratch_q);
    }
}

impl RightPreconditioner for TwoBlockStreamingPreconditioner<'_> {
    fn apply(&mut self, x: &[f64], z: &mut [f64]) {
        // Copy non-P/Q coefficients from x, zero P/Q for GS iteration
        z.copy_from_slice(x);
        z[self.coef_start_p..self.coef_start_p + self.n_groups_p].fill(0.0);
        z[self.coef_start_q..self.coef_start_q + self.n_groups_q].fill(0.0);

        // Perform inner_iters Gauss-Seidel iterations on the 2-FE subsystem
        // Each sweep uses the original input x as the RHS
        // Forward order: p then q
        for _ in 0..self.inner_iters {
            self.sweep_p(x, z);
            self.sweep_q(x, z);
        }
    }

    fn apply_transpose(&mut self, x: &[f64], z: &mut [f64]) {
        // The transpose of GS iteration requires reversing the sweep order.
        // If forward is: sweep_p then sweep_q (inner loop)
        // Then transpose is: sweep_q then sweep_p
        //
        // This is because for Gauss-Seidel solving M z = x where M = L + D + U:
        // Forward GS uses (L + D)^{-1}
        // Backward GS (transpose) uses (D + U)^{-1} = (L^T + D)^{-1}
        // which corresponds to reversing the update order.

        // Copy non-P/Q coefficients from x, zero P/Q for GS iteration
        z.copy_from_slice(x);
        z[self.coef_start_p..self.coef_start_p + self.n_groups_p].fill(0.0);
        z[self.coef_start_q..self.coef_start_q + self.n_groups_q].fill(0.0);

        // Reverse order: q then p
        for _ in 0..self.inner_iters {
            self.sweep_q(x, z);
            self.sweep_p(x, z);
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
    /// Scratch buffers
    scratch: SparseGSScratch,
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
            scratch,
        }
    }

}

/// Sparse matrix-vector product: y = M @ x (CSR format)
#[inline]
fn sparse_csr_matvec(
    row_ptr: &[usize],
    col_idx: &[usize],
    values: &[f64],
    x: &[f64],
    y: &mut [f64],
) {
    for (i, y_i) in y.iter_mut().enumerate() {
        let start = row_ptr[i];
        let end = row_ptr[i + 1];
        let mut sum = 0.0;
        for idx in start..end {
            let j = col_idx[idx];
            sum += values[idx] * x[j];
        }
        *y_i = sum;
    }
}

impl SparseGramPreconditioner {
    /// Perform Gauss-Seidel sweep with specified order.
    ///
    /// * `p_then_q = true`: sweep p first, then q (forward)
    /// * `p_then_q = false`: sweep q first, then p (transpose)
    fn gs_sweep(&mut self, x_p: &[f64], x_q: &[f64], p_then_q: bool) {
        // Initialize z_p, z_q = 0
        self.scratch.z_p.iter_mut().for_each(|v| *v = 0.0);
        self.scratch.z_q.iter_mut().for_each(|v| *v = 0.0);

        // Perform inner_iters Gauss-Seidel iterations using sparse matvec
        // Solving: diag(count_p) @ z_p + M_pq @ z_q = x_p
        //          M_qp @ z_p + diag(count_q) @ z_q = x_q
        for _ in 0..self.inner_iters {
            if p_then_q {
                self.sweep_p(x_p);
                self.sweep_q(x_q);
            } else {
                self.sweep_q(x_q);
                self.sweep_p(x_p);
            }
        }
    }

    /// Sweep p: z_p = (x_p - M_pq @ z_q) * inv_diag_p
    #[inline]
    fn sweep_p(&mut self, x_p: &[f64]) {
        sparse_csr_matvec(
            &self.mpq_row_ptr,
            &self.mpq_col_idx,
            &self.mpq_values,
            &self.scratch.z_q,
            &mut self.scratch.temp_p,
        );
        for ((z_i, &t), (&x_i, &inv_d)) in self.scratch.z_p
            .iter_mut()
            .zip(self.scratch.temp_p.iter())
            .zip(x_p.iter().zip(self.inv_diag_p.iter()))
        {
            *z_i = (x_i - t) * inv_d;
        }
    }

    /// Sweep q: z_q = (x_q - M_qp @ z_p) * inv_diag_q
    #[inline]
    fn sweep_q(&mut self, x_q: &[f64]) {
        sparse_csr_matvec(
            &self.mqp_row_ptr,
            &self.mqp_col_idx,
            &self.mqp_values,
            &self.scratch.z_p,
            &mut self.scratch.temp_q,
        );
        for ((z_i, &t), (&x_i, &inv_d)) in self.scratch.z_q
            .iter_mut()
            .zip(self.scratch.temp_q.iter())
            .zip(x_q.iter().zip(self.inv_diag_q.iter()))
        {
            *z_i = (x_i - t) * inv_d;
        }
    }

    /// Copy sweep results back to output vector.
    #[inline]
    fn copy_results(&self, x: &[f64], z: &mut [f64]) {
        z.copy_from_slice(x);
        z[self.coef_start_p..self.coef_start_p + self.n_groups_p].copy_from_slice(&self.scratch.z_p);
        z[self.coef_start_q..self.coef_start_q + self.n_groups_q].copy_from_slice(&self.scratch.z_q);
    }
}

impl RightPreconditioner for SparseGramPreconditioner {
    fn apply(&mut self, x: &[f64], z: &mut [f64]) {
        let x_p = &x[self.coef_start_p..self.coef_start_p + self.n_groups_p];
        let x_q = &x[self.coef_start_q..self.coef_start_q + self.n_groups_q];
        self.gs_sweep(x_p, x_q, true); // p then q
        self.copy_results(x, z);
    }

    fn apply_transpose(&mut self, x: &[f64], z: &mut [f64]) {
        let x_p = &x[self.coef_start_p..self.coef_start_p + self.n_groups_p];
        let x_q = &x[self.coef_start_q..self.coef_start_q + self.n_groups_q];
        self.gs_sweep(x_p, x_q, false); // q then p (reversed for transpose)
        self.copy_results(x, z);
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
pub type BoxedPreconditioner<'a> = Box<dyn RightPreconditioner + Send + Sync + 'a>;

/// Create a preconditioner from configuration.
pub fn create_preconditioner<'a>(
    ctx: &'a DemeanContext,
    kind: PreconditionerKind,
) -> BoxedPreconditioner<'a> {
    match kind {
        PreconditionerKind::None => Box::new(IdentityPreconditioner),
        PreconditionerKind::Diagonal => Box::new(DiagonalPreconditioner::new(ctx)),
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
        let mut precond = IdentityPreconditioner;
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let mut z = vec![0.0; 5];

        precond.apply(&x, &mut z);
        assert_eq!(z, x);
    }

    #[test]
    fn test_identity_is_symmetric() {
        let mut precond = IdentityPreconditioner;
        let x = vec![1.0, 2.0, 3.0];
        let mut z1 = vec![0.0; 3];
        let mut z2 = vec![0.0; 3];

        precond.apply(&x, &mut z1);
        precond.apply_transpose(&x, &mut z2);
        assert_eq!(z1, z2);
    }
}
