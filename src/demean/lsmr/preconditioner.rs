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

    /// Diagonal + nullspace deflation (P1 + P2).
    /// Removes constant mode from each FE block.
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
}

/// Boxed preconditioner for runtime polymorphism.
pub type BoxedPreconditioner = Box<dyn RightPreconditioner + Send + Sync>;

/// Create a preconditioner from configuration.
pub fn create_preconditioner(ctx: &DemeanContext, kind: PreconditionerKind) -> BoxedPreconditioner {
    match kind {
        PreconditionerKind::None => Box::new(IdentityPreconditioner),
        PreconditionerKind::Diagonal => Box::new(DiagonalPreconditioner::new(ctx)),
        PreconditionerKind::DiagonalPlusDeflation => {
            // P1 + P2: Diagonal scaling followed by nullspace deflation
            let p1 = DiagonalPreconditioner::new(ctx);
            let p2 = NullspaceDeflator::new(ctx);
            Box::new(ComposedPreconditioner::new(p1, p2, ctx.dims.n_coef))
        }
        PreconditionerKind::TwoBlockStreaming { .. } => {
            // TODO: Implement in Phase 4 (P3 task)
            // For now, fall back to diagonal + deflation
            let p1 = DiagonalPreconditioner::new(ctx);
            let p2 = NullspaceDeflator::new(ctx);
            Box::new(ComposedPreconditioner::new(p1, p2, ctx.dims.n_coef))
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
