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
            // TODO: Implement in Phase 3 (P2 task)
            // For now, fall back to diagonal only
            Box::new(DiagonalPreconditioner::new(ctx))
        }
        PreconditionerKind::TwoBlockStreaming { .. } => {
            // TODO: Implement in Phase 4 (P3 task)
            // For now, fall back to diagonal
            Box::new(DiagonalPreconditioner::new(ctx))
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
