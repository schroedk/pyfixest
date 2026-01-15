//! Linear operator abstraction for LSMR solver.
//!
//! This module defines the `LinearOperator` trait that provides matrix-free
//! matrix-vector products needed by LSMR. The `DesignMatrixOperator` wraps
//! `DemeanContext` to provide these operations using the existing scatter/gather
//! infrastructure.
//!
//! # Weighted Least Squares
//!
//! For weighted least squares `min ||W^{1/2}(b - Dx)||²`, we transform the problem
//! to `min ||b̃ - Ãx||²` where:
//! - `Ã = W^{1/2} D` (effective matrix)
//! - `b̃ = W^{1/2} b` (effective RHS)
//!
//! This gives the correct normal equations `D^T W D x = D^T W b`.

use crate::demean::types::DemeanContext;

/// A matrix-free linear operator for LSMR.
///
/// Implementations must be allocation-free in `matvec`/`rmatvec`.
/// The caller provides all buffers.
///
/// For the fixed-effects design matrix D with shape (n_obs, n_coef):
/// - `matvec`: y = D * x (coefficient space → observation space)
/// - `rmatvec`: x = D^T * y (observation space → coefficient space)
pub trait LinearOperator {
    /// Number of rows (observation space dimension).
    fn rows(&self) -> usize;

    /// Number of columns (coefficient space dimension).
    fn cols(&self) -> usize;

    /// Compute y = A * x (forward multiplication).
    ///
    /// # Arguments
    /// * `x` - Input vector in column space (length: `cols()`)
    /// * `y` - Output vector in row space (length: `rows()`), overwritten
    fn matvec(&mut self, x: &[f64], y: &mut [f64]);

    /// Compute x = A^T * y (transpose multiplication).
    ///
    /// # Arguments
    /// * `y` - Input vector in row space (length: `rows()`)
    /// * `x` - Output vector in column space (length: `cols()`), overwritten
    fn rmatvec(&mut self, y: &[f64], x: &mut [f64]);
}

/// Design matrix operator wrapping `DemeanContext`.
///
/// For weighted problems, this implements `Ã = W^{1/2} D` where W = diag(weights).
/// For unweighted problems (ctx.weights = None), this is just D.
///
/// The operator provides:
/// - `matvec(x)`: `W^{1/2} D x` - gather coefficients, scale by sqrt(weights)
/// - `rmatvec(y)`: `D^T W^{1/2} y` - scale by sqrt(weights), scatter to coefficients
///
/// This ensures the LSMR algorithm solves the correct weighted least squares problem.
///
/// # Example
///
/// With 2 FEs (firm, year) and coefficients laid out as `[firm_0..firm_K, year_0..year_T]`:
/// - Unweighted: `y[i] = coef[firm[i]] + coef[K + year[i]]`
/// - Weighted: `y[i] = sqrt(w[i]) * (coef[firm[i]] + coef[K + year[i]])`
pub struct DesignMatrixOperator<'a> {
    ctx: &'a DemeanContext,
    /// Precomputed sqrt(weights) for weighted case, None for unweighted
    sqrt_weights: Option<Vec<f64>>,
}

impl<'a> DesignMatrixOperator<'a> {
    /// Create a new design matrix operator from a `DemeanContext`.
    ///
    /// Precomputes sqrt(weights) for efficient weighted operations.
    #[inline]
    pub fn new(ctx: &'a DemeanContext) -> Self {
        let sqrt_weights = ctx.weights.as_ref().map(|w| w.iter().map(|&wi| wi.sqrt()).collect());
        Self { ctx, sqrt_weights }
    }
}

impl LinearOperator for DesignMatrixOperator<'_> {
    #[inline]
    fn rows(&self) -> usize {
        self.ctx.dims.n_obs
    }

    #[inline]
    fn cols(&self) -> usize {
        self.ctx.dims.n_coef
    }

    #[inline]
    fn matvec(&mut self, x: &[f64], y: &mut [f64]) {
        debug_assert_eq!(x.len(), self.cols(), "x length mismatch");
        debug_assert_eq!(y.len(), self.rows(), "y length mismatch");

        // y = W^{1/2} D x: gather coefficients to observation space, scale by sqrt(weights)
        y.fill(0.0);
        // First gather: y = D * x
        self.ctx.apply_design_matrix(x, y);

        // Then scale by sqrt(weights) if weighted
        if let Some(ref sqrt_w) = self.sqrt_weights {
            for (yi, &swi) in y.iter_mut().zip(sqrt_w.iter()) {
                *yi *= swi;
            }
        }
    }

    #[inline]
    fn rmatvec(&mut self, y: &[f64], x: &mut [f64]) {
        debug_assert_eq!(y.len(), self.rows(), "y length mismatch");
        debug_assert_eq!(x.len(), self.cols(), "x length mismatch");

        // x = D^T W^{1/2} y: scale by sqrt(weights), then scatter to coefficient space
        match &self.sqrt_weights {
            None => {
                // Unweighted: x = D^T y (use the unweighted scatter)
                self.apply_design_matrix_t_unweighted(y, x);
            }
            Some(sqrt_w) => {
                // Weighted: x = D^T (W^{1/2} y) = sum_i sqrt(w_i) * y_i * indicator(group)
                self.apply_design_matrix_t_sqrt_weighted(y, sqrt_w, x);
            }
        }
    }
}

impl DesignMatrixOperator<'_> {
    /// Unweighted transpose: x = D^T y
    fn apply_design_matrix_t_unweighted(&self, y: &[f64], x: &mut [f64]) {
        x.fill(0.0);
        for fe in &self.ctx.fe_infos {
            let offset = fe.coef_start;
            for (i, &g) in fe.group_ids.iter().enumerate() {
                x[offset + g] += y[i];
            }
        }
    }

    /// Sqrt-weighted transpose: x = D^T W^{1/2} y
    fn apply_design_matrix_t_sqrt_weighted(&self, y: &[f64], sqrt_w: &[f64], x: &mut [f64]) {
        x.fill(0.0);
        for fe in &self.ctx.fe_infos {
            let offset = fe.coef_start;
            for (i, &g) in fe.group_ids.iter().enumerate() {
                x[offset + g] += y[i] * sqrt_w[i];
            }
        }
    }
}

/// Preconditioned operator: wraps A with right preconditioner M.
///
/// Represents the operator `A * M^{-1}` for right-preconditioned LSMR.
/// We solve `min ||A * M^{-1} * z - b||`, then recover `x = M^{-1} * z`.
pub struct PreconditionedOperator<'a, A> {
    operator: &'a mut A,
    preconditioner: &'a mut dyn crate::demean::lsmr::preconditioner::RightPreconditioner,
    /// Scratch buffer for M^{-1} * x (length: cols)
    scratch: Vec<f64>,
}

impl<'a, A: LinearOperator> PreconditionedOperator<'a, A> {
    /// Create a new preconditioned operator.
    pub fn new(
        operator: &'a mut A,
        preconditioner: &'a mut dyn crate::demean::lsmr::preconditioner::RightPreconditioner,
    ) -> Self {
        let scratch = vec![0.0; operator.cols()];
        Self {
            operator,
            preconditioner,
            scratch,
        }
    }
}

impl<A: LinearOperator> LinearOperator for PreconditionedOperator<'_, A> {
    #[inline]
    fn rows(&self) -> usize {
        self.operator.rows()
    }

    #[inline]
    fn cols(&self) -> usize {
        self.operator.cols()
    }

    fn matvec(&mut self, x: &[f64], y: &mut [f64]) {
        // y = A * M^{-1} * x
        // 1. scratch = M^{-1} * x
        // 2. y = A * scratch
        self.preconditioner.apply(x, &mut self.scratch);
        self.operator.matvec(&self.scratch, y);
    }

    fn rmatvec(&mut self, y: &[f64], x: &mut [f64]) {
        // x = (A * M^{-1})^T * y = M^{-T} * A^T * y
        // For symmetric M: M^{-T} = M^{-1}
        // 1. scratch = A^T * y
        // 2. x = M^{-T} * scratch
        self.operator.rmatvec(y, &mut self.scratch);
        self.preconditioner.apply_transpose(&self.scratch, x);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Simple identity-like operator for testing.
    struct IdentityOp {
        n: usize,
    }

    impl LinearOperator for IdentityOp {
        fn rows(&self) -> usize {
            self.n
        }
        fn cols(&self) -> usize {
            self.n
        }
        fn matvec(&mut self, x: &[f64], y: &mut [f64]) {
            y.copy_from_slice(x);
        }
        fn rmatvec(&mut self, y: &[f64], x: &mut [f64]) {
            x.copy_from_slice(y);
        }
    }

    #[test]
    fn test_identity_operator() {
        let mut op = IdentityOp { n: 5 };
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let mut y = vec![0.0; 5];

        op.matvec(&x, &mut y);
        assert_eq!(y, x);

        let mut x2 = vec![0.0; 5];
        op.rmatvec(&y, &mut x2);
        assert_eq!(x2, x);
    }
}
