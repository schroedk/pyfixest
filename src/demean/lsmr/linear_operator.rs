//! Linear operator abstraction for LSMR solver.
//!
//! This module defines the `LinearOperator` trait that provides matrix-free
//! matrix-vector products needed by LSMR. The `DesignMatrixOperator` wraps
//! `DemeanContext` to provide these operations using the existing scatter/gather
//! infrastructure.

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
    fn matvec(&self, x: &[f64], y: &mut [f64]);

    /// Compute x = A^T * y (transpose multiplication).
    ///
    /// # Arguments
    /// * `y` - Input vector in row space (length: `rows()`)
    /// * `x` - Output vector in column space (length: `cols()`), overwritten
    fn rmatvec(&self, y: &[f64], x: &mut [f64]);
}

/// Design matrix operator wrapping `DemeanContext`.
///
/// This adapts the existing gather/scatter operations as a `LinearOperator`.
/// The design matrix D has one-hot columns for each fixed effect group.
///
/// # Example
///
/// With 2 FEs (firm, year) and coefficients laid out as `[firm_0..firm_K, year_0..year_T]`:
/// - `matvec(coef, y)`: For each observation i, y[i] = coef[firm[i]] + coef[K + year[i]]
/// - `rmatvec(y, coef)`: Accumulate y values to their respective group coefficients
pub struct DesignMatrixOperator<'a> {
    ctx: &'a DemeanContext,
}

impl<'a> DesignMatrixOperator<'a> {
    /// Create a new design matrix operator from a `DemeanContext`.
    #[inline]
    pub fn new(ctx: &'a DemeanContext) -> Self {
        Self { ctx }
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
    fn matvec(&self, x: &[f64], y: &mut [f64]) {
        debug_assert_eq!(x.len(), self.cols(), "x length mismatch");
        debug_assert_eq!(y.len(), self.rows(), "y length mismatch");

        // y = D * x: gather coefficients to observation space
        y.fill(0.0);
        self.ctx.apply_design_matrix(x, y);
    }

    #[inline]
    fn rmatvec(&self, y: &[f64], x: &mut [f64]) {
        debug_assert_eq!(y.len(), self.rows(), "y length mismatch");
        debug_assert_eq!(x.len(), self.cols(), "x length mismatch");

        // x = D^T * y: scatter observations to coefficient space
        self.ctx.apply_design_matrix_t(y, x);
    }
}

/// Preconditioned operator: wraps A with right preconditioner M.
///
/// Represents the operator `A * M^{-1}` for right-preconditioned LSMR.
/// We solve `min ||A * M^{-1} * z - b||`, then recover `x = M^{-1} * z`.
///
/// NOTE: This struct is currently unused - it will be activated when
/// full preconditioning support is implemented in later phases (P2/P3).
#[allow(dead_code)]
pub struct PreconditionedOperator<'a, A, M> {
    operator: &'a A,
    preconditioner: &'a M,
    /// Scratch buffer for M^{-1} * x (length: cols)
    /// Uses RefCell for interior mutability since LinearOperator::matvec takes &self
    scratch: std::cell::RefCell<Vec<f64>>,
}

#[allow(dead_code)]
impl<'a, A: LinearOperator, M: crate::demean::lsmr::preconditioner::RightPreconditioner>
    PreconditionedOperator<'a, A, M>
{
    /// Create a new preconditioned operator.
    pub fn new(operator: &'a A, preconditioner: &'a M) -> Self {
        let scratch = std::cell::RefCell::new(vec![0.0; operator.cols()]);
        Self {
            operator,
            preconditioner,
            scratch,
        }
    }
}

impl<A: LinearOperator, M: crate::demean::lsmr::preconditioner::RightPreconditioner> LinearOperator
    for PreconditionedOperator<'_, A, M>
{
    #[inline]
    fn rows(&self) -> usize {
        self.operator.rows()
    }

    #[inline]
    fn cols(&self) -> usize {
        self.operator.cols()
    }

    fn matvec(&self, x: &[f64], y: &mut [f64]) {
        // y = A * M^{-1} * x
        // 1. scratch = M^{-1} * x
        // 2. y = A * scratch
        let mut scratch = self.scratch.borrow_mut();
        self.preconditioner.apply(x, &mut scratch);
        self.operator.matvec(&scratch, y);
    }

    fn rmatvec(&self, y: &[f64], x: &mut [f64]) {
        // x = (A * M^{-1})^T * y = M^{-T} * A^T * y
        // For symmetric M: M^{-T} = M^{-1}
        // 1. scratch = A^T * y
        // 2. x = M^{-T} * scratch
        let mut scratch = self.scratch.borrow_mut();
        self.operator.rmatvec(y, &mut scratch);
        self.preconditioner.apply_transpose(&scratch, x);
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
        fn matvec(&self, x: &[f64], y: &mut [f64]) {
            y.copy_from_slice(x);
        }
        fn rmatvec(&self, y: &[f64], x: &mut [f64]) {
            x.copy_from_slice(y);
        }
    }

    #[test]
    fn test_identity_operator() {
        let op = IdentityOp { n: 5 };
        let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let mut y = vec![0.0; 5];

        op.matvec(&x, &mut y);
        assert_eq!(y, x);

        let mut x2 = vec![0.0; 5];
        op.rmatvec(&y, &mut x2);
        assert_eq!(x2, x);
    }
}
