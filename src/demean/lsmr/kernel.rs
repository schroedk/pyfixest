//! LSMR algorithm implementation.
//!
//! This module implements the LSMR (Least Squares Minimal Residual) algorithm
//! by Fong and Saunders (2011). LSMR uses Golub-Kahan bidiagonalization to
//! solve least-squares problems without forming the normal equations.
//!
//! # Algorithm Overview
//!
//! LSMR solves: min ||A x - b||₂
//!
//! 1. **Bidiagonalization**: Build orthonormal bases U and V such that
//!    A V = U B where B is lower bidiagonal
//! 2. **QR factorization**: Apply Givens rotations to B
//! 3. **Solution update**: Incrementally build x from the factorization
//!
//! # References
//!
//! - Fong, D. C.-L., & Saunders, M. A. (2011). LSMR: An iterative algorithm
//!   for sparse least-squares problems. SIAM J. Sci. Comput., 33(5), 2950-2971.

use super::buffers::{LSMRBuffers, LSMRState};
use super::linear_operator::LinearOperator;

/// Result of an LSMR solve.
#[derive(Clone, Copy, Debug)]
pub struct LSMRResult {
    /// Number of iterations performed.
    pub iterations: usize,
    /// Whether the algorithm converged.
    pub converged: bool,
    /// Final residual norm ||b - Ax||.
    #[allow(dead_code)]
    pub residual_norm: f64,
    /// Final ||A^T r|| (gradient norm).
    #[allow(dead_code)]
    pub atr_norm: f64,
}

/// LSMR configuration parameters.
#[derive(Clone, Copy, Debug)]
pub struct LSMRConfig {
    /// Convergence tolerance for ||A^T r|| / (||A|| ||r||).
    pub tol: f64,
    /// Maximum number of iterations.
    pub maxiter: usize,
    /// Tolerance for condition number estimate (optional early stop).
    pub conlim: f64,
}

impl Default for LSMRConfig {
    fn default() -> Self {
        Self {
            tol: 1e-8,
            maxiter: 10_000,
            conlim: 1e8,
        }
    }
}

/// LSMR solver kernel.
///
/// Implements the core LSMR algorithm using Golub-Kahan bidiagonalization.
/// Uses pre-allocated buffers for zero-allocation iteration.
pub struct LSMRKernel<'a> {
    config: LSMRConfig,
    buffers: &'a mut LSMRBuffers,
    state: LSMRState,
}

impl<'a> LSMRKernel<'a> {
    /// Create a new LSMR kernel with the given configuration and buffers.
    pub fn new(config: LSMRConfig, buffers: &'a mut LSMRBuffers) -> Self {
        Self {
            config,
            buffers,
            state: LSMRState::new(),
        }
    }

    /// Solve min ||A x - b||₂.
    ///
    /// # Arguments
    /// * `operator` - Linear operator representing matrix A
    /// * `b` - Right-hand side vector (length: n_obs)
    /// * `x` - Output solution vector (length: n_coef), overwritten
    ///
    /// # Returns
    /// `LSMRResult` with convergence information.
    pub fn solve<A: LinearOperator>(&mut self, operator: &mut A, b: &[f64], x: &mut [f64]) -> LSMRResult {
        let n_obs = operator.rows();
        let n_coef = operator.cols();

        debug_assert_eq!(b.len(), n_obs);
        debug_assert_eq!(x.len(), n_coef);
        debug_assert_eq!(self.buffers.n_obs(), n_obs);
        debug_assert_eq!(self.buffers.n_coef(), n_coef);

        // Reset state for new solve
        self.buffers.reset();
        self.state.reset();

        // Initialize: beta_1 * u_1 = b
        self.buffers.u.copy_from_slice(b);
        self.state.beta = norm2(&self.buffers.u);

        if self.state.beta > 0.0 {
            scale_inplace(&mut self.buffers.u, 1.0 / self.state.beta);
        }

        // alpha_1 * v_1 = A^T * u_1
        operator.rmatvec(&self.buffers.u, &mut self.buffers.v);
        self.state.alpha = norm2(&self.buffers.v);

        if self.state.alpha > 0.0 {
            scale_inplace(&mut self.buffers.v, 1.0 / self.state.alpha);
        }

        // Initialize w = v
        self.buffers.w.copy_from_slice(&self.buffers.v);

        // Initialize scalars for QR factorization
        self.state.rho_bar = self.state.alpha;
        self.state.phi_bar = self.state.beta;
        self.state.rnorm = self.state.beta;
        self.state.arnorm = self.state.alpha * self.state.beta;

        // Estimate of ||A||_F
        self.state.anorm = self.state.alpha;

        // Check for trivial case (zero RHS)
        if self.state.beta == 0.0 {
            x.fill(0.0);
            return LSMRResult {
                iterations: 0,
                converged: true,
                residual_norm: 0.0,
                atr_norm: 0.0,
            };
        }

        // Check for immediate convergence
        if self.state.arnorm == 0.0 {
            x.fill(0.0);
            return LSMRResult {
                iterations: 0,
                converged: true,
                residual_norm: self.state.rnorm,
                atr_norm: 0.0,
            };
        }

        // Main iteration loop
        for iter in 1..=self.config.maxiter {
            // Bidiagonalization step
            self.bidiagonalization_step(operator);

            // Apply Givens rotation and update solution
            self.qr_and_update_step();

            // Update norms for convergence check
            self.update_norms();

            // Check convergence
            let test1 = self.state.rnorm / self.state.beta; // ||r|| / ||b||
            let test2 = if self.state.anorm * self.state.rnorm > 0.0 {
                self.state.arnorm / (self.state.anorm * self.state.rnorm)
            } else {
                0.0
            };

            if test2 <= self.config.tol || test1 <= self.config.tol {
                // Copy solution to output
                x.copy_from_slice(&self.buffers.x);
                return LSMRResult {
                    iterations: iter,
                    converged: true,
                    residual_norm: self.state.rnorm,
                    atr_norm: self.state.arnorm,
                };
            }

            // Check condition number limit
            if self.state.acond >= self.config.conlim {
                x.copy_from_slice(&self.buffers.x);
                return LSMRResult {
                    iterations: iter,
                    converged: false,
                    residual_norm: self.state.rnorm,
                    atr_norm: self.state.arnorm,
                };
            }
        }

        // Max iterations reached
        x.copy_from_slice(&self.buffers.x);
        LSMRResult {
            iterations: self.config.maxiter,
            converged: false,
            residual_norm: self.state.rnorm,
            atr_norm: self.state.arnorm,
        }
    }

    /// One step of Golub-Kahan bidiagonalization.
    #[inline]
    fn bidiagonalization_step<A: LinearOperator>(&mut self, operator: &mut A) {
        // u = A * v - alpha * u
        operator.matvec(&self.buffers.v, &mut self.buffers.matvec_scratch);
        for (u_i, &s_i) in self.buffers.u.iter_mut().zip(self.buffers.matvec_scratch.iter()) {
            *u_i = s_i - self.state.alpha * *u_i;
        }

        // beta = ||u||
        self.state.beta = norm2(&self.buffers.u);

        if self.state.beta > 0.0 {
            scale_inplace(&mut self.buffers.u, 1.0 / self.state.beta);

            // Update ||A||_F estimate
            self.state.anorm = (self.state.anorm.powi(2)
                + self.state.alpha.powi(2)
                + self.state.beta.powi(2))
            .sqrt();

            // v = A^T * u - beta * v
            operator.rmatvec(&self.buffers.u, &mut self.buffers.precond_scratch);
            for (v_i, &s_i) in self.buffers.v.iter_mut().zip(self.buffers.precond_scratch.iter()) {
                *v_i = s_i - self.state.beta * *v_i;
            }

            // alpha = ||v||
            self.state.alpha = norm2(&self.buffers.v);

            if self.state.alpha > 0.0 {
                scale_inplace(&mut self.buffers.v, 1.0 / self.state.alpha);
            }
        }
    }

    /// Apply Givens rotation and update solution estimate.
    #[inline]
    fn qr_and_update_step(&mut self) {
        // Construct and apply Givens rotation
        let rho = (self.state.rho_bar.powi(2) + self.state.beta.powi(2)).sqrt();

        let c = self.state.rho_bar / rho;
        let s = self.state.beta / rho;

        let theta_new = s * self.state.alpha;
        let rho_bar_new = -c * self.state.alpha;
        let phi = c * self.state.phi_bar;
        let phi_bar_new = s * self.state.phi_bar;

        // Fused update: x = x + (phi/rho)*w and w = v - (theta_new/rho)*w
        if rho.abs() > 1e-15 {
            let factor_x = phi / rho;
            let factor_w = theta_new / rho;
            for ((x_i, w_i), &v_i) in self.buffers.x.iter_mut()
                .zip(self.buffers.w.iter_mut())
                .zip(self.buffers.v.iter())
            {
                *x_i += factor_x * *w_i;
                *w_i = v_i - factor_w * *w_i;
            }
        } else {
            self.buffers.w.copy_from_slice(&self.buffers.v);
        }

        // Update state
        self.state.rho_bar = rho_bar_new;
        self.state.phi_bar = phi_bar_new;
        self.state.c = c;
        self.state.s = s;
    }

    /// Update norm estimates for convergence checking.
    #[inline]
    fn update_norms(&mut self) {
        // ||r|| estimate from QR factorization
        self.state.rnorm = self.state.phi_bar.abs();

        // ||A^T r|| estimate
        self.state.arnorm = (self.state.alpha * self.state.c.abs() * self.state.phi_bar).abs();

        // Condition number estimate
        if self.state.anorm > 0.0 && self.state.rho_bar.abs() > 0.0 {
            self.state.acond = self.state.anorm / self.state.rho_bar.abs();
        }
    }
}

// =============================================================================
// Utility Functions
// =============================================================================

/// Compute 2-norm of a vector.
#[inline]
fn norm2(x: &[f64]) -> f64 {
    x.iter().map(|&v| v * v).sum::<f64>().sqrt()
}

/// Scale vector in place: x = alpha * x.
#[inline]
fn scale_inplace(x: &mut [f64], alpha: f64) {
    for xi in x.iter_mut() {
        *xi *= alpha;
    }
}

/// Compute dot product of two vectors.
#[allow(dead_code)]
#[inline]
fn dot(x: &[f64], y: &[f64]) -> f64 {
    x.iter().zip(y.iter()).map(|(&a, &b)| a * b).sum()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::demean::lsmr::linear_operator::LinearOperator;

    /// Simple diagonal operator for testing: A = diag(d)
    struct DiagonalOp {
        diag: Vec<f64>,
    }

    impl LinearOperator for DiagonalOp {
        fn rows(&self) -> usize {
            self.diag.len()
        }
        fn cols(&self) -> usize {
            self.diag.len()
        }
        fn matvec(&mut self, x: &[f64], y: &mut [f64]) {
            for ((y_i, &x_i), &d_i) in y.iter_mut().zip(x.iter()).zip(self.diag.iter()) {
                *y_i = d_i * x_i;
            }
        }
        fn rmatvec(&mut self, y: &[f64], x: &mut [f64]) {
            // Diagonal is symmetric
            self.matvec(y, x);
        }
    }

    #[test]
    fn test_lsmr_identity() {
        // Solve Ix = b where I is identity
        let n = 5;
        let mut op = DiagonalOp {
            diag: vec![1.0; n],
        };
        let b = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let mut x = vec![0.0; n];

        let config = LSMRConfig::default();
        let mut buffers = LSMRBuffers::new(n, n);
        let mut kernel = LSMRKernel::new(config, &mut buffers);

        let result = kernel.solve(&mut op, &b, &mut x);

        assert!(result.converged);
        for (i, (&x_i, &b_i)) in x.iter().zip(b.iter()).enumerate() {
            assert!(
                (x_i - b_i).abs() < 1e-10,
                "x[{}] = {} != b[{}] = {}",
                i,
                x_i,
                i,
                b_i
            );
        }
    }

    #[test]
    fn test_lsmr_diagonal() {
        // Solve diag(2,3,4) x = [2, 6, 12] → x = [1, 2, 3]
        let mut op = DiagonalOp {
            diag: vec![2.0, 3.0, 4.0],
        };
        let b = vec![2.0, 6.0, 12.0];
        let expected = vec![1.0, 2.0, 3.0];
        let mut x = vec![0.0; 3];

        let config = LSMRConfig::default();
        let mut buffers = LSMRBuffers::new(3, 3);
        let mut kernel = LSMRKernel::new(config, &mut buffers);

        let result = kernel.solve(&mut op, &b, &mut x);

        assert!(result.converged);
        for (i, (&x_i, &e_i)) in x.iter().zip(expected.iter()).enumerate() {
            assert!(
                (x_i - e_i).abs() < 1e-8,
                "x[{}] = {} != expected[{}] = {}",
                i,
                x_i,
                i,
                e_i
            );
        }
    }

    #[test]
    fn test_lsmr_zero_rhs() {
        // Solve Ax = 0 → x = 0
        let mut op = DiagonalOp {
            diag: vec![1.0, 2.0, 3.0],
        };
        let b = vec![0.0, 0.0, 0.0];
        let mut x = vec![1.0, 1.0, 1.0]; // Non-zero initial

        let config = LSMRConfig::default();
        let mut buffers = LSMRBuffers::new(3, 3);
        let mut kernel = LSMRKernel::new(config, &mut buffers);

        let result = kernel.solve(&mut op, &b, &mut x);

        assert!(result.converged);
        assert_eq!(result.iterations, 0);
        for (i, &x_i) in x.iter().enumerate() {
            assert!(x_i.abs() < 1e-10, "x[{}] = {} != 0", i, x_i);
        }
    }
}
