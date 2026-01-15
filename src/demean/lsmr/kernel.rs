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
    /// Damping (regularization) parameter λ.
    /// Solves min ||Ax - b||² + λ²||x||² instead of min ||Ax - b||².
    /// Set to 0.0 for standard unregularized least squares.
    pub damp: f64,
}

impl Default for LSMRConfig {
    fn default() -> Self {
        Self {
            tol: 1e-8,
            maxiter: 10_000,
            conlim: 1e8,
            damp: 0.0,
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
        self.state.alphabar = self.state.alpha; // For damped LSMR
        self.state.rho_bar = self.state.alpha;
        self.state.phi_bar = self.state.beta;
        self.state.rnorm = self.state.beta;
        self.state.arnorm = self.state.alpha * self.state.beta;

        // Estimate of ||A||_F (includes damp for regularized case)
        let damp = self.config.damp;
        self.state.anorm = (self.state.alpha * self.state.alpha + damp * damp).sqrt();

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

    /// Apply Givens rotations and update solution estimate.
    ///
    /// For damped LSMR, we apply two rotations:
    /// 1. First rotation handles damping: eliminates damp from [alphabar; damp]
    /// 2. Second rotation: eliminates beta from the result
    ///
    /// When damp=0, the first rotation is identity (chat=1) and we recover
    /// the standard LSMR algorithm.
    #[inline]
    fn qr_and_update_step(&mut self) {
        let damp = self.config.damp;

        // First rotation (for damping): [alphabar; damp] -> [alphahat; 0]
        // When damp=0: chat=1, alphahat=alphabar (no change)
        let (chat, alphahat) = if damp == 0.0 {
            (1.0, self.state.alphabar)
        } else {
            let (c, _s, r) = sym_ortho(self.state.alphabar, damp);
            (c, r)
        };

        // Second rotation: [alphahat; beta] -> [rho; 0]
        // Use direct computation to preserve sign of alphahat
        let rho = (alphahat * alphahat + self.state.beta * self.state.beta).sqrt();
        let (c, s) = if rho > 1e-15 {
            (alphahat / rho, self.state.beta / rho)
        } else {
            (1.0, 0.0)
        };

        let theta_new = s * self.state.alpha;
        let rho_bar_new = -c * self.state.alpha;
        let phi = c * chat * self.state.phi_bar;
        let phi_bar_new = s * self.state.phi_bar;

        // Fused update: x = x + (phi/rho)*w and w = v - (theta_new/rho)*w
        if rho > 1e-15 {
            let factor_x = phi / rho;
            let factor_w = theta_new / rho;
            for ((x_i, w_i), &v_i) in self
                .buffers
                .x
                .iter_mut()
                .zip(self.buffers.w.iter_mut())
                .zip(self.buffers.v.iter())
            {
                *x_i += factor_x * *w_i;
                *w_i = v_i - factor_w * *w_i;
            }
        } else {
            self.buffers.w.copy_from_slice(&self.buffers.v);
        }

        // Update state for next iteration
        // alphabar carries across iterations; when damp=0, alphabar = rho_bar (with sign)
        self.state.alphabar = rho_bar_new;
        self.state.rho_bar = rho_bar_new;
        self.state.phi_bar = phi_bar_new;
        self.state.c = c;
        self.state.s = s;
    }

    /// Update norm estimates for convergence checking.
    #[inline]
    fn update_norms(&mut self) {
        let damp = self.config.damp;

        // Update ||A||_F estimate (includes damp for regularized case)
        self.state.anorm = (self.state.anorm.powi(2)
            + self.state.alpha.powi(2)
            + self.state.beta.powi(2)
            + damp * damp)
            .sqrt();

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

/// Numerically stable Givens rotation.
///
/// Computes (c, s, r) such that:
/// ```text
/// [ c  s ] [ a ]   [ r ]
/// [-s  c ] [ b ] = [ 0 ]
/// ```
///
/// This implementation follows S.-C. Choi's stable formulation which avoids
/// potential division by very small numbers.
///
/// # References
/// - S.-C. Choi, "Iterative Methods for Singular Linear Equations and
///   Least-Squares Problems", Stanford PhD Thesis, 2006.
#[inline]
fn sym_ortho(a: f64, b: f64) -> (f64, f64, f64) {
    if b == 0.0 {
        let c = if a >= 0.0 { 1.0 } else { -1.0 };
        (c, 0.0, a.abs())
    } else if a == 0.0 {
        let s = if b >= 0.0 { 1.0 } else { -1.0 };
        (0.0, s, b.abs())
    } else if b.abs() > a.abs() {
        let tau = a / b;
        let s = (if b >= 0.0 { 1.0 } else { -1.0 }) / (1.0 + tau * tau).sqrt();
        let c = s * tau;
        let r = b / s;
        (c, s, r)
    } else {
        let tau = b / a;
        let c = (if a >= 0.0 { 1.0 } else { -1.0 }) / (1.0 + tau * tau).sqrt();
        let s = c * tau;
        let r = a / c;
        (c, s, r)
    }
}

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

    #[test]
    fn test_sym_ortho_basic() {
        // Test basic Givens rotation: should satisfy c² + s² = 1 and r = √(a² + b²)
        let (c, s, r) = sym_ortho(3.0, 4.0);
        assert!((c * c + s * s - 1.0).abs() < 1e-14, "c²+s² should be 1");
        assert!((r - 5.0).abs() < 1e-14, "r should be 5");
        // Check rotation property: c*a + s*b = r
        assert!((c * 3.0 + s * 4.0 - r).abs() < 1e-14);
    }

    #[test]
    fn test_sym_ortho_edge_cases() {
        // b = 0 case
        let (c, s, r) = sym_ortho(5.0, 0.0);
        assert_eq!(c, 1.0);
        assert_eq!(s, 0.0);
        assert_eq!(r, 5.0);

        // a = 0 case
        let (c, s, r) = sym_ortho(0.0, 5.0);
        assert_eq!(c, 0.0);
        assert_eq!(s, 1.0);
        assert_eq!(r, 5.0);

        // Negative values
        let (c, s, r) = sym_ortho(-3.0, 4.0);
        assert!((c * c + s * s - 1.0).abs() < 1e-14);
        assert!((r - 5.0).abs() < 1e-14);

        // |b| > |a| case (triggers different branch)
        let (c, s, r) = sym_ortho(1.0, 10.0);
        assert!((c * c + s * s - 1.0).abs() < 1e-14);
        assert!((r - (101.0_f64).sqrt()).abs() < 1e-14);
    }

    #[test]
    fn test_lsmr_damped_identity() {
        // Solve (I + λ²I)x = b with damping λ
        // For diagonal A=I, the damped solution is: x = b / (1 + λ²)
        let n = 5;
        let mut op = DiagonalOp {
            diag: vec![1.0; n],
        };
        let b = vec![1.0, 2.0, 3.0, 4.0, 5.0];

        // Without damping
        let config_undamped = LSMRConfig::default();
        let mut buffers = LSMRBuffers::new(n, n);
        let mut x_undamped = vec![0.0; n];
        let mut kernel = LSMRKernel::new(config_undamped, &mut buffers);
        kernel.solve(&mut op, &b, &mut x_undamped);

        // With damping λ = 1.0
        let config_damped = LSMRConfig {
            damp: 1.0,
            ..LSMRConfig::default()
        };
        let mut buffers = LSMRBuffers::new(n, n);
        let mut x_damped = vec![0.0; n];
        let mut kernel = LSMRKernel::new(config_damped, &mut buffers);
        let result = kernel.solve(&mut op, &b, &mut x_damped);

        assert!(result.converged);

        // Damped solution should be shrunk toward zero
        // For A=I and damp=λ, the solution is x = b/(1+λ²) = b/2
        for (i, (&x_d, &b_i)) in x_damped.iter().zip(b.iter()).enumerate() {
            let expected = b_i / 2.0;
            assert!(
                (x_d - expected).abs() < 1e-6,
                "damped x[{}] = {} != expected {}",
                i,
                x_d,
                expected
            );
        }

        // Verify undamped solution is different (equals b for identity)
        for (&x_u, &x_d) in x_undamped.iter().zip(x_damped.iter()) {
            assert!(
                (x_u - x_d).abs() > 0.1,
                "Damped and undamped should differ significantly"
            );
        }
    }

    #[test]
    fn test_lsmr_damped_shrinkage() {
        // Test that increasing damp shrinks the solution more
        let mut op = DiagonalOp {
            diag: vec![2.0, 3.0, 4.0],
        };
        let b = vec![2.0, 6.0, 12.0];

        // Collect solutions for different damp values
        let damp_values = [0.0, 0.1, 1.0, 10.0];
        let mut norms: Vec<f64> = Vec::new();

        for &damp in &damp_values {
            let config = LSMRConfig {
                damp,
                ..LSMRConfig::default()
            };
            let mut buffers = LSMRBuffers::new(3, 3);
            let mut x = vec![0.0; 3];
            let mut kernel = LSMRKernel::new(config, &mut buffers);
            kernel.solve(&mut op, &b, &mut x);

            let norm = x.iter().map(|&v| v * v).sum::<f64>().sqrt();
            norms.push(norm);
        }

        // Solution norm should decrease as damp increases
        for i in 1..norms.len() {
            assert!(
                norms[i] <= norms[i - 1] + 1e-10,
                "||x|| should decrease with damp: {:?}",
                norms
            );
        }
    }

    #[test]
    fn test_lsmr_zero_damp_unchanged() {
        // Verify that damp=0 gives same result as before (backward compatibility)
        let mut op = DiagonalOp {
            diag: vec![2.0, 3.0, 4.0],
        };
        let b = vec![2.0, 6.0, 12.0];
        let expected = vec![1.0, 2.0, 3.0]; // Exact solution

        let config = LSMRConfig {
            damp: 0.0,
            ..LSMRConfig::default()
        };
        let mut buffers = LSMRBuffers::new(3, 3);
        let mut x = vec![0.0; 3];
        let mut kernel = LSMRKernel::new(config, &mut buffers);

        let result = kernel.solve(&mut op, &b, &mut x);

        assert!(result.converged);
        for (i, (&x_i, &e_i)) in x.iter().zip(expected.iter()).enumerate() {
            assert!(
                (x_i - e_i).abs() < 1e-8,
                "x[{}] = {} != expected {}",
                i,
                x_i,
                e_i
            );
        }
    }
}
