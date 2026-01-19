//! LSMR algorithm implementation.
//!
//! This module implements the LSMR (Least Squares Minimal Residual) algorithm
//! by Fong and Saunders (2011). LSMR uses Golub-Kahan bidiagonalization to
//! solve least-squares problems without forming the normal equations.
//!
//! This implementation follows the reference Python implementation (pykrylov)
//! by the original authors at Stanford.
//!
//! # Algorithm Overview
//!
//! LSMR solves: min ||A x - b||₂  (or min ||Ax - b||² + λ²||x||² with damping)
//!
//! 1. **Bidiagonalization**: Build orthonormal bases U and V such that
//!    A V = U B where B is lower bidiagonal
//! 2. **QR factorization**: Apply Givens rotations to B (two rotation chains)
//! 3. **Solution update**: Incrementally build x using two direction vectors h, hbar
//!
//! # References
//!
//! - Fong, D. C.-L., & Saunders, M. A. (2011). LSMR: An iterative algorithm
//!   for sparse least-squares problems. SIAM J. Sci. Comput., 33(5), 2950-2971.
//! - Reference implementation: https://github.com/PythonOptimizers/pykrylov

use super::buffers::LSMRBuffers;
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
    /// Absolute tolerance for ||A^T r|| / (||A|| ||r||).
    pub atol: f64,
    /// Relative tolerance for ||r|| / ||b||.
    pub btol: f64,
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
            atol: 1e-8,
            btol: 1e-8,
            maxiter: 10_000,
            conlim: 1e8,
            damp: 0.0,
        }
    }
}

// Keep the old tol field for backward compatibility
impl LSMRConfig {
    /// Create config with a single tolerance (sets both atol and btol).
    pub fn with_tol(tol: f64) -> Self {
        Self {
            atol: tol,
            btol: tol,
            ..Default::default()
        }
    }
}

/// LSMR solver kernel.
///
/// Implements the core LSMR algorithm using Golub-Kahan bidiagonalization.
/// Uses pre-allocated workspace for zero-allocation iteration.
///
/// This implementation follows the reference pykrylov implementation exactly.
pub struct LSMRKernel<'a> {
    config: LSMRConfig,
    buffers: &'a mut LSMRBuffers,
}

impl<'a> LSMRKernel<'a> {
    /// Create a new LSMR kernel with the given configuration and workspace.
    pub fn new(config: LSMRConfig, buffers: &'a mut LSMRBuffers) -> Self {
        Self { config, buffers }
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
    pub fn solve<A: LinearOperator>(
        &mut self,
        operator: &mut A,
        b: &[f64],
        x: &mut [f64],
    ) -> LSMRResult {
        let n_obs = operator.rows();
        let n_coef = operator.cols();

        debug_assert_eq!(b.len(), n_obs);
        debug_assert_eq!(x.len(), n_coef);
        debug_assert_eq!(self.buffers.n_obs(), n_obs);
        debug_assert_eq!(self.buffers.n_coef(), n_coef);

        // Reset workspace for new solve
        self.buffers.reset();

        let damp = self.config.damp;
        let ctol = if self.config.conlim > 0.0 {
            1.0 / self.config.conlim
        } else {
            0.0
        };

        // =====================================================================
        // Initialize the Golub-Kahan bidiagonalization process
        // =====================================================================

        // u = b, beta = ||u||
        self.buffers.u.copy_from_slice(b);
        self.buffers.beta = norm(&self.buffers.u);
        self.buffers.normb = self.buffers.beta;

        if self.buffers.beta > 0.0 {
            // u = u / beta
            scale(&mut self.buffers.u, 1.0 / self.buffers.beta);

            // v = A^T * u
            operator.rmatvec(&self.buffers.u, &mut self.buffers.v);
            self.buffers.alpha = norm(&self.buffers.v);
        }

        if self.buffers.alpha > 0.0 {
            // v = v / alpha
            scale(&mut self.buffers.v, 1.0 / self.buffers.alpha);
        }

        // =====================================================================
        // Initialize variables for 1st iteration
        // =====================================================================

        self.buffers.zetabar = self.buffers.alpha * self.buffers.beta;
        self.buffers.alphabar = self.buffers.alpha;
        self.buffers.rho = 1.0;
        self.buffers.rhobar = 1.0;
        self.buffers.cbar = 1.0;
        self.buffers.sbar = 0.0;

        // h = v.copy()
        self.buffers.h.copy_from_slice(&self.buffers.v);
        // hbar = zeros (already zero from reset)

        // =====================================================================
        // Initialize variables for estimation of ||r||
        // =====================================================================

        self.buffers.betadd = self.buffers.beta;
        self.buffers.betad = 0.0;
        self.buffers.rhodold = 1.0;
        self.buffers.tautildeold = 0.0;
        self.buffers.thetatilde = 0.0;
        self.buffers.zeta = 0.0;
        self.buffers.d = 0.0;

        // =====================================================================
        // Initialize variables for estimation of ||A|| and cond(A)
        // =====================================================================

        self.buffers.normA2 = self.buffers.alpha * self.buffers.alpha;
        self.buffers.maxrbar = 0.0;
        self.buffers.minrbar = 1e100;
        self.buffers.normA = self.buffers.normA2.sqrt();
        self.buffers.condA = 1.0;

        // =====================================================================
        // Items for use in stopping rules
        // =====================================================================

        self.buffers.normr = self.buffers.beta;
        self.buffers.normar = self.buffers.alpha * self.buffers.beta;

        // Check for trivial case
        if self.buffers.beta == 0.0 {
            x.fill(0.0);
            return LSMRResult {
                iterations: 0,
                converged: true,
                residual_norm: 0.0,
                atr_norm: 0.0,
            };
        }

        // Exit if A'b = 0 (x = 0 is the solution)
        if self.buffers.normar == 0.0 {
            x.fill(0.0);
            return LSMRResult {
                iterations: 0,
                converged: true,
                residual_norm: self.buffers.normr,
                atr_norm: 0.0,
            };
        }

        // =====================================================================
        // Main iteration loop
        // =====================================================================

        let mut istop = 0;
        let mut final_itn = 0;

        for itn in 1..=self.config.maxiter {
            final_itn = itn;
            // =================================================================
            // Perform the next step of the bidiagonalization
            // =================================================================

            // u = A * v - alpha * u
            operator.matvec(&self.buffers.v, &mut self.buffers.matvec_scratch);
            let alpha_old = self.buffers.alpha;
            for (u_i, &scratch_i) in self
                .buffers
                .u
                .iter_mut()
                .zip(self.buffers.matvec_scratch.iter())
            {
                *u_i = scratch_i - alpha_old * *u_i;
            }
            self.buffers.beta = norm(&self.buffers.u);

            if self.buffers.beta > 0.0 {
                scale(&mut self.buffers.u, 1.0 / self.buffers.beta);

                // v = A^T * u - beta * v
                operator.rmatvec(&self.buffers.u, &mut self.buffers.precond_scratch);
                let beta = self.buffers.beta;
                for (v_i, &scratch_i) in self
                    .buffers
                    .v
                    .iter_mut()
                    .zip(self.buffers.precond_scratch.iter())
                {
                    *v_i = scratch_i - beta * *v_i;
                }
                self.buffers.alpha = norm(&self.buffers.v);

                if self.buffers.alpha > 0.0 {
                    scale(&mut self.buffers.v, 1.0 / self.buffers.alpha);
                }
            }

            // At this point, beta = beta_{k+1}, alpha = alpha_{k+1}

            // =================================================================
            // Construct rotation Qhat_{k,2k+1} (for damping)
            // =================================================================

            let (chat, shat, alphahat) = sym_ortho(self.buffers.alphabar, damp);

            // =================================================================
            // Use a plane rotation (Q_i) to turn B_i to R_i
            // =================================================================

            self.buffers.rhoold = self.buffers.rho;
            let (c, s, rho) = sym_ortho(alphahat, self.buffers.beta);
            self.buffers.c = c;
            self.buffers.s = s;
            self.buffers.rho = rho;
            let thetanew = s * self.buffers.alpha;
            self.buffers.alphabar = c * self.buffers.alpha;

            // =================================================================
            // Use a plane rotation (Qbar_i) to turn R_i^T to R_i^bar
            // =================================================================

            self.buffers.rhobarold = self.buffers.rhobar;
            let zetaold = self.buffers.zeta;
            let thetabar = self.buffers.sbar * self.buffers.rho;
            let rhotemp = self.buffers.cbar * self.buffers.rho;
            let (cbar, sbar, rhobar) = sym_ortho(self.buffers.cbar * self.buffers.rho, thetanew);
            self.buffers.cbar = cbar;
            self.buffers.sbar = sbar;
            self.buffers.rhobar = rhobar;
            self.buffers.zeta = self.buffers.cbar * self.buffers.zetabar;
            self.buffers.zetabar = -self.buffers.sbar * self.buffers.zetabar;

            // =================================================================
            // Update h, hbar, x
            // =================================================================

            // hbar = h - (thetabar * rho / (rhoold * rhobarold)) * hbar
            let hbar_factor = thetabar * self.buffers.rho
                / (self.buffers.rhoold * self.buffers.rhobarold);
            for (hbar_i, &h_i) in self.buffers.hbar.iter_mut().zip(self.buffers.h.iter()) {
                *hbar_i = h_i - hbar_factor * *hbar_i;
            }

            // x = x + (zeta / (rho * rhobar)) * hbar
            let x_factor = self.buffers.zeta / (self.buffers.rho * self.buffers.rhobar);
            for (x_i, &hbar_i) in self.buffers.x.iter_mut().zip(self.buffers.hbar.iter()) {
                *x_i += x_factor * hbar_i;
            }

            // h = v - (thetanew / rho) * h
            let h_factor = thetanew / self.buffers.rho;
            for (h_i, &v_i) in self.buffers.h.iter_mut().zip(self.buffers.v.iter()) {
                *h_i = v_i - h_factor * *h_i;
            }

            // =================================================================
            // Estimate of ||r||
            // =================================================================

            // Apply rotation Qhat_{k,2k+1}
            let betaacute = chat * self.buffers.betadd;
            let betacheck = -shat * self.buffers.betadd;

            // Apply rotation Q_{k,k+1}
            let betahat = self.buffers.c * betaacute;
            self.buffers.betadd = -self.buffers.s * betaacute;

            // Apply rotation Qtilde_{k-1}
            let thetatildeold = self.buffers.thetatilde;
            let (ctildeold, stildeold, rhotildeold) =
                sym_ortho(self.buffers.rhodold, thetabar);
            self.buffers.thetatilde = stildeold * self.buffers.rhobar;
            self.buffers.rhodold = ctildeold * self.buffers.rhobar;
            self.buffers.betad = -stildeold * self.buffers.betad + ctildeold * betahat;

            self.buffers.tautildeold =
                (zetaold - thetatildeold * self.buffers.tautildeold) / rhotildeold;
            let taud =
                (self.buffers.zeta - self.buffers.thetatilde * self.buffers.tautildeold)
                    / self.buffers.rhodold;
            self.buffers.d += betacheck * betacheck;
            self.buffers.normr = (self.buffers.d
                + (self.buffers.betad - taud).powi(2)
                + self.buffers.betadd.powi(2))
            .sqrt();

            // =================================================================
            // Estimate ||A||
            // =================================================================

            self.buffers.normA2 += self.buffers.beta * self.buffers.beta;
            self.buffers.normA = self.buffers.normA2.sqrt();
            self.buffers.normA2 += self.buffers.alpha * self.buffers.alpha;

            // =================================================================
            // Estimate cond(A)
            // =================================================================

            self.buffers.maxrbar = self.buffers.maxrbar.max(self.buffers.rhobarold);
            if itn > 1 {
                self.buffers.minrbar = self.buffers.minrbar.min(self.buffers.rhobarold);
            }
            self.buffers.condA = self.buffers.maxrbar.max(rhotemp)
                / self.buffers.minrbar.min(rhotemp);

            // =================================================================
            // Test for convergence
            // =================================================================

            // Compute norms for convergence testing
            self.buffers.normar = self.buffers.zetabar.abs();
            let normx = norm(&self.buffers.x);

            // Now use these norms to estimate certain other quantities
            let test1 = self.buffers.normr / self.buffers.normb;
            let test2 = self.buffers.normar / (self.buffers.normA * self.buffers.normr);
            let test3 = 1.0 / self.buffers.condA;
            let t1 = test1 / (1.0 + self.buffers.normA * normx / self.buffers.normb);
            let rtol =
                self.config.btol + self.config.atol * self.buffers.normA * normx / self.buffers.normb;

            // The following tests guard against extremely small values of
            // atol, btol or ctol. (The user may have set any or all of
            // the parameters atol, btol, conlim to 0.)
            if itn >= self.config.maxiter {
                istop = 7;
            }
            if 1.0 + test3 <= 1.0 {
                istop = 6;
            }
            if 1.0 + test2 <= 1.0 {
                istop = 5;
            }
            if 1.0 + t1 <= 1.0 {
                istop = 4;
            }

            // Allow for tolerances set by the user
            if test3 <= ctol {
                istop = 3;
            }
            if test2 <= self.config.atol {
                istop = 2;
            }
            if test1 <= rtol {
                istop = 1;
            }

            if istop > 0 {
                break;
            }
        }

        // Copy solution to output
        x.copy_from_slice(&self.buffers.x);

        let converged = istop != 3 && istop != 6 && istop != 7;

        LSMRResult {
            iterations: final_itn,
            converged,
            residual_norm: self.buffers.normr,
            atr_norm: self.buffers.normar,
        }
    }
}

// =============================================================================
// Utility Functions
// =============================================================================

/// Numerically stable Givens rotation (symOrtho).
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
fn norm(x: &[f64]) -> f64 {
    x.iter().map(|&v| v * v).sum::<f64>().sqrt()
}

/// Scale vector in place: x = alpha * x.
#[inline]
fn scale(x: &mut [f64], alpha: f64) {
    for xi in x.iter_mut() {
        *xi *= alpha;
    }
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

    /// Rectangular matrix operator for overdetermined least-squares testing
    struct RectangularOp {
        rows: usize,
        cols: usize,
        // Store as row-major: data[i * cols + j] = A[i, j]
        data: Vec<f64>,
    }

    impl RectangularOp {
        fn new(rows: usize, cols: usize, data: Vec<f64>) -> Self {
            assert_eq!(data.len(), rows * cols);
            Self { rows, cols, data }
        }
    }

    impl LinearOperator for RectangularOp {
        fn rows(&self) -> usize {
            self.rows
        }
        fn cols(&self) -> usize {
            self.cols
        }
        fn matvec(&mut self, x: &[f64], y: &mut [f64]) {
            // y = A * x
            for i in 0..self.rows {
                y[i] = 0.0;
                for j in 0..self.cols {
                    y[i] += self.data[i * self.cols + j] * x[j];
                }
            }
        }
        fn rmatvec(&mut self, y: &[f64], x: &mut [f64]) {
            // x = A^T * y
            for j in 0..self.cols {
                x[j] = 0.0;
                for i in 0..self.rows {
                    x[j] += self.data[i * self.cols + j] * y[i];
                }
            }
        }
    }

    #[test]
    fn test_lsmr_overdetermined() {
        // Test overdetermined system (more rows than columns)
        // A = [[1, 0], [0, 1], [1, 1]], b = [1, 2, 4]
        // Least squares solution minimizes ||Ax - b||²
        let mut op = RectangularOp::new(
            3,
            2,
            vec![
                1.0, 0.0, // row 0
                0.0, 1.0, // row 1
                1.0, 1.0, // row 2
            ],
        );
        let b = vec![1.0, 2.0, 4.0];
        let mut x = vec![0.0; 2];

        let config = LSMRConfig::default();
        let mut buffers = LSMRBuffers::new(3, 2);
        let mut kernel = LSMRKernel::new(config, &mut buffers);

        let result = kernel.solve(&mut op, &b, &mut x);

        assert!(result.converged);

        // Verify solution satisfies normal equations: A^T A x = A^T b
        // A^T A = [[2, 1], [1, 2]], A^T b = [5, 6]
        // Solution: x = [4/3, 7/3] ≈ [1.333, 2.333]
        let expected = [4.0 / 3.0, 7.0 / 3.0];
        for (i, (&x_i, &e_i)) in x.iter().zip(expected.iter()).enumerate() {
            assert!(
                (x_i - e_i).abs() < 1e-6,
                "x[{}] = {} != expected {}",
                i,
                x_i,
                e_i
            );
        }
    }
}
