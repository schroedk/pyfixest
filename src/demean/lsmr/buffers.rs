//! Pre-allocated workspace for LSMR iteration.
//!
//! This module provides `LSMRBuffers` which holds all working vectors and scalar
//! state needed by the LSMR algorithm. Workspace is allocated once and reused
//! across iterations and multiple solve() calls to avoid allocation overhead.

/// Pre-allocated workspace for LSMR iteration.
///
/// Contains both vector buffers and scalar state for the algorithm:
///
/// **Vectors** (Golub-Kahan bidiagonalization):
/// - `u`: Bidiagonalization vector in observation space
/// - `v`: Bidiagonalization vector in coefficient space
/// - `w`: Update direction vector
/// - `x`: Current solution estimate
///
/// **Scalars** (iteration state):
/// - Bidiagonalization: `alpha`, `beta`
/// - QR factorization: `alphabar`, `rho_bar`, `phi_bar`, `c`, `s`
/// - Convergence monitoring: `anorm`, `acond`, `rnorm`, `arnorm`
///
/// All state is reset via `reset()` between solves.
#[derive(Clone)]
pub struct LSMRBuffers {
    // Bidiagonalization vectors
    /// u vector in observation space (length: n_obs)
    pub(crate) u: Vec<f64>,
    /// v vector in coefficient space (length: n_coef)
    pub(crate) v: Vec<f64>,

    // Solution update vectors
    /// w vector for solution updates (length: n_coef)
    pub(crate) w: Vec<f64>,
    /// Current solution estimate (length: n_coef)
    pub(crate) x: Vec<f64>,

    // Scratch buffers
    /// Scratch for preconditioner application (length: n_coef)
    pub(crate) precond_scratch: Vec<f64>,
    /// Scratch for matvec operations (length: n_obs)
    pub(crate) matvec_scratch: Vec<f64>,

    // Bidiagonalization scalars
    /// alpha_k from bidiagonalization
    pub(crate) alpha: f64,
    /// beta_k from bidiagonalization
    pub(crate) beta: f64,

    // QR factorization state (Givens rotations)
    /// alpha_bar for damped LSMR (carries across iterations)
    pub(crate) alphabar: f64,
    /// rho_bar from QR
    pub(crate) rho_bar: f64,
    /// phi_bar from QR
    pub(crate) phi_bar: f64,
    /// c (cosine) from previous Givens rotation
    pub(crate) c: f64,
    /// s (sine) from previous Givens rotation
    pub(crate) s: f64,

    // Convergence monitoring
    /// Estimate of ||A||_F
    pub(crate) anorm: f64,
    /// Estimate of cond(A)
    pub(crate) acond: f64,
    /// ||r|| = ||b - Ax||
    pub(crate) rnorm: f64,
    /// ||A^T r||
    pub(crate) arnorm: f64,
}

impl LSMRBuffers {
    /// Create new LSMR workspace for the given problem dimensions.
    ///
    /// # Arguments
    /// * `n_obs` - Number of observations (row space dimension)
    /// * `n_coef` - Number of coefficients (column space dimension)
    pub fn new(n_obs: usize, n_coef: usize) -> Self {
        Self {
            // Vectors
            u: vec![0.0; n_obs],
            v: vec![0.0; n_coef],
            w: vec![0.0; n_coef],
            x: vec![0.0; n_coef],
            precond_scratch: vec![0.0; n_coef],
            matvec_scratch: vec![0.0; n_obs],
            // Scalars (all zero-initialized)
            alpha: 0.0,
            beta: 0.0,
            alphabar: 0.0,
            rho_bar: 0.0,
            phi_bar: 0.0,
            c: 0.0,
            s: 0.0,
            anorm: 0.0,
            acond: 0.0,
            rnorm: 0.0,
            arnorm: 0.0,
        }
    }

    /// Reset all workspace state for a new solve.
    ///
    /// This is more efficient than reallocating for repeated solves
    /// on problems with the same dimensions.
    pub fn reset(&mut self) {
        // Reset vectors
        self.u.fill(0.0);
        self.v.fill(0.0);
        self.w.fill(0.0);
        self.x.fill(0.0);
        // precond_scratch and matvec_scratch don't need reset as they're
        // always overwritten before use

        // Reset scalars
        self.alpha = 0.0;
        self.beta = 0.0;
        self.alphabar = 0.0;
        self.rho_bar = 0.0;
        self.phi_bar = 0.0;
        self.c = 0.0;
        self.s = 0.0;
        self.anorm = 0.0;
        self.acond = 0.0;
        self.rnorm = 0.0;
        self.arnorm = 0.0;
    }

    /// Get the observation space dimension.
    #[inline]
    pub fn n_obs(&self) -> usize {
        self.u.len()
    }

    /// Get the coefficient space dimension.
    #[inline]
    pub fn n_coef(&self) -> usize {
        self.v.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_buffers_dimensions() {
        let buffers = LSMRBuffers::new(1000, 50);
        assert_eq!(buffers.n_obs(), 1000);
        assert_eq!(buffers.n_coef(), 50);
    }

    #[test]
    fn test_buffers_reset() {
        let mut buffers = LSMRBuffers::new(10, 5);

        // Modify some values
        buffers.x[0] = 1.0;
        buffers.u[0] = 2.0;

        // Reset and verify
        buffers.reset();
        assert_eq!(buffers.x[0], 0.0);
        assert_eq!(buffers.u[0], 0.0);
    }
}
