//! Pre-allocated buffers for LSMR iteration.
//!
//! This module provides `LSMRBuffers` which holds all working vectors needed
//! by the LSMR algorithm. Buffers are allocated once and reused across
//! iterations and multiple solve() calls to avoid allocation overhead.

/// Pre-allocated buffers for LSMR iteration.
///
/// LSMR uses Golub-Kahan bidiagonalization which maintains several vectors:
/// - `u`: Bidiagonalization vector in observation space
/// - `v`: Bidiagonalization vector in coefficient space
/// - `w`: Update direction vector
/// - `x`: Current solution estimate
///
/// All vectors are allocated at construction and reused via `reset()`.
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
}

impl LSMRBuffers {
    /// Create new LSMR buffers for the given problem dimensions.
    ///
    /// # Arguments
    /// * `n_obs` - Number of observations (row space dimension)
    /// * `n_coef` - Number of coefficients (column space dimension)
    pub fn new(n_obs: usize, n_coef: usize) -> Self {
        Self {
            u: vec![0.0; n_obs],
            v: vec![0.0; n_coef],
            w: vec![0.0; n_coef],
            x: vec![0.0; n_coef],
            precond_scratch: vec![0.0; n_coef],
            matvec_scratch: vec![0.0; n_obs],
        }
    }

    /// Reset all buffers to zero for a new solve.
    ///
    /// This is more efficient than reallocating for repeated solves
    /// on problems with the same dimensions.
    pub fn reset(&mut self) {
        self.u.fill(0.0);
        self.v.fill(0.0);
        self.w.fill(0.0);
        self.x.fill(0.0);
        // precond_scratch and matvec_scratch don't need reset as they're
        // always overwritten before use
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

/// Scalar state for LSMR iteration.
///
/// These scalars track the bidiagonalization and QR factorization state
/// across iterations. Separated from buffers for clarity.
#[derive(Clone, Copy, Default)]
pub struct LSMRState {
    // Bidiagonalization scalars
    /// alpha_k from bidiagonalization
    pub(crate) alpha: f64,
    /// beta_k from bidiagonalization
    pub(crate) beta: f64,

    // QR factorization state (Givens rotations)
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

impl LSMRState {
    /// Create a new LSMR state with default values.
    pub fn new() -> Self {
        Self::default()
    }

    /// Reset state for a new solve.
    pub fn reset(&mut self) {
        *self = Self::default();
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
