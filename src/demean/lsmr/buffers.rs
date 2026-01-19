//! Pre-allocated workspace for LSMR iteration.
//!
//! This module provides `LSMRBuffers` which holds all working vectors and scalar
//! state needed by the LSMR algorithm. Workspace is allocated once and reused
//! across iterations and multiple solve() calls to avoid allocation overhead.

/// Pre-allocated workspace for LSMR iteration.
///
/// Contains both vector buffers and scalar state for the algorithm.
/// Follows the reference implementation by Fong & Saunders (pykrylov).
///
/// **Vectors** (Golub-Kahan bidiagonalization):
/// - `u`: Bidiagonalization vector in observation space
/// - `v`: Bidiagonalization vector in coefficient space
/// - `h`, `hbar`: Update direction vectors (LSMR uses two)
/// - `x`: Current solution estimate
///
/// **Scalars** (iteration state):
/// - Bidiagonalization: `alpha`, `beta`
/// - QR factorization: `alphabar`, `rho`, `rhobar`, `c`, `s`
/// - QR-bar factorization: `cbar`, `sbar`, `zetabar`, `zeta`
/// - Residual estimation: `betadd`, `betad`, `rhodold`, `tautildeold`, `thetatilde`, `d`
/// - Convergence monitoring: `normA2`, `normA`, `condA`, `normr`, `normar`
///
/// All state is reset via `reset()` between solves.
#[derive(Clone)]
pub struct LSMRBuffers {
    // Bidiagonalization vectors
    /// u vector in observation space (length: n_obs)
    pub(crate) u: Vec<f64>,
    /// v vector in coefficient space (length: n_coef)
    pub(crate) v: Vec<f64>,

    // Solution update vectors (LSMR uses two direction vectors)
    /// h vector for solution updates (length: n_coef)
    pub(crate) h: Vec<f64>,
    /// hbar vector for solution updates (length: n_coef)
    pub(crate) hbar: Vec<f64>,
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

    // QR factorization state (first Givens rotation chain)
    /// alphabar for QR factorization
    pub(crate) alphabar: f64,
    /// rho from current iteration
    pub(crate) rho: f64,
    /// rho from previous iteration
    pub(crate) rhoold: f64,
    /// rhobar from QR-bar factorization
    pub(crate) rhobar: f64,
    /// rhobar from previous iteration
    pub(crate) rhobarold: f64,
    /// c (cosine) from Q rotation
    pub(crate) c: f64,
    /// s (sine) from Q rotation
    pub(crate) s: f64,

    // QR-bar factorization state (second Givens rotation chain)
    /// cbar (cosine) from Qbar rotation
    pub(crate) cbar: f64,
    /// sbar (sine) from Qbar rotation
    pub(crate) sbar: f64,
    /// zetabar for solution update
    pub(crate) zetabar: f64,
    /// zeta for solution update
    pub(crate) zeta: f64,

    // Residual norm estimation (||r|| without computing r explicitly)
    /// betadd for residual estimation
    pub(crate) betadd: f64,
    /// betad for residual estimation
    pub(crate) betad: f64,
    /// rhodold for residual estimation
    pub(crate) rhodold: f64,
    /// tautildeold for residual estimation
    pub(crate) tautildeold: f64,
    /// thetatilde for residual estimation
    pub(crate) thetatilde: f64,
    /// d (accumulated term) for residual estimation
    pub(crate) d: f64,

    // Convergence monitoring
    /// Squared Frobenius norm estimate ||A||_F^2
    pub(crate) normA2: f64,
    /// Frobenius norm estimate ||A||_F
    pub(crate) normA: f64,
    /// Condition number estimate cond(A)
    pub(crate) condA: f64,
    /// maxrbar for condition estimation
    pub(crate) maxrbar: f64,
    /// minrbar for condition estimation
    pub(crate) minrbar: f64,
    /// ||r|| = ||b - Ax||
    pub(crate) normr: f64,
    /// ||A^T r|| = |zetabar|
    pub(crate) normar: f64,
    /// ||b|| for relative tolerance
    pub(crate) normb: f64,
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
            h: vec![0.0; n_coef],
            hbar: vec![0.0; n_coef],
            x: vec![0.0; n_coef],
            precond_scratch: vec![0.0; n_coef],
            matvec_scratch: vec![0.0; n_obs],
            // Scalars initialized to match pykrylov
            alpha: 0.0,
            beta: 0.0,
            alphabar: 0.0,
            rho: 1.0,
            rhoold: 1.0,
            rhobar: 1.0,
            rhobarold: 1.0,
            c: 1.0,
            s: 0.0,
            cbar: 1.0,
            sbar: 0.0,
            zetabar: 0.0,
            zeta: 0.0,
            betadd: 0.0,
            betad: 0.0,
            rhodold: 1.0,
            tautildeold: 0.0,
            thetatilde: 0.0,
            d: 0.0,
            normA2: 0.0,
            normA: 0.0,
            condA: 1.0,
            maxrbar: 0.0,
            minrbar: 1e100,
            normr: 0.0,
            normar: 0.0,
            normb: 0.0,
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
        self.h.fill(0.0);
        self.hbar.fill(0.0);
        self.x.fill(0.0);
        // precond_scratch and matvec_scratch don't need reset as they're
        // always overwritten before use

        // Reset scalars to initial values (matching pykrylov)
        self.alpha = 0.0;
        self.beta = 0.0;
        self.alphabar = 0.0;
        self.rho = 1.0;
        self.rhoold = 1.0;
        self.rhobar = 1.0;
        self.rhobarold = 1.0;
        self.c = 1.0;
        self.s = 0.0;
        self.cbar = 1.0;
        self.sbar = 0.0;
        self.zetabar = 0.0;
        self.zeta = 0.0;
        self.betadd = 0.0;
        self.betad = 0.0;
        self.rhodold = 1.0;
        self.tautildeold = 0.0;
        self.thetatilde = 0.0;
        self.d = 0.0;
        self.normA2 = 0.0;
        self.normA = 0.0;
        self.condA = 1.0;
        self.maxrbar = 0.0;
        self.minrbar = 1e100;
        self.normr = 0.0;
        self.normar = 0.0;
        self.normb = 0.0;
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
        buffers.rho = 5.0;
        buffers.cbar = 0.5;

        // Reset and verify
        buffers.reset();
        assert_eq!(buffers.x[0], 0.0);
        assert_eq!(buffers.u[0], 0.0);
        assert_eq!(buffers.rho, 1.0);
        assert_eq!(buffers.cbar, 1.0);
    }

    #[test]
    fn test_buffers_initial_values() {
        let buffers = LSMRBuffers::new(10, 5);
        // Check initial values match pykrylov
        assert_eq!(buffers.rho, 1.0);
        assert_eq!(buffers.rhobar, 1.0);
        assert_eq!(buffers.cbar, 1.0);
        assert_eq!(buffers.sbar, 0.0);
        assert_eq!(buffers.rhodold, 1.0);
        assert_eq!(buffers.minrbar, 1e100);
    }
}
