//! LSMR-based fixed effects demeaning.
//!
//! This module provides an alternative solver for fixed effects demeaning
//! using the LSMR (Least Squares Minimal Residual) algorithm instead of
//! the default accelerated Gauss-Seidel approach.
//!
//! # When to Use LSMR
//!
//! LSMR may converge faster than Gauss-Seidel for:
//! - Problems with many fixed effects (3+)
//! - Unbalanced group sizes
//! - Highly correlated fixed effects
//!
//! # Usage
//!
//! LSMR is selected via the `solver` parameter in the demean API:
//! ```python
//! demean(x, flist, weights, solver="lsmr")
//! ```

pub mod buffers;
pub mod kernel;
pub mod linear_operator;
pub mod preconditioner;

use crate::demean::demeaner::Demeaner;
use crate::demean::types::{ConvergenceState, DemeanContext, DemeanResult, LSMRConfig};
use buffers::LSMRBuffers;
use kernel::{LSMRConfig as KernelConfig, LSMRKernel};
use linear_operator::{DesignMatrixOperator, PreconditionedOperator};
use preconditioner::{create_preconditioner, BoxedPreconditioner};

/// LSMR-based demeaner.
///
/// Uses LSMR to solve the least-squares problem for FE coefficients:
/// ```text
/// min_coef ||input - D * coef||₂
/// ```
/// where D is the fixed-effects design matrix.
///
/// The demeaned output is then: `output = input - D * coef`.
pub struct LSMRDemeaner<'a> {
    ctx: &'a DemeanContext,
    config: LSMRConfig,
    buffers: LSMRBuffers,
    preconditioner: BoxedPreconditioner<'a>,
}

impl<'a> LSMRDemeaner<'a> {
    /// Create a new LSMR demeaner.
    ///
    /// # Arguments
    /// * `ctx` - The demeaning context with FE information
    /// * `config` - LSMR configuration parameters
    pub fn new(ctx: &'a DemeanContext, config: LSMRConfig) -> Self {
        let buffers = LSMRBuffers::new(ctx.dims.n_obs, ctx.dims.n_coef);
        let preconditioner = create_preconditioner(ctx, config.preconditioner);

        Self {
            ctx,
            config,
            buffers,
            preconditioner,
        }
    }
}

impl Demeaner for LSMRDemeaner<'_> {
    fn solve(&mut self, input: &[f64]) -> DemeanResult {
        let n_obs = self.ctx.dims.n_obs;
        let n_coef = self.ctx.dims.n_coef;

        debug_assert_eq!(input.len(), n_obs);

        // Create the linear operator (handles sqrt(weights) scaling internally)
        let mut base_operator = DesignMatrixOperator::new(self.ctx);

        // Create kernel config from our config
        let kernel_config = KernelConfig {
            tol: self.config.tol,
            maxiter: self.config.maxiter,
            conlim: 1e8,
        };

        // For weighted least squares, we solve min ||W^{1/2}(b - Dx)||
        // which transforms to min ||b̃ - Ãx|| where:
        // - Ã = W^{1/2} D (handled by DesignMatrixOperator)
        // - b̃ = W^{1/2} b (we scale input here)
        let scaled_input: Vec<f64> = match &self.ctx.weights {
            Some(w) => input
                .iter()
                .zip(w.iter())
                .map(|(&inp, &wi)| inp * wi.sqrt())
                .collect(),
            None => input.to_vec(),
        };

        // Solve using right-preconditioned LSMR:
        // min ||Ã * M^{-1} * z - b̃||, then recover coef = M^{-1} * z
        let mut precond_operator =
            PreconditionedOperator::new(&mut base_operator, self.preconditioner.as_mut());

        let mut kernel = LSMRKernel::new(kernel_config, &mut self.buffers);
        let mut z = vec![0.0; n_coef]; // Preconditioned variable
        let result = kernel.solve(&mut precond_operator, &scaled_input, &mut z);

        // Recover actual coefficients: coef = M^{-1} * z
        let mut coef = vec![0.0; n_coef];
        self.preconditioner.apply(&z, &mut coef);

        // Compute demeaned output: demeaned = input - D * coef
        // Note: We use the original input here, not scaled_input
        let mut demeaned = input.to_vec();
        for fe in &self.ctx.fe_infos {
            let offset = fe.coef_start;
            for (i, &g) in fe.group_ids.iter().enumerate() {
                demeaned[i] -= coef[offset + g];
            }
        }

        // Reorder coefficients to original FE order
        let fe_coefficients = self.ctx.reorder_coef_to_original(&coef);

        let convergence = if result.converged {
            ConvergenceState::Converged
        } else {
            ConvergenceState::NotConverged
        };

        DemeanResult {
            demeaned,
            fe_coefficients,
            convergence,
            iterations: result.iterations,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;

    #[test]
    fn test_lsmr_demeaner_single_fe() {
        // Simple test with 1 FE: 6 observations, 2 groups
        // Group 0: obs 0, 1, 2 with values 10, 20, 30 → mean = 20
        // Group 1: obs 3, 4, 5 with values 40, 50, 60 → mean = 50
        // Demeaned: -10, 0, 10, -10, 0, 10

        let flist = Array2::from_shape_vec((6, 1), vec![0, 0, 0, 1, 1, 1]).unwrap();
        let flist_view = flist.view();

        let ctx = DemeanContext::new(&flist_view, None, false);
        let config = LSMRConfig::default();
        let mut demeaner = LSMRDemeaner::new(&ctx, config);

        let input = vec![10.0, 20.0, 30.0, 40.0, 50.0, 60.0];

        let result = demeaner.solve(&input);

        assert_eq!(result.convergence, ConvergenceState::Converged);

        // Check demeaned values are approximately correct
        let expected = vec![-10.0, 0.0, 10.0, -10.0, 0.0, 10.0];
        for (i, (&o, &e)) in result.demeaned.iter().zip(expected.iter()).enumerate() {
            assert!(
                (o - e).abs() < 1e-6,
                "output[{}] = {} != expected {}",
                i,
                o,
                e
            );
        }
    }

    #[test]
    fn test_lsmr_demeaner_two_fe() {
        // 2 FEs with unbalanced groups - diagonal preconditioning should help
        // 100 observations with 10 groups in FE1 and 5 groups in FE2
        use crate::demean::lsmr::preconditioner::PreconditionerKind;

        let n = 100;
        let n_groups_1 = 10;
        let n_groups_2 = 5;

        // Create FE assignments
        let mut flist_data = Vec::with_capacity(n * 2);
        for i in 0..n {
            flist_data.push(i % n_groups_1); // FE1: cyclic 0-9
            flist_data.push((i / 20) % n_groups_2); // FE2: blocks of 20
        }
        let flist = Array2::from_shape_vec((n, 2), flist_data).unwrap();
        let flist_view = flist.view();

        // Create input data
        let input: Vec<f64> = (0..n).map(|i| i as f64).collect();

        let ctx = DemeanContext::new(&flist_view, None, false);

        // Test with identity (no preconditioning)
        let config_identity = LSMRConfig {
            preconditioner: PreconditionerKind::None,
            ..LSMRConfig::default()
        };
        let mut demeaner_identity = LSMRDemeaner::new(&ctx, config_identity);
        let result_identity = demeaner_identity.solve(&input);

        // Test with diagonal preconditioning
        let config_diagonal = LSMRConfig {
            preconditioner: PreconditionerKind::Diagonal,
            ..LSMRConfig::default()
        };
        let mut demeaner_diagonal = LSMRDemeaner::new(&ctx, config_diagonal);
        let result_diagonal = demeaner_diagonal.solve(&input);

        // Both should converge
        assert_eq!(result_identity.convergence, ConvergenceState::Converged);
        assert_eq!(result_diagonal.convergence, ConvergenceState::Converged);

        // Both should produce the same demeaned output
        for (i, (&a, &b)) in result_identity
            .demeaned
            .iter()
            .zip(result_diagonal.demeaned.iter())
            .enumerate()
        {
            assert!(
                (a - b).abs() < 1e-6,
                "Mismatch at {}: identity={} diagonal={}",
                i,
                a,
                b
            );
        }

        println!(
            "Identity iterations: {}, Diagonal iterations: {}",
            result_identity.iterations, result_diagonal.iterations
        );

        // Diagonal preconditioning should use fewer or equal iterations
        // (or at least not significantly more)
        assert!(
            result_diagonal.iterations <= result_identity.iterations + 5,
            "Diagonal ({}) used many more iterations than identity ({})",
            result_diagonal.iterations,
            result_identity.iterations
        );
    }
}
