//! Fixed-iteration Preconditioned Conjugate Gradient for Gram systems.
//!
//! Solves `G * x = r` where G is the bipartite Gram matrix (Deg + Adj),
//! using the tree Gram `G_T` as preconditioner.
//!
//! Uses exactly `k_inner` iterations (no early stopping) to maintain linearity,
//! which is required for use as a preconditioner in LSMR.
//!
//! # Gram vs Laplacian
//!
//! For a bipartite graph with P (left) and Q (right) nodes:
//! - Laplacian L = Deg - Adj (negative off-diagonal)
//! - Gram G = Deg + Adj (positive off-diagonal, used for normal equations)
//!
//! The Gram matrix is the correct operator for preconditioning LSMR because
//! it matches the structure of the FE design matrix's normal equations A^T A.
//!
//! # Sign Transformation
//!
//! The tree Laplacian solve L_T^{-1} can be transformed to a tree Gram solve:
//! G_T^{-1} r = S L_T^{-1} (S r)
//!
//! where S = diag(+1 for P nodes, -1 for Q nodes) flips signs of Q nodes.

use super::forest_builder::SpanningForest;
use super::solver::TreeSolverOwned;
use crate::graph::components::ComponentsCSR;
use crate::laplacian::operator::LaplacianOp;

/// Fixed-iteration Preconditioned Conjugate Gradient solver.
///
/// Solves `G * x ≈ r` where `G` is the bipartite Gram matrix and uses
/// `G_T^{-1}` (tree Gram solve via sign transformation) as preconditioner.
///
/// # Linearity Requirement
/// This solver uses exactly `k_inner` iterations regardless of convergence.
/// This ensures the mapping `r -> x` is linear, which is required for
/// use as a preconditioner in LSMR (which assumes M^{-1} is linear).
///
/// # Algorithm
/// Standard PCG with Gram matrix:
/// 1. Project RHS to Range(G)
/// 2. r₀ = rhs (since x₀ = 0)
/// 3. z₀ = G_T^{-1} r₀ (via sign transformation)
/// 4. p₀ = z₀
/// 5. For k = 1..k_inner:
///    - αₖ = (rₖ, zₖ) / (pₖ, G pₖ)
///    - xₖ₊₁ = xₖ + αₖ pₖ
///    - rₖ₊₁ = rₖ - αₖ G pₖ
///    - zₖ₊₁ = G_T^{-1} rₖ₊₁
///    - βₖ = (rₖ₊₁, zₖ₊₁) / (rₖ, zₖ)
///    - pₖ₊₁ = zₖ₊₁ + βₖ pₖ
pub struct FixedPCG {
    /// Number of inner iterations (fixed)
    k_inner: usize,
    /// Number of P nodes (left nodes), used for sign transformation
    kp: usize,
    /// Tree solver for preconditioner application
    tree_solver: TreeSolverOwned,
    /// Residual: r = rhs - G*x
    r: Vec<f64>,
    /// Preconditioned residual: z = G_T^{-1} * r
    z: Vec<f64>,
    /// Search direction
    p: Vec<f64>,
    /// G * p
    ap: Vec<f64>,
    /// Scratch for sign-transformed input to tree solve
    sign_scratch: Vec<f64>,
}

impl FixedPCG {
    /// Create a new fixed-iteration PCG solver.
    ///
    /// # Arguments
    /// * `v` - Number of nodes (problem size)
    /// * `kp` - Number of P (left) nodes, for sign transformation
    /// * `k_inner` - Number of PCG iterations (typically 5-15)
    pub fn new(v: usize, kp: usize, k_inner: usize) -> Self {
        Self {
            k_inner,
            kp,
            tree_solver: TreeSolverOwned::new(v),
            r: vec![0.0; v],
            z: vec![0.0; v],
            p: vec![0.0; v],
            ap: vec![0.0; v],
            sign_scratch: vec![0.0; v],
        }
    }

    /// Apply sign transformation: flip signs of Q nodes (indices kp..v).
    ///
    /// This transforms between Laplacian and Gram representations:
    /// S L S = G, so G^{-1} = S L^{-1} S
    #[inline]
    fn apply_sign_transform(kp: usize, x: &mut [f64]) {
        for xi in x[kp..].iter_mut() {
            *xi = -*xi;
        }
    }


    /// Solve `G * x ≈ P_G * rhs` using k_inner PCG iterations.
    ///
    /// # Arguments
    /// * `gram` - The Gram operator (provides gram_matvec)
    /// * `forest` - Spanning forest for tree preconditioner
    /// * `components` - Component structure for projection
    /// * `rhs` - Right-hand side (arbitrary input)
    /// * `x` - Output solution (overwritten, starts from zero)
    ///
    /// # Note
    /// The RHS is projected to Range(G) before solving because the Gram matrix G
    /// is singular (has alternating-sign nullspace). The projection P_G enforces
    /// sum_P = sum_Q per component.
    ///
    /// **WARNING for LSMR preconditioning**: This projection means the solver
    /// effectively computes x = G^+ * P_G * rhs (pseudoinverse), not G^{-1} * rhs.
    /// For vectors with non-zero nullspace components, this produces incorrect
    /// results as a preconditioner. Consider using `sparse_gram` instead.
    ///
    /// # Returns
    /// The final residual norm (for diagnostics only).
    pub fn solve(
        &mut self,
        gram: &LaplacianOp,
        forest: &SpanningForest,
        components: &ComponentsCSR,
        rhs: &[f64],
        x: &mut [f64],
    ) -> f64 {
        let v = gram.v;
        let kp = self.kp;
        debug_assert_eq!(rhs.len(), v);
        debug_assert_eq!(x.len(), v);

        // Initialize x = 0
        x.fill(0.0);

        // Copy RHS and project to Range(G) for the Gram system.
        // The projection P_G enforces sum_P = sum_Q per component.
        //
        // NOTE: This projection makes the preconditioner NOT a true inverse for
        // vectors outside Range(G). For LSMR preconditioning, this can cause
        // convergence issues. Consider using sparse_gram preconditioner instead.
        self.r.copy_from_slice(rhs);
        components.project_to_gram_range(&mut self.r, kp);

        // Check for zero RHS
        let rhs_norm_sq: f64 = self.r.iter().map(|&ri| ri * ri).sum();
        if rhs_norm_sq < 1e-30 {
            return 0.0;
        }

        // z = G_T^{-1} * P_G * r (tree Gram pseudoinverse)
        // The key identity: G_T^+ = S L_T^+ S where L_T^+ = L_T^{-1} * P_L
        // After projecting r to Range(G), sign-transforming gives S*r in Range(L).
        // So we can use apply_raw (no additional projection needed in L-space).
        self.sign_scratch.copy_from_slice(&self.r);
        Self::apply_sign_transform(kp, &mut self.sign_scratch);
        self.tree_solver.apply_raw(forest, &self.sign_scratch, &mut self.z);
        Self::apply_sign_transform(kp, &mut self.z);

        // p = z (initial search direction)
        self.p.copy_from_slice(&self.z);

        // rz = (r, z)
        let mut rz: f64 = dot(&self.r, &self.z);

        // PCG iterations
        for _ in 0..self.k_inner {
            // ap = G * p (Gram matvec, not Laplacian!)
            gram.gram_matvec(&self.p, &mut self.ap);

            // alpha = rz / (p, ap)
            let pap: f64 = dot(&self.p, &self.ap);
            if pap.abs() < 1e-30 {
                // p is in nullspace, stop
                break;
            }
            let alpha = rz / pap;

            // x = x + alpha * p
            axpy(alpha, &self.p, x);

            // r = r - alpha * ap
            axpy(-alpha, &self.ap, &mut self.r);

            // z = G_T^{-1} * r (via sign transformation)
            // Note: r stays in Range(G) throughout PCG (G maps to Range(G)),
            // so we can use apply_raw without projection.
            self.sign_scratch.copy_from_slice(&self.r);
            Self::apply_sign_transform(kp, &mut self.sign_scratch);
            self.tree_solver.apply_raw(forest, &self.sign_scratch, &mut self.z);
            Self::apply_sign_transform(kp, &mut self.z);

            // rz_new = (r, z)
            let rz_new = dot(&self.r, &self.z);

            // beta = rz_new / rz
            let beta = if rz.abs() > 1e-30 { rz_new / rz } else { 0.0 };

            // p = z + beta * p
            for (pi, &zi) in self.p.iter_mut().zip(self.z.iter()) {
                *pi = zi + beta * *pi;
            }

            rz = rz_new;
        }

        // Return final residual norm
        self.r.iter().map(|&ri| ri * ri).sum::<f64>().sqrt()
    }

    /// Get the number of inner iterations.
    #[inline]
    pub fn k_inner(&self) -> usize {
        self.k_inner
    }
}

/// Dot product of two vectors.
#[inline]
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b.iter()).map(|(&ai, &bi)| ai * bi).sum()
}


/// AXPY: y = y + alpha * x
#[inline]
fn axpy(alpha: f64, x: &[f64], y: &mut [f64]) {
    for (yi, &xi) in y.iter_mut().zip(x.iter()) {
        *yi += alpha * xi;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::union_find::UnionFind;

    fn make_test_problem() -> (LaplacianOp, SpanningForest, ComponentsCSR, usize) {
        // Simple connected graph
        let group_ids_p = vec![0, 0, 1, 1];
        let group_ids_q = vec![0, 1, 1, 2];
        let kp = 2;
        let kq = 3;
        let v = kp + kq;

        // Build operator (provides both Laplacian and Gram matvec)
        let op = LaplacianOp::new_streaming(
            group_ids_p.clone(),
            group_ids_q.clone(),
            kp,
            kq,
            None,
        );

        // Build components
        let mut uf = UnionFind::new(v);
        for (&gp, &gq) in group_ids_p.iter().zip(group_ids_q.iter()) {
            uf.union(gp, kp + gq);
        }
        let components = ComponentsCSR::from_union_find(&mut uf, v);

        // Build forest
        let forest = SpanningForest::build_simple(
            &group_ids_p,
            &group_ids_q,
            kp,
            kq,
            None,
            &components,
        );

        (op, forest, components, kp)
    }

    #[test]
    fn test_pcg_zero_rhs() {
        let (op, forest, components, kp) = make_test_problem();
        let v = op.v;

        let mut pcg = FixedPCG::new(v, kp, 10);
        let rhs = vec![0.0; v];
        let mut x = vec![1.0; v];

        pcg.solve(&op, &forest, &components, &rhs, &mut x);

        // With zero RHS, solution should be zero
        for &xi in &x {
            assert!(xi.abs() < 1e-10);
        }
    }

    #[test]
    fn test_pcg_convergence() {
        let (op, forest, components, kp) = make_test_problem();
        let v = op.v;

        // Create RHS with zero sum (in Range(G))
        let mut rhs = vec![1.0, -1.0, 0.5, -0.5, 0.0];
        components.project_to_range(&mut rhs);

        let mut pcg = FixedPCG::new(v, kp, 20);
        let mut x = vec![0.0; v];

        let residual = pcg.solve(&op, &forest, &components, &rhs, &mut x);

        // Should converge well
        assert!(residual < 1e-6, "Residual too large: {}", residual);

        // Verify: G * x ≈ rhs (projected)
        let mut gx = vec![0.0; v];
        op.gram_matvec(&x, &mut gx);

        let mut rhs_proj = rhs.clone();
        components.project_to_range(&mut rhs_proj);

        for i in 0..v {
            assert!(
                (gx[i] - rhs_proj[i]).abs() < 1e-5,
                "Mismatch at {}: gx={}, rhs={}",
                i,
                gx[i],
                rhs_proj[i]
            );
        }
    }

    #[test]
    fn test_pcg_projects_rhs_to_gram_range() {
        let (op, forest, components, kp) = make_test_problem();
        let v = op.v;

        // RHS not in Range(G): sum_P != sum_Q
        let rhs = vec![1.0, 1.0, 0.0, 0.0, 0.0];

        let mut pcg = FixedPCG::new(v, kp, 20);
        let mut x = vec![0.0; v];

        let residual = pcg.solve(&op, &forest, &components, &rhs, &mut x);
        assert!(residual < 1e-6, "Residual too large: {}", residual);

        let mut gx = vec![0.0; v];
        op.gram_matvec(&x, &mut gx);

        let mut rhs_proj = rhs.clone();
        components.project_to_gram_range(&mut rhs_proj, kp);

        for i in 0..v {
            assert!(
                (gx[i] - rhs_proj[i]).abs() < 1e-5,
                "Mismatch at {}: gx={}, rhs_proj={}",
                i,
                gx[i],
                rhs_proj[i]
            );
        }
    }

    #[test]
    fn test_pcg_few_iterations() {
        let (op, forest, components, kp) = make_test_problem();
        let v = op.v;

        // Even with few iterations, should make progress
        let mut rhs = vec![1.0, -1.0, 0.5, -0.5, 0.0];
        components.project_to_range(&mut rhs);

        let mut pcg = FixedPCG::new(v, kp, 3);
        let mut x = vec![0.0; v];

        let residual = pcg.solve(&op, &forest, &components, &rhs, &mut x);

        // Should reduce residual somewhat
        let rhs_norm: f64 = rhs.iter().map(|&r| r * r).sum::<f64>().sqrt();
        assert!(residual < rhs_norm, "PCG should reduce residual");
    }

    #[test]
    fn test_pcg_linearity() {
        // Test that PCG is linear: solve(a*r1 + b*r2) ≈ a*solve(r1) + b*solve(r2)
        // Note: Fixed-iteration PCG is exactly linear because it performs the
        // same number of iterations regardless of input, and all operations are linear.
        let (op, forest, components, kp) = make_test_problem();
        let v = op.v;

        let mut r1 = vec![1.0, -0.5, 0.0, -0.5, 0.0];
        let mut r2 = vec![0.0, 0.5, -1.0, 0.5, 0.0];
        components.project_to_range(&mut r1);
        components.project_to_range(&mut r2);

        let a = 2.0;
        let b = -1.5;

        // Use different PCG instances to ensure no state interference
        let mut pcg1 = FixedPCG::new(v, kp, 10);
        let mut pcg2 = FixedPCG::new(v, kp, 10);
        let mut pcg3 = FixedPCG::new(v, kp, 10);

        // Solve individually
        let mut x1 = vec![0.0; v];
        let mut x2 = vec![0.0; v];
        pcg1.solve(&op, &forest, &components, &r1, &mut x1);
        pcg2.solve(&op, &forest, &components, &r2, &mut x2);

        // Solve combined
        let r_combined: Vec<f64> = r1.iter().zip(r2.iter()).map(|(&r1i, &r2i)| a * r1i + b * r2i).collect();
        let mut x_combined = vec![0.0; v];
        pcg3.solve(&op, &forest, &components, &r_combined, &mut x_combined);

        // Expected: a*x1 + b*x2
        let x_expected: Vec<f64> = x1.iter().zip(x2.iter()).map(|(&x1i, &x2i)| a * x1i + b * x2i).collect();

        // Check linearity (up to gauge - compare differences)
        for i in 1..v {
            let diff_combined = x_combined[i] - x_combined[0];
            let diff_expected = x_expected[i] - x_expected[0];
            assert!(
                (diff_combined - diff_expected).abs() < 1e-8,
                "Linearity violation at {}: combined_diff={}, expected_diff={}",
                i,
                diff_combined,
                diff_expected
            );
        }
    }
}
