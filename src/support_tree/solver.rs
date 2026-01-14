//! Linear-time tree Laplacian solver.
//!
//! Solves `L_T * x = r` where `L_T` is the tree Laplacian in O(V) time
//! using two passes: postorder (compute flows) and preorder (recover potentials).
//!
//! # Requirements
//! - `r` must have zero sum on each component (use `project_to_range` first)
//! - Solution `x` has gauge: `x[root] = 0` for each tree root

use super::forest_builder::SpanningForest;

/// Linear-time tree Laplacian solver.
///
/// Solves `L_T * x = r` in O(V) time where `L_T` is the Laplacian of the
/// spanning forest. Uses the classical two-pass algorithm:
///
/// 1. **Postorder (flow computation)**: For each node (children first),
///    `flow[u] = r[u] + sum_{c in children(u)} flow[c]`
///
/// 2. **Preorder (potential recovery)**: For root, `x[root] = 0`.
///    For each child, `x[child] = x[parent] + flow[child] / w(child, parent)`
///
/// # Gauge
/// The solution is unique up to an additive constant per component.
/// This solver uses the gauge `x[root] = 0`.
pub struct TreeSolver {
    /// Scratch buffer for flow computation (length V)
    flow: Vec<f64>,
}

impl TreeSolver {
    /// Create a new tree solver for a forest with V nodes.
    pub fn new(v: usize) -> Self {
        Self { flow: vec![0.0; v] }
    }

    /// Solve `L_T * x = r` where `r` has zero sum per component.
    ///
    /// # Arguments
    /// * `forest` - The spanning forest structure
    /// * `r` - Right-hand side (must have zero sum on each component)
    /// * `x` - Output solution (overwritten)
    ///
    /// # Algorithm
    /// 1. Postorder pass: `flow[u] = r[u] + sum_{c in children(u)} flow[c]`
    /// 2. Preorder pass: `x[root] = 0; x[u] = x[parent] + flow[u] / w(u, parent)`
    ///
    /// # Complexity
    /// O(V) time, O(V) space for flow buffer.
    pub fn solve(&mut self, forest: &SpanningForest, r: &[f64], x: &mut [f64]) {
        let v = forest.v;
        debug_assert_eq!(r.len(), v);
        debug_assert_eq!(x.len(), v);
        debug_assert_eq!(self.flow.len(), v);

        // Pass 1: Postorder - compute flows
        // flow[u] = r[u] + sum of children's flows
        for &u in &forest.postorder {
            let mut f = r[u];
            for &child in forest.children_of(u) {
                f += self.flow[child];
            }
            self.flow[u] = f;
        }

        // Pass 2: Preorder - recover potentials
        // x[root] = 0
        // x[child] = x[parent] + flow[child] / w(child, parent)
        for &u in &forest.preorder {
            if forest.is_root(u) {
                x[u] = 0.0;
            } else {
                let parent = forest.parent[u];
                let w = forest.parent_weight[u];
                x[u] = x[parent] + self.flow[u] / w;
            }
        }
    }

    /// Apply tree preconditioner: `z = L_T^{-1} * r`.
    ///
    /// This is a convenience wrapper that:
    /// 1. Copies `r` to internal buffer
    /// 2. Projects to range (subtracts component means)
    /// 3. Solves the tree system
    ///
    /// # Arguments
    /// * `forest` - The spanning forest
    /// * `components` - Component structure for projection
    /// * `r` - Input vector
    /// * `z` - Output vector
    pub fn apply(
        &mut self,
        forest: &SpanningForest,
        components: &crate::graph::components::ComponentsCSR,
        r: &[f64],
        z: &mut [f64],
    ) {
        // Copy r to z (we'll modify z in-place for projection)
        z.copy_from_slice(r);

        // Project to range (ensure consistency)
        components.project_to_range(z);

        // Solve (using z as both projected RHS and output)
        // We need a temporary for the projected RHS since solve reads r and writes x
        // But we can reuse flow buffer temporarily
        let v = forest.v;

        // Actually, we can solve directly since z is the projected RHS
        // and we want to overwrite z with the solution
        // We need to be careful: solve() reads r, writes x
        // So we need to copy z to a temp buffer first

        // Use a simple approach: solve with z as RHS, output to self.flow,
        // then copy back
        self.solve(forest, z, &mut self.flow.clone());

        // Hmm, that's wasteful. Let me restructure.
        // Actually the cleanest is to have a separate temp buffer for the RHS
        // For now, let's just do the allocation
        let mut rhs = z.to_vec();
        self.solve(forest, &rhs, z);
        let _ = rhs; // suppress warning
    }

    /// Reset the solver for reuse (clears flow buffer).
    #[inline]
    pub fn reset(&mut self) {
        self.flow.fill(0.0);
    }
}

/// Owned tree solver with all necessary buffers.
///
/// This version owns additional buffers to avoid allocation in `apply`.
pub struct TreeSolverOwned {
    /// Scratch buffer for flow computation
    flow: Vec<f64>,
    /// Scratch buffer for projected RHS
    rhs_scratch: Vec<f64>,
}

impl TreeSolverOwned {
    /// Create a new owned tree solver for V nodes.
    pub fn new(v: usize) -> Self {
        Self {
            flow: vec![0.0; v],
            rhs_scratch: vec![0.0; v],
        }
    }

    /// Solve `L_T * x = r` (assumes r already projected to range).
    pub fn solve(&mut self, forest: &SpanningForest, r: &[f64], x: &mut [f64]) {
        let v = forest.v;
        debug_assert_eq!(r.len(), v);
        debug_assert_eq!(x.len(), v);

        // Pass 1: Postorder - compute flows
        for &u in &forest.postorder {
            let mut f = r[u];
            for &child in forest.children_of(u) {
                f += self.flow[child];
            }
            self.flow[u] = f;
        }

        // Pass 2: Preorder - recover potentials
        for &u in &forest.preorder {
            if forest.is_root(u) {
                x[u] = 0.0;
            } else {
                let parent = forest.parent[u];
                let w = forest.parent_weight[u];
                x[u] = x[parent] + self.flow[u] / w;
            }
        }
    }

    /// Apply tree preconditioner: z = L_T^{-1} * r.
    ///
    /// Projects r to range, then solves. Allocation-free.
    pub fn apply(
        &mut self,
        forest: &SpanningForest,
        components: &crate::graph::components::ComponentsCSR,
        r: &[f64],
        z: &mut [f64],
    ) {
        // Copy r to scratch and project
        self.rhs_scratch.copy_from_slice(r);
        components.project_to_range(&mut self.rhs_scratch);

        // Solve (inlined to avoid borrow conflict)
        self.solve_inner(forest, z);
    }

    /// Apply tree solve without projection: z = L_T^{-1} * r.
    ///
    /// This version does NOT project r to Range(L_T). Use this when:
    /// - r is already known to be in Range(L_T)
    /// - You need the raw solve without any gauge fixing
    /// - You're composing with other linear operations that handle the gauge
    ///
    /// The solution still uses gauge x[root] = 0.
    pub fn apply_raw(
        &mut self,
        forest: &SpanningForest,
        r: &[f64],
        z: &mut [f64],
    ) {
        // Copy r to scratch (no projection)
        self.rhs_scratch.copy_from_slice(r);

        // Solve (inlined to avoid borrow conflict)
        self.solve_inner(forest, z);
    }

    /// Internal solve helper - assumes rhs_scratch is already populated.
    #[inline]
    fn solve_inner(&mut self, forest: &SpanningForest, z: &mut [f64]) {
        // Pass 1: Postorder - compute flows
        for &u in &forest.postorder {
            let mut f = self.rhs_scratch[u];
            for &child in forest.children_of(u) {
                f += self.flow[child];
            }
            self.flow[u] = f;
        }

        // Pass 2: Preorder - recover potentials
        for &u in &forest.preorder {
            if forest.is_root(u) {
                z[u] = 0.0;
            } else {
                let parent = forest.parent[u];
                let w = forest.parent_weight[u];
                z[u] = z[parent] + self.flow[u] / w;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::graph::components::ComponentsCSR;
    use crate::graph::union_find::UnionFind;

    fn make_simple_forest() -> (SpanningForest, ComponentsCSR) {
        // Simple connected graph: 4 nodes, 3 edges forming a tree
        // Edges: (0,2), (0,3), (1,3) with weights 1.0
        let group_ids_p = vec![0, 0, 1];
        let group_ids_q = vec![0, 1, 1];
        let kp = 2;
        let kq = 2;
        let v = kp + kq;

        let mut uf = UnionFind::new(v);
        for (&gp, &gq) in group_ids_p.iter().zip(group_ids_q.iter()) {
            uf.union(gp, kp + gq);
        }
        let components = ComponentsCSR::from_union_find(&mut uf, v);

        let forest = SpanningForest::build_simple(
            &group_ids_p,
            &group_ids_q,
            kp,
            kq,
            None,
            &components,
        );

        (forest, components)
    }

    #[test]
    fn test_tree_solver_zero_rhs() {
        let (forest, _components) = make_simple_forest();
        let mut solver = TreeSolver::new(4);

        let r = vec![0.0; 4];
        let mut x = vec![1.0; 4]; // Initialize to non-zero

        solver.solve(&forest, &r, &mut x);

        // With zero RHS, solution should be zero (gauge: x[root]=0)
        for &xi in &x {
            assert!(xi.abs() < 1e-10);
        }
    }

    #[test]
    fn test_tree_solver_roundtrip() {
        let (forest, components) = make_simple_forest();
        let mut solver = TreeSolver::new(4);

        // Create a random x, compute r = L_T * x, then solve to recover x
        // First, we need to compute L_T * x
        // For tree Laplacian: (L_T * x)[u] = sum over edges (u,v) of w*(x[u] - x[v])

        // Use a simple x
        let x_orig = vec![1.0, 2.0, 3.0, 4.0];

        // Compute r = L_T * x manually
        // For tree edges: we need to look at forest.parent and parent_weight
        let mut r = vec![0.0; 4];
        for u in 0..4 {
            if !forest.is_root(u) {
                let p = forest.parent[u];
                let w = forest.parent_weight[u];
                let d = w * (x_orig[u] - x_orig[p]);
                r[u] += d;
                r[p] -= d;
            }
        }

        // Project r to range (should already be in range if L_T * x)
        let sum: f64 = r.iter().sum();
        assert!(sum.abs() < 1e-10, "r should have zero sum");

        // Solve
        let mut x_recovered = vec![0.0; 4];
        solver.solve(&forest, &r, &mut x_recovered);

        // x_recovered should equal x_orig minus the component mean
        // (because of gauge x[root] = 0)
        let mean: f64 = x_orig.iter().sum::<f64>() / 4.0;
        for i in 0..4 {
            let expected = x_orig[i] - mean;
            // Allow for different gauge choice
            let diff = (x_recovered[i] - expected).abs();
            // The solution is unique up to constant, so check differences
            if i > 0 {
                let orig_diff = x_orig[i] - x_orig[0];
                let recov_diff = x_recovered[i] - x_recovered[0];
                assert!(
                    (orig_diff - recov_diff).abs() < 1e-10,
                    "Difference mismatch at i={}: orig_diff={}, recov_diff={}",
                    i,
                    orig_diff,
                    recov_diff
                );
            }
        }
    }

    #[test]
    fn test_tree_solver_owned_apply() {
        let (forest, components) = make_simple_forest();
        let mut solver = TreeSolverOwned::new(4);

        // r with non-zero mean - apply should project it first
        let r = vec![1.0, 2.0, 3.0, 4.0];
        let mut z = vec![0.0; 4];

        solver.apply(&forest, &components, &r, &mut z);

        // After solving L_T * z = projected(r), we should have:
        // L_T * z = r - mean(r) * ones

        // Verify by computing L_T * z
        let mut lz = vec![0.0; 4];
        for u in 0..4 {
            if !forest.is_root(u) {
                let p = forest.parent[u];
                let w = forest.parent_weight[u];
                let d = w * (z[u] - z[p]);
                lz[u] += d;
                lz[p] -= d;
            }
        }

        // lz should equal r - mean(r)
        let mean = r.iter().sum::<f64>() / 4.0;
        for i in 0..4 {
            let expected = r[i] - mean;
            assert!(
                (lz[i] - expected).abs() < 1e-10,
                "Mismatch at i={}: lz={}, expected={}",
                i,
                lz[i],
                expected
            );
        }
    }

    #[test]
    fn test_tree_solver_two_components() {
        // Two disconnected components
        let group_ids_p = vec![0, 1];
        let group_ids_q = vec![0, 1];
        let kp = 2;
        let kq = 2;
        let v = kp + kq;

        let mut uf = UnionFind::new(v);
        for (&gp, &gq) in group_ids_p.iter().zip(group_ids_q.iter()) {
            uf.union(gp, kp + gq);
        }
        let components = ComponentsCSR::from_union_find(&mut uf, v);

        let forest = SpanningForest::build_simple(
            &group_ids_p,
            &group_ids_q,
            kp,
            kq,
            None,
            &components,
        );

        let mut solver = TreeSolverOwned::new(4);

        // r with zero mean per component
        // Component 1: nodes 0, 2 -> r[0] + r[2] = 0
        // Component 2: nodes 1, 3 -> r[1] + r[3] = 0
        let r = vec![1.0, 2.0, -1.0, -2.0];
        let mut z = vec![0.0; 4];

        solver.solve(&forest, &r, &mut z);

        // Verify L_T * z = r
        let mut lz = vec![0.0; 4];
        for u in 0..4 {
            if !forest.is_root(u) {
                let p = forest.parent[u];
                let w = forest.parent_weight[u];
                let d = w * (z[u] - z[p]);
                lz[u] += d;
                lz[p] -= d;
            }
        }

        for i in 0..4 {
            assert!(
                (lz[i] - r[i]).abs() < 1e-10,
                "Mismatch at i={}: lz={}, r={}",
                i,
                lz[i],
                r[i]
            );
        }
    }
}
