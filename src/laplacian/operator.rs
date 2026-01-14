//! Laplacian operator for bipartite 2-FE graphs.
//!
//! The Laplacian matrix L for a bipartite graph has the property that
//! for any vector x: (Lx)[u] = sum over edges (u,v,w) of w*(x[u] - x[v]).
//!
//! The nullspace of L consists of vectors that are constant on each
//! connected component. For solving Lx = r, we need r to be in Range(L),
//! which means r must have zero sum on each component.

use super::edges::{AggregatedEdges, EdgeSource, StreamingEdges};
use crate::graph::components::ComponentsCSR;
use crate::graph::union_find::UnionFind;

/// Bipartite graph Laplacian operator for a 2-FE pair.
///
/// Nodes are numbered:
/// - Left nodes (FE p groups): 0..kp
/// - Right nodes (FE q groups): kp..kp+kq
/// - Total nodes: V = kp + kq
///
/// Each observation creates an edge between its FE p group and FE q group.
pub struct LaplacianOp {
    /// Total number of nodes: V = kp + kq
    pub v: usize,
    /// Number of groups in FE p (left nodes)
    pub kp: usize,
    /// Number of groups in FE q (right nodes)
    pub kq: usize,
    /// Edge representation (streaming or aggregated)
    edges: Box<dyn EdgeSource + Send + Sync>,
    /// Connected component structure for gauge projection
    pub components: ComponentsCSR,
}

impl LaplacianOp {
    /// Create a Laplacian operator using streaming edges.
    ///
    /// Streaming edges iterate directly over FE index arrays with no
    /// additional storage. Matvec cost is O(n_obs).
    ///
    /// # Arguments
    /// * `group_ids_p` - FE p group IDs per observation (values in 0..kp)
    /// * `group_ids_q` - FE q group IDs per observation (values in 0..kq)
    /// * `kp` - Number of groups in FE p
    /// * `kq` - Number of groups in FE q
    /// * `weights` - Optional observation weights
    pub fn new_streaming(
        group_ids_p: Vec<usize>,
        group_ids_q: Vec<usize>,
        kp: usize,
        kq: usize,
        weights: Option<Vec<f64>>,
    ) -> Self {
        let v = kp + kq;

        // Build connected components
        let components = Self::build_components(&group_ids_p, &group_ids_q, kp, v);

        let edges = Box::new(StreamingEdges::new(group_ids_p, group_ids_q, kp, weights));

        Self {
            v,
            kp,
            kq,
            edges,
            components,
        }
    }

    /// Create a Laplacian operator using aggregated (deduplicated) edges.
    ///
    /// Aggregated edges store unique (u,v) pairs with summed weights.
    /// Matvec cost is O(P_unique) which can be << O(n_obs) for many repeats.
    ///
    /// # Arguments
    /// * `group_ids_p` - FE p group IDs per observation (values in 0..kp)
    /// * `group_ids_q` - FE q group IDs per observation (values in 0..kq)
    /// * `kp` - Number of groups in FE p
    /// * `kq` - Number of groups in FE q
    /// * `weights` - Optional observation weights
    pub fn new_aggregated(
        group_ids_p: &[usize],
        group_ids_q: &[usize],
        kp: usize,
        kq: usize,
        weights: Option<&[f64]>,
    ) -> Self {
        let v = kp + kq;

        // Build connected components
        let components = Self::build_components(group_ids_p, group_ids_q, kp, v);

        let edges = Box::new(AggregatedEdges::from_indices(
            group_ids_p,
            group_ids_q,
            kp,
            weights,
        ));

        Self {
            v,
            kp,
            kq,
            edges,
            components,
        }
    }

    /// Build connected components from FE index arrays.
    fn build_components(
        group_ids_p: &[usize],
        group_ids_q: &[usize],
        kp: usize,
        v: usize,
    ) -> ComponentsCSR {
        let mut uf = UnionFind::new(v);

        // Each observation is an edge connecting left node u to right node v
        for (&gp, &gq) in group_ids_p.iter().zip(group_ids_q.iter()) {
            let u = gp;
            let node_v = kp + gq;
            uf.union(u, node_v);
        }

        ComponentsCSR::from_union_find(&mut uf, v)
    }

    /// Apply Laplacian: y = L * x (L = Deg - Adj).
    ///
    /// Zeros y first, then accumulates the result.
    ///
    /// # Arguments
    /// * `x` - Input vector (length V)
    /// * `y` - Output vector (length V), overwritten
    #[inline]
    pub fn laplacian_matvec(&self, x: &[f64], y: &mut [f64]) {
        debug_assert_eq!(x.len(), self.v);
        debug_assert_eq!(y.len(), self.v);

        // Zero output
        y.iter_mut().for_each(|yi| *yi = 0.0);

        // Accumulate Laplacian
        self.edges.laplacian_matvec(x, y);
    }

    /// Apply Gram matrix: y = G * x (G = Deg + Adj).
    ///
    /// This is the correct operator for preconditioning the normal equations
    /// in LSMR. Zeros y first, then accumulates the result.
    ///
    /// # Arguments
    /// * `x` - Input vector (length V)
    /// * `y` - Output vector (length V), overwritten
    #[inline]
    pub fn gram_matvec(&self, x: &[f64], y: &mut [f64]) {
        debug_assert_eq!(x.len(), self.v);
        debug_assert_eq!(y.len(), self.v);

        // Zero output
        y.iter_mut().for_each(|yi| *yi = 0.0);

        // Accumulate Gram matrix
        self.edges.gram_matvec(x, y);
    }

    /// Alias for `laplacian_matvec` (backward compatibility).
    #[inline]
    pub fn matvec(&self, x: &[f64], y: &mut [f64]) {
        self.laplacian_matvec(x, y);
    }

    /// Apply Laplacian, accumulating into y (does not zero y first).
    ///
    /// # Arguments
    /// * `x` - Input vector (length V)
    /// * `y` - Output vector (length V), accumulated into
    #[inline]
    #[allow(dead_code)]
    pub fn matvec_accumulate(&self, x: &[f64], y: &mut [f64]) {
        debug_assert_eq!(x.len(), self.v);
        debug_assert_eq!(y.len(), self.v);
        self.edges.laplacian_matvec(x, y);
    }

    /// Project vector into Range(L) by subtracting component-wise means.
    ///
    /// For the Laplacian system Lx = r to be consistent, r must be in Range(L),
    /// meaning r must have zero sum on each connected component.
    ///
    /// This method modifies r in-place to satisfy this constraint.
    ///
    /// # Arguments
    /// * `r` - Vector to project (length V), modified in-place
    #[inline]
    pub fn project_to_range(&self, r: &mut [f64]) {
        self.components.project_to_range(r);
    }

    /// Get number of edges in the underlying representation.
    ///
    /// For streaming: n_obs. For aggregated: number of unique edges.
    #[inline]
    pub fn num_edges(&self) -> usize {
        self.edges.num_edges()
    }

    /// Get the number of connected components.
    #[inline]
    pub fn num_components(&self) -> usize {
        self.components.num_components
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_simple_bipartite() -> (Vec<usize>, Vec<usize>, usize, usize) {
        // Simple bipartite graph:
        // Left nodes: 0, 1 (kp=2)
        // Right nodes: 2, 3 (kq=2, shifted to indices 2,3)
        // Edges: (0,0), (0,1), (1,1) meaning observations connect:
        //   obs 0: left[0] - right[0]
        //   obs 1: left[0] - right[1]
        //   obs 2: left[1] - right[1]
        let group_ids_p = vec![0, 0, 1];
        let group_ids_q = vec![0, 1, 1];
        let kp = 2;
        let kq = 2;
        (group_ids_p, group_ids_q, kp, kq)
    }

    #[test]
    fn test_laplacian_op_streaming_basic() {
        let (group_ids_p, group_ids_q, kp, kq) = make_simple_bipartite();
        let op = LaplacianOp::new_streaming(group_ids_p, group_ids_q, kp, kq, None);

        assert_eq!(op.v, 4);
        assert_eq!(op.kp, 2);
        assert_eq!(op.kq, 2);
        assert_eq!(op.num_edges(), 3);
        assert_eq!(op.num_components(), 1); // All connected
    }

    #[test]
    fn test_laplacian_op_aggregated_basic() {
        let (group_ids_p, group_ids_q, kp, kq) = make_simple_bipartite();
        let op = LaplacianOp::new_aggregated(&group_ids_p, &group_ids_q, kp, kq, None);

        assert_eq!(op.v, 4);
        assert_eq!(op.num_components(), 1);
        // 3 unique edges: (0,2), (0,3), (1,3)
        assert_eq!(op.num_edges(), 3);
    }

    #[test]
    fn test_laplacian_matvec() {
        let (group_ids_p, group_ids_q, kp, kq) = make_simple_bipartite();
        let op = LaplacianOp::new_streaming(group_ids_p, group_ids_q, kp, kq, None);

        // x = [1, 2, 3, 4] for nodes [p0, p1, q0, q1]
        let x = vec![1.0, 2.0, 3.0, 4.0];
        let mut y = vec![0.0; 4];
        op.matvec(&x, &mut y);

        // Edge (0,2): d = 1-3 = -2, y[0] += -2, y[2] -= -2
        // Edge (0,3): d = 1-4 = -3, y[0] += -3, y[3] -= -3
        // Edge (1,3): d = 2-4 = -2, y[1] += -2, y[3] -= -2
        // y = [-2-3, -2, 2, 3+2] = [-5, -2, 2, 5]
        assert_eq!(y, vec![-5.0, -2.0, 2.0, 5.0]);
    }

    #[test]
    fn test_laplacian_symmetry() {
        // L should be symmetric: <Lx, y> = <x, Ly>
        let (group_ids_p, group_ids_q, kp, kq) = make_simple_bipartite();
        let op = LaplacianOp::new_streaming(group_ids_p, group_ids_q, kp, kq, None);

        let x = vec![1.0, 2.0, 3.0, 4.0];
        let y = vec![5.0, 6.0, 7.0, 8.0];

        let mut lx = vec![0.0; 4];
        let mut ly = vec![0.0; 4];
        op.matvec(&x, &mut lx);
        op.matvec(&y, &mut ly);

        // <Lx, y>
        let lx_dot_y: f64 = lx.iter().zip(y.iter()).map(|(&a, &b)| a * b).sum();
        // <x, Ly>
        let x_dot_ly: f64 = x.iter().zip(ly.iter()).map(|(&a, &b)| a * b).sum();

        assert!((lx_dot_y - x_dot_ly).abs() < 1e-10);
    }

    #[test]
    fn test_laplacian_psd() {
        // L should be positive semi-definite: <Lx, x> >= 0
        let (group_ids_p, group_ids_q, kp, kq) = make_simple_bipartite();
        let op = LaplacianOp::new_streaming(group_ids_p, group_ids_q, kp, kq, None);

        let x = vec![1.0, -2.0, 3.0, -4.0];
        let mut lx = vec![0.0; 4];
        op.matvec(&x, &mut lx);

        let x_dot_lx: f64 = x.iter().zip(lx.iter()).map(|(&a, &b)| a * b).sum();
        assert!(x_dot_lx >= 0.0);
    }

    #[test]
    fn test_laplacian_nullspace() {
        // Constant vector should be in nullspace: L * 1 = 0
        let (group_ids_p, group_ids_q, kp, kq) = make_simple_bipartite();
        let op = LaplacianOp::new_streaming(group_ids_p, group_ids_q, kp, kq, None);

        let ones = vec![1.0; 4];
        let mut lones = vec![0.0; 4];
        op.matvec(&ones, &mut lones);

        for &val in &lones {
            assert!(val.abs() < 1e-10);
        }
    }

    #[test]
    fn test_project_to_range() {
        let (group_ids_p, group_ids_q, kp, kq) = make_simple_bipartite();
        let op = LaplacianOp::new_streaming(group_ids_p, group_ids_q, kp, kq, None);

        let mut r = vec![1.0, 2.0, 3.0, 4.0];
        op.project_to_range(&mut r);

        // After projection, sum should be zero (single component)
        let sum: f64 = r.iter().sum();
        assert!(sum.abs() < 1e-10);
    }

    #[test]
    fn test_multiple_components() {
        // Two disconnected components:
        // Component 1: left[0] - right[0]
        // Component 2: left[1] - right[1]
        let group_ids_p = vec![0, 1];
        let group_ids_q = vec![0, 1];
        let kp = 2;
        let kq = 2;

        let op = LaplacianOp::new_streaming(group_ids_p, group_ids_q, kp, kq, None);

        assert_eq!(op.num_components(), 2);

        // Constant on each component should be in nullspace
        // x = [1, 2, 1, 2] - constant 1 on comp1 (nodes 0,2), constant 2 on comp2 (nodes 1,3)
        let x = vec![1.0, 2.0, 1.0, 2.0];
        let mut lx = vec![0.0; 4];
        op.matvec(&x, &mut lx);

        for &val in &lx {
            assert!(val.abs() < 1e-10);
        }
    }

    #[test]
    fn test_streaming_vs_aggregated_equivalence() {
        let group_ids_p = vec![0, 1, 0, 2, 1, 0];
        let group_ids_q = vec![0, 1, 0, 2, 1, 1];
        let kp = 3;
        let kq = 3;

        let streaming = LaplacianOp::new_streaming(
            group_ids_p.clone(),
            group_ids_q.clone(),
            kp,
            kq,
            None,
        );
        let aggregated = LaplacianOp::new_aggregated(&group_ids_p, &group_ids_q, kp, kq, None);

        let x: Vec<f64> = (0..6).map(|i| i as f64 + 1.0).collect();

        let mut y_streaming = vec![0.0; 6];
        let mut y_aggregated = vec![0.0; 6];

        streaming.matvec(&x, &mut y_streaming);
        aggregated.matvec(&x, &mut y_aggregated);

        for (ys, ya) in y_streaming.iter().zip(y_aggregated.iter()) {
            assert!(
                (ys - ya).abs() < 1e-10,
                "Mismatch: streaming={}, aggregated={}",
                ys,
                ya
            );
        }
    }
}
