//! Edge representations for Laplacian operators.
//!
//! This module provides two edge representations:
//! - [`StreamingEdges`]: Iterates directly over FE index arrays (no extra storage)
//! - [`AggregatedEdges`]: Deduplicated edges with summed weights (faster for many repeats)

/// Trait for edge iteration in bipartite graph operators.
///
/// Abstracts over streaming (from indices) vs aggregated (deduplicated) edges.
///
/// Two matrix-vector products are provided:
/// - `laplacian_matvec`: L = Deg - Adj (graph Laplacian, negative off-diagonal)
/// - `gram_matvec`: G = Deg + Adj (Gram matrix, positive off-diagonal)
///
/// The Gram matrix is used for preconditioning the normal equations in LSMR.
pub trait EdgeSource: Send + Sync {
    /// Number of edges.
    ///
    /// For streaming: n_obs (one edge per observation).
    /// For aggregated: number of unique (u, v) pairs.
    fn num_edges(&self) -> usize;

    /// Apply Laplacian matvec: y += L * x (accumulates into y).
    /// L = Deg - Adj (negative off-diagonal).
    ///
    /// Caller must zero y before calling if fresh result is needed.
    fn laplacian_matvec(&self, x: &[f64], y: &mut [f64]);

    /// Apply Gram matvec: y += G * x (accumulates into y).
    /// G = Deg + Adj (positive off-diagonal).
    ///
    /// This is the correct operator for preconditioning the normal equations.
    /// Caller must zero y before calling if fresh result is needed.
    fn gram_matvec(&self, x: &[f64], y: &mut [f64]);
}

// =============================================================================
// StreamingEdges: Iterate directly over FE index arrays
// =============================================================================

/// Streaming edges directly from FE index arrays.
///
/// No additional storage beyond references to the original arrays.
/// Matvec cost is O(n_obs) regardless of how many edges are repeated.
///
/// Node numbering:
/// - Left nodes (FE p): 0..kp
/// - Right nodes (FE q): kp..kp+kq
pub struct StreamingEdges {
    /// Group IDs for FE p (left nodes): u = group_ids_p[i]
    group_ids_p: Vec<usize>,
    /// Group IDs for FE q (right nodes): v = kp + group_ids_q[i]
    group_ids_q: Vec<usize>,
    /// Number of groups in FE p (used to shift q nodes)
    kp: usize,
    /// Observation weights (None = all edges have weight 1.0)
    weights: Option<Vec<f64>>,
}

impl StreamingEdges {
    /// Create streaming edges from FE index arrays.
    ///
    /// # Arguments
    /// * `group_ids_p` - Group IDs for FE p (length n_obs, values in 0..kp)
    /// * `group_ids_q` - Group IDs for FE q (length n_obs, values in 0..kq)
    /// * `kp` - Number of groups in FE p
    /// * `weights` - Optional observation weights (length n_obs if Some)
    pub fn new(
        group_ids_p: Vec<usize>,
        group_ids_q: Vec<usize>,
        kp: usize,
        weights: Option<Vec<f64>>,
    ) -> Self {
        debug_assert_eq!(group_ids_p.len(), group_ids_q.len());
        if let Some(ref w) = weights {
            debug_assert_eq!(w.len(), group_ids_p.len());
        }

        Self {
            group_ids_p,
            group_ids_q,
            kp,
            weights,
        }
    }

    /// Get the number of observations (edges).
    #[inline]
    pub fn n_obs(&self) -> usize {
        self.group_ids_p.len()
    }
}

impl EdgeSource for StreamingEdges {
    #[inline]
    fn num_edges(&self) -> usize {
        self.group_ids_p.len()
    }

    fn laplacian_matvec(&self, x: &[f64], y: &mut [f64]) {
        let kp = self.kp;

        match &self.weights {
            None => {
                // Unweighted: all edges have weight 1.0
                for (&gp, &gq) in self.group_ids_p.iter().zip(self.group_ids_q.iter()) {
                    let u = gp;
                    let v = kp + gq;
                    let d = x[u] - x[v];
                    y[u] += d;
                    y[v] -= d;
                }
            }
            Some(weights) => {
                // Weighted edges
                for ((&gp, &gq), &w) in self
                    .group_ids_p
                    .iter()
                    .zip(self.group_ids_q.iter())
                    .zip(weights.iter())
                {
                    let u = gp;
                    let v = kp + gq;
                    let d = w * (x[u] - x[v]);
                    y[u] += d;
                    y[v] -= d;
                }
            }
        }
    }

    fn gram_matvec(&self, x: &[f64], y: &mut [f64]) {
        let kp = self.kp;

        match &self.weights {
            None => {
                // Unweighted: all edges have weight 1.0
                // G = Deg + Adj: y[u] += x[u] + x[v], y[v] += x[u] + x[v]
                for (&gp, &gq) in self.group_ids_p.iter().zip(self.group_ids_q.iter()) {
                    let u = gp;
                    let v = kp + gq;
                    let s = x[u] + x[v];
                    y[u] += s;
                    y[v] += s;
                }
            }
            Some(weights) => {
                // Weighted edges
                for ((&gp, &gq), &w) in self
                    .group_ids_p
                    .iter()
                    .zip(self.group_ids_q.iter())
                    .zip(weights.iter())
                {
                    let u = gp;
                    let v = kp + gq;
                    let s = w * (x[u] + x[v]);
                    y[u] += s;
                    y[v] += s;
                }
            }
        }
    }
}

// =============================================================================
// AggregatedEdges: Deduplicated edges with summed weights
// =============================================================================

/// Aggregated unique edges with summed weights.
///
/// Storage: Three parallel arrays (u, v, weight) of length P_unique.
/// Matvec cost: O(P_unique) which can be << O(n_obs) when many (FE_p, FE_q)
/// pairs appear multiple times in the data.
///
/// Node numbering:
/// - Left nodes (FE p): 0..kp (stored directly in u)
/// - Right nodes (FE q): kp..kp+kq (stored with shift in v)
pub struct AggregatedEdges {
    /// Left node indices (in 0..kp)
    u: Vec<usize>,
    /// Right node indices (in kp..kp+kq, already shifted)
    v: Vec<usize>,
    /// Aggregated edge weights (sum of observation weights for this (u,v) pair)
    weights: Vec<f64>,
}

impl AggregatedEdges {
    /// Build aggregated edges from FE index arrays.
    ///
    /// # Algorithm
    /// 1. Compute edge keys: `key[i] = gp + kp * gq`
    /// 2. Sort indices by key (argsort)
    /// 3. Run-length encode: for each run of same key, sum weights
    /// 4. Store (u, v, total_weight) tuples
    ///
    /// # Arguments
    /// * `group_ids_p` - Group IDs for FE p (length n_obs, values in 0..kp)
    /// * `group_ids_q` - Group IDs for FE q (length n_obs, values in 0..kq)
    /// * `kp` - Number of groups in FE p
    /// * `weights` - Optional observation weights (length n_obs if Some)
    pub fn from_indices(
        group_ids_p: &[usize],
        group_ids_q: &[usize],
        kp: usize,
        weights: Option<&[f64]>,
    ) -> Self {
        let n_obs = group_ids_p.len();
        debug_assert_eq!(group_ids_q.len(), n_obs);
        if let Some(w) = weights {
            debug_assert_eq!(w.len(), n_obs);
        }

        if n_obs == 0 {
            return Self {
                u: Vec::new(),
                v: Vec::new(),
                weights: Vec::new(),
            };
        }

        // Step 1: Compute edge keys and create index array
        // Key = gp + kp * gq ensures unique key per (gp, gq) pair
        let mut indices: Vec<usize> = (0..n_obs).collect();

        // Step 2: Sort indices by edge key
        indices.sort_by_key(|&i| {
            let gp = group_ids_p[i];
            let gq = group_ids_q[i];
            gp as u64 + (kp as u64) * (gq as u64)
        });

        // Step 3: RLE - walk sorted indices, accumulate weights for runs
        let mut edge_u = Vec::with_capacity(n_obs / 2); // Estimate: at least 50% unique
        let mut edge_v = Vec::with_capacity(n_obs / 2);
        let mut edge_w = Vec::with_capacity(n_obs / 2);

        let mut i = 0;
        while i < n_obs {
            let first_idx = indices[i];
            let gp = group_ids_p[first_idx];
            let gq = group_ids_q[first_idx];
            let mut w_sum = weights.map(|ws| ws[first_idx]).unwrap_or(1.0);

            // Accumulate run of same (gp, gq)
            let mut j = i + 1;
            while j < n_obs {
                let idx = indices[j];
                if group_ids_p[idx] != gp || group_ids_q[idx] != gq {
                    break;
                }
                w_sum += weights.map(|ws| ws[idx]).unwrap_or(1.0);
                j += 1;
            }

            edge_u.push(gp);
            edge_v.push(kp + gq); // Shift to right-node space
            edge_w.push(w_sum);

            i = j;
        }

        // Shrink to fit
        edge_u.shrink_to_fit();
        edge_v.shrink_to_fit();
        edge_w.shrink_to_fit();

        Self {
            u: edge_u,
            v: edge_v,
            weights: edge_w,
        }
    }

    /// Number of unique edges.
    #[inline]
    pub fn num_unique(&self) -> usize {
        self.u.len()
    }

    /// Get compression ratio: n_obs / num_unique.
    ///
    /// Higher ratio means more duplicate edges were merged.
    #[inline]
    pub fn compression_ratio(&self, n_obs: usize) -> f64 {
        if self.u.is_empty() {
            1.0
        } else {
            n_obs as f64 / self.u.len() as f64
        }
    }
}

impl EdgeSource for AggregatedEdges {
    #[inline]
    fn num_edges(&self) -> usize {
        self.u.len()
    }

    fn laplacian_matvec(&self, x: &[f64], y: &mut [f64]) {
        for ((&u, &v), &w) in self.u.iter().zip(self.v.iter()).zip(self.weights.iter()) {
            let d = w * (x[u] - x[v]);
            y[u] += d;
            y[v] -= d;
        }
    }

    fn gram_matvec(&self, x: &[f64], y: &mut [f64]) {
        // G = Deg + Adj: y[u] += w*(x[u] + x[v]), y[v] += w*(x[u] + x[v])
        for ((&u, &v), &w) in self.u.iter().zip(self.v.iter()).zip(self.weights.iter()) {
            let s = w * (x[u] + x[v]);
            y[u] += s;
            y[v] += s;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_streaming_edges_unweighted() {
        // Simple bipartite graph: 2 left nodes, 2 right nodes
        // Edges: (0,0), (0,1), (1,1)
        let group_ids_p = vec![0, 0, 1];
        let group_ids_q = vec![0, 1, 1];
        let kp = 2;

        let edges = StreamingEdges::new(group_ids_p, group_ids_q, kp, None);
        assert_eq!(edges.num_edges(), 3);

        // Test matvec: x = [1, 2, 3, 4] for nodes [p0, p1, q0, q1]
        let x = vec![1.0, 2.0, 3.0, 4.0];
        let mut y = vec![0.0; 4];
        edges.laplacian_matvec(&x, &mut y);

        // Edge (0,2): d = 1-3 = -2, y[0] += -2, y[2] -= -2
        // Edge (0,3): d = 1-4 = -3, y[0] += -3, y[3] -= -3
        // Edge (1,3): d = 2-4 = -2, y[1] += -2, y[3] -= -2
        // y = [-2-3, -2, 2, 3+2] = [-5, -2, 2, 5]
        assert_eq!(y, vec![-5.0, -2.0, 2.0, 5.0]);
    }

    #[test]
    fn test_streaming_edges_weighted() {
        let group_ids_p = vec![0, 0];
        let group_ids_q = vec![0, 0];
        let kp = 1;
        let weights = Some(vec![2.0, 3.0]);

        let edges = StreamingEdges::new(group_ids_p, group_ids_q, kp, weights);

        let x = vec![1.0, 5.0];
        let mut y = vec![0.0; 2];
        edges.laplacian_matvec(&x, &mut y);

        // Edge (0,1) with weight 2: d = 2*(1-5) = -8
        // Edge (0,1) with weight 3: d = 3*(1-5) = -12
        // Total: y[0] = -20, y[1] = 20
        assert_eq!(y, vec![-20.0, 20.0]);
    }

    #[test]
    fn test_aggregated_edges_basic() {
        // Same edges as streaming test but with duplicates
        let group_ids_p = vec![0, 0, 1, 0]; // Extra (0,0) edge
        let group_ids_q = vec![0, 1, 1, 0];
        let kp = 2;

        let edges = AggregatedEdges::from_indices(&group_ids_p, &group_ids_q, kp, None);

        // Should have 3 unique edges: (0,0), (0,1), (1,1)
        assert_eq!(edges.num_unique(), 3);

        // The (0,0) edge should have weight 2 (appears twice)
        // Find (0,0) edge
        let idx = edges
            .u
            .iter()
            .zip(edges.v.iter())
            .position(|(&u, &v)| u == 0 && v == 2)
            .unwrap();
        assert_eq!(edges.weights[idx], 2.0);
    }

    #[test]
    fn test_aggregated_edges_weighted() {
        let group_ids_p = vec![0, 0]; // Two (0,0) edges
        let group_ids_q = vec![0, 0];
        let kp = 1;
        let weights = vec![2.0, 3.0];

        let edges = AggregatedEdges::from_indices(&group_ids_p, &group_ids_q, kp, Some(&weights));

        // Should have 1 unique edge with weight 5
        assert_eq!(edges.num_unique(), 1);
        assert_eq!(edges.weights[0], 5.0);

        // Verify matvec gives same result as streaming
        let x = vec![1.0, 5.0];
        let mut y = vec![0.0; 2];
        edges.laplacian_matvec(&x, &mut y);
        // Single edge with weight 5: d = 5*(1-5) = -20
        assert_eq!(y, vec![-20.0, 20.0]);
    }

    #[test]
    fn test_aggregated_edges_empty() {
        let edges = AggregatedEdges::from_indices(&[], &[], 0, None);
        assert_eq!(edges.num_unique(), 0);
    }

    #[test]
    fn test_streaming_vs_aggregated_equivalence() {
        // Create data with some duplicate edges
        let group_ids_p = vec![0, 1, 0, 2, 1, 0];
        let group_ids_q = vec![0, 1, 0, 2, 1, 1];
        let kp = 3;
        let kq = 3;
        let v = kp + kq;

        let streaming = StreamingEdges::new(
            group_ids_p.clone(),
            group_ids_q.clone(),
            kp,
            None,
        );
        let aggregated = AggregatedEdges::from_indices(&group_ids_p, &group_ids_q, kp, None);

        // Test that matvec produces same result
        let x: Vec<f64> = (0..v).map(|i| i as f64 + 1.0).collect();

        let mut y_streaming = vec![0.0; v];
        let mut y_aggregated = vec![0.0; v];

        streaming.laplacian_matvec(&x, &mut y_streaming);
        aggregated.laplacian_matvec(&x, &mut y_aggregated);

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
