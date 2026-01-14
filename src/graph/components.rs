//! Connected components in CSR (Compressed Sparse Row) format.
//!
//! Provides efficient iteration over nodes within each component,
//! enabling O(V) gauge projection (component-wise mean subtraction).

use super::union_find::UnionFind;

/// CSR-like representation of connected components.
///
/// Enables O(1) lookup of component ID for any node, and efficient
/// iteration over all nodes in a component (for mean subtraction).
///
/// # Memory Layout
/// - `comp_id`: length V, maps node -> component ID
/// - `offsets`: length num_components + 1, CSR offsets
/// - `nodes`: length V, node IDs grouped by component
/// - `roots`: length num_components, one root per component
#[derive(Debug, Clone)]
pub struct ComponentsCSR {
    /// Component ID for each node: `comp_id[u]` gives the component index (0..num_components)
    pub comp_id: Vec<usize>,
    /// CSR offsets: nodes in component `cid` are `nodes[offsets[cid]..offsets[cid+1]]`
    pub offsets: Vec<usize>,
    /// Node IDs grouped by component
    pub nodes: Vec<usize>,
    /// Number of connected components
    pub num_components: usize,
    /// Root node for each component (the DSU representative)
    pub roots: Vec<usize>,
}

impl ComponentsCSR {
    /// Build component CSR from a UnionFind structure.
    ///
    /// # Arguments
    /// * `uf` - A mutable UnionFind (find() may compress paths)
    /// * `n` - Number of nodes (must match uf.len())
    ///
    /// # Algorithm
    /// 1. Find all roots and assign component IDs
    /// 2. Count nodes per component
    /// 3. Build CSR offsets
    /// 4. Fill nodes array by component
    ///
    /// # Complexity
    /// O(n * α(n)) where α is inverse Ackermann (effectively O(n))
    pub fn from_union_find(uf: &mut UnionFind, n: usize) -> Self {
        debug_assert_eq!(uf.len(), n);

        if n == 0 {
            return Self {
                comp_id: Vec::new(),
                offsets: vec![0],
                nodes: Vec::new(),
                num_components: 0,
                roots: Vec::new(),
            };
        }

        // Step 1: Find all roots and assign component IDs
        // root_to_cid maps root node -> component ID
        let mut root_to_cid = vec![usize::MAX; n];
        let mut roots = Vec::new();
        let mut num_components = 0;

        for u in 0..n {
            let root = uf.find(u);
            if root_to_cid[root] == usize::MAX {
                root_to_cid[root] = num_components;
                roots.push(root);
                num_components += 1;
            }
        }

        // Step 2: Assign component ID to each node
        let mut comp_id = vec![0; n];
        for u in 0..n {
            let root = uf.find(u);
            comp_id[u] = root_to_cid[root];
        }

        // Step 3: Count nodes per component
        let mut counts = vec![0usize; num_components];
        for &cid in &comp_id {
            counts[cid] += 1;
        }

        // Step 4: Build CSR offsets (exclusive prefix sum)
        let mut offsets = vec![0; num_components + 1];
        for (i, &count) in counts.iter().enumerate() {
            offsets[i + 1] = offsets[i] + count;
        }

        // Step 5: Fill nodes array
        let mut nodes = vec![0; n];
        let mut write_pos = offsets[..num_components].to_vec();
        for u in 0..n {
            let cid = comp_id[u];
            nodes[write_pos[cid]] = u;
            write_pos[cid] += 1;
        }

        Self {
            comp_id,
            offsets,
            nodes,
            num_components,
            roots,
        }
    }

    /// Get the slice of nodes belonging to component `cid`.
    ///
    /// # Arguments
    /// * `cid` - Component ID (must be < num_components)
    ///
    /// # Returns
    /// Slice of node indices in this component.
    #[inline]
    pub fn component_nodes(&self, cid: usize) -> &[usize] {
        debug_assert!(cid < self.num_components);
        &self.nodes[self.offsets[cid]..self.offsets[cid + 1]]
    }

    /// Get the root node for component `cid`.
    ///
    /// # Arguments
    /// * `cid` - Component ID (must be < num_components)
    #[inline]
    pub fn component_root(&self, cid: usize) -> usize {
        debug_assert!(cid < self.num_components);
        self.roots[cid]
    }

    /// Get the size (number of nodes) of component `cid`.
    #[inline]
    pub fn component_size(&self, cid: usize) -> usize {
        debug_assert!(cid < self.num_components);
        self.offsets[cid + 1] - self.offsets[cid]
    }

    /// Get the component ID for node `u`.
    #[inline]
    pub fn component_of(&self, u: usize) -> usize {
        self.comp_id[u]
    }

    /// Get the total number of nodes.
    #[inline]
    pub fn num_nodes(&self) -> usize {
        self.nodes.len()
    }

    /// Project a vector into Range(L) by subtracting component-wise means.
    ///
    /// For each component C, computes mean = sum(r[u] for u in C) / |C|,
    /// then subtracts mean from each r[u] in C.
    ///
    /// This ensures r is orthogonal to the nullspace of the Laplacian
    /// (which consists of vectors constant on each component).
    ///
    /// # Arguments
    /// * `r` - Vector to project (length must equal num_nodes)
    #[inline]
    pub fn project_to_range(&self, r: &mut [f64]) {
        debug_assert_eq!(r.len(), self.num_nodes());

        for cid in 0..self.num_components {
            let nodes = self.component_nodes(cid);
            let n = nodes.len() as f64;

            // Compute mean
            let sum: f64 = nodes.iter().map(|&u| r[u]).sum();
            let mean = sum / n;

            // Subtract mean
            for &u in nodes {
                r[u] -= mean;
            }
        }
    }

    /// Compute weighted component-wise mean and subtract.
    ///
    /// For each component C:
    /// - weighted_mean = sum(w[u] * r[u] for u in C) / sum(w[u] for u in C)
    /// - r[u] -= weighted_mean for all u in C
    ///
    /// # Arguments
    /// * `r` - Vector to project
    /// * `weights` - Node weights (length must equal num_nodes)
    #[inline]
    pub fn project_to_range_weighted(&self, r: &mut [f64], weights: &[f64]) {
        debug_assert_eq!(r.len(), self.num_nodes());
        debug_assert_eq!(weights.len(), self.num_nodes());

        for cid in 0..self.num_components {
            let nodes = self.component_nodes(cid);

            // Compute weighted mean
            let mut weighted_sum = 0.0;
            let mut weight_sum = 0.0;
            for &u in nodes {
                weighted_sum += weights[u] * r[u];
                weight_sum += weights[u];
            }
            let mean = weighted_sum / weight_sum;

            // Subtract mean
            for &u in nodes {
                r[u] -= mean;
            }
        }
    }

    /// Project a vector into Range(G) for a bipartite Gram matrix.
    ///
    /// For a bipartite graph with P (nodes 0..kp) and Q (nodes kp..v):
    /// - Null(G) = span of alternating sign vectors (+1 on P, -1 on Q) per component
    /// - Range(G) = {x : sum_P(x) = sum_Q(x)} per component
    ///
    /// To project to Range(G), we adjust so that sum_P = sum_Q for each component:
    /// - Subtract (sum_P - sum_Q) / |C| from P nodes in component
    /// - Add (sum_P - sum_Q) / |C| to Q nodes in component
    ///
    /// # Arguments
    /// * `r` - Vector to project (length must equal num_nodes)
    /// * `kp` - Number of P nodes (first kp nodes are P, rest are Q)
    #[inline]
    pub fn project_to_gram_range(&self, r: &mut [f64], kp: usize) {
        debug_assert_eq!(r.len(), self.num_nodes());

        for cid in 0..self.num_components {
            let nodes = self.component_nodes(cid);

            // Compute sum_P and sum_Q for this component
            let mut sum_p = 0.0;
            let mut sum_q = 0.0;
            let mut count_p = 0;
            let mut count_q = 0;

            for &u in nodes {
                if u < kp {
                    sum_p += r[u];
                    count_p += 1;
                } else {
                    sum_q += r[u];
                    count_q += 1;
                }
            }

            // Skip if component has only P or only Q nodes (shouldn't happen in valid bipartite graph)
            let component_size = count_p + count_q;
            if count_p == 0 || count_q == 0 || component_size == 0 {
                continue;
            }

            // Correction to make sum_P = sum_Q
            let diff = sum_p - sum_q;
            let correction = diff / (component_size as f64);

            // Apply correction: subtract from P, add to Q
            for &u in nodes {
                if u < kp {
                    r[u] -= correction;
                } else {
                    r[u] += correction;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_single_component() {
        let mut uf = UnionFind::new(5);
        uf.union(0, 1);
        uf.union(1, 2);
        uf.union(2, 3);
        uf.union(3, 4);

        let csr = ComponentsCSR::from_union_find(&mut uf, 5);

        assert_eq!(csr.num_components, 1);
        assert_eq!(csr.component_size(0), 5);

        // All nodes should be in component 0
        for u in 0..5 {
            assert_eq!(csr.component_of(u), 0);
        }

        // component_nodes should contain all 5 nodes
        let nodes = csr.component_nodes(0);
        assert_eq!(nodes.len(), 5);
        let mut sorted_nodes = nodes.to_vec();
        sorted_nodes.sort();
        assert_eq!(sorted_nodes, vec![0, 1, 2, 3, 4]);
    }

    #[test]
    fn test_multiple_components() {
        let mut uf = UnionFind::new(6);
        // Component 1: {0, 1, 2}
        uf.union(0, 1);
        uf.union(1, 2);
        // Component 2: {3, 4}
        uf.union(3, 4);
        // Component 3: {5} (singleton)

        let csr = ComponentsCSR::from_union_find(&mut uf, 6);

        assert_eq!(csr.num_components, 3);

        // Check that nodes in same union-find component have same comp_id
        assert_eq!(csr.component_of(0), csr.component_of(1));
        assert_eq!(csr.component_of(1), csr.component_of(2));
        assert_eq!(csr.component_of(3), csr.component_of(4));

        // Check different components have different comp_id
        assert_ne!(csr.component_of(0), csr.component_of(3));
        assert_ne!(csr.component_of(0), csr.component_of(5));
        assert_ne!(csr.component_of(3), csr.component_of(5));

        // Check component sizes
        let cid_012 = csr.component_of(0);
        let cid_34 = csr.component_of(3);
        let cid_5 = csr.component_of(5);

        assert_eq!(csr.component_size(cid_012), 3);
        assert_eq!(csr.component_size(cid_34), 2);
        assert_eq!(csr.component_size(cid_5), 1);
    }

    #[test]
    fn test_empty() {
        let mut uf = UnionFind::new(0);
        let csr = ComponentsCSR::from_union_find(&mut uf, 0);

        assert_eq!(csr.num_components, 0);
        assert_eq!(csr.num_nodes(), 0);
    }

    #[test]
    fn test_all_singletons() {
        let mut uf = UnionFind::new(4);
        // No unions - each node is its own component

        let csr = ComponentsCSR::from_union_find(&mut uf, 4);

        assert_eq!(csr.num_components, 4);
        for cid in 0..4 {
            assert_eq!(csr.component_size(cid), 1);
        }
    }

    #[test]
    fn test_project_to_range() {
        let mut uf = UnionFind::new(4);
        uf.union(0, 1); // Component A: {0, 1}
        uf.union(2, 3); // Component B: {2, 3}

        let csr = ComponentsCSR::from_union_find(&mut uf, 4);

        // r = [1, 3, 10, 20]
        // Component A mean = (1+3)/2 = 2, so [1-2, 3-2] = [-1, 1]
        // Component B mean = (10+20)/2 = 15, so [10-15, 20-15] = [-5, 5]
        let mut r = vec![1.0, 3.0, 10.0, 20.0];
        csr.project_to_range(&mut r);

        let cid_01 = csr.component_of(0);
        let cid_23 = csr.component_of(2);

        // Check nodes in component {0,1}
        if cid_01 == csr.component_of(0) {
            assert!((r[0] - (-1.0)).abs() < 1e-10);
            assert!((r[1] - 1.0).abs() < 1e-10);
        }

        // Check nodes in component {2,3}
        if cid_23 == csr.component_of(2) {
            assert!((r[2] - (-5.0)).abs() < 1e-10);
            assert!((r[3] - 5.0).abs() < 1e-10);
        }

        // Verify each component sums to zero
        let sum_01 = r[0] + r[1];
        let sum_23 = r[2] + r[3];
        assert!(sum_01.abs() < 1e-10);
        assert!(sum_23.abs() < 1e-10);
    }

    #[test]
    fn test_project_to_range_weighted() {
        let mut uf = UnionFind::new(3);
        uf.union(0, 1);
        uf.union(1, 2);

        let csr = ComponentsCSR::from_union_find(&mut uf, 3);

        // r = [0, 10, 20], weights = [1, 2, 2]
        // weighted_mean = (1*0 + 2*10 + 2*20) / (1+2+2) = 60/5 = 12
        // r_projected = [0-12, 10-12, 20-12] = [-12, -2, 8]
        let mut r = vec![0.0, 10.0, 20.0];
        let weights = vec![1.0, 2.0, 2.0];
        csr.project_to_range_weighted(&mut r, &weights);

        assert!((r[0] - (-12.0)).abs() < 1e-10);
        assert!((r[1] - (-2.0)).abs() < 1e-10);
        assert!((r[2] - 8.0).abs() < 1e-10);

        // Verify weighted sum is zero
        let weighted_sum: f64 = r.iter().zip(weights.iter()).map(|(&ri, &wi)| ri * wi).sum();
        assert!(weighted_sum.abs() < 1e-10);
    }

    #[test]
    fn test_component_nodes_contains_correct_nodes() {
        let mut uf = UnionFind::new(5);
        uf.union(0, 2);
        uf.union(2, 4); // Component: {0, 2, 4}
        uf.union(1, 3); // Component: {1, 3}

        let csr = ComponentsCSR::from_union_find(&mut uf, 5);

        // Find which component has {0, 2, 4}
        let cid_024 = csr.component_of(0);
        let nodes_024: Vec<usize> = csr.component_nodes(cid_024).to_vec();
        let mut sorted = nodes_024.clone();
        sorted.sort();
        assert_eq!(sorted, vec![0, 2, 4]);

        // Find which component has {1, 3}
        let cid_13 = csr.component_of(1);
        let nodes_13: Vec<usize> = csr.component_nodes(cid_13).to_vec();
        let mut sorted = nodes_13.clone();
        sorted.sort();
        assert_eq!(sorted, vec![1, 3]);
    }
}
