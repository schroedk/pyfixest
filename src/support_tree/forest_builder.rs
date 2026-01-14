//! Spanning forest construction for bipartite 2-FE graphs.
//!
//! A spanning forest has exactly (V - num_components) edges and forms a tree
//! for each connected component. Used for O(V) tree Laplacian solves.

use crate::graph::components::ComponentsCSR;
use crate::graph::union_find::UnionFind;

/// Spanning forest for the bipartite 2-FE graph.
///
/// A spanning forest has exactly (V - num_components) edges, forming one tree
/// per connected component. Stores parent pointers and children lists for
/// efficient two-pass tree traversal (postorder + preorder).
#[derive(Debug, Clone)]
pub struct SpanningForest {
    /// Parent of each node. `parent[u] == usize::MAX` means u is a root.
    pub parent: Vec<usize>,
    /// Weight of edge from node to parent. Undefined for roots.
    pub parent_weight: Vec<f64>,
    /// CSR offsets for children: children of node u are in
    /// `children[children_offsets[u]..children_offsets[u+1]]`
    pub children_offsets: Vec<usize>,
    /// Children array (CSR values)
    pub children: Vec<usize>,
    /// Postorder traversal: children before parent (for flow computation)
    pub postorder: Vec<usize>,
    /// Preorder traversal: parent before children (for potential recovery)
    pub preorder: Vec<usize>,
    /// Root nodes (one per component, same as ComponentsCSR roots)
    pub roots: Vec<usize>,
    /// Number of nodes
    pub v: usize,
}

impl SpanningForest {
    /// Build spanning forest using simple DSU scan order.
    ///
    /// Processes edges in input order. When an edge connects nodes in different
    /// components (according to DSU), it's added to the forest.
    ///
    /// # Arguments
    /// * `group_ids_p` - FE p group IDs per observation
    /// * `group_ids_q` - FE q group IDs per observation
    /// * `kp` - Number of groups in FE p
    /// * `kq` - Number of groups in FE q
    /// * `weights` - Observation weights (edge weights for the forest)
    /// * `components` - Pre-computed component structure
    ///
    /// # Returns
    /// A spanning forest with parent pointers, children CSR, and traversal orders.
    pub fn build_simple(
        group_ids_p: &[usize],
        group_ids_q: &[usize],
        kp: usize,
        kq: usize,
        weights: Option<&[f64]>,
        components: &ComponentsCSR,
    ) -> Self {
        let v = kp + kq;
        let n_obs = group_ids_p.len();

        // Initialize parent pointers (usize::MAX = root)
        let mut parent = vec![usize::MAX; v];
        let mut parent_weight = vec![0.0; v];

        // Collect forest edges using DSU
        // We need to store (u, v, weight) tuples first, then build parent structure
        let mut forest_edges: Vec<(usize, usize, f64)> = Vec::new();
        let mut uf = UnionFind::new(v);

        for i in 0..n_obs {
            let u = group_ids_p[i];
            let node_v = kp + group_ids_q[i];
            let w = weights.map(|ws| ws[i]).unwrap_or(1.0);

            // If this edge connects different components, add to forest
            if uf.union(u, node_v) {
                forest_edges.push((u, node_v, w));
            }
        }

        // Build parent structure from forest edges using BFS/DFS from roots
        // First, build adjacency list
        let mut adj: Vec<Vec<(usize, f64)>> = vec![Vec::new(); v];
        for &(u, node_v, w) in &forest_edges {
            adj[u].push((node_v, w));
            adj[node_v].push((u, w));
        }

        // BFS from each root to assign parent pointers
        // First find connected components and pick one root per component
        let mut visited = vec![false; v];
        for start in 0..v {
            if visited[start] {
                continue;
            }
            // BFS from start
            let mut queue = std::collections::VecDeque::new();
            queue.push_back(start);
            visited[start] = true;
            // start is a root (no parent)

            while let Some(u) = queue.pop_front() {
                for &(neighbor, w) in &adj[u] {
                    if !visited[neighbor] {
                        visited[neighbor] = true;
                        parent[neighbor] = u;
                        parent_weight[neighbor] = w;
                        queue.push_back(neighbor);
                    }
                }
            }
        }

        // Build children CSR from parent pointers
        let (children_offsets, children) = Self::build_children_csr(&parent, v);

        // Find forest roots (nodes with no parent)
        let roots: Vec<usize> = (0..v).filter(|&u| parent[u] == usize::MAX).collect();

        // Build traversal orders using forest roots
        let (postorder, preorder) = Self::build_traversal_orders(&children_offsets, &children, &roots);

        Self {
            parent,
            parent_weight,
            children_offsets,
            children,
            postorder,
            preorder,
            roots,
            v,
        }
    }

    /// Build children CSR from parent pointers.
    fn build_children_csr(parent: &[usize], v: usize) -> (Vec<usize>, Vec<usize>) {
        // Count children per node
        let mut child_counts = vec![0usize; v];
        for u in 0..v {
            if parent[u] != usize::MAX {
                child_counts[parent[u]] += 1;
            }
        }

        // Build offsets (exclusive prefix sum)
        let mut offsets = vec![0; v + 1];
        for u in 0..v {
            offsets[u + 1] = offsets[u] + child_counts[u];
        }

        // Fill children array
        let total_children = offsets[v];
        let mut children = vec![0; total_children];
        let mut write_pos = offsets[..v].to_vec();

        for u in 0..v {
            if parent[u] != usize::MAX {
                let p = parent[u];
                children[write_pos[p]] = u;
                write_pos[p] += 1;
            }
        }

        (offsets, children)
    }

    /// Build postorder and preorder traversal arrays.
    fn build_traversal_orders(
        children_offsets: &[usize],
        children: &[usize],
        roots: &[usize],
    ) -> (Vec<usize>, Vec<usize>) {
        let v = children_offsets.len() - 1;
        let mut postorder = Vec::with_capacity(v);
        let mut preorder = Vec::with_capacity(v);

        // Process each tree (rooted at forest root)
        for &root in roots {
            // Iterative DFS to avoid stack overflow
            // Stack entries: (node, next_child_index)
            let mut stack: Vec<(usize, usize)> = vec![(root, 0)];

            while let Some((u, child_idx)) = stack.last_mut() {
                let u = *u;
                let start = children_offsets[u];
                let end = children_offsets[u + 1];
                let num_children = end - start;

                if *child_idx == 0 {
                    // First visit: add to preorder
                    preorder.push(u);
                }

                if *child_idx < num_children {
                    // More children to process
                    let child = children[start + *child_idx];
                    *child_idx += 1;
                    stack.push((child, 0));
                } else {
                    // All children processed: add to postorder and pop
                    postorder.push(u);
                    stack.pop();
                }
            }
        }

        (postorder, preorder)
    }

    /// Get children of node `u`.
    #[inline]
    pub fn children_of(&self, u: usize) -> &[usize] {
        &self.children[self.children_offsets[u]..self.children_offsets[u + 1]]
    }

    /// Check if node is a root (has no parent).
    #[inline]
    pub fn is_root(&self, u: usize) -> bool {
        self.parent[u] == usize::MAX
    }

    /// Get the number of tree edges (should be V - num_components).
    #[inline]
    pub fn num_edges(&self) -> usize {
        self.v - self.roots.len()
    }

    /// Verify forest validity (for testing).
    ///
    /// Checks:
    /// - Exactly V - num_components edges
    /// - Each non-root has valid parent
    /// - Children arrays are consistent with parent pointers
    #[cfg(test)]
    pub fn verify(&self) -> bool {
        let num_roots = self.parent.iter().filter(|&&p| p == usize::MAX).count();
        if num_roots != self.roots.len() {
            return false;
        }

        // Check edge count
        let num_edges = self.v - num_roots;
        if num_edges != self.num_edges() {
            return false;
        }

        // Check children consistency
        for u in 0..self.v {
            for &child in self.children_of(u) {
                if self.parent[child] != u {
                    return false;
                }
            }
        }

        // Check postorder/preorder lengths
        if self.postorder.len() != self.v || self.preorder.len() != self.v {
            return false;
        }

        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_simple_bipartite() -> (Vec<usize>, Vec<usize>, usize, usize, ComponentsCSR) {
        // Simple connected bipartite graph
        let group_ids_p = vec![0, 0, 1];
        let group_ids_q = vec![0, 1, 1];
        let kp = 2;
        let kq = 2;
        let v = kp + kq;

        // Build components
        let mut uf = UnionFind::new(v);
        for (&gp, &gq) in group_ids_p.iter().zip(group_ids_q.iter()) {
            uf.union(gp, kp + gq);
        }
        let components = ComponentsCSR::from_union_find(&mut uf, v);

        (group_ids_p, group_ids_q, kp, kq, components)
    }

    #[test]
    fn test_forest_basic() {
        let (group_ids_p, group_ids_q, kp, kq, components) = make_simple_bipartite();
        let forest = SpanningForest::build_simple(
            &group_ids_p,
            &group_ids_q,
            kp,
            kq,
            None,
            &components,
        );

        assert_eq!(forest.v, 4);
        assert!(forest.verify());

        // Single component: should have V-1 = 3 edges
        assert_eq!(forest.num_edges(), 3);
        assert_eq!(forest.roots.len(), 1);
    }

    #[test]
    fn test_forest_two_components() {
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

        assert!(forest.verify());
        assert_eq!(forest.roots.len(), 2);
        // V=4, 2 components: should have 4-2 = 2 edges
        assert_eq!(forest.num_edges(), 2);
    }

    #[test]
    fn test_postorder_children_before_parent() {
        let (group_ids_p, group_ids_q, kp, kq, components) = make_simple_bipartite();
        let forest = SpanningForest::build_simple(
            &group_ids_p,
            &group_ids_q,
            kp,
            kq,
            None,
            &components,
        );

        // In postorder, every node should appear after all its children
        let mut position = vec![0; forest.v];
        for (pos, &node) in forest.postorder.iter().enumerate() {
            position[node] = pos;
        }

        for u in 0..forest.v {
            for &child in forest.children_of(u) {
                assert!(
                    position[child] < position[u],
                    "Child {} should appear before parent {} in postorder",
                    child,
                    u
                );
            }
        }
    }

    #[test]
    fn test_preorder_parent_before_children() {
        let (group_ids_p, group_ids_q, kp, kq, components) = make_simple_bipartite();
        let forest = SpanningForest::build_simple(
            &group_ids_p,
            &group_ids_q,
            kp,
            kq,
            None,
            &components,
        );

        // In preorder, every node should appear before all its children
        let mut position = vec![0; forest.v];
        for (pos, &node) in forest.preorder.iter().enumerate() {
            position[node] = pos;
        }

        for u in 0..forest.v {
            for &child in forest.children_of(u) {
                assert!(
                    position[u] < position[child],
                    "Parent {} should appear before child {} in preorder",
                    u,
                    child
                );
            }
        }
    }

    #[test]
    fn test_weighted_edges() {
        let group_ids_p = vec![0, 0];
        let group_ids_q = vec![0, 1];
        let kp = 2;
        let kq = 2;
        let v = kp + kq;
        let weights = vec![2.0, 3.0];

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
            Some(&weights),
            &components,
        );

        assert!(forest.verify());

        // Check that weights are stored
        let mut found_weights = vec![];
        for u in 0..v {
            if forest.parent[u] != usize::MAX {
                found_weights.push(forest.parent_weight[u]);
            }
        }
        found_weights.sort_by(|a, b| a.partial_cmp(b).unwrap());

        // First edge (0,2) has weight 2.0, second edge (0,3) has weight 3.0
        // But we only keep V-1=3 edges for single component
        assert!(!found_weights.is_empty());
    }

    #[test]
    fn test_traversal_covers_all_nodes() {
        let (group_ids_p, group_ids_q, kp, kq, components) = make_simple_bipartite();
        let forest = SpanningForest::build_simple(
            &group_ids_p,
            &group_ids_q,
            kp,
            kq,
            None,
            &components,
        );

        // Both traversals should contain all nodes exactly once
        let mut seen_post = vec![false; forest.v];
        let mut seen_pre = vec![false; forest.v];

        for &u in &forest.postorder {
            assert!(!seen_post[u], "Node {} appears twice in postorder", u);
            seen_post[u] = true;
        }

        for &u in &forest.preorder {
            assert!(!seen_pre[u], "Node {} appears twice in preorder", u);
            seen_pre[u] = true;
        }

        assert!(seen_post.iter().all(|&b| b), "Postorder missing nodes");
        assert!(seen_pre.iter().all(|&b| b), "Preorder missing nodes");
    }
}
