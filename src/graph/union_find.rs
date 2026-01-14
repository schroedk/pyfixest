//! Disjoint Set Union (Union-Find) data structure.
//!
//! Provides efficient union and find operations for connected component detection.
//! Uses path compression and union by rank for near-O(1) amortized operations.

/// Disjoint Set Union (DSU) for connected component detection.
///
/// Uses path compression and union by rank for near-O(α(n)) amortized operations,
/// where α is the inverse Ackermann function (effectively constant for practical inputs).
#[derive(Debug, Clone)]
pub struct UnionFind {
    /// Parent pointer for each node. `parent[u] == u` means u is a root.
    parent: Vec<usize>,
    /// Rank (upper bound on tree height) for union by rank.
    rank: Vec<usize>,
}

impl UnionFind {
    /// Create a new DSU with `n` nodes, each in its own component.
    ///
    /// # Arguments
    /// * `n` - Number of nodes (0..n)
    ///
    /// # Example
    /// ```
    /// use pyfixest::graph::union_find::UnionFind;
    /// let mut uf = UnionFind::new(5);
    /// assert!(!uf.same(0, 1));
    /// uf.union(0, 1);
    /// assert!(uf.same(0, 1));
    /// ```
    #[inline]
    pub fn new(n: usize) -> Self {
        Self {
            parent: (0..n).collect(),
            rank: vec![0; n],
        }
    }

    /// Find the representative (root) of the component containing `u`.
    ///
    /// Uses path compression: all nodes on the path to root are updated
    /// to point directly to the root.
    ///
    /// # Arguments
    /// * `u` - Node index (must be < n)
    ///
    /// # Returns
    /// The root node of the component containing `u`.
    #[inline]
    pub fn find(&mut self, u: usize) -> usize {
        if self.parent[u] != u {
            self.parent[u] = self.find(self.parent[u]);
        }
        self.parent[u]
    }

    /// Union the components containing `u` and `v`.
    ///
    /// Uses union by rank to keep trees balanced.
    ///
    /// # Arguments
    /// * `u` - First node index
    /// * `v` - Second node index
    ///
    /// # Returns
    /// * `true` if they were in different components (union performed)
    /// * `false` if already in the same component (no change)
    #[inline]
    pub fn union(&mut self, u: usize, v: usize) -> bool {
        let root_u = self.find(u);
        let root_v = self.find(v);

        if root_u == root_v {
            return false;
        }

        // Union by rank: attach smaller tree under larger tree
        match self.rank[root_u].cmp(&self.rank[root_v]) {
            std::cmp::Ordering::Less => {
                self.parent[root_u] = root_v;
            }
            std::cmp::Ordering::Greater => {
                self.parent[root_v] = root_u;
            }
            std::cmp::Ordering::Equal => {
                self.parent[root_v] = root_u;
                self.rank[root_u] += 1;
            }
        }

        true
    }

    /// Check if `u` and `v` are in the same component.
    ///
    /// # Arguments
    /// * `u` - First node index
    /// * `v` - Second node index
    ///
    /// # Returns
    /// `true` if both nodes are in the same component.
    #[inline]
    pub fn same(&mut self, u: usize, v: usize) -> bool {
        self.find(u) == self.find(v)
    }

    /// Get the number of nodes in this DSU.
    #[inline]
    pub fn len(&self) -> usize {
        self.parent.len()
    }

    /// Check if this DSU is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.parent.is_empty()
    }

    /// Count the number of connected components.
    ///
    /// Note: This requires iterating over all nodes and is O(n).
    pub fn count_components(&mut self) -> usize {
        let n = self.parent.len();
        (0..n).filter(|&u| self.find(u) == u).count()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let uf = UnionFind::new(5);
        assert_eq!(uf.len(), 5);
        assert!(!uf.is_empty());
    }

    #[test]
    fn test_find_initial() {
        let mut uf = UnionFind::new(5);
        for i in 0..5 {
            assert_eq!(uf.find(i), i);
        }
    }

    #[test]
    fn test_union_and_same() {
        let mut uf = UnionFind::new(5);

        // Initially all separate
        assert!(!uf.same(0, 1));
        assert!(!uf.same(1, 2));

        // Union 0 and 1
        assert!(uf.union(0, 1)); // Should return true (different components)
        assert!(uf.same(0, 1));
        assert!(!uf.same(0, 2));

        // Union 1 and 2 (transitively connects 0, 1, 2)
        assert!(uf.union(1, 2));
        assert!(uf.same(0, 2));

        // Union already-connected nodes
        assert!(!uf.union(0, 2)); // Should return false (same component)
    }

    #[test]
    fn test_count_components() {
        let mut uf = UnionFind::new(6);
        assert_eq!(uf.count_components(), 6);

        uf.union(0, 1);
        assert_eq!(uf.count_components(), 5);

        uf.union(2, 3);
        assert_eq!(uf.count_components(), 4);

        uf.union(0, 2); // Merges {0,1} with {2,3}
        assert_eq!(uf.count_components(), 3);
    }

    #[test]
    fn test_path_compression() {
        let mut uf = UnionFind::new(5);

        // Create a chain: 0 -> 1 -> 2 -> 3 -> 4
        uf.union(3, 4);
        uf.union(2, 3);
        uf.union(1, 2);
        uf.union(0, 1);

        // All should be in the same component
        let root = uf.find(0);
        assert_eq!(uf.find(1), root);
        assert_eq!(uf.find(2), root);
        assert_eq!(uf.find(3), root);
        assert_eq!(uf.find(4), root);
    }

    #[test]
    fn test_empty() {
        let uf = UnionFind::new(0);
        assert!(uf.is_empty());
        assert_eq!(uf.len(), 0);
    }
}
