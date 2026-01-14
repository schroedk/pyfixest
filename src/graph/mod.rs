//! Graph algorithms for the Laplacian preconditioner.
//!
//! This module provides:
//! - [`UnionFind`]: Disjoint Set Union for connected component detection
//! - [`ComponentsCSR`]: CSR representation of connected components

pub mod components;
pub mod union_find;

pub use components::ComponentsCSR;
pub use union_find::UnionFind;
