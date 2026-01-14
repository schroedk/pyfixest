//! Laplacian operator for bipartite 2-FE graphs.
//!
//! This module provides:
//! - [`EdgeSource`]: Trait for edge representations
//! - [`StreamingEdges`]: Direct iteration over FE indices (no extra storage)
//! - [`AggregatedEdges`]: Deduplicated edges with summed weights
//! - [`LaplacianOp`]: Full Laplacian operator with matvec and gauge projection

pub mod edges;
pub mod operator;

pub use edges::{AggregatedEdges, EdgeSource, StreamingEdges};
pub use operator::LaplacianOp;
