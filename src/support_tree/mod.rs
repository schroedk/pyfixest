//! Support tree (spanning forest) structures for Laplacian preconditioning.
//!
//! This module provides:
//! - [`SpanningForest`]: Spanning forest construction via DSU
//! - [`TreeSolver`]: Linear-time tree Laplacian solve
//! - [`FixedPCG`]: Fixed-iteration PCG with tree preconditioner

pub mod forest_builder;
pub mod pcg;
pub mod solver;

pub use forest_builder::SpanningForest;
pub use pcg::FixedPCG;
pub use solver::{TreeSolver, TreeSolverOwned};
