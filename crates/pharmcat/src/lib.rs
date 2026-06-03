//! Rust port of PharmCAT.
//!
//! The Java implementation in `repos/PharmCAT` is the behavioral reference
//! while this crate is ported module by module.

pub mod cli;
pub mod common;
pub mod definition;
pub mod matcher;
pub mod phenotype;
pub mod pipeline;
pub mod report;
pub mod vcf;
