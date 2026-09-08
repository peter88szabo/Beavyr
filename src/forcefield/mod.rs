//! Native force fields, for geometry work that must not wait on an external program.
//!
//! Beavyr can drive xTB and Behemoth, but both are separate processes whose startup cost is
//! measured in seconds. Cleaning up a structure while building it has to feel instantaneous, so
//! that job is done in-process here.
//!
//! Currently one force field: [`dreiding`].

pub mod dreiding;
