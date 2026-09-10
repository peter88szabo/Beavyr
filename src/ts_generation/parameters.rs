//! User-facing reaction definitions translated into the existing RDA options.

use anyhow::{bail, Result};

use crate::optimizer::rda::{RdaDistanceMode, RdaOptions};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Algorithm {
    #[default]
    Rda,
    PoorMansNeb,
}

impl Algorithm {
    pub fn label(self) -> &'static str {
        match self {
            Self::Rda => "RDA",
            Self::PoorMansNeb => "Poor Man's NEB",
        }
    }

    pub fn file_stem(self) -> &'static str {
        match self {
            Self::Rda => "rda",
            Self::PoorMansNeb => "poormans_neb",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ReactionKind {
    ActiveAtoms,
    #[default]
    Bond,
    Transfer,
    Angle,
    Dihedral,
    Custom,
}

impl ReactionKind {
    pub const ALL: [Self; 6] = [
        Self::ActiveAtoms,
        Self::Bond,
        Self::Transfer,
        Self::Angle,
        Self::Dihedral,
        Self::Custom,
    ];

    pub fn label(self) -> &'static str {
        match self {
            Self::ActiveAtoms => "Active atoms / Cartesian motion",
            Self::Bond => "Bond change (2 atoms)",
            Self::Transfer => "Atom transfer (3 atoms)",
            Self::Angle => "Angle change (3 atoms)",
            Self::Dihedral => "Dihedral change (4 atoms)",
            Self::Custom => "Custom reactive coordinates",
        }
    }

    pub fn hint(self) -> &'static str {
        match self {
            Self::ActiveAtoms => {
                "Enter the atoms whose Cartesian motion defines the reaction; blank uses all atoms."
            }
            Self::Bond => "Atom numbers i j: the bond being formed, broken, or stretched.",
            Self::Transfer => "Donor, transferred atom, acceptor: i j k defines bonds i-j and j-k.",
            Self::Angle => "Atom numbers i j k; j is the vertex.",
            Self::Dihedral => "Atom numbers i j k l; j-k is the central bond.",
            Self::Custom => {
                "One coordinate per line or separated by semicolons, e.g. bonds: 1 2; 2 3."
            }
        }
    }
}

#[derive(Clone)]
pub struct Parameters {
    pub algorithm: Algorithm,
    pub reaction: ReactionKind,
    pub indices: String,
    pub active_atoms: String,
    pub align_atoms: String,
    pub bonds: String,
    pub angles: String,
    pub dihedrals: String,
    pub weights: String,
    pub extra_bonds: String,
    pub extra_angles: String,
    pub extra_dihedrals: String,
    pub gamma_values: String,
    pub options: RdaOptions,
}

impl Default for Parameters {
    fn default() -> Self {
        Self {
            algorithm: Algorithm::Rda,
            reaction: ReactionKind::Bond,
            indices: "1 2".into(),
            active_atoms: String::new(),
            align_atoms: String::new(),
            bonds: String::new(),
            angles: String::new(),
            dihedrals: String::new(),
            weights: String::new(),
            extra_bonds: String::new(),
            extra_angles: String::new(),
            extra_dihedrals: String::new(),
            gamma_values: "0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9".into(),
            options: RdaOptions {
                distance_mode: RdaDistanceMode::ReactiveInternal,
                verbosity: 0,
                rda_trajectory_file: None,
                rda_chain_file: None,
                poormans_neb_file: None,
                poormans_neb_profile_file: None,
                ..Default::default()
            },
        }
    }
}

fn tokens(text: &str) -> impl Iterator<Item = &str> {
    text.split(|c: char| c.is_whitespace() || c == ',')
        .filter(|s| !s.is_empty())
}

/// UI atom numbers are one-based; the optimizer's indices are zero-based.
pub fn atom_indices(text: &str, natoms: usize) -> Result<Vec<usize>> {
    let mut out = Vec::new();
    for token in tokens(text) {
        let number = token.parse::<usize>().map_err(|_| {
            anyhow::anyhow!("Invalid atom number {token:?}; use integers starting at 1.")
        })?;
        if number == 0 || number > natoms {
            bail!("Atom {number} is outside the endpoint's atom range 1..={natoms}.");
        }
        let index = number - 1;
        if out.contains(&index) {
            bail!("Atom {number} is repeated in the same selection.");
        }
        out.push(index);
    }
    Ok(out)
}

fn coordinate<const N: usize>(text: &str, natoms: usize) -> Result<[usize; N]> {
    let indices = atom_indices(text, natoms)?;
    indices.try_into().map_err(|v: Vec<usize>| {
        anyhow::anyhow!(
            "This coordinate needs {N} atom numbers, but {} were entered.",
            v.len()
        )
    })
}

fn coordinates<const N: usize>(text: &str, natoms: usize) -> Result<Vec<[usize; N]>> {
    text.split([';', '\n'])
        .filter(|line| !line.trim().is_empty())
        .map(|line| coordinate(line, natoms))
        .collect()
}

fn numbers(text: &str, label: &str) -> Result<Vec<f64>> {
    tokens(text)
        .map(|token| {
            let value = token
                .parse::<f64>()
                .map_err(|_| anyhow::anyhow!("Invalid {label} value {token:?}."))?;
            if !value.is_finite() {
                bail!("{label} values must be finite.");
            }
            Ok(value)
        })
        .collect()
}

impl Parameters {
    pub fn rda_options(&self, natoms: usize) -> Result<RdaOptions> {
        let mut opts = self.options.clone();
        opts.poormans_neb = self.algorithm == Algorithm::PoorMansNeb;
        opts.reactive_bonds.clear();
        opts.reactive_angles.clear();
        opts.reactive_dihedrals.clear();
        match self.reaction {
            ReactionKind::ActiveAtoms => opts.distance_mode = RdaDistanceMode::Cartesian,
            ReactionKind::Bond => opts.reactive_bonds.push(coordinate(&self.indices, natoms)?),
            ReactionKind::Transfer => {
                let [i, j, k] = coordinate(&self.indices, natoms)?;
                opts.reactive_bonds = vec![[i, j], [j, k]];
            }
            ReactionKind::Angle => opts
                .reactive_angles
                .push(coordinate(&self.indices, natoms)?),
            ReactionKind::Dihedral => opts
                .reactive_dihedrals
                .push(coordinate(&self.indices, natoms)?),
            ReactionKind::Custom => {
                opts.reactive_bonds = coordinates(&self.bonds, natoms)?;
                opts.reactive_angles = coordinates(&self.angles, natoms)?;
                opts.reactive_dihedrals = coordinates(&self.dihedrals, natoms)?;
            }
        }
        let count =
            opts.reactive_bonds.len() + opts.reactive_angles.len() + opts.reactive_dihedrals.len();
        if count == 0
            && (opts.poormans_neb || opts.distance_mode == RdaDistanceMode::ReactiveInternal)
        {
            bail!("Define at least one reactive bond, angle, or dihedral for this calculation.");
        }
        let active = atom_indices(&self.active_atoms, natoms)?;
        opts.active_atoms = (!active.is_empty()).then_some(active);
        let align = if opts.align {
            atom_indices(&self.align_atoms, natoms)?
        } else {
            Vec::new()
        };
        opts.align_atoms = (!align.is_empty()).then_some(align);
        let weights = if opts.distance_mode == RdaDistanceMode::ReactiveInternal {
            numbers(&self.weights, "coordinate weight")?
        } else {
            Vec::new()
        };
        if !weights.is_empty() {
            if weights.len() != count || weights.iter().any(|&w| w <= 0.0) {
                bail!("Enter one positive weight per reactive coordinate ({count} total), in bond/angle/dihedral order.");
            }
            opts.reactive_coordinate_weights = Some(weights);
        } else {
            opts.reactive_coordinate_weights = None;
        }
        if opts.poormans_neb {
            opts.extra_bonds = coordinates(&self.extra_bonds, natoms)?;
            opts.extra_angles = coordinates(&self.extra_angles, natoms)?;
            opts.extra_dihedrals = coordinates(&self.extra_dihedrals, natoms)?;
        } else {
            opts.gamma_values = numbers(&self.gamma_values, "gamma")?;
            if opts.gamma_values.is_empty()
                || opts.gamma_values.iter().any(|&g| g <= 0.0 || g >= 1.0)
            {
                bail!("Gamma values must lie strictly between 0 and 1.");
            }
        }
        for (label, value) in [
            ("First energy tolerance", opts.first_energy_tol_ev),
            ("Later energy tolerance", opts.later_energy_tol_ev),
            ("Distance tolerance", opts.distance_tol),
            ("Cartesian step limit", opts.cartesian_step_max_component),
            ("Internal step limit", opts.internal_step_max_component),
            ("Beta step", opts.beta_step),
        ] {
            if !value.is_finite() || value <= 0.0 {
                bail!("{label} must be positive and finite.");
            }
        }
        if !opts.beta_start.is_finite() || !(0.0..=1.0).contains(&opts.beta_start) {
            bail!("Starting beta must lie between 0 and 1.");
        }
        if opts.max_conditional_steps == 0
            || opts.max_bracket_attempts == 0
            || opts.internal_best_fit_iters == 0
        {
            bail!("Iteration limits must be at least 1.");
        }
        if opts.poormans_neb && (opts.poormans_neb_images < 3 || opts.poormans_neb_maxiter == 0) {
            bail!(
                "Poor Man's NEB needs at least 3 images and at least 1 relaxation step per image."
            );
        }
        Ok(opts)
    }
}
