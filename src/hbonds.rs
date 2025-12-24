use bevy::prelude::*;

use crate::molecule::{Molecule, covalent_radius_angstrom};
use crate::settings::MolSettings;

#[derive(Default, Reflect, GizmoConfigGroup)]
pub struct HbondGizmos;

/// Configure dashed H-bond gizmo style (runs every frame; no conflict with drawing)
pub fn configure_hbond_gizmos(
    mut cfg_store: ResMut<GizmoConfigStore>,
    settings: Res<MolSettings>,
) {
    let (cfg, _) = cfg_store.config_mut::<HbondGizmos>();
    cfg.enabled = true;
    cfg.line.width = settings.hbond_thickness.max(1.0).min(50.0);
    cfg.line.perspective = true;
    cfg.line.style = GizmoLineStyle::Dashed {
        gap_scale:  settings.hbond_gap_scale,
        line_scale: settings.hbond_line_scale,
    };
}

/// Draw dashed hydrogen bonds as gizmo lines
pub fn draw_hydrogen_bonds_dashed(
    mut gizmos: Gizmos<HbondGizmos>,
    mol: Res<Molecule>,
    settings: Res<MolSettings>,
) {
    let hb_color = settings.hbond_color;

    let n_pos   = mol.pos.len();
    let n_atoms = mol.atoms.len();

    for &(i, j, hdist) in &mol.hydrogen_bonds {
        // Skip degenerate / “off” hbonds
        if hdist <= 0.0001 {
            continue;
        }

        // Defensive: skip bonds whose indices are no longer valid
        if i >= n_pos || j >= n_pos {
            continue;
        }
        if i >= n_atoms || j >= n_atoms {
            continue;
        }

        let p0 = mol.pos[i];
        let p1 = mol.pos[j];

        let d = p1 - p0;
        let len = d.length();
        if len <= 1e-4 {
            continue;
        }
        let dn = d / len;

        // Use a small fraction of each atom's visual radius so it doesn't depend on line width
        let ri = covalent_radius_angstrom(&mol.atoms[i]) * settings.atom_scale;
        let rj = covalent_radius_angstrom(&mol.atoms[j]) * settings.atom_scale;

        let surface_inset_frac = 0.12; // ~12% inside each sphere; tweak if you like
        let inset_i = (ri * surface_inset_frac).min(ri * 0.9);
        let inset_j = (rj * surface_inset_frac).min(rj * 0.9);

        let a = p0 + dn * (ri - inset_i); // start just inside sphere i
        let b = p1 - dn * (rj - inset_j); // end   just inside sphere j

        if (b - a).length() <= 1e-4 {
            continue;
        }

        gizmos.line(a, b, hb_color);
    }
}

