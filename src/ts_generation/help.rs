//! In-app instructions for defining the reaction and interpreting TS controls.

use bevy_egui::egui;

fn bullet(ui: &mut egui::Ui, text: &str) {
    ui.horizontal_top(|ui| {
        ui.label("•");
        ui.add(egui::Label::new(text).wrap());
    });
}

fn heading(ui: &mut egui::Ui, text: &str) {
    ui.add_space(10.0);
    ui.separator();
    ui.label(egui::RichText::new(text).strong().size(17.0));
}

fn definition(ui: &mut egui::Ui, name: &str, example: &str, explanation: &str) {
    ui.group(|ui| {
        ui.strong(name);
        ui.monospace(example);
        ui.add(egui::Label::new(explanation).wrap());
    });
}

fn rda_options(ui: &mut egui::Ui) -> egui::CollapsingResponse<()> {
    egui::CollapsingHeader::new("RDA options — what each control means")
        .id_salt("ts_help_rda_options")
        .default_open(false)
        .show(ui, |ui| {
            ui.label("Start with the defaults. These settings control the short relaxations and the search between structures that relax in different directions.");
            for (name, explanation) in [
                ("Conditional relaxation: internal / Cartesian", "Chooses the coordinates used to relax each trial structure. Internal uses bond lengths, angles and dihedrals; Cartesian uses x, y and z. This is separate from the distance metric in Reaction definition."),
                ("Maximum conditional steps", "Maximum relaxation steps for each RDA trial structure. This is a per-trial limit, not a limit on the whole search."),
                ("First energy tolerance (eV)", "Stops the first midpoint relaxation when the absolute energy change between steps is below this value. A smaller tolerance asks for a smaller energy change."),
                ("Later energy tolerance (eV)", "Energy-change stopping tolerance for later RDA trial relaxations. In Poor Man's NEB, it is used for each constrained image relaxation."),
                ("Starting beta", "After relaxing the midpoint, RDA moves towards the opposite endpoint to look for a trial that relaxes in the other direction. Beta is the fraction along that move: 0 is the relaxed midpoint, 1 is the chosen endpoint."),
                ("Beta increment", "Amount added to beta for the next attempt when the search needs another trial. Trial beta is limited to the range 0–1."),
                ("Maximum bracket attempts", "Limits how many beta trials are attempted to find structures that relax in opposite directions. If no bracket is found, the result message reports this and the last candidate is returned."),
                ("Gamma values", "Fractions between the two bracket structures, each strictly between 0 and 1, separated by spaces or commas. RDA selects the trial whose distances to A and B are most nearly equal, then relaxes it. More values sample this interval more finely."),
                ("Distance tolerance (bohr / rad)", "Threshold for deciding whether relaxation changed the distances to A and B enough to assign a direction. It applies to the selected distance metric: Cartesian/bond distances use bohr; angles and dihedrals use radians. It is not an energy or force threshold."),
                ("Allow Cartesian relaxation if internal coordinates fail", "Allows an RDA trial to retry in Cartesian coordinates when its internal-coordinate relaxation fails. The energy backend and electronic method stay the same."),
                ("Cartesian step limit (bohr)", "Caps the largest x/y/z displacement component in an RDA Cartesian relaxation step. Smaller values make steps more conservative."),
                ("Internal step limit (bohr / rad)", "Caps the largest bond-length or angular change in an internal-coordinate step. Used by internal RDA relaxation and Poor Man's NEB image relaxation."),
                ("Internal-coordinate fit iterations", "Maximum fitting iterations for converting an internal-coordinate step back into Cartesian coordinates. This is separate from the number of optimization steps."),
            ] {
                ui.add_space(6.0);
                ui.strong(name);
                bullet(ui, explanation);
            }
        })
}

fn contents(ui: &mut egui::Ui) {
    ui.heading("How to use TS Generation");
    ui.label("Define reactant A, product B, and the motion that connects them. RDA generates a TS guess; Poor Man's NEB generates a constrained-relaxed path and selects its highest-energy interior image as the guess.");

    heading(ui, "1. Prepare and store A and B");
    bullet(ui, "Load an XYZ file for each endpoint, or draw/pre-optimize a structure in the viewer and press Use current as A. Prepare the product and press Use current as B.");
    bullet(ui, "A and B must have the same elements in the same atom order. Keep the identity of each atom, including any transferred H, unchanged between endpoints. Both structures must contain at least two atoms and must differ in geometry.");
    bullet(ui, "The green loaded marker confirms a stored endpoint. Show A and Show B display those copies. Editing the viewer does not update a stored endpoint: capture it again when you want to keep an edit. Pause trajectory playback before capturing a frame.");

    heading(ui, "2. Reaction definition — choose what changes");
    ui.label("Use the 1-based row numbers in Structure (XYZ): the first atom is 1. The numbers must refer to the same atoms in A and B. The examples below are placeholders; replace them with your own atom numbers. Spaces or commas separate atoms.");
    definition(ui, "Bond change (2 atoms)", "Atom numbers: 4 7", "Use for forming, breaking or stretching the 4–7 bond. Draw/load A and B with the intended initial and final bond lengths. This selects one bond-length coordinate.");
    definition(ui, "Atom transfer (3 atoms)", "Atom numbers: 4 7 12", "Enter donor, transferred atom, acceptor — in that order. This selects two bonds: 4–7 and 7–12. For H transfer, atom 7 is the same H in both structures: bonded towards donor 4 in A and towards acceptor 12 in B.");
    definition(ui, "Angle change (3 atoms)", "Atom numbers: 1 2 3", "Atom 2 is the vertex of angle 1–2–3. Use endpoints with the intended starting and finishing angles.");
    definition(ui, "Dihedral change (4 atoms)", "Atom numbers: 1 2 3 4", "Atoms 2–3 define the central bond of the 1–2–3–4 torsion. Use this for a change in torsion about that bond.");
    definition(ui, "Custom reactive coordinates", "Bonds: 4 7; 7 12\nAngles: 4 7 12\nDihedrals: 1 4 7 12", "Combine the coordinates that describe the reaction. Separate coordinates within a field with semicolons (or newlines). Leave unused fields blank. Each bond needs 2 atoms, each angle 3, and each dihedral 4; do not repeat an atom within one coordinate.");
    definition(ui, "Active atoms / Cartesian motion", "Active atoms: 4 7 12", "Use RDA when the motion is easier to describe with a group of atoms than with specific bonds or angles. The distance comparison uses the Cartesian displacement of these atoms. Blank uses all atoms. Poor Man's NEB requires at least one bond, angle or dihedral, so use one of the other definitions for it.");

    heading(ui, "3. Distance, active atoms and alignment");
    bullet(ui, "Distance in reactive coordinates compares A/B distances using the bonds, angles and dihedrals you entered. Cartesian distance on active atoms compares their x/y/z displacements. For RDA, this comparison determines whether a short relaxation moves towards A or B.");
    bullet(ui, "Active atoms selects the atoms measured by Cartesian distance, and supplies the default alignment subset. It does not freeze the other atoms. With reactive-coordinate distance, every listed reactive coordinate is measured, even if its atoms are not in Active atoms.");
    bullet(ui, "Coordinate weights: leave blank for equal weights, or enter one positive number per reactive coordinate, ordered as bonds, then angles, then dihedrals. For the transfer example, 1 1 weights its two bonds equally. Larger weights give a coordinate more influence in RDA's distance comparison. They do not change NEB's interpolated targets or constraint strength.");
    bullet(ui, "Align product to reactant removes overall translation/rotation before generating trials. Alignment atoms can select a relatively unchanged molecular framework; three or more non-collinear atoms help define its orientation. Blank uses Active atoms, or all atoms if that is also blank. Turn alignment off when the absolute relative placement of the endpoints is intentional.");

    heading(ui, "4. Choose RDA or Poor Man's NEB and generate");
    bullet(ui, "Choose the energy backend (xTB or Behemoth), executable and electronic method. Set the charge and spin multiplicity for the reaction before starting.");
    bullet(ui, "RDA uses short relaxations and the selected distance comparison to search for a TS guess between A and B. The detailed RDA controls are explained in the collapsed section below.");
    bullet(ui, "Poor Man's NEB interpolates the selected reactive coordinates between A and B. Each interior image relaxes while those coordinates are constrained; the endpoints are retained. Images includes A and B (minimum 3); more images sample the path more densely and cost more calculations.");
    bullet(ui, "Maximum relaxation steps per image limits each NEB image's optimization. Under Relaxation controls, Extra bonds / angles / dihedrals add connectivity to the internal-coordinate set; they are not additional reactive constraints. Use the same atom-list syntax as Custom reactive coordinates.");
    bullet(ui, "Press Generate RDA TS guess or Generate NEB path and TS guess. Progress shows energy/gradient calls and elapsed time. Cancel stops the calculation.");

    heading(ui, "5. Inspect and save the result");
    bullet(ui, "Show TS guess displays the candidate structure. Show path in Trajectory lets you inspect the path. Plot NEB energy shows all image energies relative to A in kJ/mol, with the selected TS guess marked; hover over a point for its exact energy.");
    bullet(ui, "Save TS guess XYZ writes one structure. Save NEB trajectory XYZ writes the complete path, including both endpoints. Save NEB energy profile exports the numerical energies; RDA has a separate search-structure export.");
    bullet(ui, "Read the result message and image relaxation messages. The output is a guess: a saddle-point optimization and a frequency calculation are needed to confirm the intended transition state.");
    ui.add_space(10.0);
    rda_options(ui);
}

/// The help window is independent of the TS panel and never changes a run.
pub fn window(ctx: &egui::Context, open: &mut bool) -> Option<egui::Rect> {
    if !*open {
        return None;
    }
    egui::Window::new("TS Generation — How To Use")
        .id(egui::Id::new("ts_generation_how_to_use"))
        .open(open)
        .resizable(true)
        .default_size([680.0, 650.0])
        .min_size([420.0, 300.0])
        .show(ctx, |ui| {
            egui::ScrollArea::vertical()
                .id_salt("ts_help_scroll")
                .show(ui, contents);
        })
        .map(|window| window.response.rect)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn detailed_rda_explanations_start_collapsed() {
        let ctx = egui::Context::default();
        let mut output = ctx.run_ui(egui::RawInput::default(), |ui| {
            let response = rda_options(ui);
            assert!(response.body_response.is_none());
        });
        output.textures_delta.clear();
    }

    #[test]
    fn help_window_opens_without_endpoints_or_a_calculation() {
        let ctx = egui::Context::default();
        let mut open = true;
        let mut rect = None;
        for _ in 0..2 {
            let mut output = ctx.run_ui(egui::RawInput::default(), |ui| {
                rect = window(ui.ctx(), &mut open);
            });
            assert!(!output.shapes.is_empty());
            output.textures_delta.clear();
        }
        assert!(open && rect.unwrap().is_finite());
    }
}
