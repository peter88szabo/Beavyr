use super::{rotors::*, state::Partner, ts_mode};
use crate::ui_style::ResponseExt;
use bevy_egui::egui;

pub fn controls(ui: &mut egui::Ui, p: &mut Partner) {
    let Some(f) = p.fragment.as_ref() else {
        return;
    };
    if f.imaginary_cm1.len() == 1 {
        egui::CollapsingHeader::new(format!("TS mode ({:.1} cm⁻¹)", f.imaginary_cm1[0]))
            .default_open(true).show(ui, |ui| {
            // A rigid partner is never displaced, so the TS controls are inert.
            ui.add_enabled_ui(!p.rigid, |ui| {
                ui.checkbox(&mut p.ts_mode.enabled,"Sample TS direction");
            });
            ui.add_enabled_ui(p.ts_mode.enabled && !p.rigid, |ui| {
                ui.horizontal_wrapped(|ui| {
                    for side in [ts_mode::Side::Both,ts_mode::Side::Positive,ts_mode::Side::Negative] {
                        ui.selectable_value(&mut p.ts_mode.side,side,side.label());
                    }
                });
                ui.horizontal_wrapped(|ui| {
                    ui.label("Maximum atom displacement (Å)");
                    ui.add(egui::DragValue::new(&mut p.ts_mode.max_displacement_angstrom).speed(0.01).range(0.0001..=10.));
                });
                ui.small("Uniform displacement along the projected TS mode; no thermal/ZPE assignment.")
                    .on_hover_text("The limit applies to the TS contribution of the most displaced atom. Positive/negative refer to the eigenvector sign, not a known reactant/product. Vibrational and rotor displacements are additional.");
            });
        });
    }
    if f.atoms.len() < 4 {
        return;
    }
    egui::CollapsingHeader::new("Internal rotors (free-rotor sampling)").show(ui, |ui| {
        ui.small("Partner atom numbers start at 1. Angles are uniform over 0–360°.");
        ui.add_enabled_ui(!p.rigid, |ui| {
            let mut remove=None;
            for (index,r) in p.rotors.iter_mut().enumerate() {
                ui.push_id(index, |ui| {
                    ui.group(|ui| {
                        ui.horizontal(|ui| {
                            ui.checkbox(&mut r.enabled,format!("Rotor {}",index+1));
                            if ui.small_button("Remove").clicked_once() { remove=Some(index); }
                        });
                        ui.add_enabled_ui(r.enabled, |ui| {
                            ui.horizontal_wrapped(|ui| {
                                ui.selectable_value(&mut r.dihedral,true,"Dihedral");
                                ui.selectable_value(&mut r.dihedral,false,"Bond");
                                for (i,label) in ["A","B","C","D"].iter().enumerate() {
                                    if !r.dihedral && (i==0 || i==3) { continue; }
                                    ui.label(*label);
                                    ui.add(egui::DragValue::new(&mut r.atoms[i]).range(1..=f.atoms.len()));
                                }
                            });
                            egui::ComboBox::from_id_salt("fragment").selected_text(r.side.label()).show_ui(ui,|ui| {
                                for side in [FragmentSide::Smaller,FragmentSide::First,FragmentSide::Second,FragmentSide::Explicit] {
                                    ui.selectable_value(&mut r.side,side,side.label());
                                }
                            });
                            if r.side==FragmentSide::Explicit {
                                ui.add(egui::TextEdit::singleline(&mut r.atom_list).hint_text("Rotating atoms: 4,5,7-10").desired_width(220.));
                            }
                            match prepare(f,r) {
                                Ok(prepared) => { ui.small(format!("Rotating atoms: {}",prepared.moving.iter().map(|i|(i+1).to_string()).collect::<Vec<_>>().join(", "))); }
                                Err(e) => { ui.colored_label(egui::Color32::YELLOW,e); }
                            }
                            ui.horizontal_wrapped(|ui| {
                                ui.label("Associated vibration");
                                let text=if r.mode==0 { "Unassigned".into() } else { format!("Mode {}",r.mode) };
                                egui::ComboBox::from_id_salt("mode").selected_text(text).show_ui(ui,|ui| {
                                    ui.selectable_value(&mut r.mode,0,"Unassigned");
                                    for (i,w) in f.frequencies.iter().enumerate() {
                                        ui.selectable_value(&mut r.mode,i+1,format!("Mode {} · {:.1} cm⁻¹",i+1,w/super::input::CM_TO_AU));
                                    }
                                });
                            });
                            if r.mode>0 {
                                ui.horizontal_wrapped(|ui| {
                                    ui.selectable_value(&mut r.treatment,ModeTreatment::ZpeOnly,"ZPE only (Pebble)");
                                    ui.selectable_value(&mut r.treatment,ModeTreatment::Replace,"Replace harmonic mode");
                                });
                            } else {
                                ui.small("Assign the torsional vibration to control its harmonic contribution.");
                            }
                            egui::CollapsingHeader::new("Bond detection").show(ui,|ui| {
                                ui.horizontal_wrapped(|ui| {
                                    ui.label("H–X (Å)"); ui.add(egui::DragValue::new(&mut r.bond_hx_angstrom).speed(0.01).range(0.1..=5.));
                                    ui.label("X–X (Å)"); ui.add(egui::DragValue::new(&mut r.bond_xx_angstrom).speed(0.01).range(0.1..=5.));
                                });
                            });
                        });
                    });
                });
            }
            if let Some(i)=remove { p.rotors.remove(i); }
            if ui.button("Add rotor").clicked_once() { p.rotors.push(Rotor::default()); }
        });
    });
}
