//! Raster annotations on captured frames, using egui's bundled font.
use super::{extras::Extras, Request};
use ab_glyph::{Font, ScaleFont};
use bevy::prelude::*;
pub struct Annotation {
    pub text: String,
    pub position: Vec2,
    pub color: [u8; 3],
}
pub struct Overlay {
    pub text: Vec<Annotation>,
    pub lines: Vec<(Vec2, Vec2, [u8; 3], f32)>,
}
pub fn prepare(request: &Request, index: usize, size: UVec2) -> Overlay {
    let frame = &request.frames[if request.rotation { 0 } else { index }];
    let Extras {
        measurements,
        labels,
        ..
    } = &request.extras;
    let mut orbit = request.orbit.clone();
    if request.rotation {
        orbit.orientation =
            Quat::from_rotation_y(std::f32::consts::TAU * index as f32 / request.count as f32)
                * orbit.orientation;
    }
    let mut projection = request.projection.clone();
    projection.update(size.x as f32, size.y as f32);
    let clip = projection.get_clip_from_view()
        * crate::isolated_scene::camera_transform(&orbit)
            .to_matrix()
            .inverse();
    let project = |p: Vec3| -> Option<Vec2> {
        let p = clip * p.extend(1.);
        if p.w <= 0. {
            return None;
        }
        let p = p.truncate() / p.w;
        if p.z < 0. || p.z > 1. {
            return None;
        }
        Some(Vec2::new(
            (p.x + 1.) * 0.5 * size.x as f32,
            (1. - p.y) * 0.5 * size.y as f32,
        ))
    };
    let atom = |i: usize| {
        if request.hidden.contains(&i) {
            None
        } else {
            frame.pos.get(i).copied()
        }
    };
    let mut overlay = Overlay {
        text: vec![],
        lines: vec![],
    };
    if request.extras.options.caption {
        let source = request.source_indices[if request.rotation { 0 } else { index }] + 1;
        let energy = frame
            .metadata
            .energy_hartree
            .map(|e| {
                format!(
                    " · {}E = {e:.8} Eh",
                    if frame.metadata.energy_is_relative {
                        "Δ"
                    } else {
                        ""
                    }
                )
            })
            .unwrap_or_default();
        overlay.text.push(Annotation {
            text: format!("Frame {source}{energy}"),
            position: Vec2::new(12., 12.),
            color: [255, 255, 255],
        });
    }
    if let Some(m) = measurements {
        let c = m.color.to_srgba();
        let color = [
            (c.red * 255.) as u8,
            (c.green * 255.) as u8,
            (c.blue * 255.) as u8,
        ];
        for pair in &m.pairs {
            if !pair.visible {
                continue;
            }
            if let (Some(a), Some(b)) = (atom(pair.a), atom(pair.b)) {
                if let (Some(pa), Some(pb)) = (project(a), project(b)) {
                    overlay.lines.push((pa, pb, color, m.line_width));
                    if m.show_labels {
                        overlay.text.push(Annotation {
                            text: format!("{:.3} Å", a.distance(b)),
                            position: (pa + pb) * 0.5,
                            color,
                        });
                    }
                }
            }
        }
        for angle in &m.angles {
            if let (Some(a), Some(b), Some(c)) = (atom(angle.a), atom(angle.b), atom(angle.c)) {
                if let (Some(x), Some(y), Some(p)) =
                    ((a - b).try_normalize(), (c - b).try_normalize(), project(b))
                {
                    overlay.text.push(Annotation {
                        text: format!("{:.2}°", x.dot(y).clamp(-1., 1.).acos().to_degrees()),
                        position: p,
                        color,
                    });
                }
            }
        }
        for d in &m.dihedrals {
            let indices = [d.a + 1, d.b + 1, d.c + 1, d.d + 1];
            if indices.iter().any(|i| request.hidden.contains(&(i - 1))) {
                continue;
            }
            if let Some(v) =
                crate::polarizability::data::Axis::Dihedral.value(frame, 0, &indices, 1., 0.)
            {
                if let Some(p) = atom(d.b)
                    .zip(atom(d.c))
                    .and_then(|(b, c)| project((b + c) * 0.5))
                {
                    overlay.text.push(Annotation {
                        text: format!("{v:.2}°"),
                        position: p,
                        color,
                    });
                }
            }
        }
    }
    if let Some(labels) = labels {
        if labels.atoms == frame.atoms && labels.show_labels {
            for label in &labels.labels {
                if let Some(p) = atom(label.atom).and_then(project) {
                    overlay.text.push(Annotation {
                        text: label.text.clone(),
                        position: p + Vec2::new(8., -16.),
                        color: [255, 255, 255],
                    });
                }
            }
        }
    }
    if request.extras.polar_labels {
        if let Some(values) = request.extras.values(frame) {
            for (i, a) in values.alpha.iter().enumerate() {
                if let Some(p) = atom(i).and_then(project) {
                    let factor = if request.extras.polar_angstrom {
                        crate::polarizability::data::BOHR3_TO_ANGSTROM3
                    } else {
                        1.
                    };
                    let unit = if request.extras.polar_angstrom {
                        "Å³"
                    } else {
                        "bohr³"
                    };
                    overlay.text.push(Annotation {
                        text: format!("{:.3} {unit}", a * factor),
                        position: p + Vec2::new(8., -16.),
                        color: [255, 255, 255],
                    });
                }
            }
        }
    }
    if let Some(r) = &request.extras.reactivity {
        if r.labels {
            if let Some(f) = r.results.frames.iter().find(|f| {
                crate::trajectory_wavefunctions::job::geometry_matches(&f.geometry, frame)
            }) {
                if let Some(values) = r.descriptor.atoms(f) {
                    for (i, value) in values.iter().enumerate() {
                        if let Some(p) = atom(i).and_then(project) {
                            overlay.text.push(Annotation {
                                text: format!("{value:+.4}"),
                                position: p + Vec2::new(8., 8.),
                                color: [255, 255, 255],
                            });
                        }
                    }
                }
            }
        }
    }
    overlay
}
fn blend(pixels: &mut [u8], size: UVec2, x: i32, y: i32, color: [u8; 3], alpha: f32) {
    if x < 0 || y < 0 || x >= size.x as i32 || y >= size.y as i32 {
        return;
    }
    let i = (y as usize * size.x as usize + x as usize) * 4;
    for c in 0..3 {
        pixels[i + c] =
            (pixels[i + c] as f32 * (1. - alpha) + color[c] as f32 * alpha).round() as u8;
    }
}
pub fn paint(overlay: &Overlay, pixels: &mut [u8], size: UVec2) -> Result<(), String> {
    if pixels.len() != size.x as usize * size.y as usize * 4 {
        return Err("Invalid animation pixel buffer.".into());
    }
    for &(a, b, color, width) in &overlay.lines {
        let steps = a.distance(b).ceil().min(20000.) as usize;
        let radius = (width / 2.).round().clamp(0., 8.) as i32;
        for n in 0..=steps {
            if n % 12 >= 7 {
                continue;
            }
            let p = a.lerp(b, n as f32 / steps.max(1) as f32);
            for dx in -radius..=radius {
                for dy in -radius..=radius {
                    blend(pixels, size, p.x as i32 + dx, p.y as i32 + dy, color, 0.9);
                }
            }
        }
    }
    if overlay.text.is_empty() {
        return Ok(());
    }
    let fonts = bevy_egui::egui::FontDefinitions::default();
    let data = fonts
        .font_data
        .get("Ubuntu-Light")
        .or_else(|| fonts.font_data.values().next())
        .ok_or("Bundled font missing")?;
    let font = ab_glyph::FontRef::try_from_slice(&data.font).map_err(|e| e.to_string())?;
    let scaled = font.as_scaled(16.);
    for annotation in &overlay.text {
        let width = annotation
            .text
            .chars()
            .map(|c| scaled.h_advance(font.glyph_id(c)))
            .sum::<f32>();
        let pos = Vec2::new(
            annotation
                .position
                .x
                .clamp(2., (size.x as f32 - width - 2.).max(2.)),
            annotation
                .position
                .y
                .clamp(2., (size.y as f32 - 20.).max(2.)),
        );
        for (offset, color) in [(Vec2::ONE, [0, 0, 0]), (Vec2::ZERO, annotation.color)] {
            let mut x = pos.x + offset.x;
            for ch in annotation.text.chars() {
                let id = font.glyph_id(ch);
                let glyph = id.with_scale_and_position(
                    16.,
                    ab_glyph::point(x, pos.y + offset.y + scaled.ascent()),
                );
                if let Some(g) = font.outline_glyph(glyph) {
                    let bounds = g.px_bounds();
                    g.draw(|gx, gy, a| {
                        blend(
                            pixels,
                            size,
                            bounds.min.x as i32 + gx as i32,
                            bounds.min.y as i32 + gy as i32,
                            color,
                            a,
                        )
                    });
                }
                x += scaled.h_advance(id);
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn exported_distance_labels_follow_geometry_and_caption_keeps_source_frame() {
        let mut mol = crate::molecule::Molecule::from_xyz("2\nframe\nH 0 0 0\nH 1 0 0\n");
        let mut m = crate::measurements::Measurements::default();
        m.new_pair(0, 1, 999.);
        let mut extras = Extras::default();
        extras.options.caption = true;
        extras.measurements = Some(m);
        let a = crate::trajectory::TrajectoryFrame {
            atoms: mol.atoms.clone(),
            pos: mol.pos.clone(),
            metadata: Default::default(),
        };
        mol.pos[1].x = 2.;
        let b = crate::trajectory::TrajectoryFrame {
            atoms: mol.atoms.clone(),
            pos: mol.pos.clone(),
            metadata: Default::default(),
        };
        let request = Request {
            extras,
            source_indices: vec![3, 7],
            frames: vec![a, b],
            count: 2,
            rotation: false,
            template: mol,
            settings: Default::default(),
            hidden: Default::default(),
            orbit: Default::default(),
            projection: Default::default(),
            path: "unused.gif".into(),
            fixed_bonds: false,
        };
        let size = UVec2::new(640, 480);
        let first = prepare(&request, 0, size);
        let last = prepare(&request, 1, size);
        assert!(first.text.iter().any(|t| t.text == "1.000 Å"));
        assert!(last.text.iter().any(|t| t.text == "2.000 Å"));
        assert!(last.text.iter().any(|t| t.text == "Frame 8"));
        let mut rgba = vec![0; 640 * 480 * 4];
        paint(&last, &mut rgba, size).unwrap();
        assert!(rgba
            .chunks_exact(4)
            .any(|p| p[0] > 0 || p[1] > 0 || p[2] > 0));
    }
}
