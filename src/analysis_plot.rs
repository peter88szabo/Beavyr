//! Shared publication plot primitives: the screen, SVG and PNG use the same layout.
use crate::ui_style::ResponseExt;
use bevy_egui::egui::{self, Color32, Pos2};
#[derive(Clone)]
pub struct Style {
    pub font: f32,
    pub line: u8,
    pub points: u8,
    pub width: f32,
    pub radius: f32,
    pub color: Color32,
    pub export_width: u32,
    pub export_height: u32,
}
impl Default for Style {
    fn default() -> Self {
        Self {
            font: 13.,
            line: 0,
            points: 0,
            width: 1.5,
            radius: 2.5,
            color: Color32::from_rgb(55, 145, 205),
            export_width: 1600,
            export_height: 1000,
        }
    }
}
pub struct Plot<'a> {
    pub points: &'a [Option<(f64, f64)>],
    pub x: &'a str,
    pub y: &'a str,
    pub note: &'a str,
    pub reverse_x: bool,
    pub point_ids: Option<&'a [usize]>,
}
#[derive(Clone)]
enum Primitive {
    Line(Pos2, Pos2, f32, Color32),
    Point(Pos2, f32, bool, Color32),
    Text(Pos2, f32, f32, String, Color32),
}
fn limits(points: &[Option<(f64, f64)>], x: bool) -> (f64, f64) {
    let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
    for (a, b) in points
        .iter()
        .flatten()
        .filter(|(a, b)| a.is_finite() && b.is_finite())
    {
        let v = if x { *a } else { *b };
        lo = lo.min(v);
        hi = hi.max(v);
    }
    if !lo.is_finite() {
        return (0., 1.);
    }
    let pad = ((hi - lo) * 0.05).max(1e-6_f64.max(lo.abs() * 1e-8));
    (lo - pad, hi + pad)
}
fn tick(x: f64) -> String {
    if x != 0. && (x.abs() < 0.001 || x.abs() >= 1e5) {
        format!("{x:.2e}")
    } else {
        format!("{x:.3}")
    }
}
fn scene(
    plot: &Plot<'_>,
    style: &Style,
    size: egui::Vec2,
    ink: Color32,
) -> (Vec<Primitive>, Vec<Option<Pos2>>) {
    let f = style.font;
    let r = egui::Rect::from_min_max(
        egui::pos2(7. * f, 3.5 * f),
        egui::pos2(size.x - 3. * f, size.y - 5. * f),
    );
    let (xl, xh) = limits(plot.points, true);
    let (yl, yh) = limits(plot.points, false);
    let map = |x: f64, y: f64| {
        let t = ((x - xl) / (xh - xl)) as f32;
        egui::pos2(
            r.left() + r.width() * if plot.reverse_x { 1. - t } else { t },
            r.bottom() - r.height() * ((y - yl) / (yh - yl)) as f32,
        )
    };
    let mut out = vec![
        Primitive::Text(egui::pos2(12., 4.), 0., f, plot.y.into(), ink),
        Primitive::Text(
            egui::pos2(size.x / 2., size.y - 3. * f),
            0.5,
            f,
            plot.x.into(),
            ink,
        ),
        Primitive::Text(
            egui::pos2(size.x / 2., size.y - 1.5 * f),
            0.5,
            f * 0.85,
            plot.note.into(),
            ink,
        ),
    ];
    for i in 0..=4 {
        let t = i as f64 / 4.;
        let y = yl + (yh - yl) * t;
        let x = xl + (xh - xl) * t;
        let p = map(x, y);
        out.push(Primitive::Line(
            egui::pos2(r.left(), p.y),
            egui::pos2(r.right(), p.y),
            0.6,
            Color32::from_gray(170),
        ));
        out.push(Primitive::Text(
            egui::pos2(r.left() - 7., p.y - f / 2.),
            1.,
            f * 0.85,
            tick(y),
            ink,
        ));
        out.push(Primitive::Text(
            egui::pos2(p.x, r.bottom() + 6.),
            0.5,
            f * 0.85,
            tick(x),
            ink,
        ));
    }
    out.push(Primitive::Line(r.left_top(), r.left_bottom(), 1., ink));
    out.push(Primitive::Line(r.left_bottom(), r.right_bottom(), 1., ink));
    let positions: Vec<_> = plot
        .points
        .iter()
        .map(|p| {
            p.filter(|(x, y)| x.is_finite() && y.is_finite())
                .map(|(x, y)| map(x, y))
        })
        .collect();
    let mut previous: Option<Pos2> = None;
    for p in &positions {
        if let Some(p) = p {
            if let Some(old) = previous {
                if style.line == 0 {
                    out.push(Primitive::Line(old, *p, style.width, style.color));
                }
                if style.line == 1 {
                    let length = old.distance(*p);
                    let mut a = 0.;
                    while a < length {
                        out.push(Primitive::Line(
                            old.lerp(*p, a / length),
                            old.lerp(*p, (a + 6.).min(length) / length),
                            style.width,
                            style.color,
                        ));
                        a += 10.;
                    }
                }
            }
            if style.points < 2 {
                out.push(Primitive::Point(
                    *p,
                    style.radius,
                    style.points == 1,
                    style.color,
                ));
            }
            previous = Some(*p);
        } else {
            previous = None;
        }
    }
    (out, positions)
}
pub fn draw(ui: &mut egui::Ui, id: &str, plot: Plot<'_>, selected: Option<usize>) -> Option<usize> {
    let key = ui.make_persistent_id(id);
    let mut style = ui
        .ctx()
        .data_mut(|d| d.get_temp::<Style>(key))
        .unwrap_or_default();
    let mut export = None;
    egui::CollapsingHeader::new("Plot appearance / export")
        .id_salt(id)
        .default_open(false)
        .show(ui, |ui| {
            ui.horizontal_wrapped(|ui| {
                ui.label("Font");
                ui.add(egui::DragValue::new(&mut style.font).range(10.0..=24.));
                ui.label("Line");
                egui::ComboBox::from_id_salt((id, "line"))
                    .selected_text(["Solid", "Dashed", "None"][style.line as usize])
                    .show_ui(ui, |ui| {
                        for (i, s) in ["Solid", "Dashed", "None"].into_iter().enumerate() {
                            ui.selectable_value(&mut style.line, i as u8, s);
                        }
                    });
                ui.label("Points");
                egui::ComboBox::from_id_salt((id, "points"))
                    .selected_text(["Circle", "Square", "None"][style.points as usize])
                    .show_ui(ui, |ui| {
                        for (i, s) in ["Circle", "Square", "None"].into_iter().enumerate() {
                            ui.selectable_value(&mut style.points, i as u8, s);
                        }
                    });
                ui.color_edit_button_srgba(&mut style.color);
            });
            ui.horizontal_wrapped(|ui| {
                ui.label("Line width");
                ui.add(
                    egui::DragValue::new(&mut style.width)
                        .range(0.5..=5.)
                        .speed(0.1),
                );
                ui.label("Point size");
                ui.add(
                    egui::DragValue::new(&mut style.radius)
                        .range(1.0..=8.)
                        .speed(0.2),
                );
                ui.label("Image size");
                ui.add(egui::DragValue::new(&mut style.export_width).range(640..=4000));
                ui.label("×");
                ui.add(egui::DragValue::new(&mut style.export_height).range(480..=3000));
                for ext in ["svg", "png"] {
                    if ui
                        .button(format!("Save {}…", ext.to_uppercase()))
                        .clicked_once()
                    {
                        export = Some(ext);
                    }
                }
            });
        });
    ui.ctx().data_mut(|d| d.insert_temp(key, style.clone()));
    if let Some(ext) = export {
        if let Some(path) = rfd::FileDialog::new()
            .set_file_name(format!("analysis.{ext}"))
            .add_filter(ext, &[ext])
            .save_file()
        {
            if let Err(e) = save(&path.with_extension(ext), &plot, &style, ext) {
                crate::runtime_warning::push(ui.ctx(), "Plot export", &e);
            }
        }
    }
    let (response, painter) = ui.allocate_painter(
        egui::vec2(
            ui.available_width().max(330.),
            ui.available_height().max(280.),
        ),
        egui::Sense::click(),
    );
    let (primitives, positions) = scene(
        &plot,
        &style,
        response.rect.size(),
        ui.visuals().text_color(),
    );
    let offset = response.rect.min.to_vec2();
    for p in primitives {
        match p {
            Primitive::Line(a, b, w, c) => {
                painter.line_segment([a + offset, b + offset], egui::Stroke::new(w, c));
            }
            Primitive::Point(p, r, square, c) => {
                if square {
                    painter.rect_filled(
                        egui::Rect::from_center_size(p + offset, egui::vec2(2. * r, 2. * r)),
                        0.,
                        c,
                    );
                } else {
                    painter.circle_filled(p + offset, r, c);
                }
            }
            Primitive::Text(p, a, f, t, c) => {
                painter.text(
                    p + offset,
                    if a == 0. {
                        egui::Align2::LEFT_TOP
                    } else if a == 1. {
                        egui::Align2::RIGHT_TOP
                    } else {
                        egui::Align2::CENTER_TOP
                    },
                    t,
                    egui::FontId::proportional(f),
                    c,
                );
            }
        }
    }
    if let Some(p) = selected.and_then(|i| positions.get(i)).and_then(|p| *p) {
        painter.circle_stroke(
            p + offset,
            style.radius + 3.,
            egui::Stroke::new(2., Color32::YELLOW),
        );
    }
    let nearest = response.hover_pos().and_then(|mouse| {
        positions
            .iter()
            .enumerate()
            .filter_map(|(i, p)| p.map(|p| (i, mouse.distance_sq(p + offset))))
            .filter(|(_, d)| *d < 225.)
            .min_by(|a, b| a.1.total_cmp(&b.1))
            .map(|(i, _)| i)
    });
    if let Some(i) = nearest {
        let (x, y) = plot.points[i].unwrap();
        response.clone().on_hover_text(format!(
            "Point {}\n{}: {x:.8}\n{}: {y:.8}\nClick to display",
            plot.point_ids
                .and_then(|ids| ids.get(i))
                .copied()
                .unwrap_or(i)
                + 1,
            plot.x,
            plot.y
        ));
        if response.clicked_once() {
            return Some(i);
        }
    }
    None
}
fn escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}
fn rgb(c: Color32) -> String {
    format!("#{:02x}{:02x}{:02x}", c.r(), c.g(), c.b())
}
pub fn save(
    path: &std::path::Path,
    plot: &Plot<'_>,
    style: &Style,
    ext: &str,
) -> Result<(), String> {
    let w = style.export_width;
    let h = style.export_height;
    let scale = (w as f32 / 800.).min(h as f32 / 500.);
    let mut s = style.clone();
    s.font *= scale;
    s.width *= scale;
    s.radius *= scale;
    let (scene, _) = scene(plot, &s, egui::vec2(w as f32, h as f32), Color32::BLACK);
    if ext == "svg" {
        let mut svg=format!("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{w}\" height=\"{h}\" viewBox=\"0 0 {w} {h}\"><rect width=\"100%\" height=\"100%\" fill=\"white\"/>");
        for p in scene {
            match p{
            Primitive::Line(a,b,width,c)=>svg.push_str(&format!("<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"{}\" stroke-width=\"{width}\"/>",a.x,a.y,b.x,b.y,rgb(c))),
            Primitive::Point(p,r,square,c)=>{if square{svg.push_str(&format!("<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" fill=\"{}\"/>",p.x-r,p.y-r,r*2.,r*2.,rgb(c)));}else{svg.push_str(&format!("<circle cx=\"{}\" cy=\"{}\" r=\"{r}\" fill=\"{}\"/>",p.x,p.y,rgb(c)));}},
            Primitive::Text(p,a,f,t,c)=>svg.push_str(&format!("<text x=\"{}\" y=\"{}\" text-anchor=\"{}\" font-family=\"sans-serif\" font-size=\"{f}\" fill=\"{}\">{}</text>",p.x,p.y+f*0.85,if a==0.{"start"}else if a==1.{"end"}else{"middle"},rgb(c),escape(&t)))
        }
        }
        svg.push_str("</svg>");
        return std::fs::write(path, svg).map_err(|e| e.to_string());
    }
    raster(path, w, h, scene)
}
fn raster(path: &std::path::Path, w: u32, h: u32, scene: Vec<Primitive>) -> Result<(), String> {
    use ab_glyph::{Font, ScaleFont};
    use bevy::{
        asset::RenderAssetUsages,
        render::render_resource::{Extent3d, TextureDimension, TextureFormat},
    };
    let fonts = egui::FontDefinitions::default();
    let data = fonts
        .font_data
        .get("Ubuntu-Light")
        .or_else(|| fonts.font_data.values().next())
        .ok_or("Bundled font missing")?;
    let font = ab_glyph::FontRef::try_from_slice(&data.font).map_err(|e| e.to_string())?;
    let mut pixels = vec![255u8; w as usize * h as usize * 4];
    let mut paint = |x: i32, y: i32, c: Color32, a: f32| {
        if x >= 0 && y >= 0 && x < w as i32 && y < h as i32 {
            let i = (y as usize * w as usize + x as usize) * 4;
            for (k, v) in [c.r(), c.g(), c.b()].into_iter().enumerate() {
                pixels[i + k] = (pixels[i + k] as f32 * (1. - a) + v as f32 * a) as u8;
            }
        }
    };
    for p in scene {
        match p {
            Primitive::Text(p, align, size, text, c) => {
                let scaled = font.as_scaled(size);
                let width = text
                    .chars()
                    .map(|c| scaled.h_advance(font.glyph_id(c)))
                    .sum::<f32>();
                let mut x = p.x - align * width;
                for ch in text.chars() {
                    let id = font.glyph_id(ch);
                    if let Some(g) = font.outline_glyph(
                        id.with_scale_and_position(size, ab_glyph::point(x, p.y + scaled.ascent())),
                    ) {
                        let b = g.px_bounds();
                        g.draw(|a, bx, v| {
                            paint(b.min.x as i32 + a as i32, b.min.y as i32 + bx as i32, c, v)
                        });
                    }
                    x += scaled.h_advance(id);
                }
            }
            Primitive::Point(p, r, square, c) => {
                for y in (p.y - r - 1.) as i32..=(p.y + r + 1.) as i32 {
                    for x in (p.x - r - 1.) as i32..=(p.x + r + 1.) as i32 {
                        let d = if square {
                            (x as f32 - p.x).abs().max((y as f32 - p.y).abs())
                        } else {
                            egui::pos2(x as f32, y as f32).distance(p)
                        };
                        paint(x, y, c, (r + 0.5 - d).clamp(0., 1.));
                    }
                }
            }
            Primitive::Line(a, b, width, c) => {
                let d = b - a;
                let n = d.length_sq();
                for y in (a.y.min(b.y) - width) as i32..=(a.y.max(b.y) + width) as i32 {
                    for x in (a.x.min(b.x) - width) as i32..=(a.x.max(b.x) + width) as i32 {
                        let p = egui::pos2(x as f32, y as f32);
                        let t = ((p - a).dot(d) / n.max(1e-12)).clamp(0., 1.);
                        let dist = p.distance(a + d * t);
                        paint(x, y, c, (width / 2. + 0.5 - dist).clamp(0., 1.));
                    }
                }
            }
        }
    }
    bevy::prelude::Image::new(
        Extent3d {
            width: w,
            height: h,
            depth_or_array_layers: 1,
        },
        TextureDimension::D2,
        pixels,
        TextureFormat::Rgba8UnormSrgb,
        RenderAssetUsages::default(),
    )
    .try_into_dynamic()
    .map_err(|e| e.to_string())?
    .save(path)
    .map_err(|e| e.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn svg_png_preserve_units_reference_and_gaps() {
        let root = std::env::temp_dir().join(format!("beavyr-plot-{}", std::process::id()));
        std::fs::create_dir_all(&root).unwrap();
        let plot = Plot {
            points: &[Some((3., 2.)), None, Some((1., -1.))],
            x: "Distance (Å)",
            y: "ΔE (kcal/mol)",
            note: "Reference: point 3 & test",
            reverse_x: true,
            point_ids: None,
        };
        let style = Style {
            export_width: 800,
            export_height: 600,
            ..Default::default()
        };
        save(&root.join("plot.svg"), &plot, &style, "svg").unwrap();
        save(&root.join("plot.png"), &plot, &style, "png").unwrap();
        let svg = std::fs::read_to_string(root.join("plot.svg")).unwrap();
        assert!(svg.contains("kcal/mol"));
        assert!(svg.contains("Reference: point 3 &amp; test"));
        assert_eq!(svg.matches("<circle").count(), 2);
        let png = std::fs::read(root.join("plot.png")).unwrap();
        assert_eq!(&png[..8], b"\x89PNG\r\n\x1a\n");
        let (primitives, positions) = scene(&plot, &style, egui::vec2(800., 600.), Color32::BLACK);
        assert!(positions[0].unwrap().x < positions[2].unwrap().x);
        assert!(positions[1].is_none());
        assert_eq!(
            primitives
                .iter()
                .filter(|p| matches!(p,Primitive::Line(_,_,_,c) if *c==style.color))
                .count(),
            0
        );
        // Keep fixtures for independent image inspection in the test report.
    }
}
