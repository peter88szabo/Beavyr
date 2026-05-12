// src/trajectory.rs
use bevy::prelude::*;
use std::path::PathBuf;

use crate::events::MoleculeChanged;
use crate::molecule::{covalent_radius_angstrom, Molecule};
use crate::settings::MolSettings;
use crate::scene::LAYER_MAIN;
use bevy::camera::visibility::RenderLayers;

#[derive(Component)]
pub struct TrajGhostAtom;
#[derive(Component)]
pub struct TrajGhostBond;
#[derive(Component)]
pub struct TrajFullAtom;
#[derive(Component)]
pub struct TrajFullBond;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PlaybackMode {
    Loop,
    Once,
    PingPong,
}

#[derive(Clone)]
pub struct TrajectoryFrame {
    pub atoms: Vec<String>,
    pub pos: Vec<Vec3>,
}

#[derive(Resource, Clone)]
pub struct TrajectoryState {
    pub frames: Vec<TrajectoryFrame>,
    pub current_frame: usize,
    pub playing: bool,
    pub direction: i32,
    pub mode: PlaybackMode,
    pub overlay_enabled: bool,
    pub overlay_dirty: bool,
    pub overlay_full_stride: usize,
    pub overlay_ghost_stride: usize,
    pub overlay_full_last: bool,
    pub ghost_alpha: f32,
    pub fixed_bonds: bool,
    pub fps: f32,
    pub accum: f32,
    pub last_applied: Option<usize>,
    pub last_dir: Option<PathBuf>,
    pub current_file: Option<PathBuf>,
}

impl Default for TrajectoryState {
    fn default() -> Self {
        Self {
            frames: Vec::new(),
            current_frame: 0,
            playing: false,
            direction: 1,
            mode: PlaybackMode::Loop,
            overlay_enabled: false,
            overlay_dirty: false,
            overlay_full_stride: 1,
            overlay_ghost_stride: 1,
            overlay_full_last: false,
            ghost_alpha: 0.2,
            fixed_bonds: false,
            fps: 12.0,
            accum: 0.0,
            last_applied: None,
            last_dir: None,
            current_file: None,
        }
    }
}

pub fn parse_multi_xyz(text: &str) -> Result<Vec<TrajectoryFrame>, String> {
    let mut lines = text.lines().enumerate().peekable();
    let mut frames: Vec<TrajectoryFrame> = Vec::new();
    let mut frame_idx: usize = 0;

    let looks_like_coord = |line: &str| -> bool {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 4 {
            return false;
        }
        parts[1].parse::<f32>().is_ok()
            && parts[2].parse::<f32>().is_ok()
            && parts[3].parse::<f32>().is_ok()
    };

    let mut expected_count: Option<usize> = None;
    loop {
        // Skip blank lines before the count line.
        let mut count_line = None;
        while let Some((line_no, line)) = lines.next() {
            if !line.trim().is_empty() {
                count_line = Some((line_no + 1, line));
                break;
            }
        }
        let Some(count_line) = count_line else { break };

        let (count_line_no, count_line_text) = count_line;
        let count = count_line_text
            .trim()
            .parse::<usize>()
            .map_err(|_| {
                format!(
                    "Invalid atom count line in trajectory (frame {}, line {}): {}",
                    frame_idx + 1,
                    count_line_no,
                    count_line_text.trim()
                )
            })?;
        if let Some(expected) = expected_count {
            if count != expected {
                return Err(format!(
                    "Trajectory frames have inconsistent atom counts (frame {}, line {}, expected {}, got {}): {}",
                    frame_idx + 1,
                    count_line_no,
                    expected,
                    count,
                    count_line_text.trim()
                ));
            }
        } else {
            expected_count = Some(count);
        }

        // Optional comment line: consume one line unless it looks like a coordinate line.
        if let Some((_, next)) = lines.peek() {
            if next.trim().is_empty() || !looks_like_coord(next) {
                let _ = lines.next();
            }
        }

        let mut frame_atoms: Vec<String> = Vec::with_capacity(count);
        let mut frame_pos: Vec<Vec3> = Vec::with_capacity(count);
        while frame_atoms.len() < count {
            let Some((line_no, line)) = lines.next() else {
                return Err(format!(
                    "Unexpected end of frame data (frame {}, line {})",
                    frame_idx + 1,
                    count_line_no
                ));
            };
            if line.trim().is_empty() {
                continue;
            }
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() < 4 {
                return Err(format!(
                    "Invalid coordinate line in trajectory (frame {}, line {}): {}",
                    frame_idx + 1,
                    line_no + 1,
                    line.trim()
                ));
            }
            let sym = parts[0].to_string();
            let x = parts[1].parse::<f32>().map_err(|_| {
                format!(
                    "Invalid X coordinate in trajectory (frame {}, line {}): {}",
                    frame_idx + 1,
                    line_no + 1,
                    line.trim()
                )
            })?;
            let y = parts[2].parse::<f32>().map_err(|_| {
                format!(
                    "Invalid Y coordinate in trajectory (frame {}, line {}): {}",
                    frame_idx + 1,
                    line_no + 1,
                    line.trim()
                )
            })?;
            let z = parts[3].parse::<f32>().map_err(|_| {
                format!(
                    "Invalid Z coordinate in trajectory (frame {}, line {}): {}",
                    frame_idx + 1,
                    line_no + 1,
                    line.trim()
                )
            })?;
            frame_atoms.push(sym);
            frame_pos.push(Vec3::new(x, y, z));
        }

        frames.push(TrajectoryFrame {
            atoms: frame_atoms,
            pos: frame_pos,
        });
        frame_idx += 1;
    }

    if frames.is_empty() {
        return Err("No frames found in trajectory".to_string());
    }

    Ok(frames)
}

pub fn apply_current_frame(
    traj: &mut TrajectoryState,
    mol: &mut Molecule,
    settings: &mut MolSettings,
    ev_changed: &mut MessageWriter<MoleculeChanged>,
    recenter: bool,
    cam: Option<&mut crate::camera::OrbitCamera>,
) {
    if traj.frames.is_empty() {
        return;
    }
    let idx = traj.current_frame.min(traj.frames.len() - 1);
    if traj.last_applied == Some(idx) && !recenter {
        return;
    }

    let frame = &traj.frames[idx];
    let old_bond_pairs: Vec<(usize, usize)> = mol.bonds.iter().map(|&(i, j, _)| (i, j)).collect();
    let topology_changed = mol.atoms != frame.atoms || mol.pos.len() != frame.pos.len();
    mol.atoms = frame.atoms.clone();
    mol.set_pos(frame.pos.clone());

    let mut bond_topology_changed = false;
    if topology_changed || !traj.fixed_bonds {
        mol.recompute_bonds(settings.bond_thresh_scale, settings.hbond_cutoff);
        bond_topology_changed = old_bond_pairs.len() != mol.bonds.len()
            || old_bond_pairs
                .iter()
                .zip(mol.bonds.iter())
                .any(|(&(oi, oj), &(ni, nj, _))| oi != ni || oj != nj);
    }

    if topology_changed || bond_topology_changed {
        settings.geometry_dirty = true;
        settings.bond_topology_dirty = true;
    } else {
        settings.coords_dirty = true;
    }
    ev_changed.write(MoleculeChanged::set_pos());

    if recenter {
        if let Some(c) = compute_centroid(&mol.pos) {
            if let Some(cam) = cam {
                cam.target = c;
            }
        }
    }

    traj.last_applied = Some(idx);
}

pub fn advance_trajectory(
    time: Res<Time>,
    mut traj: ResMut<TrajectoryState>,
    mut mol: ResMut<Molecule>,
    mut settings: ResMut<MolSettings>,
    mut ev_changed: MessageWriter<MoleculeChanged>,
) {
    if !traj.playing || traj.frames.is_empty() || traj.fps <= 0.0 {
        return;
    }

    traj.accum += time.delta_secs();
    let step_count = (traj.accum * traj.fps).floor() as usize;
    if step_count == 0 {
        return;
    }
    traj.accum -= step_count as f32 / traj.fps;

    for _ in 0..step_count {
        step_frame(&mut traj);
    }

    apply_current_frame(
        &mut traj,
        &mut mol,
        &mut settings,
        &mut ev_changed,
        false,
        None,
    );
}

fn step_frame(traj: &mut TrajectoryState) {
    let len = traj.frames.len() as i32;
    if len == 0 {
        return;
    }
    let mut next = traj.current_frame as i32 + traj.direction;

    match traj.mode {
        PlaybackMode::Loop => {
            if next >= len {
                next = 0;
            } else if next < 0 {
                next = len - 1;
            }
        }
        PlaybackMode::Once => {
            if next >= len {
                next = len - 1;
                traj.playing = false;
            } else if next < 0 {
                next = 0;
                traj.playing = false;
            }
        }
        PlaybackMode::PingPong => {
            if next >= len {
                traj.direction = -1;
                next = len - 2;
            } else if next < 0 {
                traj.direction = 1;
                next = 1;
            }
        }
    }

    traj.current_frame = next as usize;
}

pub fn rebuild_trajectory_overlay(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut traj: ResMut<TrajectoryState>,
    settings: Res<MolSettings>,
    q_ghost_atoms: Query<Entity, With<TrajGhostAtom>>,
    q_ghost_bonds: Query<Entity, With<TrajGhostBond>>,
    q_full_atoms: Query<Entity, With<TrajFullAtom>>,
    q_full_bonds: Query<Entity, With<TrajFullBond>>,
) {
    if !traj.overlay_dirty {
        return;
    }
    traj.overlay_dirty = false;

    for e in q_ghost_atoms.iter() {
        commands.entity(e).despawn();
    }
    for e in q_ghost_bonds.iter() {
        commands.entity(e).despawn();
    }
    for e in q_full_atoms.iter() {
        commands.entity(e).despawn();
    }
    for e in q_full_bonds.iter() {
        commands.entity(e).despawn();
    }

    if traj.frames.is_empty() || !traj.overlay_enabled {
        return;
    }

    let metallic = settings.metallic.clamp(0.0, 1.0);
    let rough = settings.roughness.clamp(0.02, 1.0);
    let refl = settings.reflectance.clamp(0.0, 1.0);

    let emissive = if settings.use_emissive {
        let lr = settings.emissive_color.to_linear();
        LinearRgba::new(
            lr.red * settings.emissive_strength.max(0.0),
            lr.green * settings.emissive_strength.max(0.0),
            lr.blue * settings.emissive_strength.max(0.0),
            1.0,
        )
    } else {
        LinearRgba::BLACK
    };

    if traj.overlay_enabled && traj.overlay_full_stride > 0 {
        let stride = traj.overlay_full_stride.max(1);
        for frame in traj.frames.iter().step_by(stride) {
            spawn_overlay_frame(
                &mut commands,
                &mut meshes,
                &mut materials,
                &frame.atoms,
                &frame.pos,
                &settings,
                metallic,
                rough,
                refl,
                emissive,
                1.0,
                false,
            );
        }
    }
    if traj.overlay_enabled && traj.overlay_full_last && traj.frames.len() > 1 {
        let last = traj.frames.len() - 1;
        spawn_overlay_frame(
            &mut commands,
            &mut meshes,
            &mut materials,
            &traj.frames[last].atoms,
            &traj.frames[last].pos,
            &settings,
            metallic,
            rough,
            refl,
            emissive,
            1.0,
            false,
        );
    }
    if traj.overlay_enabled && traj.overlay_ghost_stride > 0 {
        let stride = traj.overlay_ghost_stride.max(1);
        for frame in traj.frames.iter().step_by(stride) {
            spawn_overlay_frame(
                &mut commands,
                &mut meshes,
                &mut materials,
                &frame.atoms,
                &frame.pos,
                &settings,
                metallic,
                rough,
                refl,
                emissive,
                traj.ghost_alpha,
                true,
            );
        }
    }
}

fn spawn_overlay_frame(
    commands: &mut Commands,
    meshes: &mut Assets<Mesh>,
    materials: &mut Assets<StandardMaterial>,
    atoms: &[String],
    frame: &[Vec3],
    settings: &MolSettings,
    metallic: f32,
    rough: f32,
    refl: f32,
    emissive: LinearRgba,
    alpha: f32,
    is_ghost: bool,
) {
    for (sym, pos) in atoms.iter().zip(frame.iter()) {
        let r_cov = covalent_radius_angstrom(sym) * settings.atom_scale;
        let mut color = settings
            .element_colors
            .get(sym)
            .copied()
            .unwrap_or(Color::srgb(0.7, 0.7, 0.7));
        if alpha < 1.0 {
            let s = color.to_srgba();
            color = Color::srgba(s.red, s.green, s.blue, alpha);
        }

        let sphere_mesh = Sphere::new(r_cov)
            .mesh()
            .ico(settings.atom_resolution)
            .unwrap();

        let mat = materials.add(StandardMaterial {
            base_color: color,
            metallic,
            perceptual_roughness: rough,
            reflectance: refl,
            emissive,
            alpha_mode: if alpha < 1.0 { AlphaMode::Blend } else { AlphaMode::Opaque },
            ..default()
        });

        if is_ghost {
            commands.spawn((
                Mesh3d(meshes.add(sphere_mesh)),
                MeshMaterial3d(mat),
                Transform::from_translation(*pos),
                TrajGhostAtom,
                RenderLayers::layer(LAYER_MAIN),
            ));
        } else {
            commands.spawn((
                Mesh3d(meshes.add(sphere_mesh)),
                MeshMaterial3d(mat),
                Transform::from_translation(*pos),
                TrajFullAtom,
                RenderLayers::layer(LAYER_MAIN),
            ));
        }
    }

    let bonds = compute_bonds(atoms, frame, settings.bond_thresh_scale);
    for (i, j, len_cc) in bonds {
        if len_cc <= 0.0001 {
            continue;
        }
        let p0 = frame[i];
        let p1 = frame[j];
        let dir = p1 - p0;
        let dir_n = dir.normalize();
        let ri_sphere = covalent_radius_angstrom(&atoms[i]) * settings.atom_scale;
        let rj_sphere = covalent_radius_angstrom(&atoms[j]) * settings.atom_scale;
        let base = ri_sphere.min(rj_sphere);
        let radius = (base * settings.bond_radius_pct).clamp(0.01, base);
        let rot = Quat::from_rotation_arc(Vec3::Y, dir_n);
        let center = (p0 + p1) * 0.5;
        let half_length = (len_cc * 0.5 - radius).max(0.0);
        let cap = Mesh::from(Capsule3d { radius, half_length });

        let mut color = settings.uniform_bond_color;
        if alpha < 1.0 {
            let s = color.to_srgba();
            color = Color::srgba(s.red, s.green, s.blue, alpha);
        }

        let mat = materials.add(StandardMaterial {
            base_color: color,
            metallic,
            perceptual_roughness: rough,
            reflectance: refl,
            emissive,
            alpha_mode: if alpha < 1.0 { AlphaMode::Blend } else { AlphaMode::Opaque },
            ..default()
        });

        if is_ghost {
            commands.spawn((
                Mesh3d(meshes.add(cap)),
                MeshMaterial3d(mat),
                Transform {
                    translation: center,
                    rotation: rot,
                    ..default()
                },
                TrajGhostBond,
                RenderLayers::layer(LAYER_MAIN),
            ));
        } else {
            commands.spawn((
                Mesh3d(meshes.add(cap)),
                MeshMaterial3d(mat),
                Transform {
                    translation: center,
                    rotation: rot,
                    ..default()
                },
                TrajFullBond,
                RenderLayers::layer(LAYER_MAIN),
            ));
        }
    }
}

fn compute_bonds(atoms: &[String], pos: &[Vec3], thresh_scale: f32) -> Vec<(usize, usize, f32)> {
    let n = atoms.len();
    let mut bonds = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            let ri = covalent_radius_angstrom(&atoms[i]);
            let rj = covalent_radius_angstrom(&atoms[j]);
            let cutoff = (ri + rj) * thresh_scale;
            let d = pos[i].distance(pos[j]);
            if d > 0.01 && d <= cutoff {
                bonds.push((i, j, d));
            }
        }
    }
    bonds
}

fn compute_centroid(points: &[Vec3]) -> Option<Vec3> {
    if points.is_empty() {
        return None;
    }
    let mut acc = Vec3::ZERO;
    for &p in points {
        acc += p;
    }
    Some(acc / (points.len() as f32))
}
