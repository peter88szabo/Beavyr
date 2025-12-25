use bevy::prelude::Vec3;

#[derive(Clone, Copy, Debug)]
struct Vec3d {
    x: f64,
    y: f64,
    z: f64,
}

impl Vec3d {
    fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    fn add(self, other: Vec3d) -> Vec3d {
        Vec3d::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }

    fn sub(self, other: Vec3d) -> Vec3d {
        Vec3d::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }

    fn scale(self, s: f64) -> Vec3d {
        Vec3d::new(self.x * s, self.y * s, self.z * s)
    }

    fn dot(self, other: Vec3d) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    fn cross(self, other: Vec3d) -> Vec3d {
        Vec3d::new(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )
    }

    fn norm(self) -> f64 {
        self.dot(self).sqrt()
    }

    fn normalize(self) -> Vec3d {
        let n = self.norm();
        if n == 0.0 {
            self
        } else {
            self.scale(1.0 / n)
        }
    }
}

impl From<Vec3> for Vec3d {
    fn from(v: Vec3) -> Self {
        Self::new(v.x as f64, v.y as f64, v.z as f64)
    }
}

impl From<Vec3d> for Vec3 {
    fn from(v: Vec3d) -> Self {
        Vec3::new(v.x as f32, v.y as f32, v.z as f32)
    }
}

#[derive(Clone, Debug)]
pub struct ZAtom {
    pub symbol: String,
    pub bond_ref: Option<usize>,
    pub bond_len: f64,
    pub angle_ref: Option<usize>,
    pub angle_deg: f64,
    pub dihedral_ref: Option<usize>,
    pub dihedral_deg: f64,
}

fn clamp_cos(c: f64) -> f64 {
    if c > 1.0 {
        1.0
    } else if c < -1.0 {
        -1.0
    } else {
        c
    }
}

fn calc_angle(a: Vec3d, b: Vec3d, c: Vec3d) -> f64 {
    let ab = b.sub(a);
    let cb = b.sub(c);
    let cos_angle = clamp_cos(ab.dot(cb) / (ab.norm() * cb.norm()));
    cos_angle.acos().to_degrees()
}

fn pick_perpendicular(v: Vec3d) -> Vec3d {
    let axis = if v.x.abs() < 0.9 {
        Vec3d::new(1.0, 0.0, 0.0)
    } else {
        Vec3d::new(0.0, 1.0, 0.0)
    };
    v.cross(axis).normalize()
}

fn calc_dihedral_zmat(p1: Vec3d, p2: Vec3d, p3: Vec3d, p4: Vec3d) -> f64 {
    let e1 = p1.sub(p2).normalize();
    let e2 = p2.sub(p3).normalize();
    let mut n = e1.cross(e2);
    if n.norm() < 1.0e-12 {
        n = pick_perpendicular(e1);
    } else {
        n = n.normalize();
    }
    let m = n.cross(e1);
    let v = p4.sub(p1);
    v.dot(n).atan2(v.dot(m)).to_degrees()
}

pub fn xyz_to_zmat(symbols: &[String], coords: &[Vec3]) -> Vec<ZAtom> {
    let mut zmat = Vec::with_capacity(coords.len());
    let coords_d: Vec<Vec3d> = coords.iter().copied().map(Vec3d::from).collect();

    for i in 0..coords_d.len() {
        let symbol = symbols[i].clone();
        let entry = if i == 0 {
            ZAtom {
                symbol,
                bond_ref: None,
                bond_len: 0.0,
                angle_ref: None,
                angle_deg: 0.0,
                dihedral_ref: None,
                dihedral_deg: 0.0,
            }
        } else if i == 1 {
            let bond = coords_d[i].sub(coords_d[0]).norm();
            ZAtom {
                symbol,
                bond_ref: Some(1),
                bond_len: bond,
                angle_ref: None,
                angle_deg: 0.0,
                dihedral_ref: None,
                dihedral_deg: 0.0,
            }
        } else if i == 2 {
            let bond = coords_d[i].sub(coords_d[1]).norm();
            let angle = calc_angle(coords_d[i], coords_d[1], coords_d[0]);
            ZAtom {
                symbol,
                bond_ref: Some(2),
                bond_len: bond,
                angle_ref: Some(1),
                angle_deg: angle,
                dihedral_ref: None,
                dihedral_deg: 0.0,
            }
        } else {
            let bond = coords_d[i].sub(coords_d[i - 1]).norm();
            let angle = calc_angle(coords_d[i], coords_d[i - 1], coords_d[i - 2]);
            let dihedral = calc_dihedral_zmat(
                coords_d[i - 1],
                coords_d[i - 2],
                coords_d[i - 3],
                coords_d[i],
            );
            ZAtom {
                symbol,
                bond_ref: Some(i),
                bond_len: bond,
                angle_ref: Some(i - 1),
                angle_deg: angle,
                dihedral_ref: Some(i - 2),
                dihedral_deg: dihedral,
            }
        };
        zmat.push(entry);
    }

    zmat
}

pub fn zmat_to_xyz(zmat: &[ZAtom]) -> Vec<Vec3> {
    let mut coords: Vec<Vec3d> = Vec::with_capacity(zmat.len());

    for (i, atom) in zmat.iter().enumerate() {
        let clamp_ref = |value: Option<usize>, max_ref: usize| -> usize {
            let v = value.unwrap_or(1).max(1).min(max_ref);
            v - 1
        };
        let pos = match i {
            0 => Vec3d::new(0.0, 0.0, 0.0),
            1 => {
                let r1 = clamp_ref(atom.bond_ref, i);
                let bond = atom.bond_len;
                let base = coords[r1];
                base.add(Vec3d::new(bond, 0.0, 0.0))
            }
            2 => {
                let r1 = clamp_ref(atom.bond_ref, i);
                let r2 = clamp_ref(atom.angle_ref, i);
                let bond = atom.bond_len;
                let theta = atom.angle_deg.to_radians();

                let p1 = coords[r1];
                let p2 = coords[r2];
                let e1 = p1.sub(p2).normalize();
                let e2 = pick_perpendicular(e1);
                let e3 = e2.cross(e1).normalize();

                p1.add(e1.scale(-bond * theta.cos()))
                    .add(e3.scale(bond * theta.sin()))
            }
            _ => {
                let r1 = clamp_ref(atom.bond_ref, i);
                let r2 = clamp_ref(atom.angle_ref, i);
                let r3 = clamp_ref(atom.dihedral_ref, i);
                let bond = atom.bond_len;
                let theta = atom.angle_deg.to_radians();
                let phi = atom.dihedral_deg.to_radians();

                let p1 = coords[r1];
                let p2 = coords[r2];
                let p3 = coords[r3];

                let e1 = p1.sub(p2).normalize();
                let e2 = p2.sub(p3).normalize();
                let mut n = e1.cross(e2);
                if n.norm() < 1.0e-12 {
                    n = pick_perpendicular(e1);
                } else {
                    n = n.normalize();
                }
                let m = n.cross(e1);

                p1.add(e1.scale(-bond * theta.cos()))
                    .add(m.scale(bond * theta.sin() * phi.cos()))
                    .add(n.scale(bond * theta.sin() * phi.sin()))
            }
        };
        coords.push(pos);
    }

    coords.into_iter().map(Vec3::from).collect()
}
