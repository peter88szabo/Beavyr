#![allow(dead_code)]
use std::f64::consts::PI; // Import PI from the standard library
const AV: f64 = 6.022045e0;
const PH: f64 = 6.626176e0;
const SL: f64 = 2.99792458e0;
const TOCM1: f64 = 1.0e2 * (PH * AV) / (8.0 * PI * PI * SL);
const TOMHZ: f64 = 1.0e6 * (PH * AV) / (8.0 * PI * PI);

use crate::normalmode::linalg_shim::{Backend, LinAlg};
use ndarray::Array2;

pub fn get_brot(xyz: &Vec<[f64; 3]>, mass: &Vec<f64>) -> [f64; 3] {
    let mut ixx = 0.0;
    let mut iyy = 0.0;
    let mut izz = 0.0;
    let mut ixy = 0.0;
    let mut ixz = 0.0;
    let mut iyz = 0.0;

    // Calculate elements of the inertia tensor
    for (k, &mass) in mass.iter().enumerate() {
        let x = xyz[k][0];
        let y = xyz[k][1];
        let z = xyz[k][2];

        let x2 = x * x;
        let y2 = y * y;
        let z2 = z * z;

        ixx += mass * (y2 + z2);
        iyy += mass * (x2 + z2);
        izz += mass * (x2 + y2);
        ixy += mass * x * y;
        ixz += mass * x * z;
        iyz += mass * y * z;
    }

    let la = LinAlg::new(Backend::Auto).expect("inertia: failed to init linalg");
    let mut it = Array2::<f64>::zeros((3, 3));
    it[(0, 0)] = ixx;
    it[(0, 1)] = ixy;
    it[(0, 2)] = ixz;
    it[(1, 0)] = ixy;
    it[(1, 1)] = iyy;
    it[(1, 2)] = iyz;
    it[(2, 0)] = ixz;
    it[(2, 1)] = iyz;
    it[(2, 2)] = izz;
    let (eigval, _eigvec) = la.eigh(&it).expect("inertia: eigh failed");

    let mut inertia_amuang2 = [0.0; 3]; // Initialize with a size of 3
    let mut brot_cm1 = [0.0; 3]; // Initialize with a size of 3
    let mut brot_mhz = [0.0; 3]; // Initialize with a size of 3

    for i in 0..3 {
        brot_cm1[i] = TOCM1 / eigval[i];
        brot_mhz[i] = TOMHZ / eigval[i];
        inertia_amuang2[i] = eigval[i] / AV * 10.0;
    }

    //if you need brot or intertia in other units, then just print them
    return brot_cm1;
}
