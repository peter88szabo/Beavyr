const ANGSTROM_TO_BOHR: f64 = 1.0 / 0.529_177_210_903;

fn atom(x: &[f64], index: usize) -> [f64; 3] {
    let offset = 3 * index;
    [x[offset], x[offset + 1], x[offset + 2]]
}

fn largest_eigenvalue_symmetric_4(mut matrix: [[f64; 4]; 4]) -> f64 {
    for _ in 0..64 {
        let mut p = 0;
        let mut q = 1;
        let mut largest = matrix[p][q].abs();
        for i in 0..4 {
            for j in (i + 1)..4 {
                if matrix[i][j].abs() > largest {
                    largest = matrix[i][j].abs();
                    p = i;
                    q = j;
                }
            }
        }
        if largest < 1.0e-13 {
            break;
        }
        let phi = 0.5 * (2.0 * matrix[p][q]).atan2(matrix[q][q] - matrix[p][p]);
        let cosine = phi.cos();
        let sine = phi.sin();
        for k in 0..4 {
            if k == p || k == q {
                continue;
            }
            let mkp = matrix[k][p];
            let mkq = matrix[k][q];
            matrix[k][p] = cosine * mkp - sine * mkq;
            matrix[p][k] = matrix[k][p];
            matrix[k][q] = sine * mkp + cosine * mkq;
            matrix[q][k] = matrix[k][q];
        }
        let app = matrix[p][p];
        let aqq = matrix[q][q];
        let apq = matrix[p][q];
        matrix[p][p] = cosine * cosine * app - 2.0 * sine * cosine * apq + sine * sine * aqq;
        matrix[q][q] = sine * sine * app + 2.0 * sine * cosine * apq + cosine * cosine * aqq;
        matrix[p][q] = 0.0;
        matrix[q][p] = 0.0;
    }
    matrix
        .iter()
        .enumerate()
        .map(|(i, row)| row[i])
        .fold(f64::NEG_INFINITY, f64::max)
}

fn aligned_rmsd(reference: &[f64], candidate: &[f64], atoms: &[usize], mirror: bool) -> f64 {
    if atoms.is_empty() {
        return f64::INFINITY;
    }
    let n = atoms.len() as f64;
    let mut reference_center = [0.0; 3];
    let mut candidate_center = [0.0; 3];
    for &index in atoms {
        let p = atom(reference, index);
        let q = atom(candidate, index);
        for axis in 0..3 {
            reference_center[axis] += p[axis] / n;
            candidate_center[axis] += q[axis] / n;
        }
    }

    let mut covariance = [[0.0; 3]; 3];
    let mut squared_norms = 0.0;
    for &index in atoms {
        let mut p = atom(reference, index);
        let mut q = atom(candidate, index);
        for axis in 0..3 {
            p[axis] -= reference_center[axis];
            q[axis] -= candidate_center[axis];
            if mirror {
                q[axis] = -q[axis];
            }
            squared_norms += p[axis] * p[axis] + q[axis] * q[axis];
        }
        for row in 0..3 {
            for column in 0..3 {
                covariance[row][column] += p[row] * q[column];
            }
        }
    }
    let s = covariance;
    let trace = s[0][0] + s[1][1] + s[2][2];
    let horn = [
        [
            trace,
            s[1][2] - s[2][1],
            s[2][0] - s[0][2],
            s[0][1] - s[1][0],
        ],
        [
            s[1][2] - s[2][1],
            s[0][0] - s[1][1] - s[2][2],
            s[0][1] + s[1][0],
            s[0][2] + s[2][0],
        ],
        [
            s[2][0] - s[0][2],
            s[0][1] + s[1][0],
            -s[0][0] + s[1][1] - s[2][2],
            s[1][2] + s[2][1],
        ],
        [
            s[0][1] - s[1][0],
            s[0][2] + s[2][0],
            s[1][2] + s[2][1],
            -s[0][0] - s[1][1] + s[2][2],
        ],
    ];
    let best_trace = largest_eigenvalue_symmetric_4(horn);
    ((squared_norms - 2.0 * best_trace).max(0.0) / n).sqrt()
}

pub(crate) fn heavy_atom_indices(elements: &[String]) -> Vec<usize> {
    let heavy = elements
        .iter()
        .enumerate()
        .filter_map(|(index, element)| (!element.eq_ignore_ascii_case("H")).then_some(index))
        .collect::<Vec<_>>();
    if heavy.is_empty() {
        (0..elements.len()).collect()
    } else {
        heavy
    }
}

pub(crate) fn unique_against_blacklist(
    candidate: &[f64],
    blacklist: &[Vec<f64>],
    rmsd_atoms: &[usize],
    cutoff_angstrom: f64,
    chiral: bool,
) -> bool {
    let cutoff_bohr = cutoff_angstrom * ANGSTROM_TO_BOHR;
    blacklist.iter().all(|known| {
        if aligned_rmsd(known, candidate, rmsd_atoms, false) < cutoff_bohr {
            return false;
        }
        chiral || aligned_rmsd(known, candidate, rmsd_atoms, true) >= cutoff_bohr
    })
}

#[cfg(test)]
pub(crate) fn test_aligned_rmsd(reference: &[f64], candidate: &[f64], atoms: &[usize]) -> f64 {
    aligned_rmsd(reference, candidate, atoms, false)
}
