#![allow(dead_code)]
use anyhow::{bail, Result};
use ndarray::Array2;

/// Pretty-print a labeled square matrix in column blocks (AO-style).
///
/// - `title`: printed above the matrix
/// - `mat`: square matrix (n x n) of f64
/// - `labels`: AO labels, length n
/// - `block_size`: number of columns per block (e.g., 8)
/// - `val_width`: field width for numbers (e.g., 10)
/// - `val_prec`: digits after decimal (e.g., 5)
pub fn print_ao_matrix<L: AsRef<str>>(
    title: &str,
    mat: &Array2<f64>,
    labels: &[L],
    block_size: usize,
    val_width: usize,
    val_prec: usize,
) -> Result<()> {
    let n = labels.len();
    if n == 0 {
        bail!("labels are empty");
    }
    let (nr, nc) = mat.dim();
    if nr != n || nc != n {
        bail!("matrix must be n x n (n={n}); got {nr} x {nc}");
    }

    // Row/column label field width = max label length or val_width (so headers align nicely).
    let label_width = labels
        .iter()
        .map(|l| l.as_ref().len())
        .max()
        .unwrap_or(0)
        .max(val_width);

    println!("\n{title}");

    let mut col = 0usize;
    while col < n {
        let end = (col + block_size).min(n);

        // Column header
        print!("{:>width$}", "", width = label_width); // empty corner
        for j in col..end {
            print!(
                "{:>w$}",
                labels[j].as_ref(),
                w = val_width.max(labels[j].as_ref().len())
            );
        }
        println!();

        // Rows
        for i in 0..n {
            print!("{:>width$}", labels[i].as_ref(), width = label_width);
            for j in col..end {
                let v = mat[(i, j)];
                // Right-aligned fixed precision
                print!("{val:>w$.p$}", val = v, w = val_width, p = val_prec);
            }
            println!();
        }

        println!(); // space between blocks
        col = end;
    }

    Ok(())
}
