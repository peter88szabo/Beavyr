//! Interactive scan plot using the same order, units and reference as the table.
use super::*;
fn coordinate(report: &ScanReport, index: usize) -> f64 {
    let p = &report.points[index];
    if report.options.coordinate == Coordinate::Dihedral {
        p.target + (p.actual - p.target + 180.).rem_euclid(360.) - 180.
    } else {
        p.actual
    }
}
fn points(report: &ScanReport, view: ScanView) -> Vec<(usize, [f64; 2])> {
    view.indices(report)
        .into_iter()
        .map(|i| (i, [coordinate(report, i), view.energy(report, i)]))
        .collect()
}
#[cfg(test)]
fn hover_text(report: &ScanReport, index: usize, view: ScanView) -> String {
    format!(
        "Point {}\n{}: {:.6} {}\nΔE: {:.6} {}",
        index + 1,
        report.options.coordinate.label(),
        coordinate(report, index),
        report.options.coordinate.unit(),
        view.energy(report, index),
        view.unit.label()
    )
}
pub(super) fn draw(
    ui: &mut egui::Ui,
    report: &ScanReport,
    selected: usize,
    view: ScanView,
) -> Option<usize> {
    let values=points(report,view);
    let ids:Vec<_>=values.iter().map(|(i,_)|*i).collect();
    let xy:Vec<_>=values.iter().map(|(_,p)|Some((p[0],p[1]))).collect();
    let x=format!("{} ({})",report.options.coordinate.label(),report.options.coordinate.unit());
    let y=format!("ΔE ({})",view.unit.label());
    let reference=view.reference_index(report);
    let note=format!("Reference: point {} · E = {:.10} Eh",reference+1,report.points[reference].energy);
    let selected=values.iter().position(|(i,_)|*i==selected);
    let descending=values.first().zip(values.last()).is_some_and(|(a,b)|a.1[0]>b.1[0]);
    crate::analysis_plot::draw(ui,"scan_energy",crate::analysis_plot::Plot{points:&xy,x:&x,y:&y,note:&note,reverse_x:descending,point_ids:Some(&ids)},selected).map(|i|values[i].0)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn scan_plot_values_hover_and_order_match_display_settings() {
        let report = super::super::tests::output().scan.unwrap();
        let mut view = ScanView::default();
        assert_eq!(points(&report, view)[0].0, 0);
        assert!(hover_text(&report, 1, view).contains("kcal/mol"));
        view.unit = EnergyUnit::Kj;
        view.reference = Reference::Last;
        view.reversed = true;
        let p = points(&report, view);
        assert_eq!(p[0].0, 2);
        assert_eq!(p.last().unwrap().0, 0);
        assert_eq!(p.last().unwrap().1[1], 0.);
        assert!(hover_text(&report, 2, view).contains("kJ/mol"));
    }
}
