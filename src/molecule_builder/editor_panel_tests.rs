use super::*;

struct EditorFixture {
    ctx: egui::Context,
    editor: EditorRotateState,
    zmat: ZMatrixBuilderState,
    mol: Molecule,
    settings: MolSettings,
    history: StructureHistory,
    traj: TrajectoryState,
}

impl EditorFixture {
    fn new(xyz: &str) -> Self {
        let mol = Molecule::from_xyz(xyz);
        let zmat = ZMatrixBuilderState {
            zmat: zmat2xyz::xyz_to_zmat(&mol.atoms, &mol.pos),
            ..Default::default()
        };
        let ctx = egui::Context::default();
        ctx.style_mut_of(ctx.theme(), |style| {
            style.text_styles = [
                (egui::TextStyle::Heading, egui::FontId::proportional(22.0)),
                (egui::TextStyle::Body, egui::FontId::proportional(16.0)),
                (egui::TextStyle::Monospace, egui::FontId::monospace(15.0)),
                (egui::TextStyle::Button, egui::FontId::proportional(16.0)),
                (egui::TextStyle::Small, egui::FontId::proportional(13.0)),
            ]
            .into();
        });
        let mut history = StructureHistory::detached();
        // Match a normal session without writing persistent history in a UI test.
        history.is_error = false;
        history.message = None;
        Self {
            ctx,
            editor: Default::default(),
            zmat,
            mol,
            settings: Default::default(),
            history,
            traj: Default::default(),
        }
    }

    fn act(&mut self, action: Action) -> Result<Option<(BuilderChange, String)>, String> {
        perform(
            action,
            &mut self.editor,
            &mut self.zmat,
            &mut self.mol,
            &mut self.settings,
            &mut self.traj,
        )
    }

    fn frame(&mut self, events: Vec<egui::Event>, size: [f32; 2]) -> egui::FullOutput {
        let input = egui::RawInput {
            screen_rect: Some(egui::Rect::from_min_size(egui::Pos2::ZERO, size.into())),
            events,
            ..Default::default()
        };
        let mut output = self.ctx.run_ui(input, |ui| {
            let draw = |ui: &mut egui::Ui| {
                builder_ui_contents(
                    ui,
                    &mut self.editor,
                    &mut self.zmat,
                    &mut self.mol,
                    &mut self.settings,
                    &mut self.history,
                    &mut self.traj,
                );
            };
            if editor_docked(ui.ctx()) {
                egui::Panel::right("test_editor_panel")
                    .default_size(420.0)
                    .min_size(380.0)
                    .show(ui, draw);
            } else {
                crate::ui::molecule_editor_window(ui.ctx(), &mut true, draw);
            }
        });
        output.textures_delta.clear();
        output
    }

    fn click(&mut self, label: &str) {
        self.frame(vec![], [1280.0, 720.0]);
        let output = self.frame(vec![], [1280.0, 720.0]);
        let position = label_rect(&output, label)
            .unwrap_or_else(|| panic!("No control labelled {label:?}"))
            .center();
        self.frame(
            vec![
                egui::Event::PointerMoved(position),
                egui::Event::PointerButton {
                    pos: position,
                    button: egui::PointerButton::Primary,
                    pressed: true,
                    modifiers: egui::Modifiers::NONE,
                },
            ],
            [1280.0, 720.0],
        );
        self.frame(
            vec![
                egui::Event::PointerMoved(position),
                egui::Event::PointerButton {
                    pos: position,
                    button: egui::PointerButton::Primary,
                    pressed: false,
                    modifiers: egui::Modifiers::NONE,
                },
            ],
            [1280.0, 720.0],
        );
    }
}

fn label_rect(output: &egui::FullOutput, label: &str) -> Option<egui::Rect> {
    fn find(shape: &egui::Shape, label: &str) -> Option<egui::Rect> {
        match shape {
            egui::Shape::Text(text) if text.galley.text() == label => {
                Some(egui::Rect::from_min_size(text.pos, text.galley.size()))
            }
            egui::Shape::Vec(shapes) => shapes.iter().find_map(|shape| find(shape, label)),
            _ => None,
        }
    }
    output
        .shapes
        .iter()
        .find_map(|shape| find(&shape.shape, label))
}

#[test]
fn named_custom_fragment_button_reuses_placement_and_undo() {
    let mut fixture = EditorFixture::new("1\ncarbon\nC 0 0 0\n");
    let source = Molecule::from_xyz("2\nhydroxyl\nO 0 0 0\nH 0.96 0 0\n");
    let fragment = crate::molecule_builder::custom_fragments::Fragment::from_selection(
        "My hydroxyl",
        &source,
        &[0, 1],
        0,
    )
    .unwrap();
    fixture.zmat.custom_fragments.entries.clear();
    fixture.zmat.custom_fragments.entries.push(fragment);
    set_selected_atom(&mut fixture.zmat, &fixture.mol.atoms, Some(0));
    fixture.click("My fragments");
    fixture.click("My hydroxyl");
    fixture.click("Add to C1");
    assert_eq!(fixture.mol.atoms, ["C", "O", "H"]);
    assert_eq!(
        fixture.zmat.placed_editor.moving.as_deref(),
        Some(&[1, 2][..])
    );
    assert!((fixture.mol.pos[1].distance(fixture.mol.pos[2]) - 0.96).abs() < 1e-5);
    fixture.act(Action::Undo).unwrap();
    assert_eq!(fixture.mol.atoms, ["C"]);
    fixture.act(Action::Redo).unwrap();
    assert_eq!(fixture.mol.atoms, ["C", "O", "H"]);
}

#[test]
fn clicking_fragment_then_add_targets_the_existing_fragment_editor() {
    let mut fixture = EditorFixture::new("1\ncarbon\nC 0 0 0\n");
    set_selected_atom(&mut fixture.zmat, &fixture.mol.atoms, Some(0));
    fixture.click("OH");
    fixture.click("Add to C1");
    assert_eq!(fixture.mol.atoms, ["C", "O", "H"]);
    let moving = fixture
        .zmat
        .placed_editor
        .moving
        .clone()
        .expect("placement must prepare the existing tool");
    assert_eq!(moving, [1, 2]);
    assert_eq!(fixture.zmat.placed_editor.picks, [0, 1, 2]);
    assert_eq!(fixture.zmat.history.depth(), 1);

    fixture.click("Fragment Editor ↗");
    assert!(fixture.editor.window_open);
    assert_eq!(fixture.editor.moving.as_deref(), Some(&[1, 2][..]));
    assert_eq!(fixture.editor.picks, [0, 1, 2]);
    assert!(
        fixture.zmat.placed_editor.moving.is_none(),
        "only one editor should own the preview"
    );

    // Use the original Fragment Editor's preparation and live geometry path.
    let mut placed = fixture.editor.clone();
    assert!(fragment_editor_prepare(
        &mut placed,
        &mut fixture.mol,
        &mut fixture.settings
    ));
    let host = fixture.mol.pos[0];
    let before = fixture.mol.pos.clone();
    placed.axis_len_target += 0.5;
    placed.angle_deg = 30.0;
    recompute_preview(&mut placed, &mut fixture.mol, &mut fixture.settings);
    assert_eq!(fixture.mol.pos[0], host);
    assert_ne!(fixture.mol.pos[1], before[1]);
    assert_ne!(fixture.mol.pos[2], before[2]);
    assert!(placed.status.is_none(), "{:?}", placed.status);
}

#[test]
fn replacement_prepares_the_new_fragment_without_reference_picking() {
    let mut fixture = EditorFixture::new("2\nCH\nC 0 0 0\nH 1.09 0 0\n");
    set_selected_atom(&mut fixture.zmat, &fixture.mol.atoms, Some(1));
    fixture.zmat.frag_name = "-OH".into();
    fixture.act(Action::Replace).unwrap();
    assert_eq!(fixture.mol.atoms, ["C", "O", "H"]);
    assert_eq!(
        fixture.zmat.placed_editor.moving.as_deref(),
        Some(&[1, 2][..])
    );
    assert_eq!(fixture.zmat.placed_editor.picks, [0, 1, 2]);
    fixture.act(Action::Undo).unwrap();
    assert_eq!(fixture.mol.atoms, ["C", "H"]);
    assert!(fixture.zmat.placed_editor.moving.is_none());
}

#[test]
fn hydrogen_completion_and_deletion_share_one_chronological_undo_history() {
    let mut fixture = EditorFixture::new("1\ncarbon\nC 0 0 0\n");
    let original = fixture.mol.pos.clone();
    fixture.act(Action::Hydrogens).unwrap();
    assert_eq!(fixture.mol.atoms, ["C", "H", "H", "H", "H"]);
    let methane = fixture.mol.pos.clone();
    fixture.act(Action::Delete(1)).unwrap();
    assert_eq!(fixture.mol.atoms.len(), 4);
    fixture.act(Action::Undo).unwrap();
    assert_eq!(fixture.mol.pos, methane);
    fixture.act(Action::Undo).unwrap();
    assert_eq!(fixture.mol.atoms, ["C"]);
    assert_eq!(fixture.mol.pos, original);
    fixture.act(Action::Redo).unwrap();
    assert_eq!(fixture.mol.pos, methane);
}

#[test]
fn failed_replacement_keeps_geometry_and_redo_history() {
    let mut fixture = EditorFixture::new("1\ncarbon\nC 0 0 0\n");
    fixture.act(Action::ChangeElement(0, "N".into())).unwrap();
    fixture.act(Action::Undo).unwrap();
    let before = fixture.mol.pos.clone();
    let redo = fixture.zmat.history.next_redo().map(str::to_owned);
    set_selected_atom(&mut fixture.zmat, &fixture.mol.atoms, Some(0));
    fixture.zmat.frag_name = "-OH".into();
    assert!(fixture.act(Action::Replace).is_err());
    assert_eq!(fixture.mol.atoms, ["C"]);
    assert_eq!(fixture.mol.pos, before);
    assert_eq!(fixture.zmat.history.next_redo(), redo.as_deref());
    assert_eq!(fixture.zmat.history.depth(), 0);
}

#[test]
fn changing_element_preserves_coordinates_and_can_be_undone() {
    let mut fixture = EditorFixture::new("2\nCH\nC 2 3 4\nH 3.09 3 4\n");
    let original = fixture.mol.pos.clone();
    fixture.act(Action::ChangeElement(1, "F".into())).unwrap();
    assert_eq!(fixture.mol.atoms, ["C", "F"]);
    assert_eq!(fixture.mol.pos, original);
    fixture.act(Action::Undo).unwrap();
    assert_eq!(fixture.mol.atoms, ["C", "H"]);
    assert_eq!(fixture.mol.pos, original);
}

#[test]
fn live_fragment_adjustment_blocks_destructive_main_panel_edits() {
    let mut fixture = EditorFixture::new("1\ncarbon\nC 0 0 0\n");
    set_selected_atom(&mut fixture.zmat, &fixture.mol.atoms, Some(0));
    fixture.zmat.frag_name = "-OH".into();
    fixture.act(Action::Add).unwrap();
    assert!(fragment_editor_prepare(
        &mut fixture.zmat.placed_editor,
        &mut fixture.mol,
        &mut fixture.settings
    ));
    fixture.zmat.placed_editor.axis_len_target += 0.5;
    recompute_preview(
        &mut fixture.zmat.placed_editor,
        &mut fixture.mol,
        &mut fixture.settings,
    );
    let preview = fixture.mol.pos.clone();
    assert!(fixture.act(Action::Delete(0)).is_err());
    assert_eq!(fixture.mol.pos, preview);
    assert_eq!(fixture.zmat.history.depth(), 1);
}

#[test]
fn editor_starts_floating_and_user_docking_preserves_geometry_history_and_fragment_target() {
    let mut fixture = EditorFixture::new("1\ncarbon\nC 0 0 0\n");
    assert!(!editor_docked(&fixture.ctx));
    set_selected_atom(&mut fixture.zmat, &fixture.mol.atoms, Some(0));
    fixture.zmat.frag_name = "-OH".into();
    fixture.act(Action::Add).unwrap();
    let positions = fixture.mol.pos.clone();
    let moving = fixture.zmat.placed_editor.moving.clone();
    fixture.click("Dock");
    assert!(editor_docked(&fixture.ctx));
    let docked = fixture.frame(vec![], [1280.0, 720.0]);
    assert!(label_rect(&docked, "Float").is_some());
    fixture.click("Float");
    assert!(!editor_docked(&fixture.ctx));
    assert_eq!(fixture.mol.pos, positions);
    assert_eq!(fixture.zmat.placed_editor.moving, moving);
    assert_eq!(fixture.zmat.history.depth(), 1);
}

#[test]
fn main_tools_stay_visible_at_laptop_height_with_all_fragments_and_placed_editor() {
    let mut fixture = EditorFixture::new("1\ncarbon\nC 0 0 0\n");
    set_selected_atom(&mut fixture.zmat, &fixture.mol.atoms, Some(0));
    fixture.zmat.frag_name = "-Phenyl".into();
    fixture.act(Action::Add).unwrap();
    let mut state = panel_state(&fixture.ctx);
    state.all_fragments = true;
    store_panel_state(&fixture.ctx, state);
    for height in [720.0, 600.0] {
        fixture.frame(vec![], [1280.0, height]);
        let output = fixture.frame(vec![], [1280.0, height]);
        for label in [
            "Undo",
            "Add missing H",
            "Clean up geometry",
            "Z-matrix ↗",
            "Fragment Editor ↗",
        ] {
            let rect = label_rect(&output, label).unwrap_or_else(|| panic!("Missing {label}"));
            assert!(
                rect.top() >= 0.0 && rect.bottom() <= height,
                "{label} is offscreen at height {height}: {rect:?}"
            );
            assert!(
                rect.right() <= 1280.0,
                "{label} overflows the editor: {rect:?}"
            );
        }
    }
}
