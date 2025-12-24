# Repository Guidelines

## Project Structure & Module Organization
This is a Rust/Bevy application. The crate root is `src/main.rs`, with most functionality split into Rust modules in `src/`. Key modules include `camera.rs`, `scene.rs`, `ui.rs`, and domain features like `molecule.rs`, `hbonds.rs`, and `measurements.rs`. Subdirectories such as `src/exporting/`, `src/molecule_builder/`, and `src/picking/` group related systems. The `src/old_version_that_worked/` directory holds legacy code for reference only.

## Build, Test, and Development Commands
- `cargo run` builds and launches the Bevy viewer.
- `cargo build` compiles the project without running it.
- `cargo test` runs tests (currently none are defined).
- `cargo fmt` formats code with rustfmt (recommended before commits).
- `cargo clippy` runs lint checks (recommended for warnings).

## Coding Style & Naming Conventions
- Follow Rust 2021 conventions with rustfmt defaults (4-space indentation).
- Module/file names use `snake_case` (for example, `ui_measurements.rs`).
- Types use `CamelCase`; functions, variables, and systems use `snake_case`.
- Prefer small, focused Bevy systems grouped by module to keep files readable.

## Testing Guidelines
There are no existing tests or test framework setup beyond standard Rust tooling. When adding tests, use Rust’s built-in test harness (`#[test]`) and run with `cargo test`. For integration tests, add files under `tests/` and name them after the feature (for example, `tests/measurements.rs`).

## Commit & Pull Request Guidelines
Recent Git history uses short, descriptive messages without strict prefixes. Keep commits concise and specific (for example, “add bond rotation” or “fix measurement UI”). For PRs, include a short summary, testing notes (commands and results), and screenshots or GIFs for UI changes. Link related issues when applicable.

## Configuration & Dependencies
Dependencies are managed in `Cargo.toml`; Bevy is configured with the `dynamic_linking` feature. Avoid adding new dependencies without justification and document any new runtime configuration in `README.md`.
