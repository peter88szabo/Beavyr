# Vendored `dreid-typer`

DREIDING atom typing and molecular topology perception. **Copied into Beavyr's
source tree, not declared as a dependency**, so the project builds from its own
checkout with nothing to fetch or link.

## Upstream

| | |
|---|---|
| project | `dreid-typer` v0.4.0 |
| upstream | https://github.com/caltechmsc/dreid-typer |
| fork used | https://github.com/peter88szabo/dreid-typer |
| licence | MIT -- see `LICENSE-MIT` |
| copyright | © 2025 California Institute of Technology, Materials and Process Simulation Center (MSC) |
| authors | Tony Kan, William A. Goddard III |

Goddard is a co-author of the DREIDING paper itself (Mayo, Olafson & Goddard,
*J. Phys. Chem.* **1990**, *94*, 8897), which is the reason this code was
adopted rather than written fresh: typing failures are silent, producing a
plausible but wrong geometry rather than an error.

## Changes from upstream

Deliberately minimal, and confined to removing the three dependencies. Beavyr
has five dependencies on purpose, and needed none of these.

The entanglement was shallow: `serde`, `toml` and `thiserror` appeared **only**
in rule loading. The perception engine -- `aromaticity`, `electrons`,
`hybridization`, `kekulize`, `resonance`, `rings`, ~4,600 lines and the
substance of the library -- used none of them and is byte-for-byte unchanged
apart from module paths.

1. **`typing/rules_data.rs` is new, and replaces TOML parsing.**
   Upstream embeds `resources/default.rules.toml` via `include_str!` and parses
   it at runtime with `serde` + `toml`. That 408-line file, 62 rules, is
   translated here into a static Rust array built once behind a `OnceLock`.
   `get_default_rules()` keeps its signature.

   The translation was produced **mechanically** by a throwaway generator that
   parsed the TOML and emitted the literals, asserting that every rule key and
   every condition key it encountered was handled -- so a schema field could not
   be silently dropped. Fidelity is then *proved* by the ported upstream test
   suite: if the deck differs in any way that matters, `dreiding_paper.rs`
   fails.

2. **`typing/rules.rs` keeps only the schema.** `Rule` and `Conditions` stay
   public, with their `Deserialize` derives and `#[serde(...)]` attributes
   removed; `Ruleset`, `parse_rules`, `deserialize_str_keyed_map` and the
   `OnceLock` cache are gone. Its TOML-parsing unit tests are replaced by tests
   over the generated deck: the rule count, that the cache is not rebuilt, that
   no rule has empty conditions (which would shadow the whole deck), and that no
   two rules share both a priority and a condition set (which would make typing
   depend on deck order rather than on chemistry).

   A custom deck is still possible via `assign_topology_with_rules`; it is built
   in Rust rather than parsed from a file.

3. **`core/error.rs`:** `thiserror` derives replaced by hand-written `Display`,
   `Error::source` and `From` impls. Messages are unchanged. The
   `TyperError::RuleParse(toml::de::Error)` variant is removed, there being no
   TOML parser.

4. **`core/properties.rs`:** `thiserror` derives on the three parse errors
   replaced by hand-written impls, messages unchanged; the `Deserialize` impls
   for `Element` and `Hybridization` removed. Those impls only delegated to
   `FromStr`, which is untouched and still public.

5. **`lib.rs` became `mod.rs`**, and `crate::{core,perception,typing,builder}::`
   paths were rewritten to `crate::forcefield::dreiding::typer::...`. The
   doctest on the entry point is marked `ignore` (Beavyr is a binary crate, so
   doctests do not run) and rewritten against water rather than ethanol.

Nothing else was touched. No perception logic, no rule content, no atom type.

## Tests

The upstream suite is ported into Beavyr's tree under
`src/forcefield/dreiding/typer/tests.rs` -- `dreiding_paper.rs` (the molecules
from the paper), `amino_acids.rs`, `nucleic_acids.rs` and the harness. Those
tests are what carry the "validated" property across the vendoring; without
them, the rule translation in (1) would be an assumption rather than a fact.

## Updating

Re-vendor from upstream, then re-apply the five changes above. The generator
for (1) lives in the design spec's history rather than in the tree, since it
needs `toml` and must not be a Beavyr dependency; it is a dozen lines and
re-deriving it is easier than maintaining it.

See `docs/superpowers/specs/2026-09-09-dreiding-force-field-design.md`.
