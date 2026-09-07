# Rendering Roadmap

How Beavyr's atoms and bonds could look better, and what the engine already
offers for it. Written 2026-09-07 against Bevy 0.19.1.

Every API named here was checked against the crate sources in
`~/.cargo/registry/src/*/bevy_*-0.19*`, not recalled — the feature set moved a
lot between Bevy versions, so re-verify before acting on this if the engine has
been upgraded since.

## Why the materials feel unfinished

The main camera is spawned as a bare `Camera3d::default()` (`src/scene.rs:379`)
with nothing attached to it. Four consequences, in order of how much they cost
the picture:

### 1. No environment map, so two material sliders cannot work

`MolSettings` exposes `metallic`, `roughness` and `reflectance`
(`src/settings.rs:84-86`), and they are handed to `StandardMaterial`
(`src/scene.rs:350`). But with no `EnvironmentMapLight` on the camera there is
nothing for a specular surface to reflect: `metallic > 0` makes an atom go dark
rather than shiny, and `reflectance` barely registers.

This is not an incomplete feature. It is two controls that are physically
unable to do their job until the scene has an environment to reflect. It is the
first thing to fix and the largest single improvement available.

### 2. Default tonemapping desaturates the element colours

Bevy defaults to `Tonemapping::TonyMcMapface`, whose own documentation says
"brights desaturate across the spectrum". That is why CPK colours read flatter
here than in Molden or PyMOL. `Tonemapping::None` gives literal sRGB.

### 3. The light rig is world-fixed and distance-dependent

Three `PointLight`s at fixed world positions with `shadow_maps_enabled: false`
and `range: 200` (`src/scene.rs:388-433`):

* Point lights fall off with distance, so the same `light_intensity` and
  `light_distance` expose a water molecule and a 200-atom system differently --
  the sliders need re-tuning per system, which they should not.
* Being fixed in world space, orbiting the camera swings the molecule into
  gloom. Every molecular viewer attaches at least a headlight to the camera.

### 4. No ambient occlusion

The single biggest legibility gain available, and no amount of light tuning
substitutes for it: AO is what makes the crevices where atoms touch read as
depth rather than as overlapping flat discs.

## Tier 1 -- hours of work, large payoff

| change | effect |
|---|---|
| `EnvironmentMapLight` on the main camera | makes glossy and metallic atoms work at all. `GeneratedEnvironmentMapLight` also exists in 0.19 (`bevy_light-0.19.1/src/probe.rs`), which may avoid shipping an HDRI asset entirely -- check its requirements first, since keeping the dependency and asset count low is deliberate here |
| Parent the key/fill/rim rig to the camera | lighting stops depending on the orbit angle |
| A `Tonemapping` selector in Appearance | `None`, `Reinhard`, `ReinhardLuminance`, `AcesFitted`, `AgX`, `TonyMcMapface`. Fixes colour fidelity for one dropdown. `AgX` needs the `tonemapping_luts` cargo feature |
| `DistanceFog` for depth cueing | one component, and the classic VMD trick. Very effective on anything large |

## Tier 2 -- a real quality jump

* **`ScreenSpaceAmbientOcclusion`** (`bevy_pbr::ssao`). Requires a
  `DepthPrepass`.
* **`ContactShadows`** (`bevy_pbr/src/contact_shadows.rs`, new in this Bevy).
  Its own doc: "small-scale shadows in areas where traditional shadow maps may
  lack detail, such as where objects touch the ground" -- which is exactly the
  ball-and-stick case. Screen-space, needs no shadow maps, and is declared
  `#[require(DepthPrepass)]`.

  **Do these two together.** They share the depth prepass, so the second is
  nearly free once the first has paid for it.

* **Swap `PointLight` for `DirectionalLight`** on the key/fill/rim rig, for
  scale invariance, and enable shadows on the key light only.
* **An anti-aliasing selector.** `bevy_anti_alias-0.19.0` ships `fxaa`, `smaa`,
  `taa`, `dlss` and contrast-adaptive sharpening. TAA is markedly better than
  MSAA for thin bonds and for the isosurface wireframe; MSAA remains better for
  a still frame with no camera motion, so both are worth offering.
* **`OrderIndependentTransparencySettings`** (`bevy_core_pipeline::oit`). This
  one is a correctness issue rather than polish: the Surface tool already has an
  opacity slider, so transparency is currently sorted per object, which will
  produce visible artefacts wherever a translucent isosurface overlaps atoms.

## Tier 3 -- a signature look, and scale

* **Silhouette outlines.** Bevy has no built-in; it needs a post-process edge
  pass over the depth and normal buffers, or the inverted-hull trick. This is
  the illustrative "Goodsell / Mol\*" look -- very legible, and the thing that
  would make a Beavyr figure recognisable at a glance.
* **Sphere and cylinder impostors.** A custom `Material` whose fragment shader
  writes depth, drawing each atom as one camera-facing quad. Gives mathematically
  perfect spheres at any zoom -- the current `Sphere::new(r).mesh().ico(resolution)`
  (`src/scene.rs:309`) facets when zoomed in -- and is the route to 100k atoms.
  It would also make the atom-resolution slider unnecessary. The largest job
  here, and the right answer for big trajectories.
* **Supersampled export.** Render the Export tool's output at 3-4x and
  downsample. A cheap way to close part of the gap to PyMOL for publication
  figures, using the export path that already exists.

## Deliberately not worth it here

`atmosphere`, `volumetric_fog`, `ssr` (a molecule has nothing to reflect),
`meshlet` (the geometry is not dense enough to need virtual geometry), motion
blur, auto exposure, `decal`.

## Where this sits against other viewers

Beavyr already looks better than the tools chemists actually use for
quantum-chemistry output -- Molden, GabEdit, and Avogadro's older renderer --
essentially for free, because Bevy/wgpu gives it PBR materials, MSAA and 60 fps
on a large trajectory without effort spent.

It does not look as good as the tools built for figures: PyMOL's ray-traced
output remains the publication standard, and Mol\*/NGL are handsome. Tier 1 and
Tier 2 would close most of the interactive gap; only the outline work and
supersampled export address the figure-quality gap, and neither makes Beavyr a
renderer, which is not what it is for.
