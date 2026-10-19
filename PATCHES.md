# LunCoSim 0.27 convex-query fixes

This branch preserves the Parry 0.27 public API required by Avian 0.7.
It changes the geometry kernel at its owner; applications continue using the
normal native contact, intersection and point-query interfaces.

- GJK requires a separating certificate before classifying a stalled simplex
  as disjoint. Deterministic bounded perturbations resolve uncertified stalls.
  Backports: upstream `e32424193b885f4116ade337ec57b4b8d0a03152` and
  `08c43eae18537c694c4fdb4681dddc004e25b887`.
- Exactly touching EPA simplexes retain their real material witnesses and
  acquire a normal by expansion from a distinct support point. A singleton CSO
  with no distinct support point reports no normal. Adapted from upstream
  `261de840521de600580dc568068a76afebbaffac`.
- A new tetrahedron must have height above the rounding resolution of its
  support coordinates. This prevents an almost-flat simplex from being
  classified as containing the origin when all its vertices are outside it.
- Cylinder interior projection uses a deterministic nearest-cap tie break.
  Equal cap distances must not select a farther radial boundary.

`crates/parry3d-f64/tests/convex_query_contracts.rs` covers nested and separated
cylinders, cap ties, real touching witnesses in both shape orders, and a large
triangle/cylinder edge query. The application also qualifies these changes
through its continuous-collision integration and authored scene gates.

The dependency feature set is unchanged. In particular, this branch does not
enable global scalar math or change other dependency sources.
The interval trigonometry test compares endpoints using the configured simba
backend, preserving exact assertions under enhanced determinism. Its previous
host-intrinsic expectation also fails on unmodified 0.27 by one rounding bit;
the production interval implementation is unchanged.

## Thin tetrahedron separation (Roadster ramp)

GJK reduces a tetrahedron by testing its four face half-spaces and projecting
exterior native triangle features. The previous edge Voronoi-region predicates
could disagree on a thin tetrahedron and return an interior classification for
separated shapes. EPA then received a non-enclosing simplex and returned no
contact geometry. The correction retains native support witnesses and does not
retry with different seeds or synthesize a contact normal.

Regression coverage includes the recorded offset-terrain and two ramp wheel
poses, both shape orders, equivalent local translations, a prediction smaller
than the measured gap, and the exact thin simplex. The original implementation
fails the ramp regression; the corrected implementation passes it.
