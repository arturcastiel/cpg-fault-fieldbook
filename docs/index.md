# CPG + Fault Fieldbook

Two connected libraries, two kinds of question:

- **How does this work?** Follow the [numbered lessons](learn/grids.md)
  from a five-by-five toy deck to fault projection, adaptation and four
  independent checks.
- **How exactly do I call a method?** Use the
  [MATLAB-style function reference](functions/index.md): syntax, array
  shapes, defaults, variants, outputs and failure modes.
- **Which option should I choose?** Start with the
  [conditional settings guide](decisions/settings.md), then reproduce the
  comparisons rather than trusting a historical winning score.
- **What does a pass or failure look like?** Read the
  [passing case](cases/one-trace.md), [quality trap](cases/quality-trap.md)
  and [sheared-outline trap](cases/outline-trap.md).

`cpg2python` provides the `cpg2unstructured` package for reading,
topology, morphology, validity and writing grid geometry.
`fault-projections` provides `fault_projections` for authored fault
geometry, face projection, targets, adaptation and fault output.

```text
authored trace/polygon → targets (fault_projections)
                             ↓
GRDECL/EGRID → pillar grid → morph_pillars (cpg2unstructured)
                             ↓
       validity + quality + corridor + chain + rim
                             ↓
                   optional grid INCLUDE
```

**Start:** [install](start/install-linux.md) → [read a grid](learn/grids.md)
→ [project a fault](learn/faults.md) → [deform](learn/deformation.md)
→ [judge and export](learn/validation.md) → [slide a fault in Z](learn/z-slide.md).

The [previous standalone HTML wiki](legacy-field-guide.html) remains an
additional dated snapshot. Its file links assume the two upstream
checkouts are siblings of this repository; for online browsing, use
the linked upstream repositories instead. No workflow here proves
simulator acceptance or safe flow/transmissibility behavior.
