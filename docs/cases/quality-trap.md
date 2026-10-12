# Case: `ok=True` can hide nearly flat cells

**Same fixture and fault as [the successful case](one-trace.md)**,
but turn off barrier polish while keeping `boundary="pinned"` and
`w0=0.05`. Run `python -m fieldbook.cli run 08_settings`.

| Setting | Mean target residual | Top min SJ | New cells below 0.3 | Upstream `ok` |
|---|---:|---:|---:|---|
| Unpolished elastic | ~0.0063 | ~0.058 | 3 | True |
| Barrier elastic | ~0.0621 | ~0.627 | 0 | True |
| Stiffer elastic, `w0=2.0` | ~0.128 | ~0.825 | 0 | True |

**Verdict:** reject the first candidate under an explicitly chosen
0.3 scaled-Jacobian floor despite its much better trace fit. The
toy input had no sub-0.3 columns before morphing. A chosen floor
is a user policy, **not** an upstream guarantee. Look at both top
and bottom planes, target distance and the corridor before choosing
the barrier or stiffer candidate. Units are those of the toy deck.
