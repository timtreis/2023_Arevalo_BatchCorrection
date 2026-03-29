# TODO: Harmony v1 vs v2 comparison

## Context

The paper "Integration of large, complex single-cell datasets with Harmony2"
(Patikas, Yao, Madhu, Raychaudhuri, Hemberg & Korsunsky, bioRxiv March 2026,
DOI: 10.64898/2026.03.16.711825) describes a second-generation Harmony algorithm
that scales to >100M cells and prevents overintegration.

"Harmony2" is **not a separate package** — the improvements are shipped in the
existing `harmony` R package (≥1.0.1) and `harmonypy` Python package (≥0.1.0).
Manuscript companion repo: github.com/korsunskylab/harmonymanuscript2025.

## Version split

| | **Harmony v1** | **Harmony v2** |
|---|---|---|
| R | `harmony` ≤ 0.1.1 (Nov 2022) | `harmony` ≥ 1.0.1 (Sep 2023) |
| Python | `harmonypy` ≤ 0.0.10 (Jul 2024) | `harmonypy` ≥ 0.1.0 (Jan 2025) |

Key algorithmic changes (hardcoded, not toggleable via parameters):
- Soft-assignment formula: `(E+1)/(O+1)` → `E/(O+E)`
- Objective function: `log((O+E)/E)` with `2000/N` normalization
- Dynamic lambda estimation via new `alpha` parameter (default 0.2)
- Simplified Z_cos normalization (basic L2 instead of max-scaling)
- Removed dense BxN matrix bottleneck in regression step

**Cannot toggle v1/v2 behavior via parameters — must use different package versions.**

## Action items

- [ ] Verify `harmonypy` version inside `lightweight.sif` (`pip show harmonypy`)
- [ ] To compare v1 vs v2, create a second container with `harmonypy==0.0.10`
      (or add a `harmony_v1` method that pins the old version)
- [ ] Add `harmony_v2` as explicit method name in Snakefile if running both
- [ ] Add `alpha` to Optuna search space in `scripts/optimise_harmony.py`
      (new v2 parameter, default 0.2, controls lambda estimation)
- [ ] Consider enabling PyTorch GPU backend (`device` param in harmonypy 0.2.0)
- [ ] Re-run benchmarks for both versions and compare scores
