# Signal Detection vs Biological Interpretation

**Date:** 2026-06-29  
**Status:** Implemented

## Problem

CHROMER plot annotations mixed UV280 signal observations with biological claims. Labels such as `NO / LOW EXPRESSION`, bare fraction ranges (`A1-A2`), and `MULTI` were either over-interpretive or lacked context.

## Solution

Split chromatogram annotation into two stages:

1. **Signal detection** — `detect_uv280_signal()` finds peaks and, for pool-candidate methods, computes fraction ranges via existing `get_fraction_ranges()`.
2. **Interpretation** — `format_interpretation_label()` maps method + signal + optional `brain.json` metadata to a defensible summary string rendered by `annotate_summary()`.

## Method categories

| Constant | Methods | Fraction ranges | Default interpretation |
| --- | --- | --- | --- |
| `POOL_CANDIDATE_METHODS` | PROA, LEC | Yes | Candidate elution / binding fraction |
| `VERIFY_METHODS` | IMAC | No | Verify target-containing fractions |
| `SIGNAL_ONLY_METHODS` | SEC | No | Signal only; verify identity separately |

## Label matrix

| Condition | Label |
| --- | --- |
| No peak | `No clear UV280 elution peak` |
| PROA, 1 peak | `Candidate Protein A elution fraction: {range}` |
| PROA, 2+ peaks | `MULTI — Candidate Protein A elution fractions: {ranges}` |
| LEC, 1 peak | `Candidate lectin-binding fraction: {range}` |
| LEC, 2+ peaks | `MULTI — Candidate lectin-binding fractions: {ranges}` |
| IMAC, 1+ peaks | `IMAC UV280 peak detected: verify target-containing fractions` |
| SEC, 1+ peaks | `SEC UV280 peak detected: signal only; verify identity/oligomeric state separately` |
| Other, 1+ peaks | `UV280 peak detected: verify target-containing fractions` |

Ranges are comma-separated. Per-peak `V:` / `A:` markers are unchanged.

## Metadata upgrades

Optional `brain.json` field `TargetClass`:

- `IgG` or `Fc` on PROA runs → `Candidate IgG/Fc capture fraction(s): …`
- `glycoprotein` on LEC runs → `Candidate glycoprotein pool: …` (explicit only; no construct-name inference)

**Heuristic fallback (PROA only):** when `TargetClass` is absent, `ConstructID` is scanned for `IgG1`, ` IgG`, `_IgG`, `CHIg`, `CLIg`, `_Fc`, `huFc`.

Protein Queue.xlsx is not read at runtime.

## Plotting changes

- Removed diagonal `NO / LOW EXPRESSION` overlay.
- Replaced `annotate_no_expression`, `annotate_multiple_peaks`, and `annotate_fraction_ranges` with unified `annotate_summary()`.
- Font size 18 (16 when label length > 60 characters); `wrap=True` for long multi-peak strings.

## Data flow

`process_chrom()` returns `construct_meta` from `brain.json`. `chromeunicorns()` calls `detect_uv280_signal()` then `format_interpretation_label()` then `annotate_summary()`.

## Out of scope

- Reading Protein Queue.xlsx at runtime
- SEC oligomer assignment
- IMAC metal-type inference
- Peak-detection threshold changes

## Tests

`tests/test_interpretation_labels.py` covers method × peak scenarios, metadata upgrades, heuristics, and `detect_uv280_signal()` fraction behavior.
