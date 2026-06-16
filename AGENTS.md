# AGENTS.md

## Cursor Cloud specific instructions

CHROMER (Project AKTAmercy) is a single Python 3.11 CLI (`CHROMER.py`) that parses AKTA
UNICORN chromatography exports and renders annotated chromatogram plots. There is no web
server, database, or build step.

### Scope for local development
- Focus is **local compilation** and **local processing of `.UFol`/`.Result` input files**.
- The **Google layer (Sheets/Drive upload) is out of scope** for local dev.
- `.UFol`/`.Result` inputs are **in scope to process but are proprietary — never commit them**.
  Drop them in `data/DROP-OFF/` (git-ignored); they are not stored in the repo.

### Environment
- Dependencies live in a Python 3.11 virtualenv at `.venv` (created by `uv`; the startup
  update script provisions it). Run things with `.venv/bin/python ...` or `source .venv/bin/activate`.
- Pip dependencies are the single source of truth in the `pip:` block of `environment.yml`.
  Do NOT run `conda env create -f environment.yml`: its conda packages are pinned to macOS
  arm64 builds and will not resolve on Linux. The update script installs only the pip block.

### Lint / build / test / run
- Local compilation (the "build"): `.venv/bin/python -m compileall CHROMER.py` and
  `.venv/bin/python -m py_compile CHROMER.py "dev/debug/?INSPECT" dev/debug/parseLOG`.
- Tests: none exist in this repo.
- Run the app: `./CHROMER.py` (or `.venv/bin/python CHROMER.py`) from the repo root.
- `dev/debug/parseLOG <logfile>` is fully standalone (stdlib only, no Google) — runs against
  any CHROMER log in `dev/debug/logs/`.

### Important non-obvious caveat (Google auth at import time)
`CHROMER.py` authenticates with Google AND opens the `"Protein Queue"` and `"_CHROMATRIX_"`
spreadsheets at IMPORT time (module-level, ~lines 240-252), and `save_and_upload_plot()`
calls `upload_file()` on every plot. So `./CHROMER.py` (and `dev/debug/?INSPECT`, which
imports `CHROMER`) cannot run without a real Google service-account key at
`./dev/.rsc/AKTAmercy.json` plus network egress.

To exercise the **local** parse/plot pipeline without Google (the in-scope path), inject stub
modules into `sys.modules` **before** importing `CHROMER`, so the module-level Google calls
become no-ops, then call `CHROMER.process_chrom()` / `CHROMER.chromeunicorns()` directly on a
`.Result` in `data/DROP-OFF/`. Set `MPLBACKEND=Agg` for headless plotting. The modules to stub
are `gspread`, `pydrive.auth`, `pydrive.drive`, and `oauth2client.service_account`. This needs
a `brain` dict mapping the run's batch/SampleID to a `ConstructID` (normally synced from Sheets
into `data/brain.json`).

### Misc
- `config.json`'s `"mode"` key is read but never used; uploads/Sheet writes always run in `CHROMER.py`.
- Runtime dirs (safe to pre-create): `data/DROP-OFF/`, `data/DONE/`, `dev/debug/logs/`, `dev/.rsc/`.
