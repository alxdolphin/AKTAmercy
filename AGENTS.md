# AGENTS.md

## Cursor Cloud specific instructions

CHROMER (Project AKTAmercy) is a single Python 3.11 CLI (`CHROMER.py`) that parses AKTA
UNICORN chromatography exports and generates annotated chromatogram plots, syncing an
index/results to Google Sheets + Drive. There is no web server, database, or build step.

### Environment
- Dependencies live in a Python 3.11 virtualenv at `.venv` (created by `uv`; the startup
  update script provisions it). Run things with `.venv/bin/python ...` or
  `source .venv/bin/activate`.
- Pip dependencies are the single source of truth in the `pip:` block of `environment.yml`.
  Do NOT run `conda env create -f environment.yml`: its conda packages are pinned to macOS
  arm64 builds and will not resolve on Linux. The update script installs only the pip block.

### Lint / build / test / run
- Lint/build (syntax): `.venv/bin/python -m py_compile CHROMER.py "dev/debug/?INSPECT" dev/debug/parseLOG`
- Tests: none exist in this repo.
- Run the app: `./CHROMER.py` (or `.venv/bin/python CHROMER.py`) from the repo root.

### Important non-obvious caveats
- `CHROMER.py` authenticates with Google AND opens the `"Protein Queue"` and `"_CHROMATRIX_"`
  spreadsheets at IMPORT time (module-level, ~lines 240-252). So merely importing the module
  exits immediately unless a real Google service-account key exists at `./dev/.rsc/AKTAmercy.json`
  (with access to those two sheets) and network egress to Google is available. The
  `dev/debug/?INSPECT` tool imports `CHROMER`, so it has the same requirement.
- `dev/debug/parseLOG` is fully standalone (stdlib only, no Google) and runs against any
  CHROMER log file in `dev/debug/logs/` — useful for verifying the toolchain without credentials.
- A full end-to-end run also needs `.UFol`/`.Result` UNICORN input files in `data/DROP-OFF/`
  (proprietary instrument exports; none are committed to the repo).
- `config.json`'s `"mode"` key is read but never used; uploads and Sheet writes always run.
- The app expects these dirs (created on first run, but safe to pre-create):
  `data/DROP-OFF/`, `data/DONE/`, `dev/debug/logs/`, `dev/.rsc/`.
