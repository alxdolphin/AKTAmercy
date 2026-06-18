# AGENTS.md

## Cursor Cloud specific instructions

CHROMER is a single-file Python CLI (`CHROMER.py`) that post-processes AKTA UNICORN
chromatography exports (`.UFol` / `.Result` zip archives) into annotated chromatogram
JPGs. There is no server, web UI, or test/lint suite in the repo.

### Environment
- Dependencies live in a virtualenv at `.venv` (created by the startup update script).
  The repo README suggests conda, but conda is not installed here; use `.venv` + pip.
- Run commands with the venv interpreter, e.g. `.venv/bin/python ./CHROMER.py`.
- `requirements.txt` still lists `gspread`/`oauth2client`/`pydrive` from the pre-local
  (Google Drive) era; the current code does not import them. The only runtime imports are
  `matplotlib`, `scipy`, and `seaborn` (which pull in `numpy`/`pandas`).

### Running the app
- The app is headless-friendly: set `MPLBACKEND=Agg` to render without a display.
- It reads inputs from `./data/DROP-OFF/` and writes JPGs to
  `./data/DONE/<timestamp>/<METHOD>_<BATCH> | <ConstructID>.jpg`. Logs go to
  `./dev/debug/logs/`.
- Required-but-gitignored inputs (not in the repo): `brain.json` (batch->construct index;
  see README "Construct Index" for schema) and `.Result`/`.UFol` files under
  `./data/DROP-OFF/`. Without them the run does nothing useful.
- `process_file()` calls `os.remove()` on each `.Result` after processing, so a
  `.Result` placed in `DROP-OFF` is consumed on the first run — re-create it to re-run.

### Gotchas
- Duplicate log lines and repeated `SampleID is BLANK ... Skipping` WARNINGs are expected
  (the batch regex runs against every Run Log event line); they are not failures. The
  `dev/debug/parseLOG` helper summarizes warnings from a log file.
- The `dev/debug/?INSPECT` helper changes to the repo root before importing `CHROMER`;
  relative `.Result` arguments are still resolved from the launch directory.
