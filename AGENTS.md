# AGENTS.md

## Cursor Cloud specific instructions

CHROMER is a single-file Python CLI (`CHROMER.py`) that post-processes AKTA UNICORN
chromatography exports (`.UFol` / `.Result` zip archives) into annotated chromatogram
JPGs. There is no server, web UI, or test/lint suite in the repo.

### Environment
- Dependencies live in a virtualenv at `.venv` (created by the startup update script).
  The repo README suggests conda, but conda is not installed here; use `.venv` + pip.
- Run commands with the venv interpreter, e.g. `.venv/bin/python ./CHROMER.py`.
- The only runtime dependencies in `requirements.txt` are `matplotlib`, `scipy`, and
  `seaborn` (which pull in `numpy`/`pandas`).

### Running the app
- The app is headless-friendly: set `MPLBACKEND=Agg` to render without a display.
- It reads inputs from `./data/DROP-OFF/` and writes JPGs to
  `./data/DONE/<timestamp>/<METHOD>/`. Filenames are `{METHOD}_{SAMPLE}.jpg` in
  generic mode, or `{METHOD}_{SAMPLE} | {ConstructID}.jpg` when indexed enrichment
  matches. Logs go to `./dev/debug/logs/`.
- Optional gitignored input: `brain.json` (batch→construct index; see README
  "Construct Index"). `.Result`/`.UFol` files under `./data/DROP-OFF/` are also
  gitignored. Without export files the run does nothing useful; without
  `brain.json` CHROMER still runs in generic mode.
- `process_file()` calls `os.remove()` on each `.Result` after processing, so a
  `.Result` placed in `DROP-OFF` is consumed on the first run — re-create it to re-run.

### Gotchas
- Duplicate log lines and repeated `SampleID is BLANK ... Skipping` WARNINGs are expected
  (the batch regex runs against every Run Log event line); they are not failures. The
  `dev/debug/parseLOG` helper summarizes warnings from a log file.
- The `dev/debug/?INSPECT` helper changes to the repo root before importing `CHROMER`;
  relative `.Result` arguments are still resolved from the launch directory.
