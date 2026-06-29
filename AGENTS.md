# AGENTS.md

## Cursor Cloud specific instructions

CHROMER is an installable Python package (`src/chromer/`) that post-processes AKTA UNICORN
chromatography exports (`.UFol` / `.Result` zip archives) into annotated chromatogram
JPGs. There is no server or web UI.

### Environment
- Install with `pip install -e .` from the repo root (uses `pyproject.toml`).
- Python 3.11+ required. Dependencies: `matplotlib`, `scipy`, `seaborn`.
- Run commands with the venv interpreter, e.g. `.venv/bin/chromer process`.
- Set `MPLBACKEND=Agg` for headless rendering.

### Running the app
- `chromer process` reads inputs from `./data/DROP-OFF/` and writes JPGs to
  `./data/DONE/<timestamp>/<METHOD>/`. Filenames are `{METHOD}_{SAMPLE}.jpg` in
  generic mode, or `{METHOD}_{SAMPLE} | {ConstructID}.jpg` when indexed enrichment
  matches. Logs go to `./dev/debug/logs/`.
- `chromer inspect FILE` prints parsed metadata as JSON.
- `chromer parse-log LOG_FILE` extracts WARNING lines into a summary file.
- Optional gitignored input: `brain.json` (batch→construct index; see README
  "Construct Index"). `.Result`/`.UFol` files under `./data/DROP-OFF/` are also
  gitignored. Without export files the run does nothing useful; without
  `brain.json` CHROMER still runs in generic mode.
- `process_file()` calls `os.remove()` on each `.Result` after processing, so a
  `.Result` placed in `DROP-OFF` is consumed on the first run — re-create it to re-run.

### Tests
- `MPLBACKEND=Agg python -m unittest discover -s tests -v`

### Gotchas
- Duplicate log lines and repeated `SampleID is BLANK ... Skipping` WARNINGs are expected
  (the batch regex runs against every Run Log event line); they are not failures. Use
  `chromer parse-log` to summarize warnings from a log file.
- Paths in `config.json` resolve relative to the config file directory, not cwd.
