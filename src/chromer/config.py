import json
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Config:
    unicorns: str
    run_folder_base: str
    brain: str
    plotting_params: dict
    base_dir: Path


def _resolve_path(base: Path, value: str) -> str:
    path = Path(value)
    if path.is_absolute():
        return str(path)
    return str((base / path).resolve())


def load_config(path="config.json") -> Config:
    config_path = Path(path).resolve()
    data = json.loads(config_path.read_text())
    base = config_path.parent
    paths = data["paths"]
    return Config(
        unicorns=_resolve_path(base, paths["unicorns"]),
        run_folder_base=_resolve_path(base, paths["run_folder_base"]),
        brain=_resolve_path(base, paths["brain"]),
        plotting_params=data["plotting"]["params"],
        base_dir=base,
    )
