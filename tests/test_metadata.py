import json
import tempfile
import unittest
from pathlib import Path

from chromer.config import load_config
from chromer.metadata import format_title


class FormatTitleTests(unittest.TestCase):
    def test_with_construct(self):
        self.assertEqual("SEC_B01 | Foo-Bar", format_title("SEC", "B01", "Foo-Bar"))

    def test_without_construct(self):
        self.assertEqual("IMAC_B01", format_title("IMAC", "B01"))


class ConfigTests(unittest.TestCase):
    def test_resolves_relative_paths_from_config_dir(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            dropoff = root / "inputs"
            done = root / "outputs"
            brain = root / "brain.json"
            dropoff.mkdir()
            done.mkdir()
            brain.write_text("{}")

            config_path = root / "config.json"
            config_path.write_text(json.dumps({
                "paths": {
                    "unicorns": "./inputs/",
                    "run_folder_base": "./outputs/",
                    "brain": "./brain.json",
                },
                "plotting": {"params": {"axes.grid": True}},
            }))

            config = load_config(config_path)

            self.assertEqual(str(dropoff.resolve()), config.unicorns)
            self.assertEqual(str(done.resolve()), config.run_folder_base)
            self.assertEqual(str(brain.resolve()), config.brain)
            self.assertTrue(config.plotting_params["axes.grid"])


if __name__ == "__main__":
    unittest.main()
