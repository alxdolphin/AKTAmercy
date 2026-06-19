import argparse
import datetime
import json
import logging
import os
import sys
import tarfile
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

from chromer.config import load_config
from chromer.logging_utils import setup_logger, summarize_warnings
from chromer.metadata import enable_cognition, inspect_summary, process_chrom
from chromer.plotting import apply_plotting_params, chromeunicorns


def _is_within_directory(directory, target):
    abs_directory = os.path.abspath(directory)
    abs_target = os.path.abspath(target)
    return os.path.commonpath([abs_directory, abs_target]) == abs_directory


def extractall_compat(tar, path):
    try:
        tar.extractall(path=path, filter='data')
    except TypeError:
        for member in tar.getmembers():
            member_path = os.path.join(path, member.name)
            if not _is_within_directory(path, member_path):
                raise ValueError(f"Unsafe tar member path: {member.name}")
        tar.extractall(path=path)


def process_result_file(result_file, run_folder, brain):
    logging.info(f"CHROMER: Processing {os.path.basename(result_file)}")
    udata_info = process_chrom(result_file, brain)
    if not udata_info:
        return

    method = udata_info['method']
    batch = udata_info['batch']
    title = udata_info['title']

    if batch is None or method is None or batch == "" or method == "" or "None" in title:
        logging.warning(
            f"CHROMER: Batch# or Method is missing or invalid. Skipping... {os.path.basename(result_file)}"
        )
        return

    method_folder = os.path.join(run_folder, method)
    os.makedirs(method_folder, exist_ok=True)
    processed = chromeunicorns(result_file, udata_info, method_folder)
    if processed:
        os.remove(result_file)


def process_file(root, run_folder, brain, file):
    if file.endswith(".Result"):
        process_result_file(os.path.join(root, file), run_folder, brain)
    elif file.endswith(".UFol"):
        tar_file = os.path.join(root, file)
        with tarfile.open(tar_file, 'r') as tar:
            extractall_compat(tar, run_folder)
            for walk_root, _dirs, files in os.walk(run_folder):
                for nested_file in files:
                    if nested_file.endswith(".Result"):
                        process_result_file(os.path.join(walk_root, nested_file), run_folder, brain)


def cmd_process(args):
    config = load_config(args.config)
    run_folder = os.path.join(
        config.run_folder_base,
        datetime.datetime.now().strftime("%Y%m%d_%H%M%S/"),
    )
    os.makedirs(run_folder, exist_ok=True)

    setup_logger()
    apply_plotting_params(config.plotting_params)

    brain = enable_cognition(config.brain)
    if brain:
        logging.info("CHROMER: Indexed mode (brain.json loaded).")
    else:
        logging.info("CHROMER: Generic mode (no brain.json). Construct enrichment disabled.")

    all_files = []
    processed_files = set()
    for root, _dirs, files in os.walk(config.unicorns):
        for file in files:
            if file.endswith(".UFol") or file.endswith(".Result"):
                full_path = os.path.join(root, file)
                if full_path not in processed_files:
                    all_files.append((root, run_folder, brain, file))
                    processed_files.add(full_path)

    for file_args in all_files:
        process_file(*file_args)

    logging.info(
        "\n\nCHROMER: Finished processing %d files at %s",
        len(all_files),
        datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
    )


def cmd_inspect(args):
    config = load_config(args.config)
    zip_file = args.file
    if not os.path.isabs(zip_file):
        zip_file = str(Path.cwd() / zip_file)

    logging.basicConfig(level=logging.INFO)
    brain = enable_cognition(config.brain)
    result = process_chrom(zip_file, brain)
    print(json.dumps(inspect_summary(result), indent=2))


def cmd_parse_log(args):
    print(summarize_warnings(args.log_file))


def main(argv=None):
    parser = argparse.ArgumentParser(prog="chromer")
    sub = parser.add_subparsers(dest="command", required=True)

    process_parser = sub.add_parser("process", help="Process DROP-OFF exports into JPG chromatograms")
    process_parser.add_argument(
        "--config",
        default="config.json",
        help="Path to config.json (default: config.json)",
    )
    process_parser.set_defaults(func=cmd_process)

    inspect_parser = sub.add_parser("inspect", help="Inspect a .Result file and print parsed metadata")
    inspect_parser.add_argument("file", help="Path to a .Result file")
    inspect_parser.add_argument(
        "--config",
        default="config.json",
        help="Path to config.json (default: config.json)",
    )
    inspect_parser.set_defaults(func=cmd_inspect)

    parse_log_parser = sub.add_parser("parse-log", help="Extract and summarize warnings from a log file")
    parse_log_parser.add_argument("log_file", help="Path to a chromatics log file")
    parse_log_parser.set_defaults(func=cmd_parse_log)

    args = parser.parse_args(argv)
    args.func(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
