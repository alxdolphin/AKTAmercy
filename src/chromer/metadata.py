import json
import logging
import os
import re

from chromer.unicorn import pc_uni7


def enable_cognition(file_path):
    try:
        with open(file_path) as f:
            brain = json.load(f)
        logging.info("CHROMER: Brain loaded successfully. Chromatogram recognition enabled.")
        return brain
    except Exception as e:
        logging.error(f"CHROMER: Brain empty. Unable to recognize constructs.\n    REASON: {e}")
        return None


_FILENAME_RE = re.compile(r'[bB]#([^-\s]*)', re.IGNORECASE)


def _resolve_sample_from_filename(zip_file):
    filename = os.path.basename(zip_file)
    stem, _ = os.path.splitext(filename)
    filename_match = _FILENAME_RE.search(stem) or _FILENAME_RE.search(filename)
    if filename_match:
        return filename_match.group(1).upper()
    return re.sub(r'[^\w.-]', '_', stem) or None


def mine_run_metadata(udata, zip_file):
    method = sample = date = None
    batch_re = re.compile(r'Set mark "([^-\s]*)"', re.IGNORECASE)
    method_re = re.compile(r'Method:(.*)', re.IGNORECASE)
    date_re = re.compile(r'Method Run ((\d{1,2}/\d{1,2}/\d{4}))')

    for key, value in udata.items():
        if key == "Run Log":
            for data in value['data']:
                batch_match = batch_re.search(data[1])
                if batch_match and batch_match.group(1):
                    sample = batch_match.group(1).upper()
                    logging.info(f"CHROMER: SampleID is {sample}.")
                else:
                    logging.warning(f"CHROMER: SampleID is BLANK in run log line. {zip_file}")

                method_match = method_re.search(data[1])
                if method_match:
                    method = method_match.group(1)
                    if "SEC" in method:
                        method = "SEC"
                    elif "Lectin" in method:
                        method = "LEC"
                    elif "Protein A" in method:
                        method = "PROA"
                    elif "IMAC" in method:
                        method = "IMAC"
                    else:
                        method = method.split(' ')[0]
                    logging.info(f"CHROMER: Method is {method_match}. Abbreviating... {method}.")

                date_match = date_re.search(data[1])
                if date_match:
                    date = date_match.group(1)
                    logging.info(f"CHROMER: Run date is {date}.")

    if not sample:
        sample = _resolve_sample_from_filename(zip_file)
        if sample:
            logging.info(f"CHROMER: SampleID resolved from filename: {sample}.")

    return {"method": method, "sample": sample, "date": date}


def lookup_construct(sample, brain, zip_file):
    if not sample or not brain:
        return None, sample

    sample_upper = sample.upper()
    if sample_upper not in brain:
        return None, sample

    construct_name = brain[sample_upper]['ConstructID']
    construct_re = re.compile(f"{construct_name}[-_]?pVAX1?", re.IGNORECASE)
    construct_re2 = re.compile(f"{construct_name}[-_]?mc?", re.IGNORECASE)
    if construct_name is None:
        logging.warning(f"CHROMER: SampleID ({sample}) is NOT a batch#. Checking if it's a constructID... (FAILSAFE-1)")
        construct_name = sample
        for key in brain:
            if construct_re.match(key):
                sample = brain[key]['ConstructID']
                break
            elif construct_re2.match(key):
                sample = brain[key]['ConstructID']
                break
        if sample is None:
            logging.warning(f"CHROMER: SampleID ({construct_name}) is NOT a construct either. Checking if batch# is in the filename... (FAILSAFE-2)")
            filename_match = _FILENAME_RE.search(os.path.basename(zip_file))
            if filename_match:
                sample = filename_match.group(1).upper()
                logging.info(f"CHROMER: Batch# matched from filename. Linking to... {sample}.")
            else:
                logging.warning("CHROMER: Construct unrecognizable.")
                return None, sample
    else:
        construct_name = construct_name.replace("/", "-")
        logging.info(f"CHROMER: Construct recognized as {construct_name}.")

    return construct_name, sample


def format_title(method, sample, construct_name=None):
    if construct_name:
        return f"{method}_{sample} | {construct_name}"
    return f"{method}_{sample}"


def process_chrom(zip_file, brain=None):
    udata = pc_uni7(zip_file)
    udata.load(show=False)
    udata.xml_parse(show=False)
    udata.clean_up()

    meta = mine_run_metadata(udata, zip_file)
    if not meta["method"] or not meta["sample"]:
        logging.warning(f"CHROMER: Could not resolve method/sample for {zip_file}. Skipping...")
        return None

    construct_name, sample = (None, meta["sample"])
    if brain:
        construct_name, sample = lookup_construct(meta["sample"], brain, zip_file)
        if construct_name is None:
            logging.warning(f"CHROMER: SampleID ({sample}) not in brain.json. Using generic title.")

    title = format_title(meta["method"], sample, construct_name)
    logging.info(f"CHROMER: Success. Assigning index: {title}")

    return {
        "udata": udata,
        "method": meta["method"],
        "batch": sample,
        "date": meta["date"],
        "construct_name": construct_name,
        "title": title,
    }


def inspect_summary(result):
    if result is None:
        return None
    return {
        "method": result["method"],
        "batch": result["batch"],
        "date": result["date"],
        "construct_name": result["construct_name"],
        "title": result["title"],
        "curves": sorted(result["udata"].keys()),
    }
