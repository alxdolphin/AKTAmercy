import datetime
import logging
import os
import re


def setup_logger():
    logging.getLogger('matplotlib').setLevel(logging.WARNING)
    logging.getLogger('pc_uni7').setLevel(logging.WARNING)
    log_filename = datetime.datetime.now().strftime("%Y%m%d_%H%M%S") + '_chromatics.log'
    log_dir = './dev/debug/logs'
    os.makedirs(log_dir, exist_ok=True)
    log_filepath = os.path.join(log_dir, log_filename)
    logging.basicConfig(
        filename=log_filepath,
        level=logging.DEBUG,
        format='%(asctime)s [%(levelname)s] %(message)s\n\n',
        datefmt='%Y-%m-%d %H:%M:%S',
    )
    logging.info(
        "CHROMER: Initialized at %s\n\n",
        datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
    )


def summarize_warnings(log_file_path):
    pattern = re.compile(r'.*\[WARNING\].*')
    warning_logs = []
    unprocessed_chromatograms = []
    count = 0
    output_file_path = os.path.splitext(log_file_path)[0] + '_warnings.log'

    with open(log_file_path, 'r') as log_file:
        with open(output_file_path, 'w') as output_file:
            for line in log_file:
                if pattern.match(line):
                    warning_logs.append(line)
                    count += 1
                    match = re.search(r'Skipping... (.+)$', line)
                    if match:
                        unprocessed_chromatograms.append(match.group(1))

                if len(warning_logs) >= 1000:
                    output_file.writelines(warning_logs)
                    warning_logs.clear()

            if warning_logs:
                output_file.writelines(warning_logs)

            output_file.write("\n\n--- Summary ---\n")
            output_file.write("Total WARN messages: {}\n".format(count))
            output_file.write("Unprocessed Chromatograms:\n")
            for identifier in unprocessed_chromatograms:
                output_file.write(f"{identifier}\n")

    return f"Extracted {count} warning logs and saved to {output_file_path}"
