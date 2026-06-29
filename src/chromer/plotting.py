import logging
import os

import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import mpl_toolkits.axisartist as AA
import seaborn as sns
from mpl_toolkits.axes_grid1 import host_subplot
from scipy.signal import find_peaks

from chromer.fractions import (
    SIGNAL_ONLY_METHODS,
    detect_uv280_signal,
    format_interpretation_label,
)


def annotate_fractions(host, frac_data, injection_time=0):
    for i in range(len(frac_data)):
        adjusted_x = frac_data[i][0] - injection_time
        host.axvline(x=adjusted_x, ymin=0.065, ymax=0.0, color='crimson', linewidth=0.5)
        mid_x = adjusted_x if i == len(frac_data) - 1 else (adjusted_x + frac_data[i + 1][0] - injection_time) / 2
        host.annotate(
            str(frac_data[i][1]),
            xy=(mid_x, 0),
            xytext=(0, -5),
            textcoords='offset points',
            ha='center',
            fontsize=12,
            rotation=90,
            fontweight='bold',
        )


def annotate_peaks(host, x_values, y_values, peaks):
    for i, peak in enumerate(peaks):
        offset = (i % 2) * 10
        host.annotate(
            f"V: {x_values[peak]:.2f}\nA: {y_values[peak]:.2f}",
            (x_values[peak], y_values[peak]),
            textcoords="offset points",
            xytext=(0, 5 + offset),
            ha='center',
            fontweight='bold',
            fontsize=10,
            color='red',
            bbox=dict(facecolor='none', edgecolor='red', boxstyle='round,pad=0.2'),
        )


def plot_data(host, x_values, y_values, color='blue', linewidth=5, label='UV 280nm'):
    host.plot(x_values, y_values, color=color, linewidth=linewidth, label=label)
    return host


def save_plot(title, method_folder):
    filename = os.path.join(method_folder, f"{title}.jpg")
    plt.savefig(filename, bbox_inches='tight', dpi=300)
    logging.info(f"CHROMER: Saved chromatogram to {filename}")


def annotate_summary(host, label):
    fontsize = 16 if len(label) > 60 else 18
    host.text(
        1,
        1.05,
        label,
        ha="right",
        va="top",
        transform=host.transAxes,
        fontsize=fontsize,
        fontweight="bold",
        color="black",
        wrap=True,
    )


def apply_plotting_params(params):
    pylab.rcParams.update(params)


def chromeunicorns(zip_file, udata_info, method_folder):
    udata = udata_info['udata']
    method = udata_info['method']
    batch = udata_info['batch']
    title = udata_info['title']
    date = udata_info['date']

    if not udata or 'Fractions' not in udata or None in [batch, method, title] or "" in [batch, method, title]:
        print(
            f"CHROMER: Data for {zip_file} is invalid. "
            f"Invalid: FRACTIONS - {not udata or udata.get('Fractions') is None}, "
            f"BATCH - {batch is None}, METHOD - {method is None}. Skipping..."
        )
        return False

    plt.figure(figsize=(15, 10), edgecolor='black')
    sns.set_style("whitegrid")
    host = host_subplot(111, axes_class=AA.Axes)

    injection_time = 0
    if method == 'SEC':
        injection_time = udata.get('Injection', {}).get('data', [[0, None]])[0][0]

    x_values = []
    y_values = []
    for key, value in udata.items():
        if key in ["UV 1_280", "UV"]:
            x_values = [x[0] - injection_time for x in value['data']]
            y_values = [y[1] for y in value['data']]
            plot_data(host, x_values, y_values)

    if 'Fractions' in udata and method not in SIGNAL_ONLY_METHODS:
        xlim_min = udata['Fractions']['data'][0][0]
        xlim_max = udata['Fractions']['data'][-1][0]
    elif method in SIGNAL_ONLY_METHODS:
        xlim_min = 0
        xlim_max = 30

    host.set_xlim(xlim_min, xlim_max)
    y_values_in_xlim = [y for x, y in zip(x_values, y_values) if xlim_min <= x <= xlim_max]
    if y_values_in_xlim:
        ymax = max(y_values_in_xlim)
        host.set_ylim(-0.5, ymax * 1.05)
    else:
        logging.warning("No y_values within specified x limits. Skipping...")

    peaks, _properties = find_peaks(y_values, prominence=5, width=1, height=10)
    peaks = [peak for peak in peaks if xlim_min <= x_values[peak] <= xlim_max]
    annotate_peaks(host, x_values, y_values, peaks)

    frac_data = udata['Fractions']['data'] if 'Fractions' in udata else None
    if frac_data:
        annotate_fractions(host, frac_data, injection_time)

    signal = detect_uv280_signal(x_values, y_values, peaks, frac_data, method)
    label = format_interpretation_label(method, signal, udata_info.get("construct_meta"))
    annotate_summary(host, label)

    host.set_title(title, loc='left', fontsize=24, fontweight='bold', color='black', pad=20)
    host.legend(fancybox=True, framealpha=1, shadow=True, borderpad=1, fontsize=22, loc='upper right', edgecolor='black')
    host.set_xlabel('Elution volume (ml)')
    host.set_ylabel('Absorbance (mAU)')

    host.axis["bottom"].major_ticklabels.set_pad(20)
    host.axis["bottom"].label.set_weight('bold')
    host.axis["left"].label.set_weight('bold')

    if date is not None:
        plt.text(0, 1.015, f"{date}", ha='left', va='top', transform=plt.gca().transAxes, style='italic', fontsize=10)

    plt.minorticks_on()
    plt.subplots_adjust(bottom=0.2)
    plt.tight_layout(pad=2)
    save_plot(title, method_folder)
    plt.close()
    return True
