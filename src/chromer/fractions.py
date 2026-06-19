from scipy.signal import peak_widths


def get_fraction_ranges(x_values, y_values, peaks, frac_data):
    frac_ranges = []
    for peak in peaks:
        widths = peak_widths(y_values, [peak])[0][0] * 1.5
        start = int(max(0, peak - widths / 2))
        end = int(min(len(y_values) - 1, peak + widths / 2))
        frac_start = min(frac_data, key=lambda x: abs(x[0] - x_values[start]))[1]
        frac_end = min(frac_data, key=lambda x: abs(x[0] - x_values[end]))[1]
        if frac_start and frac_end:
            frac_ranges.append((f"{frac_start}-{frac_end}", end))
    return sorted(frac_ranges, key=lambda x: x[1])


def annotate_fraction_ranges(host, frac_ranges):
    for _frac_range, _end_frac_index in frac_ranges:
        host.text(1, 1.05, _frac_range, ha='right', va='top', transform=host.transAxes, fontsize=24, fontweight='bold', color='black')
