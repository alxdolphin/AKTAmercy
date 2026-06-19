from scipy.signal import peak_widths


AFFINITY_METHODS = {"LEC", "PROA", "IMAC"}


def get_fraction_ranges(x_values, y_values, peaks, frac_data):
    if not frac_data:
        return []

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


def get_pooling_recommendation(method, x_values, y_values, peaks, frac_data):
    method_label = method or "Unknown"
    if not peaks:
        peak_type = "affinity peak" if method in AFFINITY_METHODS else "signal peak"
        return {
            "kind": "no_signal",
            "label": f"{method_label} no clear {peak_type}",
            "ranges": [],
        }

    if method == "SEC":
        return {
            "kind": "signal_only",
            "label": "SEC peak detected: signal only",
            "ranges": [],
        }

    if method not in AFFINITY_METHODS:
        return {
            "kind": "signal_only",
            "label": f"{method_label} peak detected: no method-specific pool",
            "ranges": [],
        }

    if len(peaks) > 1:
        return {
            "kind": "manual_review",
            "label": f"{method_label} multiple peaks: review manually",
            "ranges": [],
        }

    frac_ranges = get_fraction_ranges(x_values, y_values, peaks, frac_data)
    if not frac_ranges:
        return {
            "kind": "manual_review",
            "label": f"{method_label} peak detected: review manually",
            "ranges": [],
        }

    frac_range = frac_ranges[0][0]
    return {
        "kind": "pool",
        "label": f"Conservative {method_label} pool: {frac_range}",
        "ranges": frac_ranges,
    }


def annotate_fraction_ranges(host, frac_ranges):
    for _frac_range, _end_frac_index in frac_ranges:
        host.text(1, 1.05, _frac_range, ha='right', va='top', transform=host.transAxes, fontsize=24, fontweight='bold', color='black')
