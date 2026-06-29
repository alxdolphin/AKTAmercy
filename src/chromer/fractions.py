import re
from dataclasses import dataclass

from scipy.signal import peak_widths


POOL_CANDIDATE_METHODS = frozenset({"LEC", "PROA"})
VERIFY_METHODS = frozenset({"IMAC"})
SIGNAL_ONLY_METHODS = frozenset({"SEC"})
KNOWN_METHODS = POOL_CANDIDATE_METHODS | VERIFY_METHODS | SIGNAL_ONLY_METHODS

IGG_FC_PATTERN = re.compile(r"IgG1| IgG|_IgG|CHIg|CLIg|_Fc|huFc", re.IGNORECASE)


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


@dataclass
class SignalResult:
    has_peak: bool
    peak_count: int
    frac_ranges: list
    is_multi: bool


def detect_uv280_signal(x_values, y_values, peaks, frac_data, method):
    peak_count = len(peaks)
    has_peak = peak_count > 0
    frac_ranges = []
    if has_peak and method in POOL_CANDIDATE_METHODS and frac_data:
        frac_ranges = [
            label for label, _ in get_fraction_ranges(x_values, y_values, peaks, frac_data)
        ]
    is_multi = has_peak and peak_count > 1 and method in POOL_CANDIDATE_METHODS
    return SignalResult(
        has_peak=has_peak,
        peak_count=peak_count,
        frac_ranges=frac_ranges,
        is_multi=is_multi,
    )


def _is_igg_fc_target(construct_meta):
    target_class = construct_meta.get("TargetClass")
    if target_class in ("IgG", "Fc"):
        return True
    construct_id = construct_meta.get("ConstructID", "")
    return bool(IGG_FC_PATTERN.search(construct_id))


def _format_pool_candidate_label(
    method,
    signal,
    construct_meta,
    default_single,
    default_multi,
    upgraded_single,
    upgraded_multi,
):
    ranges_text = ", ".join(signal.frac_ranges)
    use_upgrade = False
    if method == "PROA":
        use_upgrade = _is_igg_fc_target(construct_meta) and signal.frac_ranges
    elif method == "LEC":
        use_upgrade = construct_meta.get("TargetClass") == "glycoprotein" and signal.frac_ranges

    if use_upgrade:
        if signal.is_multi:
            return f"MULTI — {upgraded_multi}: {ranges_text}"
        return f"{upgraded_single}: {signal.frac_ranges[0]}"

    if signal.is_multi:
        return f"MULTI — {default_multi}: {ranges_text}"
    if signal.frac_ranges:
        return f"{default_single}: {signal.frac_ranges[0]}"
    return "UV280 peak detected: verify target-containing fractions"


def format_interpretation_label(method, signal, construct_meta=None):
    construct_meta = construct_meta or {}
    if not signal.has_peak:
        return "No clear UV280 elution peak"

    if method == "PROA":
        return _format_pool_candidate_label(
            method,
            signal,
            construct_meta,
            "Candidate Protein A elution fraction",
            "Candidate Protein A elution fractions",
            "Candidate IgG/Fc capture fraction",
            "Candidate IgG/Fc capture fractions",
        )
    if method == "LEC":
        return _format_pool_candidate_label(
            method,
            signal,
            construct_meta,
            "Candidate lectin-binding fraction",
            "Candidate lectin-binding fractions",
            "Candidate glycoprotein pool",
            "Candidate glycoprotein pool",
        )
    if method in VERIFY_METHODS:
        return "IMAC UV280 peak detected: verify target-containing fractions"
    if method in SIGNAL_ONLY_METHODS:
        return "SEC UV280 peak detected: signal only; verify identity/oligomeric state separately"
    return "UV280 peak detected: verify target-containing fractions"
