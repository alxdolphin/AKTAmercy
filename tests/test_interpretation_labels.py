import unittest

from chromer.fractions import SignalResult, detect_uv280_signal, format_interpretation_label


def _signal(has_peak, peak_count=0, frac_ranges=None, is_multi=False):
    return SignalResult(
        has_peak=has_peak,
        peak_count=peak_count,
        frac_ranges=frac_ranges or [],
        is_multi=is_multi,
    )


class InterpretationLabelTests(unittest.TestCase):
    def test_no_peak_any_method(self):
        signal = _signal(False)
        for method in ("PROA", "LEC", "IMAC", "SEC", "Nickel"):
            self.assertEqual(
                format_interpretation_label(method, signal),
                "No clear UV280 elution peak",
            )

    def test_proa_single_peak(self):
        signal = _signal(True, 1, ["A1-A2"])
        self.assertEqual(
            format_interpretation_label("PROA", signal),
            "Candidate Protein A elution fraction: A1-A2",
        )

    def test_proa_multi_peak(self):
        signal = _signal(True, 2, ["A1-A2", "A3-A4"], is_multi=True)
        self.assertEqual(
            format_interpretation_label("PROA", signal),
            "MULTI — Candidate Protein A elution fractions: A1-A2, A3-A4",
        )

    def test_lec_single_peak(self):
        signal = _signal(True, 1, ["B1-B2"])
        self.assertEqual(
            format_interpretation_label("LEC", signal),
            "Candidate lectin-binding fraction: B1-B2",
        )

    def test_lec_multi_peak(self):
        signal = _signal(True, 2, ["B1-B2", "B3-B3"], is_multi=True)
        self.assertEqual(
            format_interpretation_label("LEC", signal),
            "MULTI — Candidate lectin-binding fractions: B1-B2, B3-B3",
        )

    def test_imac_peak(self):
        signal = _signal(True, 1)
        self.assertEqual(
            format_interpretation_label("IMAC", signal),
            "IMAC UV280 peak detected: verify target-containing fractions",
        )

    def test_sec_peak(self):
        signal = _signal(True, 1)
        self.assertEqual(
            format_interpretation_label("SEC", signal),
            "SEC UV280 peak detected: signal only; verify identity/oligomeric state separately",
        )

    def test_other_method_peak(self):
        signal = _signal(True, 1)
        self.assertEqual(
            format_interpretation_label("Nickel", signal),
            "UV280 peak detected: verify target-containing fractions",
        )

    def test_proa_target_class_igg_upgrade(self):
        signal = _signal(True, 1, ["A1-A2"])
        meta = {"ConstructID": "SomeConstruct", "TargetClass": "IgG"}
        self.assertEqual(
            format_interpretation_label("PROA", signal, meta),
            "Candidate IgG/Fc capture fraction: A1-A2",
        )

    def test_proa_construct_heuristic_upgrade(self):
        signal = _signal(True, 1, ["A1-A2"])
        meta = {"ConstructID": "PGT128- IgG1"}
        self.assertEqual(
            format_interpretation_label("PROA", signal, meta),
            "Candidate IgG/Fc capture fraction: A1-A2",
        )

    def test_lec_gp140_name_stays_cautious(self):
        signal = _signal(True, 1, ["A1-A2"])
        meta = {"ConstructID": "BG505_gp140_foldon_PDGFR_pVAX"}
        self.assertEqual(
            format_interpretation_label("LEC", signal, meta),
            "Candidate lectin-binding fraction: A1-A2",
        )

    def test_lec_glycoprotein_target_class_upgrade(self):
        signal = _signal(True, 2, ["A1-A2", "A3-A4"], is_multi=True)
        meta = {"ConstructID": "BG505_gp140", "TargetClass": "glycoprotein"}
        self.assertEqual(
            format_interpretation_label("LEC", signal, meta),
            "MULTI — Candidate glycoprotein pool: A1-A2, A3-A4",
        )


class DetectUv280SignalTests(unittest.TestCase):
    def test_pool_candidate_multi_peak_collects_fractions(self):
        x_values = list(range(9))
        y_values = [0, 1, 9, 1, 0, 1, 12, 1, 0]
        frac_data = [(0, "A1"), (3, "A2"), (6, "A3"), (8, "A4")]

        signal = detect_uv280_signal(x_values, y_values, [2, 6], frac_data, "PROA")

        self.assertTrue(signal.has_peak)
        self.assertEqual(2, signal.peak_count)
        self.assertTrue(signal.is_multi)
        self.assertEqual(["A1-A2", "A3-A3"], signal.frac_ranges)

    def test_imac_does_not_compute_fraction_ranges(self):
        x_values = list(range(9))
        y_values = [0, 1, 9, 1, 0, 1, 12, 1, 0]
        frac_data = [(0, "A1"), (3, "A2"), (6, "A3"), (8, "A4")]

        signal = detect_uv280_signal(x_values, y_values, [2], frac_data, "IMAC")

        self.assertTrue(signal.has_peak)
        self.assertEqual([], signal.frac_ranges)
        self.assertFalse(signal.is_multi)


if __name__ == "__main__":
    unittest.main()
