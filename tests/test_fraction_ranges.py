import unittest

from chromer.fractions import get_fraction_ranges, get_pooling_recommendation


class FractionRangeTests(unittest.TestCase):
    def test_returns_ranges_for_each_peak(self):
        x_values = list(range(9))
        y_values = [0, 1, 9, 1, 0, 1, 12, 1, 0]
        frac_data = [(0, "A1"), (3, "A2"), (6, "A3"), (8, "A4")]

        ranges = get_fraction_ranges(x_values, y_values, [2, 6], frac_data)

        self.assertEqual(["A1-A2", "A3-A3"], [label for label, _ in ranges])


class PoolingRecommendationTests(unittest.TestCase):
    def setUp(self):
        self.x_values = list(range(9))
        self.y_values = [0, 1, 9, 1, 0, 1, 12, 1, 0]
        self.frac_data = [(0, "A1"), (3, "A2"), (6, "A3"), (8, "A4")]

    def test_affinity_method_recommends_single_conservative_pool(self):
        recommendation = get_pooling_recommendation(
            "PROA",
            self.x_values,
            self.y_values,
            [2],
            self.frac_data,
        )

        self.assertEqual("pool", recommendation["kind"])
        self.assertEqual("Conservative PROA pool: A1-A2", recommendation["label"])

    def test_affinity_method_defers_multiple_peaks_to_manual_review(self):
        recommendation = get_pooling_recommendation(
            "LEC",
            self.x_values,
            self.y_values,
            [2, 6],
            self.frac_data,
        )

        self.assertEqual("manual_review", recommendation["kind"])
        self.assertEqual("LEC multiple peaks: review manually", recommendation["label"])

    def test_sec_reports_signal_detection_without_pool_recommendation(self):
        recommendation = get_pooling_recommendation(
            "SEC",
            self.x_values,
            self.y_values,
            [2],
            self.frac_data,
        )

        self.assertEqual("signal_only", recommendation["kind"])
        self.assertEqual("SEC peak detected: signal only", recommendation["label"])

    def test_unknown_method_reports_generic_signal_detection(self):
        recommendation = get_pooling_recommendation(
            "CUSTOM",
            self.x_values,
            self.y_values,
            [2],
            self.frac_data,
        )

        self.assertEqual("signal_only", recommendation["kind"])
        self.assertEqual("CUSTOM peak detected: no method-specific pool", recommendation["label"])

    def test_affinity_no_peak_reports_affinity_specific_interpretation(self):
        recommendation = get_pooling_recommendation(
            "IMAC",
            self.x_values,
            self.y_values,
            [],
            self.frac_data,
        )

        self.assertEqual("no_signal", recommendation["kind"])
        self.assertEqual("IMAC no clear affinity peak", recommendation["label"])

    def test_affinity_single_peak_without_fractions_requires_manual_review(self):
        recommendation = get_pooling_recommendation(
            "PROA",
            self.x_values,
            self.y_values,
            [2],
            [],
        )

        self.assertEqual("manual_review", recommendation["kind"])
        self.assertEqual("PROA peak detected: review manually", recommendation["label"])


if __name__ == "__main__":
    unittest.main()
