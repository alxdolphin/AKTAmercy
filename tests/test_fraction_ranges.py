import unittest

from CHROMER import get_fraction_ranges


class FractionRangeTests(unittest.TestCase):
    def test_returns_ranges_for_each_peak(self):
        x_values = list(range(9))
        y_values = [0, 1, 9, 1, 0, 1, 12, 1, 0]
        frac_data = [(0, "A1"), (3, "A2"), (6, "A3"), (8, "A4")]

        ranges = get_fraction_ranges(x_values, y_values, [2, 6], frac_data)

        self.assertEqual(["A1-A2", "A3-A3"], [label for label, _ in ranges])


if __name__ == "__main__":
    unittest.main()
