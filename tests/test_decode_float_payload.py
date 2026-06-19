import struct
import unittest
from unittest.mock import patch

from chromer.unicorn import (
    MARKER_BYTE,
    decode_float_payload,
    validate_chrom_pair_lengths,
)


def _legacy_blob(floats):
    payload = b"".join(struct.pack("<f", value) for value in floats)
    return b"\x00" * 47 + payload + b"\x00" * 48


def _marker_blob(floats, prefix=b"\x00" * 8):
    count = len(floats)
    frame = (
        struct.pack("<i", count)
        + bytes([MARKER_BYTE])
        + b"".join(struct.pack("<f", value) for value in floats)
    )
    return prefix + frame


class DecodeFloatPayloadTests(unittest.TestCase):
    def test_legacy_path_roundtrip(self):
        expected = [0.0, 1.5, 3.25, 42.0]
        values, decoder = decode_float_payload(_legacy_blob(expected))

        self.assertEqual("legacy", decoder)
        self.assertEqual(expected, values)

    def test_marker_path_roundtrip(self):
        expected = [0.0, 2.5, 5.0, 10.0]
        values, decoder = decode_float_payload(_marker_blob(expected))

        self.assertEqual("marker", decoder)
        self.assertEqual(expected, values)

    def test_marker_preferred_over_legacy(self):
        expected = [1.0, 2.0, 3.0]
        marker_frame = (
            struct.pack("<i", len(expected))
            + bytes([MARKER_BYTE])
            + b"".join(struct.pack("<f", value) for value in expected)
        )
        blob = _legacy_blob([99.0, 88.0])[:47] + marker_frame + b"\x00" * 48

        values, decoder = decode_float_payload(blob)

        self.assertEqual("marker", decoder)
        self.assertEqual(expected, values)

    def test_invalid_marker_falls_back(self):
        expected = [0.5, 1.5, 2.5]
        legacy = _legacy_blob(expected)
        header = bytearray(legacy[:47])
        header[0:5] = struct.pack("<i", 999999) + bytes([MARKER_BYTE])
        blob = bytes(header) + legacy[47:]

        values, decoder = decode_float_payload(blob)

        self.assertEqual("legacy", decoder)
        self.assertEqual(expected, values)

    @patch("chromer.unicorn.logging.warning")
    def test_length_validation_helper(self, mock_warning):
        validate_chrom_pair_lengths(
            "Chrom.1_2_True",
            {
                "CoordinateData.Volumes": [0.0, 1.0, 2.0],
                "CoordinateData.Amplitudes": [0.0, 1.0],
            },
        )

        mock_warning.assert_called_once_with(
            "CHROMER: length mismatch in %s: volumes=%d amplitudes=%d",
            "Chrom.1_2_True",
            3,
            2,
        )


if __name__ == "__main__":
    unittest.main()
