import io
import logging
import struct
import xml.etree.ElementTree as ET
from collections import OrderedDict
from zipfile import ZipFile, is_zipfile

MARKER_BYTE = 0x0B
LEGACY_DATA_OFFSET = 47
LEGACY_TAIL_BYTES = 48
MAX_REASONABLE_FLOAT_COUNT = 10_000_000  # ponytail: reject garbage int32 matches


def _decode_legacy_offset(blob):
    read_size = len(blob) - LEGACY_TAIL_BYTES
    return [
        struct.unpack("<f", blob[i:i + 4])[0]
        for i in range(LEGACY_DATA_OFFSET, read_size, 4)
    ]


def _try_marker_decode(blob):
    idx = blob.find(bytes([MARKER_BYTE]))
    while idx != -1:
        if idx >= 4:
            count = struct.unpack("<i", blob[idx - 4:idx])[0]
            payload_start = idx + 1
            payload_end = payload_start + count * 4
            if (
                0 < count <= MAX_REASONABLE_FLOAT_COUNT
                and payload_end <= len(blob)
            ):
                return [
                    struct.unpack("<f", blob[i:i + 4])[0]
                    for i in range(payload_start, payload_end, 4)
                ]
        idx = blob.find(bytes([MARKER_BYTE]), idx + 1)
    return None


def decode_float_payload(blob):
    """Return (floats, decoder_name). decoder_name is 'marker' or 'legacy'."""
    marker_values = _try_marker_decode(blob)
    if marker_values is not None:
        logging.info(
            "CHROMER: decode_float_payload marker decoder (count=%d, blob_len=%d)",
            len(marker_values),
            len(blob),
        )
        return marker_values, "marker"

    legacy_values = _decode_legacy_offset(blob)
    logging.info(
        "CHROMER: decode_float_payload legacy decoder (count=%d, blob_len=%d)",
        len(legacy_values),
        len(blob),
    )
    return legacy_values, "legacy"


def validate_chrom_pair_lengths(chrom_key, chrom_dict):
    volumes = chrom_dict.get("CoordinateData.Volumes")
    amplitudes = chrom_dict.get("CoordinateData.Amplitudes")
    if not isinstance(volumes, list) or not isinstance(amplitudes, list):
        return
    if len(volumes) != len(amplitudes):
        logging.warning(
            "CHROMER: length mismatch in %s: volumes=%d amplitudes=%d",
            chrom_key,
            len(volumes),
            len(amplitudes),
        )


class pc_uni7(OrderedDict):
    zip_magic_start = b'\x50\x4B\x03\x04\x2D\x00\x00\x00\x08'
    zip_magic_end = b'\x50\x4B\x05\x06\x00\x00\x00\x00'

    SensData_id = 0
    SensData_id2 = 0
    Fractions_id = 0
    Fractions_id2 = 0

    def __init__(self, udata_file):
        OrderedDict.__init__(self)
        self.file_name = udata_file
        self.inject_vol = 0.0
        self.run_name = 'blank'

    def load(self, show=False):
        with open(self.file_name, 'rb') as f:
            udataut_zip = ZipFile(f)
            zip_data = self.zip2dict(udataut_zip)
            self.update(zip_data)
            proc_yes = []
            proc_no = []
            for i in self.keys():
                tmp_raw = io.BytesIO(udataut_zip.read(i))
                f_header = tmp_raw.read(9)
                if f_header == self.zip_magic_start:
                    proper_zip = tmp_raw.getvalue()
                    f_end = proper_zip.rindex(self.zip_magic_end) + 22
                    tmp_raw = io.BytesIO(proper_zip[0:f_end])
                if is_zipfile(tmp_raw):
                    tmp_zip = ZipFile(tmp_raw)
                    x = {i: self.zip2dict(tmp_zip)}
                    self.update(x)
                    proc_yes.append(i)
                else:
                    proc_no.append(i)
            if show:
                print("Loaded " + self.file_name + " into memory")
                print("\n-Supported-")
                for i in proc_yes:
                    print(" " + i)
                print("\n-Not supported-")
                for i in proc_no:
                    print(" " + i)
        to_process = []
        for i in self.keys():
            if "Chrom" in i and "Xml" not in i:
                to_process.append(i)
        if show:
            print("\nFiles to process:")
            for i in to_process:
                print(" " + i)
        for i in to_process:
            for n in self[i].keys():
                if "DataType" in n:
                    a = self[i][n]
                    b = a.decode('utf-8')
                    x = b.strip("\r\n")
                else:
                    x = self.unpacker(self[i][n])
                tmp_dict = {n: x}
                self[i].update(tmp_dict)
            validate_chrom_pair_lengths(i, self[i])
        if show:
            print("Finished decoding x/y-data!")

    @staticmethod
    def zip2dict(udata):
        mydict = {}
        for i in udata.NameToInfo:
            tmp_dict = {i: udata.read(i)}
            mydict.update(tmp_dict)
        return mydict

    @staticmethod
    def unpacker(udata):
        values, _decoder = decode_float_payload(udata)
        return values

    def xml_parse(self, show=False):
        tree = ET.fromstring(self['Chrom.1.Xml'])
        mc = tree.find('Curves')
        me = tree.find('EventCurves')
        event_dict = {}
        for i in range(len(me)):
            magic_id = self.SensData_id
            e_name = me[i].find('Name').text
            if e_name == 'Fraction':
                e_name = 'Fractions'
            e_orig = me[i].find('IsOriginalData').text
            e_list = me[i].find('Events')
            e_data = []
            for e in range(len(e_list)):
                e_vol = float(e_list[e].find('EventVolume').text)
                e_txt = e_list[e].find('EventText').text
                e_data.append((e_vol, e_txt))
            if e_orig == "false":
                print("not added - not orig data")
            if e_orig == "true":
                x = {'run_name': "Blank", 'data': e_data, 'data_name': e_name, 'magic_id': magic_id}
                event_dict.update({e_name: x})
        self.update(event_dict)
        chrom_dict = {}
        for i in range(len(mc)):
            d_type = mc[i].attrib['CurveDataType']
            d_name = mc[i].find('Name').text
            d_fname = mc[i].find('CurvePoints')[0][1].text
            d_unit = mc[i].find('AmplitudeUnit').text
            magic_id = self.SensData_id
            try:
                x_dat = self[d_fname]['CoordinateData.Volumes']
                y_dat = self[d_fname]['CoordinateData.Amplitudes']
                zdata = list(zip(x_dat, y_dat))
                if d_name == "UV cell path length":
                    d_name = "xUV cell path length"
                x = {'run_name': "Blank", 'data': zdata, 'unit': d_unit, 'data_name': d_name, 'data_type': d_type, 'magic_id': magic_id}
                chrom_dict.update({d_name: x})
            except KeyError:
                pass
            if show:
                print("---")
                print(d_type)
                print(d_name)
                print(d_fname)
                print(d_unit)
        self.update(chrom_dict)

    def clean_up(self):
        manifest = ET.fromstring(self['Manifest.xml'])
        for i in range(len(manifest)):
            file_name = manifest[i][0].text
            self.pop(file_name)
        self.pop('Manifest.xml')
