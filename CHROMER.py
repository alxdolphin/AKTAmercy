#!/usr/bin/env python3

'''
PROBLEMS
- AFFINITY Chromatograms are numerically offset by an unknown factor (sample flow?)
- Log files are filled with duplicate statements, and require a complementary script to clean up
- BIG: Google's API requires re-authentication every 1 hour, which is a problem for long runs
- Cleaning up of plot annotations (overlap, etc.)
- Matrix system is hacky and needs to be reworked, preferably directly into the queue
- Output organization is messy and needs to be reworked
'''

## IMPORTS - GENERAL ##
import datetime # For timestamping log files
import io # For reading zip files, and other file IO
import json # For user configuration, construct-chromatogram pairing, and chromatogram parameter X/Y dictionaries
import logging
import os # For file IO
import re # For parsing dictionary keys for chomatogram metadata (batch, method, etc.)
import struct 
import tarfile # For extracting UFol files
# import tempfile # For creating temporary directories | Not used - may be useful in the future, but max coverage is goal for now

## IMPORTS - PyCORN-SPECIFIC ##
import xml.etree.ElementTree as ET 
from collections import OrderedDict
from zipfile import ZipFile, is_zipfile
# import numpy as np # For peak detection | Not used - Scipy is more efficient

## IMPORTS - GOOGLE API ##
import gspread # For accessing the Queue and Matrix spreadsheets
from pydrive.auth import GoogleAuth # For upload of chromatograms to Google Drive
from pydrive.drive import GoogleDrive 
from oauth2client.service_account import ServiceAccountCredentials
import time # For exponential back-off when uploading to Google Drive | TODO: See if this is still necessary

## IMPORTS - PLOTTING ##
import matplotlib.pylab as pylab # For plotting chromatograms
import matplotlib.pyplot as plt 
import mpl_toolkits.axisartist as AA 
import seaborn as sns # For some additional plotting functionality
from mpl_toolkits.axes_grid1 import host_subplot
from scipy.signal import find_peaks, peak_widths # For peak detection

## PyCORN - UNI7 HACK ##
class pc_uni7(OrderedDict):
    '''
    A class for holding the pycorn/RESv6 data
    A subclass of `dict`, with the form `data_name`: `data`.
    '''
    # for manual zip-detection
    zip_magic_start = b'\x50\x4B\x03\x04\x2D\x00\x00\x00\x08'
    zip_magic_end = b'\x50\x4B\x05\x06\x00\x00\x00\x00'
    
    # hack to get pycorn-bin to move on
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
        '''
        zip-files inside the zip-bundle are replaced by dicts, again with dicts with filename:content
        Chrom.#_#_True (=zip-files) files are unpacked from binary to floats by unpacker()
        To access x/y-value of Chrom.1_2:
        udata = pc_uni7("mybundle.zip")
        udata.load()
        x = udata['Chrom.1_2_True']['CoordinateData.Volumes']
        y = udata['Chrom.1_2_True']['CoordinateData.Amplitudes']
        '''
        with open(self.file_name, 'rb') as f:
            udataut_zip = ZipFile(f)
            zip_data = self.zip2dict(udataut_zip)
            self.update(zip_data)
            proc_yes = []
            proc_no = []
            for i in self.keys():
                tmp_raw = io.BytesIO(udataut_zip.read(i))
                f_header = tmp_raw.read(9)
                # tmp_raw.seek(0)
                # the following if block is to fix the non-standard zip files
                # by stripping out all the null-bytes at the end
                # see https://bugs.python.org/issue24621
                if f_header == self.zip_magic_start:
                    proper_zip = tmp_raw.getvalue()
                    f_end = proper_zip.rindex(self.zip_magic_end) + 22
                    tmp_raw = io.BytesIO(proper_zip[0:f_end])
                if is_zipfile(tmp_raw):
                    tmp_zip = ZipFile(tmp_raw)
                    x = {i:self.zip2dict(tmp_zip)}
                    self.update(x)
                    proc_yes.append(i)
                else:
                    pass
                    proc_no.append(i)
            if show:
                print("Loaded " + self.file_name + " into memory")
                print("\n-Supported-")
                for i in proc_yes:
                    print(" " + i)
                print("\n-Not supported-")
                for i in proc_no:
                    print(" " + i)
        # filter out data we dont deal with atm
        to_process = []
        for i in self.keys():
            if "Chrom" in i and not "Xml" in i:
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
                tmp_dict = {n:x}
                self[i].update(tmp_dict)
        if show:
            print("Finished decoding x/y-data!")

    @staticmethod
    def zip2dict(udata):
        '''
        udataut = zip object
        outout = dict with filename:file-object pairs
        '''
        mydict = {}
        for i in udata.NameToInfo:
            tmp_dict = {i:udata.read(i)}
            mydict.update(tmp_dict)
        return(mydict)
    
    @staticmethod
    def unpacker(udata):
        '''
        udataut = data block
        output = list of values
        '''
        read_size = len(udata) - 48
        values = []
        for i in range(47, read_size, 4):
            x = struct.unpack("<f", udata[i:i+4])
            x = x[0]
            values.append(x)
        return(values)
   
    def xml_parse(self,show=False):
        '''
        parses parts of the Chrom.1.Xml and creates a res3-like dict
        '''
        tree = ET.fromstring(self['Chrom.1.Xml'])
        mc = tree.find('Curves')
        me = tree.find('EventCurves')
        #print(tree.tag) # NOTE: I commented these out because they were printing to the console, and I don't think they're necessary (adalton)
        #print(tree.attrib) # NOTE: See above (adalton)
        event_dict = {}
        for i in range(len(me)):
            magic_id = self.SensData_id
            e_type = me[i].attrib['EventCurveType']
            e_name = me[i].find('Name').text
            if e_name == 'Fraction':
                e_name = 'Fractions' # another hack for pycorn-bin
            e_orig = me[i].find('IsOriginalData').text
            e_list = me[i].find('Events')
            e_data = []
            for e in range(len(e_list)):
                e_vol = float(e_list[e].find('EventVolume').text)
                e_txt = e_list[e].find('EventText').text
                e_data.append((e_vol,e_txt))
            if e_orig == "false":
                print("not added - not orig data")
            if e_orig == "true":
                #print("added - orig data") # NOTE: See above (adalton)
                x = {'run_name':"Blank", 'data': e_data, 'data_name':e_name, 'magic_id':magic_id}
                event_dict.update({e_name:x})
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
                zdata = list(zip(x_dat,y_dat))
                if d_name == "UV cell path length":
                    d_name = "xUV cell path length" # hack to prevent pycorn-bin from picking this up
                x = {'run_name':"Blank", 'data': zdata, 'unit': d_unit, 'data_name':d_name, 'data_type':d_type, 'magic_id':magic_id}
                chrom_dict.update({d_name:x})
            except:
                KeyError
                # don't deal with data that does not make sense atm
                # orig2.zip contains UV-blocks that are (edited) copies of
                # original UV-trace but they dont have the volume data
            if show:
                print("---")
                print(d_type)
                print(d_name)
                print(d_fname)
                print(d_unit)
        self.update(chrom_dict)
    
    def clean_up(self):
        '''
        deletes everything and just keeps relevant run-date
        resulting dict is more like res3
        '''
        manifest = ET.fromstring(self['Manifest.xml'])
        for i in range(len(manifest)):
            file_name = manifest[i][0].text
            self.pop(file_name)
        self.pop('Manifest.xml')

## GLOBALS / CONFIG ##
with open('config.json') as f:
    config = json.load(f)
    
unicorns_path = config['paths']['unicorns']
run_folder_base_path = config['paths']['run_folder_base']
credentials_path = config['paths']['credentials']

scopes = ["https://spreadsheets.google.com/feeds", 
          "https://www.googleapis.com/auth/spreadsheets",
          "https://www.googleapis.com/auth/drive.file", 
          "https://www.googleapis.com/auth/drive"]

creds = ServiceAccountCredentials.from_json_keyfile_name(credentials_path, scopes)

gauth = GoogleAuth()
gauth.credentials = creds
drive = GoogleDrive(gauth)
gauth = GoogleAuth()
gauth.credentials = creds
drive = GoogleDrive(gauth)
client = gspread.authorize(creds)
queue_spreadsheet = config['spreadsheets']['queue']
queue = client.open(queue_spreadsheet)
matrix_spreadsheet = config['spreadsheets']['matrix'] # TODO: Write version without matrix
matrix = client.open(matrix_spreadsheet)
brain_path = config['paths']['brain']
mode = config['mode']
params = config['plotting']['params']

## LOGGER ##
def setup_logger():
    logging.getLogger('matplotlib').setLevel(logging.WARNING)
    logging.getLogger('oauth2client').setLevel(logging.WARNING)
    logging.getLogger('googleapiclient').setLevel(logging.WARNING)
    logging.getLogger('pc_uni7').setLevel(logging.WARNING)
    log_filename = datetime.datetime.now().strftime("%Y%m%d_%H%M%S") + '_chromatics.log'
    log_filepath = os.path.join('./dev/debug/logs', log_filename)
    logging.basicConfig(filename=log_filepath, level=logging.DEBUG,
                        format='%(asctime)s [%(levelname)s] %(message)s\n\n',
                        datefmt='%Y-%m-%d %H:%M:%S')
    logging.info(f"CHROMER: Initialized at {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")    

## CHROMATOGRAM-CONSTRUCT PAIRING ##
def construct_recognition(queue, output_path):
    worksheet_queue = queue.get_worksheet(0)
    worksheet_antibodies = queue.get_worksheet(1)
    values_queue = worksheet_queue.get_all_values()
    values_antibodies = worksheet_antibodies.get_all_values()

    filtered_queue = {
        row[3].upper(): {"RowIndex": f"Queue!A{index+2}", "ConstructID": row[4], "Ext/1000": row[13], "MW/1000": row[15]} # TODO: Add SEC Column and Pooled fractionS
        for index, row in enumerate(values_queue[1:]) if len(row) >= 16
    }
    filtered_antibodies = {
        row[3].upper(): {"RowIndex": f"Antibodies!A{index+2}", "ConstructID": row[4], "Ext/1000": row[14], "MW/1000": row[16]} 
        for index, row in enumerate(values_antibodies[1:]) if len(row) >= 16
    }
    consolidated_dict = {**filtered_queue, **filtered_antibodies}
    consolidated_json = json.dumps(consolidated_dict, indent=4, sort_keys=True, ensure_ascii=False)
    
    try:
        with open(output_path, 'w') as f:
            f.write(consolidated_json)
    except Exception as e:
        logging.error(f"CHROMER: Failed to fetch updated Queue. Falling back to recent memory.\n    REASON: {e}")
        return None
def enable_cognition(file_path): # Kept as a seperate function in case internet access is not available
    try:
        with open(file_path) as f:
            brain = json.load(f)
        logging.info("CHROMER: Brain loaded successfully. Chromatogram recognition enabled.")
        return brain
    except Exception as e:
        logging.error(f"CHROMER: Brain empty. Unable to recognize constructs.\n    REASON: {e}")
        return None 


## CHROMATOGRAM META-MINING ##
def process_chrom(zip_file, brain):
    udata = pc_uni7(zip_file)
    udata.load(show=False)
    udata.xml_parse(show=False)
    udata.clean_up()

    method = batch = date = None
    batch_re = re.compile(r'Set mark "([^-\s]*)"', re.IGNORECASE)
    method_re = re.compile(r'Method:(.*)', re.IGNORECASE)
    date_re = re.compile(r'Method Run ((\d{1,2}/\d{1,2}/\d{4}))')
    filename_re = re.compile(r'[bB]#([^-\s]*)', re.IGNORECASE)

    for key, value in udata.items():
        if key == "Run Log":
            for data in value['data']:
               
                # Determine batch
                batch_match = batch_re.search(data[1])
                if batch_match and batch_match.group(1):
                    batch = batch_match.group(1).upper()
                    logging.info(f"CHROMER: SampleID is {batch}.")
                else:
                    logging.warning(f"CHROMER: SampleID is BLANK. Skipping... {zip_file}.")
               
                # Determine method
                method_match = method_re.search(data[1])
                if method_match:
                    method = method_match.group(1)
                    if "SEC" in method:
                        method = "SEC"
                    elif "Lectin" in method:
                        method = "LEC"
                    elif "Protein A" in method:
                        method = "PROA"
                    elif "IMAC" in method:
                        method = "IMAC"
                    else:
                        method = method.split(' ')[0]
                    logging.info(f"CHROMER: Method is {method_match}. Abbreviating... {method}.")

                # Determine date
                date_match = date_re.search(data[1])
                if date_match:
                    date = date_match.group(1)
                    logging.info(f"CHROMER: Run date is {date}.")
                    
            construct_name = None
            if batch:
                batch_upper = batch.upper()
                if batch_upper in brain:
                    construct_name = brain[batch_upper]['ConstructID']
                    construct_re = re.compile(f"{construct_name}[-_]?pVAX1?", re.IGNORECASE)
                    construct_re2 = re.compile(f"{construct_name}[-_]?mc?", re.IGNORECASE)
                    if construct_name is None:
                        logging.warning(f"CHROMER: SampleID ({batch}) is NOT a batch#. Checking if it's a constructID... (FAILSAFE-1) ")
                        construct_name = batch
                        for key in brain:
                            if construct_re.match(key):
                                batch = brain[key]['ConstructID']
                                break
                            elif construct_re2.match(key):
                                batch = brain[key]['ConstructID']
                                break
                        if batch is None:
                            logging.warning(f"CHROMER: SampleID ({construct_name}) is NOT a construct either. Checking if batch# is in the filename... (FAILSAFE-2)")
                            filename = os.path.basename(zip_file)
                            filename_match = filename_re.search(filename)
                            if filename_match:
                                batch = filename_match.group(1).upper()
                                logging.info(f"CHROMER: Batch# matched from filename. Linking to... {batch}.")
                            else:
                                logging.warning(f"CHROMER: Construct unrecognizable. Skipping...")
                                return None, None, None
                    else:
                        construct_name = construct_name.replace("/", "-")
                        logging.info(f"CHROMER: Construct recognized as {construct_name}.")
                
    Title = (f"{method}_{batch} | {construct_name}")
    logging.info(f"CHROMER: Success. Assigning index: {Title}")

    if 'udata' not in locals():
        udata = None

    return {"udata": udata, "method": method, "batch": batch, "date": date, "construct_name": construct_name, "title": Title}

def upload_file(drive, filepath, method, batch, brain, max_retries=5):
    logging.info(f"CHROMER: Sending to the CHROMATRIX as {filepath}")

    for n in range(max_retries):
        try:
            file = drive.CreateFile({'title': os.path.basename(filepath), 'parents': [{'id':'1-f-RQYJ_w_YAlcddoEfFmZzxAsyvLZEQ'}]})
            file.SetContentFile(filepath)
            file.Upload()

            file_id = file['id']
            file_link = f"https://drive.google.com/file/d/{file_id}/view?usp=sharing"

            if batch and batch.upper() in (key.upper() for key in brain):
                try:    
                    if 'Antibodies' in brain[batch.upper()]['RowIndex']:
                        logging.info(f"CHROMER: {batch} is an antibody. Updating Antibodies sheet... {filepath}")
                        worksheet = matrix.get_worksheet(1)
                    else:
                        logging.info(f"CHROMER: {batch} is a non-antibody protein. Updating Queue sheet... {filepath}")
                        worksheet = matrix.get_worksheet(0)
                
                    column = 'E' if method != 'SEC' else 'F'
                    cell = brain[batch.upper()]['RowIndex'].replace('A', column)
                    worksheet.update_acell(cell.split('!')[1], file_link)
                    
                except Exception as e:
                    logging.error(f"CHROMER: Unable to access CHROMATRIX.\n    REASON: {e}")
                    raise e
            break
        except Exception as e:
            if n < max_retries - 1: # if it's not the last attempt
                sleep_time = 2 ** n  # Exponential back-off
                time.sleep(sleep_time)
                logging.warning(f"CHROMER: Exhausted API allowance. Sleeping {sleep_time} seconds before recalling...")
            else:
                logging.error("CHROMER: Max retries exceeded. Failing...")
                raise e

def annotate_fractions(host, frac_data, injection_time=0):
    for i in range(len(frac_data)):
        adjusted_x = frac_data[i][0] - injection_time
        host.axvline(x=adjusted_x, ymin=0.065, ymax=0.0, color='crimson', linewidth=0.5)
        mid_x = adjusted_x if i == len(frac_data) - 1 else (adjusted_x + frac_data[i+1][0] - injection_time) / 2
        host.annotate(str(frac_data[i][1]), xy=(mid_x, 0), xytext=(0, -5), textcoords='offset points', 
                      ha='center', fontsize=12, rotation=90, fontweight='bold')
            
def annotate_peaks(host, x_values, y_values, peaks):
    for i, peak in enumerate(peaks):
        offset = (i % 2) * 10  # prevent overlap of annotations
        host.annotate(f"V: {x_values[peak]:.2f}\nA: {y_values[peak]:.2f}", 
                      (x_values[peak], y_values[peak]), textcoords="offset points", 
                      xytext=(0,5 + offset), ha='center', fontweight='bold', fontsize=10, color='red', 
                      bbox=dict(facecolor='none', edgecolor='red', boxstyle='round,pad=0.2'))

def plot_data(host, x_values, y_values, color='blue', linewidth=5, label='UV 280nm'):
    host.plot(x_values, y_values, color=color, linewidth=linewidth, label=label)
    return host

def save_and_upload_plot(Title, method_folder, method, batch, brain):
    filename = os.path.join(method_folder, f"{Title}.jpg")
    plt.savefig(filename, bbox_inches='tight', dpi=300)
    upload_file(drive, filename, method, batch, brain)

    
def annotate_no_expression(host):
    host.text(0.5, 0.5, 'NO / LOW EXPRESSION', fontsize=50, fontweight='bold', color='red', ha='center', va='center', alpha=0.5, rotation=45, transform=host.transAxes)

def annotate_multiple_peaks(host):
    host.text(1, 1.05, 'MULTI', ha='right', va='top', transform=host.transAxes, fontsize=24, fontweight='bold', color='black')

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
    for i, (frac_range, end_frac_index) in enumerate(frac_ranges):
        host.text(1, 1.05, frac_range, ha='right', va='top', transform=host.transAxes, fontsize=24, fontweight='bold', color='black')

def chromeunicorns(zip_file, udata, run_folder, brain):    
    udata_INFO = process_chrom(zip_file, brain)
    udata = udata_INFO['udata']
    method = udata_INFO['method']
    batch = udata_INFO['batch']
    Title = udata_INFO['title']
    date = udata_INFO['date']
    
    
    if 'Fractions' not in udata or None in [batch, method, Title] or "" in [batch, method, Title]:
        print(f"CHROMER: Data for {zip_file} is invalid. Invalid: FRACTIONS - {udata.get('Fractions') is None}, BATCH - {batch is None}, METHOD - {method is None}. Skipping...")
        return


    plt.figure(figsize=(15, 10), edgecolor='black')
    sns.set_style("whitegrid")
    host = host_subplot(111, axes_class=AA.Axes)

    injection_time = 0
    if method == 'SEC':
        injection_time = udata.get('Injection', {}).get('data', [[0, None]])[0][0]

    for key, value in udata.items():
        if key in ["UV 1_280", "UV"]:
            x_values = [x[0] - injection_time for x in value['data']] 
            y_values = [y[1] for y in value['data']]
            plot_data(host, x_values, y_values)
    
    
    if 'Fractions' in udata and method != "SEC":
        xlim_min = udata['Fractions']['data'][0][0]
        xlim_max = udata['Fractions']['data'][-1][0]
    elif method == "SEC":
        xlim_min = 0
        xlim_max = 30
            
    host.set_xlim(xlim_min, xlim_max)        
    y_values_in_xlim = [y for x, y in zip(x_values, y_values) if xlim_min <= x <= xlim_max]
    if y_values_in_xlim:
        ymax = max(y_values_in_xlim)
        host.set_ylim(-0.5, ymax * 1.05)
    else:
        logging.warning("No y_values within specified x limits. Skipping...")
    
    peaks, properties = find_peaks(y_values, prominence=5, width=1, height=10)
    peaks = [peak for peak in peaks if xlim_min <= x_values[peak] <= xlim_max]
    annotate_peaks(host, x_values, y_values, peaks)

    if 'Fractions' in udata:
        annotate_fractions(host, udata['Fractions']['data'], injection_time)

    if not peaks:
        annotate_no_expression(host)
    elif len(peaks) > 1 and method != "SEC" and method != "IMAC":
        annotate_multiple_peaks(host)
    else:
        if 'Fractions' in udata and method != "SEC" and method != "IMAC":
            frac_ranges = get_fraction_ranges(x_values, y_values, peaks, udata['Fractions']['data'])
            annotate_fraction_ranges(host, frac_ranges)
            
    host.set_title(Title, loc='left', fontsize=24, fontweight='bold', color='black', pad=20)
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
    save_and_upload_plot(Title, run_folder, method, batch, brain)
    plt.close()
    
def process_file(root, run_folder, brain, file):
    if file.endswith(".Result"):
        logging.info(f"CHROMER: Processing {file}")
        result_file = os.path.join(root, file)
        udata_INFO = process_chrom(result_file, brain)
        udata = udata_INFO['udata']
        method = udata_INFO['method']
        batch = udata_INFO['batch']
        Title = udata_INFO['title']
        
        if batch is None or method is None or batch == "" or method == "" or "None" in Title:
            logging.warning(f"CHROMER: Batch# or Method is missing or invalid. Skipping... {file}")
        else:
            method_folder = os.path.join(run_folder, method)
            os.makedirs(method_folder, exist_ok=True)
            chromeunicorns(result_file, udata_INFO['udata'], run_folder, brain)
            os.remove(result_file)  # Delete the processed file | TODO: This doesn't actually work?
    elif file.endswith(".UFol"):
        tar_file = os.path.join(root, file)
        with tarfile.open(tar_file, 'r') as tar:
            tar.extractall(path=run_folder)
            for root, dirs, files in os.walk(run_folder):
                for file in files:
                    if file.endswith(".Result"):
                        logging.info(f"CHROMER: Processing {file}")
                        result_file = os.path.join(root, file)
                        udata_INFO = process_chrom(result_file, brain)
                        udata = udata_INFO['udata']
                        method = udata_INFO['method']
                        batch = udata_INFO['batch']
                        Title = udata_INFO['title']
                        
                        if batch is None or method is None or batch == "" or method == "" or "None" in Title:
                            logging.warning(f"CHROMER: Batch# or Method is missing or invalid. Skipping... {file}")
                        else:
                            method_folder = os.path.join(run_folder, method)
                            os.makedirs(method_folder, exist_ok=True)
                            chromeunicorns(result_file, udata_INFO['udata'], run_folder, brain)
                            os.remove(result_file)  # Delete the processed file
        
if __name__ == "__main__":
    unicorns = config['paths']['unicorns']
    run_folder = os.path.join(config['paths']['run_folder_base'], datetime.datetime.now().strftime("%Y%m%d_%H%M%S/"))
    os.makedirs(run_folder, exist_ok=True)

    setup_logger()
    pylab.rcParams.update(params)

    try:
        construct_recognition(queue, config['paths']['brain'])
        brain = enable_cognition(config['paths']['brain'])
    except Exception as e:
        logging.warning(f"[UPDATES NOT FETCHED. ATTEMPTING TO LOAD MOST RECENT RECORD.]\n    REASON:{e}.", exc_INFO=True)
        try:
            brain = enable_cognition(config['paths']['brain'])
        except Exception as e:
            logging.error(f"[INDEX FAILED TO UPDATE AND FALLBACK WAS MISSING OR INVALID.]\n    REASON:{e}", exc_INFO=True)
            
    all_files = []
    processed_files = set()  # set to keep track of processed files
    for root, dirs, files in os.walk(unicorns):
        for file in files:
            if file.endswith(".UFol") or file.endswith(".Result"):
                full_path = os.path.join(root, file)  # get the full path of the file
                if full_path not in processed_files:  # check if the full path of the file has already been processed
                    all_files.append((root, run_folder, brain, file))
                    processed_files.add(full_path)  # add the full path of the file to processed_files set | TODO: This doesn't actually work either?

    for file in all_files:
        process_file(*file)
    
    logging.info(f"\n\nCHROMER: Finished processing {len(all_files)} files at {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
