# **Project AKTAmercy - "CHROMER"**

[![VER](https://img.shields.io/badge/VERSION-1.0.0-darkgreen.svg)] |
[![DEV](https://img.shields.io/badge/DEVELOPMENT-MAINTENANCE-blue.svg)] |
[![STATUS](https://img.shields.io/badge/STATUS-STABLE-darkgreen.svg)]

CHROMER is a Python automaton that post-processes chromatographic data from AKTA Pure HPLC systems running UNICORN (7.0+). It compiles `.UFol` and `.Result` files into annotated JPG chromatograms locally.

## Features

- **PyCORN-powered Parsing**: Parses proprietary `.UFol` (*UNICORNS*) archives from UNICORN Evaluation software, using a PyCORN-inspired parser by [Yasar L. Ahmed](https://github.com/pyahmed)
- **Chromatogram Recognition**: Parses sample name, purification method, and run date from UNICORN exports. Optionally enriches titles with construct metadata via a local index (`brain.json`)
- **Peak Detection**: Detects peaks in the chromatogram and determines pool fractions based on peak area
- **Local Output**: Writes compiled chromatograms to timestamped folders under `./data/DONE/`

## Developmental Roadmap

<details>
<summary>
<b> 1.0.0 </b> :white_check_mark: 
</summary> 

```diff
Advancements
+ SEC Chromatograms generated are now numerically accurate, accounting for flow rate and injection point.
+ Index algorithm now covers a majority of targets
+ Debug tools <b>?INSPECT</b> and <b>parseLOG</b> added to identify deviants easily.
+ CHROMER -> CHROMER

Regressions
- Multiprocessing removed for simplicity, potential for reimplementation at later release
- AFFINITY chromatograms are visually congruent, but numerically (x values / volume) incorrect.
```
</details>

## Known Issues and Limitations

- **Sample metadata**: Sample_ID, method, and date are read from the UNICORN run log. Without `brain.json`, CHROMER runs in generic mode and titles chromatograms as `{METHOD}_{SAMPLE}`. With `brain.json`, indexed mode adds construct names to titles when the sample matches the index; unmatched samples still render with a generic title.
- **Duplicate Logging**: Even without multiprocessing, the log generates duplicate lines for each sample. `parseLOG` is required to process logs when looking for deviants.
- **Pooling Fractions (AFFINITY)**: The implementation is incomplete, and lacks the ability to enumerate multiple peaks. Some of the pooling ranges are a bit too conservative.
- **Plotting**: Revisions to the plotting are likely, to include more information and to make the plots more readable.

Please report any issues you encounter [here](https://github.com/alxdolphin/AKTAmercy/issues), and feel free to contribute to the project by submitting a [pull request](https://github.com/alxdolphin/AKTAmercy/pulls).

## Usage

To use CHROMER, you need Python installed on your system. A virtual environment via [Anaconda](https://www.anaconda.com/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) is recommended.

```bash
# Create and activate the environment
conda env create -f environment.yml
conda activate CHROMER

# Clone the repository
git clone https://github.com/alxdolphin/AKTAmercy.git
cd AKTAmercy

# Place UNICORN files in the drop-off directory
mkdir -p ./data/DROP-OFF

# Run the script
./CHROMER.py
```

The script processes all `.UFol` and `.Result` files in `./data/DROP-OFF/` and writes compiled chromatograms to `./data/DONE/<timestamp>/`, organized by purification method.

## Construct Index (`brain.json`, optional)

CHROMER works without `brain.json`. When the file is present at `./brain.json` (or the path in `config.json`), CHROMER runs in **indexed mode** and maps batch numbers to construct metadata. When absent, **generic mode** still parses exports and renders chromatograms.

Each index entry uses the batch number (uppercase) as the key:

```json
{
    "BATCH01": {
        "ConstructID": "MyConstruct_pVax",
        "Ext/1000": "210",
        "MW/1000": "150.00"
    }
}
```

| Field | Description |
| --- | --- |
| `ConstructID` | Construct name used in chromatogram titles |
| `Ext/1000` | Extinction coefficient (×1000) |
| `MW/1000` | Molecular weight (×1000) |

Maintain `brain.json` manually or export from your lab's tracking system. CHROMER loads this file at startup when present and does not fetch remote updates. Samples not found in the index still produce chromatograms with generic titles.

## Configuration

`config.json` controls local paths and plotting parameters:

```json
{
    "paths": {
        "unicorns": "./data/DROP-OFF/",
        "run_folder_base": "./data/DONE/",
        "brain": "./brain.json"
    },
    "plotting": { "params": { ... } }
}
```

## State of the Art - UNICORN vs CHROMER

`v1.0.0-rc1`

| | UNICORN | CHROMAUTOGRAM |
|:---:|:---:|:---:|
| **Lectin** | <sub><img src="data/DONE/_sample/UNICORN_LEC.png" width="400"></sub> | <sub><img src="data/DONE/_sample/CHROMER_LEC.jpg" width="400"></sub> |

<center>
The X values are slightly off in the CHROMAUTOGRAM, but the plot is otherwise identical to the UNICORN plot, with the peak falling within the same fraction range and reaching the same height.
</center>

| | UNICORN | CHROMAUTOGRAM |
|:---:|:---:|:---:|
| **SEC** | <sub><img src="data/DONE/_sample/UNICORN_SEC.png" width="400"></sub> | <sub><img src="data/DONE/_sample/CHROMER_SEC.jpg" width="400"> </sub> |

<center>
Features like peak area shading, an overview plot, and more are planned in future releases.

## Resources and Acknowledgements

</center>

| Category | Libraries/APIs |
| --- | --- |
| **Data Handling** | [PyCORN](https://github.com/pyahmed/PyCORN) |
| **Python Standard Library** | [datetime](https://docs.python.org/3/library/datetime.html), [io](https://docs.python.org/3/library/io.html), [json](https://docs.python.org/3/library/json.html), [logging](https://docs.python.org/3/library/logging.html), [os](https://docs.python.org/3/library/os.html), [re](https://docs.python.org/3/library/re.html), [struct](https://docs.python.org/3/library/struct.html), [tarfile](https://docs.python.org/3/library/tarfile.html), [xml.etree.ElementTree](https://docs.python.org/3/library/xml.etree.elementtree.html), [collections.OrderedDict](https://docs.python.org/3/library/collections.html#collections.OrderedDict), [zipfile](https://docs.python.org/3/library/zipfile.html) |
| **Data Visualization** | [numpy](https://numpy.org/), [matplotlib](https://matplotlib.org/stable/api/pyplot_summary.html), [mpl_toolkits](https://matplotlib.org/stable/api/artist_api.html), [seaborn](https://seaborn.pydata.org/), [scipy.signal](https://docs.scipy.org/doc/scipy/reference/signal.html) |
