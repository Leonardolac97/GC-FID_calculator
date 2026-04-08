# GC-FID HDO Calculator GUI

A PyQt6 desktop application for extracting peak tables from Agilent-style GC-FID PDF reports, filtering compounds by retention time, estimating reactant conversion, applying external calibration files, and generating time-on-stream plots and CSV summaries for hydrodeoxygenation (HDO) datasets.

This tool reads the **Front Signal Results** section from GC report PDFs, filters peaks against a user-defined retention-time table, and computes derived quantities such as:

- reactant wt.% and conversion
- calibrated wt.%
- mass flow (`g/min`)
- molar flow (`mol/min`)
- mol percentages
- grouped selectivity trends
- Pearson correlation maps
- carbon recovery over time

The current GUI contains these tabs:

1. **Retention time filter**
2. **Product trend**
3. **Product calibration**
4. **Mol percentage over time**
5. **Selectivity**
6. **Pearson correlation chart**
7. **Carbon balance**

The workflow, dependencies, default reaction settings, and file-loading behavior are defined directly in `GC_FID_HDO_calculator.py`. The script imports `pandas`, `PyMuPDF` (`fitz`), `numpy`, `periodictable`, `seaborn`, `matplotlib`, `scipy`, and `PyQt6`. The GUI also expects reaction settings to be loadable from JSON, retention tables from CSV/TXT, and calibration tables from CSV/TXT. fileciteturn2file0

---

## Main features

- Import multiple GC-FID PDF reports
- Extract the **Front Signal Results** peak table from each PDF
- Filter peaks by a custom retention-time table and tolerance
- Estimate reactant wt.% and conversion from a linear calibration
- Apply external calibration for multiple products
- Compute `wt%`, `MW`, `g/min`, `mol/min`, `mol%`, `mol%_without_solv`, `mol%_normalized`, and `gC/min`
- Generate grouped product trends and grouped selectivity plots
- Build Pearson correlation heatmaps from normalized molar fractions
- Compute carbon recovery and export summaries automatically

The script creates `CSV_Files` and `Plot_Summary` folders next to the imported PDF files, and writes outputs such as `complete.csv`, `filtered.csv`, `conversion.csv`, `calibrated.csv`, `selectivity_group_results.csv`, `carbon_recovery.csv`, plus several PNG figures. fileciteturn2file0

---

## Installation

### Option 1 — Conda / Mamba

```bash
conda env create -f environment.yml
conda activate gc-fid-hdo-gui
python GC_FID_HDO_calculator.py
```

### Option 2 — pip + virtual environment

```bash
python -m venv .venv
```

#### Windows

```bash
.venv\Scripts\activate
pip install --upgrade pip
pip install -r requirements.txt
python GC_FID_HDO_calculator.py
```

#### macOS / Linux

```bash
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
python GC_FID_HDO_calculator.py
```

---

## Requirements

- Python 3.10 or newer recommended
- Desktop environment capable of running PyQt6
- GC-FID PDF reports that contain a parsable **Front Signal Results** section

The script uses these libraries directly in the source code: `pandas`, `fitz` / PyMuPDF, `numpy`, `periodictable`, `seaborn`, `json`, `matplotlib`, `scipy.stats`, and `PyQt6`. fileciteturn2file0

---

## Running the GUI

```bash
python GC_FID_HDO_calculator.py
```

When the GUI opens, the usual order is:

1. Load a **retention time file**
2. Load a **reaction configuration**
3. Import your **GC-FID PDFs**
4. Click **Save CSVs**
5. Load a **calibration file**
6. Click **Apply Calibration and Save CSV**
7. Use the other tabs to generate plots and summaries

The script starts a `QApplication`, opens the `RetentionTimeFilterApp` main window, and runs the GUI event loop directly from `main()`. fileciteturn2file0

---

## Input files

This is the part that must be very clear for users. The GUI depends on three user-provided reference files plus the GC-FID PDF reports.

### 1. GC-FID PDF reports

The application imports one or more PDF reports and extracts the **Front Signal Results** section. The parser expects a table with five values per row:

- `Retention Time`
- `Area`
- `Area %`
- `Height`
- `Height %`

This is hard-coded in the PDF parser and the script builds a dataframe with exactly those columns after extraction. fileciteturn2file0

A sample uploaded PDF shows the expected structure under **Front Signal Results**, with rows such as retention time, area, area percent, height, and height percent. fileciteturn1file6

**Important assumptions about the PDF filenames:**

- the script strips `.pdf` to build `Sample Name`
- it later converts `Sample Name` to `float` for time-on-stream plots
- therefore, filenames should ideally be numeric, for example:

```text
0.pdf
1.pdf
2.pdf
4.pdf
8.pdf
24.pdf
```

or at least contain a clean numeric part that can be interpreted as the time on stream. The import routine sorts files using the numeric part extracted from the filename. fileciteturn2file0

If filenames are not numeric enough, plots that depend on `Sample Name.astype(float)` may fail or behave unpredictably. That limitation comes from the current implementation. fileciteturn2file0

### 2. Retention time file

The GUI can load a retention table from `.csv` or `.txt`.

#### Accepted formats

**TXT format**
- tab-separated
- no header required
- two columns only:
  1. retention time
  2. compound name

Example:

```text
9.386	Formaldehyde
9.520	Methanol
20.736	Guaiacol
21.208	Decalin cis
```

The current example file follows exactly this format. fileciteturn1file10

**CSV format**
- must contain the columns:
  - `Retention Time`
  - `Compound`

If those column names are missing, the GUI raises an error while loading the table. fileciteturn2file0

### 3. Calibration file

The calibration file is used in the **Product calibration** tab.

#### Accepted formats

**TXT format**
- tab-separated
- no header required
- five columns, in this order:
  1. `RT`
  2. `Compound`
  3. `Formula`
  4. `Slope`
  5. `Intercept`

Example rows from the uploaded file:

```text
9.386	Formaldehyde	CH2O	263521	-239
20.736	Guaiacol	C7H8O2	1064289.148	136790.506
30.608	2,3,4,5,6-Pentamethylphenol	C11H16O	1596360.094	20749.8884
```

The example calibration file follows that exact structure. fileciteturn1file9

**CSV format**
- must contain the columns:
  - `RT`
  - `Compound`
  - `Formula`
  - `Slope`
  - `Intercept`

If any of those required column names are missing, calibration loading fails. fileciteturn2file0

### 4. Reaction configuration JSON

The reaction configuration is loaded from a JSON file and overwrites the default values stored in the script.

The current JSON structure supports:

- `tolerance`
- `calib_flow`
- `experiment_flow`
- `slope_calibration`
- `linear_co_calibration`
- `reactant_name`
- `reactant_rt`
- `reactant_formula`
- `reactant_molar_mass`
- `reactant_nC`
- `reactant_wt`
- `solvents` as a list of dictionaries
- `key_products`
- `groups`

The uploaded example JSON includes two solvent entries, each with `name`, `rt`, `molar_mass`, and `wt_fraction`. The code reads `solvents` directly from JSON, so multiple solvents are supported as long as each entry follows that dictionary structure. fileciteturn1file4

Example:

```json
{
  "tolerance": 0.09,
  "calib_flow": 0.06,
  "experiment_flow": 0.06,
  "reactant_name": "Guaiacol",
  "reactant_rt": 20.736,
  "reactant_formula": "C7H8O2",
  "reactant_molar_mass": 124.139,
  "reactant_nC": 7,
  "reactant_wt": 10.0,
  "solvents": [
    {
      "name": "Decalin trans",
      "rt": 20.062,
      "molar_mass": 138.23,
      "wt_fraction": 0.36
    },
    {
      "name": "Decalin cis",
      "rt": 21.208,
      "molar_mass": 138.23,
      "wt_fraction": 0.54
    }
  ]
}
```

The current implementation also checks that:

```text
reactant_wt/100 + sum(solvent wt_fraction) ≈ 1.0
```

If that total is inconsistent, calibration/export is stopped with a warning. fileciteturn2file0

---

## Output structure

After importing PDFs, the script creates two folders next to the PDF files:

```text
<your_pdf_folder>/
├── 0.pdf
├── 1.pdf
├── 2.pdf
├── ...
├── CSV_Files/
└── Plot_Summary/
```

This behavior is hard-coded in `import_pdfs()`, which creates `CSV_Files` and `Plot_Summary` inside the folder containing the imported PDFs. fileciteturn2file0

### Folder diagram after a typical full run

```text
project_folder/
├── 0.pdf
├── 1.pdf
├── 2.pdf
├── reaction_config.json
├── retention_time.txt
├── calibration_file.txt
├── CSV_Files/
│   ├── complete.csv
│   ├── filtered.csv
│   ├── conversion.csv
│   ├── calibrated.csv
│   ├── selectivity_group_results.csv
│   ├── selectivity_group_avg_results.csv
│   └── carbon_recovery.csv
└── Plot_Summary/
    ├── conversion_plot.png
    ├── conversion_zoom_plot.png
    ├── mol_percentage.png
    ├── trend_selectivity_group_plot.png
    ├── avg_selectivity_group_plot.png
    ├── pearson_plot.png
    └── carbon_plot.png
```

These filenames are explicitly written by the script when saving CSV and PNG outputs. fileciteturn2file0

### What each output contains

#### `CSV_Files/complete.csv`
All extracted peaks from the imported PDFs before retention-time filtering. fileciteturn2file0

#### `CSV_Files/filtered.csv`
Only peaks matching the retention table within the configured tolerance, plus reactant wt.% and reactant conversion columns where applicable. fileciteturn2file0

#### `CSV_Files/conversion.csv`
Mean reactant conversion by sample name / time on stream. fileciteturn2file0

#### `CSV_Files/calibrated.csv`
Main processed dataset containing original filtered data plus calibration-derived columns such as `wt%`, `Molecular formula`, `MW`, `g/min`, `mol/min`, `mol%`, `mol%_without_solv`, `mol%_normalized`, `gC/min`, and `Carbon recovery (%)`. fileciteturn2file0

#### `CSV_Files/selectivity_group_results.csv`
Grouped selectivity values over time on stream. fileciteturn2file0

#### `CSV_Files/selectivity_group_avg_results.csv`
Average grouped selectivity values over the full run. fileciteturn2file0

#### `CSV_Files/carbon_recovery.csv`
Carbon recovery summary by sample name / time on stream. fileciteturn2file0

#### `Plot_Summary/*.png`
Automatically exported figures generated from the GUI tabs. Current hard-coded filenames are:

- `conversion_plot.png`
- `conversion_zoom_plot.png`
- `mol_percentage.png`
- `trend_selectivity_group_plot.png`
- `avg_selectivity_group_plot.png`
- `pearson_plot.png`
- `carbon_plot.png` fileciteturn2file0

---

## Repository layout

A clean repository structure for publication could look like this:

```text
gc-fid-hdo-calculator/
├── GC_FID_HDO_calculator.py
├── README.md
├── LICENSE
├── CITATION.cff
├── requirements.txt
├── environment.yml
├── .gitignore
├── examples/
│   ├── retention_time_example.txt
│   ├── calibration_file_example.txt
│   ├── reaction_config_example.json
│   └── sample_report_excerpt.pdf
└── docs/
    └── screenshots/
```

---

## Known limitations

- The parser currently reads only the **Front Signal Results** block from the PDF. It does not use Back Signal or Aux Detector data. fileciteturn2file0
- The application assumes sample names can be converted to numeric time-on-stream values in several plotting functions. fileciteturn2file0
- The code is currently a single-file GUI application with some experiment-specific grouping logic defined directly in the script through `self.key_products` and `self.groups`. fileciteturn2file0
- Pearson plots require a successful calibration step because they depend on `latest_calibrated_df` and `mol%_normalized`. fileciteturn2file0
- The script imports `seaborn` specifically for the Pearson heatmap implementation. fileciteturn2file0

---

## Citation

Please cite the software using the metadata in `CITATION.cff`.

A suggested manuscript-style citation is:

> Almeida de Campos, L. *GC-FID HDO Calculator GUI: a PyQt-based workflow for extraction, calibration, and time-on-stream analysis of GC-FID HDO datasets*.

Update the DOI, version, repository URL, and release date in `CITATION.cff` once the repository is public and a permanent archive has been created.

---

## License

This repository currently includes an **MIT License** as a simple permissive default for code sharing.

If you prefer a more attribution-focused academic license or need institutional constraints, revise `LICENSE` before publishing the repository.

---

## Acknowledgment

This project was developed by **Leonardo Almeida de Campos** at **Karlsruhe Institute of Technology (KIT) – ITCP**. The script header also notes AI-assisted development support. fileciteturn2file0
