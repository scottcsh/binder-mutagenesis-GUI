# Binder Mutagenesis GUI

A lightweight PySide6 desktop GUI for ProteinMPNN-based binder mutagenesis screening with structural scoring.

The GUI is designed for Rocky Linux workstation/server environments and focuses on generating residue-range variants from one input `.pdb` file, then scoring every generated variant with:

- **ProteinMPNN**
- **EvoEF2 ComputeBinding**
- **PDBe-KB PISA**
- **Mutation count summary for improved variants**

---

## Table of Contents
- [Installation](#installation)
- [System Requirements](#system-requirements)
- [First-run Configuration](#first-run-configuration)
- [Running the App](#running-the-app)
- [Input and Output](#input-and-output)
- [Mutagenesis Workflow](#mutagenesis-workflow)
- [Result Tables](#result-tables)
- [Project Layout](#project-layout)
- [Acknowledgements](#acknowledgements)

---

## Installation

Clone or upload this project directory to your Rocky Linux machine.

Create a virtual environment and install Python dependencies:

```bash
cd binder_mutagenesis_GUI
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

If your system does not already provide virtualenv support:

```bash
sudo dnf install -y python3 python3-pip python3-virtualenv
```

---

## System Requirements

This project depends on external tools that are not installed through `pip`.

Required external programs and services:

- ProteinMPNN repository
- Bash
- EvoEF2 executable
- PDBe-KB PISA REST access

Python packages are listed in:

```text
requirements.txt
```

Current Python dependencies:

- PySide6
- Biopython
- requests

---

## First-run Configuration

Open the GUI and click **Settings**.

Configure the following paths and parameters:

### ProteinMPNN

- `ProteinMPNN directory`
  - Path to the ProteinMPNN repository.
  - The GUI expects these scripts under that directory:
    - `helper_scripts/parse_multiple_chains.py`
    - `helper_scripts/assign_fixed_chains.py`
    - `helper_scripts/make_fixed_positions_dict.py`
    - `protein_mpnn_run.py`
- `ProteinMPNN env setup`
  - Optional shell setup commands run before ProteinMPNN.
  - Example:

```bash
module load miniforge3/25.3.0-3
conda activate proteinmpnn
```

- `ProteinMPNN seed`
  - Random seed passed to ProteinMPNN.
- `ProteinMPNN batch size`
  - Batch size passed to ProteinMPNN.

### EvoEF2

- `EvoEF2 executable`
  - Path or command name for EvoEF2.
- `EvoEF2 threads`
  - Value passed through `OMP_NUM_THREADS`.

Settings are stored locally in:

```text
resources/settings.json
```

This file is intentionally ignored by Git because it contains machine-specific paths.

---

## Running the App

From the project root:

```bash
source .venv/bin/activate
python main.py
```

If Qt reports DBus warnings in SSH/X11 or non-GNOME sessions, the app may still work normally. If the GUI does not open, try:

```bash
dbus-run-session -- python main.py
```

---

## Input and Output

Input file:

- `.pdb`

The GUI requires:

- Input PDB
- Target chain
- Design chain
- Residue range to screen
- ProteinMPNN replicate count

Output directory is optional. If it is blank, the GUI writes output under:

```text
<input_pdb_directory>/_binder_mutagenesis_gui/
```

Each run creates:

```text
<output_directory>/<input_pdb_stem>_mutagenesis/
```

Important output files:

- `mutagenesis_results.csv`
  - One row per generated variant PDB.
- `mutation_counts.csv`
  - Mutation count summary for variants that improve at least one score.
- `variants/`
  - Coordinate-preserving variant PDB files generated from ProteinMPNN sequences.
- `input_pisa/`
  - PISA input and output files for the original input PDB.
- `variant_pisa/`
  - PISA input and output files for generated variants.
- `proteinmpnn/`
  - ProteinMPNN outputs and generated run script.

---

## Mutagenesis Workflow

### Input scoring

The GUI first reports an **Input Result** table.

For the original input PDB, it reports:

- Input PDB filename
- Requested design-chain residue range
- Native sequence for that residue range
- EvoEF2 binding score
- PISA interface area
- PISA interface solvation energy
- PISA p-value

### ProteinMPNN

ProteinMPNN is run once for the whole input PDB.

The requested design-chain residue range is treated as the mutable region. Other design-chain residues are fixed through ProteinMPNN fixed positions.

ProteinMPNN is called with sampling temperature:

```text
0.4
```

The replicate count is controlled by the GUI input field. The default is:

```text
5000
```

### Variant PDB generation

For every ProteinMPNN replicate sequence, the GUI writes one coordinate-preserving variant PDB.

Only residue names in the requested design-chain range are changed. Coordinates and all other PDB text are preserved from the input PDB.

### EvoEF2

EvoEF2 is run for the input PDB and every variant PDB:

```bash
EvoEF2 --command=ComputeBinding --pdb=<variant.pdb>
```

The reported value is parsed from the `Total = ...` line and keeps the original sign.

### PISA

PISA uses **target + design chains**.

For each structure:

1. Extract target and design chains into a prepared PDB.
2. Submit the prepared PDB to the PDBe-KB PISA REST workflow.
3. Parse `PISA_interface.csv`.
4. Use the interface with the lowest `int_solv_energy`.

Reported columns:

- `Interface area`
- `Interface solvation E`
- `p-value`

PISA jobs are submitted with concurrency 5.

### Mutation summary

After all replicate scores are complete, the GUI compares each variant to the input PDB.

A variant is considered improved if at least one of these conditions is true:

- EvoEF2 binding is lower than the input value.
- Interface area is higher than the input value.
- Interface solvation energy is lower than the input value.

For improved variants, mutations in the requested residue range are counted in this format:

```text
A123V
```

Mutation counts are sorted from highest to lowest.

---

## Result Tables

### Input Result

One row of input baseline values:

- `Input PDB`
- `Residues`
- `Sequence`
- `EvoEF2`
- `PISA Area`
- `PISA Solv E`
- `PISA p-value`

### Replication Output

One row per generated variant PDB:

- `variant pdb`
- `sequence`
- `EvoEF2 binding`
- `Interface area`
- `Interface solvation E`
- `p-value`
- `status`

The GUI displays only the variant PDB filename. The full path is still written to `mutagenesis_results.csv`.

### Mutagenesis Result

Mutation count summary:

- `mutation`
- `count`
- `EvoEF2 improved`
- `Interface area improved`
- `Solvation E improved`

---

## Project Layout

```text
binder_mutagenesis_GUI/
  app/
    main_window.py
    services/
      processor.py
  resources/
    settings.json
  scripts/
    pisa_batch_dir_to_csv.py
  main.py
  README.md
  requirements.txt
```

`resources/settings.json` is generated locally and ignored by Git.

---

## Acknowledgements

This GUI is a thin interface layer around external tools and services including:

- ProteinMPNN
- EvoEF2
- PDBe-KB PISA

Please follow the licenses, citation policies, and usage requirements of each upstream project.

</br>
</br>

[Return to top](#table-of-contents)
