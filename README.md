# PprA Interface Analyzer

This directory contains scripts to analyze residue-residue contacts in PprA interfaces from various crystal structures and visualize them as schematic diagrams.

## Workflow

The analysis consists of two steps:

1.  **Identify Contacts:** Run the PyMOL script to find contacts across different PDB models and generate a CSV report.
2.  **Generate Schematics:** Run the visualization script to draw vector graphics (SVG) based on the CSV data.

## 1. Identify Contacts

The script `find_pprA_interface_contacts.py` loads multiple PprA crystal structures, aligns them, and identifies ionic, hydrogen bond, hydrophobic, and other interactions across P, F, S, and M interfaces.

### Requirements
- PyMOL (must be run within PyMOL or with a python environment that has `pymol` installed)
- pandas

### Usage
From within this directory:
```bash
pymol find_pprA_interface_contacts.py
```
Or from an interactive PyMOL session:
```python
run find_pprA_interface_contacts.py
```

**Output:**
- `pprA_residue_contacts.csv`: A CSV file containing all identified contacts.

*Note: You can edit parameters like model sources, distance cutoffs, and output filenames directly within the Python script.*

## 2. Generate Schematics

The script `draw_pprA_contact_schemiatic.py` reads the CSV file generated in the previous step and creates schematic diagrams of the interactions.

### Requirements
- python 3
- pycairo (`pip install pycairo`)
- pandas

### Usage
From within this directory:
```bash
python draw_pprA_contact_schemiatic.py -i pprA_residue_contacts.csv -o ./output_plots
```
- `-i`: Path to the input CSV file (required).
- `-o`: Path to the output directory (optional, defaults to `./`).

**Output:**
- SVG files for different interfaces (e.g., `intra_contacts.svg`, `P_contacts.svg`, `F_contacts.svg`, `S_contacts.svg`) saved to the specified output directory.

*Note: Visualization settings such as colors, dimensions, and residue numbering references can still be modified directly within the Python script.*