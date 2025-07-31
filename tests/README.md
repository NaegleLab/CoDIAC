# CoDIAC Pipeline – Master Test Runner

The testing suite ensures each step of the pipeline runs correctly — validating input parsing, file generation, sequence alignment, and feature extraction.

## How to Run All Tests

To run the full test suite:

```bash
python master_test_runner.py
```

## Step 1 – Testing the Setup _(run_setup_task_tests.py)_
This step verifies that:

Required packages and modules are available

Input files and directories are accessible and correctly formatted

The PDB reference file (temp_pdb_filtered.csv) is correctly read

## Step 2 – Binary Adjacency File Generation _(run_adjfile_tests.py)_
In this step:

Adjacency files (.txt, _BF.txt) and mmCIF (.cif) files are generated from Arpeggio-generated .json files located in Data/

Tests verify:

- Proper file and directory creation

- Polymer entity consistency between the adjacency files, PDB reference file and the original structural .cif file

## Step 3 – Contact Map Feature Extraction _(run_contact_maps_tests.py)_
This stage validates both intra-entity (domain–domain) and inter-entity (domain–ligand) contact maps.

Tests include:

- Generating contact maps from the binary adjacency files

- Validating .fasta and feature output files

- Ensuring sequence alignment and offsets are correctly applied

- Identifying unmodeled regions in the PDB structure

- Confirming that domain boundaries and entity mapping from the PDB reference file are accurate

- Checking whether PTMs are properly located and used to crop sequence windows and that these sequences appear correctly in .fasta files

- Checking contact dictionaries for ligand–domain interactions fall within aligned structure sequence ranges


