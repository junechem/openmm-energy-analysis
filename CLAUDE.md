# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This repository contains `omm_energy.py`, an OpenMM Energy Analysis Tool that mimics the functionality of GROMACS's `gmx energy` command for analyzing energy data from OpenMM molecular dynamics simulations.

## Core Architecture

The tool is structured as a single Python script with several key functional components:

### Data Processing Pipeline
1. **File Reading** (`read_openmm_log`): Parses OpenMM log files with tab-separated data and quoted column headers
2. **Time Filtering** (`filter_by_time`): Applies time range constraints to the dataset
3. **Interactive Selection** (`display_column_options`, `get_user_selection`): Provides GROMACS-style interactive column selection
4. **Statistical Analysis** (`calculate_statistics`, `print_statistics`): Computes averages, error estimates, RMSD, and drift
5. **Fluctuation Properties** (`calculate_isothermal_compressibility`, `print_fluctuation_properties`): Calculates thermodynamic properties from fluctuations
6. **Output Generation** (`generate_xvg_files`): Creates GROMACS-compatible .xvg time-series files

### Key Features
- **GROMACS Interface Compatibility**: Mimics `gmx energy` user experience with numbered column selection
- **Time Window Analysis**: `-b` and `-e` flags for begin/end time filtering
- **Per-molecule Scaling**: `-nmol` flag for normalizing energy quantities
- **Fluctuation Properties**: `-fluct_props` flag calculates isothermal compressibility using `Kt = (<V²> - <V>²)/(kbT*<V>)`
- **XVG File Generation**: Automatically creates time-series data files for each selected property

### Input/Output Format
- **Input**: OpenMM log files (default: `md.log`) with tab-separated columns and quoted headers starting with `#`
- **Output**: Statistical summaries and `.xvg` files with GROMACS-compatible format

## Dependencies

- Python 3.x
- pandas: For data manipulation and CSV parsing
- numpy: For statistical calculations and array operations
- Standard library: argparse, sys, re, os, typing

## Usage Commands

```bash
# Run with conda environment (pandas dependency)
source ~/miniconda/etc/profile.d/conda.sh && conda activate base && python3 omm_energy.py

# Basic analysis
python3 omm_energy.py

# Time window analysis
python3 omm_energy.py -b 100 -e 500

# Per-molecule quantities
python3 omm_energy.py -nmol 343

# Fluctuation properties calculation
python3 omm_energy.py -fluct_props

# Combined analysis
python3 omm_energy.py -b 100 -e 500 -nmol 343 -fluct_props
```

## Testing

The script can be tested interactively by running with sample input:
```bash
echo -e "5\n6\n0" | python3 omm_energy.py -b 100 -e 200 -fluct_props
```
This selects Temperature (column 5) and Box Volume (column 6), enabling fluctuation property calculations.

## Development Notes

- **Time Column Detection**: The script automatically detects time columns by searching for keywords ('time', 'ps', 'picosecond', 'step', 'frame') in column names. Time is treated like any other column - no special treatment in the interface.
- **Time Filtering**: When no time column is found, `-b` and `-e` options are ignored with a warning message.
- **XVG File Generation**: Uses actual time values when available, or frame numbers (starting from 1) when no time column exists.
- **Energy Scaling**: Energy columns are automatically scaled when `-nmol` is provided.
- **Fluctuation Properties**: Require both temperature and volume columns to be selected.
- **Filename Generation**: XVG files use cleaned property names (spaces → underscores, units removed).
- **Error Handling**: Includes file not found, parsing errors, and missing time column warnings.