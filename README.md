# OpenMM Energy Analysis Tool

A Python tool that mimics the functionality of GROMACS's `gmx energy` command for analyzing energy data from OpenMM molecular dynamics simulations.

## Overview

This tool provides a familiar GROMACS-style interface for analyzing energy trajectories from OpenMM simulations. It reads OpenMM log files with tab-separated data and quoted column headers, allowing interactive selection of energy terms for statistical analysis and visualization.

## Features

- **GROMACS Interface Compatibility**: Interactive numbered column selection mimicking `gmx energy`
- **Time Window Analysis**: Filter data by time range using `-b` (begin) and `-e` (end) flags
- **Per-molecule Scaling**: Normalize energy quantities using the `-nmol` flag
- **Statistical Analysis**: Calculate averages, error estimates, RMSD, and drift for selected properties
- **Fluctuation Properties**: Calculate isothermal compressibility from volume and temperature fluctuations
- **XVG File Generation**: Generate GROMACS-compatible time-series data files
- **Automatic Time Detection**: Smart detection of time columns in various formats

## Installation

### Requirements

- Python 3.x
- pandas
- numpy

### Setup

1. Clone this repository:
```bash
git clone https://github.com/yourusername/openmm-energy-analysis.git
cd openmm-energy-analysis
```

2. Install dependencies:
```bash
# Using conda (recommended)
conda install pandas numpy

# Or using pip
pip install pandas numpy
```

## Usage

### Basic Usage

```bash
# Basic analysis with interactive column selection
python3 omm_energy.py

# Analyze specific log file
python3 omm_energy.py -f simulation.log

# With conda environment
source ~/miniconda/etc/profile.d/conda.sh && conda activate base && python3 omm_energy.py
```

### Advanced Options

```bash
# Time window analysis (analyze data from 100 to 500 ps)
python3 omm_energy.py -b 100 -e 500

# Per-molecule quantities (normalize by number of molecules)
python3 omm_energy.py -nmol 343

# Calculate fluctuation properties (isothermal compressibility)
python3 omm_energy.py -fluct_props

# Combined analysis
python3 omm_energy.py -b 100 -e 500 -nmol 343 -fluct_props
```

### Interactive Testing

For automated testing with predefined selections:
```bash
# Select Temperature (column 5) and Box Volume (column 6), then exit
echo -e "5\n6\n0" | python3 omm_energy.py -b 100 -e 200 -fluct_props
```

## Input Format

The tool expects OpenMM log files with:
- Tab-separated values
- Column headers starting with `#` and enclosed in quotes
- Time data (automatically detected from column names containing: 'time', 'ps', 'picosecond', 'step', 'frame')

Example input format:
```
#"Time (ps)"	"Potential Energy (kJ/mol)"	"Temperature (K)"	"Box Volume (nm^3)"
0.000	-12345.67	298.15	125.678
1.000	-12340.23	299.82	125.701
...
```

## Output

### Statistical Analysis
- **Average**: Mean value over the selected time range
- **Err.Est.**: Standard error of the mean
- **RMSD**: Root mean square deviation (standard deviation)
- **Tot-Drift**: Total drift (difference between first and last values)

### Fluctuation Properties
When both temperature and volume are selected with `-fluct_props`:
- **Isothermal Compressibility**: κ = (⟨V²⟩ - ⟨V⟩²)/(k_BT⟨V⟩)
- **Bulk Modulus**: Reciprocal of compressibility
- **Volume per Mole**: Average system volume

### XVG Files
Time-series data files compatible with GROMACS plotting tools:
- Automatic filename generation from property names
- GROMACS-style headers with metadata
- Two-column format: time/frame vs. property value

## Examples

### Example 1: Basic Energy Analysis
```bash
python3 omm_energy.py -f md.log
# Select columns interactively:
# 2  (Potential Energy)
# 3  (Kinetic Energy)
# 0  (exit)
```

### Example 2: Temperature and Pressure Analysis with Time Window
```bash
python3 omm_energy.py -f md.log -b 50 -e 200
# Analyze only data from 50 ps to 200 ps
```

### Example 3: Fluctuation Properties Analysis
```bash
python3 omm_energy.py -f md.log -fluct_props
# Select both temperature and volume columns to calculate compressibility
```

## File Structure

```
├── omm_energy.py          # Main analysis script
├── CLAUDE.md             # Development guidance for Claude Code
├── README.md             # This file
└── *.xvg                 # Generated time-series files (after running)
```

## Core Functions

- `read_openmm_log()`: Parse OpenMM log files with proper header handling
- `find_time_column()`: Automatically detect time-based columns
- `filter_by_time()`: Apply time range filtering
- `display_column_options()`: Show interactive column selection menu
- `calculate_statistics()`: Compute statistical properties
- `calculate_isothermal_compressibility()`: Calculate thermodynamic properties
- `generate_xvg_files()`: Create GROMACS-compatible output files

## Technical Details

### Time Column Detection
The tool automatically identifies time columns by searching for keywords in column names:
- 'time', 'ps', 'picosecond'
- 'step', 'frame'

### Energy Scaling
When using `-nmol`, energy-related columns (containing 'Energy' in the name) are automatically divided by the number of molecules for per-molecule quantities.

### Fluctuation Properties
Isothermal compressibility calculation requires:
1. Volume data (column name containing 'volume' and 'nm^3')
2. Temperature data (column name containing 'temperature')
3. Both columns must be selected for calculation

## Limitations

- Currently supports only tab-separated OpenMM log files
- Fluctuation properties require specific column naming conventions
- Error analysis is basic (standard error only)
- No block averaging for error estimation yet

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

## License

This project is open source. Please see the LICENSE file for details.

## Acknowledgments

- Inspired by GROMACS `gmx energy` tool
- Built for OpenMM simulation analysis
- Compatible with standard molecular dynamics analysis workflows

## Support

For issues, questions, or contributions, please use the GitHub issue tracker.