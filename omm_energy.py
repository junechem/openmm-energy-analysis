#!/usr/bin/env python3
"""
OpenMM Energy Analysis Tool
A tool to analyze energy data from OpenMM simulations, similar to gmx energy.
"""

import argparse
import sys
import numpy as np
import pandas as pd
import re
import os
from typing import List, Tuple, Optional

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Analyze energy data from OpenMM simulations',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -f md.log
  %(prog)s -f md.log -b 100 -e 500
  %(prog)s -f md.log -nmol 343 -b 100
  %(prog)s -f md.log -fluct_props
        """
    )

    parser.add_argument('-f', '--file', default='md.log',
                        help='Input log file (default: md.log)')
    parser.add_argument('-b', '--begin', type=float, default=None,
                        help='Begin analysis at time (ps)')
    parser.add_argument('-e', '--end', type=float, default=None,
                        help='End analysis at time (ps)')
    parser.add_argument('-nmol', '--nmol', type=int, default=None,
                        help='Number of molecules for per-molecule quantities')
    parser.add_argument('-fluct_props', '--fluct_props', action='store_true',
                        help='Calculate fluctuation properties (isothermal compressibility, etc.)')

    return parser.parse_args()

def read_openmm_log(filename: str) -> pd.DataFrame:
    """Read OpenMM log file and return pandas DataFrame"""
    try:
        # Read the first line to get headers
        with open(filename, 'r') as f:
            first_line = f.readline().strip()
            if first_line.startswith('#'):
                first_line = first_line[1:]  # Remove the # character

        # Read the file, handling the tab-separated format
        df = pd.read_csv(filename, sep='\t', skiprows=1)

        # Set the column names from the first line
        column_names = [col.strip('"') for col in first_line.split('\t')]
        df.columns = column_names

        return df
    except FileNotFoundError:
        print(f"Error: Could not find file '{filename}'", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file '{filename}': {e}", file=sys.stderr)
        sys.exit(1)

def find_time_column(df: pd.DataFrame) -> Optional[str]:
    """Find the time column by examining column names"""
    columns = df.columns.tolist()

    # Common time column identifiers
    time_indicators = ['time', 'ps', 'picosecond', 'step', 'frame']

    # Look for exact matches or columns containing time indicators
    for col in columns:
        col_lower = col.lower()
        if any(indicator in col_lower for indicator in time_indicators):
            return col

    # If no time column found, return None
    return None

def filter_by_time(df: pd.DataFrame, begin: Optional[float] = None,
                   end: Optional[float] = None) -> pd.DataFrame:
    """Filter dataframe by time range"""
    time_col = find_time_column(df)

    if time_col is None:
        if begin is not None or end is not None:
            print("Warning: No time column found. Time filtering options (-b, -e) will be ignored.", file=sys.stderr)
        return df

    if begin is not None:
        df = df[df[time_col] >= begin]
    if end is not None:
        df = df[df[time_col] <= end]

    return df

def display_column_options(df: pd.DataFrame) -> None:
    """Display available columns for selection"""
    print("\nSelect the terms you want from the following list by", file=sys.stderr)
    print("selecting either (part of) the name or the number or a combination.", file=sys.stderr)
    print("End your selection with an empty line or a zero.", file=sys.stderr)
    print("-" * 71, file=sys.stderr)

    columns = df.columns.tolist()

    # Display columns in a nice format, 4 per row
    for i in range(0, len(columns), 4):
        row_cols = columns[i:i+4]
        row_str = ""
        for j, col in enumerate(row_cols):
            col_num = i + j + 1
            # Truncate long column names for display
            display_name = col[:15] + "..." if len(col) > 18 else col
            row_str += f"{col_num:3d}  {display_name:<18} "
        print(row_str, file=sys.stderr)

def get_user_selection(df: pd.DataFrame) -> List[int]:
    """Get user's column selection"""
    columns = df.columns.tolist()
    selected_indices = []

    while True:
        try:
            user_input = input().strip()

            # Empty input or 0 ends selection
            if not user_input or user_input == '0':
                break

            # Try to parse as number
            try:
                col_num = int(user_input)
                if 1 <= col_num <= len(columns):
                    if col_num - 1 not in selected_indices:
                        selected_indices.append(col_num - 1)
                else:
                    print(f"Invalid column number: {col_num}", file=sys.stderr)
                    continue
            except ValueError:
                # Try to match by name
                matches = []
                for i, col in enumerate(columns):
                    if user_input.lower() in col.lower():
                        matches.append(i)

                if len(matches) == 1:
                    if matches[0] not in selected_indices:
                        selected_indices.append(matches[0])
                elif len(matches) > 1:
                    print(f"Multiple matches found for '{user_input}':", file=sys.stderr)
                    for match_idx in matches:
                        print(f"  {match_idx + 1}: {columns[match_idx]}", file=sys.stderr)
                    continue
                else:
                    print(f"No matches found for '{user_input}'", file=sys.stderr)
                    continue

        except (EOFError, KeyboardInterrupt):
            print("\nExiting...", file=sys.stderr)
            sys.exit(0)

    return selected_indices

def calculate_statistics(data: np.ndarray) -> Tuple[float, float, float, float]:
    """Calculate mean, error estimate, RMSD, and drift"""
    mean_val = np.mean(data)

    # Simple error estimate (standard error)
    error_est = np.std(data) / np.sqrt(len(data))

    # RMSD (standard deviation)
    rmsd = np.std(data)

    # Total drift (difference between first and last values)
    drift = data[-1] - data[0] if len(data) > 1 else 0.0

    return mean_val, error_est, rmsd, drift

def get_unit_from_column_name(col_name: str) -> str:
    """Extract unit from column name"""
    if '(' in col_name and ')' in col_name:
        start = col_name.rfind('(')
        end = col_name.rfind(')')
        if start < end:
            return col_name[start+1:end]
    return ""

def clean_property_name_for_filename(col_name: str) -> str:
    """Clean property name to create a valid filename"""
    # Remove units in parentheses
    clean_name = re.sub(r'\s*\([^)]*\)', '', col_name)
    # Replace spaces with underscores
    clean_name = re.sub(r'\s+', '_', clean_name.strip())
    # Remove or replace other problematic characters
    clean_name = re.sub(r'[<>:"/\\|?*]', '', clean_name)
    # Remove leading/trailing underscores and convert to lowercase
    clean_name = clean_name.strip('_').lower()
    return clean_name

def generate_xvg_files(df: pd.DataFrame, selected_indices: List[int],
                      nmol: Optional[int] = None) -> None:
    """Generate .xvg files for each selected property"""
    columns = df.columns.tolist()
    time_col = find_time_column(df)

    # Use time data if available, otherwise use row indices starting from 1
    if time_col is not None:
        time_data = df[time_col].values
        x_label = "Time (ps)"
    else:
        time_data = np.arange(1, len(df) + 1)
        x_label = "Frame"

    for idx in selected_indices:
        col_name = columns[idx]
        data = df[col_name].values.copy()

        # Apply per-molecule scaling if requested
        if nmol is not None and 'Energy' in col_name:
            data = data / nmol

        # Create filename
        clean_name = clean_property_name_for_filename(col_name)
        filename = f"{clean_name}.xvg"

        # Create .xvg file content
        unit = get_unit_from_column_name(col_name)
        unit_str = f" ({unit})" if unit else ""

        with open(filename, 'w') as f:
            # Write GROMACS-style header
            f.write("# This file was created by OpenMM Energy Analysis Tool\n")
            f.write("# Generated from OpenMM simulation data\n")
            f.write(f"# Property: {col_name}\n")
            f.write("@    title \"OpenMM Energy Data\"\n")
            f.write(f"@    xaxis  label \"{x_label}\"\n")
            f.write(f"@    yaxis  label \"{col_name.split('(')[0].strip()}{unit_str}\"\n")
            f.write("@TYPE xy\n")

            # Write data
            for t, val in zip(time_data, data):
                f.write(f"{t:12.6f}  {val:12.6f}\n")

        print(f"Generated {filename}", file=sys.stderr)

def calculate_isothermal_compressibility(volume_data: np.ndarray,
                                       temperature: float) -> float:
    """
    Calculate isothermal compressibility using the formula:
    Kt = (<V^2> - <V>^2) / (kbT * <V>)

    Args:
        volume_data: Array of volume values in nm^3
        temperature: Temperature in K

    Returns:
        Isothermal compressibility in m^3/J
    """
    # Convert volume from nm^3 to m^3
    volume_m3 = volume_data * 1e-27  # 1 nm^3 = 1e-27 m^3

    # Calculate averages
    v_avg = np.mean(volume_m3)
    v2_avg = np.mean(volume_m3**2)

    # Boltzmann constant in J/K
    kb = 1.380649e-23

    # Calculate isothermal compressibility
    kt = (v2_avg - v_avg**2) / (kb * temperature * v_avg)

    return kt

def detect_fluctuation_properties(df: pd.DataFrame, selected_indices: List[int]) -> dict:
    """
    Detect which fluctuation properties can be calculated based on selected columns

    Returns:
        Dictionary with available fluctuation properties and their column indices
    """
    columns = df.columns.tolist()
    available_props = {}

    # Check for volume column (needed for isothermal compressibility)
    volume_idx = None
    for idx in selected_indices:
        col_name = columns[idx].lower()
        if 'volume' in col_name and 'nm^3' in col_name:
            volume_idx = idx
            break

    if volume_idx is not None:
        available_props['isothermal_compressibility'] = volume_idx

    # Check for temperature column (needed for compressibility calculation)
    temp_idx = None
    for idx in selected_indices:
        col_name = columns[idx].lower()
        if 'temperature' in col_name:
            temp_idx = idx
            break

    if temp_idx is not None:
        available_props['temperature'] = temp_idx

    return available_props

def print_fluctuation_properties(df: pd.DataFrame, selected_indices: List[int],
                                nmol: Optional[int] = None) -> None:
    """Print fluctuation properties calculations"""
    available_props = detect_fluctuation_properties(df, selected_indices)
    columns = df.columns.tolist()

    # Need both volume and temperature for isothermal compressibility
    if 'isothermal_compressibility' in available_props and 'temperature' in available_props:
        volume_idx = available_props['isothermal_compressibility']
        temp_idx = available_props['temperature']

        volume_data = df[columns[volume_idx]].values
        temp_data = df[columns[temp_idx]].values
        avg_temp = np.mean(temp_data)

        print(f"\nTemperature dependent fluctuation properties at T = {avg_temp:.3f}.")
        print()
        print("Heat capacities obtained from fluctuations do *not* include")
        print("quantum corrections. If you want to get a more accurate estimate")
        print("please use a quantum correction program.")
        print()
        print("WARNING: Please verify that your simulations are converged and perform")
        print("a block-averaging error analysis (not implemented in omm_energy yet)")

        # Calculate isothermal compressibility
        kt = calculate_isothermal_compressibility(volume_data, avg_temp)

        # Convert volume to m^3/mol using Avogadro's number
        avg_volume_nm3 = np.mean(volume_data)
        avg_volume_m3_per_mol = avg_volume_nm3 * 1e-27 * 6.02214076e23  # nm^3 to m^3/mol

        # Bulk modulus is reciprocal of compressibility
        bulk_modulus = 1.0 / kt if kt > 0 else 0.0

        print(f"Volume                                   = {avg_volume_m3_per_mol:.5e} m^3/mol")
        print(f"Isothermal Compressibility Kappa         = {kt:.4e} (m^3/J)")
        print(f"Adiabatic bulk modulus                   = {bulk_modulus:.5e} (J/m^3)")
        print()
    else:
        missing = []
        if 'isothermal_compressibility' not in available_props:
            missing.append("volume")
        if 'temperature' not in available_props:
            missing.append("temperature")

        print(f"\nCannot calculate fluctuation properties.", file=sys.stderr)
        print(f"Missing required properties: {', '.join(missing)}", file=sys.stderr)
        print("Please select volume and temperature columns to enable fluctuation property calculations.", file=sys.stderr)

def print_statistics(df: pd.DataFrame, selected_indices: List[int],
                    nmol: Optional[int] = None) -> None:
    """Print statistics for selected columns"""
    columns = df.columns.tolist()
    time_col = find_time_column(df)
    n_points = len(df)

    # Get time range info if time column exists
    if time_col is not None:
        time_data = df[time_col].values
        time_start = time_data[0]
        time_end = time_data[-1]
        print(f"\nLast energy frame read {n_points} time {time_end:8.3f}", file=sys.stderr)
        print(f"\nStatistics over {n_points} steps [ {time_start:8.4f} through {time_end:8.4f} ps ], {len(selected_indices)} data sets")
    else:
        print(f"\nRead {n_points} data frames", file=sys.stderr)
        print(f"\nStatistics over {n_points} data points, {len(selected_indices)} data sets")

    print(f"All statistics are over {n_points} points")
    print()
    print("Energy                      Average   Err.Est.       RMSD  Tot-Drift")
    print("-" * 79)

    for idx in selected_indices:
        col_name = columns[idx]
        data = df[col_name].values

        # Apply per-molecule scaling if requested
        if nmol is not None and 'Energy' in col_name:
            data = data / nmol

        mean_val, error_est, rmsd, drift = calculate_statistics(data)
        unit = get_unit_from_column_name(col_name)

        # Clean column name for display
        display_name = col_name.split('(')[0].strip()
        if len(display_name) > 27:
            display_name = display_name[:24] + "..."

        print(f"{display_name:<27} {mean_val:10.3f} {error_est:10.4g} {rmsd:10.4g} {drift:10.4g}  ({unit})")

def print_header():
    """Print program header information to stderr"""
    print("                :-) OpenMM Energy Analysis Tool (-:", file=sys.stderr)
    print("", file=sys.stderr)
    print("                          OpenMM Energy is written by:", file=sys.stderr)
    print("     Robert J. Weldon      Peter Eastman        Yutong Zhao       John D. Chodera", file=sys.stderr)
    print("    and the OpenMM development team", file=sys.stderr)
    print("", file=sys.stderr)
    print("               Based on the GROMACS gmx energy interface", file=sys.stderr)
    print("         Check out http://openmm.org for more information.", file=sys.stderr)
    print("", file=sys.stderr)
    print("OpenMM Energy:      omm_energy.py, version 1.0", file=sys.stderr)
    print("Working dir:        " + os.getcwd(), file=sys.stderr)
    print("Command line:", file=sys.stderr)
    print("  " + " ".join(sys.argv), file=sys.stderr)
    print("", file=sys.stderr)

def main():
    """Main function"""
    args = parse_arguments()

    # Print header information
    print_header()

    # Read the log file
    df = read_openmm_log(args.file)

    # Apply time filtering
    if args.begin is not None or args.end is not None:
        df = filter_by_time(df, args.begin, args.end)

        if len(df) == 0:
            print("Error: No data points in the specified time range", file=sys.stderr)
            sys.exit(1)

    # Display file info
    print(f"\nOpened {args.file} as OpenMM energy file", file=sys.stderr)
    print(file=sys.stderr)

    # Display column options
    display_column_options(df)

    # Get user selection
    selected_indices = get_user_selection(df)

    if not selected_indices:
        print("No columns selected. Exiting.", file=sys.stderr)
        sys.exit(0)

    # Print statistics
    print_statistics(df, selected_indices, args.nmol)

    # Print fluctuation properties if requested
    if args.fluct_props:
        print_fluctuation_properties(df, selected_indices, args.nmol)

    # Generate .xvg files for selected properties
    print(file=sys.stderr)
    generate_xvg_files(df, selected_indices, args.nmol)

    print("\nOMM Energy reminds you: \"Energize your simulations!\" (OpenMM)", file=sys.stderr)

if __name__ == "__main__":
    main()