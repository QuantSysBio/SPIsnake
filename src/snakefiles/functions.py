""" Functions for writing SPIsnake input mass/rt lists.
"""
import glob
import json
import sys

import pandas as pd
from pyteomics import mgf

PROTON = 1.007276466622

def extract_mgf_spectrum_data(spectrum, input_filename):
    """ Function to extract useful information from a single scan in an mgf file.
    """
    # Dictionary for all data.
    spec_data = {}

    # Add source file and scan number identifiers.
    spec_data['source'] = input_filename.split('/')[-1].split('.mgf')[0]
    spec_data['scan'] = int(spectrum['params']['title'].split('scan=')[-1].strip('"'))

    # Calculate precursor mass from m/z.
    pep_mz = float(spectrum['params']['pepmass'][0])
    charge = int(spectrum['params']['charge'][0])
    pep_mass = (pep_mz*charge) - (PROTON*charge)
    spec_data['precursorMass'] = pep_mass

    # Add other important scan results.
    spec_data['charge'] = charge
    spec_data['retentionTime'] = float(spectrum['params']['rtinseconds'])
    spec_data['fragmentMZs'] = json.dumps(list(spectrum['m/z array']))

    return spec_data

def generate_mass_rt_list_single_file(input_filename):
    # Fill list of scan data.
    collected_scan_data = []
    with mgf.read(input_filename) as reader:
        for spectrum in reader:
            spectrum_data = extract_mgf_spectrum_data(spectrum, input_filename)
            collected_scan_data.append(spectrum_data)


    # Create DataFrame out of the scan data.
    scan_df = pd.DataFrame(collected_scan_data)

    return scan_df
