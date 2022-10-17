#!/usr/bin/env python

import io
import os
import sys

import h5py
import numpy as np
import yaml

HDF5_OPTS = {'dtype': 'float32', 'compression': 'gzip'}
"""dict: options used when creating non-scalar HDF5 datasets

default is store as half-precision and compress with gzip
"""


class AesopusTable:
    """Storage for a single AESOPUS table

    Attributes:
        params (dict): composition parameters of this table
        data (array): tabulated opacities (log10); [logTs, logRs]
    """
    
    def __init__(self, description, table):
        """Create table from relevant parts of AESOPUS output

        Args:
            description (str): line beginning # TABLE
            tbl (str): uncommented lines
        """
        
        self._description = description
        self._parse_description()
        self._read_table(table)
    
    def _parse_description(self):
        """Extract composition parameters"""
        
        # split and throw away # TABLE N part
        items = self._description.split()[3:]
        
        # extract parameter names, removing =
        keys = [k.strip('=') for k in items[0::2]]
        
        # extract parameter values, converting to floats
        values = []
        for v in items[1::2]:
            try:
                fv = float(v)
            except ValueError:
                if v == 'NEGATIVE':
                    fv = -1
            finally:
                values.append(fv)
        
        self.params = dict(zip(keys, values))
    
    def _read_table(self, table):
        """Extract opacities into array"""
        
        stream = io.StringIO(''.join(table))
        data = np.loadtxt(stream, comments='#', dtype='float32')
        
        # data has temperatures in decreasing order; reverse
        self.data = np.flipud(data[:, 1:])


def read_AESOPUS_tables(filename, nT):
    """Extract individual tables from AESOPUS output

    Args:
        filename (str): name of AESOPUS data file
        nTs (int): number of AESOPUS temperatures

    Returns:
        list of AesopusTable objects
    """
    
    # first, we read the entire file
    with open(filename) as f:
        data = f.readlines()
    
    # then, loop through and identify individual tables
    tables = {}
    summary_line = ''
    for i, line in enumerate(data):
        if "TABLE" in line:
            # the first time we hit TABLE, it is in a overall summary like
            # TOTAL NUMBER OF TABLES:   mh x mco x mc x mn =   100
            # each time thereafter it marks the start of a new table
            # and contains information about the composition parameters
            # TABLE      1     X=  0.500000     Y=  0.480000 ...
            if not summary_line:
                summary_line = line
            else:
                # save the line numbers and descriptions
                tables[i] = line
    
    # now, go through the tables that were found and extract each one
    extracted_tables = []
    for l, description in tables.items():
        # there are 103 + nT lines between the TABLE line and the end
        nl = 103 + nT
        table = data[l:l + nl]
        t = AesopusTable(description, table)
        extracted_tables.append(t)
    
    return extracted_tables


def write_AESOPUS_tables(tables, h5group, nT, nR):
    """Write a set of AesopusTable into an HDF5 group

    Args:
        filename (list of AesopusTable): name of AESOPUS data file
        h5group (h5.Group): HDF5 group to store data

    Returns:
        None
    """
    
    # extract unique parameters
    def unique_params(name, tables):
        vals = sorted(set([t.params[name] for t in tables]))
        return np.array(vals, dtype='float32')
    
    # make dataset
    def save_dataset(name, array):
        h5group.create_dataset(name, data=array, **HDF5_OPTS)
    
    # get the unique compositions parameters
    # it is assumed there will be a table for combination Xs x fCOs x fCs x fNs
    AESOPUS_Xs = unique_params('X', tables)
    AESOPUS_fCOs = unique_params('FCO', tables)
    AESOPUS_fCs = unique_params('FC', tables)
    AESOPUS_fNs = unique_params('FN', tables)
    
    # save compositions to file
    save_dataset("Xs", AESOPUS_Xs)
    save_dataset("fCOs", AESOPUS_fCOs)
    save_dataset("fCs", AESOPUS_fCs)
    save_dataset("fNs", AESOPUS_fNs)
    
    # make a big dataset that can hold all opacity tables
    nX = AESOPUS_Xs.size
    nCO = AESOPUS_fCOs.size
    nC = AESOPUS_fCs.size
    nN = AESOPUS_fNs.size
    
    dset = h5group.create_dataset("kap", (nT, nR, nX, nCO, nC, nN),
                                  **HDF5_OPTS)
    
    # put the tables in an iterator
    # loop through the parameters in the same order used in the AESOPUS file
    # Xs vary most slowly; fNs most rapidly
    ts = iter(tables)
    for iX, fX in enumerate(AESOPUS_Xs):
        for iCO, fCO in enumerate(AESOPUS_fCOs):
            for iC, fC in enumerate(AESOPUS_fCs):
                for iN, fN in enumerate(AESOPUS_fNs):
                    t = next(ts)
                    print(t._description.strip())
                    data = t.data
                    dset[:, :, iX, iCO, iC, iN] = data


def main(config):
    # one must manually extract these values from the AESOPUS file
    # they depend on the solar abundance pattern
    Zsun = config['Zsun']
    C_div_Z_sun = config['C_div_Z_sun']
    N_div_Z_sun = config['N_div_Z_sun']
    O_div_Z_sun = config['O_div_Z_sun']
    
    # the AESOPUS web form returns a file for each Z
    # one must manually provide the list of files and their base metallicities
    # these should be provided in increasing order
    AESOPUS_Zs = np.array(config['Zs'], dtype='float32')
    
    AESOPUS_logRs = np.array(config['logRs'], dtype='float32')
    AESOPUS_logTs = np.array(config['logTs'], dtype='float32')
    
    # open HDF5 file
    with h5py.File(config['output'], "w") as h5file:
        
        # calculate and store reference composition values
        fCO_ref = np.log10(C_div_Z_sun / O_div_Z_sun)
        fC_ref = np.log10(C_div_Z_sun)
        fN_ref = np.log10(N_div_Z_sun)
        
        h5file.create_dataset("Zsun", data=Zsun)
        h5file.create_dataset("fCO_ref", data=fCO_ref)
        h5file.create_dataset("fC_ref", data=fC_ref)
        h5file.create_dataset("fN_ref", data=fN_ref)
        
        # store logTs and logRs
        h5file.create_dataset("logTs", data=AESOPUS_logTs, **HDF5_OPTS)
        h5file.create_dataset("logRs", data=AESOPUS_logRs, **HDF5_OPTS)
        
        # store Zs
        h5file.create_dataset("Zs", data=AESOPUS_Zs, **HDF5_OPTS)
        
        # now, parse and store the tables associated with each Z
        for Z, table, in zip(AESOPUS_Zs, config['files']):
            
            # the Z_id string must be such that a lexicographic sort
            # gives the same order as a numerical sort in Z
            Z_id = f'{Z:8.6f}'
            Z_grp = h5file.create_group(Z_id)
            
            filename = os.path.join(config['directory'], table)
            tables = read_AESOPUS_tables(filename, AESOPUS_logTs.size)
            write_AESOPUS_tables(tables, Z_grp,
                                 AESOPUS_logTs.size, AESOPUS_logRs.size)


if __name__ == '__main__':
    
    with open(sys.argv[1]) as f:
        y = yaml.safe_load(f.read())
    
    main(y)
