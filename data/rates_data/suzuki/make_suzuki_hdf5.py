#!/usr/bin/python

import io
import lzma
import re

import h5py
import numpy as np


def make_mesa_rxn_id(isos, wk_str):
    mesa_isos = []
    for iso in isos:
        m = re.match('(?P<A>[0-9]{1,3})(?P<Z>[A-Z][a-z]?)', iso)
        mesa_isos.append('{Z}{A}'.format(**m.groupdict()).lower())
    return '_'.join((mesa_isos[0], wk_str, mesa_isos[-1]))


class SuzukiTable:
    _columns = ('mu', 'dQ', 'Vs', 'rate', 'nu', 'gamma')
    _compression_opts = {'compression': 'gzip', }
    
    def __init__(self, description, data):
        
        self.description = description
        
        data1d = np.loadtxt(io.BytesIO(data))
        
        self.logRhoYs = np.unique(data1d[:, 0])
        self.logTs = np.unique(data1d[:, 1])
        
        nRhoYs = len(self.logRhoYs)
        nTs = len(self.logTs)
        
        self.data = {}
        for i, col in enumerate(self._columns):
            self.data[col] = data1d[:, i + 2].reshape((nRhoYs, nTs,))
    
    def add_1D_datasets(self, hdf5_group):
        
        hdf5_group.create_dataset("lYeRhos", data=self.logRhoYs,
                                  **self._compression_opts)
        hdf5_group.create_dataset("logTs", data=self.logTs,
                                  **self._compression_opts)
    
    def add_2D_datasets(self, hdf5_group):
        
        for k, v in self.data.items():
            hdf5_group.create_dataset(k, data=v, **self._compression_opts)


class SuzukiData:
    
    def __init__(self, filename):
        
        self.filename = filename
        
        self._read()
    
    def _read(self):
        
        # dictionary to store data tables
        self._tables = {}
        
        with lzma.open(self.filename, 'rt') as f:
            block = []
            description = None
            for line in f:
                if not line.strip():
                    pass
                elif line.startswith('*'):
                    pass
                elif line.strip().startswith('(MeV)'):
                    # lines like this ought to be commented, but aren't
                    pass
                elif line.startswith('!'):
                    if block:
                        data = ''.join(block).encode()
                        
                        print(f'..read table {rxn_id} ({process})')
                        t = SuzukiTable(description, data)
                        
                        if rxn_id in self._tables:
                            self._tables[rxn_id][process] = t
                        else:
                            self._tables[rxn_id] = {process: t}
                        
                        block = []
                        description = None
                    if description is None:
                        isos = re.findall('[0-9]{1,3}[A-Z][a-z]?', line)
                        if 'e-capture' in line:
                            process = 'capture'
                            wk_str = 'wk'
                        elif 'betap' in line:
                            process = 'decay'
                            wk_str = 'wk'
                        elif 'beta-decay' in line:
                            process = 'decay'
                            wk_str = 'wk-minus'
                        elif 'p(e+)-capture' in line:
                            process = 'capture'
                            wk_str = 'wk-minus'
                        rxn_id = make_mesa_rxn_id(isos, wk_str)
                        description = process
                else:
                    block.append(line)
        return
    
    def reactions(self):
        return self._tables.keys()
    
    def get_reaction_data(self, rxn):
        return self._tables[rxn]


def add_data(h5file, filename):
    # read suzuki datafile
    print(f'reading {filename}')
    data = SuzukiData(filename)
    
    for rxn_id in data.reactions():
        
        # make one group per MESA reaction
        print(f'..created sub-group {rxn_id}')
        rxn_grp = h5file.create_group(rxn_id)
        
        # get the data
        rxn_data = data.get_reaction_data(rxn_id)
        
        for process, table in rxn_data.items():
            
            # make one group per process
            print(f'....created sub-sub-group {process}')
            process_grp = rxn_grp.create_group(process)
            
            print(f'......added 2D datasets')
            table.add_2D_datasets(process_grp)
        
        else:
            
            print(f'....added 1D datasets')
            table.add_1D_datasets(rxn_grp)


suzuki_files = (
    'A17_FO_ScrExp.odat.xz',
    'A18_OFNe_ScrExp.odat.xz',
    'A19_OFNeNa_ScrExp.odat.xz',
    'A20_OFNeNaMg_ScrExp.odat.xz',
    'A21_OFNeNaMg_ScrExp.odat.xz',
    'A22_FNeNaMg_ScrExp.odat.xz',
    'A23_FNeNaMgAl_ScrExp.odat.xz',
    'A24_NeNaMgAlSi_ScrExp.odat.xz',
    'A25_NeNaMgAlSi_ScrExp.odat.xz',
    'A26_NaMgAlSi_ScrExp.odat.xz',
    'A27_NaMgAlSiP_ScrExp.odat.xz',
    'A28_NaMgAlSiPS_ScrExp.odat.xz',
)

if __name__ == '__main__':
    
    # open hdf5 file
    with h5py.File("Suzuki2016.h5", "w") as h5file:
        
        # add the data
        for filename in suzuki_files:
            add_data(h5file, filename)
