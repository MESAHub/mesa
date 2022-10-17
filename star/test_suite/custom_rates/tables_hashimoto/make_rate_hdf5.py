#!/usr/bin/python

import h5py
import numpy as np


class GMPTable:
    _columns = ('uf', 'lbeta+', 'leps-', 'lrnu', 'lbeta-', 'leps+', 'lrnubar')
    
    def __init__(self, filename):
        
        data1d = np.loadtxt(filename, skiprows=4)
        
        self.T9s = np.unique(data1d[:, 0])
        self.lYeRhos = np.unique(data1d[:, 1])
        
        nYeRhos = len(self.lYeRhos)
        nTs = len(self.T9s)
        
        self.data = {}
        for i, col in enumerate(self._columns):
            self.data[col] = data1d[:, i + 2].reshape((nYeRhos, nTs))


if __name__ == '__main__':
    
    filename = "n14c-hr.dat"
    table = GMPTable(filename)
    
    # open hdf5 file
    with h5py.File("GMP_r_n14_wk_c14.h5", "w") as h5file:
        
        # the 1D data sets
        h5file.create_dataset("lYeRhos", data=table.lYeRhos)
        h5file.create_dataset("T9s", data=table.T9s)
        
        # the 2D data sets
        h5file.create_dataset("ldecay", data=table.data['lbeta+'])
        h5file.create_dataset("lcapture", data=table.data['leps-'])
        h5file.create_dataset("lneutrino", data=table.data['lrnu'])
    
    # open hdf5 file
    with h5py.File("GMP_r_c14_wk-minus_n14.h5", "w") as h5file:
        
        # the 1D data sets
        h5file.create_dataset("lYeRhos", data=table.lYeRhos)
        h5file.create_dataset("T9s", data=table.T9s)
        
        # the 2D data sets
        h5file.create_dataset("ldecay", data=table.data['lbeta-'])
        h5file.create_dataset("lcapture", data=table.data['leps+'])
        h5file.create_dataset("lneutrino", data=table.data['lrnubar'])
