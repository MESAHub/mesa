'''This python file was made exlusively for calulcating bi-cubic splines of LANL ATOMIC OPLIB opacity tables, not OPAL or OP (which contain 9.999 values in different table locations)'''
'''lasted updated by EbF 12-2023'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import pandas as pd
import os
import subprocess
import scipy as sp
import glob
import shutil
import sys
import time

start = time.time()

def np_array(shape):
    return np.array(shape, dtype=object)


def smooth_opac(fname):
    # importing files from history.data in each model folder
    file_name = fname.split('/')
    file_name1 = file_name[0]
    file_name2 = file_name[1]
    location = fname #file_name1 + ".data"
    #print ('calculating bi-cubic spline of table:',fname)
    #file_path1 = main_file + "/" + file_name1 + "/LOGS/history.data"
    model = pd.read_csv(location,skiprows=5,delim_whitespace=True)

    logt_vals = np.array(model.index)
    logr_vals = np.array(model.columns, dtype = float)
    rows = []
    for i in range (0,74):
        row_val = np.array(model.iloc[i,:])
        rows.append(row_val)
        
    rows = np.array(rows)


    '''now we remove the 9.9999e+0 values from opacity files'''
    for i in range (0,len(logr_vals)):
        interp_number = 0
        for j in range(0,len(logt_vals)):
            if ((rows[j,i] >= 9.99e0)  and (interp_number ==0)):
                interp_number = 1
                
                dx_to_interp = (logt_vals[j-1]- logt_vals[j-2])
                dy_to_interp =rows[j-1,i]- rows[j-2,i]
                rows[j,i] = (dy_to_interp/dx_to_interp)*(logt_vals[j]- logt_vals[j-interp_number]) + rows[j-interp_number,i]
            elif ((rows[j,i] >= 9.99e0) and (interp_number >=1)):
                interp_number = interp_number + 1 # add 1 to interpolation number
                rows[j,i] = (dy_to_interp/dx_to_interp)*(logt_vals[j]- logt_vals[j-interp_number]) + rows[j-interp_number,i]
        



    '''Generate a spline'''
    spline = sp.interpolate.RectBivariateSpline(logt_vals, logr_vals, rows, kx=3, bbox=[3.750, None, None, None], ky=3, s=0)
    # bbox=[None, None, None, None]

    '''spline dimensions'''
    #tarr1 = np.linspace(3.75, 8, 426)
    #tarr2 = np.linspace(8.025,9.05,42)
    #tarr = np.concatenate((tarr1,tarr2))
    tarr = np.linspace(3.75,9.05,213)
    rarr = np.linspace(-8, 1.5, 39)
    tgrid, rgrid = np.meshgrid(tarr, rarr, indexing="ij")

    kap_data_interp = spline(tgrid, rgrid, grid=False)



    '''output to new spline file'''
    
    input = open(location, "r")
    content = input.readlines()

    np.set_printoptions(linewidth=1000)

    output_name = output_file_dir + "/" + file_name2
    #print ('output spline name: ', output_name)
    with open(output_name, "w") as output:
        # write header
        for i in range (0,5):
            output.write(content[i])

        # write log R header
        index = np.array2string(rarr, formatter={'float_kind':lambda x: "{: 8.2f}".format(x)})
        index = index.replace(' [', '').replace('[', '').replace(']', '')
        #index = index.replace(' ', '     ')
        index = "       " + index
        output.write(index)
        output.write("\n")
        output.write("\n")
        # write logt + kap_vals
        for i in range(0,len(tarr)):
                kap_row_val = kap_data_interp[i,:]
                row = np.array2string(kap_row_val, formatter={'float_kind':lambda x: "{: 8.4f}".format(x)})
                row = row.replace(' [', '').replace('[', '').replace(']', '')
                row = f'{tarr[i]:.3f}' + '  ' + row
                #row = row.replace(' ', '  ')
                output.write(row)
                output.write("\n")

    return







'''Import opacity file'''
  
opacity_file_dir = sys.argv[1]
 #'MESA_OPLIB_GS98_standard_tables'#input('enter name of input opacity file directory: ')
output_file_dir = sys.argv[1] + "_spline"
#'spline_MESA_OPLIB_GS98_standard_tables' #input('enter name of output opacity file directory: ')

finame = []
for name in glob.glob(opacity_file_dir+"/*.data"):
    finame.append(name)
#print (finame)
# remove old output_file_dir if it exists.
shutil.rmtree(output_file_dir, ignore_errors=True)
# remake empty output file dir.
os.mkdir(output_file_dir)

# loop through all opacity file in input dir and spline.
for i in range(0,len(finame)):
    smooth_opac(finame[i])

# fix table limits and number of columns
sed_command = "sed -i 's/ 74    3.764000    9.065000/213    3.750000    9.050000/g' " + output_file_dir + "/*.data"
subprocess.run([sed_command], capture_output=True, shell=True, check=False)

# print final run time
end = time.time()
print(end - start)
