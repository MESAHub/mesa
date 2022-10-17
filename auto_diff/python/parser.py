from os import listdir
from os.path import isfile, join

from yaml import load, Loader

from auto_diff_type import AutoDiffType
from functions import *
from make_auto_diff_type import make_auto_diff_type
from partial import Partial
from utils import py_to_fort, tab

# Get config files
config_path = '../config'
config_files = [f for f in listdir(config_path) if
                isfile(join(config_path, f)) and '.config' in f]
config_files = [join(config_path, f) for f in config_files]

# compilation_list stores a list of all the fortran files that will need
# compiling.
# This is used in the makefile.
compilation_list = []
compilation_list.append('support_functions.f90')

# use_list stores a list of all private auto_diff modules that need
# importing into the public auto_diff module.
use_list = []
use_list.append(tab + 'use support_functions')

# Loop over config files, building the relevant module for each.
for f in config_files:
    print(f)
    with open(f, 'r') as fi:
        data = load(fi, Loader=Loader)
        
        # gfortran does not (as of September 2021) support variable-length
        # arrays in parameterized-derived-types. So stick with fixed-length
        # arrays. If this changes in the future you can set fixed_length
        # to False and use variable-length arrays as desired.
        if data['array'] and data['fixed_length']:
            array_length = data['array_length']
        else:
            array_length = None
        
        # Read desired highest-order partial derivatives
        partials = list(
            Partial(orders, data['array']) for orders in data['orders'])
        
        # Build auto_diff type with those and all lower-order derivatives.
        adr = AutoDiffType(data['name'], data['array'], array_length, partials)
        
        out_fi = open('../private/' + data['name'] + '_module.f90', 'w+')
        out_fi.write(py_to_fort(
            make_auto_diff_type(adr, unary_operators, binary_operators,
                                comparison_operators, intrinsics)))
        out_fi.close()
        compilation_list.append(data['name'] + '.f90')
        use_list.append(tab + 'use ' + data['name'] + '_module')

# Write list of files to compile.
compilation_list.append('auto_diff.f90')
fi = open('../make/source_list.txt', 'w+')
fi.write(' '.join(compilation_list))
fi.close()

# Write public module that just imports all the private ones.
fi = open('../public/auto_diff.f90', 'w+')
fi.write('module auto_diff\n')
fi.write('\n'.join(use_list))
fi.write('\n')
fi.write('end module auto_diff')
fi.close()
