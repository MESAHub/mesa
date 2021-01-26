from yaml import load, Loader
from os import listdir
from os.path import isfile, join
from partial import Partial
from auto_diff_type import AutoDiffType
from make_auto_diff_type import make_auto_diff_type
from utils import py_to_fort, tab
from functions import *

# Get config files
config_path = '../config'
config_files = [f for f in listdir(config_path) if isfile(join(config_path, f)) and '.config' in f]
config_files = [join(config_path, f) for f in config_files]

compilation_list = []
use_list = []
compilation_list.append('support_functions.f90')
use_list.append(tab + 'use support_functions')

for f in config_files:
	print(f)
	if '_array' not in f:
		continue
	with open(f, 'r') as fi:
		data = load(fi, Loader=Loader)

		if data['array'] and data['fixed_length']:
			array_length = data['array_length']
		else:
			array_length = None

		partials = list(Partial(orders, data['array']) for orders in data['orders'])
		adr = AutoDiffType(data['name'], data['array'], array_length, partials)

		out_fi = open('../private/' + data['name'] + '_module.f90', 'w+')
		out_fi.write(py_to_fort(make_auto_diff_type(adr, unary_operators, binary_operators, comparison_operators, intrinsics)))
		out_fi.close()
		compilation_list.append(data['name'] + '.f90')
		use_list.append(tab + 'use ' + data['name'] + '_module')

compilation_list.append('auto_diff.f90')
fi = open('../make/source_list.txt', 'w+')
fi.write(' '.join(compilation_list))
fi.close()

fi = open('../public/auto_diff.f90', 'w+')
fi.write('module auto_diff\n')
fi.write('\n'.join(use_list))
fi.write('\n')
fi.write('end module auto_diff')
fi.close()
