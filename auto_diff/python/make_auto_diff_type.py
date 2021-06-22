from utils import tab, indent_list, indent_list_of_str

def make_auto_diff_type(auto_diff_type, unary_operators, binary_operators, comparison_operators, intrinsics):
	# Starting boilerplate
	header = ['module ' + auto_diff_type.name + '_module']
	begin = [tab + 'use const_def, only: dp, ln10, pi', tab + 'use utils_lib', tab + 'use support_functions', tab + 'use math_lib']
	begin = begin + ['', tab + 'implicit none', '', tab + 'private', '']
	contains = [tab + 'contains', '']

	# Export interfaces
	begin.append('public :: ' + auto_diff_type.name + ', &')

	interfaces = []
	functions = []

	# Type
	type_lines = auto_diff_type.type() + ['']

	# Assignment
	interfaces.append('interface assignment(=)')
	begin.append(tab + 'assignment(=), &')

	assignment = auto_diff_type.assignment_routine_from_self()
	functions.append(str(assignment))
	functions.append('')
	interfaces.append(tab + 'module procedure ' + assignment.name)

	assignment = auto_diff_type.assignment_routine_from_real_dp()
	functions.append(str(assignment))
	functions.append('')
	interfaces.append(tab + 'module procedure ' + assignment.name)

	assignment = auto_diff_type.assignment_routine_from_int()
	functions.append(str(assignment))
	functions.append('')
	interfaces.append(tab + 'module procedure ' + assignment.name)

	interfaces.append('end interface assignment(=)')
	interfaces.append('')

	# Comparison
	for op, opname in comparison_operators:
		interfaces.append('interface operator(' + op + ')')
		begin.append(tab + 'operator(' + op + '), &')

		# (auto_diff, auto_diff)
		self_comparison = auto_diff_type.comparison_function_with_self(opname, op)
		functions.append(str(self_comparison))
		functions.append('')
		interfaces.append(tab + 'module procedure ' + self_comparison.name)

		# (auto_diff, real(dp)) and vice-versa
		real_comparisons = auto_diff_type.comparison_function_with_real_dp(opname, op)
		functions.append(str(real_comparisons[0]))
		functions.append('')
		functions.append(str(real_comparisons[1]))
		functions.append('')
		interfaces.append(tab + 'module procedure ' + real_comparisons[0].name)
		interfaces.append(tab + 'module procedure ' + real_comparisons[1].name)

		# (auto_diff, integer) and vice-versa
		int_comparisons = auto_diff_type.comparison_function_with_int(opname, op)
		functions.append(str(int_comparisons[0]))
		functions.append('')
		functions.append(str(int_comparisons[1]))
		functions.append('')
		interfaces.append(tab + 'module procedure ' + int_comparisons[0].name)
		interfaces.append(tab + 'module procedure ' + int_comparisons[1].name)

		interfaces.append('end interface operator(' + op + ')')
		interfaces.append('')

	# Generic Operators
	unary_operator = auto_diff_type.generic_unary_operator_function()
	functions.append(str(unary_operator))
	functions.append('')

	begin.append(tab + 'make_unop, &')
	interfaces.append('interface make_unop')
	interfaces.append(tab + 'module procedure ' + unary_operator.name)
	interfaces.append('end interface make_unop')
	interfaces.append('')

	binary_operator = auto_diff_type.generic_binary_operator_function()
	functions.append(str(binary_operator))
	functions.append('')

	begin.append(tab + 'make_binop, &')
	interfaces.append('interface make_binop')
	interfaces.append(tab + 'module procedure ' + binary_operator.name)
	interfaces.append('end interface make_binop')
	interfaces.append('')

	# Specific Operators
	for op, opname in unary_operators:

		function = auto_diff_type.specific_unary_operator_function(opname, op)
		functions.append(str(function))
		functions.append('')

		if opname in intrinsics:
			interopname = 'operator(' + intrinsics[opname] + ')'
		else:
			interopname = opname

		begin.append(tab + interopname + ', &')
		interfaces.append('interface ' + interopname)
		interfaces.append(tab + 'module procedure ' + function.name)
		interfaces.append('end interface ' + interopname)
		interfaces.append('')

	for op, opname in binary_operators:
		# (auto_diff, auto_diff)
		function = auto_diff_type.specific_binary_operator_function_self(opname, op)
		functions.append(str(function))
		functions.append('')

		if opname in intrinsics:
			interopname = 'operator(' + intrinsics[opname] + ')'
		else:
			interopname = opname

		begin.append(tab + interopname + ', &')
		interfaces.append('interface ' + interopname)
		interfaces.append(tab + 'module procedure ' + function.name)

		# (auto_diff, real(dp)) and vice-versa
		f1, f2 = auto_diff_type.specific_binary_operator_function_real_dp(opname, op)
		functions.append(str(f1))
		functions.append('')
		functions.append(str(f2))
		functions.append('')

		interfaces.append(tab + 'module procedure ' + f1.name)
		interfaces.append(tab + 'module procedure ' + f2.name)

		# (auto_diff, integer) and vice-versa
		f1, f2 = auto_diff_type.specific_binary_operator_function_int(opname, op)
		functions.append(str(f1))
		functions.append('')
		functions.append(str(f2))
		functions.append('')

		interfaces.append(tab + 'module procedure ' + f1.name)
		interfaces.append(tab + 'module procedure ' + f2.name)
		interfaces.append('end interface ' + interopname)
		interfaces.append('')

	# Derivative functions
	if auto_diff_type.array:
		num_independent_non_array = auto_diff_type.num_independent - 1
	else:
		num_independent_non_array = auto_diff_type.num_independent

	for i in range(num_independent_non_array):
		derivative_function = auto_diff_type.derivative_function(i)
		functions.append(str(derivative_function))
		functions.append('')

		begin.append(tab + 'differentiate_' + str(i+1) + ', &')
		interfaces.append('interface differentiate_' + str(i+1))
		interfaces.append(tab + 'module procedure ' + derivative_function.name)
		interfaces.append('end interface differentiate_' + str(i+1))
		interfaces.append('')

	# Clean up begin
	begin[-1] = begin[-1][:-3] # Remove ending ', &'
	begin.append('')

	# De-duplicate begin. This is necessary because the unary_minus operator and subtraction have the same interface, but get added twice.
	begin_de_duplicated = []
	for b in begin:
		if b not in begin_de_duplicated:
			begin_de_duplicated.append(b)
	begin = begin_de_duplicated
	# Ending boilerplate
	ender = ['end module ' + auto_diff_type.name + '_module']

	# Put it all together
	lines = header + indent_list(begin + type_lines + interfaces) + contains + indent_list_of_str(functions) + ender

	return '\n'.join(lines)
