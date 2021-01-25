from sympy import symbols, factorial, diff, DiracDelta, simplify
from utils import wrap_element, fortran_substitutions, ddelta, zero_function
from sympy.simplify.cse_main import cse
from sympy.utilities.iterables import numbered_symbols

def unary_generic_chain_rule(auto_diff_type, fixed_length=None):
	'''
	Produces a list of lines of Fortran code to compute all specified partial derivatives of a generic unary function z(x),
	where x is itself a function of independent variables (val1, val2, val3, ...).
	
	This is done by constructing a series representation of z out to the desired order,
	then passing that to the unary_specific_chain_rule method.
	The result is given in terms of the derivatives of z which appear in the series expansion.
	'''

	# Construct the symbols corresponding to the partial derivatives of z.
	# These are d1z, d2z, ..., dNz, giving dz/dx, d^2 / dx^2, and so on.
	z_symbol_strs = ['z_' + str(p).replace('val1','x') for p in auto_diff_type.unary_partials]
	z_symbol_str = ' '.join(z_symbol_strs)
	z_syms = wrap_element(symbols(z_symbol_str, real=True))

	def operator(x):
		# Construct z as a power series in terms of its partial derivatives (z_syms) with
		# respect to the x.
		z = sum(sym * x**p.orders[0] / factorial(p.orders[0]) for sym,p in zip(*(z_syms, auto_diff_type.unary_partials)))
		return z

	return unary_specific_chain_rule(auto_diff_type, operator, xval=0, fixed_length=fixed_length)


def binary_generic_chain_rule(auto_diff_type, fixed_length=None):
	'''
	Produces a list of lines of Fortran code to compute all specified partial derivatives of a generic unary function z(x),
	where x is itself a function of independent variables (val1, val2, val3, ...).
	
	This is done by constructing a series representation of z out to the desired order,
	then passing that to the unary_specific_chain_rule method.
	The result is given in terms of the derivatives of z which appear in the series expansion.
	'''

	# Construct the symbols corresponding to the partial derivatives of z.
	# These are d1z, d2z, ..., dNz, giving dz/dx, d^2 / dx^2, and so on.
	z_symbol_strs = ['z_' + str(p).replace('val1','x').replace('val2','y') for p in auto_diff_type.binary_partials]
	z_symbol_str = ' '.join(z_symbol_strs)
	z_syms = wrap_element(symbols(z_symbol_str, real=True))

	def operator(x, y):
		# Construct z as a power series in terms of its partial derivatives (z_syms) with
		# respect to the x and y.
		z = sum(sym * x**p.orders[0] * y**p.orders[1] / (factorial(p.orders[0]) * factorial(p.orders[1])) for sym,p in zip(*(z_syms, auto_diff_type.binary_partials)))
		return z

	return binary_specific_chain_rule(auto_diff_type, operator, xval=0, yval=0, fixed_length=fixed_length)


def unary_specific_chain_rule(auto_diff_type, operator, xval=None, fixed_length=None):
	'''
	Produces a list of lines of Fortran code to compute all specified partial derivatives of a specified unary function z(x),
	where x is itself a function of independent variables (val1, val2, val3, ...) and z(x) === operator(x).

	This is done via the power series method.
		1. x is first constructed as a power series in terms of the desired partial derivatives.
		2. Then z is constructed as a function of x.
		3. Finally, we use sympy to take symbolic derivatives and extract the chain rule result.

	The result is evaluated at the specified xval, if given. Otherwise the value of x is left symbolic.
	'''

	partials = auto_diff_type.partials

	# Construct sympy variables corresponding to the various independent variables.
	# These never appear on the Fortran side, but we keep the naming consistent to correspond to the
	# names in partial_orders.
	# So these are called val1, val2, ..., valN.
	indep_symbol_str = ' '.join(auto_diff_type.partials[0].val_name(i) for i in range(auto_diff_type.num_independent))
	indep_syms = wrap_element(symbols(indep_symbol_str, real=True))

	# Construct sympy variables corresponding to the various derivatives dx/d(...).
	# Note that these variable names correspond to the names we'll use on the Fortran side, so
	# we can just directly map sympy expressions to strings and get valid Fortran :-)
	# Hence these are called x%d1val1, x%d2val1, ..., x%d1val2, x%d2val2, ..., x%d1val1_d1val2, ...
	# The first integer in each 'd_val_' block is the number of derivatives,
	# the second is the independent variable those derivatives are with respect to.
	x_symbol_str = ' '.join(auto_diff_type.partial_str_in_instance('x', p).replace(':','colon') for p in partials)
	x_syms = wrap_element(symbols(x_symbol_str, real=True))

	# Construct x as a power series in terms of its partial derivatives (sym) with respect to the independent
	# variables (indep).
	x = 0
	for p,sym in zip(*(partials, x_syms)):
		term = sym
		for order, indep in zip(*(p.orders, indep_syms)):
			term = term * indep ** order / factorial(order)
		x = x + term

	z = operator(x)

	# Extract chain rule expressions
	expressions = []
	left_hand_names = []
	derivatives = []
	for p in partials:
		# Beginning of Fortran line, for writing out the answer
		unary_symbol_str = auto_diff_type.partial_str_in_instance('unary', p).replace(':','colon')
		
		# Construct partial derivative specified by orders
		# of the power series, evaluated at (val1, val2, ...) = 0.
		d = z
		for order, indep in zip(*(p.orders, indep_syms)):
			d = diff(d, indep, order)
			d = d.replace(DiracDelta, zero_function) # Diract Delta is only non-zero in a set of measure zero, so we set it to zero.
			d = d.subs(indep, 0)

		if xval is not None:
			d = d.subs(x_syms[0], 0)

		d = simplify(d)

		expressions.append(d)
		left_hand_names.append(unary_symbol_str)

	# Optimize expressions by eliminating common sub-expressions
	intermediate_vars_to_declare = []
	common_sub_expressions = cse(expressions, symbols=numbered_symbols('q'))
	expressions = common_sub_expressions[1] # After substitution

	# Determine which sub-expressions are array-valued.
	is_array = {}
	for sym, sub_expr in common_sub_expressions[0]:
		is_array[sym] = False

	for i in range(len(common_sub_expressions[0])):
		for sym, sub_expr in common_sub_expressions[0]:
			if 'colon' in str(sub_expr): # Intermediate is an array type.
				is_array[sym] = True
			else:
				for s in sub_expr.args:
					if s in is_array:
						if is_array[s]:
							is_array[sym] = True # Intermediate is defined in terms of an array type and so is an array type.

	# Write the lines defining the sub-expressions
	for sym, sub_expr in common_sub_expressions[0][::-1]: # Have to go in reverse order because we're prepending.
		if is_array[sym]: # # Intermediate is an array type.
			if fixed_length is None:
				sym_str = str(sym) + '(1:x%n)' # All of these expressions have x defined as an auto_diff type, so it has an %n property.
			else:
				sym_str = str(sym) + '(1:' + str(fixed_length) + ')'
		else:
			sym_str = str(sym)

		left_hand_names.insert(0, sym_str)
		expressions.insert(0, sub_expr.subs(sym, symbols(sym_str)))
		intermediate_vars_to_declare.append(sym_str)

	for d,unary_symbol_str in zip(*(expressions, left_hand_names)):
		derivatives.append(unary_symbol_str + ' = ' + str(fortran_substitutions(d)))

	declarations = []
	for var in intermediate_vars_to_declare:
		declarations.append('real(dp) :: ' + var)

	return derivatives, declarations


def binary_specific_chain_rule(auto_diff_type, operator, xval=None, yval=None, fixed_length=None):
	'''
	Produces a list of lines of Fortran code to compute all specified partial derivatives of a specific
	binary function z(x,y), where x and y are functions of independent variables (val1, val2, val3, ...)
	and z === operator(x,y).

	This is done via the power series method.
		1. x and y are first constructed as a power series in terms of the desired partial derivatives.
		2. Then z is constructed as a power series in x and y in terms of the same partial derivatives.
		3. Finally, we use sympy to take symbolic derivatives and extract the chain rule result.

	The result is evaluated at the specified xval and yval, if given. Otherwise x and y are left symbolic.
	'''

	partials = auto_diff_type.partials

	# Construct sympy variables corresponding to the various independent variables.
	# These never appear on the Fortran side, but we keep the naming consistent to correspond to the
	# names in partial_orders.
	# So these are called val1, val2, ..., valN.
	indep_symbol_str = ' '.join(auto_diff_type.partials[0].val_name(i) for i in range(auto_diff_type.num_independent))
	indep_syms = wrap_element(symbols(indep_symbol_str, real=True))

	# Construct sympy variables corresponding to the various derivatives dx/d(...).
	# Note that these variable names correspond to the names we'll use on the Fortran side, so
	# we can just directly map sympy expressions to strings and get valid Fortran :-)
	# Hence these are called x%d1val1, x%d2val1, ..., x%d1val2, x%d2val2, ..., x%d1val1_d1val2, ...
	# The first integer in each 'd_val_' block is the number of derivatives,
	# the second is the independent variable those derivatives are with respect to.
	# Note that unlike in the generic form we include an offset which is the value at which the derivative is evaluated.
	x_symbol_str = ' '.join(auto_diff_type.partial_str_in_instance('x', p).replace(':','colon') for p in partials)
	x_syms = wrap_element(symbols(x_symbol_str, real=True))
	y_symbol_str = ' '.join(auto_diff_type.partial_str_in_instance('y', p).replace(':','colon') for p in partials)
	y_syms = wrap_element(symbols(y_symbol_str, real=True))

	# Construct x and y as power series in terms of their partial derivatives (sym)
	# with respect to the independent variables (indep).
	x = 0
	for p,sym in zip(*(partials, x_syms)):
		term = sym
		for order, indep in zip(*(p.orders, indep_syms)):
			term = term * indep ** order / factorial(order)
		x = x + term

	y = 0
	for p,sym in zip(*(partials, y_syms)):
		term = sym
		for order, indep in zip(*(p.orders, indep_syms)):
			term = term * indep ** order / factorial(order)
		y = y + term

	z = operator(x,y)

	# Extract chain rule expressions
	expressions = []
	left_hand_names = []
	derivatives = []
	for p in partials:
		# Beginning of Fortran line, for writing out the answer
		binary_symbol_str = auto_diff_type.partial_str_in_instance('binary', p).replace(':','colon')
		
		# Construct partial derivative specified by orders
		# of the power series, evaluated at (val1, val2, ...) = 0.
		d = z

		for order, indep in zip(*(p.orders, indep_syms)):
			d = diff(d, indep, order)
			d = d.replace(DiracDelta, zero_function) # Diract Delta is only non-zero in a set of measure zero, so we set it to zero.
			d = d.subs(indep, 0)

		if xval is not None:
			d = d.subs(x_syms[0], xval)
		if yval is not None:
			d = d.subs(y_syms[0], yval)

		expressions.append(d)
		left_hand_names.append(binary_symbol_str)

	# Optimize expressions by eliminating common sub-expressions
	intermediate_vars_to_declare = []
	common_sub_expressions = cse(expressions, symbols=numbered_symbols('q'))
	expressions = common_sub_expressions[1] # After substitution

	# Determine which sub-expressions are array-valued.
	is_array = {}
	for sym, sub_expr in common_sub_expressions[0]:
		is_array[sym] = False

	for i in range(len(common_sub_expressions[0])):
		for sym, sub_expr in common_sub_expressions[0]:
			if str(sym) == 'q34':
				print(sym,sub_expr)
			if 'colon' in str(sub_expr): # Intermediate is an array type.
				is_array[sym] = True
			else:
				for s in sub_expr.free_symbols:
					if s in is_array:
						if str(sym) == 'q34':
							print(sym, s, is_array[s])
						if is_array[s]:
							is_array[sym] = True # Intermediate is defined in terms of an array type and so is an array type.

	# Write the lines defining the sub-expressions
	for sym, sub_expr in common_sub_expressions[0][::-1]: # Have to go in reverse order because we're prepending.
		if is_array[sym]: # # Intermediate is an array type.
			if fixed_length is None:
				sym_str = str(sym) + '(1:x%n)' # All of these expressions have x defined as an auto_diff type, so it has an %n property.
			else:
				sym_str = str(sym) + '(1:' + str(fixed_length) + ')'
		else:
			sym_str = str(sym)

		left_hand_names.insert(0, sym_str)
		expressions.insert(0, sub_expr.subs(sym, symbols(sym_str)))
		intermediate_vars_to_declare.append(sym_str)

	for d,unary_symbol_str in zip(*(expressions, left_hand_names)):
		derivatives.append(unary_symbol_str + ' = ' + str(fortran_substitutions(d)))

	declarations = []
	for var in intermediate_vars_to_declare:
		declarations.append('real(dp) :: ' + var)

	return derivatives, declarations
