from sympy import symbols, Function, preorder_traversal, Integer, Rational, Pow
from sympy import diff as diff_sym
from sympy import DiracDelta, Pow, Integer, symbols, Rational, Float, srepr
import re
from collections.abc import Iterable

tab = '   '
x,y,z = symbols('x y z', real=True)
ppow = Function('PPow')
ddelta = Function('ddelta')
ssqrt = Function('ssqrt')
powm1 = Function('powm1')
powN = list(Function('pow' + str(i)) for i in range(9))

def zero_function(*args):
	return 0

def substitute_powm1(expr):
	# Substitutes pow(x, -i) for Div(1, x) to avoid log calls.

	# Search for Pow
	to_sub = [p for p in preorder_traversal(expr) if p.func == Pow]

	# Search for -1 arguments
	to_sub = [p for p in to_sub if p.args[1] == Integer(-1)]

	# Build substitutions
	new_ex = [() for p in to_sub]

def fortran_substitutions(deriv):
	# Replace calls to half-integer powers with calls to sqrt**n
	deriv = substitute_pow(deriv)

	# Next replace rational numbers with floating-point approximations
	# to avoid Fortran integer division and issues in pow calls.
	deriv = substitute_rational(deriv)

	# Put in _dp's for the remaining floats.
	deriv = substitute_dp(deriv)

	# Now we replace ** with PPow so we can later do a string replacement to get just the fortran pow
	deriv = deriv.replace(Pow, ppow)

	return deriv

def substitute_rational(expr):
	# Substitutes i/j for float(i)/float(j) for all integers i,j.
	# This avoids Fortran integer division from ruining our expressions through rounding.

	# Search for Rational
	to_sub = [p for p in preorder_traversal(expr) if isinstance(p, Rational)]

	# Build substitutions
	new_ex = [symbols(str(float(p)) + '_dp') for p in to_sub]

	# Substitute. We use xreplace because this does a direct replacement
	# of only exactly matching subexpressions, rather than what subs does
	# (reason to try to shoehorn our replacement in) or what replace does
	# (wildcard matching).
	for (old,new) in zip(to_sub, new_ex):
		expr = expr.xreplace({old:new})

	return expr

def substitute_dp(expr):
	# Substitutes x_dp for x whenever x is a type Float.

	# Search for Float
	to_sub = [p for p in preorder_traversal(expr) if isinstance(p, Float)]

	# Build substitutions
	new_ex = [symbols(str(float(p)) + '_dp') for p in to_sub]

	# Substitute. We use xreplace because this does a direct replacement
	# of only exactly matching subexpressions, rather than what subs does
	# (reason to try to shoehorn our replacement in) or what replace does
	# (wildcard matching).
	for (old,new) in zip(to_sub, new_ex):
		expr = expr.xreplace({old:new})

	return expr

def substitute_pow(expr): 
	# Substitutes Pow(x, i/2) for a call to Pow(sqrt(x), i) for all integers i.
	# Then substitutes Pow(x, i) for a call to powI(x) for integers 8 >= i > 0.
	# This is done for speed.

	### Handle half-integer powers

	# Search for Pow
	to_sub = [p for p in preorder_traversal(expr) if p.func == Pow]

	# Search for Rational arguments
	to_sub = [p for p in to_sub if isinstance(p.args[1], Rational)]

	# Search for half-integer arguments
	to_sub = [p for p in to_sub if p.args[1].q == Integer(2)]

	# Build substitutions
	# We use ssqrt to prevent sympy from simplifying them away
	new_ex = [Pow(ssqrt(p.args[0]), p.args[1].p) for p in to_sub]

	# Substitute. We use xreplace because this does a direct replacement
	# of only exactly matching subexpressions, rather than what subs does
	# (reason to try to shoehorn our replacement in) or what replace does
	# (wildcard matching).
	for (old,new) in zip(to_sub, new_ex):
		expr = expr.xreplace({old:new})

	### Turn negative integer powers into 1/x**(positive integer)

	# Search for instances of Pow
	to_sub = [p for p in preorder_traversal(expr) if p.func == Pow]

	# Search for Integer arguments
	to_sub = [p for p in to_sub if isinstance(p.args[1], Integer)]

	# Search for negative integer arguments less than 0
	to_sub = [p for p in to_sub if p.args[1] < 0]

	# Build substitutions
	new_ex = [powm1(Pow(p.args[0],-p.args[1])) for p in to_sub]

	# Substitute. We use xreplace because this does a direct replacement
	# of only exactly matching subexpressions, rather than what subs does
	# (reason to try to shoehorn our replacement in) or what replace does
	# (wildcard matching).
	for (old,new) in zip(to_sub, new_ex):
		expr = expr.xreplace({old:new})

	### Turn positive integer powers into powN calls

	# Search for instances of Pow
	to_sub = [p for p in preorder_traversal(expr) if p.func == Pow]

	# Search for Integer arguments
	to_sub = [p for p in to_sub if isinstance(p.args[1], Integer)]

	# Search for positive integer arguments less than 9
	to_sub = [p for p in to_sub if p.args[1] > 0 and p.args[1] <= 8]


	# Build substitutions
	new_ex = [powN[p.args[1]](p.args[0]) for p in to_sub]

	# Substitute. We use xreplace because this does a direct replacement
	# of only exactly matching subexpressions, rather than what subs does
	# (reason to try to shoehorn our replacement in) or what replace does
	# (wildcard matching).
	for (old,new) in zip(to_sub, new_ex):
		expr = expr.xreplace({old:new})

	return expr

def py_to_fort(expr):

	# Replace all instances of a**b with pow(a,b)
	# We made the sympy output use PPow, so now we just replace PPow -> pow
	expr = expr.replace('PPow', 'pow')

	# We made the sympy output use ssqrt, so now we just replace ssqrt -> pow
	expr = expr.replace('ssqrt', 'sqrt')

	# The sign function used by sympy isn't native to Fortran. We've written a
	# sgn function that does the same job, located in private/support_functions.f90.
	# So we divert calls from sign -> sgn.
	expr = expr.replace('sign(', 'sgn(')

	# Next we take advantage of MESA/const having pre-computing log10
	expr = expr.replace('safe_log(10)', 'ln10')
	expr = expr.replace('log(10)', 'ln10')

	# Next we replace 'colon' with ':' because sympy doesn't like colons in the middle of variable names.
	expr = expr.replace('colon', ':')

	return expr

def wrap_element(x):
	# If x is not iterable, wrap it in a list.
	# Otherwise, return as-is.

	if isinstance(x, Iterable):
		return x
	else:
		return [x]

def indent_list(lines):
	return [tab + l for l in lines]

def indent_list_of_str(strings):
	ret = []
	for s in strings:
		s = s.split('\n')
		s = '\n'.join(tab + l for l in s)
		ret.append(s)
	return ret
