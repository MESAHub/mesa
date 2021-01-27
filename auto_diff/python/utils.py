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
sgn = Function('sgn')
powN = list(Function('pow' + str(i)) for i in range(9))

def zero_function(*args):
	return 0

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
	to_sub = [p for p in preorder_traversal(expr) if isinstance(p, Rational) and p != 1 and p != -1]

	# Build substitutions
	new_ex = [symbols(str(float(p)) + '_dp') for p in to_sub]

	# Substitute. We use xreplace because this does a direct replacement
	# of only exactly matching subexpressions, rather than what subs does
	# (reason to try to shoehorn our replacement in) or what replace does
	# (wildcard matching).
	expr = expr.xreplace(dict({old:new for old,new in zip(*(to_sub,new_ex))}))


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
	expr = expr.xreplace(dict({old:new for old,new in zip(*(to_sub,new_ex))}))


	return expr

def substitute_pow(expr): 
	expr_prev = None
	while expr_prev != expr:
		# Replacements can interfere with each other.
		# For instance in changing (x**3+y**2)**6 into
		# pow6(pow3(x) + pow2(y)) we would identify the substitutions
		# (x**3+y**2)**6  -> pow6(x**3+y**2)
		# x**3 -> pow3(x)
		# y**2 -> pow2(y)
		# If we make the second substitution first then the first substitution
		# no longer matches (because it's looking for x**3+y**2 not pow3(x)+y**2).
		# There may be a clever way to handle this, but I'm just doing an 'iterate-till-it-stops-changing'
		# approach because I know that'll work.
		expr_prev = expr

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
		expr = expr.xreplace(dict({old:new for old,new in zip(*(to_sub,new_ex))}))


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
		expr = expr.xreplace(dict({old:new for old,new in zip(*(to_sub,new_ex))}))


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
		expr = expr.xreplace(dict({old:new for old,new in zip(*(to_sub,new_ex))}))

	return expr

def py_to_fort(expr):

	# Replace all instances of a**b with pow(a,b)
	# We made the sympy output use PPow, so now we just replace PPow -> pow
	expr = expr.replace('PPow', 'pow')

	# We made the sympy output use ssqrt, so now we just replace ssqrt -> pow
	expr = expr.replace('ssqrt', 'sqrt')

	# Next we take advantage of MESA/const having pre-computing log10
	expr = expr.replace('safe_log(10.0_dp)', 'ln10')
	expr = expr.replace('log(10.0_dp)', 'ln10')

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
