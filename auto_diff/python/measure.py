'''
Sympy allows us to simplify expressions.
Instead of having it simplify using its built-in notion of what
a simple expression is, we choose to simplify to minimize computational cost.
'''

from sympy import *
from sympy.abc import x,y
from utils import substitute_pow, powN
import numpy as np

# 'basic' here means roughly a one-cycle op.
# 'div' is division, which takes ~30 cycles.
# 'special' is a special function, which takes ~1000 cycles.
# DIRACDELTA and DERIVATIVE get eliminated in post-processing and so are free.
special = 1000
div = 30
basic = 1
weights = {
	'SIN': special,
	'COS': special,
	'TAN': special,
	'TANH': special,
	'COSH': special,
	'SINH': special,
	'ASIN': special,
	'ACOS': special,
	'ATAN': special,
	'ATANH': special,
	'ACOSH': special,
	'ASINH': special,
	'EXP': special,
	'LOG': special,
	'POW': special,
	'ADD': basic,
	'MUL': basic,
	'NEG': basic,
	'SUB': basic,
	'HEAVISIDE': basic,
	'ABS': basic,
	'DIV': div,
	'SGN': basic,
	'POWM1': div,
	'SSQRT': special,
	'DIRACDELTA': 0,
	'DERIVATIVE': 0
}

# The cost of powN is a multiply times ~lnN because it's done by building up x^{1,2,4,...}
# and combining these.
for i,p in enumerate(powN):
	weights['POW' + str(i)] = basic * np.log2(i+1)

ln10 = Symbol('Q')
ln2 = Symbol('R')

def weighted_count_ops(expr_original, verbose=False):
	expr = substitute_pow(expr_original)

	# In post-processing we'll replace log(10) and log(2) with constants, so
	# don't count them towards the cost.
	expr = expr.xreplace({log(10):ln10, log(2):ln2}) 

	# Break up the expression into components we can count.
	vis = expr.count_ops(visual=True)
	components = list(vis.args)

	if len(components) == 0:
		return 0

	if isinstance(components[0], Rational): # Means there's just one operation and it occurs repeatedly.
		parts = {str(components[1]):components[0]} # First argument is the number of times the op was used. Second is op name.
	else:
		parts = {}
		for a in components:
			if len(a.args) == 0:
				parts[str(a)] = 1
			else:
				parts[str(a.args[1])] = a.args[0] # First argument is the number of times the op was used. Second is op name.

	if verbose:
		print(expr_original)
		print(expr)
		print(srepr(expr))
		print(components)
		print(parts)

	measure = 0
	for key,val in parts.items():
		measure += weights[key] * val

	return measure

