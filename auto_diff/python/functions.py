from sympy import *
from sympy.codegen.cfunctions import *

pi = symbols('pi')

# Supported unary functions
unary_operators = [
	(lambda x: sign(x), 'sign'),
	(lambda x: sqrt(x * Heaviside(x)), 'safe_sqrt'),
	(lambda x: -1*x, 'unary_minus'),
	(lambda x: exp(x), 'exp'),
	(lambda x: expm1(x), 'expm1'),
	(lambda x: 10**x, 'exp10'),
	(lambda x: 1/x, 'powm1'),
	(lambda x: log(x), 'log'),
	(lambda x: log1p(x), 'log1p'),
	(lambda x: log(x), 'safe_log'),
	(lambda x: log(x,10), 'log10'),
	(lambda x: log(x,10), 'safe_log10'),
	(lambda x: log(x,2), 'log2'),
	(lambda x: sin(x), 'sin'),
	(lambda x: cos(x), 'cos'),
	(lambda x: tan(x), 'tan'),
	(lambda x: sin(pi*x), 'sinpi'),
	(lambda x: cos(pi*x), 'cospi'),
	(lambda x: tan(pi*x), 'tanpi'),
	(lambda x: sinh(x), 'sinh'),
	(lambda x: cosh(x), 'cosh'),
	(lambda x: tanh(x), 'tanh'),
	(lambda x: asin(x), 'asin'),
	(lambda x: acos(x), 'acos'),
	(lambda x: atan(x), 'atan'),
	(lambda x: asin(x)/pi, 'asinpi'),
	(lambda x: acos(x)/pi, 'acospi'),
	(lambda x: atan(x)/pi, 'atanpi'),
	(lambda x: asinh(x), 'asinh'),
	(lambda x: acosh(x), 'acosh'),
	(lambda x: atanh(x), 'atanh'),
	(lambda x: sqrt(x), 'sqrt'),
	(lambda x: x**2, 'pow2'),
	(lambda x: x**3, 'pow3'),
	(lambda x: x**4, 'pow4'),
	(lambda x: x**5, 'pow5'),
	(lambda x: x**6, 'pow6'),
	(lambda x: x**7, 'pow7'),
	(lambda x: x**8, 'pow8'),
	(lambda x: abs(x), 'abs')
]

# Supported binary functions

def Dim(x,y):
	return (x-y+abs(x-y))/Float(2)

def Sub(x,y):
	return x-y

def Div(x,y):
	return x/y

binary_operators = [
	(Add, 'add'),
	(Sub, 'sub'),
	(Mul, 'mul'),
	(Div, 'div'),
	(Pow, 'pow'),
	(Max, 'max'),
	(Min, 'min'),
	(Dim, 'dim')
]

# Supported comparison operators
comparison_operators = [
	('.eq.', 'equal'),
	('.ne.', 'neq'),
	('.gt.', 'greater'),
	('.lt.', 'less'),
	('.le.', 'leq'),
	('.ge.', 'geq')
]

# Names of functions that map onto intrinsic operators
intrinsics = {'add':'+', 'sub':'-', 'mul':'*', 'div':'/', 'unary_minus':'-'}
