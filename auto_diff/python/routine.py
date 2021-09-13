from utils import py_to_fort, tab

class FortranRoutine:
	def __init__(self, name, arguments, body):
		self.name = name # Routine name
		self.arguments = arguments # List of tuples of the form (argument_name, argument_declaration_string, argument_intent)
		self.body = body # List of lines forming the body of the routine.

	def __str__(self):
		'''
		Returns a string representing a valid fortran implementation of this routine.
		This involves constructing some boilerplate (the routine declaration and ending), as well as all
		variable declarations, and putting it together in the right order.
		'''

		# Construct the header
		header = 'subroutine ' + self.name
		args = ', '.join(arg[0] for arg in self.arguments)
		header = header + '(' + args + ')'

		# Construct declarations
		declarations = [tab + arg[1] + ', intent(' + arg[2] + ') :: ' + arg[0] for arg in self.arguments]

		# Indent body
		body = [tab + b for b in self.body]

		# Construct ender
		ender = 'end subroutine ' + self.name

		lines = [header] + declarations + body + [ender]
		return py_to_fort('\n'.join(lines))

class FortranFunction:
	def __init__(self, name, arguments, result, body):
		self.name = name # Function name
		self.arguments = arguments # List of tuples of the form (argument_name, argument_declaration_string, argument_intent)
		self.result = result # Tuple (argument_name, argument_declaration_string).
		self.body = body # List of lines forming the body of the function.

	def __str__(self):
		'''
		Returns a string representing a valid fortran implementation of this routine.
		This involves constructing some boilerplate (the routine declaration and ending), as well as all
		variable declarations, and putting it together in the right order.
		Note that because this is a function there is an explicitly listed return argument, which forms part of the beginning boilerplate.
		'''

		# Construct the header
		header = 'function ' + self.name
		args = ', '.join(arg[0] for arg in self.arguments)
		header = header + '(' + args + ')'
		header = header + ' result(' + self.result[0] + ')'

		# Construct declarations
		declarations = [tab + arg[1] + ', intent(' + arg[2] + ') :: ' + arg[0] for arg in self.arguments]
		declarations.append(tab + self.result[1] + ' :: ' + self.result[0])

		# Indent body
		body = [tab + b for b in self.body]

		# Construct ender
		ender = 'end function ' + self.name

		lines = [header] + declarations + body + [ender]
		return py_to_fort('\n'.join(lines))