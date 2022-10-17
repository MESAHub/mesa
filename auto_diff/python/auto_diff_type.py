from sympy import symbols

from chain_rule import unary_generic_chain_rule, binary_generic_chain_rule, \
    unary_specific_chain_rule, binary_specific_chain_rule
from partial import Partial
from routine import FortranRoutine, FortranFunction
from utils import tab


class AutoDiffType:
    def __init__(self, name, array, array_length, partials):
        '''
        Stores a list of partials that is complete, in the sense that there
        is enough information
        to compute the chain rule within that set of partials, and sorted by
        total order.
        '''
        
        self.name = name
        self.num_independent = list(partials)[0].num_independent
        self.array = list(partials)[0].array
        self.array_length = array_length
        
        # Complete the partials list
        partials = set(partials)
        complete = False
        while not complete:
            ps = list(partials)
            for p in ps:
                partials.update(p.completion_partials())
            if len(partials) == len(ps):
                complete = True
        
        self.partials = sorted(list(partials), key=lambda p: [p.net_order,
                                                              tuple(-o for o in
                                                                    p.orders)])
        self.max_order = max(p.net_order for p in self.partials)
        
        self.unary_partials = sorted(
            list(Partial((i,), False) for i in range(self.max_order + 1)),
            key=lambda p: [p.net_order, tuple(-o for o in p.orders)])
        self.binary_partials = sorted(list(
            Partial((i, j), False) for i in range(self.max_order + 1) for j in
            range(self.max_order + 1) if i + j <= self.max_order),
                                      key=lambda p: [p.net_order, tuple(
                                          -o for o in p.orders)])
    
    def declare_name(self, ref=None):
        if (not self.array) or (
        not self.array_length is None):  # Either no array type or else it's
            # a fixed length.
            return 'type(' + str(self.name) + ')'
        else:
            # For variable-length array types we have to know the length at
            # declare-type or else declare the length as deferred.
            # To declare it as referred pass ref='*'. Otherwise pass 'x' (or
            # 'y', 'z', etc.) to specify which
            # variable the length is derived from.
            if ref != '*':
                ref = ref + '%n'
            return 'type(' + str(self.name) + '(n=' + str(ref) + '))'
    
    def type(self):
        '''
        Returns the fortran type declaration for this type.
        '''
        
        # Construct header
        if self.array and self.array_length is None:
            header = 'type :: ' + self.name + '(n)'
        else:
            header = 'type :: ' + self.name
        
        # Construct variable declarations
        variables = []
        if self.array and self.array_length is None:
            variables.append(
                'integer, len :: n ! Length of array of independent variables')
        
        for p in self.partials:
            p_name = str(p)
            if 'Array' in p_name:
                if self.array_length is None:
                    variables.append('real(dp) :: ' + p_name + '(n)')
                else:
                    variables.append('real(dp) :: ' + p_name + '(' + str(
                        self.array_length) + ')')
            else:
                variables.append('real(dp) :: ' + p_name)
        
        # Construct ender
        ender = 'end type ' + self.name
        
        # Put it all together
        lines = [header] + [tab + v for v in variables] + [ender]
        
        return lines
    
    def partial_str_in_instance(self, instance_name, p):
        '''
        Returns a string representing the specified partial derivative as a
        member of an instance of this auto_diff type.
        '''
        s = str(p)
        s = instance_name + '%' + s
        if 'Array' in s:
            if self.array_length is None:
                s = s + '(1:' + instance_name + '%n)'
            else:
                s = s + '(1:' + str(self.array_length) + ')'
        return s
    
    def assignment_routine_from_self(self):
        '''
        Constructs a routine for assigning one object of this type to another.
        '''
        routine_name = 'assign_from_self'
        routine_arguments = [('this', self.declare_name(ref='*'), 'out'),
                             ('other', self.declare_name(ref='this'), 'in')]
        routine_body = ['this%' + str(p) + ' = other%' + str(p) for p in
                        self.partials]
        
        return FortranRoutine(routine_name, routine_arguments, routine_body)
    
    def assignment_routine_from_real_dp(self):
        '''
        Construct a routine for assigning a real(dp) to this type.
        '''
        
        routine_name = 'assign_from_real_dp'
        routine_arguments = [('this', self.declare_name(ref='*'), 'out'),
                             ('other', 'real(dp)', 'in')]
        routine_body = ['this%val = other'] + ['this%' + str(p) + ' = 0_dp' for
                                               p in self.partials if
                                               p.net_order > 0]
        
        return FortranRoutine(routine_name, routine_arguments, routine_body)
    
    def assignment_routine_from_int(self):
        '''
        Construct a routine for assigning an integer to this type.
        '''
        
        # Sanitize function name in case the other type is a real(dp). Note
        # similiar special case handling will be needed for quad precision.
        routine_name = 'assign_from_int'
        routine_arguments = [('this', self.declare_name(ref='*'), 'out'),
                             ('other', 'integer', 'in')]
        routine_body = ['this%val = other'] + ['this%' + str(p) + ' = 0_dp' for
                                               p in self.partials if
                                               p.net_order > 0]
        
        return FortranRoutine(routine_name, routine_arguments, routine_body)
    
    def comparison_function_with_self(self, operator_name, operator):
        '''
        Makes a function for comparing two objects of this type using the
        specified operator.
        '''
        
        function_name = operator_name + '_self'
        function_arguments = [('this', self.declare_name(ref='*'), 'in'),
                              ('other', self.declare_name(ref='this'), 'in')]
        function_result = ('z', 'logical')
        function_body = ['z = (this%val ' + operator + ' other%val)']
        
        return FortranFunction(function_name, function_arguments,
                               function_result, function_body)
    
    def comparison_function_with_real_dp(self, operator_name, operator):
        '''
        Returns functions for comparing an object of this type with a real(
        dp) and vice-versa.
        
        '''
        
        function_name = operator_name + '_' + self.name + '_' + 'real_dp'
        function_arguments = [('this', self.declare_name(ref='*'), 'in'),
                              ('other', 'real(dp)', 'in')]
        function_result = ('z', 'logical')
        function_body = ['z = (this%val ' + operator + ' other)']
        
        first_function = FortranFunction(function_name, function_arguments,
                                         function_result, function_body)
        
        function_name = operator_name + '_' + 'real_dp' + '_' + self.name
        function_arguments = [('this', 'real(dp)', 'in'),
                              ('other', self.declare_name(ref='*'), 'in')]
        function_result = ('z', 'logical')
        function_body = ['z = (this ' + operator + ' other%val)']
        
        second_function = FortranFunction(function_name, function_arguments,
                                          function_result, function_body)
        
        return first_function, second_function
    
    def comparison_function_with_int(self, operator_name, operator):
        '''
        Returns functions for comparing an object of this type with an
        integer and vice-versa.
        
        '''
        
        function_name = operator_name + '_' + self.name + '_' + 'int'
        function_arguments = [('this', self.declare_name(ref='*'), 'in'),
                              ('other', 'integer', 'in')]
        function_result = ('z', 'logical')
        function_body = ['z = (this%val ' + operator + ' other)']
        
        first_function = FortranFunction(function_name, function_arguments,
                                         function_result, function_body)
        
        function_name = operator_name + '_' + 'int' + '_' + self.name
        function_arguments = [('this', 'integer', 'in'),
                              ('other', self.declare_name(ref='*'), 'in')]
        function_result = ('z', 'logical')
        function_body = ['z = (this ' + operator + ' other%val)']
        
        second_function = FortranFunction(function_name, function_arguments,
                                          function_result, function_body)
        
        return first_function, second_function
    
    def derivative_function(self, variable_index):
        '''
        Returns a function which takes a derivative with respect to the
        specified variable.
        '''
        
        function_name = 'differentiate_' + self.name + '_' + str(
            variable_index + 1)
        function_arguments = [('this', self.declare_name(ref='*'), 'in')]
        function_result = ('derivative', self.declare_name(ref='this'))
        
        function_body = []
        for p in self.partials:
            p_increment = p.increment_order(variable_index)
            
            if p_increment in self.partials:  # The higher-order derivative
                # exists.
                line = 'derivative%' + str(p) + ' = this%' + str(p_increment)
            else:  # Fill with zeros otherwise.
                line = 'derivative%' + str(p) + ' = 0_dp'
            
            function_body.append(line)
        
        return FortranFunction(function_name, function_arguments,
                               function_result, function_body)
    
    def generic_unary_operator_function(self):
        '''
        Returns a function which implements the chain rule on a general
        unary operator.
        '''
        
        function_name = 'make_unary_operator'
        function_arguments = [('x', self.declare_name(ref='*'), 'in')]
        function_arguments = function_arguments + [
            ('z_' + str(p).replace('val1', 'x'), 'real(dp)', 'in') for p in
            self.unary_partials]
        function_result = ('unary', self.declare_name(ref='x'))
        function_body, function_declarations = unary_generic_chain_rule(self,
                                                                        fixed_length=self.array_length)
        function_body = function_declarations + function_body
        
        return FortranFunction(function_name, function_arguments,
                               function_result, function_body)
    
    def generic_binary_operator_function(self):
        '''
        Returns a function which implements the chain rule on a general
        unary operator.
        '''
        
        function_name = 'make_binary_operator'
        function_arguments = [('x', self.declare_name(ref='*'), 'in'),
                              ('y', self.declare_name(ref='*'), 'in')]
        function_arguments = function_arguments + [('z_' + str(p).replace(
            'val1', 'x').replace('val2', 'y'), 'real(dp)', 'in') for p in
                                                   self.binary_partials]
        function_result = ('binary', self.declare_name(ref='x'))
        function_body, function_declarations = binary_generic_chain_rule(self,
                                                                         fixed_length=self.array_length)
        function_body = function_declarations + function_body
        
        return FortranFunction(function_name, function_arguments,
                               function_result, function_body)
    
    def specific_unary_operator_function(self, operator_name, operator):
        '''
        Returns a function which implements the specified unary operator.
        '''
        
        function_name = operator_name + '_self'
        function_arguments = [('x', self.declare_name(ref='*'), 'in')]
        function_result = ('unary', self.declare_name(ref='x'))
        function_body, function_declarations = unary_specific_chain_rule(self,
                                                                         operator,
                                                                         fixed_length=self.array_length)
        function_body = function_declarations + function_body
        
        # Special case handling for safe_log
        if 'safe' in operator_name:
            for i in range(len(function_body)):
                function_body[i] = function_body[i].replace('log', 'safe_log')
        
        return FortranFunction(function_name, function_arguments,
                               function_result, function_body)
    
    def specific_binary_operator_function_self(self, operator_name, operator):
        '''
        Returns a function which implements the specified binary operator.
        '''
        
        function_name = operator_name + '_self'
        function_arguments = [('x', self.declare_name(ref='*'), 'in'),
                              ('y', self.declare_name(ref='*'), 'in')]
        function_result = ('binary', self.declare_name(ref='x'))
        function_body, function_declarations = binary_specific_chain_rule(self,
                                                                          operator,
                                                                          fixed_length=self.array_length)
        function_body = function_declarations + function_body
        
        return FortranFunction(function_name, function_arguments,
                               function_result, function_body)
    
    def specific_binary_operator_function_real_dp(self, operator_name,
                                                  operator):
        '''
        Returns functions which implements the specified binary operator with a real(dp) object
        in the two possible orders (this,real), (real,this).
        This is done by producing a unary operator that implements operator(this type, real), treating the
        real as a constant symbol, and then using the specific operator unary chain rule on that.
        '''
        
        # (this, real)
        function_name = operator_name + '_self_real'
        function_arguments = [('x', self.declare_name(ref='*'), 'in'),
                              ('y', 'real(dp)', 'in')]
        
        y = symbols('y', real=True)
        unary_operator = lambda x: operator(x, y)
        
        function_result = ('unary', self.declare_name(ref='x'))
        function_body, function_declarations = unary_specific_chain_rule(self,
                                                                         unary_operator,
                                                                         fixed_length=self.array_length)
        function_body = function_declarations + function_body
        function_1 = FortranFunction(function_name, function_arguments,
                                     function_result, function_body)
        
        # (real, this)
        function_name = operator_name + '_real_self'
        function_arguments = [('z', 'real(dp)', 'in'),
                              ('x', self.declare_name(ref='*'), 'in')]
        
        y = symbols('z', real=True)
        unary_operator = lambda x: operator(y, x)
        
        function_result = ('unary', self.declare_name(ref='x'))
        function_body, function_declarations = unary_specific_chain_rule(self,
                                                                         unary_operator,
                                                                         fixed_length=self.array_length)
        function_body = function_declarations + function_body
        function_2 = FortranFunction(function_name, function_arguments,
                                     function_result, function_body)
        
        return function_1, function_2
    
    def specific_binary_operator_function_int(self, operator_name, operator):
        '''
        Returns functions which implements the specified binary operator with an integer.
        in the two possible orders (this,int), (int,this).
        This is done by producing a unary operator that implements operator(this type, int), treating the
        int as a constant symbol, and then using the specific operator unary chain rule on that.
        '''
        
        # (this, int)
        function_name = operator_name + '_self_int'
        function_arguments = [('x', self.declare_name(ref='*'), 'in'),
                              ('y', 'integer', 'in')]
        
        y = symbols('y_dp', real=True)
        unary_operator = lambda x: operator(x, y)
        
        function_result = ('unary', self.declare_name(ref='x'))
        function_body, function_declarations = unary_specific_chain_rule(self,
                                                                         unary_operator,
                                                                         fixed_length=self.array_length)
        function_body = ['real(dp) :: y_dp'] + function_declarations + [
            'y_dp = y'] + function_body
        function_1 = FortranFunction(function_name, function_arguments,
                                     function_result, function_body)
        
        # (int, this)
        function_name = operator_name + '_int_self'
        function_arguments = [('z', 'integer', 'in'),
                              ('x', self.declare_name(ref='*'), 'in')]
        
        y = symbols('y_dp', real=True)
        unary_operator = lambda x: operator(y, x)
        
        function_result = ('unary', self.declare_name(ref='x'))
        function_body, function_declarations = unary_specific_chain_rule(self,
                                                                         unary_operator,
                                                                         fixed_length=self.array_length)
        function_body = ['real(dp) :: y_dp'] + function_declarations + [
            'y_dp = z'] + function_body
        
        function_2 = FortranFunction(function_name, function_arguments,
                                     function_result, function_body)
        
        return function_1, function_2
