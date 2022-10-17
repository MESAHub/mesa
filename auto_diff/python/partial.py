class Partial:
    def __init__(self, orders, array):
        '''
        A partial is a list of the orders of derivatives with respect to the
        independent variables.
        For instance, orders=(0,1,0,2) corresponds to d^3/(dx_1 dx_3^2).

        If array is true, the last independent variable is actually an array
        of independent variables
        (i.e. a vector-valued independent variable).
        The variables the array contains are restricted to first-order
        derivatives.
        '''
        self.array = array
        self.num_independent = len(orders)
        self.orders = tuple(orders)
        self.net_order = sum(orders)
    
    def __hash__(self):
        return hash(self.orders) + hash(self.num_independent) + hash(
            self.array)
    
    def __str__(self):
        if self.net_order == 0:
            return 'val'
        
        parts = []
        for i, o in enumerate(self.orders):
            if o != 0:
                parts.append('d' + str(o) + self.val_name(i))
        name = '_'.join(parts)
        return name
    
    def __eq__(self, other):
        return (self.orders == other.orders and self.array == other.array)
    
    def val_name(self, index):
        '''
        Returns the name of the specified independent variable.
        '''
        if self.array and index == self.num_independent - 1:
            return 'Array'
        else:
            return 'val' + str(index + 1)
    
    def completion_partials(self):
        '''
        In order to apply the chain rule with this partial, we need all
        partials that are one-lower order (min zero).
        '''
        
        # Go through each independent variable and construct the partial
        # that's one-lower order.
        # Skip those where this pushes the order of any independent variable
        # below zero.
        completion = []
        for i in range(self.num_independent):
            if self.orders[i] > 0:
                orders = list(self.orders)
                orders[i] -= 1
                completion.append(Partial(orders, self.array))
        
        return completion
    
    def increment_order(self, variable_index):
        '''
        Returns the partial with the derivative order of the specified
        variable incremented.
        '''
        
        orders = list(self.orders)
        orders[variable_index] += 1
        orders = tuple(orders)
        return Partial(orders, self.array)
