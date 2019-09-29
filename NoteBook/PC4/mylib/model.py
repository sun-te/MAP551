import numpy as np


class vanderpol_model(object):

    def __init__(self, eps):
        self.eps = eps

    def fcn(self, t, y):
        y1, y2 = y
        eps = self.eps
        y1_dot = y2
        y2_dot = eps*(1-y1*y1)*y2 - y1
        return [y1_dot, y2_dot]

    def jac(self, y):
        y1, y2 = y
        eps = self.eps
        return np.array([[0, 1], [-1 -2*eps*y1*y2, eps*(1-y1*y1)]])

class brusselator_model(object): 
 
    def __init__(self, a, b) : 
        self.a = a 
        self.b = b 
 
    def fcn(self, t, y):
        y1, y2 = y
        a = self.a
        b = self.b 
        y1_dot = 1 - (b+1)*y1 + a*y1*y1*y2 
        y2_dot = b*y1 - a*y1*y1*y2  
        return (y1_dot, y2_dot)

    def jac(self, y):
        y1, y2 = y
        a = self.a
        b = self.b 
        return np.array([[2*a*y1*y2 -(b+1), a*y1*y1], [-2*a*y1*y2 + b, -a*y1*y1]])
