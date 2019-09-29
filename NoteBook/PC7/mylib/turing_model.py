import numpy as np


class turing_model(object):

    def __init__(self, a, b, xmin, xmax, nx) :
        self.a = a
        self.b = b
        self.du = 1.0
        self.dv = 1.5
        self.delta = 8.0
        self.xmin = xmin
        self.xmax = xmax
        self.nx = nx
        self.dx = (xmax-xmin)/(nx-1)

    def fcn_radau(self, n, t, y, ydot, rpar, ipar):
        a = self.a
        b = self.b
        du = self.du
        dv = self.dv
        delta = self.delta
        nx = self.nx
        dx = self.dx
        duoverdxdx = du/(dx*dx)
        dvoverdxdx = dv/(dx*dx)

        ui   = y[0]
        vi   = y[1]
        uip1 = y[2]
        vip1 = y[3]
        ydot[0] =           duoverdxdx*(-ui + uip1) + (a - ui - (4*ui*vi)/(1+ui*ui))
        ydot[1] = delta * ( dvoverdxdx*(-vi + vip1) +  b*(ui - (ui*vi)/(1+ui*ui)) )

        for ix in range(1, nx-1):
            irow = ix*2
            
            uim1 = y[irow-2]  
            vim1 = y[irow-1]  
            ui   = y[irow]  
            vi   = y[irow+1]  
            uip1 = y[irow+2]  
            vip1 = y[irow+3]  
            ydot[irow]   =           duoverdxdx*(uim1 -2*ui + uip1) + (a - ui - (4*ui*vi)/(1+ui*ui))
            ydot[irow+1] = delta * ( dvoverdxdx*(vim1 -2*vi + vip1) +  b*(ui - (ui*vi)/(1+ui*ui)) )

        uim1 = y[2*nx -4]
        vim1 = y[2*nx -3]
        ui   = y[2*nx -2]
        vi   = y[2*nx -1]
        ydot[2*nx-2] =           duoverdxdx*(uim1 -ui) + (a - ui - (4*ui*vi)/(1+ui*ui))
        ydot[2*nx-1] = delta * ( dvoverdxdx*(vim1 -vi) +  b*(ui - (ui*vi)/(1+ui*ui)) )
