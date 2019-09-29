import numpy as np

class nagumo_model(object): 
 
    def __init__(self, k, d, xmin, xmax, nx) : 
        self.k = k
        self.d = d
        self.xmin = xmin 
        self.xmax = xmax
        self.nx = nx
        #self.nx = nx-2 
        self.dx = (xmax-xmin)/(nx-1) 
 
    def fcn(self, t, y):
        k = self.k
        d = self.d
        nx = self.nx
        dx = self.dx
        doverdxdx = d/(dx*dx)

        ydot = np.zeros(y.size)

        ydot[0] = doverdxdx*(-y[0] + y[1]) + k*y[0]*y[0]*(1 - y[0])
        for ix in range(1, nx-1):
            ydot[ix] = doverdxdx*(y[ix-1] - 2*y[ix] + y[ix+1]) + k*y[ix]*y[ix]*(1 - y[ix])
        ydot[nx-1] = doverdxdx*(y[nx-2] - y[nx-1]) + k*y[nx-1]*y[nx-1]*(1 - y[nx-1])

        return ydot

    def fcn_radau(self, n, t, y, ydot, rpar, ipar):
        k = self.k
        d = self.d
        nx = self.nx
        dx = self.dx
        doverdxdx = d/(dx*dx)

        ydot[0] = doverdxdx*(-y[0] + y[1]) + k*y[0]*y[0]*(1 - y[0])
        for ix in range(1, nx-1):
            ydot[ix] = doverdxdx*(y[ix-1] - 2*y[ix] + y[ix+1]) + k*y[ix]*y[ix]*(1 - y[ix])
        ydot[nx-1] = doverdxdx*(y[nx-2] - y[nx-1]) + k*y[nx-1]*y[nx-1]*(1 - y[nx-1])

    def fcn_reac_radau(self, n, t, y, ydot, rpar, ipar):
        k = self.k
        for ix in range(n[0]):
            ydot[ix] = k*y[ix]*y[ix]*(1 - y[ix])

    def fcn_diff_radau(self, n, t, y, ydot, rpar, ipar):
        k = self.k
        d = self.d
        nx = self.nx
        dx = self.dx
        doverdxdx = d/(dx*dx)

        ydot[0] = doverdxdx*(-y[0] + y[1])
        for ix in range(1, nx-1):
            ydot[ix] = doverdxdx*(y[ix-1] - 2*y[ix] + y[ix+1])
        ydot[nx-1] = doverdxdx*(y[nx-2] - y[nx-1])

    def fcn_rock(self, n, t, y, ydot):
        k = self.k
        d = self.d
        nx = self.nx
        dx = self.dx
        doverdxdx = d/(dx*dx)

        ydot[0] = doverdxdx*(-y[0] + y[1]) + k*y[0]*y[0]*(1 - y[0])
        for ix in range(1, nx-1):
            ydot[ix] = doverdxdx*(y[ix-1] - 2*y[ix] + y[ix+1]) + k*y[ix]*y[ix]*(1 - y[ix])
        ydot[nx-1] = doverdxdx*(y[nx-2] - y[nx-1]) + k*y[nx-1]*y[nx-1]*(1 - y[nx-1])

    def fcn_diff_rock(self, n, t, y, ydot):
        d = self.d
        nx = self.nx
        dx = self.dx
        doverdxdx = d/(dx*dx)

        ydot[0] = doverdxdx*(-y[0] + y[1]) 
        for ix in range(1, nx-1):
            ydot[ix] = doverdxdx*(y[ix-1] - 2*y[ix] + y[ix+1])
        ydot[nx-1] = doverdxdx*(y[nx-2] - y[nx-1])

    def fcn_exact(self, t):
        k = self.k
        d = self.d
        xmin = self.xmin
        xmax = self.xmax
        nx = self.nx
        dx = self.dx
        x0 = -20.

        v = (1./np.sqrt(2.))*(np.sqrt(k*d))
        cst  = -(1./np.sqrt(2.))*(np.sqrt(k/d));

        ##x = np.linspace(xmin, xmax, nx)
        x = np.linspace(xmin+dx, xmax-dx, nx)
        y = np.exp(cst*(x-x0-v*t)) / (1. + np.exp(cst*(x-x0-v*t)))
        return y
