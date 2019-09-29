import numpy as np

class thermal_explosion_model(object): 
 
    def __init__(self, lamb, xmin, xmax, nx) : 
        self.lamb = lamb
        self.xmin = xmin 
        self.xmax = xmax
        self.nx = nx-2 
        self.dx = (xmax-xmin)/(nx-1) 
 
    def fcn(self, t, y):
        lamb = self.lamb
        nx = self.nx
        dx = self.dx
        oneoverlambdxdx = lamb/(dx*dx)

        ydot = np.zeros(y.size)

        ydot[0] = oneoverlambdxdx*(-y[0] + y[1]) + np.exp(y[0])
        for ix in range(1, nx-1):
            ydot[ix] = oneoverlambdxdx*(y[ix-1] - 2*y[ix] + y[ix+1]) + np.exp(y[ix])
        ydot[nx-1] = oneoverlambdxdx*(y[nx-2] - y[nx-1]) + np.exp(y[nx-1])

        return ydot

    def fcn_radau(self, n, t, y, ydot, rpar, ipar):
        lamb = self.lamb
        nx = self.nx
        dx = self.dx
        oneoverlambdxdx = 1/(lamb*dx*dx)

        ydot[0] = oneoverlambdxdx*(-2*y[0] + y[1]) + np.exp(y[0])
        for ix in range(1, nx-1):
            ydot[ix] = oneoverlambdxdx*(y[ix-1] - 2*y[ix] + y[ix+1]) + np.exp(y[ix])
        ydot[nx-1] = oneoverlambdxdx*(y[nx-2] - 2*y[nx-1]) + np.exp(y[nx-1])

    def fcn_reac_radau(self, n, t, y, ydot, rpar, ipar):
        for ix in range(n[0]):
            ydot[ix] = np.exp(y[0])

    def fcn_rock(self, n, t, y, ydot):
        lamb = self.lamb
        nx = self.nx
        dx = self.dx
        oneoverlambdxdx = 1/(lamb*dx*dx)

        ydot[0] = oneoverlambdxdx*(-2*y[0] + y[1]) + + np.exp(y[0])
        for ix in range(1, nx-1):
            ydot[ix] = oneoverlambdxdx*(y[ix-1] - 2*y[ix] + y[ix+1]) + np.exp(y[ix])
        ydot[nx-1] = oneoverlambdxdx*(y[nx-2] - 2*y[nx-1]) + np.exp(y[nx-1])

    def fcn_diff_rock(self, n, t, y, ydot):
        lamb = self.lamb
        nx = self.nx
        dx = self.dx
        oneoverlambdxdx = 1/(lamb*dx*dx)

        ydot[0] = oneoverlambdxdx*(-y[0] + y[1]) 
        for ix in range(1, nx-1):
            ydot[ix] = oneoverlambdxdx*(y[ix-1] - 2*y[ix] + y[ix+1])
        ydot[nx-1] = oneoverlambdxdx*(y[nx-2] - y[nx-1]) 

class fuel_thermal_explosion_model(object):

    def __init__(self, lamb, eps, xmin, xmax, nx) :
        self.lamb = lamb
        self.eps = eps
        self.xmin = xmin
        self.xmax = xmax
        self.nx = nx-2
        self.dx = (xmax-xmin)/(nx-1)

    def fcn_radau(self, n, t, y, ydot, rpar, ipar):
        lamb = self.lamb
        eps = self.eps
        nx = self.nx
        dx = self.dx
        oneonlambdxdx = 1./(lamb*dx*dx)

        theta_i = y[0]
        Y_i = y[1]
        theta_ip1 = y[2]
        Y_ip1 = y[3]
        ydot[0] = oneonlambdxdx*(-2*theta_i + theta_ip1) + np.exp(theta_i)*Y_i
        ydot[1] = oneonlambdxdx*(-Y_i + Y_ip1) - eps * np.exp(theta_i)*Y_i

        for ix in range(1, nx-1):
            irow = ix*2
            
            theta_im1 = y[irow-2]  
            Y_im1 = y[irow-1]  
            theta_i = y[irow]  
            Y_i = y[irow+1]  
            theta_ip1 = y[irow+2]  
            Y_ip1 = y[irow+3]  

            ydot[irow]   = oneonlambdxdx*(theta_im1 -2*theta_i + theta_ip1) + np.exp(theta_i)*Y_i
            ydot[irow+1] = oneonlambdxdx*(Y_im1 -2*Y_i + Y_ip1) - eps * np.exp(theta_i)*Y_i


        theta_im1 = y[2*nx -4]
        Y_im1 = y[2*nx - 3]
        theta_i = y[2*nx -2]
        Y_i = y[2*nx - 1]
        ydot[2*nx-2] = oneonlambdxdx*(theta_im1 - 2*theta_i) + np.exp(theta_i)*Y_i
        ydot[2*nx-1] = oneonlambdxdx*(Y_im1 - Y_i) - eps * np.exp(theta_i)*Y_i
