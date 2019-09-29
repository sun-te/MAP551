import numpy as np


class vanderpol_model:

    def __init__(self, eps):
        self.eps = eps

    def fcn(self, t, y):
        y1, y2 = y
        eps = self.eps
        y1_dot = y2
        y2_dot = eps*(1-y1*y1)*y2 - y1
        return np.array([y1_dot, y2_dot])

    def fcn_radau(self, n, t, y, f, rpar, ipar):
        eps = self.eps
        y1 = y[0]
        y2 = y[1]
        f[0] = y2
        f[1] = eps*(1-y1*y1)*y2 - y1

    def jac(self, t, y):
        y1, y2 = y
        eps = self.eps
        return np.array([[0, 1], [-1 -2*eps*y1*y2, eps*(1-y1*y1)]])

class brusselator_model: 
 
    def __init__(self, a, b) : 
        self.a = a 
        self.b = b 
 
    def fcn(self, t, y):
        y1, y2 = y
        a = self.a
        b = self.b 
        y1_dot = 1 - (b+1)*y1 + a*y1*y1*y2 
        y2_dot = b*y1 - a*y1*y1*y2  
        return np.array([y1_dot, y2_dot])

    def fcn_radau(self, n, t, y, f, rpar, ipar):
        a = self.a
        b = self.b 
        y1 = y[0]
        y2 = y[1]
        f[0] = 1 - (b+1)*y1 + a*y1*y1*y2
        f[1] = b*y1 - a*y1*y1*y2

    def jac(self, t, y):
        y1, y2 = y
        a = self.a
        b = self.b 
        return np.array([[2*a*y1*y2 -(b+1), a*y1*y1], [-2*a*y1*y2 + b, -a*y1*y1]])

class oregonator_model:

    def __init__(self, eps, mu, f, q, y10=0., y20=0., y30=0.):
        self.eps = eps
        self.mu = mu
        self.f = f
        self.q = q
        self.y10 = y10
        self.y20 = y20
        self.y30 = y30

    def jac(self, t, y):
        y1, y2, y3 = y
        eps = self.eps
        mu = self.mu
        f = self.f
        q = self.q
        return np.array([[-1   , 1                   , 0             ],
                         [0    , (1/eps)*(-y3+1-2*y2), (1/eps)*(q-y2)],
                         [ f/mu, -y3/mu              , (1/mu)*(-q-y2)]])

    def fcn(self, t, y):
        y1, y2, y3 = y
        eps = self.eps
        mu = self.mu
        f = self.f
        q = self.q
        y1_dot = y2 - y1
        y2_dot = (1/eps)*(q*y3 - y3*y2 + y2*(1-y2))
        y3_dot = (1/mu)*(-q*y3 - y3*y2 + f*y1)
        return np.array([y1_dot, y2_dot, y3_dot])

    def fcn_ps(self, t, y):
        y1, y2 = y
        eps = self.eps
        f = self.f
        q = self.q
        y1_dot = y2 - y1
        y2_dot = (1/eps)*(f*y1*(q -y2)/(q+y2) + y2*(1-y2))
        return np.array([y1_dot, y2_dot])

    def fcn_radau(self, n, t, y, ydot, rpar, ipar):
        eps = self.eps
        mu = self.mu
        f = self.f
        q = self.q
        y1 = y[0]
        y2 = y[1]
        y3 = y[2]
        ydot[0] = y2 - y1
        ydot[1] = (1/eps)*(q*y3 - y3*y2 + y2*(1-y2))
        ydot[2] = (1/mu)*(-q*y3 - y3*y2 + f*y1)

    def fcn_ps_radau(self, n, t, y, ydot, rpar, ipar):
        eps = self.eps
        f = self.f
        q = self.q
        y1 = y[0]
        y2 = y[1]
        ydot[0] = y2 - y1
        ydot[1] = (1/eps)*(f*y1*(q -y2)/(q+y2) + y2*(1-y2))

    def fcn_a1(self, t, y):
        y1, y2 = y
        eps = self.eps
        mu = self.mu
        f = self.f
        q = self.q
        y30 = self.y30
        y1_dot = y2 - y1
        y2_dot = (1/eps)*(q*y30 - y30*y2 + y2*(1-y2))
        return np.array([y1_dot, y2_dot])

    def fcn_a1_radau(self, n, t, y, ydot, rpar, ipar):
        eps = self.eps
        mu = self.mu
        f = self.f
        q = self.q
        y30 = self.y30
        y1 = y[0]
        y2 = y[1]
        ydot[0] = y2 - y1
        ydot[1] = (1/eps)*(q*y30 - y30*y2 + y2*(1-y2))

    def fcn_b1(self, t, y):
        y3 = y
        eps = self.eps
        mu = self.mu
        f = self.f
        q = self.q
        y10 = self.y10
        y20 = self.y20
        y3_dot = (1/mu)*(-q*y3 - y3*y20 + f*y10)
        return y3_dot

    def fcn_b1_radau(self, n, t, y, ydot, rpar, ipar):
        eps = self.eps
        mu = self.mu
        f = self.f
        q = self.q
        y10 = self.y10
        y20 = self.y20
        y30 = self.y30
        y3 = y[0]
        ydot[0] = (1/mu)*(-q*y3 - y3*y20 + f*y10)

    def fcn_a2(self, t, y):
        y1, y2 = y
        eps = self.eps
        mu = self.mu
        f = self.f
        q = self.q
        y1_dot = y2 - y1
        y2_dot = (1/eps)*(f*y1*((q-y2)/(q+y2)) + y2*(1-y2))
        return (y1_dot, y2_dot)

    def fcn_a2_radau(self, n, t, y, ydot, rpar, ipar):
        eps = self.eps
        mu = self.mu
        f = self.f
        q = self.q
        y1 = y[0]
        y2 = y[1]
        ydot[0] = y2 - y1
        ydot[1] = (1/eps)*(f*y1*((q-y2)/(q+y2)) + y2*(1-y2))

    def fcn_b2(self, t, y):
        y2, y3 = y
        eps = self.eps
        mu = self.mu
        f = self.f
        q = self.q
        y10 = self.y10
        y2_dot = (1/eps)*(q-y2)*(y3 - (f*y10)/(q+y2))
        y3_dot = (1/mu)*( (-q-y2)*y3 + f*y10 )
        return (y2_dot, y3_dot)

    def fcn_b2_radau(self, n, t, y, ydot, rpar, ipar):
        eps = self.eps
        mu = self.mu
        f = self.f
        q = self.q
        y10 = self.y10
        y2 = y[0]
        y3 = y[1]
        ydot[0] = (1/eps)*(q-y2)*(y3 - (f*y10)/(q+y2))
        ydot[1] = (1/mu)*( (-q-y2)*y3 + f*y10 )

class heat_model:

    def __init__(self, xmin, xmax, nx) :
        self.xmin = xmin
        self.xmax = xmax
        self.nx = nx
        self.dx = (xmax-xmin)/(nx-1)

    def fcn(self, t, y):
        nx = self.nx
        dx = self.dx
        oneondxdx = 1/(dx*dx)

        ydot = np.empty(nx)
        ydot[0] = oneondxdx * (-y[0] + y[1])
        for ix in range(1, nx-1):
            ydot[ix] = oneondxdx * (y[ix-1] - 2*y[ix] + y[ix+1])
        ydot[nx-1] = oneondxdx * (y[nx-2] - y[nx-1])
        return ydot

    def fcn_radau(self, n, t, y, f, rpar, ipar):
        nx = self.nx
        dx = self.dx
        oneondxdx = 1/(dx*dx)

        f[0] = oneondxdx * (-y[0] + y[1])
        for ix in range(1, nx-1):
            f[ix] = oneondxdx * (y[ix-1] - 2*y[ix] + y[ix+1])
        f[nx-1] = oneondxdx * (y[nx-2] - y[nx-1])

    def fcn_rock(self, n, t, y, dy):
        nx = self.nx
        dx = self.dx
        oneondxdx = 1/(dx*dx)

        dy[0] = oneondxdx * (-y[0] + y[1])
        for ix in range(1, nx-1):
            dy[ix] = oneondxdx * (y[ix-1] - 2*y[ix] + y[ix+1])
        dy[nx-1] = oneondxdx * (y[nx-2] - y[nx-1])

    def fcn_exact(self, t):
        xmin = self.xmin
        xmax = self.xmax
        nx = self.nx
        x = np.linspace(xmin, xmax, nx)
        y = (1/(2*np.sqrt(np.pi*t))) * np.exp(-(x*x)/(4.*t))
        return y
