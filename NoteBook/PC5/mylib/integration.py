import numpy as np
import ctypes as ct
from scipy.integrate import solve_ivp, ode
from scipy.optimize import fsolve
import copy

#############################################################################
class ode_result:
    def __init__(self, y, t, nfev):
        self.y = y
        self.t = t
        self.nfev = nfev 

#############################################################################
def backward_euler(tini, tend, nt, yini, fcn): 
 
    dt = (tend-tini) / (nt-1)
    t = np.linspace(tini, tend, nt)

    yini_array = np.array(yini)
    neq = yini_array.size

    y = np.zeros((neq, nt), order='F')
    y[:,0] = yini_array
 
    def g(uip1, *args): 
        uip, tip1 = args 
        return uip1 - uip - dt*fcn(tip1, uip1)
 
    for it, tn  in enumerate(t[:-1]):
        yn = y[:,it]
        y0 = yn + dt*fcn(tn, yn)
        # solve y[:,it+1] - y[:,it] - dt * fcn(tini + (it+1)*dt, y[:,it+1]) = 0 
        y[:,it+1] = fsolve(g, y0, (yn, tn+dt))

    nfev = nt-1
         
    return ode_result(y, t, nfev) 

#############################################################################
def forward_euler(tini, tend, nt, yini, fcn):

    dt = (tend-tini) / (nt-1)
    t = np.linspace(tini, tend, nt)

    yini_array = np.array(yini)
    neq = yini_array.size

    y = np.zeros((neq, nt), order='FORTRAN')
    y[:,0] = yini_array

    for it, tn  in enumerate(t[:-1]):
        yn = y[:,it]
        y[:,it+1] = yn + dt*fcn(tn, yn)

    nfev = nt-1

    return ode_result(y, t, nfev)

#############################################################################
def rk2(tini, tend, nt, yini, fcn):

    dt = (tend-tini) / (nt-1)
    t = np.linspace(tini, tend, nt)

    yini_array = np.array(yini)
    neq = yini_array.size

    y = np.zeros((neq, nt), order='FORTRAN')
    y[:,0] = yini_array

    for it, tn  in enumerate(t[:-1]):
        yn = y[:,it]
        k1 = fcn(tn, yn)
        k2 = fcn(tn + 0.5*dt, yn + dt*(0.5*k1))
        y[:,it+1] = yn + dt*k2

    nfev = 2*(nt-1)

    return ode_result(y, t, nfev)

#############################################################################
def rk3(tini, tend, nt, yini, fcn):

    dt = (tend-tini) / (nt-1)
    t = np.linspace(tini, tend, nt)

    yini_array = np.array(yini)
    neq = yini_array.size

    y = np.zeros((neq, nt), order='FORTRAN')
    y[:,0] = yini_array

    for it, tn  in enumerate(t[:-1]):
        yn = y[:,it]
        k1 = fcn(tn, yn)
        k2 = fcn(tn + 0.5*dt, yn + dt*(0.5*k1))
        k3 = fcn(tn + dt, yn + dt*(-k1 + 2*k2))
        y[:,it+1] = yn + (dt/6)*(k1+4*k2+k3)

    nfev = 3*(nt-1)

    return ode_result(y, t, nfev)

#############################################################################
def rk4(tini, tend, nt, yini, fcn):

    dt = (tend-tini) / (nt-1)
    t = np.linspace(tini, tend, nt)

    yini_array = np.array(yini)
    neq = yini_array.size

    y = np.zeros((neq, nt), order='FORTRAN')
    y[:,0] = yini_array

    for it, tn  in enumerate(t[:-1]):
        yn = y[:,it]
        k1 = fcn(tn, yn)
        k2 = fcn(tn + 0.5*dt, yn + dt*(0.5*k1))
        k3 = fcn(tn + 0.5*dt, yn + dt*(0.5*k2))
        k4 = fcn(tn + dt, yn + dt*k3)
        y[:,it+1] = yn + (dt/6)*(k1+2*k2+2*k3+k4)

    nfev = 4*(nt-1)

    return ode_result(y, t, nfev)

#############################################################################
def rk38(tini, tend, nt, yini, fcn):

    # butcher table of the "3/8 rule" 
    c2 = 1/3
    c3 = 2/3
    c4 = 1

    a21 = 1/3
    a31 = -1/3
    a32 = 1
    a41 = 1
    a42 = -1
    a43 = 1

    b1 = 1/8
    b2 = 3/8
    b3 = 3/8
    b4 = 1/8

    dt = (tend-tini) / (nt-1)
    t = np.linspace(tini, tend, nt)

    yini_array = np.array(yini)
    neq = yini_array.size

    y = np.zeros((neq, nt), order='F')
    y[:,0] = yini_array

    for it, tn  in enumerate(t[:-1]):
        yn = y[:,it]
        k1 = fcn(tn, yn)
        k2 = fcn(tn + c2*dt, yn + dt*(a21*k1))
        k3 = fcn(tn + c3*dt, yn + dt*(a31*k1 + a32*k2))
        k4 = fcn(tn + c4*dt, yn + dt*(a41*k1 + a42*k2 + a43*k3))
        y[:,it+1] = yn + dt*(b1*k1+b2*k2+b3*k3+b4*k4)

    nfev = 4*(nt-1)

    return ode_result(y, t, nfev)

#############################################################################
def rk_embedded(tini, tend, yini, fcn, tol, dt_ini=1.e-2):

    # butcher table of the "3/8 rule" 
    c2 = 1/3
    c3 = 2/3
    c4 = 1

    a21 = 1/3
    a31 = -1/3
    a32 = 1
    a41 = 1
    a42 = -1
    a43 = 1

    b1 = 1/8
    b2 = 3/8
    b3 = 3/8
    b4 = 1/8

    # coefficients bhat of a third order method 
    b1hat = 2*b1 - 1/6 
    b2hat = 2*(1-c2)*b2 
    b3hat = 2*(1-c3)*b3
    b4hat = 0.
    b5hat = 1/6
     
    dt = dt_ini

    yini_array = np.array(yini)
    neq = yini_array.size

    t = [tini]
    y = [yini_array] 

    tn = tini
    yn = yini_array

    k5 = np.array(fcn(tn, yn))
    nfev = 1

    while (tn < tend):

        k1 = k5
        k2 = fcn(tn + c2*dt, yn + dt*a21*k1)
        k3 = fcn(tn + c3*dt, yn + dt*(a31*k1 + a32*k2))
        k4 = fcn(tn + c4*dt, yn + dt*(a41*k1 + a42*k2 + a43*k3))

        ynp1 = yn + dt*(b1*k1 + b2*k2 + b3*k3 + b4*k4)

        k5 = np.array(fcn(tn+dt, ynp1))
        nfev += 4

        yhatnp1 = yn + dt*(b1hat*k1 + b2hat*k2 + b3hat*k3 + b4hat*k4 + b5hat*k5)

        err = np.sqrt(1./neq) * np.linalg.norm((ynp1 - yhatnp1)/(1+np.maximum(np.abs(yn), np.abs(ynp1))))
        dtopt = dt * min(5, max(0.2, 0.9*np.power(tol/err, 1./4.)))

        if (err < tol):

            # update for next step
            tn = tn + dt
            yn = ynp1

            # save time and solution 
            t.append(tn)
            y.append(yn)

            # compute next dt
            dt = min(dtopt, tend-tn)

        else:
            dt = dtopt
            k5 = k1

    t = np.array(t)
    y = np.array(np.transpose(np.array(y)), order='F')

    return ode_result(y, t, nfev)

#############################################################################
class radau_result:
    def __init__(self, y, t, nfev, njev, nstep, naccpt, nrejct, ndec, nsol):
        self.y = y
        self.t = t
        self.nfev = nfev 
        self.njev = njev 
        self.nstep = nstep
        self.naccpt = naccpt
        self.nrejct = nrejct
        self.ndec = ndec
        self.nsol = nsol

#############################################################################
def radau5(tini, tend, yini, fcn, njac, rtol, atol, iout=1):

    c_integration = ct.CDLL("./mylib/lib_radau_rock.so")

    tsol=[]
    ysol=[]

    def solout(nr, told, t, y, cont, lrc, n, rpar, ipar, irtrn):
        ##import ipdb; ipdb.set_trace()
        tsol.append(t[0])
        tmp = []
        for i in range(n[0]):
            tmp.append(y[i])
        ysol.append(np.array(tmp))

    fcn_type = ct.CFUNCTYPE(None, ct.POINTER(ct.c_int), ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),
                            ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.POINTER(ct.c_int))
    solout_type = ct.CFUNCTYPE(None, ct.POINTER(ct.c_int), ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),
                               ct.POINTER(ct.c_double), ct.POINTER(ct.c_double), ct.POINTER(ct.c_int),
                               ct.POINTER(ct.c_int), ct.POINTER(ct.c_double), ct.POINTER(ct.c_int), ct.POINTER(ct.c_int))

    c_radau5 = c_integration.radau5_integration
    c_radau5.argtypes = [ct.c_double, ct.c_double, ct.c_int,
                         np.ctypeslib.ndpointer(dtype = np.float64),
                         np.ctypeslib.ndpointer(dtype = np.float64),
                         fcn_type, solout_type, ct.c_double, ct.c_double, ct.c_int, 
                         ct.c_int, np.ctypeslib.ndpointer(dtype = np.int32)]
    c_radau5.restype = None

    callable_fcn = fcn_type(fcn)
    callable_solout = solout_type(solout)

    yini_array = np.array(yini)
    neq = yini_array.size
    yn = np.zeros(neq)
    info= np.zeros(7, dtype=np.int32)
    c_radau5(tini, tend, neq, yini_array, yn, callable_fcn, callable_solout, rtol, atol, njac, iout, info)

    if iout == 1:
        tsol = np.array(tsol)
        ysol = np.array(np.transpose(np.array(ysol)), order='F')
    else:
        tsol = tend
        ysol = yn

    nfev   = info[0]  # number of function evaluations
    njev   = info[1]  # number of jacobian evaluations
    nstep  = info[2]  # number of computed steps
    naccpt = info[3]  # number of accepted steps
    nrejct = info[4]  # number of rejected steps
    ndec   = info[5]  # number of lu-decompositions
    nsol   = info[6]  # number of forward-backward substitutions

    return radau_result(ysol, tsol, nfev, njev, nstep, naccpt, nrejct, ndec, nsol)

#############################################################################
class rock_result:
    def __init__(self, y, nfev, nstep, naccpt, nrejct, nfevrho, nstage):
        self.y = y
        self.nfev = nfev 
        self.nstep = nstep
        self.naccpt = naccpt
        self.nrejct = nrejct
        self.nfevrho = nfevrho
        self.nstage = nstage

#############################################################################
def rock4(tini, tend, yini, fcn, tol):

    c_integration = ct.CDLL("./mylib/lib_radau_rock.so")

    fcn_type = ct.CFUNCTYPE(None, ct.POINTER(ct.c_int), ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),
                            ct.POINTER(ct.c_double))

    c_rock4 = c_integration.rock4_integration
    c_rock4.argtypes = [ct.c_double, ct.c_double, ct.c_int,
                        np.ctypeslib.ndpointer(dtype=np.float64),
                        np.ctypeslib.ndpointer(dtype=np.float64),
                        fcn_type, ct.c_double, np.ctypeslib.ndpointer(dtype=np.int32)]
    c_rock4.restype = None

    callable_fcn = fcn_type(fcn)

    neq = yini.size
    y = np.zeros(neq)
    info = np.zeros(8, dtype=np.int32)
    c_rock4(tini, tend, neq, yini, y, callable_fcn, tol, info)

    nfev    = info[0]   # number of function evaluations.
    nstep   = info[1]   # number of computed steps
    naccpt  = info[2]   # number of accepted steps
    nrejct  = info[3]   # number of rejected steps
    nfevrho = info[4]   # number of evaluations of f used to estimate the spectral radius
    nstage  = info[5]   # maximum number of stages used

    return rock_result(y, nfev, nstep, naccpt, nrejct, nfevrho, nstage)

#############################################################################
def rock4_dt(tini, tend, yini, fcn, tol):

    c_integration = ct.CDLL("./mylib/lib_radau_rock.so")

    fcn_type = ct.CFUNCTYPE(None, ct.POINTER(ct.c_int), ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),
                            ct.POINTER(ct.c_double))

    c_rock4 = c_integration.rock4_integration_dt
    c_rock4.argtypes = [ct.c_double, ct.c_double, ct.c_int,
                        np.ctypeslib.ndpointer(dtype=np.float64),
                        np.ctypeslib.ndpointer(dtype=np.float64),
                        fcn_type, ct.c_double, np.ctypeslib.ndpointer(dtype=np.float64), np.ctypeslib.ndpointer(dtype=np.int32)]
    c_rock4.restype = None

    callable_fcn = fcn_type(fcn)

    neq = yini.size
    y = np.zeros(neq)
    tsol = np.zeros(10000)
    info = np.zeros(8, dtype=np.int32)
    c_rock4(tini, tend, neq, yini, y, callable_fcn, tol, tsol, info)

    nfev    = info[0]   # number of function evaluations.
    nstep   = info[1]   # number of computed steps
    naccpt  = info[2]   # number of accepted steps
    nrejct  = info[3]   # number of rejected steps
    nfevrho = info[4]   # number of evaluations of f used to estimate the spectral radius
    nstage  = info[5]   # maximum number of stages used

    return (tsol)
