import numpy as np
import ctypes as ct

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
def radau5(tini, tend, yini, fcn, njac, rtol=1.e-6, atol=1.e-6, iout=0):

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
def rock4(tini, tend, yini, fcn, tol=1.e-6):

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
def strang(tini, tend, nt, yini, fcn_diff, fcn_reac):

    t = np.linspace(tini, tend, nt)
    dt = (tend-tini) / (nt-1)
    #nt-1 = (tend-tini) 

    ysol = yini 

    for it, ti in enumerate(t[:-1]):
        #print(ti, ti+dt)
        sol = radau5(ti, ti+dt/2, ysol, fcn_reac, 0)
        ysol = sol.y
        ##print(np.linalg.norm(ysol))
        ##sol = radau5(ti, ti+dt, ysol, fcn_diff, 1)
        sol = rock4(ti, ti+dt, ysol, fcn_diff)
        ysol = sol.y
        sol = radau5(ti+dt/2, ti+dt, ysol, fcn_reac, 0)
        ysol = sol.y
        #print(np.linalg.norm(ysol))

    return ysol
