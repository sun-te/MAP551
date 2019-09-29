import numpy as np
from scipy.integrate import solve_ivp, ode
from scipy.optimize import fsolve
import copy

class ode_result:
    def __init__(self, y, t, nfev):
        self.y = y
        self.t = t
        self.nfev = nfev 

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

    y = np.zeros((neq, nt), order='FORTRAN')
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

class embedded_ode_result:

    def __init__(self, y, t, nfev, t_rej, dt, dt_rej, loc_err_est, loc_err):
        self.y = y
        self.t = t
        self.nfev = nfev 
        self.t_rej = t_rej
        self.dt = dt
        self.dt_rej = dt_rej
        self.loc_err_est = loc_err_est
        self.loc_err = loc_err

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

    t_acc = [tini]
    t_rej = []
    dt_acc = []
    dt_rej = []
    y = [yini_array] 
    loc_err_est = [0.]
    loc_err = [0.]

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

            # compute and save exact local error
            sol = solve_ivp(fcn, (tn, tn+dt), yn, method="RK45", rtol=1.e-12, atol=1.e-12)
            loc_err.append(np.sqrt(1./neq) * np.linalg.norm(ynp1 - sol.y[:,-1]))
            
            # update for next step
            tn = tn + dt
            yn = ynp1

            # save time, solution and estimate local error
            t_acc.append(tn)
            dt_acc.append(dt) 
            y.append(yn)
            loc_err_est.append(err)

            # compute next dt
            dt = min(dtopt, tend-tn)

        else:
            t_rej.append(tn)
            dt_rej.append(dt) 
            dt = dtopt
            k5 = k1

    y = np.array(np.transpose(np.array(y)), order='F')

    return embedded_ode_result(y, np.array(t_acc), nfev,  np.array(t_rej), np.array(dt_acc), np.array(dt_rej),
                               np.array(loc_err_est), np.array(loc_err))

def dopri853(tini, tend, uini, fcn, tol):

    time=[]
    usol=[]

    def solout(t, y):
        #print("  t = {}, u = {}".format(t, y))
        time.append(t)
        usol.append(copy.copy(y))

    solver = ode(fcn)
    solver.set_integrator('dop853', atol=tol, rtol=tol, nsteps=10000)
    solver.set_initial_value(uini, tini)
    solver.set_solout(solout)

    solver.integrate(tend)

    return np.array(time), np.array(usol)

def quasi_exact(fcn, yini, t_eval):

    solver = ode(fcn)
    solver.set_integrator('dop853', atol=1.e-12, rtol=1.e-12, nsteps=10000)
    solver.set_initial_value(yini, t_eval[0])

    y=[yini]

    for t in t_eval[1::]:
       solver.integrate(t)
       y.append(solver.y)

    return np.array(y)
