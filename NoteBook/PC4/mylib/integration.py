import numpy as np
from scipy.optimize import fsolve, root
from scipy.integrate import solve_ivp, ode

import copy

def backward_euler(tini, tend, nt, uini, fcn): 
 
    dt = (tend-tini) / (nt-1) 
 
    uini_array = np.array(uini)
    neq = uini_array.size 
 
    u = np.zeros((nt, neq))
    u[0] = uini_array 
 
    def g(uip1, *args): 
        uip, tip1 = args 
        return uip1 - uip - dt*np.array(fcn(tip1, uip1))
 
    for it in range(nt-1) : 
        # solve u[it+1] - u[it] - dt * fcn(tini + (it+1)*dt, u[it+1]) = 0 
        u0 = u[it] + dt * np.array(fcn(tini+ it*dt, u[it])) 
        u[it+1] = fsolve(g, u0, (u[it], tini + (it+1)*dt)) 
         
    return u 

def forward_euler(tini, tend, nt, uini, fcn):

    dt = (tend-tini) / (nt-1)

    uini_array = np.array(uini)
    neq = uini_array.size

    u = np.zeros((nt, neq))
    u[0] = uini_array

    for it in range(nt-1) :
        u[it+1] = u[it] + dt * np.array(fcn(tini+ it*dt, u[it]))

    return u

def rk2(tini, tend, nt, uini, fcn):

    dt = (tend-tini) / (nt-1)

    uini_array = np.array(uini)
    neq = uini_array.size

    u = np.zeros((nt, neq))
    u[0] = uini_array

    for it in range(nt-1):

        un = u[it]
        tn = tini + it*dt 

        k1 = np.array(fcn(tn, un))
        k2 = np.array(fcn(tn + 0.5*dt, un + dt*(0.5*k1)))

        u[it+1] = u[it] + dt*k2

    return u

def rk3(tini, tend, nt, uini, fcn):

    dt = (tend-tini) / (nt-1)

    uini_array = np.array(uini)
    neq = uini_array.size

    u = np.zeros((nt, neq))
    u[0] = uini_array

    for it in range(nt-1):

        un = u[it]
        tn = tini + it*dt 

        k1 = np.array(fcn(tn, un))
        k2 = np.array(fcn(tn + 0.5*dt, un + dt*(0.5*k1)))
        k3 = np.array(fcn(tn + dt, un + dt*(-k1 + 2*k2)))

        u[it+1] = u[it] + (dt/6)*(k1+4*k2+k3)

    return u

def rk4(tini, tend, nt, uini, fcn):

    dt = (tend-tini) / (nt-1)

    uini_array = np.array(uini)
    neq = uini_array.size

    u = np.zeros((nt, neq))
    u[0] = uini_array

    for it in range(nt-1):

        un = u[it]
        tn = tini + it*dt 

        k1 = np.array(fcn(tn, un))
        k2 = np.array(fcn(tn + 0.5*dt, un + dt*(0.5*k1)))
        k3 = np.array(fcn(tn + 0.5*dt, un + dt*(0.5*k2)))
        k4 = np.array(fcn(tn + dt, un + dt*k3))

        u[it+1] = u[it] + (dt/6)*(k1+2*k2+2*k3+k4)

    return u

def rk38(tini, tend, nt, uini, fcn):

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

    uini_array = np.array(uini)
    neq = uini_array.size

    u = np.zeros((nt, neq))
    u[0] = uini_array

    for it in range(nt-1):

        un = u[it]
        tn = tini + it*dt 

        k1 = np.array(fcn(tn, un))
        k2 = np.array(fcn(tn + c2*dt, un + dt*a21*k1))
        k3 = np.array(fcn(tn + c3*dt, un + dt*(a31*k1 + a32*k2)))
        k4 = np.array(fcn(tn + c4*dt, un + dt*(a41*k1 + a42*k2 + a43*k3)))

        u[it+1] = u[it] + dt*(b1*k1+b2*k2+b3*k3+b4*k4)

    return u

class embedded_sol(object):

    def __init__(self, y, t, t_rej, dt, dt_rej, loc_err, loc_err_exa, glob_err) :
        self.y = y
        self.t = t
        self.t_rej = t_rej
        self.dt = dt
        self.dt_rej = dt_rej
        self.loc_err = loc_err
        self.loc_err_exa = loc_err_exa
        self.glob_err = glob_err 

def rk43_embedded(tini, tend, nt, yini, fcn, tol):

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

    b1hat = 2*b1 - 1/6 
    b2hat = 2*(1-c2)*b2 
    b3hat = 2*(1-c3)*b3
    b4hat = 0.
    b5hat = 1/6
     
    dt = (tend-tini) / (nt-1)

    yini_array = np.array(yini)
    neq = yini_array.size

    t_accepted = [tini]
    t_rejected = []
    dt_accepted = []
    dt_rejected = []
    y = [yini_array] 
    loc_err = [0.]
    loc_err_exa = [0.]

    tn = tini
    yn = yini_array

    while (tn < tend):

        k1 = np.array(fcn(tn, yn))
        k2 = np.array(fcn(tn + c2*dt, yn + dt*a21*k1))
        k3 = np.array(fcn(tn + c3*dt, yn + dt*(a31*k1 + a32*k2)))
        k4 = np.array(fcn(tn + c4*dt, yn + dt*(a41*k1 + a42*k2 + a43*k3)))

        ynp1 = yn + dt*(b1*k1 + b2*k2 + b3*k3 + b4*k4)
        k5 = np.array(fcn(tn+dt, ynp1))

        yhatnp1 = yn + dt*(b1hat*k1 + b2hat*k2 + b3hat*k3 + b4hat*k4 + b5hat*k5)

        err = np.sqrt(1./neq) * np.linalg.norm((ynp1 - yhatnp1)/(1+np.maximum(np.abs(yn), np.abs(ynp1))))
        dtopt = dt * min(5, max(0.2, 0.9*np.power(tol/err, 1./4.)))

        if (err < tol):

            # compute and save exact local error
            sol = solve_ivp(fcn, (tn, tn+dt), yn, rtol=1.e-12, atol=1.e-12)
            loc_err_exa.append(np.sqrt(1./neq) * np.linalg.norm(ynp1 - sol.y[:,-1]))
            
            # update for next step
            tn = tn + dt
            yn = ynp1

            # save time, solution and estimate local error
            t_accepted.append(tn)
            dt_accepted.append(dt) 
            y.append(yn)
            loc_err.append(err)

            # compute next dt
            dt = min(dtopt, tend-tn)

        else:
            t_rejected.append(tn)
            dt_rejected.append(dt) 
            dt = dtopt


    y = np.array(y)
    sol = solve_ivp(fcn, (tini, tend), yini, rtol=1.e-12, atol=1.e-12, t_eval=t_accepted)
    glob_err = np.empty(sol.t.size)
    for i in range(sol.t.size):
        glob_err[i] = np.linalg.norm(y[i] - sol.y[:,i])
        
    return embedded_sol(y, np.array(t_accepted), np.array(t_rejected), np.array(dt_accepted), np.array(dt_rejected),
                        np.array(loc_err), np.array(loc_err_exa), glob_err)

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

