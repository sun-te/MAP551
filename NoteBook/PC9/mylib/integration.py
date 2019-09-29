import numpy as np
import copy
from scipy.optimize import fsolve
from scipy.integrate import ode

#################################################################
class ode_result:
    def __init__(self, y, t):
        self.y = y
        self.t = t

#################################################################
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

    return ode_result(y, t)

#################################################################
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

    return ode_result(y, t)

#################################################################
def symplectic_euler(tini, tend, nt, yini, H_p, H_q):
    
    dt = (tend-tini) / (nt-1) 
    t = np.linspace(tini, tend, nt)

    yini_array = np.array(yini)
    neq = yini_array.size

    y = np.zeros((neq, nt), order='FORTRAN')
    y[:,0] = yini_array

    q, p = np.split(y, 2)
    for it in range(nt-1):
        q_n = q[:,it]
        p_n = p[:,it]
        q[:,it+1] = q_n + dt*H_p(p_n)
        p[:,it+1] = p_n - dt*H_q(q[:,it+1])

    return ode_result(y,t)

#################################################################
def stormer_verlet(tini, tend, nt, yini, H_p, H_q):
    
    dt = (tend-tini) / (nt-1) 
    t = np.linspace(tini, tend, nt)

    yini_array = np.array(yini)
    neq = yini_array.size

    y = np.zeros((neq, nt), order='FORTRAN')
    y[:,0] = yini_array
    
    q, p = np.split(y, 2)
    for it in range(nt-1):
        q_n = q[:,it]
        p_n = p[:,it]
        p_np05 = p_n - (dt/2)*H_q(q_n)
        q[:,it+1] = q_n + dt*H_p(p_np05)
        p[:,it+1] = p_np05 - (dt/2)*H_q(q[:,it+1])
        
    return ode_result(y, t)

#################################################################
def stormer_verlet_step(dt, yini, H_p, H_q):
   
    q_n, p_n = np.split(yini,2)

    p_np05 = p_n - (dt/2)*H_q(q_n)
    q_np1  = q_n + dt*H_p(p_np05)
    p_np1  = p_np05 - (dt/2)*H_q(q_np1)

    y = np.concatenate((q_np1, p_np1))
 
    return y

#################################################################
def optimized_815(tini, tend, nt, yini, H_p, H_q):
    
    dt = (tend-tini) / (nt-1) 
    t = np.linspace(tini, tend, nt)

    yini_array = np.array(yini)
    neq = yini_array.size

    y = np.zeros((neq, nt), order='FORTRAN')
    y[:,0] = yini_array

    nstep = 15
    gamma = np.zeros(nstep+1)
    gamma[0]  =  0.
    gamma[1]  =  0.74167036435061295344822780
    gamma[2]  = -0.40910082580003159399730010
    gamma[3]  =  0.19075471029623837995387626
    gamma[4]  = -0.57386247111608226665638773
    gamma[5]  =  0.29906418130365592384446354
    gamma[6]  =  0.33462491824529818378495798
    gamma[7]  =  0.31529309239676659663205666
    gamma[8]  = -0.79688793935291635401978884
    gamma[9]  = gamma[7]
    gamma[10] = gamma[6]
    gamma[11] = gamma[5]
    gamma[12] = gamma[4]
    gamma[13] = gamma[3]
    gamma[14] = gamma[2]
    gamma[15] = gamma[1]
                          
    ytmp = y[:,0]

    for it in range(nt-1):
 
        for istep in range(nstep):

            ti = 0.0
            te = gamma[istep+1]*dt
            ytmp = stormer_verlet_step(te-ti, ytmp, H_p, H_q)

        y[:,it+1] = ytmp
        
    return ode_result(y, t)

#################################################################
def optimized_815_scov(tini, tend, nt, yini, V_q, Fp, expp): 
    
    dt = (tend-tini) / (nt-1) 

    yini_array = np.array(yini)
    neq = yini_array.size

    y = np.zeros((nt, neq))
    y[0] = yini_array

    nstep = 15
    gamma = np.zeros(nstep+1)
    gamma[0]  =  0.
    gamma[1]  =  0.74167036435061295344822780
    gamma[2]  = -0.40910082580003159399730010
    gamma[3]  =  0.19075471029623837995387626
    gamma[4]  = -0.57386247111608226665638773
    gamma[5]  =  0.29906418130365592384446354
    gamma[6]  =  0.33462491824529818378495798
    gamma[7]  =  0.31529309239676659663205666
    gamma[8]  = -0.79688793935291635401978884
    gamma[9]  = gamma[7]
    gamma[10] = gamma[6]
    gamma[11] = gamma[5]
    gamma[12] = gamma[4]
    gamma[13] = gamma[3]
    gamma[14] = gamma[2]
    gamma[15] = gamma[1]
                          
    ytmp = y[0]
    for it in range(nt-1):
 
        for istep in range(nstep):

            ti = 0.0
            te = gamma[istep+1]*dt
            ytmp = scov(ti, te, 2, ytmp, V_q, Fp, expp)

        y[it+1] = ytmp

    return y

#################################################################
def scov(tini, tend, nt, yini, V_q, Fp, expp):

    dt = (tend-tini) / (nt-1)

    yini_array = np.array(yini)
    neq = yini_array.size

    y = yini_array

    for it in range(nt-1):
        q_n = y[0:neq//2]
        p_n = y[neq//2:neq]

        p_np05 = p_n - (dt/2)*V_q(q_n)
        q_np1  = q_n + Fp(dt, p_np05)
        p_np1  = expp(dt, p_np05) - (dt/2)*V_q(q_np1)

        y[0:neq//2] = q_np1
        y[neq//2:neq] = p_np1

    return y

#################################################################
def scovel(tini, tend, nt, yini, V_q, Fp, expp):
    
    dt = (tend-tini) / (nt-1) 

    yini_array = np.array(yini)
    neq = yini_array.size

    y = np.zeros((nt, neq))
    y[0] = yini_array

    for it in range(nt-1):
        q_n = y[it, 0:neq//2]
        p_n = y[it, neq//2:neq]

        p_np05 = p_n - (dt/2)*V_q(q_n)
        q_np1  = q_n + Fp(dt, p_np05)
        p_np1  = expp(dt, p_np05) - (dt/2)*V_q(q_np1)

        y[it+1, 0:neq//2] = q_np1
        y[it+1, neq//2:neq] = p_np1

    return y

#################################################################
def dopri853(tini, tend, nt, yini, fcn, tol=1.e-6):

    time=[]
    ysol=[] 

    def solout(t, y): 
        time.append(t)
        ysol.append(copy.copy(y))
 
    solver = ode(fcn)
    solver.set_integrator('dop853', atol=tol, rtol=tol, nsteps=10000)
    solver.set_initial_value(yini, tini)
    solver.set_solout(solout)
    solver.integrate(tend)


    return ode_result(np.transpose(np.array(ysol)), np.array(time))
