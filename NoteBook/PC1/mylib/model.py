import numpy as np

class first_model(object) : 
 
    def __init__(self) :
        pass
        
    def fcn(self, t, theta) : 
        theta_dot = np.exp(theta)
        return theta_dot
    def Y_t(self, t, Y) : 
        B  =1
        beta=30
        T0 = 400.
        Tb = 3000.
        dY =-B*np.exp(-(beta*T0)/(Tb-(Tb-T0)*Y)) * Y
        return dY
    
    
class second_model(object) :

    def __init__(self, alpha) :
        self.alpha = alpha

    def fcn(self, t, theta) :
        alpha = self.alpha
        theta_dot = np.exp(theta) - alpha*theta 
        return theta_dot

class second_bis_model(object) :

    def __init__(self, alpha, Tr, beta) :
        self.alpha = alpha
        self.Tr = Tr
        self.beta = beta 

    def fcn(self, t, y) :
        theta, Y = y
        alpha = self.alpha
        beta = self.beta
        Tr = self.Tr
        theta_dot = np.exp(theta/(1+(theta/beta))) * Y - alpha*theta 
        Y_dot     = (-1/Tr) * np.exp(theta/(1+(theta/beta))) * Y
        return (theta_dot, Y_dot)

class second_bis_fk_model(object) :

    def __init__(self, alpha, Tr) :
        self.alpha = alpha
        self.Tr = Tr

    def fcn(self, t, y) :
        theta, Y = y
        alpha = self.alpha
        Tr = self.Tr
        theta_dot = np.exp(theta) * Y - alpha*theta 
        Y_dot     = (-1./Tr) * np.exp(theta) * Y
        return (theta_dot, Y_dot)

class third_model(object) :

    def __init__(self, alpha0, mu, a, thetac) :
        self.alpha0 = alpha0
        self.mu = mu
        self.a = a
        self.thetac = thetac

    def fcn(self, t, y) :
        theta, psi = y
        alpha0 = self.alpha0
        mu = self.mu
        a = self.a
        thetac = self.thetac
        theta_dot = np.exp(theta) - alpha0*(1 + mu*psi*psi) * theta
        psi_dot = -a*psi*(psi*psi+thetac-theta)
        return (theta_dot, psi_dot)
