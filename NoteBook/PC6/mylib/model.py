import numpy as np

class brusselator_model: 
 
    def __init__(self, a, b=3): 
        self.a = a 
        self.b = b 
 
    def fcn(self, t, y):
        y1, y2 = y
        a = self.a
        b = self.b 
        y1_dot = a - (b+1)*y1 + y1*y1*y2 
        y2_dot = b*y1 - y1*y1*y2  
        return np.array([y1_dot, y2_dot])

    def jac(self, t, y):
        y1, y2 = y
        a = self.a
        b = self.b 
        return np.array([[-(b+1) + 2*y1*y2, y1*y1], [b - 2*y1*y2, -y1*y1]])

    def dfoverda(self, y):
        return np.array([1, 0])

    def jac_eq(self):
        a = self.a
        b = self.b 
        return np.array([[b-1, a*a], [-b, -a*a]])

class thermal_explosion_model:

    def __init__(self, fk):
        self.fk = fk

    def fcn(self, t, y):
        fk = self.fk
        y1 = y[0]
        y1_dot = fk * np.exp(y1) - y1
        return np.array([y1_dot])

    def jac(self, t, y) :
        fk = self.fk
        y1 = y[0]
        return np.array([[fk * np.exp(y1) - 1]])

    def dfoverda(self, y):
        return np.exp(y) 

    def dyoverds(self, y):
        fk = self.fk
        tmp = np.exp(-y) - fk
        return np.sqrt(1/(1+tmp*tmp)) 

    def daoverds(self, y):
        fk = self.fk
        tmp = np.exp(-y) - fk
        return tmp * np.sqrt(1/(1+tmp*tmp)) 

    #def jac_eq(self):
    #    fk = self.fk
    #    return np.array([], [-b, -a*a]])

class bead_hoop_model:

    def __init__(self, omega):
        self.omega = omega

    def fcn(self, t, y):
        y1, y2 = y
        omega = self.omega
        omega_c = np.sqrt(9.81)
        y1_dot = y2
        y2_dot = np.sin(y1)*(omega*omega*np.cos(y1) - omega_c*omega_c)
        return np.array([y1_dot, y2_dot])

    def jac(self, t, y):
        y1, y2 = y
        omega = self.omega
        omega_c = np.sqrt(9.81)
        j21 = omega*omega*(np.cos(y1)*np.cos(y1) - np.sin(y1)*np.sin(y1)) - omega_c*omega_c*np.cos(y1)
        return np.array([[0., 1.], [j21 , 0]])

    def dfoverda(self, y):
        y1, y2 = y
        omega = self.omega
        return np.array([0., 2*omega*np.sin(y1)*np.cos(y1)])

    def dyoverds(self, y):
        y1, y2 = y
        omega = self.omega
        omega_c = np.sqrt(9.81)
        return np.array([0., 0.])

    def daoverds(self, y):
        y1, y2 = y
        omega = self.omega
        omega_c = np.sqrt(9.81)
        return np.array([1.])

    def jac_eq_1(self):
        omega = self.omega
        omega_c = np.sqrt(9.81)
        j21 = omega*omega - omega_c*omega_c
        return np.array([[0., 1.], [j21 , 0]])

    def jac_eq_2(self):
        omega = self.omega
        omega_c = np.sqrt(9.81)
        
        if (omega<=omega_c):
            j21 = omega*omega - omega_c*omega_c
            return np.array([[0., 1.], [j21 , 0]])
        else:
            siny1siny1 = 1 - (omega_c*omega_c)/(omega*omega)
            j21 = omega_c*omega_c - omega*omega* siny1siny1 - (omega_c*omega_c*omega_c*omega_c)/(omega*omega)
            return np.array([[0., 1.], [j21 , 0]])


class bead_hoop_friction_model:


    def __init__(self, omega, alpha):
        self.omega = omega
        self.alpha = alpha

    def fcn(self, t, y):
        y1, y2 = y
        omega = self.omega
        alpha = self.alpha
        omega_c = np.sqrt(9.81)
        y1_dot = y2
        y2_dot = np.sin(y1)*(omega*omega*np.cos(y1) - omega_c*omega_c) - alpha*y2
        return np.array([y1_dot, y2_dot])

    def jac_eq_1(self):
        omega = self.omega
        alpha = self.alpha
        omega_c = np.sqrt(9.81)
        j21 = omega*omega - omega_c*omega_c
        return np.array([[0., 1.], [j21 , -alpha]])

    def jac_eq_2(self):
        omega = self.omega
        omega_c = np.sqrt(9.81)
        alpha = self.alpha
        
        if (omega<=omega_c):
            j21 = omega*omega - omega_c*omega_c
            return np.array([[0., 1.], [j21 , -alpha]])
        else:
            siny1siny1 = 1 - (omega_c*omega_c)/(omega*omega)
            j21 = omega_c*omega_c - omega*omega* siny1siny1 - (omega_c*omega_c*omega_c*omega_c)/(omega*omega)
            return np.array([[0., 1.], [j21 , -alpha]])
