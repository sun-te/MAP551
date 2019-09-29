import numpy as np

#################################################################
class double_well_potential_model:

    def __init__(self):
        pass
    
    def fcn(self, t, y):
        q, p = y
        q_dot = self.H_p(p)
        p_dot = -self.H_q(q)        
        return np.array([q_dot, p_dot])

    def H_p(self, p):
        return p
    
    def H_q(self, q):
        return 4*(q*q*q - q)
    
    def hamiltonian(self, y):
        q, p = np.split(y,2)
        return (0.5*p*p) + (q*q-1)*(q*q-1)


#################################################################
class bead_hoop_model :

    def __init__(self, omega):
        self.omega = omega

    def fcn(self, t, y):
        q, p = y
        q_dot = self.H_p(p)
        p_dot = -self.H_q(q)
        return np.array([q_dot, p_dot])

    def H_p(self, p):
        return p
    
    def H_q(self, q):
        omega = self.omega
        omega_c = np.sqrt(9.81)
        return -np.sin(q)*(omega*omega*np.cos(q) - omega_c*omega_c)

    def hamiltonian(self, y):
        q, p = np.split(y,2)
        omega = self.omega 
        oo = omega*omega
        return 0.5*(-oo*np.sin(q)*np.sin(q) + p*p) - 9.81*np.cos(q)

#################################################################
class three_body_model:

    def __init__(self):
        self.m0 = 1.00000597682
        self.m1 = 9.54786104043e-4
        self.m2 = 2.85583733151e-4
        self.G = 2.95912208286e-4

    def fcn(self, t, y):
        q = y[0:9]
        p = y[9:18]
        q_dot = self.H_p(p)
        p_dot = -self.H_q(q)
        return np.concatenate((q_dot, p_dot))

    def H_p(self, p):
        m0 = self.m0
        m1 = self.m1
        m2 = self.m2
        Hp = np.zeros(9)
        Hp[0:3] = p[0:3]/m0
        Hp[3:6] = p[3:6]/m1
        Hp[6:9] = p[6:9]/m2
        return Hp
    
    def H_q(self, q):
        m0 = self.m0
        m1 = self.m1
        m2 = self.m2
        G = self.G
        Hq = np.zeros(9)
        Hq[0:3] = ( G*m0*m1*((q[0:3]-q[3:6])/np.power(np.linalg.norm(q[0:3]-q[3:6]),3)) 
                   +G*m0*m2*((q[0:3]-q[6:9])/np.power(np.linalg.norm(q[0:3]-q[6:9]),3)))
        Hq[3:6] = ( G*m1*m0*((q[3:6]-q[0:3])/np.power(np.linalg.norm(q[3:6]-q[0:3]),3))
                   +G*m1*m2*((q[3:6]-q[6:9])/np.power(np.linalg.norm(q[3:6]-q[6:9]),3)))
        Hq[6:9] = ( G*m2*m0*((q[6:9]-q[0:3])/np.power(np.linalg.norm(q[6:9]-q[0:3]),3)) 
                   +G*m2*m1*((q[6:9]-q[3:6])/np.power(np.linalg.norm(q[6:9]-q[3:6]),3)))
        return Hq

    def hamiltonian(self, y):
    
        m0 = self.m0
        m1 = self.m1
        m2 = self.m2
        G = self.G
        nt = y.shape[1]
        neq = y.shape[0]
        q = y[0:neq//2]
        p = y[neq//2:neq]
        ham = np.zeros(nt)

        for i in range(nt):
            ham[i] = ( 0.5*( (1/m0)*np.dot(p[0:3,i],p[0:3,i]) 
                            +(1/m1)*np.dot(p[3:6,i],p[3:6,i]) 
                            +(1/m2)*np.dot(p[6:9,i],p[6:9,i]) ) 
                      - G *( (m0*m1)/np.linalg.norm(q[0:3,i]-q[3:6,i]) 
                            +(m0*m2)/np.linalg.norm(q[0:3,i]-q[6:9,i])
                            +(m1*m2)/np.linalg.norm(q[3:6,i]-q[6:9,i]) ) )

        return ham

#################################################################
class arenstorf_model: 
 
    def __init__(self, mu): 
        self.mu = mu 
 
    def fcn(self, t, y) : 
        y1,y2,y3,y4 = y 
        mu = self.mu 
        r1 = np.sqrt((y1+mu)*(y1+mu) + y2*y2) 
        r2 = np.sqrt((y1-1+mu)*(y1-1+mu) + y2*y2) 
        y1_dot = y3 
        y2_dot = y4 
        y3_dot = y1 + 2*y4 - (1-mu)*(y1+mu)/(r1*r1*r1) - mu*(y1 - 1 + mu)/(r2*r2*r2) 
        y4_dot = y2 - 2*y3 - (1-mu)*y2/(r1*r1*r1) - mu*y2/(r2*r2*r2) 
        return np.array([y1_dot, y2_dot, y3_dot, y4_dot])

    def V_q(self, q): 
        mu = self.mu 
        q1 = q[0]  
        q2 = q[1]  
        
        # not rhe same r1, r2 than fcn
        r1 = np.sqrt((q1-1+mu)*(q1-1+mu) + q2*q2) 
        r2 = np.sqrt((q1+mu)*(q1+mu) + q2*q2) 
         
        ##Hq1 = q1 + 2*p2 - (1-mu)*(q1+mu)/(r1*r1*r1) - mu*(q1 - 1 + mu)/(r2*r2*r2) 
        ##Hq2 = q2 - 2*p1 - (1-mu)*q2/(r1*r1*r1) - mu*q2/(r2*r2*r2) 
        V_q1 = -q1 + (1-mu)*(q1+mu)/(r2*r2*r2) + mu*(q1-1+ mu)/(r1*r1*r1) 
        V_q2 = -q2 + (1-mu)*q2/(r2*r2*r2) + mu*q2/(r1*r1*r1) 
 
        return np.array((V_q1, V_q2)) 
 
    def Fp(self, t, p): 
        p1 = p[0]  
        p2 = p[1] 
 
        w = 2. 
         
        Fp1 = t*p1 + ((1-np.cos(w*t))/(w*w))*2*p2 + ((np.sin(w*t) - w*t)/(w*w*w))*4*p1 
        Fp2 = t*p2 - ((1-np.cos(w*t))/(w*w))*2*p1 + ((np.sin(w*t) - w*t)/(w*w*w))*4*p2 
         
        return np.array((Fp1, Fp2)) 
 
    def expp(self, t, p): 
        p1 = p[0]  
        p2 = p[1]  
 
        w = 2. 
         
        expp1 = p1 + (np.sin(w*t)/w)*2*p2 - 2*(np.sin(0.5*w*t)/w)*(np.sin(0.5*w*t)/w)*4*p1 
        expp2 = p2 - (np.sin(w*t)/w)*2*p1 - 2*(np.sin(0.5*w*t)/w)*(np.sin(0.5*w*t)/w)*4*p2 
         
        return np.array((expp1, expp2)) 

    def get_hamiltonian(self, y):
    
        mu = self.mu
        nt = y.shape[0]
        neq = y.shape[1]
        q = y[:,0:neq//2]
        p = y[:,neq//2:neq]
        ham = np.zeros(nt)

        for i in range(nt):
            
            q1 = q[i,0]
            q2 = q[i,0]
            r1 = np.sqrt((q1-1+mu)*(q1-1+mu) + q2*q2) 
            r2 = np.sqrt((q1+mu)*(q1+mu) + q2*q2) 
            ham[i] = 0.5*np.dot(p[i],p[i]) - 0.5*(q1*q1+q2*q2) - mu/r1 - (1-mu)/r2

        return ham
