import numpy as np

def natural_continuation_order_0(a_ini, a_end, na, yeq_ini, model, b,tol=1.e-8, nmax_newton=100):
    
    a = np.linspace(a_ini, a_end, na)
        
    yeq = np.zeros((len(yeq_ini), na), order='F')
    yeq[:,0] = yeq_ini
        
    for it, ai in enumerate(a[1:]):
        
        mod = model(ai,b)
        fcn = mod.fcn 
        jac = mod.jac
        
        # Newton iteration
        yn = yeq[:, it]
        #print("++++++++++++++++++++++++++++++++++++++")
        #print("Initialization ||f(y)|| = ",  np.linalg.norm(fcn(0, yn)))
        for it_newton in range(nmax_newton):
            yn += np.linalg.solve(jac(0, yn), -fcn(0, yn))
            #print(yn)
            ##print("Iter nb ", it_newton+1, " ||f(y)|| = ",  np.linalg.norm(fcn(0, yn)))
            if (np.linalg.norm(fcn(0, yn)) < tol):
                #print(it_newton)
                break
        yeq[:, it+1] = yn

    return(yeq)

def natural_continuation_order_1(a_ini, a_end, na, yeq_ini, model,b, tol=1.e-8, nmax_newton=100):
    
    a = np.linspace(a_ini, a_end, na)
    da = (a_end - a_ini)/(na-1)
        
    yeq = np.zeros((len(yeq_ini), na), order='F')
    yeq[:,0] = yeq_ini 

    for it, ai in enumerate(a[1:]):
        
        mod = model(ai,b)
        fcn = mod.fcn 
        jac = mod.jac
        dfoverda = mod.dfoverda
        
        # Newton iteration
        yn = yeq[:, it] + np.linalg.solve(jac(0, yeq_ini), -dfoverda(yeq[:, it]))*da
        #print("++++++++++++++++++++++++++++++++++++++")
        #print("Initialization ||f(y)|| = ",  np.linalg.norm(fcn(0, yn)))
        for it_newton in range(nmax_newton):
            yn += np.linalg.solve(jac(0, yn), -fcn(0, yn))
            ##print("Iter nb ", it_newton+1, " ||f(y)|| = ",  np.linalg.norm(fcn(0, yn)))
            if (np.linalg.norm(fcn(0, yn)) < tol):
                #print(it_newton)
                break
        yeq[:, it+1] = yn

    return(yeq)


def pseudo_arclength_continuation(s_ini, s_end, ns, a_ini, yeq_ini, model, tol=1.e-8, nmax_newton=10):

    ds = (s_end - s_ini) / (ns-1)

    a = np.zeros(ns)

    ny = len(yeq_ini)     
    yeq = np.zeros((ny, ns), order='F')
    yeq[:,0] = yeq_ini 
    
    for it in range(ns-1):

        ##print("Continuation it ", it) 

        mod = model(a[it])
        fcn = mod.fcn 
        jac = mod.jac
        dyoverds = mod.dyoverds
        daoverds = mod.daoverds
        dfoverda = mod.dfoverda

        yn = np.array(yeq[:, it])
        an = a[it]

        dyoverdsn = dyoverds(yn)
        daoverdsn = daoverds(yn)

        for it_newton in range(nmax_newton):

            ##print("  Newton it ", it_newton) 

            jac_aug = np.block([ [ jac(0, yn), dfoverda(yn).reshape(ny,1)],
                                 [ dyoverdsn,  daoverdsn    ]])
            #print("J_aug ") 
            #print(jac_aug) 

            n = np.dot(yn - yeq[:, it], dyoverdsn) + np.dot(an - a[it], daoverdsn) - ds  
            rhs = np.block([fcn(0, yn), n])
            #print("rhs")
            #print(-rhs)

            y_aug = np.linalg.solve(jac_aug, -rhs)
            
            yn += y_aug[:-1]
            an += y_aug[-1]
            #print("y_aug")
            print(y_aug)

            mod = model(an)
            fcn = mod.fcn 
            jac = mod.jac

            if (np.linalg.norm(y_aug) < tol):
                break
            
        yeq[:, it+1] = np.array(yn)
        a[it+1] = an

    return a, yeq




def pseudo_arclength_continuation_bis(s_ini, s_end, ns, a_ini, yeq_ini, model, tol=1.e-8, nmax_newton=1):

    ds = (s_end - s_ini) / (ns-1)

    a = np.zeros(ns)

    ny = len(yeq_ini)     
    yeq = np.zeros((ny, ns), order='F')
    yeq[:,0] = yeq_ini 
    
    for it in range(ns-1):

        print("Continuation it ", it) 

        mod = model(a[it])
        fcn = mod.fcn 
        jac = mod.jac
        dyoverds = mod.dyoverds
        daoverds = mod.daoverds
        dfoverda = mod.dfoverda

        yn = np.array(yeq[:, it])
        an = a[it]

        dyoverdsn = dyoverds(yn)
        daoverdsn = daoverds(yn)

        for it_newton in range(nmax_newton):

            print("  Newton it ", it_newton) 

            jac_aug = np.block([ [ jac(0, yn), dfoverda(yn).reshape(ny,1)],
                                 [ dyoverdsn,  daoverdsn    ]])
            print("J_aug ") 
            print(jac_aug) 

            #print("dyoverds :", dyoverdsn) 
            #print("daoverds :", daoverdsn) 
            #print("ds :", ds) 
            #n = (yn - yeq[:, it])*dyoverdsn + (an - a[it])*daoverdsn - ds  
            n = np.dot(yn - yeq[:, it], dyoverdsn) + np.dot(an - a[it], daoverdsn) - ds  
            #print("n") 
            #print(n) 
            rhs = np.block([fcn(0, yn), n])
            print("rhs")
            print(-rhs)

            y_aug = np.linalg.solve(jac_aug, -rhs)
            print("y_aug")
            print(y_aug)
            yn += y_aug[:-1]
            an += y_aug[-1]

            #print("fcn", np.linalg.norm(fcn(0,yn)))
            mod = model(an)
            fcn = mod.fcn 
            jac = mod.jac
            #n = np.dot(yn - yeq[:, it], dyoverdsn) + np.dot(an - a[it], daoverdsn) - ds  
            #rhs = np.block([fcn(0, yn), n])
            #print(np.linalg.norm(fcn(0,yn)))
            #print(np.linalg.norm(n))
            #print("rhs", np.linalg.norm(rhs))
            

            if (np.linalg.norm(rhs) < tol):
                break
            if (np.linalg.norm(y_aug) < tol):
                break
            
        yeq[:, it+1] = np.array(yn)
        a[it+1] = an
        print("yeq ", yeq[:, it+1]) 
        print("an", an) 

    return a, yeq



def pseudo_arclength_continuation_for_hoop(s_ini, s_end, ns, omega_ini, y1_ini, tol=1.e-6):

    ds = (s_end - s_ini) / (ns-1)
    
    y1 = np.zeros((ns, 3))
    omega = np.zeros((ns,3))
    
    y1[0] = y1_ini 
    omega[0] = omega_ini 

    #print("y1_ini", y1_ini) 
    #print("omega_ini", omega_ini) 
    
    omega_c = np.sqrt(9.81)
    
    for it_cont in range(ns-1):
       
        ##print("it cont :", it_cont)

        # Unstable branch          
        y1_n = y1[it_cont,0]
        omega_n = omega[it_cont,0]
        dy1onds = 0.
        domegaonds = 1.
        for it_newton in range(101):
         
            R_n = np.sin(y1_n)*(omega_n*omega_n*np.cos(y1_n) - omega_c*omega_c)
            n_n = (y1_n - y1[it_cont,0])*dy1onds + (omega_n - omega[it_cont,0])*domegaonds - ds
            rhs = np.array([-R_n, -n_n])    

            dRondy1 = (omega_n*omega_n)*np.cos(2*y1_n) - omega_c*omega_c*np.cos(y1_n)
            dRondomega = 2*omega_n*np.sin(y1_n)*np.cos(y1_n)
            J_aug = np.array([ [dRondy1, dRondomega], 
                               [dy1onds, domegaonds] ])
           
            y_aug = np.linalg.solve(J_aug, rhs)
            y1_n = y1_n + y_aug[0]
            omega_n = omega_n + y_aug[1]
            
            if (np.linalg.norm(y_aug) < tol):
                break

        y1[it_cont+1,0] = y1_n
        omega[it_cont+1,0] = omega_n
        
        # Stable branch 1
        y1_n = y1[it_cont,1]
        omega_n = omega[it_cont,1]
        if omega[it_cont,1] > omega_c: 
            tmp = (omega[it_cont,1]/2)*np.tan(y1[it_cont,1])
            dy1onds = np.sqrt(1/(1+tmp*tmp))
            domegaonds = tmp*dy1onds

        for it_newton in range(101):
            
            if omega[it_cont,1] < omega_c:
                R_n = np.sin(y1_n)*(omega_n*omega_n*np.cos(y1_n) - omega_c*omega_c)
                dRondy1 = (omega_n*omega_n)*np.cos(2*y1_n) - omega_c*omega_c*np.cos(y1_n)
                dRondomega = 2*omega_n*np.sin(y1_n)*np.cos(y1_n)
            else: 
                R_n = (omega_n*omega_n*np.cos(y1_n) - omega_c*omega_c)
                dRondy1 = -(omega_n*omega_n)*np.sin(y1_n)
                dRondomega = 2*omega_n*np.cos(y1_n)

            n_n = (y1_n - y1[it_cont,1])*dy1onds + (omega_n - omega[it_cont,1])*domegaonds - ds
            rhs = np.array([-R_n, -n_n])    

            J_aug = np.array([ [dRondy1, dRondomega], 
                               [dy1onds, domegaonds] ])
           
            y_aug = np.linalg.solve(J_aug, rhs)
            y1_n = y1_n + y_aug[0]
            omega_n = omega_n + y_aug[1]

            if (np.linalg.norm(y_aug) < tol):
                break

        y1[it_cont+1,1] = y1_n
        omega[it_cont+1,1] = omega_n
        

        # Stable branch 2
        y1_n = y1[it_cont,2]
        omega_n = omega[it_cont,2]
        if omega[it_cont,1] > omega_c: 
            tmp = (omega[it_cont,2]/2)*np.tan(y1[it_cont,2])
            dy1onds = -np.sqrt(1/(1+tmp*tmp))
            domegaonds = tmp*dy1onds

        for it_newton in range(101):
            
            if omega[it_cont,2] < omega_c:
                R_n = np.sin(y1_n)*(omega_n*omega_n*np.cos(y1_n) - omega_c*omega_c)
                dRondy1 = (omega_n*omega_n)*np.cos(2*y1_n) - omega_c*omega_c*np.cos(y1_n)
                dRondomega = 2*omega_n*np.sin(y1_n)*np.cos(y1_n)
            else: 
                R_n = (omega_n*omega_n*np.cos(y1_n) - omega_c*omega_c)
                dRondy1 = -(omega_n*omega_n)*np.sin(y1_n)
                dRondomega = 2*omega_n*np.cos(y1_n)

            n_n = (y1_n - y1[it_cont,2])*dy1onds + (omega_n - omega[it_cont,2])*domegaonds - ds
            rhs = np.array([-R_n, -n_n])    

            J_aug = np.array([ [dRondy1, dRondomega], 
                               [dy1onds, domegaonds] ])
           
            y_aug = np.linalg.solve(J_aug, rhs)
            y1_n = y1_n + y_aug[0]
            omega_n = omega_n + y_aug[1]

            if (np.linalg.norm(y_aug) < tol):
                break

        y1[it_cont+1,2] = y1_n
        omega[it_cont+1,2] = omega_n
         
    return omega, y1
