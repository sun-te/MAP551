import numpy as np

from scipy.integrate import solve_ivp

from bokeh.io import push_notebook, show, output_notebook
from bokeh.plotting import figure
from bokeh.layouts import column
from bokeh.models import PrintfTickFormatter

from ipywidgets import interact, IntSlider, Dropdown, FloatSlider, Text

from model import brusselator_model
import integration as integration

def plot_embedded_cost():

    yini = (1.5 , 3.)
    tini = 0.
    tend = 20.
    
    bm = brusselator_model(a=1, b=3)
    fcn = bm.fcn
            
    l_nt = [101, 501, 1001, 5001, 10001, 20001]
    fe_rk38=[]
    norm_rk38=[]
    
    #for nt in l_nt:
    #    sol = solve_ivp(fcn, (tini, tend), yini, rtol=1.e-12, atol=1.e-12, t_eval=np.linspace(tini, tend, nt))
    #    yexa = np.transpose(sol.y)
    #    yrk38 = integration.rk38(tini, tend, nt, yini, fcn)
    #    fe_rk38.append(4*(nt-1))
    #    norm_rk38.append(np.linalg.norm(yexa-yrk38) / np.sqrt(nt))
        
    #l_tol = [1.e-2, 1.e-4, 1.e-6, 1.e-8, 1.e-10, 1.e-12]
    #l_tol = [1.e-10, 1.e-12]
    l_tol = [1.e-12]

    fe_rk43_emb=[]
    norm_rk43_emb=[]
    
    for tol in l_tol:
        sol_rk43_emb = integration.rk43_embedded(tini, tend, 1000, yini, fcn, tol)
        print(sol_rk43_emb.t.size)  
        print(sol_rk43_emb.dt)
        sol = solve_ivp(fcn, (tini, tend), yini, rtol=1.e-12, atol=1.e-12, t_eval=sol_rk43_emb.t)
        fe_rk43_emb.append(4*sol_rk43_emb.dt.size + 4*sol_rk43_emb.dt_rej.size)
        err = 0 
        for i in range(sol_rk43_emb.dt.size):
            ldt = sol_rk43_emb.dt[i]
            err += ldt * ((sol.y[0,i+1]-sol_rk43_emb.y[i+1,0])*(sol.y[0,i+1]-sol_rk43_emb.y[i+1,0]) + \
                          (sol.y[1,i+1]-sol_rk43_emb.y[i+1,1])*(sol.y[1,i+1]-sol_rk43_emb.y[i+1,1])  )
        err = np.sqrt(err)    
        norm_rk43_emb.append(err)
    

    fig = figure(x_axis_type="log", y_axis_type="log", plot_height=500, plot_width=900, \
                 title= "Number of function evaluations vs norm of global error" )
    fig.x(norm_rk38, fe_rk38, legend="rk38", color="Indigo")
    fig.line(norm_rk38, fe_rk38, legend="rk38", color="Indigo")
    fig.x(norm_rk43_emb, fe_rk43_emb, legend="rk embedded", color="Green")
    fig.line(norm_rk43_emb, fe_rk43_emb, legend="rk embedded", color="Green")

    ##show(fig)
    
plot_embedded_cost()
