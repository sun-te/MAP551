#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  9 17:07:05 2018

@author: te.sun
"""

import numpy as np
import pandas as pd
import os
os.chdir('/users/eleves-a/2017/te.sun/MAP551/bz/bz_1d/run')
from outil import write_input, read_data
from bokeh.io import  output_notebook, push_notebook, show
from bokeh.plotting import figure
from bokeh.layouts import column, row
from bokeh.palettes import Category10,Category20
from bokeh.models import PrintfTickFormatter

from ipywidgets import FloatProgress, IntProgress
from IPython.display import display
#%%
'''Question2.1'''
#type of model (nb of equations) : 2 or 3 
neq=2
#tini: intial time
tini=0.0
#xmin: left limit of domain
xmin=0.0
#xmax: right limit of domain
xmax=4.0
#nx: nb of spatial discretization points
nx=4001
#method: integration method (radau5, strang)
method='radau5'
#tol: tolerance for radau5 and rock4
tol='1.d-12'
#%%
sols=[]
t_plot=[0]

for tend in np.linspace(0.5,5.5,6):
    nt=int(tend*100)+1
    list_input=[neq,tini,tend,nt,xmin,xmax,nx,method,tol]
    write_input(list_input)
    os.system('./bz_1d')
    sol_final,num_final=read_data("sol_num.dat")
    sols.append(sol_final)
    t_plot.append(tend)
np.save('traveling_waves_neq'+str(neq)+'.npy',np.array(sols))
#%%
sols=np.load('traveling_waves_neq'+str(neq)+'.npy')
sol_init,num_init=read_data("sol_ini.dat")
x=sol_init[0]
t_plot=[0]
for tend in np.linspace(0.5,5.5,11):t_plot.append(tend)
if(neq==2):
    fig_solb = figure(x_range=(xmin,xmax), plot_height=400, plot_width=950, title="Solution")
    fig_solc = figure(x_range=(xmin,xmax), plot_height=400, plot_width=950, title="Solution")
    it=0
    fig_solc.line(x, sol_init[2],line_width=2,color="Green",legend=f'sol c(x) at t={t_plot[it]:.2f}')
    fig_solb.line(x, sol_init[1],line_width=2,color="Green",legend=f'sol b(x) at t={t_plot[it]:.2f}')
    for sol in sols:
        it+=1
        fig_solc.line(x, sol[2],line_width=2,color=Category20[20][it],legend=f'sol c(x) at t={t_plot[it]:.2f}')
        fig_solb.line(x, sol[1],line_width=2,color=Category20[20][it],legend=f'sol b(x) at t={t_plot[it]:.2f}')
    show(column(fig_solb,fig_solc))
if(neq==3):
    fig_sola = figure(x_range=(xmin,xmax+.5), plot_height=400, plot_width=950, title="Solution")
    fig_solb = figure(x_range=(xmin,xmax+.5), plot_height=400, plot_width=950, title="Solution")
    fig_solc = figure(x_range=(xmin,xmax+.5), plot_height=400, plot_width=950, title="Solution")
    it=0
    fig_solc.line(x, sol_init[3],line_width=2,color="Green",legend=f'sol c(x) at t={t_plot[it]}')
    fig_solb.line(x, sol_init[2],line_width=2,color="Green",legend=f'sol b(x) at t={t_plot[it]}')
    fig_sola.line(x, sol_init[1],line_width=2,color="Green",legend=f'sol a(x) at t={t_plot[it]}')
    for sol in sols:
        it+=1
        fig_solc.line(x, sol[3],line_width=2,color=Category20[20][it],legend=f'sol c(x) at t={t_plot[it]}')
        fig_solb.line(x, sol[2],line_width=2,color=Category20[20][it],legend=f'sol b(x) at t={t_plot[it]}')
        fig_sola.line(x, sol[1],line_width=2,color=Category20[20][it],legend=f'sol a(x) at t={t_plot[it]}')
    show(column(fig_sola,fig_solb,fig_solc))