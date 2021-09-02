#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 02:00:59 2021

@author: mac
"""
import numpy as np

def grid_xt_to_omegak(grid_history):
    # FFT grid_history(x,t) to obtain grid_omegak(k,omega)
    grid_omegak=np.fft.rfft2(grid_history)
    grid_omegak=grid_omegak[:grid_omegak.shape[0]//2,:]
    grid_omegak=grid_omegak[::-1,:]
    return grid_omegak

def grid_xt_to_kt(grid_history):
    grid_kt=np.fft.rfft(grid_history,axis=0)
    return grid_kt

def distribution_function_grid(x,v_grid,dv,r,v,NG,Nv):
    distribution_function=np.histogram2d(r,v,bins=[NG,Nv])[0]
    distribution_function/=dv
    x_space_distribution=np.sum(distribution_function,axis=1)
    v_space_distribution=np.sum(distribution_function,axis=0)
    X,V=np.meshgrid(x,v_grid)
    
    return distribution_function, X, V, x_space_distribution, v_space_distribution

def calculate_field_grid_quantities(v_grid,distribution_function):
    
    return n_grid, u_grid, P_grid, T_grid