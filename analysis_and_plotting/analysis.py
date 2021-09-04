#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 02:00:59 2021

@author: mac
"""
import numpy as np

def grid_xt_to_omegak(k,grid_history):
    '''
    Use 2DFFT to obtain grid_omegak(k,omega) from grid_history(x,t).
    Because The 2DFFT of a real matrix has symmetry, 
    only the upper left part of the 2DFFT result 
    corresponding to (0<=k<=np.pi/dx,0<=omega<=np.pi/dt) is returned.
    
    Parameters
    ----------
    grid_history : numpy.ndarray
        The field quantities in real space at grid points with shape=(NG,NT). 
        grid_history[i,j] = field(x=i*dx,j*dt).

    Returns
    -------
    grid_omegak : numpy.ndarray
        The transformed complex field(k,omega) of shape (len(k),len(omega)).
        grid_omegak[i,j]=FFTfield(k=i*dk=i*2.*np.pi/L,omega=j*domega=j*2.*np.pi/(DT*NT))

    '''
    grid_omegak=np.fft.rfft2(grid_history)
    grid_omegak=grid_omegak[:len(k),:]
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

# def calculate_field_grid_quantities(v_grid,distribution_function):
    
#     return n_grid, u_grid, P_grid, T_grid

# def mode_energy():
#     return E