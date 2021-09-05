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

def generate_units(UnitSystem):
    units={'t':'(arb. unit)',
       'm':'(arb. unit)',
       'q':'(arb. unit)',
       'r':'(arb. unit)',
       'v':'(arb. unit)',
       'rho_grid':'(arb. unit)',
       'phi_grid':'(arb. unit)',
       'E_grid':'(arb. unit)',
       'k':'(arb. unit)',
       'omega':'(arb. unit)',
       'Momentum':'(arb. unit)',
       'Energy':'(arb. unit)',
       'arb. unit':'(arb. unit)'
       }

    if UnitSystem=='SI':
        units['t']='(s)'
        units['m']='(kg)'
        units['q']='(C)'
        units['r']='(m)'
        units['v']='(m/s)'
        units['rho_grid']='(C/m)'
        units['phi_grid']='(V)'
        units['E_grid']='(V/m)'
        units['k']='(1/m)'
        units['omega']='(1/s)'
        units['Momentum']='(kg*m/s)'
        units['Energy']='(J)'
        
    if UnitSystem=='Normalized':
        units['t']='(1/ω_p)'
        units['m']='(kg)'
        units['q']='(C)'
        units['r']='(c/ω_p)'
        units['v']='(c)'
        units['rho_grid']='(C*ω_p/c)'
        units['phi_grid']='(V)'
        units['E_grid']='(V*ω_p/c)'
        units['k']='(ω_p/c)'
        units['omega']='(ω_p)'
        units['Momentum']='(kg*c)'
        units['Energy']='(J)'
        
    return units
# def calculate_field_grid_quantities(v_grid,distribution_function):
    
#     return n_grid, u_grid, P_grid, T_grid

# def mode_energy():
#     return E

# def calculate_omega_plasma():
#     return omega_plasma

# def calculate_lambda_Debye():
#     return lambda_Debye