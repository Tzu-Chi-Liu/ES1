#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 02:00:59 2021

@author: mac
"""
import numpy as np

def grid_xt_to_omegak(k,grid_history):
    '''
    Use 2DFFT with forward normalization to obtain grid_omegak(k,omega) from grid_history(x,t).
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
        Normalization is chosen such that 
            the transformed grid_omegak has been scaled by 1./(NG*NT).
        grid_omegak[i,j]=FFTfield(k=i*dk=i*2.*np.pi/L,omega=j*domega=j*2.*np.pi/(DT*NT)).

    '''
    grid_omegak=np.fft.rfft2(grid_history,norm='forward') # rfft2 only keeps the left part 
    grid_omegak=grid_omegak[:len(k),:] # Keep only the upper left corner
    return grid_omegak

def grid_xt_to_kt(grid_history):
    grid_kt=np.fft.rfft(grid_history,axis=0,norm='forward')
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
        units['m']='(m_e)'
        units['q']='(e)'
        units['r']='(c/ω_p)'
        units['v']='(c)'
        units['rho_grid']='(C*ω_p/c)'
        units['phi_grid']='(m_e*c^2/e)'
        units['E_grid']='(m_e*c*ω_p/e)'
        units['k']='(ω_p/c)'
        units['omega']='(ω_p)'
        units['Momentum']='(m_e*c)'
        units['Energy']='(m_e*c^2)'
        
    return units

# def calculate_field_grid_quantities(r,v_grid,distribution_function,IW):
#     n_grid=np.zeros(NG)
#     u_grid=np.zeros(NG)
#     P_grid=np.zeros(NG)
#     T_grid=np.zeros(NG)
#     for i in range(N):
#         x=
#     return n_grid, u_grid, P_grid, T_grid

def energy(dx, NG, phi_grid, rho_grid, species_parameters, m, v, v_old,):
    E_F, E_D, E_T = 0.0, 0.0, 0.0
    # Field energy history
    # E_F=0.5*np.dot(phi_grid,rho_grid)*dx
    E_F=0.5*np.real(np.dot(np.fft.fft(phi_grid),
                                 np.conjugate(np.fft.fft(rho_grid))))*dx/NG
    # E_F[step]=0.5*epsilon0*np.dot(E_grid,E_grid)*dx
    
    # Particle drift and thermal energy histories
    particle_counter=0
    for species in range(len(species_parameters)):
        N_species=species_parameters[species]['N']
        m_species=species_parameters[species]['m']
        avg_v_drift=np.mean(0.5*(v[particle_counter:particle_counter+N_species]
                                  +v_old[particle_counter:particle_counter+N_species]))
        E_D_species=N_species*0.5*m_species*avg_v_drift**2.
        E_T_species=np.sum(0.5*m[particle_counter:particle_counter+N_species]\
                          *v[particle_counter:particle_counter+N_species]
                          *v_old[particle_counter:particle_counter+N_species])-E_D_species
            
        E_D+=E_D_species
        E_T+=E_T_species
        particle_counter+=N_species
    
    E_total = E_F + E_D + E_T
    
    return E_F, E_D, E_T, E_total

# def mode_energy():
#     return E

def plasma_parameters(input_txt_parameters,species_parameters, v, v_old):
    '''
    Calculate plasma frequency and Debye length for each species.

    Parameters
    ----------
    input_txt_parameters : TYPE
        DESCRIPTION.
    species_parameters : TYPE
        DESCRIPTION.

    Returns
    -------
    omega_plasma_species : numpy.ndarray
        np array of shape(NSP,), where the ith entry is the plasma frequency for the ith species.
    lambda_Debye_species : numpy.ndarray
        np array of shape(NSP,), where the ith entry is the Debye length for the ith species.

    '''
    epsilon0            = input_txt_parameters['epsilon0']        # vacuum permittivity
    k_B                 = input_txt_parameters['k_B']
    L                   = input_txt_parameters['L']        # length of box
    
    omega_plasma_species = np.zeros(len(species_parameters))
    lambda_Debye_species = np.zeros(len(species_parameters))
    
    particle_counter = 0
    for species in range(len(species_parameters)):
        N_species = species_parameters[species]['N']
        m_species = species_parameters[species]['m']
        q_species = species_parameters[species]['q']
        
        n_species = N_species/L
        T_species = m_species\
                    *np.std(0.5*(v[particle_counter:particle_counter+N_species]\
                                 +v_old[particle_counter:particle_counter+N_species]))\
                    /k_B
        
        omega_plasma_species[species] = np.sqrt((n_species*q_species**2.)/(epsilon0*m_species))
        lambda_Debye_species[species] = np.sqrt((epsilon0*k_B*T_species)/(n_species*q_species**2.))
        
        particle_counter += N_species
        
    return omega_plasma_species, lambda_Debye_species

# def calculate_omega_plasma_grid(n_grid,species_parameters,input_txt_parameters):
#     omega_plasma_squared_sum=np.zeros(len(n_grid))
    
#     for species in range(len(species_parameters)):
#         epsilon0=
#         m=
#         q=
#         for i in range(NG):
            
#             omega_plasma_squared_sum[i]+=n_grid[i]*q**2./(epsilon0*m)
    
#     omega_plasma_grid=np.sqrt(omega_plasma_squared_sum)
#     return omega_plasma_grid

# def calculate_lambda_Debye_grid():
#     return lambda_Debye_grid