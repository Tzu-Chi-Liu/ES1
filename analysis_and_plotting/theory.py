#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  6 01:47:10 2021

@author: mac
"""
import numpy as np

def dispersion_relation(input_txt_parameters,species_parameters):
    epsilon0            = input_txt_parameters['epsilon0']        # vacuum permittivity
    k_B                 = input_txt_parameters['k_B']
    
    L                   = input_txt_parameters['L']        # length of box
    NG                  = input_txt_parameters['NG']         # number of grids
    dx                  = L/NG # grid spacing
    k                   = 2.*np.pi*np.fft.rfftfreq(NG,dx)
    
    InitialCondition    = input_txt_parameters['InitialCondition']
    
    # Calculate plasma frequency for each species
    omega_plasma_species=np.zeros(len(species_parameters))
    for species in range(len(species_parameters)):
        N=species_parameters[species]['N']
        
        m_species=species_parameters[species]['m']
        q_species=species_parameters[species]['q']
        
        n          = N/L
        omega_plasma_species[species]    = np.sqrt((n*q_species**2.)/(epsilon0*m_species))
    
    if InitialCondition=='Plasma_Oscillation':
        v0         = species_parameters[0]['v0']
        v_sigma    = species_parameters[0]['v_sigma']        
        T          = species_parameters[0]['T']   
        
        omega_plasma=np.sqrt(np.sum(omega_plasma_species**2.))
        
        # theoretical dispersion relation
        theoretical_omega_of_k=[np.zeros_like(k),np.zeros_like(k)]
        if v_sigma==0:
            theoretical_omega_of_k[0]=k*v0+omega_plasma*np.cos(0.5*dx*k)
            theoretical_omega_of_k[1]=k*v0-omega_plasma*np.cos(0.5*dx*k)
        else:
            theoretical_omega_of_k[0]=k*v0+np.sqrt(omega_plasma**2+3*v_sigma**2*k**2)
            theoretical_omega_of_k[1]=k*v0-np.sqrt(omega_plasma**2+3*v_sigma**2*k**2)
            
        # theoretical growth rate
        theoretical_growth_rate=np.zeros_like(k)
        
    elif InitialCondition=='Two_Stream_Instability':
        v0         = species_parameters[0]['v0']        # velocity
        # v0_sigma   = species_parameters[0]['v0_sigma']        # velocity width
        # charge     = parameters['charge']
        omega_plasma=omega_plasma_species[0]
        
        # theoretical dispersion relation & growth rate
        theoretical_omega_of_k=[np.zeros_like(k),np.zeros_like(k)]
        theoretical_growth_rate=np.zeros_like(k)
        for mode in range(len(k)):
            if omega_plasma*np.sqrt(4.*k[mode]**2*v0**2+omega_plasma**2)\
                >k[mode]**2*v0**2+omega_plasma**2:
                theoretical_growth_rate[mode]=np.sqrt(-k[mode]**2*v0**2-omega_plasma**2
                                                      +omega_plasma*np.sqrt(4.*k[mode]**2*v0**2
                                                                            +omega_plasma**2))
            else:
                theoretical_omega_of_k[1][mode]=np.sqrt(k[mode]**2*v0**2+omega_plasma**2
                                                        -omega_plasma*np.sqrt(4.*k[mode]**2*v0**2
                                                                              +omega_plasma**2))
            theoretical_omega_of_k[0][mode]=np.sqrt(k[mode]**2*v0**2+omega_plasma**2
                                                        +omega_plasma*np.sqrt(4.*k[mode]**2*v0**2
                                                                              +omega_plasma**2))
        
    else:
        theoretical_omega_of_k=np.zeros_like(k)
        theoretical_growth_rate=np.zeros_like(k)
        
    return theoretical_omega_of_k,theoretical_growth_rate