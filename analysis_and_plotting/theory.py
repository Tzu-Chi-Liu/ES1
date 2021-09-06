#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  6 01:47:10 2021

@author: mac
"""
import numpy as np

def dispersion_relation(input_txt_parameters,species_parameters):
    e                   = input_txt_parameters['e']     # elementary charge
    epsilon0            = input_txt_parameters['epsilon0']        # vacuum permittivity
    m_e                 = input_txt_parameters['m_e']        # electron mass
    m_p                 = input_txt_parameters['m_p']     # proton mass
    k_B                 = input_txt_parameters['k_B']
    
    L                   = input_txt_parameters['L']        # length of box
    NG                  = input_txt_parameters['NG']         # number of grids
    dx                  = L/NG # grid spacing
    k                   = 2.*np.pi*np.fft.rfftfreq(NG,dx)
    
    InitialCondition    = input_txt_parameters['InitialCondition']
    
    if InitialCondition=='Plasma_Oscillation':
        omega_plasma_sum_squared=0
        for species in range(len(species_parameters)):
            N=species_parameters[species]['N']
            
            if species_parameters[species]['m']=='m_e':
                m_species=m_e
            elif species_parameters[species]['m']=='m_p':
                m_species=m_p
            else:
                m_species=species_parameters[species]['m']
                
            if species_parameters[species]['q']=='e':
                q_species=e
            elif species_parameters[species]['q']=='-e':
                q_species=-e
            else:
                q_species=species_parameters[species]['q']
            
            # Assign mass and charge to each particle
            v0         = species_parameters[species]['v0']
            v_sigma    = species_parameters[species]['v_sigma']        # Velocity width (v0_sigma**2 ~ Temperature)
            T          = species_parameters[species]['T']        # Temperature for a Maxwellian velocity distribution
            
            # n0         = 0.5*parameters['N']/parameters['L']
        
            omega_plasma_species    = np.sqrt(((N/L)*q_species**2.)/(epsilon0*m_species)) # plasma frequency
            omega_plasma_sum_squared+=omega_plasma_species**2.
        
        omega_plasma=np.sqrt(omega_plasma_sum_squared)
        # theoretical dispersion relation
        if v_sigma==0:
            theoretical_omega_of_k=omega_plasma*np.cos(0.5*dx*k)
        else:
            theoretical_omega_of_k=np.sqrt(omega_plasma**2+3*v_sigma**2*k**2)
        
        # theoretical growth rate
        theoretical_growth_rate=np.zeros_like(k)
        
    # if InitialCondition=='Two_Stream_Instability':
    #     v0         = parameters['v0']        # velocity
    #     v0_sigma   = parameters['v0_sigma']        # velocity width
    #     charge     = parameters['charge']
        
    #     # theoretical dispersion relation & growth rate
    #     # theoretical_omeag_of_k=[np.zeros_like(k),np.zeros_like(k)]
    #     # theoretical_growth_rate=np.zeros_like(k)
    #     # for mode in range(len(k)):
    #     #     if omega_plasma*np.sqrt(4.*k[mode]**2*v0**2+omega_plasma**2)>k[mode]**2*v0**2+omega_plasma**2:
    #     #         theoretical_growth_rate[mode]=np.sqrt(-k**2*v0**2-omega_plasma**2
    #     #                                               +omega_plasma*np.sqrt(4.*k**2*v0**2+omega_plasma**2))
    #     #     else:
    #     #         theoretical_omega_of_k[0][mode]=np.sqrt(k**2*v0**2+omega_plasma**2
    #     #                                                 +omega_plasma*np.sqrt(4.*k**2*v0**2+omega_plasma**2))
    #     #         theoretical_omega_of_k[1][mode]=np.sqrt(k**2*v0**2+omega_plasma**2
    #     #                                                 -omega_plasma*np.sqrt(4.*k**2*v0**2+omega_plasma**2))
    #     theoretical_omega_of_k=[np.sqrt(k**2*v0**2+omega_plasma**2
    #                                     +omega_plasma*np.sqrt(4.*k**2*v0**2+omega_plasma**2)),
    #                             -np.sqrt(k**2*v0**2+omega_plasma**2
    #                                     +omega_plasma*np.sqrt(4.*k**2*v0**2+omega_plasma**2)),
    #                             np.sqrt(k**2*v0**2+omega_plasma**2
    #                                     -omega_plasma*np.sqrt(4.*k**2*v0**2+omega_plasma**2)),
    #                             -np.sqrt(k**2*v0**2+omega_plasma**2
    #                                     -omega_plasma*np.sqrt(4.*k**2*v0**2+omega_plasma**2))]
    #     for omegak in theoretical_omega_of_k:
    #         omegak[np.isnan(omegak)] = 0.
        
    #     theoretical_growth_rate=np.sqrt(-k**2*v0**2-omega_plasma**2
    #                                     +omega_plasma*np.sqrt(4.*k**2*v0**2+omega_plasma**2))
    #     theoretical_growth_rate[np.isnan(theoretical_growth_rate)]=0.
        
    else:
        theoretical_omega_of_k=np.zeros_like(k)
        theoretical_growth_rate=np.zeros_like(k)
    return theoretical_omega_of_k,theoretical_growth_rate