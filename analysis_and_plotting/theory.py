#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  6 01:47:10 2021

@author: mac
"""
import numpy as np

def dispersion_relation(parameters):
    e                   = parameters['e']     # elementary charge
    epsilon0            = parameters['epsilon0']        # vacuum permittivity
    m_e                 = parameters['m_e']        # electron mass
    m_p                 = parameters['m_p']     # proton mass
    k_B                 = parameters['k_B']
    
    InitialCondition    = parameters['InitialCondition']
    
    N                   = parameters['N']      # number of particles
    L                   = parameters['L']        # length of box
    NG                  = parameters['NG']         # number of grids
    dx                  = L/NG # grid spacing
    k                   = 2.*np.pi*np.fft.rfftfreq(NG,dx)
    
    omega_plasma        = np.sqrt(((N/(2*L))*e**2)/(epsilon0*m_e)) # plasma frequency
    
    if InitialCondition=='Plasma_Oscillation':
        v0         = parameters['v0']        # velocity
        v_sigma   = parameters['v_sigma']
        T         = parameters['T']
        n0         = 0.5*parameters['N']/parameters['L']
        
        # first half electrons, second half protons

        
        # theoretical dispersion relation
        if v_sigma==0:
            theoretical_omega_of_k=omega_plasma*np.cos(0.5*dx*k)
        else:
            theoretical_omega_of_k=np.sqrt(omega_plasma**2+3*v_sigma**2*k**2)
        
        # theoretical growth rate
        theoretical_growth_rate=np.zeros_like(k)
        
    if InitialCondition=='Two_Stream_Instability':
        v0         = parameters['v0']        # velocity
        v0_sigma   = parameters['v0_sigma']        # velocity width
        charge     = parameters['charge']
        
        # theoretical dispersion relation & growth rate
        # theoretical_omeag_of_k=[np.zeros_like(k),np.zeros_like(k)]
        # theoretical_growth_rate=np.zeros_like(k)
        # for mode in range(len(k)):
        #     if omega_plasma*np.sqrt(4.*k[mode]**2*v0**2+omega_plasma**2)>k[mode]**2*v0**2+omega_plasma**2:
        #         theoretical_growth_rate[mode]=np.sqrt(-k**2*v0**2-omega_plasma**2
        #                                               +omega_plasma*np.sqrt(4.*k**2*v0**2+omega_plasma**2))
        #     else:
        #         theoretical_omega_of_k[0][mode]=np.sqrt(k**2*v0**2+omega_plasma**2
        #                                                 +omega_plasma*np.sqrt(4.*k**2*v0**2+omega_plasma**2))
        #         theoretical_omega_of_k[1][mode]=np.sqrt(k**2*v0**2+omega_plasma**2
        #                                                 -omega_plasma*np.sqrt(4.*k**2*v0**2+omega_plasma**2))
        theoretical_omega_of_k=[np.sqrt(k**2*v0**2+omega_plasma**2
                                        +omega_plasma*np.sqrt(4.*k**2*v0**2+omega_plasma**2)),
                                -np.sqrt(k**2*v0**2+omega_plasma**2
                                        +omega_plasma*np.sqrt(4.*k**2*v0**2+omega_plasma**2)),
                                np.sqrt(k**2*v0**2+omega_plasma**2
                                        -omega_plasma*np.sqrt(4.*k**2*v0**2+omega_plasma**2)),
                                -np.sqrt(k**2*v0**2+omega_plasma**2
                                        -omega_plasma*np.sqrt(4.*k**2*v0**2+omega_plasma**2))]
        for omegak in theoretical_omega_of_k:
            omegak[np.isnan(omegak)] = 0.
        
        theoretical_growth_rate=np.sqrt(-k**2*v0**2-omega_plasma**2
                                        +omega_plasma*np.sqrt(4.*k**2*v0**2+omega_plasma**2))
        theoretical_growth_rate[np.isnan(theoretical_growth_rate)]=0.
        
    return theoretical_omega_of_k,theoretical_growth_rate