#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 11:22:18 2021

@author: mac
"""
import numpy as np

# =============================================================================
# Build m, q, r, v from parameters
# =============================================================================

def load_particles(parameters):
    '''
    For generating the particle m, q, r, v np arrays from the parameters dictionary.

    Parameters
    ----------
    parameters : dict
        Parameters read from `input.txt` using input_output.load_input(input_file_loc).

    Returns
    -------
    m : numpy.ndarray
        Mass of each particle.
    q : numpy.ndarray
        Charge of each particle.
    r : numpy.ndarray
        Initial position of each particle.
    v : numpy.ndarray
        Initial velocity of each particle.

    '''
    
    # Species 1
    # N                   = 2000       # Number of particles
    # density             = 10.0       # not sure to input or calculate
    # WP
    # WC
    # QM
    
    # # distribution      = Maxwellian # Used for choosing different initial velocity distribution function
    # v0                  = 
    # v_sigma             = 0.0        # Velocity width (v0_sigma**2 ~ Temperature)
    # T                   = 0.1        # Temperature for a Maxwellian velocity distribution
    
    # Modes               = 1:50:1     # Excited modes (seperate by comma or i:j:k for modes = list(range(i,j,k)))
    # X1                  = 
    # V1                  = 
    # THETAX              = 
    # THETAV              = 
    
    
    
    k=2.*np.pi/parameters['L']
    density=parameters['N']/parameters['L']
    
    # if parameters['InitialCondition']=='Particle Pair Oscillation':
    #     N=2
    #     m=np.array([parameters['m_e'],parameters['m_e']])
    #     q=np.array([parameters['e'],-parameters['e']])
    #     r=np.array([parameters['L']/4.,2*parameters['L']/4.])
    #     v=np.array([0.,0.])
        
    if parameters['InitialCondition']=='Plasma_Oscillation':
        A          = parameters['A']
        modes      = parameters['Modes']
        v_sigma   = parameters['v_sigma']
        n0         = 0.5*parameters['N']/parameters['L']
        
        # first half electrons, second half protons
        m=np.ones(parameters['N'])*parameters['m_e']
        m[parameters['N']//2:]=parameters['m_p']
        q=-1.*np.ones(parameters['N'])*parameters['e']
        q[parameters['N']//2:]*=-1.
        r=np.append(np.linspace(0,parameters['L'],parameters['N'])//2,np.linspace(0,parameters['L'],parameters['N']//2))
        v=np.append(np.random.normal(0,v_sigma,size=(1,parameters['N']//2)),np.zeros(parameters['N']//2))
        
        
        # excite mode in modes list
        for mode in modes:
            r[:parameters['N']//2]+=A*mode/n0*np.sin(mode*k*r[:parameters['N']//2])
        r=r%parameters['L']
    
    if parameters['InitialCondition']=='Two_Stream_Instability':
        v0         = parameters['v0']     # velocity
        v0_sigma   = parameters['v0_sigma']       # velocity width
        charge     = parameters['charge']
        A          = parameters['A']
        modes      = parameters['Modes']
        
        m=np.ones(parameters['N'])*parameters['m_e'] 
        q=-1*np.ones(parameters['N'])*parameters['e']
        if charge == 'opposite':
            q[parameters['N']//2:]*=-1.
        r=np.random.random(parameters['N'])*parameters['L']
        v=np.append(np.random.normal(v0,v0_sigma,size=(1,parameters['N']//2)),
                    np.random.normal(-v0,v0_sigma,size=(1,parameters['N']//2)))
        
        for mode in modes:
            r[:parameters['N']//2]+=A*mode/density*np.sin(mode*k*r[:parameters['N']//2])
        r=r%parameters['L']
        
    # if parameters['InitialCondition']=='Single Beam':
    #     v0         = 1.0        # velocity
    #     v0_sigma   = 0.0        # velocity width
    #     A          = 0.01       # sinusoidal perturbation amplitude
        
    #     m=np.ones(parameters['N']))*parameters['m_e']
    #     q=-1*np.ones(parameters['N']))*parameters['e']
    #     # q[parameters['N'])//2:]*=-1.
    #     r=np.random.random(parameters['N']))*parameters['L']
    #     # v=np.random.normal(v0,v0_sigma,size=(1,parameters['N'])))[0]
    #     # v[parameters['N'])//2:]*=-1.
    #     v=np.random.normal(v0,v0_sigma,size=(1,parameters['N'])))[0]
        
    return m,q,r,v

