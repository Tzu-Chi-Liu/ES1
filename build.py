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

def build(parameters):
    '''
    For generating the particle m, q, r, v np arrays from the parameters dictionary.

    Parameters
    ----------
    parameters : dict
        Parameters read from `input.txt`.

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
    # if parameters['InitialCondition']=='Particle Pair Oscillation':
    #     N=2
    #     m=np.array([m_e,m_e])
    #     q=np.array([e_e,-e_e])
    #     r=np.array([L/4.,2*L/4.])
    #     v=np.array([0.,0.])
        
    if parameters['InitialCondition']=='Plasma_Oscillation':
        A          = float(parameters['A'])
        if ':' in parameters['Modes']:
            modes  = list(range(parameters['Modes'].split(':')[0],parameters['Modes'].split(':')[1]))
        else:
            modes  = [int(mode) for mode in parameters['Modes'].split(',')]
        v_sigma   = float(parameters['v_sigma'])
        v_gt          = float(parameters['v_gt'])         # velocity
        n0         = 0.5*float(parameters['N'])/float(parameters['L'])
        
        # first half electrons, second half protons
        m=np.ones(N)*m_e
        m[N//2:]=m_p
        q=-1.*np.ones(N)*e_e
        q[N//2:]*=-1.
        r=np.append(np.linspace(0,L,N//2),np.linspace(0,L,N//2))
        v=np.append(np.random.normal(0,v_sigma,size=(1,N//2)),np.zeros(N//2))
        
        # Galilean transformation
        v+=v_gt
        
        # excite mode in modes list
        for mode in modes:
            r[:N//2]+=A*mode/n0*np.sin(mode*k*r[:N//2])
        r=r%L
    
    if parameters['InitialCondition']=='Two_Stream_Instability':
        v0         = float(parameters['v0'])        # velocity
        v0_sigma   = float(parameters['v0_sigma'])        # velocity width
        charge     = parameters['charge']
        A          = float(parameters['A'])
        if ':' in parameters['Modes']:
            modes  = list(range(parameters['Modes'].split(':')[0],parameters['Modes'].split(':')[1]))
        else:
            modes  = [int(mode) for mode in parameters['Modes'].split(',')]
        
        m=np.ones(N)*m_e 
        q=-1*np.ones(N)*e_e
        if charge == 'opposite':
            q[N//2:]*=-1.
        r=np.random.random(N)*L
        v=np.append(np.random.normal(v0,v0_sigma,size=(1,N//2)),
                    np.random.normal(-v0,v0_sigma,size=(1,N//2)))
        
        for mode in modes:
            r[:N//2]+=A*mode/density*np.sin(mode*k*r[:N//2])
        r=r%L
        
    # if parameters['InitialCondition']=='Single Beam':
    #     v0         = 1.0        # velocity
    #     v0_sigma   = 0.0        # velocity width
    #     A          = 0.01       # sinusoidal perturbation amplitude
        
    #     m=np.ones(N)*m_e
    #     q=-1*np.ones(N)*e_e
    #     # q[N//2:]*=-1.
    #     r=np.random.random(N)*L
    #     # v=np.random.normal(v0,v0_sigma,size=(1,N))[0]
    #     # v[N//2:]*=-1.
    #     v=np.random.normal(v0,v0_sigma,size=(1,N))[0]
        
    return m,q,r,v