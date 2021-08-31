#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 11:22:18 2021

@author: mac
"""
import numpy as np

def build(parameters):
    if InitialCondition=='Particle Pair Oscillation':
        N=2
        m=np.array([m_e,m_e])
        q=np.array([e_e,-e_e])
        r=np.array([L/4.,2*L/4.])
        v=np.array([0.,0.])
        
    if InitialCondition=='Plasma Oscillation':
        A          = 1e-5
        v0         = 1.         # velocity
        k          = 2.*np.pi/L # velocity wavenumber 
        n0         = 0.5*N/L
        v0_sigma   = 0.0
        
        # first half electrons, second half protons
        m=np.ones(N)*m_e
        m[N//2:]=m_p
        q=-1.*np.ones(N)*e_e
        q[N//2:]*=-1.
        r=np.append(np.linspace(0,L,N//2),np.linspace(0,L,N//2))
        v=np.append(np.random.normal(0,v0_sigma,size=(1,N//2)),np.zeros(N//2))
        
        # sinusoidal perturbation
        r[:N//2]+=A/(n0*k)*np.sin(k*r[:N//2])
        r=r%L
        
        # excite all modes
        # for i in range(1,NG//2+1):
        #     r[:N//2]+=A*i/(n0)*np.sin(i*k*r[:N//2])
        # r=r%L
        
        # excite first few modes
        # for i in range(1,51):
        #     r[:N//2]+=A/(n0)*np.sin(i*k*r[:N//2])
        # r=r%L
        
        # r[15*N//64:17*N//64]+=A/(n0)
    
    if InitialCondition=='Two Stream Instability':
        v0         = 1.0        # velocity
        v0_sigma   = 0.1        # velocity width
        A          = 0.01       # sinusoidal perturbation amplitude
        k          = 2.*np.pi/L # sinusoidal perturbation wavevector 
        
        m=np.ones(N)*m_e
        q=-1*np.ones(N)*e_e
        # q[N//2:]*=-1.
        r=np.random.random(N)*L
        # v=np.random.normal(v0,v0_sigma,size=(1,N))[0]
        # v[N//2:]*=-1.
        v=np.append(np.random.normal(v0,v0_sigma,size=(1,N//2)),
                    np.random.normal(-v0,v0_sigma,size=(1,N//2)))
        # v+=A*np.sin(k*r) # add perturbation
    
    if InitialCondition=='Single Beam':
        v0         = 1.0        # velocity
        v0_sigma   = 0.0        # velocity width
        A          = 0.01       # sinusoidal perturbation amplitude
        k          = 2.*np.pi/L # sinusoidal perturbation wavevector 
        
        m=np.ones(N)*m_e
        q=-1*np.ones(N)*e_e
        # q[N//2:]*=-1.
        r=np.random.random(N)*L
        # v=np.random.normal(v0,v0_sigma,size=(1,N))[0]
        # v[N//2:]*=-1.
        v=np.random.normal(v0,v0_sigma,size=(1,N))[0]
        
    return m,q,r,v