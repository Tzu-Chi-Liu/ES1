#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 03:01:12 2021

@author: mac
"""
import numpy as np

# =============================================================================
# source weighting
# =============================================================================
def weightrho(r,e,IW,N,NG,dx):
    # Weighting particle positions r to grid density rho_grid
    rho_grid=np.zeros(NG)
    for i in range(N):
        x=int(r[i]//dx)
        if IW=='NGP' or IW=='1':
            if np.abs(r[i]%dx)<0.5*dx:
                rho_grid[x%NG]+=e[i]/dx
            else:
                rho_grid[(x+1)%NG]+=e[i]/dx
        if IW=='CIC' or IW=='2' or IW=='energy conserving' or IW=='3':
            rho_grid[x%NG]+=e[i]/dx*(dx-r[i]%dx)/dx
            rho_grid[(x+1)%NG]+=e[i]/dx*(r[i]%dx)/dx
            
    return rho_grid

# =============================================================================
# Poisson solver
# =============================================================================
def solve_poisson(r,e,IW,NG,dx,epsilon0,solver='FFT',operator='local'):
    rho_grid=weightrho(r,e,IW)
    
    if solver=='FFT':
        # Use FFT to solve Poisson equation with periodic boundary condition
        rho_k=np.fft.rfft(rho_grid)
        k=2*np.pi*np.fft.rfftfreq(NG,dx)
        k[0]=1.     # to avoid division by 0
        rho_k[0]=0. # set DC component (average) of rho to 0 (average charge density = 0)
        
        # laplacian
        if operator=='nonlocal':
            phi_k=rho_k/(epsilon0*k**2)
        if operator=='local':
            K=k*(np.sin(0.5*k*dx)/(0.5*k*dx)) # Fourier transform of discrete laplacian
            phi_k=rho_k/(epsilon0*K**2)
        phi_grid=np.fft.irfft(phi_k)
    
    # Differentiate potential phi to obtain field E 
    gradient=(np.diagflat(np.ones(NG-1),1)-np.diagflat(np.ones(NG-1),-1))*1/(2.*dx)
    # gradient operator also has periodic BC
    gradient[0,NG-1]=-1./(2.*dx)
    gradient[NG-1,0]=1./(2.*dx)
    E_grid=-gradient@phi_grid
    
    return phi_grid,rho_grid,E_grid,k,rho_k,phi_k

# =============================================================================
# force weighting
# =============================================================================
def weightField(m,q,r,v,E_grid,IW,NG,dx):
    # Weighting grid field E to obtain field E_p at particle position r
    X=int(r//dx)
    if IW=='1' or IW=='NGP' or IW=='energy conserving' or IW=='3':
        if np.abs(r%dx)<(0.5)*dx:    
            E_par=E_grid[X%NG]
        else:
            E_par=E_grid[(X+1)%NG]
    if IW=='2' or IW=='CIC':
        E_par=E_grid[X%NG]*(dx-r%dx)/dx+E_grid[(X+1)%NG]*(r%dx)/dx
        
    return E_par

# =============================================================================
# integration of EOM
# =============================================================================
def update_motion(m,q,r,v,E_grid,scheme,IW,DT,L):
    global step 
    E_par=weightField(m,q,r,v,E_grid,IW)
    a=q*E_par/m
    
    if scheme=='leapfrog':
        if step==0:
            v-=a*0.5*DT
            
        v+=a*DT
        
        r+=v*DT
        r=r%L
        
    if scheme=='KDK':
        v+=a*0.5*DT
        
        r+=v*DT
        r=r%L
        
        E_par=weightField(m,q,r,v,E_grid,IW)
        a=q*E_par/m
        v+=a*0.5*DT
        
    if scheme=='DKD':
        r+=v*0.5*DT
        r=r%L
        
        v+=a*DT
        
        r+=v*0.5*DT
        r=r%L
        