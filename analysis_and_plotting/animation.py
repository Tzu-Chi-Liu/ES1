#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 00:51:02 2021

@author: mac
"""
import numpy as np
import matplotlib.pyplot as plt

def diagnostics_animation(t,x,r,v,phi_grid,E_grid,rho_grid,step,N,NG,units,save_dir):
    global animation_fig
    if step==0:
        animation_fig=plt.figure(figsize=(16,8))
    
    # clear figure
    plt.clf()
    
    # phase space
    ax1=animation_fig.add_subplot(221)
    ax1.scatter(r[:N//2],v[:N//2],fc=(1,0,0,0.3),s=1)
    ax1.scatter(r[N//2:],v[N//2:],fc=(0,0,1,0.3),s=1)
    # ax1.hexbin(r,v,gridsize=200)
    # tracker particles 
    ax1.scatter(r[N//4],v[N//4],c='black',s=30,lw=1,edgecolor='white') # tracker 1
    ax1.scatter(r[-1],v[-1],c='white',s=30,lw=1,edgecolor='black') # tracker 2
    ax1.set_xlabel('position x '+units['r'])
    ax1.set_ylabel('velocity v '+units['v'])
    ax1.set_title('Phase Space')
    # ax1.set_xlim(0.,L)
    ax1.set_ylim(-1.1*np.max(np.abs(v)),1.1*np.max(np.abs(v)))
    # ax1.set_ylim(-7e-10,7e-10)
    
    # velocity distribution
    ax2=animation_fig.add_subplot(222)
    ax2.hist(v[:N//2],bins=100,range=((np.min(v),np.max(v))),
              color='r',orientation='horizontal')
    ax2.hist(v[N//2:],bins=100,range=((np.min(v),np.max(v))),
              color='b',alpha=0.5,orientation='horizontal')
    ax2.hist(v,bins=100,range=((np.min(v),np.max(v))),color='white',
              edgecolor='black',lw=1.,orientation='horizontal',histtype='step')
    ax2.set_xlabel('number of particles f(v)dv')
    ax2.set_ylabel('velocity v '+units['v'])
    ax2.set_title('Velocity Distribution')
    ax2.set_ylim(-1.1*np.max(np.abs(v)),1.1*np.max(np.abs(v)))
    # ax2.set_ylim(-7e-10,7e-10)
    
    # density at grid positions
    ax3=animation_fig.add_subplot(223)
    ax3.plot(x,rho_grid-np.mean(rho_grid))
    # ax3.hist(r,bins=NG)
    ax3.set_xlabel('position x '+units['r'])
    ax3.set_ylabel(r'charge density $\rho$ '+units['rho_grid'])
    ax3.set_title('Charge Density')
    ax3.set_ylim(-1.1*np.max(np.abs(rho_grid-np.mean(rho_grid))),
                  1.1*np.max(np.abs(rho_grid-np.mean(rho_grid))))
    # ax3.set_ylim(-6e-7,6e-7)
    # ax3.set_title('Density')
    
    # potential and field at grid positions
    ax4=animation_fig.add_subplot(224)
    ax8=ax4.twinx()
    ax4.plot(x,phi_grid,'blue')
    ax8.plot(x,E_grid,'orange')
    ax4.plot(x,np.zeros(NG),'black',lw=1)
    ax4.set_xlabel('position x '+units['r'])
    ax4.set_ylabel(r'potential $\phi$ '+units['phi_grid'],color='blue')
    ax8.set_ylabel('electric field E '+units['E_grid'],color='orange')
    ax4.set_title('Potential and Electric Fields')
    ax4.set_ylim(-1.1*np.max(np.abs(phi_grid)),1.1*np.max(np.abs(phi_grid)))
    ax8.set_ylim(-1.1*np.max(np.abs(E_grid)),1.1*np.max(np.abs(E_grid)))
    # ax4.set_ylim(-1.5e-8,1.5e-8)
    # ax8.set_ylim(-1e-7,1e-7)
    # ax4.set_title('Potential')
    
    animation_fig.suptitle('Snapshot at t = %.3f (s)'%t[step])
    plt.tight_layout()
    plt.pause(0.1)
        
    return animation_fig