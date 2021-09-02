#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 02:00:29 2021

@author: mac
"""
import numpy as np
import matplotlib.pyplot as plt

def energy_history_plot(t,E_D,E_T,E_F,units,tlim=None,Elim=None):
    fig,ax=plt.subplots()
    ax.plot(t,E_D,label='Drift')
    ax.plot(t,E_T,label='Thermal')
    ax.plot(t,E_F,label='Field')
    ax.plot(t,E_D+E_T+E_F,label='Total')
    
    if tlim is not None:
        ax.set_xlim(tlim)
    if Elim is not None:
        ax.set_ylim(Elim)
    ax.set_xlabel('Time t '+units['t'])
    ax.set_ylabel('Energy E '+units['Energy'])
    ax.legend(bbox_to_anchor=(0,1.01,1,0.2),loc='lower left'
               ,mode='expand',borderaxespad=0.,ncol=4)
    ax.set_title('Energy History',y=1.07)
    
    return fig

def momentum_change_history_plot(t,P,units,tlim=None,Plim=None):
    fig,ax=plt.subplots()
    ax.plot(t,P-P[0])
    
    if tlim is not None:
        ax.set_xlim(tlim)
    if Plim is not None:
        ax.set_ylim(Plim)
    ax.set_xlabel('Time t '+units['t'])
    ax.set_ylabel('Total momentum P '+units['Momentum'])
    ax.set_title('Total Momentum Change History')
    return fig

def grid_history_plot(x,t,dx,DT,grid_history,units,projection='2d',xlim=None,tlim=None):
    fig=plt.figure()
    if projection=='2d':
        ax=fig.add_subplots()
        pos=ax.imshow(grid_history.T[xlim[0]:xlim[1],tlim[0]:tlim[1]],
                      extent=(x[xlim[0]]-0.5*dx,x[xlim[1]]-0.5*dx,
                              t[tlim[0]]-0.5*DT,t[tlim[1]]-0.5*DT))
        
        cbar=fig.colorbar(pos,ax=ax)
        cbar.set_label('electric field E '+units['E_grid'],rotation=90)
        
    elif projection=='3d':
        ax=fig.add_subplot(projection='3d')
        x_mesh,T_mesh=np.meshgrid(x,t)
        ax.plot_surface(T_mesh,x_mesh,grid_history.T)
    
        if xlim is not None:
            ax.set_xlim(xlim)
        if tlim is not None:
            ax.set_ylim(tlim)
    ax.set_xlabel('Time t '+units['t'])
    ax.set_ylabel('Position x '+units['r'])
    ax.set_zlabel(r'Potential $\phi$ '+units['phi_grid'])
    ax.set_title('Field history')
    
    return fig

def selected_mode_history_plot(k,t,grid_kt,selected_modes,
                               plot_theoretical_growth_rate,theoretical_growth_rate,units,
                               scale='linear',klim=None,tlim=None):
    fig=plt.figure(figsize=(8,8))
    for index in range(len(selected_modes)):
        ax=fig.add_subplot(len(selected_modes),1,index+1)
        ax.plot(t[:],np.abs(grid_kt[selected_modes[index],:]),label='Simulation result')
        if plot_theoretical_growth_rate:
            growth_rate_label='Theoretical growth rate = %.3f'\
                               %theoretical_growth_rate[selected_modes[index]]\
                               +' '+units['omega']
            ax.plot(t[:],np.abs(grid_kt[selected_modes[index],0])
                    *np.exp(theoretical_growth_rate[selected_modes[index]]*t[:]),
                    label=growth_rate_label)
           
            
        title='Mode '+str(selected_modes[index])
        if klim is not None:
            ax.set_xlim(klim)
        if tlim is not None:
            ax.set_ylim(tlim)
        ax.set_xlabel('Time t '+units['t'])
        ax.set_ylabel('Amplitude '+units['arb. unit'])
        ax.set_yscale(scale)
        ax.legend(loc='upper left')
        ax.set_title(title)
        
    plt.tight_layout()
    return fig

def all_mode_history_plot(k,t,dk,DT,grid_kt,units,klim=None,tlim=None):
    fig,ax=plt.subplots()
    pos=ax.imshow(np.abs(grid_kt)[klim[0]:klim[1],tlim[0]:tlim[1]],
                  extent=(k[klim[0]]-0.5*dk,k[klim[1]]-0.5*dk,
                          t[tlim[0]]-0.5*DT,t[tlim[1]]-0.5*DT))
    
    if klim is not None:
        ax.set_xlim(klim)
    if tlim is not None:
        ax.set_ylim(tlim)
    # ax11.set_xlabel(r'Time $t$ '+units['t'])
    # ax11.set_ylabel(r'Wave number $k$ '+units['k'])
    cbar=fig.colorbar(pos,ax=ax)
    ax.set_xlabel('Step #')
    ax.set_ylabel('Mode number')
    ax.set_title('')
    cbar.set_label('Amplitude '+units['arb. unit'],rotation=90)
    return fig

def tracker_particle_trajectory_plot(t,R,V,Tracker,units,
                                     space=['R-t'],tlim=None,Rlim=None,Vlim=None):
    fig=plt.figure()
    if 'R-t' in space:
        ax=fig.add_subplot(len(space),1,space.index('R-t')+1)
        for particle in range(Tracker):
            ax.plot(t,R[particle,:],label='Particle '+str(particle+1))
        
        if tlim is not None:
            ax.set_xlim(tlim)
        if Rlim is not None:
            ax.set_ylim(Rlim)
        ax.set_xlabel('Time t '+units['t'])
        ax.set_ylabel('Position r '+units['r'])
        ax.legend(loc='upper left')
        ax.set_title('R-t Plot')
        
    if 'V-t' in space:
        ax=fig.add_subplot(len(space),1,space.index('V-t')+1)
        for particle in range(Tracker):
            ax.plot(t,V[particle,:],label='Particle '+str(particle+1))
        
        if tlim is not None:
            ax.set_xlim(tlim)
        if Vlim is not None:
            ax.set_ylim(Vlim)
        ax.set_xlabel('Time t '+units['t'])
        ax.set_ylabel('Velocity v '+units['v'])
        ax.legend(loc='upper left')
        ax.set_title('V-t Plot')
        
    if 'R-V' in space:
        ax=fig.add_subplot(len(space),1,space.index('R-V')+1)
        for particle in range(Tracker):
            ax.plot(R[particle,:],V[particle,:],label='Particle '+str(particle+1))
        
        if Rlim is not None:
            ax.set_xlim(Rlim)
        if Vlim is not None:
            ax.set_ylim(Vlim)
        ax.set_xlabel('Position r '+units['r'])
        ax.set_ylabel('Velocity v'+units['v'])
        ax.legend(loc='upper left')
        ax.set_title('R-V Plot')
        
    return fig

# =============================================================================
# dispersion relation
# =============================================================================
def dispersion_relation_plot(k,omega,dk,domega,grid_omegak,units,
                             plot_theoretical_dispersion_relation,theoretical_omega_of_k,
                             klim=None,omegalim=None):
    fig,ax=plt.subplots()
    pos=ax.imshow(np.abs(grid_omegak)[klim[0]:klim[1],omegalim[0]:omegalim[1]],
                  extent=(k[klim[0]]-0.5*dk,k[klim[1]]-0.5*dk,
                          omega[omegalim[0]]-0.5*domega,omega[omegalim[1]]-0.5*domega))
    cbar=fig.colorbar(pos,ax=ax)
    ax.set_xlabel(r'wave number $k$ '+units['k'])
    ax.set_ylabel(r'angular frequency $\omega$ '+units['omega'])
    ax.set_title('Dispersion Relation')
    cbar.set_label('Amplitude '+units['arb. unit'],rotation=90)
    
    # plot theoretical dispersion relation    
    if plot_theoretical_dispersion_relation:
        ax.plot(k,theoretical_omega_of_k,'red',lw=0.5,label='theoretical')
        ax.legend(loc='upper right')

    return fig 

def phase_space_plot(r,v,N,units):
    fig,ax=plt.subplots()
    ax.scatter(r[:N//2],v[:N//2],fc=(1,0,0,0.3),s=1)
    ax.scatter(r[N//2:],v[N//2:],fc=(0,0,1,0.3),s=1)
    # ax1.hexbin(r,v,gridsize=200)
    ax.scatter(r[N//4],v[N//4],c='black',s=30,lw=1,edgecolor='white') # tracker 1
    ax.scatter(r[-1],v[-1],c='white',s=30,lw=1,edgecolor='black') # tracker 2
    
    ax.set_xlabel('position x '+units['r'])
    ax.set_ylabel('velocity v '+units['v'])
    ax.set_title('Phase Space')
    ax.set_xlim()
    ax.set_ylim()
    
    return fig

def distribution_function_grid_plot(X,V,distribution_function,projection,units):
    fig=plt.figure()
    if projection=='2d':
        ax=fig.add_subplots()
        ax.imshow(distribution_function)
    if projection=='3d':
        ax=fig.add_subplots(projection=='3d')
        ax.plot_surface(X,V,distribution_function.T)
    ax.set_xlabel('position x '+units['r'])
    ax.set_ylabel('velocity v '+units['v'])
    ax.set_title('Phase Space')
    
    return fig

def grid_plot(x,grid,units):
    fig,ax=plt.subplots()
    ax.plot(x,grid)
    ax.set_xlabel('position x '+units['r'])
    ax.set_ylabel(r'charge density $\rho$ '+units['rho_grid'])
    ax.set_title('Charge Density')
    return fig

