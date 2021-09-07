#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 02:00:29 2021

@author: mac
"""
import numpy as np
import matplotlib.pyplot as plt

def energy_history_plot(t,E_D,E_T,E_F,units,tlim='default',Elim='default'):
    fig,ax=plt.subplots()
    ax.plot(t,E_D,label='Drift')
    ax.plot(t,E_T,label='Thermal')
    ax.plot(t,E_F,label='Field')
    ax.plot(t,E_D+E_T+E_F,label='Total')
    
    if tlim != 'default':
        ax.set_xlim(tlim)
    if Elim != 'default':
        ax.set_ylim(Elim)
    ax.set_xlabel('Time t '+units['t'])
    ax.set_ylabel('Energy E '+units['Energy'])
    ax.legend(loc='best')
    # ax.legend(bbox_to_anchor=(0,1.01,1,0.2),loc='best'
    #            ,mode='expand',borderaxespad=0.,ncol=4)
    ax.set_title('Energy History',y=1.07)
    
    return fig, ax

def momentum_change_history_plot(t,P,units,tlim='default',Plim='default'):
    fig,ax=plt.subplots()
    ax.plot(t,P)
    
    if tlim != 'default':
        ax.set_xlim(tlim)
    if Plim != 'default':
        ax.set_ylim(Plim)
    ax.set_xlabel('Time t '+units['t'])
    ax.set_ylabel(r'Total momentum change fraction $\frac{\Delta P}{\sum \left|P\right|}$ '
                  +units['Momentum'])
    ax.set_title('Total Momentum Change History')
    return fig, ax

def grid_history_plot(x_grid,t,dx,DT,grid_history,units,
                      projection='2d',tlim='default',xlim='default'):
    fig=plt.figure()
    if projection=='2d':
        ax=fig.add_subplot()
        grid_history=grid_history[::-1,:] # So that (x = 0, t = 0) is at bottom left corner
        
        pos=ax.imshow(grid_history,
                      extent=(t[0]-0.5*DT,t[-1]+0.5*DT,
                              x_grid[0]-0.5*dx,x_grid[-1]+0.5*dx))
        
        if tlim != 'default':
            ax.set_xlim(tlim)
        if xlim != 'default':
            ax.set_ylim(xlim)
        cbar=fig.colorbar(pos,ax=ax)
        cbar.set_label('Electric field E '+units['E_grid'],rotation=90)
        ax.set_xlabel('Time t '+units['t'])
        ax.set_ylabel('Position x '+units['r'])
        
    elif projection=='3d':
        ax=fig.add_subplot(projection='3d')
        x_mesh,T_mesh=np.meshgrid(x_grid,t)
        ax.plot_surface(T_mesh,x_mesh,grid_history.T)
    
        if xlim != 'default':
            ax.set_xlim(xlim)
        if tlim != 'default':
            ax.set_ylim(tlim)
        ax.set_xlabel('Time t '+units['t'])
        ax.set_ylabel('Position x '+units['r'])
        ax.set_zlabel(r'Potential $\phi$ '+units['phi_grid'])
        
    ax.set_title('Field history')
    
    return fig, ax

def selected_modes_history_plot(k,t,grid_kt,selected_modes,
                               plot_theoretical_growth_rate,theoretical_growth_rate,units,
                               part='abs',scale='linear',tlim='default',Alim='default'):
    fig=plt.figure(figsize=(8,8))
    if part=='abs':
        grid_kt=2.*np.abs(grid_kt) # Because c_n = 0.5*(a_n-i*b_n)
    elif part=='real':
        grid_kt=2.*np.real(grid_kt) # Because c_n = 0.5*(a_n-i*b_n)
    elif part=='imag':
        grid_kt=2.*np.imag(grid_kt) # Because c_n = 0.5*(a_n-i*b_n)
        
    for index in range(len(selected_modes)):
        ax=fig.add_subplot(len(selected_modes),1,index+1)
        ax.plot(t,grid_kt[selected_modes[index],:],label='Simulation result')
        if plot_theoretical_growth_rate:
            growth_rate_label='Theoretical growth rate = %.3f'\
                               %theoretical_growth_rate[selected_modes[index]]\
                               +' '+units['omega']
            ax.plot(t[:],grid_kt[selected_modes[index],0]
                    *np.exp(theoretical_growth_rate[selected_modes[index]]*t[:]),
                    label=growth_rate_label)
           
            
        title='Mode '+str(selected_modes[index])\
                +', k = '+str(k[selected_modes[index]])+' '+units['k']
        if tlim != 'default':
            ax.set_xlim(tlim)
        if Alim != 'default':
            ax.set_ylim(Alim)
        ax.set_xlabel('Time t '+units['t'])
        ax.set_ylabel('Amplitude '+units['arb. unit'])
        ax.set_yscale(scale)
        ax.legend(loc='best')
        ax.set_title(title)
        
    plt.suptitle('Mode Amplitude History')
    plt.tight_layout()
    return fig, ax

def all_modes_history_plot(k,t,dk,DT,grid_kt,units,
                           part='abs',scale='linear',tlim='default',klim='default'):
    fig,ax=plt.subplots()
    grid_kt=grid_kt[::-1,:] # So that k=0 in in the lowest row
    grid_kt=grid_kt[1:-1,:] # Don't draw k=0 and k= np.pi/dx because these two rows==0,
                            # and taking log of 0 is -inf
    
    if part=='abs':
        grid_kt=2.*np.abs(grid_kt) # Because c_n = 0.5*(a_n-i*b_n)
    elif part=='real':
        grid_kt=2.*np.real(grid_kt) # Because c_n = 0.5*(a_n-i*b_n)
    elif part=='imag':
        grid_kt=2.*np.imag(grid_kt) # Because c_n = 0.5*(a_n-i*b_n)
    if scale=='log':
        grid_kt=np.log10(2.*np.abs(grid_kt))
    
    pos=ax.imshow(grid_kt,extent=(t[0]-0.5*DT,t[-1]+0.5*DT,k[1]-0.5*dk,k[-2]+0.5*dk))
        
    if tlim != 'default':
        ax.set_xlim(tlim)
    if klim != 'default':
        ax.set_ylim(klim)
        
    cbar=fig.colorbar(pos,ax=ax)
    ax.set_xlabel(r'Time $t$ '+units['t'])
    ax.set_ylabel(r'Wave number $k$ '+units['k'])
    
    # ax.set_xlabel('Step #')
    # ax.set_ylabel('Mode number')
    ax.set_title('Mode Amplitue History')
    if scale=='linear':
        cbar.set_label('Amplitude '+units['arb. unit'],rotation=90)
    elif scale=='log':
        cbar.set_label(r'$\log_{10}(Amplitude)$ '+units['arb. unit'],rotation=90)
    plt.tight_layout()
    
    return fig, ax

def tracker_particle_trajectory_plot(t,R,V,NTracker,units,
                                     space=['R-t'],tlim='default',Rlim='default',Vlim='default'):
    fig=plt.figure()
    if 'R-t' in space:
        ax=fig.add_subplot(len(space),1,space.index('R-t')+1)
        for particle in range(NTracker):
            ax.plot(t,R[particle,:],label='Particle '+str(particle+1))
        
        if tlim != 'default':
            ax.set_xlim(tlim)
        if Rlim != 'default':
            ax.set_ylim(Rlim)
        ax.set_xlabel('Time t '+units['t'])
        ax.set_ylabel('Position r '+units['r'])
        ax.legend(loc='best')
        ax.set_title('R-t Plot')
        
    if 'V-t' in space:
        ax=fig.add_subplot(len(space),1,space.index('V-t')+1)
        for particle in range(NTracker):
            ax.plot(t,V[particle,:],label='Particle '+str(particle+1))
        
        if tlim != 'default':
            ax.set_xlim(tlim)
        if Vlim != 'default':
            ax.set_ylim(Vlim)
        ax.set_xlabel('Time t '+units['t'])
        ax.set_ylabel('Velocity v '+units['v'])
        ax.legend(loc='best')
        ax.set_title('V-t Plot')
        
    if 'R-V' in space:
        ax=fig.add_subplot(len(space),1,space.index('R-V')+1)
        for particle in range(NTracker):
            ax.plot(R[particle,:],V[particle,:],label='Particle '+str(particle+1))
        
        if Rlim != 'default':
            ax.set_xlim(Rlim)
        if Vlim != 'default':
            ax.set_ylim(Vlim)
        ax.set_xlabel('Position r '+units['r'])
        ax.set_ylabel('Velocity v'+units['v'])
        ax.legend(loc='best')
        ax.set_title('R-V Plot')
        
    plt.tight_layout()
    
    return fig, ax

# =============================================================================
# dispersion relation
# =============================================================================
def dispersion_relation_plot(k,omega,dk,domega,grid_omegak,units,
                             plot_theoretical_dispersion_relation,theoretical_omega_of_k,
                             part='abs',scale='linear',klim='default',omegalim='default'):
    fig,ax=plt.subplots()
    grid_omegak=grid_omegak.T # So that the horizontal axis is k and the vertical axis is omega
    grid_omegak=grid_omegak[::-1,:] # So that (k = 0, omega = 0) is at bottom left corner
    grid_omegak=grid_omegak[:,1:-1] # Don't draw k=0 and k= np.pi/dx because these two rows==0,
                                    # and taking log of 0 is -inf

    if part=='abs':
        grid_omegak=2.*np.abs(grid_omegak) # Because c_n = 0.5*(a_n-i*b_n)
    elif part=='real':
        grid_omegak=2.*np.real(grid_omegak) # Because c_n = 0.5*(a_n-i*b_n)
    elif part=='imag':
        grid_omegak=2.*np.imag(grid_omegak) # Because c_n = 0.5*(a_n-i*b_n)
    if scale=='log':
        grid_omegak=np.log10(2.*np.abs(grid_omegak))
        
    pos=ax.imshow(grid_omegak,extent=(k[1]-0.5*dk,k[-2]+0.5*dk,
                                      omega[0]-0.5*domega,omega[-1]+0.5*domega))
    cbar=fig.colorbar(pos,ax=ax)
    
    # plot theoretical dispersion relation    
    if plot_theoretical_dispersion_relation:
        for i in range(len(theoretical_omega_of_k)):
            ax.plot(k,theoretical_omega_of_k[i],'red',lw=0.5,label='Theoretical')
        ax.legend(loc='best')
    
    if klim != 'default':
        ax.set_xlim(klim)
    if omegalim != 'default':
        ax.set_ylim(omegalim)
    ax.set_xlabel(r'Wave number $k$ '+units['k'])
    ax.set_ylabel(r'Angular frequency $\omega$ '+units['omega'])
    ax.set_title('Dispersion Relation')
    if scale=='linear':
        cbar.set_label('Amplitude '+units['arb. unit'],rotation=90)
    elif scale=='log':
        cbar.set_label(r'$\log_{10}(Amplitude)$ '+units['arb. unit'],rotation=90)

    return fig, ax 

def phase_space_plot(species_parameters,r,v,units):
    fig,ax=plt.subplots()
    particle_counter=0

    for species in range(len(species_parameters)):
        ax.scatter(r[particle_counter:particle_counter+species_parameters[species]['N']],
                   v[particle_counter:particle_counter+species_parameters[species]['N']],
                   s=1,alpha=0.5,label='species '+str(species+1))
        particle_counter+=species_parameters[species]['N']
    # ax1.hexbin(r,v,gridsize=200)
    
    # Tracker particles
    # ax.scatter(r[N//4],v[N//4],c='black',s=30,lw=1,edgecolor='white') # tracker 1
    # ax.scatter(r[-1],v[-1],c='white',s=30,lw=1,edgecolor='black') # tracker 2
    
    ax.set_xlabel('Position x '+units['r'])
    ax.set_ylabel('Velocity v '+units['v'])
    ax.legend(loc='best')
    ax.set_title('Phase Space')
    ax.set_xlim()
    ax.set_ylim()
    
    return fig, ax

def distribution_function_grid_plot(X,V,distribution_function,projection,units):
    fig=plt.figure()
    if projection=='2d':
        ax=fig.add_subplot()
        ax.imshow(distribution_function)
    if projection=='3d':
        ax=fig.add_subplot(projection=='3d')
        ax.plot_surface(X,V,distribution_function.T)
    ax.set_xlabel('position x '+units['r'])
    ax.set_ylabel('velocity v '+units['v'])
    ax.set_title('Phase Space')
    
    return fig, ax

def grid_plot(x,grid,units):
    fig,ax=plt.subplots()
    ax.plot(x,grid)
    ax.set_xlabel('position x '+units['r'])
    ax.set_ylabel(r' $\$ '+units[''])
    ax.set_title('')
    return fig, ax

def history_frequency_plot(omega,FFTspectrum,units):
    fig,ax=plt.subplots()
    ax.plot(omega,FFTspectrum)
    ax.set_xlabel(r'Frequency $omega$ '+units['omega'])
    ax.set_ylabel('Amplitude '+units['arb. units'])
    ax.set_title('')
    return fig,ax

# =============================================================================
# 
# =============================================================================
if __name__=='__main__':
    import read_output