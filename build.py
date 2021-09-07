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

def load_particles(input_txt_parameters,species_parameters):
    '''
    For generating the particle m, q, r, v np arrays from the parameters dictionary.

    Parameters
    ----------
    input_txt_parameters : dict
        Parameters read from `input.txt` using input_output.load_input(input_file_loc).
    species_parameters: list
        List of dicts, where each dict contains parameters for each particle species

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
    m_e=input_txt_parameters['m_e']
    m_p=input_txt_parameters['m_p']
    e=input_txt_parameters['e']
    k_B=input_txt_parameters['k_B']
    
    L=input_txt_parameters['L']
    
    InitialCondition=input_txt_parameters['InitialCondition']
    
    k=2.*np.pi/L  # First (lowest) mode in the system (smallest k)
    # density=N/L
    
    m=np.array([])
    q=np.array([])
    r=np.array([])
    v=np.array([])
    for species in range(len(species_parameters)):
        # Assign mass and charge to each particle
        N=species_parameters[species]['N']
        if species_parameters[species]['m']=='m_e':
            m_species=m_e*np.ones(N)
        elif species_parameters[species]['m']=='m_p':
            m_species=m_p*np.ones(N)
        else:
            m_species=species_parameters[species]['m']*np.ones(N)
            
        if species_parameters[species]['q']=='e':
            q_species=e*np.ones(N)
        elif species_parameters[species]['q']=='-e':
            q_species=-e*np.ones(N)
        else:
            q_species=species_parameters[species]['q']*np.ones(N)
        # place particles in x-v phase space
        # distribution      = Maxwellian # Used for choosing different initial velocity distribution function
        v0         = species_parameters[species]['v0']
        v_sigma    = species_parameters[species]['v_sigma']        # Velocity width (v0_sigma**2 ~ Temperature)
        T          = species_parameters[species]['T']        # Temperature for a Maxwellian velocity distribution
        
        r_species=np.linspace(0.,L,N,endpoint=False)
        v_species=np.random.normal(v0,v_sigma,size=N)
        
        # excite mode in modes list
        modes               = species_parameters[species]['Modes']     
        X1                  = species_parameters[species]['X1']
        V1                  = species_parameters[species]['V1']
        THETAX              = species_parameters[species]['THETAX']
        THETAV              = species_parameters[species]['THETAV']
        
        for mode in modes:
            r_species+=X1*np.cos(mode*k*r_species+THETAX)
            v_species+=V1*np.sin(mode*k*v_species+THETAV)
            
        # Periodic boundary condition
        r_species=r_species%L
        
        # Append to the final m, q, r, v array
        m=np.append(m,m_species)
        q=np.append(q,q_species)
        r=np.append(r,r_species)
        v=np.append(v,v_species)
 
    # if InitialCondition=='Particle Pair Oscillation':
    #     N=2
    #     m=np.array([m_e,m_e])
    #     q=np.array([e,-e])
    #     r=np.array([L/4.,2*L/4.])
    #     v=np.array([0.,0.])
    
    # if InitialCondition=='Two_Stream_Instability':
    #     v0         = parameters['v0']     # velocity
    #     v0_sigma   = parameters['v0_sigma']       # velocity width
    #     charge     = parameters['charge']
    #     A          = parameters['A']
    #     modes      = parameters['Modes']
        
    #     m=np.ones(N)*m_e 
    #     q=-1*np.ones(N)*e
    #     if charge == 'opposite':
    #         q[N//2:]*=-1.
    #     r=np.random.random(N)*L
    #     v=np.append(np.random.normal(v0,v0_sigma,size=(1,N//2)),
    #                 np.random.normal(-v0,v0_sigma,size=(1,N//2)))
        
    #     for mode in modes:
    #         r[:N//2]+=A*mode/density*np.sin(mode*k*r[:N//2])
    #     r=r%L
        
    # if InitialCondition=='Single Beam':
    #     v0         = 1.0        # velocity
    #     v0_sigma   = 0.0        # velocity width
    #     A          = 0.01       # sinusoidal perturbation amplitude
        
    #     m=np.ones(N)*m_e
    #     q=-1*np.ones(N)*e
    #     # q[N//2:]*=-1.
    #     r=np.random.random(N)*L
    #     # v=np.random.normal(v0,v0_sigma,size=(1,N))[0]
    #     # v[N//2:]*=-1.
    #     v=np.random.normal(v0,v0_sigma,size=(1,N))[0]
    
    # if InitialCondition=='Self Defined':
    #     v0         = parameters['v0']        # velocity
    #     v0_sigma   = parameters['v0_sigma']        # velocity width
    #     charge     = parameters['charge']
    #     A          = parameters['A']
    #     modes  = parameters['Modes']
        
    #     m=np.ones(N)*m_e 
    #     q=-1*np.ones(N)*e
    #     if charge == 'opposite':
    #         q[N//2:]*=-1.
    #     r=np.random.random(N)*L
    #     v=np.append(np.random.normal(v0,v0_sigma,size=(1,N//2)),
    #                 np.random.normal(-v0,v0_sigma,size=(1,N//2)))
        
    #     for mode in modes:
    #         r[:N//2]+=A*mode/density*np.sin(k[mode]*r[:N//2])
    #     r=r%L
        
    return m,q,r,v

# =============================================================================
# 
# =============================================================================
if __name__=='__main__':
    import input_output 
    # input_folder_loc=sys.argv[1] # For running in terminal
    input_folder_loc='simulation_results/EXAMPLE/inputs' # for running in IDE (ex:Spyder)
    input_txt_parameters, save_dir, species_parameters = input_output.load_input(input_folder_loc)
    
    m,q,r,v=load_particles(input_txt_parameters,species_parameters)
    
    from analysis_and_plotting import analysis
    UnitSystem          = input_txt_parameters['UnitSystem']
    units=analysis.generate_units(UnitSystem)
    
    from analysis_and_plotting import plotting
    plotting.phase_space_plot(species_parameters,r,v,units)
    