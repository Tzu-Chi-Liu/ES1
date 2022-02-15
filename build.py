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
    For generating the t=0. particle m, q, r, v np arrays from the species_parameters.

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
    # Load constants
    m_e=input_txt_parameters['m_e']
    m_p=input_txt_parameters['m_p']
    e=input_txt_parameters['e']
    k_B=input_txt_parameters['k_B']
    
    L=input_txt_parameters['L']
    
    k=2.*np.pi/L  # First (lowest) mode in the system (smallest k)
    # density=N/L
    
    # Output arrays
    m=np.array([])
    q=np.array([])
    r=np.array([])
    v=np.array([])
    
    for species in range(len(species_parameters)):
        # Mass and charge of each species
        N=species_parameters[species]['N'] # Number of particles in the species
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
            
        # Velocity distribution function parameters for each species
        # distribution      = Maxwellian 
        v0         = species_parameters[species]['v0']
        v_sigma    = species_parameters[species]['v_sigma']      
        # T          = species_parameters[species]['T']       
        # v_sigma    = np.sqrt(k_B*T/m_species)
        
        # place particles in x-v phase space
        r_species=np.linspace(0.,L,N,endpoint=False)
        # r_species=np.random.rand(N)*L
        v_species=np.random.normal(v0,v_sigma,size=N)
        
        # excite mode in 'modes' list
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
        m=np.append(m,m_species*np.ones(N))
        q=np.append(q,q_species*np.ones(N))
        r=np.append(r,r_species)
        v=np.append(v,v_species)
    
    return m,q,r,v

# =============================================================================
# 
# =============================================================================
if __name__=='__main__':
    import input_output 
    # input_folder_loc=sys.argv[1] # For running in terminal
    input_folder_loc='simulation_results/Numeric_Instability/Single_Cold_Stream/untitled folder 4/inputs' # for running in IDE (ex:Spyder)
    input_txt_parameters, save_dir, species_parameters = input_output.load_input(input_folder_loc)
    
    m,q,r,v=load_particles(input_txt_parameters,species_parameters)
    
    import PIC_functions
    N = len(m) # Total number of particles
    epsilon0            = input_txt_parameters['epsilon0']  
    L                   = input_txt_parameters['L']      
    NG                  = input_txt_parameters['NG']    
    dx                  = L/NG # grid spacing
    IW                  = input_txt_parameters['IW']
    rho_grid = PIC_functions.weightrho(r, q, IW, N, NG, dx)
    phi_grid, E_grid = PIC_functions.solve_poisson(rho_grid, NG, dx, epsilon0,
                                               solver='FFT', operator='local')
    
    from analysis_and_plotting import analysis
    UnitSystem          = input_txt_parameters['UnitSystem']
    units=analysis.generate_units(UnitSystem)
    
    from analysis_and_plotting import plotting
    x_grid              = np.arange(0., L, dx) # np array for x axis
    plotting.phase_space_plot(species_parameters,r,v,units)
    plotting.grid_plot(x_grid,rho_grid,units,'density (arb. units)', 'density')
    plotting.grid_plot(x_grid,phi_grid,units,'potential (arb. units)', 'potential')
    plotting.grid_plot(x_grid,E_grid,units,'field (arb. units)', 'field')
    
    import matplotlib.pyplot as plt
    plt.figure()
    plt.plot(np.abs(np.fft.rfft(phi_grid)))
    