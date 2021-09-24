#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 02:31:43 2021

@author: mac
"""
import numpy as np
import matplotlib.pyplot as plt
import sys

import input_output
import build
import PIC_functions

from analysis_and_plotting import theory
from analysis_and_plotting import animation
from analysis_and_plotting import analysis
from analysis_and_plotting import plotting

# =============================================================================
# Load input file into parameters dict & save to same dir as input file
# =============================================================================
# input_folder_loc=sys.argv[1] # For running in terminal
input_folder_loc='example_problems/cold_plasma_oscillation/inputs' # for running in IDE (ex:Spyder)
input_txt_parameters, save_dir, species_parameters = input_output.load_input(input_folder_loc)

# =============================================================================
# Load from parameters dict
# =============================================================================
UnitSystem          = input_txt_parameters['UnitSystem']

e                   = input_txt_parameters['e']     
epsilon0            = input_txt_parameters['epsilon0']       
m_e                 = input_txt_parameters['m_e']       
m_p                 = input_txt_parameters['m_p']     
k_B                 = input_txt_parameters['k_B'] 

L                   = input_txt_parameters['L']       
NT                  = input_txt_parameters['NT']       

# density             = input_txt_parameters['density']       # not sure to input or calculate

InitialCondition    = input_txt_parameters['InitialCondition']

NG                  = input_txt_parameters['NG']    
DT                  = input_txt_parameters['DT']      

IW                  = input_txt_parameters['IW']
scheme              = input_txt_parameters['scheme']
solver              = input_txt_parameters['solver']

NTracker            = input_txt_parameters['NTracker']       

save_output                      = input_txt_parameters['save_output']   
output_Ndecimal_places           = input_txt_parameters['output_Ndecimal_places']   

save_animation                   = input_txt_parameters['save_animation']    
save_energy_history              = input_txt_parameters['save_energy_history']          
save_momentum_change_history     = input_txt_parameters['save_momentum_change_history']    
save_grid_history                = input_txt_parameters['save_grid_history']    
save_selected_modes_history      = input_txt_parameters['save_selected_modes_history']    
save_all_modes_history           = input_txt_parameters['save_all_modes_history']    
save_tracker_particle_trajectory = input_txt_parameters['save_tracker_particle_trajectory']    
save_omegak                      = input_txt_parameters['save_omegak']    

# snapshot animations
plot_animation                       = input_txt_parameters['plot_animation']   

# energy-time history 
plot_energy_history                  = input_txt_parameters['plot_energy_history']   
energy_history_tlim                  = input_txt_parameters['energy_history_tlim']   
energy_history_Elim                  = input_txt_parameters['energy_history_Elim']   

# momentum_change-time history 
plot_momentum_change_history         = input_txt_parameters['plot_momentum_change_history']   
momentum_change_history_tlim         = input_txt_parameters['momentum_change_history_tlim']   
momentum_change_history_Plim         = input_txt_parameters['momentum_change_history_Plim']   

# grid history 
plot_grid_history                    = input_txt_parameters['plot_grid_history']
grid_history_projection              = input_txt_parameters['grid_history_projection']
grid_history_xlim                    = input_txt_parameters['grid_history_xlim']
grid_history_tlim                    = input_txt_parameters['grid_history_tlim']

# selected_mode history 
plot_selected_modes_history          = input_txt_parameters['plot_selected_modes_history']
selected_modes                       = input_txt_parameters['selected_modes']
plot_theoretical_growth_rate         = input_txt_parameters['plot_theoretical_growth_rate']
selected_modes_part                  = input_txt_parameters['selected_modes_part']
selected_modes_scale                 = input_txt_parameters['selected_modes_scale']
selected_mode_history_tlim           = input_txt_parameters['selected_mode_history_tlim']
selected_mode_history_Alim           = input_txt_parameters['selected_mode_history_Alim']

# all_modes history 
plot_all_modes_history               = input_txt_parameters['plot_all_modes_history']
all_modes_part                       = input_txt_parameters['all_modes_part']
all_modes_scale                      = input_txt_parameters['all_modes_scale']
all_mode_history_tlim                = input_txt_parameters['all_mode_history_tlim']
all_mode_history_klim                = input_txt_parameters['all_mode_history_klim']

# tracker_particle_trajectory_history
plot_tracker_particle_trajectory     = input_txt_parameters['plot_tracker_particle_trajectory']   
tracker_particle_trajectory_space    = input_txt_parameters['tracker_particle_trajectory_space']    
tracker_particle_trajectory_tlim     = input_txt_parameters['tracker_particle_trajectory_tlim']   
tracker_particle_trajectory_Rlim     = input_txt_parameters['tracker_particle_trajectory_Rlim']   
tracker_particle_trajectory_Vlim     = input_txt_parameters['tracker_particle_trajectory_Vlim']   

# dispersion relation
plot_omegak                          = input_txt_parameters['plot_omegak']   
plot_theoretical_dispersion_relation = input_txt_parameters['plot_theoretical_dispersion_relation'] 
dispersion_relation_part             = input_txt_parameters['dispersion_relation_part']
dispersion_relation_scale            = input_txt_parameters['dispersion_relation_scale']  
dispersion_relation_klim             = input_txt_parameters['dispersion_relation_klim']   
dispersion_relation_omegalim         = input_txt_parameters['dispersion_relation_omegalim']   

# =============================================================================
# Units for physical quantities
# =============================================================================
units=analysis.generate_units(UnitSystem)

# =============================================================================
# Load particles
# =============================================================================
m,q,r,v=build.load_particles(input_txt_parameters,species_parameters)
N=len(m) # Total number of particles

# =============================================================================
# Derived parameters 
# =============================================================================
# density             = N/L      # not sure to input or calculate

dx                  = L/NG # grid spacing
x_grid                   = np.arange(0.,L,dx) # np array for x axis
k                   = 2.*np.pi*np.fft.rfftfreq(NG,dx)
dk                  = k[1]-k[0]

t                   = np.arange(0.,NT*DT,DT) # np array for time axis
omega               = 2.*np.pi*np.fft.rfftfreq(len(t),DT)
domega              = omega[1]-omega[0]

# =============================================================================
# Plasma parameters
# =============================================================================
omega_plasma        = np.sqrt(((N/L)*e**2)/(epsilon0*m_e)) # plasma frequency

print(('plasma frequency ω_p = %.4f '+units['omega']
        +',\nω_p*DT = %.4f, total steps NT = %i\n')
      %(omega_plasma,omega_plasma*DT,NT))

# lambda_debye        = np.sqrt((epsilon0*k_B*T)/((N/(2*L)*e**2))) # Deybe length
# N_D                 = N/(2*L)*lambda_debye
# print(('Debye length = %.4f '+units['r']
#       +',\nDebye length/L = %.4f, number of particles in one Debye length = %.4f\n')
#       %(lambda_debye,lambda_debye/L,N_D))

# =============================================================================
# diagnostics physical quantites
# =============================================================================
Tracker_index = np.random.choice(N,NTracker,replace=False)
R             = np.zeros((NTracker,len(t))) # tracker particle phase space trajectory
V             = np.zeros((NTracker,len(t)))

P     = np.zeros(len(t)) # total momentum
P_abs = np.zeros(len(t)) # Sum of magnitudes of momentum for each particle

E_D     = np.zeros(len(t)) # drift kinetic energy
E_T     = np.zeros(len(t)) # thermal kinetic energy
E_F     = np.zeros(len(t)) # field energy
E_total = np.zeros(len(t))

# shape=(NG,len(t)) , field(x,t) = grid_history[int(x//dx)%NG,int(t//DT)%NT]
rho_grid_history = np.zeros((NG,len(t)))
phi_grid_history = np.zeros((NG,len(t)))
E_grid_history   = np.zeros((NG,len(t)))

    
# =============================================================================
# Theory
# =============================================================================
theoretical_omega_of_k,theoretical_growth_rate=theory.dispersion_relation(input_txt_parameters,
                                                                          species_parameters)

# =============================================================================
# Analysis    
# =============================================================================
grid_kt=analysis.grid_xt_to_kt(E_grid_history)
grid_omegak=analysis.grid_xt_to_omegak(k,E_grid_history)

# =============================================================================
# Plotting
# =============================================================================
if plot_energy_history:
    fig1=plotting.energy_history_plot(t,E_D,E_T,E_F,units,
                                      tlim=energy_history_tlim,Elim=energy_history_Elim)

# momentum_change-time history 
if plot_momentum_change_history:
    fig2=plotting.momentum_change_history_plot(t, (P-P[0])/P_abs, units,
                                               tlim=momentum_change_history_tlim,
                                               Plim=momentum_change_history_Plim)

# grid history 
if plot_grid_history:
    fig3=plotting.grid_history_plot(x_grid, t, dx, DT, E_grid_history, 
                                    units,projection=grid_history_projection,
                                    xlim=grid_history_xlim,tlim=grid_history_tlim)
    
# selected_mode history 
if plot_selected_modes_history:
    fig4=plotting.selected_modes_history_plot(k, t, grid_kt, selected_modes, 
                                              plot_theoretical_growth_rate, 
                                              theoretical_growth_rate, units,
                                              part=selected_modes_part,
                                              scale=selected_modes_scale,
                                              tlim=selected_mode_history_tlim,
                                              Alim=selected_mode_history_Alim)

# all_modes history 
if plot_all_modes_history:
    fig5=plotting.all_modes_history_plot(k, t, dk, DT, grid_kt, units,
                                         part=all_modes_part,scale=all_modes_scale,
                                         tlim=all_mode_history_tlim,
                                         klim=all_mode_history_klim)
    
# tracker_particle_trajectory_history
if plot_tracker_particle_trajectory:
    fig6=plotting.tracker_particle_trajectory_plot(t, R, V, Tracker_index, units,
                                                   space=tracker_particle_trajectory_space,
                                                   tlim=tracker_particle_trajectory_tlim,
                                                   Rlim=tracker_particle_trajectory_Rlim,
                                                   Vlim=tracker_particle_trajectory_Vlim)
    
# dispersion relation
if plot_omegak:
    fig7=plotting.dispersion_relation_plot(k, omega, dk, domega, 
                                           grid_omegak, units, 
                                           plot_theoretical_dispersion_relation, 
                                           theoretical_omega_of_k,
                                           part=dispersion_relation_part,
                                           scale=dispersion_relation_scale,
                                           klim=dispersion_relation_klim,
                                           omegalim=dispersion_relation_omegalim)

plt.show()