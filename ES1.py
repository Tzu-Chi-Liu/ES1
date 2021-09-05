import numpy as np
import matplotlib.pyplot as plt
import sys
import os

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
# input_file_loc=sys.argv[1] # For running in terminal
input_file_loc='simulation_results/cold_plasma_oscillation/EXAMPLE/input.txt' # for running in IDE (ex:Spyder)
parameters, save_dir = input_output.load_input(input_file_loc)

# =============================================================================
# Load from parameters dict
# =============================================================================
UnitSystem          = parameters['UnitSystem']
e                   = parameters['e']     # elementary charge
epsilon0            = parameters['epsilon0']        # vacuum permittivity
m_e                 = parameters['m_e']        # electron mass
m_p                 = parameters['m_p']     # proton mass
k_B                 = parameters['k_B'] 

N                   = parameters['N']      # number of particles
L                   = parameters['L']        # length of box

density             = parameters['density']       # not sure to input or calculate

InitialCondition    = parameters['InitialCondition']

NG                  = parameters['NG']         # number of grids

NT                  = parameters['NT']       # Total simulation time
DT                  = parameters['DT']       # Time step

IW                  = parameters['IW']
scheme              = parameters['scheme']
solver              = parameters['solver']

NTracker            = parameters['NTracker']          # Number of tracker particles

# .txt output files
save_output         = parameters['save_output']
output_decimal_places = parameters['output_Ndecimal_places']

# 
save_animation      = parameters['save_animation']     # save animation of snapshots
save_omegak         = parameters['save_omegak']       # save dispersion relation

# snapshot animations
plot_animation      = parameters['plot_animation']  # Plot animation of snapshots

# energy-time history 
plot_energy_history = parameters['plot_energy_history']           # Plot history
# energy_history_tlim= 0.,1.      # Separate by comma
# energy_history_Elim= 0.,1.

# momentum_change-time history 
plot_momentum_change_history = parameters['plot_momentum_change_history'] # Plot history
# momentum_change_history_tlim= 0.,1.      # Separate by comma
# momentum_change_history_Elim= 0.,1.

# grid history 
plot_grid_history = parameters['plot_grid_history']           # Plot history
# grid_history_projection = 3d
# grid_history_xlim= 0.,1.      # Separate by comma
# grid_history_tlim= 0.,1.

# selected_mode history 
plot_selected_modes_history = parameters['plot_selected_modes_history']           # Plot history
selected_modes = parameters['selected_modes']
plot_theoretical_growth_rate = parameters['plot_theoretical_growth_rate']

# all_mode history 
plot_all_modes_history = parameters['plot_all_modes_history']           # Plot history

# tracker_particle_trajectory_history
plot_trajectory     = parameters['plot_trajectory']  # Trajectory of tracker particle(s)
# space='R-t','V-t','R-V'      
# tlim=None
# Rlim=None
# Vlim=None

# dispersion relation
plot_omegak         = parameters['plot_omegak']       # Plot dispersion relation
plot_theoretical_dispersion_relation = parameters['plot_theoretical_dispersion_relation']

# =============================================================================
# Units for physical quantities
# =============================================================================
units=analysis.generate_units(UnitSystem)

# =============================================================================
# Load particles
# =============================================================================
m,q,r,v=build.load_particles(parameters)
    
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

omega_plasma        = np.sqrt(((N/(2*L))*e**2)/(epsilon0*m_e)) # plasma frequency

print(('Assuming half of N are electrons, other half are protons,\n'
       +'plasma frequency = %.4f '+units['omega']
       +',\nomega_plasma*DT = %.4f, total steps = %i\n')
      %(omega_plasma,omega_plasma*DT,NT))

# lambda_debye        = np.sqrt((epsilon0*k_B*T)/((N/(2*L)*e**2))) # Deybe length
# N_D                 = N/(2*L)*lambda_debye
# print(('Assuming half of N are electrons, other half are protons,\n'+
#       'Debye length = %.4f '+units['r']
#       +',\nDebye length/L = %.4f, number of particles in one Debye length = %.4f\n')
#       %(lambda_debye,lambda_debye/L,N_D))

# =============================================================================
# diagnostics physical quantites
# =============================================================================
R     = np.zeros((NTracker,len(t))) # tracker particle phase space trajectory
V     = np.zeros((NTracker,len(t)))

P     = np.zeros(len(t)) # total momentum
E_D   = np.zeros(len(t)) # drift kinetic energy
E_T   = np.zeros(len(t)) # thermal kinetic energy
E_F   = np.zeros(len(t)) # field energy

# shape=(NG,len(t)) , field(x,t) = grid_history[int(x//dx)%NG,int(t//DT)%NT]
rho_grid_history = np.zeros((NG,len(t)))
phi_grid_history = np.zeros((NG,len(t)))
E_grid_history   = np.zeros((NG,len(t)))

# =============================================================================
# external fields
# =============================================================================
def external_field(x_grid,t):
    phiX=np.zeros(NG)
    EX=np.zeros(NG)
    
    return phiX,EX
            
# =============================================================================
# main simulation loop
# =============================================================================
for step in range(len(t)):
    # solve field
    rho_grid=PIC_functions.weightrho(r,q,IW,N,NG,dx)
    phi_grid,E_grid=PIC_functions.solve_poisson(rho_grid,NG,dx,epsilon0,
                                               solver='FFT',operator='local')
    
    # copy old position & velocity
    r_old=np.copy(r)
    v_old=np.copy(v)
    
    # update particle i velocity v[i] and position r[i] at time t
    for i in range(N):
        E_par=PIC_functions.weightField(r[i],E_grid,IW,NG,dx)
        r[i],v[i]=PIC_functions.update_motion(m[i],q[i],r[i],v[i],E_par,
                                              scheme,IW,DT,L,E_grid,step)
        
    # tracker particle trajectory & histories
    # for i in np.random.choice(N,NTracker,replace=False):
    #     R[i,step]=r[i]
    #     V[i,step]=v[i]
    P[step]=np.sum(m*0.5*(v_old+v))
    E_F[step]=0.5*np.dot(phi_grid*dx,rho_grid)
    E_D[step]=np.sum(0.5*m*v*v_old)
    rho_grid_history[:,step]=rho_grid
    phi_grid_history[:,step]=phi_grid
    E_grid_history[:,step]=E_grid
        
    if plot_animation:
        NSP=2
        pause_time=0.01
        animation.diagnostics_animation(t,NSP,N,r,0.5*(v_old+v),
                                        NG,x_grid,rho_grid,phi_grid,E_grid,
                                        step,units,pause_time,save_dir)
    
    # save output
    if save_output:
        input_output.output_to_file(step,t,N,m,q,r,v,NG,x_grid,rho_grid,phi_grid,E_grid,
                                    InitialCondition,save_dir,units,output_decimal_places)
    
    print(('step = %d '+units['t']
           +', P = %+.3e'+units['Momentum']
           +', E = %+.3e'+units['Energy'])
          %(step,P[step],E_D[step]+E_T[step]+E_F[step]))
    
# =============================================================================
# Theory
# =============================================================================
theoretical_omega_of_k,theoretical_growth_rate=theory.dispersion_relation(parameters)

# =============================================================================
# Analysis    
# =============================================================================
grid_kt=analysis.grid_xt_to_kt(E_grid_history)
grid_omegak=analysis.grid_xt_to_omegak(k,E_grid_history)

# =============================================================================
# Plotting
# =============================================================================
if plot_energy_history:
    fig1=plotting.energy_history_plot(t,E_D,E_T,E_F,units)

# momentum_change-time history 
if plot_momentum_change_history:
    fig2=plotting.momentum_change_history_plot(t, P, units)

# grid history 
if plot_grid_history:
    fig3=plotting.grid_history_plot(x_grid, t, dx, DT, E_grid_history, units)
    
# selected_mode history 
if plot_selected_modes_history:
    fig4=plotting.selected_mode_history_plot(k, t, grid_kt, selected_modes, 
                                             plot_theoretical_growth_rate, 
                                             theoretical_growth_rate, units)

# all_modes history 
if plot_all_modes_history:
    fig5=plotting.all_mode_history_plot(k, t, dk, DT, grid_kt, units)
    
# tracker_particle_trajectory_history
if plot_trajectory:
    fig6=plotting.tracker_particle_trajectory_plot(t, R, V, NTracker, units)
    
# dispersion relation
if plot_omegak:
    fig7=plotting.dispersion_relation_plot(k, omega, dk, domega, 
                                           grid_omegak, units, 
                                           plot_theoretical_dispersion_relation, 
                                           theoretical_omega_of_k)
