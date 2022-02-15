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
# input_folder_loc = sys.argv[1] # For running in terminal
input_folder_loc = 'simulation_results/Numeric_Instability/Single_Cold_Stream/untitled folder 2/inputs' # for running in IDE (ex:Spyder)
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

plot_animation                       = input_txt_parameters['plot_animation']   

plot_energy_history                  = input_txt_parameters['plot_energy_history']   
energy_history_tlim                  = input_txt_parameters['energy_history_tlim']   
energy_history_Elim                  = input_txt_parameters['energy_history_Elim']   

plot_momentum_change_history         = input_txt_parameters['plot_momentum_change_history']   
momentum_change_history_tlim         = input_txt_parameters['momentum_change_history_tlim']   
momentum_change_history_Plim         = input_txt_parameters['momentum_change_history_Plim']   

plot_grid_history                    = input_txt_parameters['plot_grid_history']
grid_history_projection              = input_txt_parameters['grid_history_projection']
grid_history_xlim                    = input_txt_parameters['grid_history_xlim']
grid_history_tlim                    = input_txt_parameters['grid_history_tlim']

plot_selected_modes_history          = input_txt_parameters['plot_selected_modes_history']
selected_modes                       = input_txt_parameters['selected_modes']
plot_theoretical_growth_rate         = input_txt_parameters['plot_theoretical_growth_rate']
selected_modes_part                  = input_txt_parameters['selected_modes_part']
selected_modes_scale                 = input_txt_parameters['selected_modes_scale']
selected_mode_history_tlim           = input_txt_parameters['selected_mode_history_tlim']
selected_mode_history_Alim           = input_txt_parameters['selected_mode_history_Alim']

plot_all_modes_history               = input_txt_parameters['plot_all_modes_history']
all_modes_part                       = input_txt_parameters['all_modes_part']
all_modes_scale                      = input_txt_parameters['all_modes_scale']
all_mode_history_tlim                = input_txt_parameters['all_mode_history_tlim']
all_mode_history_klim                = input_txt_parameters['all_mode_history_klim']

plot_tracker_particle_trajectory     = input_txt_parameters['plot_tracker_particle_trajectory']   
tracker_particle_trajectory_space    = input_txt_parameters['tracker_particle_trajectory_space']    
tracker_particle_trajectory_tlim     = input_txt_parameters['tracker_particle_trajectory_tlim']   
tracker_particle_trajectory_Rlim     = input_txt_parameters['tracker_particle_trajectory_Rlim']   
tracker_particle_trajectory_Vlim     = input_txt_parameters['tracker_particle_trajectory_Vlim']   

plot_omegak                          = input_txt_parameters['plot_omegak']   
plot_theoretical_dispersion_relation = input_txt_parameters['plot_theoretical_dispersion_relation'] 
dispersion_relation_part             = input_txt_parameters['dispersion_relation_part']
dispersion_relation_scale            = input_txt_parameters['dispersion_relation_scale']  
dispersion_relation_klim             = input_txt_parameters['dispersion_relation_klim']   
dispersion_relation_omegalim         = input_txt_parameters['dispersion_relation_omegalim']   

# =============================================================================
# Units for physical quantities
# =============================================================================
units = analysis.generate_units(UnitSystem)

# =============================================================================
# Load particles
# =============================================================================
m, q, r, v = build.load_particles(input_txt_parameters, species_parameters)
N = len(m) # Total number of particles

# =============================================================================
# Derived parameters 
# =============================================================================
dx                  = L/NG # grid spacing
x_grid              = np.arange(0., L, dx) # np array for x axis
t                   = np.linspace(0., NT*DT, NT+1, endpoint=True) # np array for time axis

# =============================================================================
# Plasma parameters
# =============================================================================
omega_plasma_species, lambda_Debye_species = analysis.plasma_parameters(input_txt_parameters,
                                                                        species_parameters, v, v)

print('In the %d species: '%len(species_parameters))
print(('Highest plasma frequency ω_p = %.4f' + units['omega'] + ', ω_p*DT = %.4f')
      %(np.max(omega_plasma_species), np.max(omega_plasma_species)*DT))
print(('Shorest Debye length  λ_D = %.4f' + units['r'] + ',  λ_D/dx = %.4f')
      %(np.max(lambda_Debye_species), np.max(lambda_Debye_species)/dx))
print(('Number of particles in one shortest Debye length = %.4f, λ_D/L = %.4f\n')
      %(N/L*np.max(lambda_Debye_species),np.max(lambda_Debye_species)/L))

# =============================================================================
# diagnostics physical quantites
# =============================================================================
Tracker_index = np.random.choice(N,NTracker, replace=False)
R             = np.zeros((NTracker, len(t))) # tracker particle phase space trajectory
V             = np.zeros((NTracker, len(t)))

P     = np.zeros(len(t)) # total momentum
P_abs = np.zeros(len(t)) # Sum of magnitudes of momentum for each particle

E_D     = np.zeros(len(t)) # drift kinetic energy
E_T     = np.zeros(len(t)) # thermal kinetic energy
E_F     = np.zeros(len(t)) # field energy
E_total = np.zeros(len(t))

# shape=(NG,len(t)) , field(x,t) = grid_history[int(x//dx)%NG,int(t//DT)%NT]
rho_grid_history = np.zeros((NG, len(t)))
phi_grid_history = np.zeros((NG, len(t)))
E_grid_history   = np.zeros((NG, len(t)))
            
omega_plasma_history = np.zeros((len(species_parameters), len(t)))
lambda_debye_history = np.zeros((len(species_parameters), len(t)))
# =============================================================================
# main simulation loop
# =============================================================================
print('Total steps NT = %d, simulation starting...'%NT)
for step in range(len(t)):
    # solve field
    rho_grid = PIC_functions.weightrho(r, q, IW, N, NG, dx)
    phi_grid, E_grid = PIC_functions.solve_poisson(rho_grid, NG, dx, epsilon0,
                                               solver='FFT', operator='local')
    
    # copy old position & velocity
    r_old = np.copy(r)
    v_old = np.copy(v)
    
    # update particle i velocity v[i] and position r[i] at time t
    for i in range(N):
        E_par = PIC_functions.weightField(r[i], E_grid, IW, NG, dx)
        r[i],v[i]=PIC_functions.update_motion(m[i], q[i], r[i], v[i], E_par,
                                              scheme, IW, DT, L, E_grid, step)
        
    # Energy histories
    E_F[step], E_D[step], E_T[step], E_total[step] = analysis.energy(dx, NG, phi_grid, rho_grid, 
                                                                     species_parameters, 
                                                                     m, v, v_old)
        

    # tracker particle trajectory & histories
    for i in range(len(Tracker_index)):
        R[i, step] = r[Tracker_index[i]]
        V[i, step] = 0.5*(v[Tracker_index[i]]+v_old[Tracker_index[i]])
        
    # Total momentum history
    P[step] = np.sum(m*0.5*(v_old+v))
    P_abs[step] = np.sum(m*0.5*np.abs(v_old+v))
    
    # Record grid histories
    rho_grid_history[:, step] = rho_grid
    phi_grid_history[:, step] = phi_grid
    E_grid_history[:, step]   = E_grid
    
    # Record plasma parameters histories for each species
    omega_plasma, lambda_debye = analysis.plasma_parameters(input_txt_parameters, 
                                                            species_parameters,
                                                            v, v_old) 
    omega_plasma_history[:, step] = omega_plasma
    lambda_debye_history[:, step] = lambda_debye
    
    # Plot snapshot animation during simulation
    if plot_animation:
        NSP=2
        pause_time=0.0
        animation.diagnostics_animation(t,NSP,N,r,0.5*(v_old+v),
                                        NG,x_grid,rho_grid,phi_grid,E_grid,
                                        step,units,pause_time,save_dir)
        
    # save output
    if save_output:
        input_output.output_to_file(step,t,N,m,q,r,0.5*(v+v_old),
                                    NG,x_grid,rho_grid,phi_grid,E_grid,
                                    InitialCondition,save_dir,units,output_Ndecimal_places)
    
    # print current step, total momentum and total energy
    print(('step = %d '+units['t']
           +', P = %+.3e'+units['Momentum']
           +', E = %+.3e'+units['Energy'])
          %(step,P[step],E_total[step]))
    
print('Simulation complete')

# =============================================================================
# Theory
# =============================================================================
theoretical_omega_of_k,theoretical_growth_rate=theory.dispersion_relation(input_txt_parameters,
                                                                          species_parameters)

# =============================================================================
# Analysis    
# =============================================================================
k                   = 2.*np.pi*np.fft.rfftfreq(NG,dx)
dk                  = k[1]-k[0]
omega               = 2.*np.pi*np.fft.rfftfreq(len(t),DT)
domega              = omega[1]-omega[0]

grid_kt=analysis.grid_xt_to_kt(E_grid_history)
grid_omegak=analysis.grid_xt_to_omegak(k,E_grid_history)

# =============================================================================
# Plotting
# =============================================================================
if plot_energy_history:
    fig1=plotting.energy_history_plot(t,E_D,E_T,E_F,units,
                                      tlim=energy_history_tlim,Elim=energy_history_Elim)
    if save_energy_history:
        fig1[0].savefig(save_dir+'/energy_history.pdf')

# momentum_change-time history 
if plot_momentum_change_history:
    fig2=plotting.momentum_change_history_plot(t, (P-P[0])/P_abs, units,
                                               tlim=momentum_change_history_tlim,
                                               Plim=momentum_change_history_Plim)
    if save_momentum_change_history:
        fig2[0].savefig(save_dir+'/momentum_change_history.pdf')

# grid history 
if plot_grid_history:
    fig3=plotting.grid_history_plot(x_grid, t, dx, DT, E_grid_history, 
                                    units,projection=grid_history_projection,
                                    xlim=grid_history_xlim,tlim=grid_history_tlim)
    if save_grid_history:
        fig3[0].savefig(save_dir+'/grid_history.pdf')
    
# selected_mode history 
if plot_selected_modes_history:
    fig4=plotting.selected_modes_history_plot(k, t, grid_kt, selected_modes, 
                                              plot_theoretical_growth_rate, 
                                              theoretical_growth_rate, units,
                                              part=selected_modes_part,
                                              scale=selected_modes_scale,
                                              tlim=selected_mode_history_tlim,
                                              Alim=selected_mode_history_Alim)
    if save_selected_modes_history:
        fig4[0].savefig(save_dir+'/selected_modes_history.pdf')

# all_modes history 
if plot_all_modes_history:
    fig5=plotting.all_modes_history_plot(k, t, dk, DT, grid_kt, units,
                                         part=all_modes_part,scale=all_modes_scale,
                                         tlim=all_mode_history_tlim,
                                         klim=all_mode_history_klim)
    if save_all_modes_history:
        fig5[0].savefig(save_dir+'/all_modes_history.pdf')
    
# tracker_particle_trajectory_history
if plot_tracker_particle_trajectory:
    fig6=plotting.tracker_particle_trajectory_plot(t, R, V, Tracker_index, units,
                                                   space=tracker_particle_trajectory_space,
                                                   tlim=tracker_particle_trajectory_tlim,
                                                   Rlim=tracker_particle_trajectory_Rlim,
                                                   Vlim=tracker_particle_trajectory_Vlim)
    if save_tracker_particle_trajectory:
        fig6[0].savefig(save_dir+'/tracker_trajectory.pdf')
        
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
    if save_omegak:
        fig7[0].savefig(save_dir+'/dispersion_relation.pdf')

plt.show()
