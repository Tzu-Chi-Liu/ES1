import numpy as np
import matplotlib.pyplot as plt
import sys
import os

import input_output
import build
import PIC_functions

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

NT               = parameters['NT']       # Total simulation time
DT                  = parameters['DT']       # Time step

IW              = parameters['IW']
scheme              = parameters['scheme']
solver              = parameters['solver']

NTracker             = parameters['NTracker']          # Number of tracker particles

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
# space='R-t','V-t','R-V'      # Choose to plot 'R-t', 'V-t', 'R-V' (separate by comma to plot mulitple subplots)
# tlim=None
# Rlim=None
# Vlim=None

# dispersion relation
plot_omegak         = parameters['plot_omegak']       # Plot dispersion relation
plot_theoretical_dispersion_relation = parameters['plot_theoretical_dispersion_relation']

# =============================================================================
# Units for physical quantities
# =============================================================================
units={'t':'(arb. unit)',
       'm':'(arb. unit)',
       'q':'(arb. unit)',
       'r':'(arb. unit)',
       'v':'(arb. unit)',
       'rho_grid':'(arb. unit)',
       'phi_grid':'(arb. unit)',
       'E_grid':'(arb. unit)',
       'k':'(arb. unit)',
       'omega':'(arb. unit)',
       'Momentum':'(arb. unit)',
       'Energy':'(arb. unit)',
       'arb. unit':'(arb. unit)'
       }

if UnitSystem=='SI':
    units['t']='(s)'
    units['m']='(kg)'
    units['q']='(C)'
    units['r']='(m)'
    units['v']='(m/s)'
    units['rho_grid']='(C/m)'
    units['phi_grid']='(V)'
    units['E_grid']='(V/m)'
    units['k']='(1/m)'
    units['omega']='(1/s)'
    units['Momentum']='(kg*m/s)'
    units['Energy']='(J)'
    
if UnitSystem=='Normalized':
    units['t']='(1/ω_p)'
    units['m']='(kg)'
    units['q']='(C)'
    units['r']='(c/ω_p)'
    units['v']='(c)'
    units['rho_grid']='(C*ω_p/c)'
    units['phi_grid']='(V)'
    units['E_grid']='(V*ω_p/c)'
    units['k']='(ω_p/c)'
    units['omega']='(ω_p)'
    units['Momentum']='(kg*c)'
    units['Energy']='(J)'
    
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

print(('Assuming half of N are electrons, other half are protons,\n'+
      'plasma frequency = %.4f '+units['omega']+',\nomega_plasma*DT = %.4f, total steps = %i\n')
      %(omega_plasma,omega_plasma*DT,NT))

# =============================================================================
# initial conditions
# =============================================================================
# if InitialCondition=='Particle Pair Oscillation':
#     N=2
#     m=np.array([m_e,m_e])
#     q=np.array([e,-e])
#     r=np.array([L/4.,2*L/4.])
#     v=np.array([0.,0.])
        
if InitialCondition=='Plasma_Oscillation':
    A          = parameters['A']
    modes  = parameters['Modes']
    v0         = parameters['v0']        # velocity
    v_sigma   = parameters['v_sigma']
    T         = parameters['T']
    n0         = 0.5*parameters['N']/parameters['L']
    
    # first half electrons, second half protons
    m=np.ones(N)*m_e
    m[N//2:]=m_p
    q=-1.*np.ones(N)*e
    q[N//2:]*=-1.
    r=np.append(np.linspace(0,L,N//2),np.linspace(0,L,N//2))
    v=np.append(np.random.normal(v0,v_sigma,size=(1,N//2)),np.zeros(N//2))
    
    
    # excite mode in modes list
    for mode in modes:
        r[:N//2]+=A/n0*np.sin(k[mode]*r[:N//2])
    r=r%L
    
    # theoretical dispersion relation
    if v_sigma==0:
        theoretical_omega_of_k=omega_plasma*np.cos(0.5*dx*k)
    else:
        theoretical_omega_of_k=np.sqrt(omega_plasma**2+3*v_sigma**2*k**2)
    
    # theoretical growth rate
    theoretical_growth_rate=np.zeros_like(k)
    
if InitialCondition=='Two_Stream_Instability':
    v0         = parameters['v0']        # velocity
    v0_sigma   = parameters['v0_sigma']        # velocity width
    charge     = parameters['charge']
    A          = parameters['A']
    modes      = parameters['Modes']

    
    m=np.ones(N)*m_e 
    q=-1*np.ones(N)*e
    if charge == 'opposite':
        q[N//2:]*=-1.
    r=np.random.random(N)*L
    v=np.append(np.random.normal(v0,v0_sigma,size=(1,N//2)),
                np.random.normal(-v0,v0_sigma,size=(1,N//2)))
    
    for mode in modes:
        r[:N//2]+=A*mode/density*np.sin(k[mode]*r[:N//2])
    r=r%L
    
    # theoretical dispersion relation & growth rate
    # theoretical_omeag_of_k=[np.zeros_like(k),np.zeros_like(k)]
    # theoretical_growth_rate=np.zeros_like(k)
    # for mode in range(len(k)):
    #     if omega_plasma*np.sqrt(4.*k[mode]**2*v0**2+omega_plasma**2)>k[mode]**2*v0**2+omega_plasma**2:
    #         theoretical_growth_rate[mode]=np.sqrt(-k**2*v0**2-omega_plasma**2
    #                                               +omega_plasma*np.sqrt(4.*k**2*v0**2+omega_plasma**2))
    #     else:
    #         theoretical_omega_of_k[0][mode]=np.sqrt(k**2*v0**2+omega_plasma**2
    #                                                 +omega_plasma*np.sqrt(4.*k**2*v0**2+omega_plasma**2))
    #         theoretical_omega_of_k[1][mode]=np.sqrt(k**2*v0**2+omega_plasma**2
    #                                                 -omega_plasma*np.sqrt(4.*k**2*v0**2+omega_plasma**2))
    theoretical_omega_of_k=[np.sqrt(k**2*v0**2+omega_plasma**2
                                    +omega_plasma*np.sqrt(4.*k**2*v0**2+omega_plasma**2)),
                            -np.sqrt(k**2*v0**2+omega_plasma**2
                                    +omega_plasma*np.sqrt(4.*k**2*v0**2+omega_plasma**2)),
                            np.sqrt(k**2*v0**2+omega_plasma**2
                                    -omega_plasma*np.sqrt(4.*k**2*v0**2+omega_plasma**2)),
                            -np.sqrt(k**2*v0**2+omega_plasma**2
                                    -omega_plasma*np.sqrt(4.*k**2*v0**2+omega_plasma**2))]
    for omegak in theoretical_omega_of_k:
        omegak[np.isnan(omegak)] = 0.
    
    theoretical_growth_rate=np.sqrt(-k**2*v0**2-omega_plasma**2
                                    +omega_plasma*np.sqrt(4.*k**2*v0**2+omega_plasma**2))
    theoretical_growth_rate[np.isnan(theoretical_growth_rate)]=0.
    
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
    v0         = parameters['v0']        # velocity
    v0_sigma   = parameters['v0_sigma']        # velocity width
    charge     = parameters['charge']
    A          = parameters['A']
    modes  = parameters['Modes']
    
    m=np.ones(N)*m_e 
    q=-1*np.ones(N)*e
    if charge == 'opposite':
        q[N//2:]*=-1.
    r=np.random.random(N)*L
    v=np.append(np.random.normal(v0,v0_sigma,size=(1,N//2)),
                np.random.normal(-v0,v0_sigma,size=(1,N//2)))
    
    for mode in modes:
        r[:N//2]+=A*mode/density*np.sin(k[mode]*r[:N//2])
    r=r%L
    
# =============================================================================
# More derived parameters
# =============================================================================
lambda_debye        = np.sqrt((epsilon0*k_B*T)/((N/(2*L)*e**2))) # Deybe length
N_D                 = N/(2*L)*lambda_debye
print(('Assuming half of N are electrons, other half are protons,\n'+
      'Debye length = %.4f '+units['r']
      +',\nDebye length/L = %.4f, number of particles in one Debye length = %.4f\n')
      %(lambda_debye,lambda_debye/L,N_D))

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
# diagnostics
# =============================================================================
def diagnostics_animation(r,v,phi_grid,E_grid,rho_grid,save_dir):
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
    ax3.plot(x_grid,rho_grid-np.mean(rho_grid))
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
    ax4.plot(x_grid,phi_grid,'blue')
    ax8.plot(x_grid,E_grid,'orange')
    ax4.plot(x_grid,np.zeros(NG),'black',lw=1)
    ax4.set_xlabel('position x '+units['r'])
    ax4.set_ylabel(r'potential $\phi$ '+units['phi_grid'],color='blue')
    ax8.set_ylabel('electric field E '+units['E_grid'],color='orange')
    ax4.set_title('Potential and Electric Fields')
    ax4.set_ylim(-1.1*np.max(np.abs(phi_grid)),1.1*np.max(np.abs(phi_grid)))
    ax8.set_ylim(-1.1*np.max(np.abs(E_grid)),1.1*np.max(np.abs(E_grid)))
    # ax4.set_ylim(-1.5e-8,1.5e-8)
    # ax8.set_ylim(-1e-7,1e-7)
    # ax4.set_title('Potential')
    
    animation_fig.suptitle(InitialCondition+'\n'+'Snapshot at t = %.3f (s)'%t[step])
    plt.tight_layout()
    plt.pause(0.1)
    
    if save_animation:
        animation_save_dir=save_dir+'/animation/'
        if not os.path.isdir(animation_save_dir):
            os.mkdir(animation_save_dir)
        animation_fig.savefig(animation_save_dir+str(step)+'.png')
        
    return animation_fig

# =============================================================================
# dispersion relation
# =============================================================================
def dispersion_relation(grid_history,save_dir):
    global x_grid,t
    # plot grid_history(x,t)
    fig4=plt.figure()
    ax7=fig4.add_subplot(projection='3d')
    x_mesh,T_mesh=np.meshgrid(x_grid,t)
    grid_history=np.transpose(grid_history)
    ax7.plot_surface(T_mesh,x_mesh,grid_history)
    ax7.set_xlabel('time t '+units['t'])
    ax7.set_ylabel('position x '+units['r'])
    ax7.set_zlabel(r'potential $\phi$ '+units['phi_grid'])
    ax7.set_title('Field history')
            
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
        r[i],v[i]=PIC_functions.update_motion(m[i],q[i],r[i],v[i],E_par,scheme,IW,DT,L,E_grid,step)
        
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
        diagnostics_animation(r,0.5*(v_old+v),phi_grid,E_grid,rho_grid,save_dir)
        
    # save output
    if save_output:
        input_output.output_to_file(step,t,N,m,q,r,v,NG,x_grid,rho_grid,phi_grid,E_grid,
                                    InitialCondition,save_dir,units,output_decimal_places)
    
    print(('step = %d '+units['t']
           +', P = %+.3e'+units['Momentum']
           +', E = %+.3e'+units['Energy'])
          %(step,P[step],E_D[step]+E_T[step]+E_F[step]))
    
# =============================================================================
# Analysis    
# =============================================================================
grid_kt=analysis.grid_xt_to_kt(E_grid_history)
grid_omegak=analysis.grid_xt_to_omegak(k,E_grid_history)

# =============================================================================
# Plotting
# =============================================================================
# # dispersion relation
# if plot_omegak:
    # grid_omegak=dispersion_relation(E_grid_history,save_dir)

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
