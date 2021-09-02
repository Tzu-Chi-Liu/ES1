import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import datetime

import build
import PIC_functions
from analysis_and_plotting import plotting

# =============================================================================
# Load input file into parameters dict & save to same dir as input file
# =============================================================================
# input_file_loc=sys.argv[1] # For running in terminal
input_file_loc='simulation_results/cold_plasma_oscillation/EXAMPLE/input.txt' # for running in IDE (ex:Spyder)
specific_title=input_file_loc.split('/')[2]
time='#'+datetime.datetime.now().strftime('%Y%m%d%H%M')
print('Input file location:',input_file_loc,'\n')

with open(input_file_loc) as f:
    text=f.read().splitlines()
text=[line for line in text if not line.startswith('#')]
text=[line for line in text if not line=='']
text=[[line.replace(' ','').split('=')[0],line.replace(' ','').split('=')[1].split('#')[0]] 
      for line in text]

parameters={}
for line in text:
    parameters[line[0]]=line[1]

save_dir=input_file_loc.split('/')[0]+'/'+input_file_loc.split('/')[1]\
         +'/'+input_file_loc.split('/')[2]+'/results'
# save_dir=input_file_loc.split('/')[0]+'/'+input_file_loc.split('/')[1]\
#          +'/'+time+'/results'
if (not os.path.isdir(save_dir)):
    os.mkdir(save_dir)

# =============================================================================
# Unit system
# =============================================================================
UnitSystem         = parameters['UnitSystem']

# =============================================================================
# physical constants
# =============================================================================
e                 = float(parameters['e'])     # elementary charge
epsilon0            = float(parameters['epsilon0'])        # vacuum permittivity
m_e                 = float(parameters['m_e'])        # electron mass
m_p                 = float(parameters['m_p'])     # proton mass
k_B                 = float(parameters['k_B']) 

# =============================================================================
# Boundary & initial (physical) conditions
# =============================================================================
N                   = int(parameters['N'])      # number of particles
L                   = float(parameters['L'])        # length of box

density             = float(parameters['density'])       # not sure to input or calculate

InitialCondition    = parameters['InitialCondition']

# =============================================================================
# Simulation (unphysical) parameters
# =============================================================================
NG                  = int(parameters['NG'])         # number of grids

NT               = int(parameters['NT'])       # Total simulation time
DT                  = float(parameters['DT'])       # Time step

IW              = parameters['IW']
scheme              = parameters['scheme']
solver              = parameters['solver']

# =============================================================================
# Analysis parameters
# =============================================================================
Tracker             = int(parameters['Tracker'])          # Number of tracker particles
save_output         = bool(int(parameters['save_output']))
output_decimal_places = int(parameters['output_decimal_places'])
plot_animation      = bool(int(parameters['plot_animation']))  # Plot animation of snapshots
save_animation      = bool(int(parameters['save_animation']))      # save animation of snapshots
plot_history        = bool(int(parameters['plot_history']))       # Plot history
save_history        = bool(int(parameters['save_history']))       # save history
selected_modes = [int(mode) for mode in parameters['selected_modes'].split(',')]
plot_theoretical_growth_rate = bool(int(parameters['plot_theoretical_growth_rate']))
plot_omegak         = bool(int(parameters['plot_omegak']))       # Plot dispersion relation
save_omegak         = bool(int(parameters['save_omegak']))       # save dispersion relation
plot_theoretical_dispersion_relation = bool(int(parameters['plot_theoretical_dispersion_relation']))
plot_trajectory     = bool(int(parameters['plot_trajectory']))  # Trajectory of tracker particle(s)

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
x                   = np.arange(0.,L,dx) # np array for x axis
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
    A          = float(parameters['A'])
    if ':' in parameters['Modes']:
        modes  = list(range(int(parameters['Modes'].split(':')[0]),
                            int(parameters['Modes'].split(':')[1]),
                            int(parameters['Modes'].split(':')[2])))
    else:
        modes  = [int(mode) for mode in parameters['Modes'].split(',')]
    v_sigma   = float(parameters['v_sigma'])
    T         = float(parameters['T'])
    v_gt          = float(parameters['v_gt'])         # velocity
    n0         = 0.5*float(parameters['N'])/float(parameters['L'])
    
    # first half electrons, second half protons
    m=np.ones(N)*m_e
    m[N//2:]=m_p
    q=-1.*np.ones(N)*e
    q[N//2:]*=-1.
    r=np.append(np.linspace(0,L,N//2),np.linspace(0,L,N//2))
    v=np.append(np.random.normal(0,v_sigma,size=(1,N//2)),np.zeros(N//2))
    
    # Galilean transformation
    v+=v_gt
    
    # excite mode in modes list
    for mode in modes:
        r[:N//2]+=A*mode/n0*np.sin(k[mode]*r[:N//2])
    r=r%L
    
    # theoretical dispersion relation
    if v_sigma==0:
        theoretical_omega_of_k=omega_plasma*np.cos(0.5*dx*k)
    else:
        theoretical_omega_of_k=np.sqrt(omega_plasma**2+3*v_sigma**2*k**2)
        
if InitialCondition=='Two_Stream_Instability':
    v0         = float(parameters['v0'])        # velocity
    v0_sigma   = float(parameters['v0_sigma'])        # velocity width
    charge     = parameters['charge']
    A          = float(parameters['A'])
    if ':' in parameters['Modes']:
        modes  = list(range(int(parameters['Modes'].split(':')[0]),
                            int(parameters['Modes'].split(':')[1]),
                            int(parameters['Modes'].split(':')[2])))
    else:
        modes  = [int(mode) for mode in parameters['Modes'].split(',')]
    
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
    v0         = float(parameters['v0'])        # velocity
    v0_sigma   = float(parameters['v0_sigma'])        # velocity width
    charge     = parameters['charge']
    A          = float(parameters['A'])
    if ':' in parameters['Modes']:
        modes  = list(range(int(parameters['Modes'].split(':')[0]),
                            int(parameters['Modes'].split(':')[1]),
                            int(parameters['Modes'].split(':')[2])))
    else:
        modes  = [int(mode) for mode in parameters['Modes'].split(',')]
    
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
R     = np.zeros((Tracker,len(t))) # tracker particle phase space trajectory
V     = np.zeros((Tracker,len(t)))

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
def external_field(x,t):
    phiX=np.zeros(NG)
    EX=np.zeros(NG)
    
    return phiX,EX

# =============================================================================
# source weighting
# =============================================================================
def weightrho(r,q,IW):
    # Weighting particle positions r to grid density rho_grid
    rho_grid=np.zeros(NG)
    for i in range(N):
        x=int(r[i]//dx)
        if IW=='NGP' or IW=='1':
            if np.abs(r[i]%dx)<0.5*dx:
                rho_grid[x%NG]+=q[i]/dx
            else:
                rho_grid[(x+1)%NG]+=q[i]/dx
        if IW=='CIC' or IW=='2' or IW=='energy conserving' or IW=='3':
            rho_grid[x%NG]+=q[i]/dx*(dx-r[i]%dx)/dx
            rho_grid[(x+1)%NG]+=q[i]/dx*(r[i]%dx)/dx
            
    return rho_grid

# =============================================================================
# Poisson solver
# =============================================================================
def solve_poisson(r,q,IW,solver='FFT',operator='local'):
    rho_grid=weightrho(r,q,IW)
    
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
def weightField(m,q,r,v,E_grid,IW):
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
def update_motion(m,q,r,v,E_grid,scheme,IW):
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
        
    return r,v

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
# history
# =============================================================================
def history(t,E_D,E_T,E_F,P,save_dir):
    fig2,ax5=plt.subplots()
    ax5.plot(t,E_D,label='drift')
    ax5.plot(t,E_T,label='thermal')
    ax5.plot(t,E_F,label='field')
    ax5.plot(t,E_D+E_T+E_F,label='total')
    ax5.set_xlabel('time t '+units['t'])
    ax5.set_ylabel('energy E '+units['Energy'])
    ax5.legend(bbox_to_anchor=(0,1.01,1,0.2),loc='lower left'
               ,mode='expand',borderaxespad=0.,ncol=4)
    ax5.set_title('Energy History',y=1.07)
    
    fig3,ax6=plt.subplots()
    ax6.plot(t,P-P[0])
    ax6.set_xlabel('time t '+units['t'])
    ax6.set_ylabel('total momentum P '+units['Momentum'])
    ax6.set_title('Total Momentum Change History')
    
    if save_history:
        fig2.savefig(save_dir+'/energy_history.png')
        fig3.savefig(save_dir+'/momentum_history.png')
        
    return True

# =============================================================================
# Time evolution of different modes
# =============================================================================
def mode_evolution(grid_history,save_dir):
     # field(k,t)=grid_kt[int(k//k[0]%NG,int(t//DT)%NT]
    grid_kt=np.fft.rfft(grid_history,axis=0)
    
    # plot amplitude v.s. t for grid_kt selected modes
    fig6=plt.figure(figsize=(8,8))
    for index in range(len(selected_modes)):
        ax=fig6.add_subplot(len(selected_modes),1,index+1)
        title='Mode '+str(selected_modes[index])
        ax.plot(t[:],np.abs(grid_kt[selected_modes[index],:]),label='Simulation result')
        if plot_theoretical_growth_rate:
            growth_rate_label='Theoretical growth rate = %.3f'\
                               %theoretical_growth_rate[selected_modes[index]]\
                               +' '+units['omega']
            ax.plot(t[:],np.abs(grid_kt[selected_modes[index],0])
                    *np.exp(theoretical_growth_rate[selected_modes[index]]*t[:]),
                    label=growth_rate_label)
            ax.legend(loc='upper left')
        ax.set_xlabel('Time t '+units['t'])
        ax.set_ylabel('Amplitude '+units['arb. unit'])
        # ax.set_yscale('log')
        ax.set_title(title)
        
    plt.tight_layout()
    
    # imshow grid_kt
    fig7,ax11=plt.subplots()
    pos=ax11.imshow(np.abs(grid_kt)[:20,:300],
                # extent=(t[0]-0.5*DT,t[100]-0.5*DT,k[10]-0.5*dk,k[0]-0.5*dk)
                )
    # ax11.set_xlabel(r'Time $t$ '+units['t'])
    # ax11.set_ylabel(r'Wave number $k$ '+units['k'])
    cbar=fig7.colorbar(pos,ax=ax11)
    ax11.set_xlabel('Step #')
    ax11.set_ylabel('Mode number')
    ax11.set_title('')
    cbar.set_label('Amplitude '+units['arb. unit'],rotation=90)
    
    return grid_kt


# =============================================================================
# dispersion relation
# =============================================================================
def dispersion_relation(grid_history,save_dir):
    global x,t
    # plot grid_history(x,t)
    fig4=plt.figure()
    ax7=fig4.add_subplot(projection='3d')
    x_mesh,T_mesh=np.meshgrid(x,t)
    grid_history=np.transpose(grid_history)
    ax7.plot_surface(T_mesh,x_mesh,grid_history)
    ax7.set_xlabel('time t '+units['t'])
    ax7.set_ylabel('position x '+units['r'])
    ax7.set_zlabel(r'potential $\phi$ '+units['phi_grid'])
    ax7.set_title('Field history')
    
    # imshow grid_history(x,t)
    fig8,ax12=plt.subplots()
    pos=ax12.imshow(grid_history.T)
    ax12.set_title('Field history')
    cbar=fig8.colorbar(pos,ax=ax12)
    cbar.set_label('electric field E '+units['E_grid'],rotation=90)
    
    # FFT grid_history(x,t) to obtain grid_omegak(k,omega)
    grid_omegak=np.fft.rfft2(grid_history)
    lenk=len(k)
    lenomega=len(omega)
    
    # 
    grid_omegak=grid_omegak[:grid_omegak.shape[0]//2,:]
    grid_omegak=grid_omegak[::-1,:]
    
    # plot grid_omegak(k,omega)
    fig5,ax9=plt.subplots()
    pos=ax9.imshow(np.abs(grid_omegak)[9*grid_omegak.shape[0]//10:,:]
               ,extent=(k[0]-0.5*dk,k[-1]-0.5*dk
                        ,omega[0]-0.5*domega,omega[lenomega//10]-0.5*domega))
    cbar=fig5.colorbar(pos,ax=ax9)
    ax9.set_xlabel(r'wave number $k$ '+units['k'])
    ax9.set_ylabel(r'angular frequency $\omega$ '+units['omega'])
    ax9.set_title(InitialCondition.replace('_',' ')+' Dispersion Relation')
    cbar.set_label('Amplitude '+units['arb. unit'],rotation=90)
    
    # plot theoretical dispersion relation    
    if plot_theoretical_dispersion_relation:
        ax9.plot(k,theoretical_omega_of_k,'red',lw=0.5,label='theoretical')
        ax9.legend(loc='upper right')
    
    # plt.tight_layout()
    plt.show()
    
    if save_omegak:
        fig4.savefig(save_dir+'/field_history.png')
        fig5.savefig(save_dir+'/dispersion_relation.png')
            
    return grid_omegak

# =============================================================================
# save output to file
# =============================================================================
def output_to_file(m,q,r,v,rho_grid,phi_grid,E_grid,
                InitialCondition,save_dir,units,output_decimal_places,time):
    # Save m,q,r(t),v(t),phi(t),rho(t),E(t)
    particle_save_dir=save_dir+'/particle/'
    if not os.path.isdir(particle_save_dir):
            os.mkdir(particle_save_dir)
    particle_save_filename=particle_save_dir+f'{step:05}'+'.txt'
    with open(particle_save_filename,'w') as f:
        f.write(('# ES1 particle output at t = %.'+str(output_decimal_places)+'f ')%t[step]+units['t']
                +'\n')
        f.write(('# number  '
               +'m '+'{:'+str(output_decimal_places+7)+'}'
               +'q '+'{:'+str(output_decimal_places+7)+'}'
               +'r '+'{:'+str(output_decimal_places+7)+'}'
               +'v '+'{:'+str(output_decimal_places+7)+'}'+'\n')
              .format(units['m'],units['q'],units['r'],units['v']))
        for particle in range(N):
            f.write(('%8d  '
                   +'%+.'+str(output_decimal_places)+'e  '
                   +'%+.'+str(output_decimal_places)+'e  '
                   +'%+.'+str(output_decimal_places)+'e  '
                   +'%+.'+str(output_decimal_places)+'e'+'\n')
                  %(particle,m[particle],q[particle],r[particle],v[particle]))
        
    field_save_dir=save_dir+'/field/'
    if not os.path.isdir(field_save_dir):
            os.mkdir(field_save_dir)
    field_save_filename=field_save_dir+f'{step:05}'+'.txt'
    with open(field_save_filename,'w') as f:
        f.write(('# ES1 field output at t = %.'+str(output_decimal_places)+'f ')%t[step]+units['t']
                +'\n')
        f.write(('# grid  '
               +'ρ '+'{:'+str(output_decimal_places+7)+'}'
               +'φ '+'{:'+str(output_decimal_places+7)+'}'
               +'E '+'{:'+str(output_decimal_places+7)+'}'+'\n')
              .format(units['rho_grid'],units['phi_grid'],units['E_grid']))
        for grid in range(NG):
            f.write(('%6d  '
                   +'%+.'+str(output_decimal_places)+'e  '
                   +'%+.'+str(output_decimal_places)+'e  '
                   +'%+.'+str(output_decimal_places)+'e'+'\n')
                  %(grid,rho_grid[grid],phi_grid[grid],E_grid[grid]))
            
# =============================================================================
# main simulation loop
# =============================================================================
for step in range(len(t)):
    #solve field
    phi_grid,rho_grid,E_grid=solve_poisson(r,q,IW)[0:3]
    k,rho_k,phi_k=solve_poisson(r,q,IW)[3:6]
    
    v_old=np.copy(v)
    if plot_animation:
        diagnostics_animation(r,0.5*(v_old+v),phi_grid,E_grid,rho_grid,save_dir)
    
    # update particle i velocity v[i] and position r[i] at time t
    for i in range(N):
        r[i],v[i]=update_motion(m[i],q[i],r[i],v[i],E_grid,scheme,IW)
        # tracker particle
        if i==N//4:
            print('r = %.6f(m), v = %.6f(m/s), i = %i, t = %.3f(s)'%(r[i],v[i],step,t[step]))
            R[0,step]=r[i]

    # history
    E_F[step]=0.5*np.dot(phi_grid*dx,rho_grid)
    P[step]=np.sum(m*0.5*(v_old+v))
    E_D[step]=np.sum(0.5*m*v*v_old)
    rho_grid_history[:,step]=rho_grid
    phi_grid_history[:,step]=phi_grid
    E_grid_history[:,step]=E_grid
        
    # save output
    if save_output:
        output_to_file(m,q,r,v,rho_grid,phi_grid,E_grid,
                InitialCondition,save_dir,units,output_decimal_places,time)
    
# history
if plot_history:
    history(t,E_D,E_T,E_F,P,save_dir)

# mode evolution
if len(selected_modes)>0:
    mode_evolution(E_grid_history,save_dir)
    
# dispersion relation
if plot_omegak:
    grid_omegak=dispersion_relation(E_grid_history,save_dir)
