import sys
import datetime
import numpy as np

# =============================================================================
# Load input file into parameters dict & save to same dir as input file
# =============================================================================
# input_file_loc=sys.argv[1] # For running in terminal
input_file_loc='simulation_results/example_input_file/input.txt' # for running in IDE (ex:Spyder)

def load_input(input_file_loc):
    '''
    For reading parameters from `input.txt` and generating save_dir for simulation outputs

    Parameters
    ----------
    input_file_loc : str
        Path to `input.txt`.

    Returns
    -------
    parameters : dict
        Parameters for simulation .
    save_dir : str
        Path to directory which saves simulation results.

    '''
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
    
    return parameters, save_dir

parameters, save_dir = load_input(input_file_loc)
# =============================================================================
# physical constants
# =============================================================================
e_e                 = float(parameters['e_e'])     # elementary charge
epsilon0            = float(parameters['epsilon0'])        # vacuum permittivity
m_e                 = float(parameters['m_e'])        # electron mass
m_p                 = float(parameters['m_p'])     # proton mass

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

T_end               = float(parameters['T_end'])       # Total simulation time
dt                  = float(parameters['dt'])       # Time step

weight              = parameters['weight']
scheme              = parameters['scheme']
solver              = parameters['solver']

# =============================================================================
# Analysis parameters
# =============================================================================
Tracker             = int(parameters['Tracker'])          # Number of tracker particles
plot_animation      = bool(parameters['plot_animation'])       # Plot animation of snapshots
save_animation      = bool(parameters['save_animation'])       # save animation of snapshots
plot_history        = bool(parameters['plot_history'])       # Plot history
save_history        = bool(parameters['save_history'])       # save history
plot_omegak         = bool(parameters['plot_omegak'])       # Plot dispersion relation
save_omegak         = bool(parameters['save_omegak'])       # save dispersion relation
plot_trajectory     = bool(parameters['plot_trajectory'])      # Plot trajectory of tracker particle(s)

# =============================================================================
# Derived parameters 
# =============================================================================
dx                  = L/NG # grid spacing
x                   = np.arange(0,L,dx) # np array for x axis
k                   = 2.*np.pi/L  # wave number for first mode 
T                   = np.arange(0,T_end,dt) # np array for time axis
# density             = N/L      # not sure to input or calculate
omega_plasma        = np.sqrt(((N/(2*L))*e_e**2)/(epsilon0*m_e)) # plasma frequency

print('Assuming half of N are electrons, other half are protons,\n'+
      'plasma frequency = %.4f (s^{-1}),\nomega_plasma*dt = %.4f, total steps = %i\n'
      %(omega_plasma,omega_plasma*dt,int(T_end/dt)))