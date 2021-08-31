import sys

# =============================================================================
# Load input file into parameters dict
# =============================================================================
# program_name=sys.argv[0]
# input_file_loc=sys.argv[1]
input_file_loc='simulation_results/Two_stream_instability/test/input.txt'
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
print(parameters)
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
v0                  = float(parameters['v0'])        # Velocity
v0_sigma            = float(parameters['v0_sigma'])        # Velocity width
A                   = float(parameters['A'])       # Sinusoidal perturbation amplitude
Mode                = int(parameters['Mode'])          # Sinusoidal perturbation wavevector 

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
