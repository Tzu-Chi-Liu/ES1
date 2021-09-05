# ES1A 1-dimensional electrostatic particle-in-cell plasma simulation code## Requirements- Python - Numpy- Matplotlib## Running the code- Terminal:     1. Uncomment `input_file_loc=sys.argv[1] # For running in terminal` in ES1.py (should be at line 10)    2. Comment `input_file_loc='simulation_results/example_input_file/input.txt' # for running in IDE (ex:Spyder)` in ES1.py (should be at line 11)    3. Run `python ES1.py <path_to_input.txt>` in terminal- IDE(ex: Spyder):     1. Uncomment `input_file_loc='simulation_results/example_input_file/input.txt' # for running in IDE (ex:Spyder)` in ES1.py (should be at line 11)    2. Change the path in `''` to <path_to_input.txt>    3. Comment `input_file_loc=sys.argv[1] # For running in terminal` in ES1.py (should be at line 10)    4. Run normally as any other scripts## Simulation problems- Cold plasma oscillations- Warm plasma waves- Two stream instability- Beam plasma instability- Landau damping## To do list - Strange oscillation for charge density in plasma oscillation- Leave fewer plotting functions and delete the more similar ones? - Output file should record variables at 0*dt, dt, 2*dt, 3*dt, ..., NT*dt (NT+1 files)- Warning for overwriting existing simulation results- Record used variable names- fft draw real() part, imag() part, abs() part- imshow() check axis extent is correct or not- imshow() check if index=0 matches axis=0- Choose between T and v_sigma- For particles, choose in input.txt to save in output files (1 file per step, m and q use only 1 file per simulation):    1. x-v space distribution function    2. x-space and v-space distribution function    3. r, v for all particles- output file: single row for m, q, r, v of same particle/ rho_grid, phi_grid, E_grid of same position - Negative k?- Add n_grid, P_grid, T_grid- Use separate `example_input_file/input.txt` for different Simulation problems- Add phi(k,t)- Add 2d imshow for field_history- Document which axis corresponds to what for field_history, field_kt, field_omegak- Set amplitude for excited modes in `input.txt`- Normalized units, cgs(Heaviside) (Need to rewrite equation?)?- Electrostatic energy for each mode- Will naming all figures axes to `fig,ax` in functions make bad things happen?- Add and plot theoretical dispersion relations for Simulation problems listed below- Set axis range for plots - Add colorbar for imshow() plots- Save particle distribution function at grid positions (1 distribution function for 1 speices)?- Separate plot_* from save_*- Input follow Birdsall ES1 (use NSP files for each species?)- One function for plotting one figure?- Use steps instead of DT?- For `input.txt` and `build.py`:    - Write prepared-for-running parameters for Simulation problems listed below- Output file frequency- Different unit systems (SI&normalized)- Use consistent variable names (capitalize, underscore, etc)- Need to explicitly add fixed ion background?- Input file use .py instead of .txt?- Rename parameters with self-explanatory names (ex:t->step)- on/off for plt.show()- Read from input file - Save output to file- Read from output file for further analysis (ex: plotting)- Add plotters- Use classes?- 1d2v ES1 with (const) B field- Benchmark program- Add serial/parallel c++ version## Stuff- Need to save `input.txt` before running?- DT affects energy conservation- FFT:    - NG and L (and dx = L/NG) affects highest possible k and k space resolution        - Mode 1 k1 = 2.*np.pi/L = dk        - Highest mode kn = 2.*np.pi/(2.*dx) = np.pi/dx = np.pi*NG/L        - Number of modes = NG        - Increase resolution in k space <=> decrease dk <=> increase L        - Similarily, for DT and NT (and T_end=DT*NT)         - Mode 1 omega1 = 2.*np.pi/(DT*NT) = domega        - Highest mode omegan = 2.*np.pi/(2.*DT) = np.pi/DT = np.pi*steps/T_end        - Number of modes = steps        - Increase resolution in omega space <=> decrease domega <=> increase T_end <=> increase NT- Don't modify .py and .cpp for input and simulation parameters- Memory requirement of simulation results:    - Assuming 1D code (ES1) with E field in 128 grid, 40000 particles        - Particle parameters: m, q, r, v        - Field parameters: rho, phi, E        - Data size for each step = (40000*4+128x3)*8(bytes/parameter)=1283072 bytes ~ 1.28 MB/step (at least)        - Saving x-v space distribution functions:            - Distribution function shape = (NG, Nv) where Nv is the number of grid points in v space            - If NG*Nv < N*4: Save distribution function, else: Save m, q, r, v        - Only save v-space distribution function and x-space distribution function instead?    - Assuming fully 3D code with both E and B fields in 128x128x128 grid, 10000 particles        - Particle parameters: m, q, r_x, r_y, r_z, v_x, v_y, v_z        - Field parameters: rho, J_x, J_y, J_z, E_x, E_y, E_z, B_x, B_y, B_z        - Data size for each step = (10000*8+128**3*10)*8(bytes/parameter)=168412160 bytes ~ 168 MB/step (at least)