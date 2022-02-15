#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 25 02:22:42 2021

@author: mac
"""
from analysis_and_plotting import read_output

file_loc='simulation_results/Numeric_Instability/Single_Cold_Stream/untitled folder 2/results/particle/00010.txt'
r,v=read_output.read_snapshot(file_loc)

import input_output
input_folder_loc = 'simulation_results/Numeric_Instability/Single_Cold_Stream/untitled folder 2/inputs' # for running in IDE (ex:Spyder)
input_txt_parameters, save_dir, species_parameters = input_output.load_input(input_folder_loc)

from analysis_and_plotting import analysis
UnitSystem          = input_txt_parameters['UnitSystem']
units = analysis.generate_units(UnitSystem)

from analysis_and_plotting import plotting
fig, ax = plotting.phase_space_plot(species_parameters, r, v, units)

file_loc = 'simulation_results/Numeric_Instability/Single_Cold_Stream/untitled folder 2/results/field/00010.txt'
x_grid,rho_grid,phi_grid,E_grid=read_output.read_snapshot(file_loc)
fig, ax = plotting.grid_plot(x_grid, E_grid, units, r'Potential $\phi$ '+units['phi_grid'], title='Potential')
