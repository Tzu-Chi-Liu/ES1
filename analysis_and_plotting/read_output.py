#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 00:20:18 2021

@author: mac
"""
import numpy as np

file_loc='../simulation_results/cold_plasma_oscillation/EXAMPLE/results/field/00000.txt'
with open(file_loc) as f:
    header=[]
    for lines in f:
        if lines.startswith('#'):
            header.append(lines.replace('\n',''))
data=np.loadtxt(file_loc)

if 'particle' in file_loc:
    m = data[:,1]
    q = data[:,2]
    r = data[:,3]
    v = data[:,4]
    
if 'field' in file_loc:
    x        = data[:,1]
    rho_grid = data[:,2]
    phi_grid = data[:,3]
    E_grid   = data[:,4]