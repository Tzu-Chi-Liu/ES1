#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 00:46:21 2021

@author: mac
"""
import datetime

input_file_loc=sys.argv[1]
input_file_loc='simulation_results/Two_stream_instability/test/input.txt'
InitialCondition='Two_stream_instability'

specific_title=input_file_loc.split('/')[2]
time='#'+datetime.datetime.now().strftime('%Y%m%d%H%M')
def save_output(m,q,r,v,phi_grid,rho_grid,E_grid,InitialCondition,specific_title,time):
# Save m,q,r(t),v(t),phi(t),rho(t),E(t)
    if specific_title!='':
        save_dir='simulation_results/'+InitialCondition+'/'+specific_title+'/results'
    else:
        save_dir='simulation_results/'+InitialCondition+'/'+time+'/results'
    
    particle_save_filename=save_dir+'/particle/'+str(t)+'.txt'
    with open(particle_save_filename,'w') as f:
        f.write('# mass\n')
        for mass in m:
            f.write(str(mass)+'\n')
        f.write('# charge\n')
        for charge in q:
            f.write(str(charge)+'\n')
        f.write('# position\n')
        for position in r[t]:
            f.write(str(position)+'\n')
        f.write('# velocity\n')
        for velocity in v[t]:
            f.write(str(velocity)+'\n')
        
    field_save_filename=save_dir+'/field/'+str(t)+'.txt'
    with open(field_save_filename,'w') as f:
        f.write('# rho_grid\m')
        for rho in rho_grid:
            f.write(str(rho)+'\n')
        f.write('# phi_grid\n')
        for phi in phi_grid:
            f.write(str(phi)+'\n')
        for E in E_grid:
            f.write(str(E)+'\n')
    
t=0
for t in range(1000):
    save_output(m,q,r,v,phi_grid,rho_grid,E_grid,InitialCondition,specific_title,time)