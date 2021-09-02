#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 00:46:21 2021

@author: mac
"""
import datetime

N=5
NG=11
t=0.01

def save_output(m,q,r,v,rho_grid,phi_grid,E_grid,
                InitialCondition,specific_title,units,output_decimal_places,time):
# Save m,q,r(t),v(t),phi(t),rho(t),E(t)
    if specific_title!='':
        save_dir='simulation_results/'+InitialCondition+'/'+specific_title+'/results'
    else:
        save_dir='simulation_results/'+InitialCondition+'/'+time+'/results'
    
    particle_save_filename=save_dir+'/particle/'+str(t)+'.txt'
    with open(particle_save_filename,'w') as f:
        f.write(('# ES1 particle output at t = %.'+str(output_decimal_places)+'f ')%t+units['t']+'\n')
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
                   +'%+.'+str(output_decimal_places)+'e  '+'\n')
                  %(particle,m[particle],q[particle],r[particle],v[particle]))
        
    field_save_filename=save_dir+'/field/'+str(t)+'.txt'
    with open(field_save_filename,'w') as f:
        f.write(('# ES1 field output at t = %.'+str(output_decimal_places)+'f ')%t+units['t']+'\n')
        f.write(('# grid  '
               +'ρ '+'{:'+str(output_decimal_places+7)+'}'
               +'φ '+'{:'+str(output_decimal_places+7)+'}'
               +'E '+'{:'+str(output_decimal_places+7)+'}'+'\n')
              .format(units['rho_grid'],units['phi_grid'],units['E_grid']))
        for grid in range(NG):
            f.write(('%6d  '
                   +'%+.'+str(output_decimal_places)+'e  '
                   +'%+.'+str(output_decimal_places)+'e  '
                   +'%+.'+str(output_decimal_places)+'e  '+'\n')
                  %(grid,rho_grid[grid],phi_grid[grid],E_grid[grid]))
    