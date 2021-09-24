#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 00:20:18 2021

@author: mac
"""
import os
import numpy as np

def read_snapshot(file_loc):
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
        
        return m,q,r,v
        
    elif 'field' in file_loc:
        x_grid   = data[:,1]
        rho_grid = data[:,2]
        phi_grid = data[:,3]
        E_grid   = data[:,4]
        
        return x_grid,rho_grid,phi_grid,E_grid
    
def read_history(folder_loc):
    filenames=[file for file in os.listdir(folder_loc) if os.path.isfile(os.path.join(folder_loc,file))]
    filenames.sort()
    
    t=[]
    if 'particle' in folder_loc:
        R=[]
        V=[]
        for filename in filenames:
            with open(os.path.join(folder_loc,filename)) as f:
                first_line=f.readline()
                t.append(float(first_line.split(' = ')[1].split(' (')[0]))
        return np.array(t),np.array(R),np.array(V)
    
    elif 'field' in folder_loc:
        x_grid=[]
        grid_history=[]
        for filename in filenames:
            with open(os.path.join(folder_loc,filename)) as f:
                first_line=f.readline()
                t.append(float(first_line.split(' = ')[1].split(' (')[0]))
            data=np.loadtxt(os.path.join(folder_loc,filename))
            
            x_grid=data[:,1]
            grid_history.append(data[:,4])  # E_grid_history
        return np.array(t),np.array(x_grid),np.array(grid_history)

# =============================================================================
# 
# =============================================================================
if __name__=='__main__':
    file_loc='/Users/mac/小皮球/school/NTU/physics department/computational physics/PIC/ES1/simulation_results/Two_Stream_Instability/Two_Stream_Instability/results/particle/00000.txt'
    m,q,r,v=read_snapshot(file_loc)
    # x_grid,rho_grid,phi_grid,E_grid=read_snapshot(file_loc)
    
    folder_loc='../simulation_results/cold_plasma_oscillation/EXAMPLE/results/field/'
    # R,V=read_history(folder_loc)
    t,x_grid,grid_history=read_history(folder_loc)
