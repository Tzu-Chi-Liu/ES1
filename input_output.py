#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  4 01:21:07 2021

@author: mac
"""
import sys
import os
import datetime

# =============================================================================
# Load input file into parameters dict & save to same dir as input file
# =============================================================================
def load_input(input_folder_loc):
    '''
    For reading parameters from `input.txt` and generating save_dir for simulation outputs

    Parameters
    ----------
    input_folder_loc : str
        Path to foler containing `input.txt`, `species_1.txt`, `species_2.txt`, ...

    Returns
    -------
    input_txt_parameters : dict
        Parameters for simulation stored in `input.txt`.
    save_dir : str
        Path to directory which saves simulation results.
    species_parameters: list
        List of dicts of length NSP, where each dict contains parameters for each particle species

    '''
    print('Input folder location:\n',input_folder_loc,'\n')
    
    input_files=[filename for filename in os.listdir(input_folder_loc) 
                 if os.path.isfile(os.path.join(input_folder_loc,filename))]
    input_files.sort()
    
    species_parameters=[]
    for filename in input_files:
        input_file_loc=os.path.join(input_folder_loc,filename)
        
        with open(input_file_loc) as f:
            text=f.read().splitlines()
        text=[line for line in text if not line.startswith('#')]
        text=[line for line in text if not line=='']
        text=[[line.replace(' ','').split('=')[0],line.replace(' ','').split('=')[1].split('#')[0]] 
              for line in text]
        
        parameters={}
        for line in text:
            parameters[line[0]]=line[1]
        
        # Convert variable types
        for keys in parameters.keys():
            if ('plot' in keys) or ('save' in keys) \
                or ('True' in parameters[keys]) or ('False' in parameters[keys]):
                if parameters[keys]=='True' or parameters[keys]=='1':
                    parameters[keys]=True
                elif parameters[keys]=='False'or parameters[keys]=='0':
                    parameters[keys]=False
            elif 'lim' in keys:
                if ('(' in parameters[keys] and ')' in parameters[keys]):
                    parameters[keys]=tuple([float(lim) for lim in parameters[keys].strip('(').strip(')').split(',')])
                elif parameters[keys]=="'default'":
                    parameters[keys]=parameters[keys].strip("'")
            elif "'" in parameters[keys]:
                parameters[keys]=parameters[keys].replace("'",'')
            elif 'N' in keys:
                parameters[keys]=int(parameters[keys])
            elif ('.' in parameters[keys]) or ('e' in parameters[keys]):
                parameters[keys]=float(parameters[keys])
            elif ('[' in parameters[keys] and ']' in parameters[keys]):
                if parameters[keys]=='[]':
                    parameters[keys]=[]
                else:
                    parameters[keys]=[int(numbers) for numbers in parameters[keys].strip('[').strip(']').split(',')]
            
        if 'input' in filename:
            input_txt_parameters=parameters
            # define save_dir
            save_dir=input_file_loc.split('/')[0]+'/'+input_file_loc.split('/')[1]\
                 +'/'+input_file_loc.split('/')[2]+'/results'
            # time='#'+datetime.datetime.now().strftime('%Y%m%d%H%M')
            # save_dir=input_file_loc.split('/')[0]+'/'+input_file_loc.split('/')[1]\
            #          +'/'+time+'/results'
            
            # make save_dir folder
            save_output         = parameters['save_output']
            if save_output:
                if (not os.path.isdir(save_dir)):
                    os.mkdir(save_dir)
                    
        elif 'species' in filename:
            if parameters['m']=='m_e':
                parameters['m']=input_txt_parameters['m_e']
            elif parameters['m']=='m_p':
                parameters['m']=input_txt_parameters['m_p']

            if parameters['q']=='e':
                parameters['q']=input_txt_parameters['e']
            elif parameters['q']=='-e':
                parameters['q']=-input_txt_parameters['e']
                
            species_parameters.append(parameters)
    
    return input_txt_parameters, save_dir, species_parameters

# =============================================================================
# Output
# =============================================================================
def output_to_file(step,t,N,m,q,r,v,NG,x_grid,rho_grid,phi_grid,E_grid,
                InitialCondition,save_dir,units,output_decimal_places):
        
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
                 +'x '+'{:'+str(output_decimal_places+7)+'}'
                 +'ρ '+'{:'+str(output_decimal_places+7)+'}'
                 +'φ '+'{:'+str(output_decimal_places+7)+'}'
                 +'E '+'{:'+str(output_decimal_places+7)+'}'+'\n')
                .format(units['r'],units['rho_grid'],units['phi_grid'],units['E_grid']))
        for grid in range(NG):
            f.write(('%6d  '
                     +'%+.'+str(output_decimal_places)+'e  '
                     +'%+.'+str(output_decimal_places)+'e  '
                     +'%+.'+str(output_decimal_places)+'e  '
                     +'%+.'+str(output_decimal_places)+'e'+'\n')
                    %(grid,x_grid[grid],rho_grid[grid],phi_grid[grid],E_grid[grid]))
            
# =============================================================================
# 
# =============================================================================
if __name__=='__main__':
    # input_folder_loc=sys.argv[1] # For running in terminal
    input_folder_loc='simulation_results/EXAMPLE/inputs' # for running in IDE (ex:Spyder)
    input_txt_parameters, save_dir, species_parameters = load_input(input_folder_loc)
    print('save_dir: ')
    print(save_dir)
    
    print('\ninput_txt_parameters: ')
    for key in input_txt_parameters.keys():
        print(type(input_txt_parameters[key]),key+':',input_txt_parameters[key])
    
    print('\nNumber of species NSP =',len(species_parameters))
    for species in range(len(species_parameters)):
        print('\nSpecies', species+1)
        for key in species_parameters[species].keys():
            print(type(species_parameters[species][key]),key+':',
                  species_parameters[species][key])
        
            