#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 00:20:18 2021

@author: mac
"""
file_loc=''
def read_from_output(file_loc):
    with open(file_loc) as f:
        file=f.read().splitlines()
    return file
file=read_from_output(file_loc)