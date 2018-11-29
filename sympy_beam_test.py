#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 13:12:58 2018

@author: ethan
"""
from sympy.physics.continuum_mechanics.beam import Beam as beam_mech
from sympy import symbols
E,L,I = symbols('W,L,I')
R1, M1 = symbols('R1, M1')
L = 9.0
beam_sym = beam_mech(L,E,I)
#point load up
beam_sym.apply_load(R1,0,-1)
beam_sym.apply_load(M1,0,-2)
beam_sym.apply_load(-8,0,0, end = 5)
beam_sym.apply_load(50,5,-2,)
beam_sym.apply_load(-12,9,-1)
beam_sym.bc_slope.append((0,0))
beam_sym.bc_deflection.append((0,0))
beam_sym.solve_for_reaction_loads(R1,M1)
print(beam_sym.reaction_loads)
