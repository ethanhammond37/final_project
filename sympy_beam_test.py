#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 13:12:58 2018

@author: ethan
"""
import sys
from PyQt5.QtWidgets import (QMainWindow, QDesktopWidget, QApplication, \
                             QWidget, QDialog, QLineEdit, QAction,\
                             QMessageBox,QFileDialog,QMenu, QSizePolicy, \
                             QComboBox, QPushButton, QGridLayout, QLabel)
from PyQt5.QtCore import Qt

from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.figure import Figure



class testtest():
    
    def __init__(self):
        from sympy.physics.continuum_mechanics.beam import Beam as beam_mech
        from sympy import symbols
        E,L,I = symbols('E,L,I')
        self.R1, self.R2, self.M1 = symbols('R1,R2, M1')
        L = 100
        E = 29000000
        I = 1/12
        beam_sym = beam_mech(L,E,I)
        #point load up
        beam_sym.apply_load(self.R1,0,-1)
        beam_sym.apply_load(self.R2,100,-1)
        beam_sym.apply_load(self.M1,0,-2)
        #beam_sym.apply_load(-80000,0,0, end = 50)
        beam_sym.apply_load(50,5,-2,)
        beam_sym.apply_load(-12,50,-1)
        beam_sym.bc_slope.append((0,0))
        beam_sym.bc_deflection.append((0,0))
        beam_sym.bc_deflection.append((100,0))
        #beam_sym.apply_load(28.4,0,0,end=100)
        print(beam_sym.applied_loads)
        beam_sym.solve_for_reaction_loads(self.R1,self.M1,self.R2)
        print(beam_sym.boundary_conditions)
        print(beam_sym.applied_loads)
        print(beam_sym.reaction_loads)
        print("test")
        print(beam_sym.max_shear_force())
        print(beam_sym.max_deflection())
        self.plot = beam_sym.plot_deflection()
        #beam_sym.plot_shear_force()
        #beam_sym.plot_bending_moment()
        print(beam_sym.load)
        #plot.show()
        