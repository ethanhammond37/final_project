###############################################################################
#                         Calculating stresses in a beam
###############################################################################
#must be running sympy 1.3
#command to change is: conda install sympy=1.3

import numpy as np
import sympy as sp
from sympy.physics.continuum_mechanics.beam import Beam as beam_mech
from sympy.core import S
from sympy.functions import Piecewise, SingularityFunction
#from sympy.physics.continuum_mechanics.beam import Beam3D as beam_mech_3D

###############################################################################
#universal units
units = "IPS"

################################# Beams #######################################
class beam:
    
    def __init__(self,L,H,W,E,dens):
        #L = length, H = height, W = width, I = moment of inertia
        #c = cenroid, E = modulus of elasticity
        #I = moment of inertia about horizonatal centroid
        self.L,self.H,self.W,self.E,self.dens = L,H,W,E,dens
        self.beam_x0, self.beam_x1 = 0,L
        
class rectangular_beam(beam):
    
    def __init__(self, L, H, W, E, dens):
        beam.__init__(self, L,H,W,E,dens)
        self.I = W*(H**3)/12
        self.c = (H/2)
        vol = L*H*W
        self.func = vol*dens/L
        
class rectangular_hallow_beam(beam):
    
    def __init__(self, L, H, W, h, w, E, dens):
        beam.__init__(self, L,H,W,E,dens)
        self.h, self.w = h,w
        self.I = (W)*(H**3)/12 - (w)*(h**3)/12
        vol = L*(H*W-h*w)
        self.func = vol*dens/L
        
        
class circular_beam(beam):
    
    def __init__(self, L, D, E, dens):
        beam.__init__(self, L,D,D,E,dens)
        self.D = D
        self.I = (np.pi*D**4)/64
        self.c = (D/2)
        vol = L*np.pi*(D/2)**2
        self.func = vol*dens/L
        
class circular_hallow_beam(beam):
    
    def __init__(self, L, D, d, E, dens):
        beam.__init__(self, L,D,D,E,dens)
        self.D, self.d = D,d
        self.I = (np.pi*(D**4-d**4))/64
        self.c = (D/2)
        vol = L*np.pi*((D**2)-(d**2))/4
        self.func = vol*dens/L
        
class i_beam(beam):
    
    def __init__(self, L, H, W, th, tw, E, dens):
        beam.__init__(self, L,H,W,E,dens)
        self.th, self.tw = th, tw
        b1 = (W-tw)/2
        d1 = (H-2*th)
        self.I = W*(H**3)/12 - b1*(d1**3)/12
        self.c = (H/2)
        vol = ((2*W*th)+(tw*(H-2*th)))*L
        self.func = vol*dens/L
        
class triangular_beam(beam):
    
    def __init__(self, L, H, W, E, dens):
        beam.__init__(self, L,H,W,E,dens)
        self. I = W*(H**3)/36
        self.c = (H/3)
        vol = L*W*H/2
        self.func = vol*dens/L
        
############################### Supports ######################################
class support:
    
    def __init__(self,x,rotational_freedom,x_freedom,y_freedom):
        self.x, self.r_f = x,rotational_freedom
        self.x_f, self.y_f = x_freedom,y_freedom
        
class fixed_support(support):
    #no axis of freedom
    def __init__(self,x):
        support.__init__(self,x,False,False,False)
        
class y_ball_support(support):
    #can rotate and move up/down
    def __init__(self,x):
        support.__init__(self,x,True,False,True)
        
class y_roller_support(support):
    #can only move up and down
    def __init__(self,x):
        support.__init__(self,x,False,False,True)
    
class x_ball_support(support):
    #can rotate and move left/right
    def __init__(self,x):
        support.__init__(self,x,True,True,False)
        
class x_roller_support(support):
    #can only 
    def __init__(self,x):
        support.__init__(self,x,False,True,False)
        
class pin_support(support):
    
    def __init__(self,x):
        support.__init__(self,x,True,False,False)
        
############################# Loads ###########################################
class load_class:
    
    def __init__(self,start,finish,load_0,load_1,order):
        self.x0, self.x1, self.l0, self.l1,self.order = start, finish, load_0, load_1, order
        
class point_load(load_class):
    
    def __init__(self,x,load):
        load_class.__init__(self,x,x,load,load,-1)
        self.func = load
        
class distributed_load(load_class):
    
    def __init__(self,x0,x1,load_per_unit):
        load_class.__init__(self,x0,x1, load_per_unit, load_per_unit,0)
        self.func = load_per_unit
        
class linear_load(load_class):
    
    def __init__(self,x0,x1,load0,load1):
        load_class.__init__(self,x0,x1,load0,load1,1)
        x = sp.Symbol("x")
        self.func = load0 + (load1-load0)/(x1-x0) * (x - x0)
        
class parabolic_load(load_class):
    
    def __init__(self,x0,x1,load0,load1):
        load_class.__init__(self,x0,x1,load0,load1,2)
        x = sp.Symbol("x")
        if load1 > load0:
            self.func = load0 + ((x-x0)**2)*(load1-load0)/((x1-x0)**2)
            
############################# creating the GUI ################################
import sys
from PyQt5.QtWidgets import (QMainWindow, QDesktopWidget, QApplication, \
                             QWidget, QDialog, QLineEdit, QAction,\
                             QFileDialog,QMenu, QSizePolicy, QPushButton,\
                             QGridLayout, QLabel)
from PyQt5.QtCore import Qt

from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.figure import Figure

class main_GUI(QMainWindow):
    
    def __init__(self,parent = None):
        QMainWindow.__init__(self,parent)
        #init beam
        main_GUI.beam_type = None
        self.beam = None
        #init supports
        self.supports = [None,None,None,None,None]
        #init loads
        self.loads = [None,None,None,None,None]
        #setting title
        self.title = "Beam loading"
        #setting menu
        self.menu_usermade(add = "init")
        #setting the geometry
        self.geometry(1800,1200,0,0)
        #the plot
        self.deflection_plot = Deflection_Plot()
        self.shear_plot = Shear_Plot()
        self.bending_plot = Bending_Plot()
        self.load_plot = Load_Plot()
        #creating buttons
        self.button_clear = QPushButton("Clear")
        
        #clearing
        self.button_update = QPushButton("Update")
        self.button_update.clicked.connect(self.update_all)
        #adding a widget
        self.widget = QDialog()
        self.layout = QGridLayout()
        self.layout.addWidget(self.load_plot,1,1,6,4)
        self.layout.addWidget(self.shear_plot,7,1,6,4)
        self.layout.addWidget(self.bending_plot,13,1,6,4)
        self.layout.addWidget(self.deflection_plot,19,1,6,4)
        for col in range(1,4):
            self.layout.setColumnStretch(col,4)
        self.layout.addWidget(self.button_clear,25,1,1,2,Qt.AlignCenter)
        self.layout.addWidget(self.button_update,25,3,1,2,Qt.AlignCenter)
        self.widget.setLayout(self.layout)
        self.setCentralWidget(self.widget)
        
    def geometry(self,x_initial,y_initial,x_plus,y_plus):
        if x_initial > 0:
            self.width = x_initial
        if y_initial > 0:
            self.height = y_initial
        if x_plus > 0:
            self.width = self.width + x_plus
        if y_plus > 0:
            self.height = self.height + y_plus
        self.setGeometry(0,0,self.width,self.height)
        #getting current screen
        qtscreen = self.frameGeometry()
        centerpoint = QDesktopWidget().availableGeometry().center()
        qtscreen.moveCenter(centerpoint)
        self.move(qtscreen.topLeft())
        self.show()
        
# =============================================================================
# Creating the menue
# =============================================================================
        
    def menu_usermade(self,add=None):
        if add == "init":
            self.menuFile = self.menuBar().addMenu("&File")
            #saveas stuff
            #
            #add menu #########################################################
            self.menuAdd = self.menuBar().addMenu("&Add")       
            #beam menu ###########################################
            self.subBeam = QMenu("&Beam",self)
            #rectangle
            self.actionRect = QAction("Rectangular",self)
            self.actionRect.triggered.connect(self.create_rect_beam)
            #rectangle(hollow)
            self.actionRect_hall = QAction("Rectangular (hallow)",self)
            self.actionRect_hall.triggered.connect(self.create_rect_hall_beam)
            #circular
            self.actionCirc = QAction("Circular",self)
            self.actionCirc.triggered.connect(self.create_circular_beam)
            #circular(hollow)
            self.actionCirc_hall = QAction("Circular (hallow)",self)
            self.actionCirc_hall.triggered.connect(self.create_circular_hall_beam)
            #I beam
            self.actionI = QAction("I-beam",self)
            self.actionI.triggered.connect(self.create_i_beam)
            #triangular beam
            self.actionTria = QAction("Triangular",self)
            self.actionTria.triggered.connect(self.create_tria_beam)
            self.subBeam.addActions([self.actionRect,self.actionRect_hall, \
                                      self.actionCirc,self.actionCirc_hall, \
                                      self.actionI,self.actionTria])
            self.menuAdd.addMenu(self.subBeam)
            #support menu ##########################################
            self.subSupport = QMenu("&Support",self)
            self.subsubLim_x = QMenu("Limit horizontal movement", self)
            self.subsubLim_y = QMenu("Limit verticle movement",self)
            self.subsubLim_both = QMenu("Limit horizontal and verticle "+ \
                                        "movement", self)
            #fixed support
            self.actionFixed = QAction("Fixed",self)
            self.actionFixed.triggered.connect(self.create_fixed_support)
            #pin support
            self.actionPin = QAction("Pin",self)
            self.actionPin.triggered.connect(self.create_pin_support)
            #yBall support
            self.actionBall_y = QAction("Ball (can rotate)", self)
            self.actionBall_y.triggered.connect(self.create_y_ball_support)
            #yRoll support
            self.actionRoll_y = QAction("Roller (cannot rotate)",self)
            self.actionRoll_y.triggered.connect(self.create_y_roller_support)
            #xBall support
            self.actionBall_x = QAction("Ball (can rotate)",self)
            self.actionBall_x.triggered.connect(self.create_x_ball_support)
            #XBall
            self.actionRoll_x = QAction("Roller (cannot rotate)",self)
            self.actionRoll_x.triggered.connect(self.create_x_roller_support)
            #sub menues
            self.subSupport.addMenu(self.subsubLim_x)
            self.subsubLim_x.addActions([self.actionBall_y,self.actionRoll_y])
            self.subSupport.addMenu(self.subsubLim_y)
            self.subsubLim_y.addActions([self.actionBall_x,self.actionRoll_x])
            self.subSupport.addMenu(self.subsubLim_both)
            self.subsubLim_both.addActions([self.actionFixed,self.actionPin])
            self.menuAdd.addMenu(self.subSupport)
            #load menu ################################################
            self.subLoad = QMenu("&Load",self)
            #point load
            self.actionPoint = QAction("Point",self)
            self.actionPoint.triggered.connect(self.create_point_load)
            #distributed load
            self.actionDist = QAction("Distributed",self)
            self.actionDist.triggered.connect(self.create_distributed_load)
            #linear load
            self.actionLinear = QAction("Linear",self)
            self.actionLinear.triggered.connect(self.create_linear_load)
            #parabolic load
            self.actionPar = QAction("Parabolic",self)
            self.actionPar.triggered.connect(self.create_parabolic_load)
            self.subLoad.addActions([self.actionPoint,self.actionDist,\
                                    self.actionLinear,self.actionPar])
            self.menuAdd.addMenu(self.subLoad)
            #remove menu ######################################################
            self.menuRemove = self.menuBar().addMenu("&Remove")
            #supports
            self.actionRemoveS1 = QAction("Support 1",self)
            self.actionRemoveS1.triggered.connect(self.remove_support_1)
            self.actionRemoveS2 = QAction("Support 2",self)
            self.actionRemoveS2.triggered.connect(self.remove_support_2)
            self.actionRemoveS3 = QAction("Support 3",self)
            self.actionRemoveS3.triggered.connect(self.remove_support_3)
            self.actionRemoveS4 = QAction("Support 4",self)
            self.actionRemoveS4.triggered.connect(self.remove_support_4)
            self.actionRemoveS5 = QAction("Support 5",self)
            self.actionRemoveS5.triggered.connect(self.remove_support_5)
            #loads
            self.actionRemoveL1 = QAction("Load 1",self)
            self.actionRemoveL1.triggered.connect(self.remove_load_1)
            self.actionRemoveL2 = QAction("Load 2",self)
            self.actionRemoveL2.triggered.connect(self.remove_load_2)
            self.actionRemoveL3 = QAction("Load 3",self)
            self.actionRemoveL3.triggered.connect(self.remove_load_3)
            self.actionRemoveL4 = QAction("Load 4",self)
            self.actionRemoveL4.triggered.connect(self.remove_load_4)
            self.actionRemoveL5 = QAction("Load 5",self)
            self.actionRemoveL5.triggered.connect(self.remove_load_5)
            self.menuRemove.addActions([self.actionRemoveS1,\
                                        self.actionRemoveS2,\
                                        self.actionRemoveS3,\
                                        self.actionRemoveS4,\
                                        self.actionRemoveS5,\
                                        self.actionRemoveL1,\
                                        self.actionRemoveL2,\
                                        self.actionRemoveL3,\
                                        self.actionRemoveL4,\
                                        self.actionRemoveL5])
        elif add != None:
            exec(add)
            
    def get_distance(self):
        return("inches")
        
    def get_force(self):
        return("lbs")
        
    def get_stress(self):
        return("psi")
            
    def update_all(self):
        #updating the beam
        if main_GUI.beam_type != None:
            exec("self.update_"+main_GUI.beam_type+"()")
            #updating the supports
            for sup in self.supports:
                if sup != None:
                    exec(sup)
            for load in self.loads:
                if load != None:
                    exec(load)
            #setting up the equations
            self.singularity_work()
            #solving for the supports and showing them on the GUI
            self.update_supports()
            #pulls singularity equations
            self.singularity_stats()
            #recreates plots
            self.deflection_plot.redraw(self.beam.L,self.deflection_sym)
            self.shear_plot.redraw(self.beam.L,self.shear_sym)
            self.bending_plot.redraw(self.beam.L,self.bending_sym)
            self.load_plot.redraw(self.beam.L,self.load_sym)

    def singularity_work(self):
        #beam
        self.beam_sym = beam_mech(self.beam.L,self.beam.E,self.beam.I)
        #its own weight
        self.beam_sym.apply_load(-self.beam.func,0,0,self.beam.L)
        #loads
        for i in range(0,len(self.loads)):
            if self.loads[i] != None:
                num = str(int(i + 1))
                if eval("self.load_"+num+".order") < 0:
                    self.beam_sym.apply_load(-eval("self.load_"+num+".l0"),\
                                             eval("self.load_"+num+".x0"),\
                                             eval("self.load_"+num+".order"))
                else:
                    self.beam_sym.apply_load(-eval("self.load_"+num+".l0"),\
                                             eval("self.load_"+num+".x0"),\
                                             eval("self.load_"+num+".order"),\
                                             end=eval("self.load_"+num+".x1"))
                    
##### add linear and parabolid loads                    
                    
                    
                    
                    
                    
                    
                    
        #supports
        self.supports_reactions = []
        self.supports_moments = []
        for i in range(0,len(self.supports)):
            if self.supports[i] != None:
                num = str(int(i+1))
                #load limiting y axis movement
                if eval("self.support_"+num+".y_f") == False:
                    exec("self.R"+num+"=sp.symbols('R"+num+"')")
                    exec("self.supports_reactions.append('self.R"+num+"')")
                    exec("self.beam_sym.apply_load((self.R"+num+\
                         "),self.support_"+num+".x,-1)")
                    exec("self.beam_sym.bc_deflection.append((self.support_"+\
                                                              num+".x,0))")
                if eval("self.support_"+num+".r_f") == False:
                    exec("self.M"+num+"=sp.symbols('M"+num+"')")
                    exec("self.supports_moments.append('self.M"+num+"')")
                    exec("self.beam_sym.apply_load((self.M"+num+\
                         "),self.support_"+num+".x,-2)")
                    exec("self.beam_sym.bc_slope.append((self.support_"+num+\
                                                         ".x,0))")
        exec("self.beam_sym.solve_for_reaction_loads("+\
             str(self.supports_reactions).replace("'","")[1:-1]+","+\
             str(self.supports_moments).replace("'","")[1:-1]+")")
        self.reactions = self.beam_sym.reaction_loads
        
    def singularity_stats(self):
        self.deflection_sym = self.beam_sym.deflection()
        self.bending_sym = self.beam_sym.bending_moment()
        self.shear_sym = self.beam_sym.shear_force()
        self.slope_sym = self.beam_sym.slope()
        self.load_sym = self.beam_sym.load
        
# =============================================================================
# creating beam widgets        
# =============================================================================
    def create_sub(self,obj,row_start,col_start,row_end,col_end,*args):
        #removes old widgets
        if self.layout.itemAtPosition(row_start,col_start) != None:
            for row in range(row_start,row_end + 1):
                for col in range(col_start,col_end + 1):
                    if self.layout.itemAtPosition(row,col) != None:
                        self.layout.itemAtPosition(row,col).widget().deleteLater()
        #adds new widgets
        i = 1
        for row in range(row_start,row_end + 1):
            for col in range(col_start, col_end + 1):
                if i < len(args):
                    if row == row_start:
                        if col == col_start:
                            exec("self.label"+str(row)+"_"+str(col)+\
                                 "=QLabel('')")
                            exec("self.layout.addWidget(self.label"+str(row)+\
                                                        "_"+str(col)+","+\
                                                        str(row)+","+str(col)+\
                                                        ",1,10,Qt.AlignCenter)")
                    else:
                        if col % 2 == 1:
                            exec("self.label"+str(row)+"_"+str(col)+"=QLabel('')")
                            exec("self.layout.addWidget(self.label"+str(row)+\
                                                        "_"+str(col)+","+\
                                                        str(row)+","+str(col)+\
                                                        ",1,1,Qt.AlignRight)")
                        else:
                            exec("self.box"+str(row)+"_"+str(col)+\
                                 "=QLineEdit('')")
                            exec("self.layout.addWidget(self.box"+str(row)+\
                                                        "_"+str(col)+","+\
                                                        str(row)+","+str(col)+\
                                                        ",1,1)")
                            i = i + 1
        #gives the widgets values and calling functions
        i = 1
        for row in range(row_start,row_end + 1):
            for col in range(col_start, col_end + 1):
                if i < len(args):
                    if row == row_start:
                        if col == col_start:
                            exec("self.label"+str(row)+"_"+str(col)+\
                                 ".setText(args[0])")
                    else:
                        if col % 2 == 1:
                            exec("self.label"+str(row)+"_"+str(col)+\
                                 ".setText(args[i]+':')")
                        else:
                            exec("self."+obj+"_"+args[i]+\
                                 "= lambda self: self.box"+str(row)+"_"+\
                                 str(col)+".text()")
                            i = i + 1
              
    def create_rect_beam(self):
        main_GUI.beam_type = "rect_beam"
        self.create_sub("beam",1,5,3,14,"Rectangular Beam","Length","Height",\
                        "Width","Youngs_Modulus","Specific_Weight")
        
    def update_rect_beam(self):
        self.beam = rectangular_beam(float(self.beam_Length(self)),\
                                     float(self.beam_Height(self)),\
                                     float(self.beam_Width(self)),\
                                     float(self.beam_Youngs_Modulus(self)),\
                                     float(self.beam_Specific_Weight(self)))
        
    def create_rect_hall_beam(self):
        main_GUI.beam_type = "rect_hall_beam"
        self.create_sub("beam",1,5,3,14,"Rectangular (Hollow) Beam","Length",\
                        "Height","Width","Interior_Width","Interior_Height",\
                        "Youngs_Modulus","Specific_Weight")    
        
    def update_rect_hall_beam(self):
        self.beam = rectangular_hallow_beam(float(self.beam_Length(self)),\
                                     float(self.beam_Height(self)),\
                                     float(self.beam_Width(self)),\
                                     float(self.beam_Interior_Height(self)),\
                                     float(self.beam_Interior_Width(self)),\
                                     float(self.beam_Youngs_Modulus(self)),\
                                     float(self.beam_Specific_Weight(self)))
        
    def create_circular_beam(self):
        main_GUI.beam_type = "circular_beam"
        self.create_sub("beam",1,5,3,14,"Circular Beam","Length","Diameter",\
                        "Youngs_Modulus","Specific_Weight")    
        
    def update_circular_beam(self):
        self.beam = circular_beam(float(self.beam_Length(self)),\
                                     float(self.beam_Diameter(self)),\
                                     float(self.beam_Youngs_Modulus(self)),\
                                     float(self.beam_Specific_Weight(self)))
        
    def create_circular_hall_beam(self):
        main_GUI.beam_type = "circular_hall_beam"
        self.create_sub("beam",1,5,3,14,"Circular (Hallow) Beam","Length",\
                        "Diameter","Interior_Diameter","Youngs_Modulus",\
                        "Specific_Weight")
        
    def update_circular_hall_beam(self):
        self.beam = circular_hallow_beam(float(self.beam_Length(self)),\
                                     float(self.beam_Diameter(self)),\
                                     float(self.beam_Interior_Diameter(self)),\
                                     float(self.beam_Youngs_Modulus(self)),\
                                     float(self.beam_Specific_Weight(self)))
        
    def create_i_beam(self):
        main_GUI.beam_type = "i_beam"
        self.create_sub("beam",1,5,3,14,"I Beam","Length","Height","Width",\
                        "Top_Bottom_Tickness","Center_Thickness",\
                        "Youngs_Modulus","Specific_Weight")
        
    def update_i_beam(self):
        self.beam = i_beam(float(self.beam_Length(self)),\
                           float(self.beam_Height(self)),\
                           float(self.beam_Width(self)),\
                           float(self.beam_Top_Bottom_Thickness(self)),\
                           float(self.beam_Center_Thickness(self)),\
                           float(self.beam_Youngs_Modulus(self)),\
                           float(self.beam_Specific_Weight(self)))
        
    def create_tria_beam(self):
        main_GUI.beam_type = "tria_beam"
        self.create_sub("beam",1,5,3,14,"Rectangular Beam","Length","Height",\
                        "Width","Youngs_Modulus","Specific_Weight")
        
    def update_tria_beam(self):
        self.beam = triangular_beam(float(self.beam_Length(self)),\
                                     float(self.beam_Height(self)),\
                                     float(self.beam_Width(self)),\
                                     float(self.beam_Youngs_Modulus(self)),\
                                     float(self.beam_Specific_Weight(self)))
        
# =============================================================================
# creating support widgets        
# =============================================================================
        
    def first_empty_row(self,row_start,row_end,col):
        for row in range(row_start,row_end):
            if self.layout.itemAtPosition(row,col) == None:
                return row
        return None
    
    def remove_support_1(self):
        self.create_sub("remove",4,5,5,14)
        self.supports[0] = None
        
    def remove_support_2(self):
        self.create_sub("remove",6,5,7,14)
        self.supports[1] = None 
    
    def remove_support_3(self):
        self.create_sub("remove",8,5,9,14)
        self.supports[2] = None
     
    def remove_support_4(self):
        self.create_sub("remove",10,5,11,14)
        self.supports[2] = None
        
    def remove_support_5(self):
        self.create_sub("remove",12,5,13,14)
        self.supports[2] = None
        
        
    def update_supports(self):
        for i in range(0,len(self.supports)):
            num = str(int(i+1))
            row = int(i*2+5)
            col = int(7)
            if self.supports[i] != None:
                if eval("self.support_"+num+".y_f") == False:
                    self.create_sub("remove",row,col,row,col+1)
                    force = eval("self.reactions[self.R"+num+"]")
                    exec("self.label"+str(row)+"_"+str(col)+"=QLabel('Reaction_Force:')")
                    exec("self.layout.addWidget(self.label"+str(row)+"_"+str(col)+","+str(row)+","+str(col)+",1,1,Qt.AlignRight)")
                    col = col + 1
                    exec("self.label"+str(row)+"_"+str(col)+"=QLabel('"+str(round(force,2))+"')")
                    exec("self.layout.addWidget(self.label"+str(row)+"_"+str(col)+","+str(row)+","+str(col)+",1,1)")
                    col = col + 1
                if eval("self.support_"+num+".r_f") == False:
                    self.create_sub("remove",row,col,row,col+1)
                    moment = eval("self.reactions[self.M"+num+"]")
                    exec("self.label"+str(row)+"_"+str(col)+"=QLabel('Reaction_Moment:')")
                    exec("self.layout.addWidget(self.label"+str(row)+"_"+str(col)+","+str(row)+","+str(col)+",1,1,Qt.AlignRight)")
                    col = col + 1
                    exec("self.label"+str(row)+"_"+str(col)+"=QLabel('"+str(round(moment,2))+"')")
                    exec("self.layout.addWidget(self.label"+str(row)+"_"+str(col)+","+str(row)+","+str(col)+",1,1)")
                    col = col + 1
        
        
    def create_fixed_support(self):
        row = self.first_empty_row(4,13,5)
        if row != None:
            num = str(int(row/2-1))
            self.create_sub("support"+"_"+num,row,5,row+1,14,\
                            "Support_"+num+" Fixed Support","Location_x")
            self.supports[int(num)-1] = ("self.support_"+num+\
                          "= fixed_support(float(self.support_"+num+\
                          "_Location_x(self)))")
            
    def create_y_ball_support(self):
        row = self.first_empty_row(4,13,5)
        if row != None:
            num = str(int(row/2-1))
            self.create_sub("support"+"_"+num,row,5,row+1,14,\
                            "Support_"+num+" Verticle Ball Support",\
                            "Location_x")
            self.supports[int(num)-1] = ("self.support_"+num+\
                          "= y_ball_support(float(self.support_"+num+\
                          "_Location_x(self)))")
            
    def create_y_roller_support(self):
        row = self.first_empty_row(4,13,5)
        if row != None:
            num = str(int(row/2-1))
            self.create_sub("support"+"_"+num,row,5,row+1,14,\
                            "Support_"+num+" Verticle Roller Support",\
                            "Location_x")
            self.supports[int(num)-1] = ("self.support_"+num+\
                          "= y_roller_support(float(self.support_"+num+\
                          "_Location_x(self)))")
            
    def create_x_ball_support(self):
        row = self.first_empty_row(4,13,5)
        if row != None:
            num = str(int(row/2-1))
            self.create_sub("support"+"_"+num,row,5,row+1,14,\
                            "Support_"+num+" Horizontal Ball Support",\
                            "Location_x")
            self.supports[int(num)-1] = ("self.support_"+num+\
                          "= x_ball_support(float(self.support_"+num+\
                          "_Location_x(self)))")
            
    def create_x_roller_support(self):
        row = self.first_empty_row(4,13,5)
        if row != None:
            num = str(int(row/2-1))
            self.create_sub("support"+"_"+num,row,5,row+1,14,\
                            "Support_"+num+" Horizontal Roller Support",\
                            "Location_x")
            self.supports[int(num)-1] = ("self.support_"+num+\
                          "= x_roller_support(float(self.support_"+num+\
                          "_Location_x(self)))")
            
    def create_pin_support(self):
        row = self.first_empty_row(4,13,5)
        if row != None:
            num = str(int(row/2-1))
            self.create_sub("support"+"_"+num,row,5,row+1,14,\
                            "Support_"+num+" Pin Support","Location_x")
            self.supports[int(num)-1] = ("self.support_"+num+\
                          "= pin_support(float(self.support_"+num+\
                          "_Location_x(self)))")
            
# =============================================================================
# creating load widgets        
# ============================================================================
    
    def remove_load_1(self):
        self.create_sub("remove",14,5,15,14)
        self.loads[0] = None
        
    def remove_load_2(self):
        self.create_sub("remove",16,5,17,14)
        self.loads[1] = None
        
    def remove_load_3(self):
        self.create_sub("remove",18,5,19,14)
        self.loads[2] = None
        
    def remove_load_4(self):
        self.create_sub("remove",20,5,21,14)
        self.loads[3] = None
        
    def remove_load_5(self):
        self.create_sub("remove",22,5,23,14)
        self.loads[0] = None
        
    def create_point_load(self):
        row = self.first_empty_row(14,23,5)
        if row != None:
            num = str(int(row/2-6))
            self.create_sub("load_"+num,row,5,row+1,14,"Load_"+num+\
                            " Point Load","Location_x","Load")
            self.loads[int(num)-1] = ("self.load_"+num+\
                       " = point_load(float(self.load_"+num+\
                       "_Location_x(self)),float(self.load_"+num+\
                       "_Load(self)))")
            
    def create_distributed_load(self):
        row = self.first_empty_row(14,23,5)
        if row != None:
            num = str(int(row/2-6))
            self.create_sub("load_"+num,row,5,row+1,14,"Load_"+num+\
                            " Distributed Load","Start_x","End_x","Load")
            self.loads[int(num)-1] = ("self.load_"+num+\
                       " = distributed_load(float(self.load_"+num+\
                       "_Start_x(self)),float(self.load_"+num+\
                       "_End_x(self)),float(self.load_"+num+\
                       "_Load(self)))")
            
    def create_linear_load(self):
        row = self.first_empty_row(14,23,5)
        if row != None:
            num = str(int(row/2-6))
            self.create_sub("load_"+num,row,5,row+1,14,"Load_"+num+\
                            " Linear Load","Start_x","End_x","Load_start","Load_end")
            self.loads[int(num)-1] = ("self.load_"+num+\
                       " = linear_load(float(self.load_"+num+\
                       "_Start_x(self)),float(self.load_"+num+\
                       "_End_x(self)),float(self.load_"+num+\
                       "_Load_start(self)),float(self.load_"+num+\
                       "_Load_end(self)))")
            
    def create_parabolic_load(self):
        row = self.first_empty_row(14,23,5)
        if row != None:
            num = str(int(row/2-6))
            self.create_sub("load_"+num,row,5,row+1,14,"Load_"+num+\
                            " Parabolic Load","Start_x","End_x","Load_start","Load_end")
            self.loads[int(num)-1] = ("self.load_"+num+\
                       " = parabolic_load(float(self.load_"+num+\
                       "_Start_x(self)),float(self.load_"+num+\
                       "_End_x(self)),float(self.load_"+num+\
                       "_Load_start(self)),float(self.load_"+num+\
                       "_Load_end(self)))")
            

class Deflection_Plot(FigureCanvas):
    
    def __init__(self):
        self.fig = Figure()
        self.deflection_plot = self.fig.add_subplot(111)
        self.deflection_plot.clear()
        distance = main_GUI.get_distance(main_GUI)
        self.deflection_plot.set_xlabel("x ("+distance+")")
        self.deflection_plot.set_ylabel("y ("+distance+")")
        self.deflection_plot.set_title("Deflection")
        FigureCanvas.__init__(self,self.fig)
        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
    
    def redraw(self, l, y):
        n = 500
        self.deflection_plot.clear()
        distance = main_GUI.get_distance(main_GUI)
        x = sp.symbols('x')
        x_array = np.linspace(0,l,n)
        x_array[n-1] = x_array[n-1] - 10**-7
        y_array = [y.subs(x,p) for p in x_array]
        y_max,y_min = max(y_array),min(y_array)
        if y_max > (-y_min): 
            deflection = y_max
        else:
            deflection = y_min
        y_zero = [0.0 for p in x_array]
        self.deflection_plot.plot(x_array,y_array,'r',x_array,y_zero,'k')
        self.deflection_plot.set_title(str(main_GUI.beam_type) +" Deflection (max deflection = "+str(round(deflection,6))+" "+distance+")")
        self.draw()
        
class Shear_Plot(FigureCanvas):
    
    def __init__(self):
        self.fig = Figure()
        self.shear_plot = self.fig.add_subplot(111)
        self.shear_plot.clear()
        distance = main_GUI.get_distance(main_GUI)
        force = main_GUI.get_force(main_GUI)
        self.shear_plot.set_xlabel("x ("+distance+")")
        self.shear_plot.set_ylabel("y ("+force+")")
        self.shear_plot.set_title("Shear")
        FigureCanvas.__init__(self,self.fig)
        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
    
    def redraw(self, l, y):
        n = 500
        self.shear_plot.clear()
        force = main_GUI.get_force(main_GUI)
        distance = main_GUI.get_distance(main_GUI)
        x = sp.symbols('x')
        x_array = np.linspace(0,l,n)
        x_array[n-1] = x_array[n-1] - 10**-7
        y_array = [y.subs(x,p) for p in x_array]
        y_max,y_min = max(y_array),min(y_array)
        if y_min == float('-inf'):
            shear = y_max
        elif y_max == float('+inf'):
            shear = y_min
        elif y_max > (-y_min):
            shear = y_max
        else:
            shear = y_min
        y_zero = [0.0 for p in x_array]
        self.shear_plot.plot(x_array,y_array,'r',x_array,y_zero,'k')
        self.shear_plot.set_title(str(main_GUI.beam_type) +" Shear (max shear = "+str(round(shear,6))+" "+force+")")
        self.shear_plot.set_xlabel("x ("+distance+")")
        self.shear_plot.set_ylabel("y ("+force+")")
        self.draw()
        
class Bending_Plot(FigureCanvas):
    
    def __init__(self):
        self.fig = Figure()
        self.bending_plot = self.fig.add_subplot(111)
        self.bending_plot.clear()
        distance = main_GUI.get_distance(main_GUI)
        force = main_GUI.get_force(main_GUI)
        self.bending_plot.set_xlabel("x ("+distance+")")
        self.bending_plot.set_ylabel("y ("+force+"-"+distance+")")
        self.bending_plot.set_title("Bending")
        FigureCanvas.__init__(self,self.fig)
        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
    
    def redraw(self, l, y):
        n = 500
        self.bending_plot.clear()
        force = main_GUI.get_force(main_GUI)
        distance = main_GUI.get_distance(main_GUI)
        x = sp.symbols('x')
        x_array = np.linspace(0,l,n)
        x_array[n-1] = x_array[n-1] - 10**-7
        y_array = [y.subs(x,p) for p in x_array]
        y_max,y_min = max(y_array),min(y_array)
        if y_min == float('-inf'):
            bending = y_max
        elif y_max == float('+inf'):
            bending = y_min
        elif y_max > (-y_min):
            bending = y_max
        else:
            bending = y_min
        y_zero = [0.0 for p in x_array]
        self.bending_plot.plot(x_array,y_array,'r',x_array,y_zero,'k')
        self.bending_plot.set_title(str(main_GUI.beam_type) +" Bending (max moment = "+str(round(bending,4))+" "+force+"-"+distance+")")
        self.bending_plot.set_xlabel("x ("+distance+")")
        self.bending_plot.set_ylabel("y ("+force+"-"+distance+")")
        self.draw()
        
class Load_Plot(FigureCanvas):
    
    def __init__(self):
        self.fig = Figure()
        self.load_plot = self.fig.add_subplot(111)
        self.load_plot.clear()
        distance = main_GUI.get_distance(main_GUI)
        force = main_GUI.get_force(main_GUI)
        self.load_plot.set_xlabel("x ("+distance+")")
        self.load_plot.set_ylabel("y ("+force+")")
        self.load_plot.set_title("Loading")
        FigureCanvas.__init__(self,self.fig)
        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
    
    def redraw(self, l, y):
        n = 500
        self.load_plot.clear()
        force = main_GUI.get_force(main_GUI)
        distance = main_GUI.get_distance(main_GUI)
        x = sp.symbols('x')
        x_array = np.linspace(0,l,n)
        y_array = [y.subs(x,p) for p in x_array]
        print(y_array[0])
        for i in range(0,len(y_array)):
            if str(y_array[i]) == 'nan':
                y_array[i] = 0.0
        y_array[0],y_array[n-1] = 0.0,0.0
        y_max,y_min = max(y_array),min(y_array)
        if y_min == float('-inf'):
            load = y_max
        elif y_max == float('+inf'):
            load = y_min
        elif y_max > (-y_min):
            load = y_max
        else:
            load = y_min
        y_zero = [0.0 for p in x_array]
        self.load_plot.plot(x_array,y_array,'r',x_array,y_zero,'k')
        self.load_plot.set_title(str(main_GUI.beam_type) +" Loading (max load = "+str(round(load,4))+" "+force+")")
        self.load_plot.set_xlabel("x ("+distance+")")
        self.load_plot.set_ylabel("y ("+force+")")
        self.draw()
        
if __name__ == "__main__":
    app = QApplication(sys.argv)
    ex = main_GUI()
    app.exec_()
    
#FigureCanvas.printJpeg
