###############################################################################
#                         Calculating stresses in a beam
###############################################################################
import numpy as np
import sympy as sp

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
        self.func = vol*dens
        
class rectangular_hallow_beam(beam):
    
    def __init__(self, L, H, W, h, w, E, dens):
        beam.__init__(self, L,H,W,E,dens)
        self.h, self.w = h,w
        self.I = (W)*(H**3)/12 - (w)*(h**3)/12
        vol = L*(H*W-h*w)
        self.func = vol*dens
        
        
class circular_beam(beam):
    
    def __init__(self, L, D, E, dens):
        beam.__init__(self, L,D,D,E,dens)
        self.D = D
        self.I = (np.pi*D**4)/64
        self.c = (D/2)
        vol = L*np.pi*(D/2)**2
        self.func = vol*dens
        
class circular_hallow_beam(beam):
    
    def __init__(self, L, D, d, E, dens):
        beam.__init__(self, L,D,D,E,dens)
        self.D, self.d = D,d
        self.I = (np.pi*(D**4-d**4))/64
        self.c = (D/2)
        vol = L*np.pi*((D**2)-(d**2))/4
        self.func = vol*dens
        
class i_beam(beam):
    
    def __init__(self, L, H, W, th, tw, E, dens):
        beam.__init__(self, L,H,W,E,dens)
        self.th, self.tw = th, tw
        b1 = (W-tw)/2
        d1 = (H-2*th)
        self.I = W*(H**3)/12 - b1*(d1**3)/12
        self.c = (H/2)
        vol = ((2*W*th)+(tw*(H-2*th)))*L
        self.func = vol*dens
        
class triangular_beam(beam):
    
    def __init__(self, L, H, W, E, dens):
        beam.__init__(self, L,H,W,E,dens)
        self. I = W*(H**3)/36
        self.c = (H/3)
        vol = L*W*H/2
        self.func = vol*dens
        
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
    
    def __init__(self,start,finish,load_0,load_1):
        self.x0, self.x1, self.l0, self.l1 = start, finish, load_0, load_1
        
class point_load(load_class):
    
    def __init__(self,x,load):
        load_class.__init__(self,x,x,load,load)
        self.func = load
        
class distributed_load(load_class):
    
    def __init__(self,x0,x1,load_per_unit):
        load_class.__init__(self,x0,x1, load_per_unit, load_per_unit)
        self.func = load_per_unit
        
class linear_load(load_class):
    
    def __init__(self,x0,x1,load0,load1):
        load_class.__init__(self,x0,x1,load0,load1)
        x = sp.Symbol("x")
        self.func = load0 + (load1-load0)/(x1-x0) * (x - x0)
        
class parabolic_load(load_class):
    
    def __init__(self,x0,x1,load0,load1):
        load_class.__init__(self,x0,x1,load0,load1)
        x = sp.Symbol("x")
        if load1 > load0:
            self.func = load0 + ((x-x0)**2)*(load1-load0)/((x1-x0)**2)
            
############################# creating the GUI ################################
import sys
from PyQt5.QtWidgets import (QMainWindow, QDesktopWidget, QApplication, \
                             QWidget, QDialog, QLineEdit, QAction,\
                             QMessageBox,QFileDialog,QMenu, QSizePolicy, \
                             QComboBox, QPushButton, QGridLayout, QLabel)
from PyQt5.QtCore import Qt

from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.figure import Figure

class main_GUI(QMainWindow):
    
    def __init__(self,parent = None):
        QMainWindow.__init__(self,parent)
        #init beam
        self.beam = None
        #init supports
        self.supports = [None,None,None]
        #init loads
        self.loads = [None,None,None,None,None]
        #setting title
        self.title = "Beam loading"
        #setting menu
        self.menu_usermade(add = "init")
        #setting the geometry
        self.geometry(1800,600,0,0)
        #the plot
        self.plot = MatplotlibCanvas()
        #creating buttons
        self.button_clear = QPushButton("Clear")
        
        #clearing
        self.button_update = QPushButton("Update")
        self.button_update.clicked.connect(self.update_all)
        #adding a widget
        self.widget = QDialog()
        self.layout = QGridLayout()
        self.layout.addWidget(self.plot,1,1,20,4)
        for col in range(1,4):
            self.layout.setColumnStretch(col,4)
        self.layout.addWidget(self.button_clear,21,1,1,2,Qt.AlignCenter)
        self.layout.addWidget(self.button_update,21,3,1,2,Qt.AlignCenter)
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
                                        self.actionRemoveL1,\
                                        self.actionRemoveL2,\
                                        self.actionRemoveL3,\
                                        self.actionRemoveL4,\
                                        self.actionRemoveL5])
        elif add != None:
            exec(add)
            
    def update_all(self):
        #updating the beam
        if self.beam != None:
            exec("self.update_"+self.beam_type+"()")
        #updating the supports
        for sup in self.supports:
            if sup != None:
                exec(sup)
        for load in self.loads:
            if load != None:
                exec(load)       
        
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
        self.beam_type = "rect_beam"
        self.create_sub("beam",1,5,3,14,"Rectangular Beam","Length","Height",\
                        "Width","Youngs_Modulus","Specific_Gravity")
        
    def update_rect_beam(self):
        self.beam = rectangular_beam(float(self.beam_Length(self)),\
                                     float(self.beam_Height(self)),\
                                     float(self.beam_Width(self)),\
                                     float(self.beam_Youngs_Modulus(self)),\
                                     float(self.beam_Specific_Gravity(self)))
        
    def create_rect_hall_beam(self):
        self.beam_type = "rect_hall_beam"
        self.create_sub("beam",1,5,3,14,"Rectangular (Hollow) Beam","Length",\
                        "Height","Width","Interior_Width","Interior_Height",\
                        "Youngs_Modulus","Specific_Gravity")    
        
    def update_rect_hall_beam(self):
        self.beam = rectangular_hallow_beam(float(self.beam_Length(self)),\
                                     float(self.beam_Height(self)),\
                                     float(self.beam_Width(self)),\
                                     float(self.beam_Interior_Height(self)),\
                                     float(self.beam_Interior_Width(self)),\
                                     float(self.beam_Youngs_Modulus(self)),\
                                     float(self.beam_Specific_Gravity(self)))
        
    def create_circular_beam(self):
        self.beam_type = "circular_beam"
        self.create_sub("beam",1,5,3,14,"Circular Beam","Length","Diameter",\
                        "Youngs_Modulus","Specific_Gravity")    
        
    def update_circular_beam(self):
        self.beam = circular_beam(float(self.beam_Length(self)),\
                                     float(self.beam_Diameter(self)),\
                                     float(self.beam_Youngs_Modulus(self)),\
                                     float(self.beam_Specific_Gravity(self)))
        
    def create_circular_hall_beam(self):
        self.beam_type = "circular_hall_beam"
        self.create_sub("beam",1,5,3,14,"Circular (Hallow) Beam","Length",\
                        "Diameter","Interior_Diameter","Youngs_Modulus",\
                        "Specific_Gravity")
        
    def update_circular_hall_beam(self):
        self.beam = circular_hallow_beam(float(self.beam_Length(self)),\
                                     float(self.beam_Diameter(self)),\
                                     float(self.beam_Interior_Diameter(self)),\
                                     float(self.beam_Youngs_Modulus(self)),\
                                     float(self.beam_Specific_Gravity(self)))
        
    def create_i_beam(self):
        self.beam_type = "i_beam"
        self.create_sub("beam",1,5,3,14,"I Beam","Length","Height","Width",\
                        "Top_Bottom_Tickness","Center_Thickness",\
                        "Youngs_Modulus","Specific_Gravity")
        
    def update_i_beam(self):
        self.beam = i_beam(float(self.beam_Length(self)),\
                           float(self.beam_Height(self)),\
                           float(self.beam_Width(self)),\
                           float(self.beam_Top_Bottom_Thickness(self)),\
                           float(self.beam_Center_Thickness(self)),\
                           float(self.beam_Youngs_Modulus(self)),\
                           float(self.beam_Specific_Gravity(self)))
        
    def create_tria_beam(self):
        self.beam_type = "tria_beam"
        self.create_sub("beam",1,5,3,14,"Rectangular Beam","Length","Height",\
                        "Width","Youngs_Modulus","Specific_Gravity")
        
    def update_tria_beam(self):
        self.beam = triangular_beam(float(self.beam_Length(self)),\
                                     float(self.beam_Height(self)),\
                                     float(self.beam_Width(self)),\
                                     float(self.beam_Youngs_Modulus(self)),\
                                     float(self.beam_Specific_Gravity(self)))
        
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
        
    def create_fixed_support(self):
        row = self.first_empty_row(4,9,5)
        if row != None:
            num = str(int(row/2-1))
            self.create_sub("support"+"_"+num,row,5,row+1,14,\
                            "Support_"+num+" Fixed Support","Location_x")
            self.supports[int(num)-1] = ("self.support_"+num+\
                          "= fixed_support(float(self.support_"+num+\
                          "_Location_x(self)))")
            
    def create_y_ball_support(self):
        row = self.first_empty_row(4,9,5)
        if row != None:
            num = str(int(row/2-1))
            self.create_sub("support"+"_"+num,row,5,row+1,14,\
                            "Support_"+num+" Verticle Ball Support",\
                            "Location_x")
            self.supports[int(num)-1] = ("self.support_"+num+\
                          "= y_ball_support(float(self.support_"+num+\
                          "_Location_x(self)))")
            
    def create_y_roller_support(self):
        row = self.first_empty_row(4,9,5)
        if row != None:
            num = str(int(row/2-1))
            self.create_sub("support"+"_"+num,row,5,row+1,14,\
                            "Support_"+num+" Verticle Roller Support",\
                            "Location_x")
            self.supports[int(num)-1] = ("self.support_"+num+\
                          "= y_roller_support(float(self.support_"+num+\
                          "_Location_x(self)))")
            
    def create_x_ball_support(self):
        row = self.first_empty_row(4,9,5)
        if row != None:
            num = str(int(row/2-1))
            self.create_sub("support"+"_"+num,row,5,row+1,14,\
                            "Support_"+num+" Horizontal Ball Support",\
                            "Location_x")
            self.supports[int(num)-1] = ("self.support_"+num+\
                          "= x_ball_support(float(self.support_"+num+\
                          "_Location_x(self)))")
            
    def create_x_roller_support(self):
        row = self.first_empty_row(4,9,5)
        if row != None:
            num = str(int(row/2-1))
            self.create_sub("support"+"_"+num,row,5,row+1,14,\
                            "Support_"+num+" Horizontal Roller Support",\
                            "Location_x")
            self.supports[int(num)-1] = ("self.support_"+num+\
                          "= x_roller_support(float(self.support_"+num+\
                          "_Location_x(self)))")
            
    def create_pin_support(self):
        row = self.first_empty_row(4,9,5)
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
        self.create_sub("remove",10,5,11,14)
        self.loads[0] = None
        
    def remove_load_2(self):
        self.create_sub("remove",12,5,13,14)
        self.loads[1] = None
        
    def remove_load_3(self):
        self.create_sub("remove",14,5,15,14)
        self.loads[2] = None
        
    def remove_load_4(self):
        self.create_sub("remove",16,5,17,14)
        self.loads[3] = None
        
    def remove_load_5(self):
        self.create_sub("remove",18,5,19,14)
        self.loads[0] = None
        
    def create_point_load(self):
        row = self.first_empty_row(10,19,5)
        if row != None:
            num = str(int(row/2-4))
            self.create_sub("load_"+num,row,5,row+1,14,"Load_"+num+\
                            " Point Load","Location_x","Load")
            self.loads[int(num)-1] = ("self.load_"+num+\
                       " = point_load(float(self.load_"+num+\
                       "_Location_x(self)),float(self.load_"+num+\
                       "_Load(self)))")
            
    def create_distributed_load(self):
        row = self.first_empty_row(10,19,5)
        if row != None:
            num = str(int(row/2-4))
            self.create_sub("load_"+num,row,5,row+1,14,"Load_"+num+\
                            " Distributed Load","Start_x","End_x","Load")
            self.loads[int(num)-1] = ("self.load_"+num+\
                       " = distributed_load(float(self.load_"+num+\
                       "_Start_x(self)),float(self.load_"+num+\
                       "_End_x(self)),float(self.load_"+num+\
                       "_Load(self)))")
            
    def create_linear_load(self):
        row = self.first_empty_row(10,19,5)
        if row != None:
            num = str(int(row/2-4))
            self.create_sub("load_"+num,row,5,row+1,14,"Load_"+num+\
                            " Linear Load","Start_x","End_x","Load_start","Load_end")
            self.loads[int(num)-1] = ("self.load_"+num+\
                       " = linear_load(float(self.load_"+num+\
                       "_Start_x(self)),float(self.load_"+num+\
                       "_End_x(self)),float(self.load_"+num+\
                       "_Load_start(self)),float(self.load_"+num+\
                       "_Load_end(self)))")
            
    def create_parabolic_load(self):
        row = self.first_empty_row(10,19,5)
        if row != None:
            num = str(int(row/2-4))
            self.create_sub("load_"+num,row,5,row+1,14,"Load_"+num+\
                            " Parabolic Load","Start_x","End_x","Load_start","Load_end")
            self.loads[int(num)-1] = ("self.load_"+num+\
                       " = parabolic_load(float(self.load_"+num+\
                       "_Start_x(self)),float(self.load_"+num+\
                       "_End_x(self)),float(self.load_"+num+\
                       "_Load_start(self)),float(self.load_"+num+\
                       "_Load_end(self)))")
            

class MatplotlibCanvas(FigureCanvas):
    
    def __init__(self):
        self.fig = Figure()
        self.axes = self.fig.add_subplot(111)
        self.axes.clear()
        units = self.get_units()
        self.axes.set_xlabel("x ("+units+")")
        self.axes.set_ylabel("y ("+units+")")
        self.axes.set_title("")
        
        FigureCanvas.__init__(self,self.fig)
        FigureCanvas.setSizePolicy(self,
                                   QSizePolicy.Expanding,
                                   QSizePolicy.Expanding)
        #FigureCanvas.setMaximumHeight(400)
        FigureCanvas.updateGeometry(self)
    
    def redraw(self):
        self.axes.clear()
        x = np.arange(1,10,1)
        y = np.arange(1,10,1)
        self.axes.plot(x,y)
        self.draw()
        
    def get_units(self):
        return "IPS"
    
    def get_density(self):
        return "lbf/in^3"
        


        
if __name__ == "__main__":
    app = QApplication(sys.argv)
    ex = main_GUI()
    app.exec_()
    
#FigureCanvas.printJpeg
