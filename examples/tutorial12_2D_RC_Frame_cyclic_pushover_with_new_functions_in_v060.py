#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 10:48:00 2024

@author: odarroyo
"""

from openseespy.opensees import *
import matplotlib.pyplot as plt
import opseestools.analisis as an
import opseestools.utilidades as ut
import opsvis as opsv
import numpy as np

#%%

wipe()
model('basic','-ndm',2,'-ndf',3)

#%% Nodos y apoyos

coordx = [0, 6.0, 10.0]
coordy = [0, 4.0, 7.0]
ut.creategrid(coordx,coordy)
fixY(0,1,1,1) 

#%% Materiales

fc = 28 # en MPa
fy = 420 # en MPa
tag_noconf, tag_conf, tag_acero = ut.col_materials(fc,fy)

#%% Secciones

Bcol = 0.3 # base de la columna
Hcol = 0.3 # altura de la columna
Bvig = 0.3 # base de la viga
Hvig = 0.4 # altura de la viga
c = 0.05  # recubrimiento de las secciones
nFibCover, nFibCore, nFibZcore  = 8, 16, 10 # numero de fibras en zona no confinada, en Z (x) en zona confinada y Y (y) zona confinada
As4 = 0.000127 # area barra #4
As5 = 0.0002 # area barra #5

col30x30 = 101 # tag de la columna
vig30x40 = 201 # tag de la viga

ut.create_rect_RC_section(col30x30, Hcol, Bcol, c, tag_conf, tag_noconf, tag_acero , 4, As5, 4, As5, 4, As5)
ut.create_rect_RC_section(vig30x40, Hvig, Bvig, c, tag_conf, tag_noconf, tag_acero , 3, As4, 4, As4)


#%% Elementos

tagcols,tagbeams = ut.create_elements(coordx,coordy,col30x30,vig30x40) 

#%% Cargas

ut.load_beams(-30,-20,tagbeams)

#%% Analisis gravedad y ver modelo

an.gravedad()
loadConst('-time',0.0)
opsv.plot_model()

#%% Analisis pushover

ut.pushover_loads(coordy)

displ = np.array([0.005,-0.005,0.01,-0.01,0.02,-0.02,0.025,-0.025,0.03,-0.03,0.035,-0.035,0.04,-0.04,0.045,-0.045])*coordy[-1]
disp = displ.tolist()
dtecho,Vbasal = an.pushover2C(disp,0.001,getNodeTags()[-1],1)









