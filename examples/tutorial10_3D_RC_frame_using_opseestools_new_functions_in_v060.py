# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 14:49:27 2024

@author: HOME
"""
#%% Cargado librerías

from openseespy.opensees import *
import matplotlib.pyplot as plt
import opsvis as opsv
import utilidades as ut
import opseestools.analisis3D as an
import opseestools.utilidades as ut2
import vfo.vfo as vfo
import numpy as np
import pandas as pd
#%% Creación del modelo

wipe()
model('basic','-ndm',3,'-ndf',6)


#%% Definicion de nodos
coordx = [0, 6, 12]
coordy = [0, 6, 12]
coordz = [0, 3, 6, 9]

masas = [90,90,90]

coords = ut.creategrid3D(coordx,coordy,coordz,1,masas)
fixZ(0.0,1,1,1,1,1,1)

#%% Definición de materiales (son según norma colombiana NSR-10, similares detallamiento especial del ACI)
fc = 28
fy = 420
noconf, conf, acero = ut2.col_materials(fc,fy)

#%% Definición de elementos

Bcol = 0.35 # base de la columna
Hcol = 0.35 # altura de la columna
Bvig = 0.25 # base de la viga
Hvig = 0.40 # altura de la viga

c = 0.05  # recubrimiento de las secciones

As4 = 0.000127 # area barra #4
As5 = 0.0002 # area barra #5
As6 = 0.000286
As7 = 0.000387 # area barra #7

col30x30 = 101 # tag de la columna
vig30x40 = 201 # tag de la viga

ut.create_rect_RC_section(col30x30, Hcol, Bcol, c, conf, noconf, acero, 3, As5, 3, As5, 4, As5)
ut.create_rect_RC_section(vig30x40, Hvig, Bvig, c, conf, noconf, acero, 3, As4, 3, As4)

#%% Creando los elementos
coltags = [col30x30,col30x30, col30x30] # one section tag per floor
# cols, vigx, vigy = ut.create_elements3D(coordx, coordy, coordz, col30x30, vig30x40, vig30x40)
cols, vigx, vigy, sectag_col = ut.create_elements3D(coordx, coordy, coordz, coltags, vig30x40, vig30x40)
opsv.plot_model(node_labels=0,gauss_points=False)
plt.show()

#%% Creando la losa
hslab = 0.10 # altura de la losa
Eslab = 1000*4400*(28)**0.5 # módulo de elasticidad de la losa
pois = 0.3 # relación de Poisson de la losa


ut.create_slabs(coordx, coordy, coordz, hslab, Eslab, pois)
vfo.plot_model(show_nodetags='yes', show_eletags='yes', show_nodes='yes')

#%% Cargando las vigas
floorx = -20  
floory = -20
roofx = -10 
roofy = -10
ut.load_beams3D(-20, -10, -20, -10, vigx, vigy, coordx, coordy)

#%% calculando modos

eig = eigen('-fullGenLapack',(len(coordz)-1)*2)
modalProperties('-print', '-unorm')

# vfo.plot_modeshape(scale=15, contour='X', modenumber=1)

#%% Analizando el modelo

an.gravedad()
plt.show()
loadConst('-time',0.0)

#%% Pushover direccion X

ut.pushover_loads3D(coordz)
dtecho,Vbasal = an.pushover2(coordz[-1]*0.03, 0.001, len(coordz)-1, 1,norm=[coordz[-1],np.sum(masas)*9,81])
plt.show()

#%% Pushover direccion Y

# ut.pushover_loads3D(coordz,pushdir='y')
# dtecho,Vbasal = an.pushover2(coordz[-1]*0.03, 0.001, len(coordz)-1, 2)
# plt.show()

#%%

# plt.plot(dtechoNL,VbasalNL,dtecho,Vbasal)
# plt.xlim([0,0.2])
# plt.ylim([0,300])
# plt.legend(['sin losa','con losa'])
# plt.xlabel('Desp. techo (m)')
# plt.ylabel('Corte basal (kN)')
# plt.show()

#%% 

# plt.plot(dtecho,Vbasal,dtecho2,Vbasal2)
# plt.legen(['12cm', '10cm'])