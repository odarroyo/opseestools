# -*- coding: utf-8 -*-
"""
Created on Fri Sep  6 14:49:27 2024

@author: HOME
"""
#%% Cargado librerías

from openseespy.opensees import *
import matplotlib.pyplot as plt
import opsvis as opsv
import opseestools.analisis3D as an
import opseestools.utilidades as ut
import vfo.vfo as vfo
import numpy as np
import pandas as pd
import time
#%% Creación del modelo

wipe()
model('basic','-ndm',3,'-ndf',6)


#%% Definicion de nodos
coordx = [0, 8.6, 12.6]
coordy = [0, 7.3, 14.6, 21.9]
coordz = [0, 3, 6, 9, 12, 15]

masas = [210,210,210,210,210]

coords = ut.creategrid3D(coordx,coordy,coordz,1,masas)
fixZ(0.0,1,1,1,1,1,1)

#%% Definición de materiales (son según norma colombiana NSR-10, similares detallamiento especial del ACI)
fc = 28
fy = 420
noconf, conf, acero = ut.col_materials(fc,fy,nps=4,tension='yes')

#%% Definición de elementos

Bcol = 0.40 # base de la columna
Hcol = 1.00 # altura de la columna
Bvig = 0.40 # base de la viga
Hvig = 0.50 # altura de la viga
Bvig2 = 0.40 # base de la viga
Hvig2 = 0.50 # altura de la viga

c = 0.05  # recubrimiento de las secciones

As4 = 0.000127 # area barra #4
As5 = 0.0002 # area barra #5
As6 = 0.000286
As7 = 0.000387 # area barra #7
As8 = 0.000508 # area barra #8

col30x30 = 101 # tag de la columna
vig30x40 = 201 # tag de la viga
vig30x40_2 = 202 # tag de la segunda viga

ut.create_rect_RC_section(col30x30, Hcol, Bcol, c, conf, noconf, acero, 4, As8, 4, As8, 6, As8)
ut.create_rect_RC_section(vig30x40, Hvig, Bvig, c, conf, noconf, acero, 3, As5, 6, As5)
ut.create_rect_RC_section(vig30x40_2, Hvig2, Bvig2, c, conf, noconf, acero, 3, As5, 3, As5)

#%% Creando los elementos
eje_y = [col30x30]*4
p1 = [eje_y]*3
coltags = [p1]*5 # one section tag per floor

vigax1 = [vig30x40]*4 # son cuatro ejes de vigas que van en la dirección X
vigap1 = [vigax1]*2 # hay dos luces en la dirección X.
vigased = [vigap1]*5

vigay1 = [vig30x40_2]*3 # son cuatro ejes de vigas que van en la dirección X
vigap1y = [vigay1]*3 # hay dos luces en la dirección X.
vigasedY = [vigap1y]*5 

# cols, vigx, vigy = ut.create_elements3D(coordx, coordy, coordz, col30x30, vig30x40, vig30x40_2)
cols, vigx, vigy, sectag_col, sectag_vigx, sectag_vigy = ut.create_elements3D2(coordx, coordy, coordz, coltags, vigased, vigasedY)
opsv.plot_model(node_labels=1,gauss_points=False)
plt.show()

#%% Creando la losa
hslab = 0.15 # altura de la losa
Eslab = 1000*4400*(28)**0.5 # módulo de elasticidad de la losa
pois = 0.3 # relación de Poisson de la losa


ut.create_slabs(coordx, coordy, coordz, hslab, Eslab, pois)
vfo.plot_model(show_nodetags='yes', show_eletags='yes', show_nodes='yes')

#%% Cargando las vigas
floorx = -30
floory = -30
roofx = -20
roofy = -20
ut.load_beams3D(floorx, floory, roofx, roofy, vigx, vigy, coordx, coordy)

#%% calculando modos

eig = eigen('-fullGenLapack',(len(coordz)-1)*2)
modalProperties('-print', '-unorm')

vfo.plot_modeshape(scale=15, contour='X', modenumber=1)

#%% Analizando el modelo
recorder('Element','-file','losa.txt','-ele',taglosa,'forces')
recorder('Element','-file','losa.txt','-ele',taglosa,'stresses')

an.gravedad()
plt.show()
loadConst('-time',0.0)

#%% Pushover direccion X
elementos = cols + vigx
# fuerzas = [eleForce(el)[2] for el in cols]
ut.pushover_loads3D(coordz)
stime = time.time()
dtecho,Vbasal = an.pushover2(coordz[-1]*0.05, 0.0025, len(coordz)-1, 1,norm=[coordz[-1],np.sum(masas)*9,81],Tol=1e-2)
etime = time.time()
ttotal = etime - stime
print('tiempo de ejecucion: ',ttotal,'segundos')
plt.show()

#%% Pushover direccion Y
# elementos = cols + vigy
# ut.pushover_loads3D(coordz,pushdir='y')
# # dtecho,Vbasal = an.pushover2(coordz[-1]*0.03, 0.001, len(coordz)-1, 2, norm=[coordz[-1],np.sum(masas)*9,81])
# # dtecho,Vbasal,fuerzas,rotaciones = an.pushover2Rot(coordz[-1]*0.03, 0.001, len(coordz)-1, 2, elements=elementos, norm=[coordz[-1],np.sum(masas)*9,81])

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