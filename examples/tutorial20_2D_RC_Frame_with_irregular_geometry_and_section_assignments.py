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
import itertools

#%%

wipe()
model('basic','-ndm',2,'-ndf',3)

#%% Nodos y apoyos

#user can enter coordinates

coordx = [0, 5.0, 12.0, 18.0]
coordy = [3*i for i in range(8)]
ut.creategrid(coordx,coordy)
fixY(0,1,1,1) 
opsv.plot_model()
#%% Materiales

fc = 28 # en MPa
fy = 420 # en MPa
tag_noconf, tag_conf, tag_acero = ut.col_materials(fc,fy,'DMO',nps=3)

#%% Secciones

Bcol = 0.55 # base de la columna
Hcol = 0.55 # altura de la columna
Bcol2 = 0.6
Hcol2 = 0.6
Bvig = 0.45 # base de la viga
Hvig = 0.5 # altura de la viga
Bvig2 = 0.45
Hvig2 = 0.6
c = 0.05  # recubrimiento de las secciones
As4 = 0.000127 # area barra #4
As5 = 0.0002 # area barra #5
As6 = 0.00028
As7 = 0.00387
As8 = 0.000508
col55x55 = 101 # tag de la columna
col60x60 = 102 # tag de la columna
vig45x50 = 201 # tag de la viga
vig45x60 = 202 # tag de la viga
ut.create_rect_RC_section(col55x55, Hcol, Bcol, c, tag_conf, tag_noconf, tag_acero , 5, As8, 5, As8, 6, As8)
ut.create_rect_RC_section(vig45x50, Hvig, Bvig, c, tag_conf, tag_noconf, tag_acero , 9, As6, 5, As4)
ut.create_rect_RC_section(col60x60, Hcol2, Bcol2, c, tag_conf, tag_noconf, tag_acero , 5, As8, 5, As8, 6, As8)
ut.create_rect_RC_section(vig45x60, Hvig2, Bvig2, c, tag_conf, tag_noconf, tag_acero , 9, As6, 5, As4)

#%% Elementos
# floorwise column definition.In this case there are as many columns as elements in coordx and floors as elements in coordy


columns_floor1 = [col55x55,col60x60,col60x60,col55x55] # colums of the first floor
columns_floor2 = [col55x55,col60x60,col60x60,col55x55]
columns_floor3 = [col55x55,col60x60,col60x60,col55x55]
columns_floor4 = [col55x55,col55x55,col55x55,col55x55]
columns_floor5 = [col55x55,col55x55,col55x55,col55x55]
columns_floor6 = [col55x55,col55x55,col55x55,col55x55]
columns_floor7 = ['None',col55x55,col55x55,'None']
building_columns = [columns_floor1,columns_floor2,columns_floor3,
                    columns_floor4,columns_floor5,columns_floor6, columns_floor7]

# floorwise beam definition.In this case there are as many beams as elements-1 in coordx and floors as elements in coordy

beams_floor1 = [vig45x50,vig45x60,vig45x60] # colums of the first floor
beams_floor2 = [vig45x50,vig45x60,vig45x60]
beams_floor3 = [vig45x50,vig45x60,vig45x60]
beams_floor4 = [vig45x50,vig45x60,vig45x60]
beams_floor5 = [vig45x50,'None',vig45x50] # second span has no beam
beams_floor6 = [vig45x50,vig45x50,vig45x50]
beams_floor7 = ['None',vig45x50,'None'] # first and third span have no beam

building_beams = [beams_floor1,beams_floor2,beams_floor3,
                    beams_floor4,beams_floor5,beams_floor6,beams_floor7]

tagcols,tagbeams,column_info,beam_info = ut.create_elements2(coordx,coordy,building_columns,building_beams,output=1) 

ut.remove_hanging_nodes(tagcols,tagbeams)
opsv.plot_model()

# model_nodes_updated = np.array(getNodeTags())
# floor_number = 7
# index_np1 = model_nodes_updated%(1000)==floor_number
# nodes_floor = list(model_nodes_updated[index_np1])

floor_diaphragms = [1,1,1,1,0,1,0]
ut.apply_diaphragms(floor_diaphragms,output=1)

#%% Cargas

# Case a: same for all floor beams and same for all roof beams
# Case b: beamwise

load_type = 'beamwise' # options: same or beamwise
tag = 1
output = 1
if load_type == 'same':
    floor_beam_loads = 70
    roof_beam_loads = 50
    ut.load_beams(-floor_beam_loads,-roof_beam_loads,tagbeams)
elif load_type == 'beamwise':
    beams_loads_floor1 = [-70,-70,-70] # beams of the first floor
    beams_loads_floor2 = [-70,-70,-70]
    beams_loads_floor3 = [-70,-70,-70]
    beams_loads_floor4 = [-70,-70,-70]
    beams_loads_floor5 = [-70,-70] # because floor 5 only has two beams
    beams_loads_floor6 = [-70,-70,-70]
    beams_loads_floor7 = [-50]
    
    beam_loads = [beams_loads_floor1,beams_loads_floor2,beams_loads_floor3,
                  beams_loads_floor4,beams_loads_floor5,beams_loads_floor6,
                  beams_loads_floor7]
    beam_loads_flatten2 = list(itertools.chain.from_iterable(beam_loads))

    if len(beam_loads_flatten2) != len(tagbeams):
        print('number of beam loads and beams dont match') # stop here don't continue
    else:
        ut.load_beams2(beam_loads,tagbeams,output=1)
        

#%% Analisis gravedad y ver modelo

an.gravedad()
loadConst('-time',0.0)
# opsvis.set_plot_props(line_width=4, point_size=2, notebook=True)
# opsvis.plot_model()1
opsv.plot_model(node_labels=1,element_labels=1,gauss_points=False)

#%% Analisis pushover
leftmost_nodes = ut.find_leftmost_nodes(coordy)  
ut.pushover_loads(coordy,nodes=leftmost_nodes)
elements = tagcols+tagbeams
nodes_control = [1000]+[int(i) for i in leftmost_nodes]
# nodes_control
dtecho,Vbasal,drifts,rotations = an.pushover2DRot(0.05*coordy[-1],0.001,getNodeTags()[-1],1,nodes_control,elements)


#%% PostProcessing

story_drifts = drifts[:len(dtecho),:]
column_rotations = rotations[:len(tagcols),:len(dtecho),[1,2]]
beam_rotations = rotations[len(tagcols):,:len(dtecho),[1,2]]



