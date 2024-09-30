"""
Created on Thu Oct 13 16:19:41 2022
@author: Orlando Arroyo y Daniela Novoa
"""
"""----------------------------------------------------------------------------
                Generador de pórticos en concreto reforzado 3D
                             ETABS 17 to OpenSees
----------------------------------------------------------------------------"""

#%% Import libraries ---->
import pandas as pd
import numpy as np
from openseespy.opensees import * 
import opsvis as opsv
import matplotlib.pyplot as plt
import opseestools.analisis3D as an
import vfo.vfo as vfo
import os as os
import time
import opseestools.utilidades as ut

# Libreria temporal
import Etabs2Op_Library as lb

#%% Parámetros de entrada ---->
ruta = os.path.abspath(os.getcwd()) # Ruta del directorio actual (En donde se encuentra almancenado el archivo .py)
directory_path = os.path.join(os.path.join(ruta, "00_Excel_Results")) # Ruta de la carpeta de los Excels
Excel_ETABS = 'Excel_ETABS2OpenSees.xlsx'
# Excel consolidado para leer información del modelo ETABS
lb.consolidate_excel(directory_path,Excel_ETABS) 

#%% Nombre de variables - Dataframes ---->
#   df_Joints: Información de las coordenadas x, y, z globales de los joints
#   df_Frames: Información Dimensiones y material de cada sección (Columnas y Vigas)
#   df_MatProp: Información propiedades de los materiales definidos en el modelo

#%% Creación del modelo ---->

wipe()
model('basic', '-ndm', 3, '-ndf', 6)

npint = 5 # Numero de puntos de integración para las secciones

#%% Lectura y creación de los nodos ---->

# Exportar tabla "Objects and Elements - Joints" --> Información de coordenadas x, y, z de cada nodo
df_Joints = pd.read_excel(Excel_ETABS,sheet_name='Objects and Elements - Joints')

# Los Tags de los nodos generados (nodes_label) corresponden al "Element Label" del Joint definido por ETABS

xcoord,ycoord,zcoord,Nodes_Label = lb.genNodes3D(df_Joints)

print("|-----------------------------------------------------|")
print("|------------------ Nodes generated ------------------|")

# opsv.plot_model(node_labels=0)
# plt.show()

#%% Definición de materiales ---->

# Exportar tabla "Material Properties - Concrete" y "Frame Sections"
df_Frames = pd.read_excel(Excel_ETABS,sheet_name='Frame Sections')
df_MatProp = pd.read_excel(Excel_ETABS,sheet_name='Material Properties - Concrete')
df_MatProp.rename(columns={'Name':'Material'}, inplace=True) # Cambiar nombre de la columna

# Incluir en df_Frames el valor de f'c [kPa]
df_Frames= pd.merge(df_Frames, df_MatProp[['Material','Fc']],on='Material',how='inner') # Dataframe con el valor de FC para cada sección

# Asignar material a las secciones
UnConf_Label,Conf_Label,Steel_Label = lb.AssMaterials(df_Frames)

# Incluir al dataframe los Tags de los materiales
df_Frames['Steel_Label'] = Steel_Label
df_Frames['UnConf_Label'] = UnConf_Label
df_Frames['Conf_Label'] = Conf_Label

print("|-----------------------------------------------------|")
print("|----------------- Materials created -----------------|")

#%% Definicion de la losa ---->

# Importar "Shell Sections - Slab"
df_SecShell = pd.read_excel(Excel_ETABS,sheet_name='Shell Sections - Slab')

# Incluir en el dataframe el falor de f'c
df_SecShell = pd.merge(df_SecShell, df_MatProp[['Material','Fc']],on='Material',how='inner') # Dataframe con el valor de FC para cada sección
df_SecShell.rename(columns={'Name':'Section'}, inplace=True)

# Importar "Objets and Elements - Shells"
df_OEShell = pd.read_excel(Excel_ETABS,sheet_name='Objects and Elements - Shells')

# Importar "Shell Assignments - Sections"
df_AsgShell = pd.read_excel(Excel_ETABS,sheet_name='Shell Assignments - Sections')
df_AsgShell.rename(columns={'Unique Name':'Element Label'}, inplace=True) # Cambiar nombre de la columna

# Mezclar los dataframes
df_Shell = pd.merge(df_OEShell, df_AsgShell[['Element Label','Section']],on='Element Label',how='inner') 
df_Shell = pd.merge(df_Shell, df_SecShell[['Section','Slab Thickness','Fc']],on='Section',how='inner') 


slabtags = lb.AssSlabs(df_Shell)
df_Shell['slabtag'] = slabtags

print("|-----------------------------------------------------|")
print("|----------------- Shells generated ------------------|")

#%% Información de las secciones y transformaciones ---->

df_Frames.rename(columns={'Name':'Frame Property'}, inplace=True)

bsect = df_Frames['t2'].to_numpy() # Base de las secciones
hsect = df_Frames['t3'].to_numpy() # Altura de las secciones

coltrans = 1000000
vigXtrans = 2000000
vigYtrans = 3000000

geomTransf('PDelta',coltrans,*[0, -1, 0])
geomTransf('Linear',vigXtrans,*[0, -1, 0])
geomTransf('Linear',vigYtrans,*[1, 0, 0])

print("|-----------------------------------------------------|")
print("|-------------- Transformations created --------------|")

#%% Información del refuerzo de las secciones ---->

# Al dataframe con la información del refuerzo hay que agregarle la información de las dimensiones
df_Col = pd.read_excel(Excel_ETABS,sheet_name='Concrete Column Rebar Data') # lee las coordenadas de los muros
df_Col = pd.merge(df_Col, df_Frames[['Frame Property','t3','t2','Steel_Label','UnConf_Label','Conf_Label']],on='Frame Property',how='inner')

df_Beam = pd.read_excel(Excel_ETABS,sheet_name='Concrete Beam Rebar Data')
df_Beam = pd.merge(df_Beam, df_Frames[['Frame Property','t3','t2','Steel_Label','UnConf_Label','Conf_Label']],on='Frame Property',how='inner')

#%% Creación de las secciones de fibras ---->

colsectags, beamsectags = lb.fibSection(df_Col,npint,df_Beam)

df_Col['eletag'] = colsectags
df_Beam['eletag'] = beamsectags

df_Cons = pd.concat([df_Col[['Frame Property', 't3', 't2', 'eletag']],df_Beam[['Frame Property', 't3', 't2', 'eletag']]])

print("|-----------------------------------------------------|")
print("|-------------- Fiber sections created ---------------|")

#%% Consolidación de Información de los elementos ---->

# Joints de los elementos (Vigas y Columnas)
df_Elem = pd.read_excel(Excel_ETABS,sheet_name='Objects and Elements - Frames')

# Secciones de los elementos
df_Sections = pd.read_excel(Excel_ETABS,sheet_name='Frame Assignments - Sections')

df_Elem.rename(columns={'Element Label':'Unique Name'}, inplace=True) # Cambiarle el nombre para hacer merge
df_Sections = pd.merge(df_Sections, df_Elem[['Unique Name','Joint I','Joint J']],on='Unique Name',how='inner') 

df_Sections.rename(columns={'Analysis Section':'Frame Property'}, inplace=True)

# Dataframe información de elementos columnas y vigas
df_Elements = pd.merge(df_Sections, df_Cons[['Frame Property','t3','t2','eletag']],on='Frame Property',how='inner') 

#%% Creación de los elementos ---->

lb.genElements(df_Elements,coltrans,Nodes_Label,xcoord,ycoord,zcoord)
print("|-----------------------------------------------------|")
print("|----------------- Elements generated ----------------|")

#%% Para aplicar las cargas verticales ---->

df_Elements.rename(columns={'Unique Name':'Element Label'}, inplace=True)
df_Elm_Beams = df_Elements[df_Elements['Design Type'] == 'Beam']

# Cargar Tabla "Tributary Area and LLRF"
df_Tributary = pd.read_excel(Excel_ETABS, sheet_name='Tributary Area and LLRF')
df_Tributary.rename(columns={'Unique Name':'Element Label'}, inplace=True)
df_Tributary = df_Tributary[df_Tributary['Design Type'] == 'Beam'].drop(columns=['Story', 'Label', 'Design Type', 'LLRF'])

# Cargar y procesar datos de cargas - Tabla "Shell Loads - Uniform
df_Shell_Loads1 = pd.read_excel(Excel_ETABS, sheet_name='Shell Loads - Uniform')
df_Shell_Loads = df_Shell_Loads1[df_Shell_Loads1['Load Pattern'] == 'CV'].rename(columns={'Load': 'Load CV'})

df_Shell2_Loads = df_Shell_Loads1[df_Shell_Loads1['Load Pattern'] == 'CMsobreimpuesta'].drop(columns=['Story', 'Label', 'Load Pattern', 'Direction']).rename(columns={'Load': 'Load CM'}).reset_index(drop=True)
df_Shell_Loads = pd.merge(df_Shell_Loads, df_Shell2_Loads[['Unique Name', 'Load CM']], on='Unique Name', how='inner')

df_Elements_Beams = lb.AppBeamLoads(df_OEShell,df_Elm_Beams,df_Shell_Loads,df_Joints,df_Tributary)

timeSeries('Linear', 1)
pattern('Plain',1,1)

nels = len(df_Elements_Beams)

for ind3 in range(nels):
    aa = df_Elements_Beams.iloc[ind3]
    eleLoad('-ele',int(aa['Element Label']),'-type','-beamUniform',-aa['Force at Start'],0.0)
    
print("|-----------------------------------------------------|")
print("|------------------ Cargas aplicadas -----------------|")
#%% Para mostrarlo con extrusión ---->
# ele_shapes = {}
# for index, row in df_Elements.iterrows():
#     ele_shapes[int(row['Unique Name'])] = ['rect', [row['t2'], row['t3']]]
# # este ciclo extrae los nodos que conforman cada una de las losas del modelo y los deja en nodosM
# opsv.plot_extruded_shapes_3d(ele_shapes)

vfo.plot_model(show_nodetags='yes',show_eletags='yes')

#%% Asignación de diafragmas y masas ---->

df_Diaf = pd.read_excel(Excel_ETABS,sheet_name='Mass Summary by Diaphragm')

alturas = np.sort(pd.unique(df_Joints['Global Z']))
altur = alturas[1::]

Xcent = np.flipud(df_Diaf['X Mass Center'])
Ycent = np.flipud(df_Diaf['Y Mass Center'])

# ****Masa debe estar en toneladas****
Masa = np.flipud(df_Diaf['Mass X'])/1000

Inercia = np.flipud(df_Diaf['Mass Moment of Inertia'])

index1 = df_Joints['Global Z'] == alturas[0]
dia1 = df_Joints[index1]['Element Label'].astype(int).to_list()

for ind,alt in enumerate(altur):   
    index1 = df_Joints['Global Z'] == alt
    dia1 = df_Joints[index1]['Element Label'].astype(int).to_list()
    # Crea un nodo en el centro de masa
    node(100000000*ind,float(Xcent[ind]),float(Ycent[ind]),float(alt))
    
    fix(100000000*ind,0,0,1,1,1,0)
    
    mass(100000000*ind,float(Masa[ind]),float(Masa[ind]),float(Masa[ind]),0.0,0.0,float(Inercia[ind]))
    rigidDiaphragm(3,100000000*ind,*dia1)

print("|-----------------------------------------------------|")
print("|----------- Diaphragm and masses assigned -----------|")
print("|-----------------------------------------------------|")

#%%
eig = eigen(len(altur))

T1 = 2*3.1416/np.sqrt(eig[0])
T2 = 2*3.1416/np.sqrt(eig[1])

# modalProperties('-file', 'modal.txt', '-unorm', '-return')
modalProperties('-print', '-unorm')
# vfo.plot_modeshape(scale=50, contour='X', modenumber=3)

an.gravedad()
