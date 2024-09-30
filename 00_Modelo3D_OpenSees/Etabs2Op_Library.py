# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 14:48:10 2024

@author: alein
"""

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

"""----------------------------------------------------------------------------
                                Consolidate Excel
----------------------------------------------------------------------------"""
def consolidate_excel(directory_path,output_file):
    
    """
    Reads Excel files located in a directory and verifies the sheets names.
    It consolidates the information into a single Excel file, assigning the data to specific sheets.
    The sheet names correspond to those assigned by ETABS version 17.
    
    Parameters
    ----------
    directory_path : String
        Directory path containing the Excel files.
    output_file : String
        Consolidate Excel file name.
   
    Returns
    -------
    Dataframe
        Dataframe to consolidate Excel file.
    """
    # Expected sheets names directory
    expected_sheets = {'Concrete Beam Rebar Data':'Concrete Beam Rebar Data',
                       'Concrete Column Rebar Data':'Concrete Column Rebar Data',
                       'Frame Assignments - Sections':'Frame Assignments - Sections',
                       'Frame Sections':'Frame Sections',
                       'Mass Summary by Diaphragm':'Mass Summary by Diaphragm',
                       'Material Properties - Concrete':'Material Properties - Concrete',
                       'Modal Participating Mass Ratios':'Modal Participating Mass Ratios',
                       'Objects and Elements - Frames':'Objects and Elements - Frames',
                       'Objects and Elements - Joints':'Objects and Elements - Joints',
                       'Objects and Elements - Shells':'Objects and Elements - Shells',
                       'Shell Assignments - Sections':'Shell Assignments - Sections',
                       'Shell Loads - Uniform':'Shell Loads - Uniform',
                       'Shell Sections - Slab':'Shell Sections - Slab',
                       'Tributary Area and LLRF':'Tributary Area and LLRF'}
    file_index = 0
    # Create a ExcelWriter file
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        file_paths = [os.path.join(directory_path,file) for file in os.listdir(directory_path) if file.endswith('.xlsx')]
        for file_path in file_paths:
            excel_file = pd.ExcelFile(file_path) 
            for sheet_name in excel_file.sheet_names:
                for var_name, expected_sheet_name in expected_sheets.items():
                    if sheet_name == expected_sheet_name:
                        file_index = file_index+1
                        df = pd.read_excel(file_path,sheet_name=sheet_name,skiprows=1).drop(index=0)
                        df.to_excel(writer,sheet_name=expected_sheet_name,index=False)
                        
        print("|-----------------------------------------------------|")
        print("|-------------- CONSOLIDATION COMPLETE ---------------|")
        print(f"|-------------- {file_index} sheets were assigned --------------|")
        print(f"|----- File saved as '{output_file}' -----|")
        
    

"""----------------------------------------------------------------------------
                                Generate 3D Nodes
----------------------------------------------------------------------------"""
def genNodes3D(df):
    
    """
    Read the x, y, and z coordinates of each joint modeled in ETABS and generate nodes based on that information.
    Nodes at the base are modeled by restraining all degrees of freedom (fixed)
    
    Parameters
    ----------
    df : DataFrame
        DataFrame with the Joint Coordinate Information
    
    Returns
    -------
    xcoord : Array
        X Coordinate of each Joint
    ycoord : Array
        Y Coordinate of each Joint
    zcoord : Array
        Z Coordinate of each Joint
    nlabel : Array
        Joints Tag
    """
    
    xcoord = df['Global X'].to_numpy() # x coordinate
    ycoord = df['Global Y'].to_numpy() # y coordinate
    zcoord = df['Global Z'].to_numpy() # z coordinate
    nlabel = df['Element Label'].to_numpy(dtype=int) # Joint Tag
    
    # Cycle to create the nodes of the model
    nnodes = len(xcoord)
    for i in range(nnodes):
        node(int(nlabel[i]),float(xcoord[i]),float(ycoord[i]),float(zcoord[i]))
    
    # Restrictions at the base
    fixZ(0.0,1,1,1,1,1,1)   # Fixed base nodes
    
    return xcoord,ycoord,zcoord,nlabel

"""----------------------------------------------------------------------------
                     Asignar materiales a las secciones
----------------------------------------------------------------------------"""
def AssMaterials(df_Frames):
    
    Unconf_Tag, Conf_Tag, Steel_Tag = [], [], []
    
    for index in range(len(df_Frames)):
        
        aa = df_Frames.iloc[index]
        fc = (aa['Fc'])/1000 # El f'c debe estar dado en KPa desde el modelo
        
        # El Tag del material viene asociado a la cantidad de secciones (Frames)
        # que el modelo tiene porgramados, aunque no los esté utilizando el modelo
        # (Puede que el diseñador sin necesidad haya creado muchos y no los haya eliminado)
        
        # Tag de los materiales:
        #   Tag del material ACERO --> comienza en 100
        #   Tag del material CONCRETO NO CONFINADO --> comienza en 200
        #   Tag del material CONCRETO CONFINADO --> comienza en 300
        
        # Comienza en 1, termina en número de frames
        Tag_Steel = int(index)
        Tag_UnConf = int(len(df_Frames)+(index))
        Tag_Conf = int(len(df_Frames)*2+(index))
            
        unctag,conftag,steeltag = col_materials(Tag_Steel,Tag_UnConf,Tag_Conf,fc,420,'DMO','tension')
        
        Unconf_Tag.append(unctag)
        Conf_Tag.append(conftag)
        Steel_Tag.append(steeltag)
    
    return Unconf_Tag,Conf_Tag,Steel_Tag

"""----------------------------------------------------------------------------
                     Asignar materiales a las losas
                             Generar losas
----------------------------------------------------------------------------"""

def AssSlabs(df_Shell):
    slabtags = []
    for index in range(len(df_Shell)):
        aa = df_Shell.iloc[index]
        fc = aa['Fc']
        hslab = aa['Slab Thickness']
        Eslab = 1000*4400*(fc/1000)**0.5
        seclosa = -(index+1)
        # Se crea la seccion de la losa
        section('ElasticMembranePlateSection', seclosa, 1.0*Eslab, 0.0, hslab, 0.0)
        
        tag = -(1000+index)
        
        n1 = int(aa['Joint 1'])
        n2 = int(aa['Joint 2'])
        n3 = int(aa['Joint 3'])
        n4 = int(aa['Joint 4'])
        nodoslosa = [n1,n2,n3,n4]
        
        element('ShellDKGQ', tag, *nodoslosa, seclosa)
        slabtags.append(tag)
    
    return slabtags
    

"""----------------------------------------------------------------------------
                            Crear secciones de fibras
----------------------------------------------------------------------------"""
def fibSection(df_Col,npint,df_Beam):
    
    colsectags = []
    for index in range(len(df_Col)):
        aa = df_Col.iloc[index]
        steel_area = aa['Corner Bar Area']/1000000 # --> Área en m^2
        
        id_tag = 1000 + index
        ut.BuildRCSection(id_tag, aa['t3'], aa['t2'], aa['Cover'], aa['Cover'], int(aa['Conf_Label']), 
                          int(aa['UnConf_Label']), int(aa['Steel_Label']), int(aa['# Long. Bars 2-axis']), 
                          steel_area, int(aa['# Long. Bars 2-axis']), steel_area, 
                          int(2*(aa['# Long. Bars 3-axis']-2)), steel_area, 12, 12, 8, 8)
        
        colsectags.append(id_tag)
        beamIntegration('Lobatto', id_tag, id_tag, npint)

    ncols = len(colsectags)
    beamsectags = []
    for ind1 in range(len(df_Beam)):
        tagbeam = 1000 + (ind1+ncols)
        aa = df_Beam.iloc[ind1]
        steel_top = aa['area top']/1000000
        steel_bot = aa['area bottom']/1000000
        
        ut.BuildRCSection(tagbeam, aa['t3'], aa['t2'], aa['Top Cover'], aa['Top Cover'], int(aa['Conf_Label']), 
                          int(aa['UnConf_Label']), int(aa['Steel_Label']), int(aa['#top']), steel_top, 
                          int(aa['#bottom']), steel_bot, 2, 1e-14, 12, 12, 8, 8)
        beamsectags.append(tagbeam)
        beamIntegration('Lobatto', tagbeam, tagbeam, npint)
        
    return colsectags,beamsectags

"""----------------------------------------------------------------------------
                     Generar elementos (Columnas y Vigas)
----------------------------------------------------------------------------"""

def genElements(df_Elements,coltrans,Nodes_Label,xcoord,ycoord,zcoord):
    
    nels = len(df_Elements)
    for ind2 in range(nels):
        # print(ind2)
        aa = df_Elements.iloc[ind2]
        nodoI = aa['Joint I']
        nodoJ = aa['Joint J']
        eleNodes = [int(nodoI), int(nodoJ)]
        if aa['Design Type'] == 'Column':
            transfTag = coltrans
        elif aa['Design Type'] == 'Beam':
            iA = np.where(Nodes_Label==nodoI)
            iB = np.where(Nodes_Label==nodoJ)
            A = np.array([xcoord[iA][0],ycoord[iA][0],zcoord[iA][0]])
            B = np.array([xcoord[iB][0],ycoord[iB][0],zcoord[iB][0]])
            AB = B-A
            zg = [0,0,1]
            cP = np.cross(AB,zg)
            vectrans = cP/np.linalg.norm(cP)
            geomTransf('Linear',int(aa['Unique Name']),*vectrans)
            transfTag = int(aa['Unique Name'])
            
        element('forceBeamColumn', int(aa['Unique Name']), *eleNodes, transfTag, int(aa['eletag']))

"""----------------------------------------------------------------------------
                     Aplicar carga distribuida a vigas
----------------------------------------------------------------------------"""

def AppBeamLoads(df_OEShell,df_Elm_Beams_Original,df_Shell_Loads,df_Joints,df_Tributary):
    
    df_Elm_Beams = df_Elm_Beams_Original.copy()
    # Procesar la matriz Shell2Beam
    Shell2Beam = np.array([
        [int(bb['Element Label']) if ((aa['Joint I'] == bb['Joint 1'] and aa['Joint J'] == bb['Joint 2']) or #Joint 1 y 2 
                                    (aa['Joint I'] == bb['Joint 2'] and aa['Joint J'] == bb['Joint 1']) or #Joint 1 y 2
                                    
                                    (aa['Joint I'] == bb['Joint 2'] and aa['Joint J'] == bb['Joint 3']) or #Joint 2 y 3
                                    (aa['Joint I'] == bb['Joint 3'] and aa['Joint J'] == bb['Joint 2']) or #Joint 2 y 3
                                    
                                    (aa['Joint I'] == bb['Joint 3'] and aa['Joint J'] == bb['Joint 4']) or #Joint 3 y 4
                                    (aa['Joint I'] == bb['Joint 4'] and aa['Joint J'] == bb['Joint 3']) or #Joint 3 y 4
                                    
                                    (aa['Joint I'] == bb['Joint 4'] and aa['Joint J'] == bb['Joint 1']) or #Joint 4 y 1
                                    (aa['Joint I'] == bb['Joint 1'] and aa['Joint J'] == bb['Joint 4'])) #Joint 4 y 1
          else 0
          for _, bb in df_OEShell.iterrows()]
        for _, aa in df_Elm_Beams.iterrows()
    ])

    # Extraer cargas asociadas a cada elemento
    Shell_2_Beam = [[num for num in row if num != 0] for row in Shell2Beam.tolist()]

    CVLoads = [np.mean([float(df_Shell_Loads[df_Shell_Loads['Unique Name'] == label]['Load CV']) for label in labels]) if labels else 0 for labels in Shell_2_Beam]
    CMLoads = [np.mean([float(df_Shell_Loads[df_Shell_Loads['Unique Name'] == label]['Load CM']) for label in labels]) if labels else 0 for labels in Shell_2_Beam]

    df_Elm_Beams['Load CM'] = CMLoads
    df_Elm_Beams['Load CV'] = CVLoads

    Length = []
    for index in range(len(df_Elm_Beams)):
        aa = df_Elm_Beams.iloc[index]
        JointI = int(aa['Joint I'])
        JointJ = int(aa['Joint J'])
        LengthI,LengthJ = [],[]
        for index2 in range(len(df_Joints)):
            bb = df_Joints.iloc[index2]
            Joint = int(bb['Element Label'])   
            if JointI == Joint:
                LengthI.append([bb['Global X'],bb['Global Y']])
            if JointJ == Joint:
                LengthJ.append([bb['Global X'],bb['Global Y']])
                
        Length.append(np.abs(LengthI[0][0]-LengthJ[0][0])+np.abs(LengthI[0][1]-LengthJ[0][1]))
        
    df_Elm_Beams['Length'] = Length

    # # Agregar y calcular la fuerza al DataFrame final
    df_Elm_Beams = df_Elm_Beams.merge(df_Tributary[['Element Label', 'Tributary Area']], on='Element Label', how='inner')

    df_Elm_Beams['Force at Start'] = (df_Elm_Beams['Tributary Area'] / df_Elm_Beams['Length']) * (1.0 * df_Elm_Beams['Load CM'] + 0.25 * df_Elm_Beams['Load CV'])

    return df_Elm_Beams

def col_materials(steeltag, unctag, conftag,fcn=28,fy=420,detailing='DES',tension = 'tension'):
    '''
    Generates materials for concrete and steel. The concrete has regularization applied
    and generates both, unconfined and confined concrete. Steel is defined based on the 
    Dhakal and Maekawa model accounting for buckling and low-cycle fatigue

    Parameters
    ----------
    zone : String, optional
        Detailing level. By default it is special detailing. User can input 'DMO' for an intermediate detailing of Colombian design code.
    fy : Float, optional
        Steel yielding stress in MPa. The default is 420.
    fcn : Float, optional
        Concrete compressive stress in MPa. The default is 28.
    tension : String, optional
        String to define if concrete has tension capacity. The default is 'tension'. Use no to ignore tension
    steeltag : Integer, optional
        Tag of the steel material. The default is 100.
    unctag : Integer, optional
        Tag of the unconfined concrete material. The default is 102.
    conftag : Integer, optional
        Tag of the confined concrete material. The default is 101.

    Returns
    -------
    list
        List with the material tags for the unconfined concrete, confined concrete and reinforcement steel.

    '''
    
    # Steel properties
    fy_1 = fy                                                              # fy del acero
    fu_1 = fy_1*1.4 # aprox considerando los valores reportados por Julian                                                       
    fh_1 = fy_1
    E1 = 200000
    ey_1 = fy_1/E1
    eh = 0.01
    eu = 0.1
    
    # concrete properties
    fc = fcn*1000
    Ec = 4400*np.sqrt(fcn)*1000
    ec = 2*fc/Ec
    fcu = 0.2*fc
    e20_c = e20Lobatto2(28, 3000, 5, 28, 24000, ec)
    ecu = e20_c
    
    if tension == 'tension':
        uniaxialMaterial('Concrete02', unctag, fc, ec, fcu, ecu)
    else:
        uniaxialMaterial('Concrete01', unctag, fc, ec, fcu, ecu)
    
    if detailing == 'DMO':
        D = 12
        L = 8*D
        s_steel1, e_steel1 = dhakal(fy_1, fu_1, ey_1, eh, eu, L, D)
        s1p1, s2p1, s3p1, s4p1, e1p1, e2p1, e3p1, e4p1 = s_steel1[0], s_steel1[1], s_steel1[2], s_steel1[3],e_steel1[0], e_steel1[1], e_steel1[2], e_steel1[3]
        s1n1, s2n1, s3n1, s4n1, e1n1, e2n1, e3n1, e4n1 = s_steel1[4], s_steel1[5], s_steel1[6], s_steel1[7],e_steel1[4], e_steel1[5], e_steel1[6], e_steel1[7]
        uniaxialMaterial('HystereticSM',steeltag,'-posEnv',s1p1,e1p1,s2p1,e2p1,s3p1,e3p1,s4p1,e4p1,'-negEnv',s1n1,e1n1,s2n1,e2n1,s3n1,e3n1,s4n1,e4n1)
        # Para el concreto confinado
        k=1.25
        fcc=fc*k
        Ec = 4400*np.sqrt(fcc/1000)*1000
        ecc= 2*fcc/Ec
        fucc=0.2*fcc
        e20_cc = e20Lobatto2(2*fc, 3000, 5, 28*k, Ec/1000, ecc)
        eucc=e20_cc
        if tension == 'tension':
            uniaxialMaterial('Concrete02', conftag, fcc, ecc, fucc, eucc)
        else:
            uniaxialMaterial('Concrete01', conftag, fcc, ecc, fucc, eucc) 
        
    elif detailing =='DES':
        D = 12
        L = 6*D
        s_steel1, e_steel1 = dhakal(fy_1, fu_1, ey_1, eh, eu, L, D)
        s1p1, s2p1, s3p1, s4p1, e1p1, e2p1, e3p1, e4p1 = s_steel1[0], s_steel1[1], s_steel1[2], s_steel1[3],e_steel1[0], e_steel1[1], e_steel1[2], e_steel1[3]
        s1n1, s2n1, s3n1, s4n1, e1n1, e2n1, e3n1, e4n1 = s_steel1[4], s_steel1[5], s_steel1[6], s_steel1[7],e_steel1[4], e_steel1[5], e_steel1[6], e_steel1[7]
        uniaxialMaterial('HystereticSM',steeltag,'-posEnv',s1p1,e1p1,s2p1,e2p1,s3p1,e3p1,s4p1,e4p1,'-negEnv',s1n1,e1n1,s2n1,e2n1,s3n1,e3n1,s4n1,e4n1)
        # Para el concreto confinado
        k=1.3
        fcc=fc*k
        Ec = 4400*np.sqrt(fcc/1000)*1000
        ecc= 2*fcc/Ec
        fucc=0.2*fcc
        e20_cc = e20Lobatto2(2*fc, 3000, 5, 28*k, Ec/1000, ecc)
        eucc=e20_cc
        if tension == 'tension':
            uniaxialMaterial('Concrete02', conftag, fcc, ecc, fucc, eucc)
        else:
            uniaxialMaterial('Concrete01', conftag, fcc, ecc, fucc, eucc)
            
    return unctag,conftag,steeltag


def e20Lobatto2(Gfc,Lel,npint,fc,E,e0):
    '''
    Calculates the ultimate strain for a concrete material applying regularization based on the constant fracture energy proposed by Coleman and Spacone

    Parameters
    ----------
    Gfc : float
        facture energy in N/mm.
    Lel : float
        element length.
    npint : int
        number of integration points.
    fc : float
        concrete compressive strength in MPa.
    E : float
        concrete modulus of elasticity strength in MPa..
    e0 : float
        strain associated to fc.

    Returns
    -------
    e20 : float
        ultimate strain corresponding to 0.2fc according to Coleman and Spacone.

    '''
    
    # TODO TIENE QUE ESTAR EN UNIDADES DE N y mm
    # Gfc entra en N/mm: Energía de fractura
    # Lel es la longitud del elemento en mm
    # npint es el número de puntos de integración
    # fc es el esfuerzo a compresión del concreto en N/mm2 (MPa)
    # E es el módulo de elasticidad del concreto en N/mm2 (MPa)
    # e0 es la deformación del concreto en fc
    
    if npint == 4:
        LIP = Lel/2*1/6
    elif npint == 5:
        LIP = Lel/2*1/10
    elif npint == 6:
        LIP = Lel/2*1/15
    elif npint == 3:
        LIP = Lel/2*1/3
    else:
        LIP = 0.1*Lel
        print('numero de puntos no soportado')
    
    e20 = Gfc/(0.6*fc*LIP)-0.8*fc/E+e0
    
    return e20

def dhakal(Fyy, Fuu, eyy, ehh, euu, Lb, Db):
    #Fy = esfuerzo de fluencia del acero [MPa]
    #Fu = esfuerzo último del acero [MPa]
    #ey, eh, eu = deformaciones unitarias del acero
    #L = Espaciamiento entre estribos (Longitud libre para pandearse) [mm]
    #D = Diámetro de la barra [mm]
    fy = Fyy
    fu = Fuu
    fh = fy+0.01
    ey = eyy
    eu = euu
    eh = ehh
    L = Lb
    D = Db
    alfa = 0.75
    efracture = 0.05
    espalling = 0.004
    sigma_u = 0.2*fy
    p1 = [eh,eu]
    p2 = [fh,fu]
    Es = fy/ey
    m = -0.02*Es
    eas = np.max([(55-2.3*np.sqrt(fy/100)*L/D)*ey,7*ey])
    sigma_l = np.interp(eas,p1,p2)
    sigma_as = np.max([alfa*(1.1-0.016*np.sqrt(fy/100)*L/D)*sigma_l,0.2*fy])
    eu_d = (sigma_u-sigma_as)/m + eas
    sigma_f = np.interp(efracture,p1,p2)
    sigma_s = np.interp(espalling,p1,p2)
    # strain = [-eu_d,-espalling,-ey,0.0,ey,efracture,eu]
    # stress = [-sigma_u,-sigma_s,-fy,0.0, fy,sigma_f, sigma_u]
    strain = [-eu_d,-eas,-espalling,-ey,0.0,ey,eh,efracture,eu]
    stress = [-sigma_u,-sigma_as,-sigma_s,-fy,0.0, fy,fh,sigma_f, sigma_u]

    s = [fy*1000,fh*1000,sigma_f*1000, sigma_u*1000,-fy*1000,-sigma_s*1000,-sigma_as*1000,-sigma_u*1000]
    e = [ey,eh,efracture,eu,-ey,-espalling,-eas,-eu_d]
 
    
    return s, e
