# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 12:50:00 2022

@author: orlandoaram
"""
from openseespy.opensees import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import argrelextrema
from scipy.stats import gmean
from scipy.fft import fft, ifft
from scipy.integrate import cumulative_trapezoid
import pandas as pd

def MomentCurvature(secTag, axialLoad, maxK, numIncr=300):
    '''
    
    The original script is available in the OpenSees Wiki. Cumputes the moment curvature of a section
    
    Parameters
    ----------
    secTag : int
        tag of the section.
    axialLoad : float
        applied axial load.
    maxK : float
        max curvature.
    numIncr : int, optional
        number of steps for the calculation. The default is 300.

    Returns
    -------
    M : list
        moments of the section.
    curv : list
        curvature.

    '''
    # Script tomado de la librería de OpenSeespy de la web
    # secTag es el tag de la sección
    # axialLoad es la carga axial de la sección
    # maxK es la curvatura
    # numIncr es el número de incrementos
    # Define two nodes at (0,0)
    
    model('basic','-ndm',2,'-ndf',3)
    a = getNodeTags()
    if not a:
        a = [0] # si no hay nodos inicia la lista en cero
    n1 = max(a)+1
    n2 = max(a)+2
    print(n1,n2)
    node(n1, 0.0, 0.0)
    node(n2, 0.0, 0.0)

    # Fix all degrees of freedom except axial and bending
    fix(n1, 1, 1, 1)
    fix(n2, 0, 1, 0)
    
    # Define element
    #                             tag ndI ndJ  secTag
    b = getEleTags()
    if not b:
        b = [0] # si no hay nodos inicia la lista en cero
    e1 = max(b)+1
    element('zeroLengthSection',  e1,   n1,   n2,  secTag)

    # Define constant axial load
    timeSeries('Constant', n1)
    pattern('Plain', n1, n1)
    load(n2, axialLoad, 0.0, 0.0)

    # Define analysis parameters
    wipeAnalysis()
    integrator('LoadControl', 1.0)
    system('SparseGeneral', '-piv')
    test('NormUnbalance', 1e-9, 10)
    numberer('Plain')
    constraints('Plain')
    algorithm('Newton')
    analysis('Static')

    # Do one analysis for constant axial load
    analyze(1)
    loadConst('-time',0.0)

    # Define reference moment
    timeSeries('Linear', n2)
    pattern('Plain',n2, n2)
    load(n2, 0.0, 0.0, 1.0)

    # Compute curvature increment
    dK = maxK / numIncr

    # Use displacement control at node 2 for section analysis
    wipeAnalysis()
    integrator('DisplacementControl', n2,3,dK,1,dK,dK)
    system('SparseGeneral', '-piv')
    test('NormUnbalance', 1e-9, 10)
    numberer('Plain')
    constraints('Plain')
    algorithm('Newton')
    analysis('Static')
    
    M = [0]
    curv = [0]
    
    # Do the section analysis
    for i in range(numIncr):
        analyze(1)
        curv.append(nodeDisp(n2,3))
        M.append(getTime())
    plt.figure()
    plt.plot(curv,M)
    plt.xlabel('Curvatura')
    plt.ylabel('Momento (kN-m)')
    # wipe()
    nodes = [n1,n2]
    return M,curv

def testMaterial(matTag,displ):
    '''
    

    Parameters
    ----------
    matTag : int
        tag of the material.
    displ : list
        list with the peaks of the displacement cycles.

    Returns
    -------
    Disp : array
        DESCRIPTION.
    F : TYPE
        DESCRIPTION.

    '''
    
    
    # wipe()
    
    model('basic','-ndm',2,'-ndf',3)
    # h = getNodeTags()

    node(100,0.0,0.0)
    node(200,0.0,0.0)
    
    fix(100,1,1,1)
    fix(200,1,1,0)
    
    controlnode = 200
    element('zeroLength',1,100,200,'-mat',matTag,'-dir',6)
    
    recorder('Node','-file','MPhi.out','-time','-node',2,'-dof',3,'disp')
    recorder('Element','-file','Moment.out','-time','-ele',1,'force')
    
    ratio = 1/1000
    
    timeSeries('Linear',1)
    pattern('Plain',1,1)
    load(200,0.0,0.0,1.0)
    
    constraints('Plain')
    numberer('Plain')
    system('BandGeneral')
    test('EnergyIncr',1e-6,1000)
    algorithm('Newton')
    
    currentDisp = 0.0
    Disp = [0]
    F = [0]
    nSteps = 1000
    
    for i in displ:
        Dincr = ratio*i/nSteps
        integrator('DisplacementControl',controlnode,3,Dincr)
        analysis('Static')
        
        if Dincr > 0:
            Dmax = Dincr*nSteps
            ok = 0
            while ok == 0 and currentDisp < Dmax:
                ok = analyze(1)
                currentDisp = nodeDisp(controlnode,3)
                F.append(getTime())
                Disp.append(currentDisp)
        elif Dincr < 0:
            Dmax = Dincr*nSteps
            ok = 0
            while ok == 0 and currentDisp > Dmax:
                ok = analyze(1)
                currentDisp = nodeDisp(controlnode,3)
                F.append(getTime())
                Disp.append(currentDisp)
    Fcurr = getTime()
    if ok != 0:
        print('Fallo la convergencia en ',Fcurr)
    else:
        print('Analisis completo')
    
    plt.figure()
    plt.plot(Disp,F)
    plt.xlabel('deformación unitaria (m/m)')
    plt.ylabel('esfuerzo (kPa)')
    return Disp,F
    
def BuildRCSection(ID,HSec,BSec,coverH,coverB,coreID,coverID,steelID,numBarsTop,barAreaTop,numBarsBot,barAreaBot,numBarsIntTot,barAreaInt,nfCoreY,nfCoreZ,nfCoverY,nfCoverZ):
    '''
    Define a procedure which generates a rectangular reinforced concrete section with one layer of steel at the top & bottom, skin reinforcement and a confined core.
	by: Silvia Mazzoni, 2006, adapted from Michael H. Scott, 2003

    Parameters
    ----------
    ID : int
        unique ID for the section.
    HSec : float
        depth of section, along local-y axis.
    BSec : float
        width of section, along local-z axis.
    coverH : float
        cover along the section height.
    coverB : float
        cover along the section width.
    coreID : int
        material tag for the core patch.
    coverID : int
        material tag for the cover patches.
    steelID : int
        material tag for the reinforcing steel.
    numBarsTop : int
        DESCRIPTION.
    barAreaTop : float
        cross-sectional area of each reinforcing bar in top layer
    numBarsBot : int
        number of reinforcing bars in the bottom layer
    barAreaBot : float
        cross-sectional area of each reinforcing bar in bottom layer
    numBarsIntTot : int
        TOTAL number of reinforcing bars on the intermediate layers, symmetric about z axis and 2 bars per layer-- needs to be an even integer
    barAreaInt : float
        cross-sectional area of each reinforcing bar in intermediate layer
    nfCoreY : int
        number of fibers in the core patch in the y direction
    nfCoreZ : int
        number of fibers in the core patch in the z direction
    nfCoverY : int
        number of fibers in the cover patches with long sides in the y direction
    nfCoverZ : int
        number of fibers in the cover patches with long sides in the z direction

    Returns
    -------
    Section defined in your model

    '''
    
    
    # Define a procedure which generates a rectangular reinforced concrete section
	# with one layer of steel at the top & bottom, skin reinforcement and a 
	# confined core.
	#		by: Silvia Mazzoni, 2006
	#			adapted from Michael H. Scott, 2003
	# 
	# Formal arguments
	#    id - tag for the section that is generated by this procedure
	#    HSec - depth of section, along local-y axis
	#    BSec - width of section, along local-z axis
	#    cH - distance from section boundary to neutral axis of reinforcement
	#    cB - distance from section boundary to side of reinforcement
	#    coreID - material tag for the core patch
	#    coverID - material tag for the cover patches
	#    steelID - material tag for the reinforcing steel
	#    numBarsTop - number of reinforcing bars in the top layer
	#    numBarsBot - number of reinforcing bars in the bottom layer
	#    numBarsIntTot - TOTAL number of reinforcing bars on the intermediate layers, symmetric about z axis and 2 bars per layer-- needs to be an even integer
	#    barAreaTop - cross-sectional area of each reinforcing bar in top layer
	#    barAreaBot - cross-sectional area of each reinforcing bar in bottom layer
	#    barAreaInt - cross-sectional area of each reinforcing bar in intermediate layer 
	#    nfCoreY - number of fibers in the core patch in the y direction
	#    nfCoreZ - number of fibers in the core patch in the z direction
	#    nfCoverY - number of fibers in the cover patches with long sides in the y direction
	#    nfCoverZ - number of fibers in the cover patches with long sides in the z direction
    
    coverY = HSec/2.0
    coverZ = BSec/2.0
    coreY = coverY - coverH
    coreZ = coverZ - coverB
    numBarsInt = int(numBarsIntTot/2)
    GJ = 1e6
    nespacios = numBarsInt + 1
    a = HSec - 2*coverH
    b = a/nespacios
        
    section('Fiber',ID,'-GJ',GJ)
    patch('quad',coreID,nfCoreZ,nfCoreY,-coreY,coreZ,-coreY,-coreZ,coreY,-coreZ,coreY,coreZ)
    patch('quad',coverID,2,nfCoverY,-coverY,coverZ,-coreY,coreZ,coreY,coreZ,coverY,coverZ)
    patch('quad',coverID,2,nfCoverY,-coreY,-coreZ,-coverY,-coverZ,coverY,-coverZ,coreY,-coreZ)
    patch('quad',coverID,nfCoverZ,2,-coverY,coverZ,-coverY,-coverZ,-coreY,-coreZ,-coreY,coreZ)
    patch('quad',coverID,nfCoverZ,2,coreY,coreZ,coreY,-coreZ,coverY,-coverZ,coverY,coverZ)    
    layer('straight',steelID,numBarsInt,barAreaInt,-coreY+b,coreZ,coreY-b,coreZ) # este
    layer('straight',steelID,numBarsInt,barAreaInt,-coreY+b,-coreZ,coreY-b,-coreZ) # y este
    layer('straight',steelID,numBarsTop,barAreaTop,coreY,coreZ,coreY,-coreZ)
    layer('straight',steelID,numBarsBot,barAreaBot,-coreY,coreZ,-coreY,-coreZ)
    
    
        

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


def nse(pred,obs):
    '''
    Calculates the normalized nash sutcliffe efficiency index

    Parameters
    ----------
    pred : numpy array
        time series predicted by the numerical model 
    obs : numpy array
        observed time series (for instance, experimental test).

    Returns
    -------
    ns : float
        normalized nash-sutcliffe index.

    '''
    # calcula el índice normalizado de Nash-Sutcliffe
    # pred es la predicción numérica
    # obs es el ensayo
    mean = np.mean(obs)
    denom = np.sum((obs-mean)**2)
    nume = np.sum((obs-pred)**2)
    ns1 = 1-(nume/denom)
    ns = 1/(2-ns1)
    return ns


def kge(pred,obs):
    '''
    Calculates the Kling-Gupta efficienty index

    Parameters
    ----------
    pred : numpy array
        time series predicted by the numerical model 
    obs : numpy array
        observed time series (for instance, experimental test).

    Returns
    -------
    kge : float
        Kling-Gupta index.

    '''
    # calcula el índice de Kling Gupta
    # pred es la predicción numérica
    # obs es el ensayo
    r = np.corrcoef(pred,obs)[0,1]
    ssim = np.std(pred)
    sobs = np.std(obs)
    msim = np.mean(pred)
    mobs = np.mean(obs)
    cr = r-1
    cs = ssim/sobs - 1
    cm = msim/mobs - 1
    kge = 1-np.sqrt(cr**2 + cs**2 + cm**2)
    return kge


def newmarkL(T,xi,GM,delta_t,betha = 1/4, gamma = 1/2 ,u0 = 0,v0 = 0,P0 = 0):
    #T: periodo de la estrutura
    #xi: porcentaje de amortiguamiento crítico
    #GM: registro en unidades consistentes
    #delta_t: delta de tiempo del registro
    #betha, gamma: parámetros del Newmark. Por defecto utiliza método lineal de interpolación
    #u0,v0,a0: condiciones iniciales de desplazamiento velocidad y aceleración
    
    w = 2*np.pi/T
    m = 1.0
    k = m*w**2
    c = 2*xi*m*w
    # Calculos iniciales
    a0 = (P0-c*v0-k*u0)/m                     # Aceleración inicial
    a1 = m/(betha*(delta_t**2))+(gamma*c)/(betha*delta_t)
    a2 = m/(betha*delta_t)+(gamma/betha-1)*c
    a3 = (1/(2*betha)-1)*m+delta_t*(gamma/(2*betha)-1)*c
    k_g = k+a1
    
    # INICIAR VARIABLES
    
    Npts = len(GM)+1
    Desplz = np.zeros((Npts,1))
    Vel = np.zeros((Npts,1))
    Acel = np.zeros((Npts,1))
    Tiempo = np.linspace(0,delta_t*(Npts-1),Npts)
    P_1 = GM*m
    P_2 = np.zeros((Npts,1))
    
    for i in range(Npts-1):
        P_2[i+1] = P_1[i] + a1*Desplz[i] + a2*Vel[i] + a3*Acel[i]
        Desplz[i+1] = P_2[i+1]/k_g
        Vel[i+1] = gamma/(betha*delta_t)*(Desplz[i+1]-Desplz[i])+(1-(gamma/betha))*Vel[i]+delta_t*(1-(gamma/(2*betha)))*Acel[i]
        Acel[i+1] = (Desplz[i+1]-Desplz[i])/(betha*(delta_t**2))-Vel[i]/(betha*delta_t)-(1/(2*betha)-1)*Acel[i]
        
    return Tiempo,Desplz,Vel,Acel


def newmarkLA(T,xi,GM,delta_t,flag = 'all',betha = 1/4, gamma = 1/2 ,u0 = 0,v0 = 0,P0 = 0):
    '''
    Calculates the response of a SDOF based on the Newmark method
    \n T: period of SDOF
    \n xi: percent of critical damping
    \n GM: ground motion acceleration. use consistent units
    \n delta_t: dt of the record
    \n flag: use 'max' to obtain maximum values of displacement, velocity and absolute acceleration
    \n betha, gamma: parameters of Newmark Method. 
    \n u0,v0,a0: initial conditions for displacement, velocity and acceleration
    '''
    
    w = 2*np.pi/T
    m = 1.0
    k = m*w**2
    c = 2*xi*m*w
    # Calculos iniciales
    a0 = (P0-c*v0-k*u0)/m                     # Aceleración inicial
    a1 = m/(betha*(delta_t**2))+(gamma*c)/(betha*delta_t)
    a2 = m/(betha*delta_t)+(gamma/betha-1)*c
    a3 = (1/(2*betha)-1)*m+delta_t*(gamma/(2*betha)-1)*c
    k_g = k+a1
    
    # INICIAR VARIABLES
    
    Npts = len(GM)
    Desplz = np.zeros((Npts))
    Vel = np.zeros((Npts))
    Acel = np.zeros((Npts))
    Tiempo = np.linspace(0,delta_t*(Npts-1),Npts)
    P_1 = GM*m
    P_2 = np.zeros((Npts))
    
    for i in range(Npts-2):
        P_2[i+1] = P_1[i] + a1*Desplz[i] + a2*Vel[i] + a3*Acel[i]
        Desplz[i+1] = P_2[i+1]/k_g
        Vel[i+1] = gamma/(betha*delta_t)*(Desplz[i+1]-Desplz[i])+(1-(gamma/betha))*Vel[i]+delta_t*(1-(gamma/(2*betha)))*Acel[i]
        Acel[i+1] = (Desplz[i+1]-Desplz[i])/(betha*(delta_t**2))-Vel[i]/(betha*delta_t)-(1/(2*betha)-1)*Acel[i]
    
    AcelAbs = Acel + GM # Aquí se calcula la aceleración absoluta
    
    if flag == 'max':
        TT = np.max(np.abs(Tiempo))
        DD = np.max(np.abs(Desplz))
        VV = np.max(np.abs(Vel))
        AA = np.max(np.abs(AcelAbs))
    else:
        TT = Tiempo
        DD = Desplz
        VV = Vel
        AA = Acel
    
    return TT,DD,VV,AA

def spectrum2(GM,delta_t,xi):
    N = 400
    T = np.linspace(0.02,3,N)
    Sa = np.zeros(N)
    U = np.zeros(N)
    V = np.zeros(N)
    
    for i,per in enumerate(T):
        w = 2*np.pi/per
        Tiempo,Desplz,Vel,Acel = newmarkLA(per,xi,GM,delta_t,'max')
        U[i] = Desplz
        V[i] = Vel
        Sa[i] = Desplz*w**2
    return T,Sa


def spectrum4(GM,dt,xi=0.05,rango=[0.02,3.0],N=300):
    ''' Calculates the Sa spectrum for a record
        \n GM: el registro en .txt. Por ejemplo 'GM01.txt'
        \n dt: dt del registro
        \n xi: porcentaje del amortiguamiento crítico
        \n rango: rango de periodos en un vector
        \n N: número de puntos
           
    '''
    m = 1
    T = np.linspace(rango[0],rango[1],N)
    w = 2*np.pi/T
    k = m*w**2
    c = 2*xi*m*w
    Sa = np.zeros(N)
    U = np.zeros(N)
    A = np.zeros(N)
    dmax,amax = [],[]
    for indx, frec in enumerate(w):
        umax,ufin,uperm,amax,tamax = sdfResponse(m,xi,k[indx],1e16,0.05,dt,GM,dt)
        U[indx] = umax
        Sa[indx] = umax*frec**2
        A[indx] = amax
    return T,Sa,U,A

# def espectroNSR(Aa,Av,Fa,Fv,I):
#     T = np.linspace(0,4,500)
#     T0 = 0.1*(Av*Fv)/(Aa*Fa)
#     Tc = 0.48*(Av*Fv)/(Aa*Fa)
#     Tl = 2.4*Fv
#     Sa = (T < T0)*2.5*Aa*Fa*I*(0.4+0.6*T/T0) + ((T0 < T) & (T < Tc))*2.5*Aa*Fa*I + ((Tc < T) & (T < Tl))*1.2*Av*Fv*I/T + (Tl < T)*1.2*Av*Fv*I*Tl/T**2
#     return T,Sa

# T,Sa = espectroNSR(0.15, 0.2, 2.1, 3.2, 1.0)

# plt.plot(T,Sa)

def creategrid(xloc,yloc):
    ''' Creates a 2D grid 
        \n xloc: coordinates of x
        \n yloc: coordinates of y
        Function returns nodes 1000,1001... for the first vertical axis, then 2000,2001... and so on
    '''
    ny = len(yloc)
    nx = len(xloc)
    # ----------------------Crear nodos de la estructura----------------------|
    for i in range(nx):
        for j in range(ny):
            nnode = 1000*(i+1)+j
            node(nnode,xloc[i],yloc[j])

def creategrid3D(xloc,yloc,zloc,dia=1,floor_mass=[1.0]):
    '''
    
    Parameters
    ----------
    xloc : list
        List of X coordinates.
    yloc : list
        List of Y coordinates.
    zloc : list
        List of Z coordinates.
    dia : int, optional
        DESCRIPTION. The default is 0. Change to 1 for a Rigid Diaphragm
    floor_mass : list, optional
        List with the masses per floor starting from the first floor to the roof. It must have one less than zloc.
    floor_inertia : list, optional
        List with the inertia per floor starting from the first floor to the roof. It must have one less than zloc.

    Returns
    -------
    df : TYPE
        DESCRIPTION.

    '''
    ny = len(yloc)
    nx = len(xloc)
    nz = len(zloc)
    coord = []
    
    # ----------------------Crear nodos de la estructura----------------------|
    for i in range(nx):
        for j in range(ny):
            for z in range(nz):
                nnode = 10000*(i+1)+100*(j)+z
                node(nnode,xloc[i],yloc[j],zloc[z])
                coord.append([nnode,xloc[i],yloc[j],zloc[z],z])

    df = pd.DataFrame(coord, columns=['nlabel','x','y','z','floor'])
    if dia == 1:
        for z in range(1,len(zloc)):
            node(z,np.max(xloc)/2,np.max(xloc)/2,zloc[z])
            nfloor = df[df['floor']==z]
            nodes_floor = nfloor.nlabel.to_list()            
            fix(z,0,0,1,1,1,0)
            mass(z,floor_mass[z-1],floor_mass[z-1],floor_mass[z-1],0.0,0.0,0.0)
            rigidDiaphragm(3,z,*nodes_floor)
    
    return df

def BuildISection(secID,matID,d,tw,bf,tf,nfdw,nftw,nfbf,nftf):
    ''' función para crear una sección en forma de I. 
        \n Recibe:
        \n secID: ID para la sección a crear
        \n matID: material de la sección
        \n d: altura total de la sección
        \n tw: espesor del alma
        \n bf: ancho de la aleta
        \n tf: espesor de la aleta
        \n nfdw, nftw, nfbf, nftf: fibras a lo largo de la longitud del alma, ancho del alma, ancho y espesor de la aleta
    '''
    
    
    dw = d-2*tf
    y1 = -d/2
    y2 = -dw/2 
    y3 = dw/2
    y4 = d/2 
    
    z1 =-bf/2
    z2 =-tw/2 
    z3 = tw/2 
    z4 = bf/2 
    GJ = 1e6
    section('Fiber',secID,'-GJ',GJ)
    patch('quad',matID,nfbf,nftf,y1,z4,y1,z1,y2,z1,y2,z4)
    patch('quad',matID,nftw,nfdw,y2,z3,y2,z2,y3,z2,y3,z3)
    patch('quad',matID,nfbf,nftf,y3,z4,y3,z1,y4,z1,y4,z4)

def plot_Wall_T_BE(matConf, matInco, bW, bF, BEU, BED, BEL, BER, Lww, LwF, nMax, nMin):
    
    cover = 0.025
    SteelTag = 4
    
    fib_sec_1 = [['section', 'Fiber', 2, '-GJ', 1e16],
         ['patch','rect',matInco, nMax, nMin, -Lww/2, -bW/2, Lww/2, bW/2], 
         ['patch','rect',matConf, int(nMax/3), nMin, Lww/2, -bW/2, Lww/2+BEU, bW/2],
         ['patch','rect',matConf, int(nMax/3), nMin, -Lww/2-BED, -bW/2, -Lww/2, bW/2],
         # ['patch','rect',matInco, nMin, int(nMax/3), Lww/2+BEU-bF, -bW/2-LwF, Lww/2+BEU, -bW/2],
         # ['patch','rect',matInco, nMin, int(nMax/3), Lww/2+BEU-bF, bW/2, Lww/2+BEU, bW/2+LwF],
         # ['patch','rect',matConf, nMin, int(nMax/3), Lww/2+BEU-bF, bW/2+LwF, Lww/2+BEU, bW/2+LwF+BEL],
         # ['patch','rect',matConf, nMin, int(nMax/3), Lww/2+BEU-bF, -bW/2-LwF-BER, Lww/2+BEU, -bW/2-LwF],
         ['layer','straight',SteelTag,3,area_bar4,2.48,-0.075+cover,2.025+cover,-0.075+cover],    ####Refuerzo en el alma confinamientos
         ['layer','straight',SteelTag,3,area_bar4,2.48,0.075-cover,2.025+cover,0.075-cover],
         ['layer','straight',SteelTag,3,area_bar4,-2.48,-0.075+cover,-2.025-cover,-0.075+cover],
         ['layer','straight',SteelTag,3,area_bar4,-2.48,0.075-cover,-2.025-cover,0.075-cover],
         
         ['layer','straight',SteelTag,27,area_malla7,-2+cover,-0.075+cover,2.0-cover,-0.075+cover],
         ['layer','straight',SteelTag,27,area_malla7,-2+cover,0.075-cover,2.0-cover,0.075-cover],
         
         # ['layer','straight',SteelTag,3,area_bar4,2.5-cover,-2.5+cover,2.5-cover,-2.5-cover+BEf],  #####Refurezo en patin confinamientos
         # ['layer','straight',SteelTag,3,area_bar4,2.35+cover,-2.5+cover,2.35+cover,-2.5-cover+BEf],
         # ['layer','straight',SteelTag,3,area_bar4,2.5-cover,2.5-cover,2.5-cover,2.5+cover-BEf],
         # ['layer','straight',SteelTag,3,area_bar4,2.35+cover,2.5-cover,2.35+cover,2.5+cover-BEf],
         
         ['layer','straight',SteelTag,13,area_malla7,2.475-cover,-2.50+cover+BEf,2.475-cover,-2.50-cover+BEf+Flange],  #superior
         ['layer','straight',SteelTag,13,area_malla7,2.375+cover,-2.50+cover+BEf,2.375+cover,-2.50-cover+BEf+Flange],    #inferior
         ['layer','straight',SteelTag,13,area_malla7,2.475-cover,2.50-cover-BEf,2.475-cover,2.50+cover-BEf-Flange],    #superior
         ['layer','straight',SteelTag,13,area_malla7,2.375+cover,2.50-cover-BEf,2.375+cover,2.50+cover-BEf-Flange]]      #inferior
 
    # plt.gca().invert_xaxis()
   
    return fib_sec_1
  
    matcolor = ['r', 'lightgrey', 'gold', 'w', 'w', 'w']
    opsv.plot_fiber_section(fib_sec_1, matcolor=matcolor)
    plt.axis('equal')
    plt.axhline(y=0, color='r', lw = 0.5)
    plt.axvline(x=0, color='r', lw = 0.5)
    # plt.ylim(0,14)
    # plt.xlim(-6,6)
    plt.gca().invert_xaxis()
    
def dackal(Fyy, Fuu, eyy, ehh, euu, Lb, Db):
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

def residual_disp(drifts,npts):
    ''' Calcula el drift residual de una estructura sometida a un terremoto
        Recibe dos entradas:
            drifts contiene los drifts a lo largo del sismo
            npts es el número de puntos hasta donde llega el registro
        
        Este algoritmo requiere que se haya corrido un periodo de vibración libre luego del final del registro
        
    '''
    freevib = drifts[npts:-1] # extrae los valores de drifts a partir de donde comenzó la vibración libre
    peaks_ind = argrelextrema(freevib, np.greater) # calcula los indices de los puntos de los máximos
    peaks_ind = peaks_ind[0]
    valleys_ind = argrelextrema(freevib, np.less) # calcula los indices de los puntos de los mínimos
    valleys_ind = valleys_ind[0]
    drifts1 = np.abs(freevib[peaks_ind]) # identifica los drifts de los picos
    drifts2 = np.abs(freevib[valleys_ind]) # identifica los drifts de los valles
    resdrift = np.mean(np.concatenate((drifts1,drifts2))) # promedia los drifts de picos y valles
    return resdrift
                       
def Sa_avg(T,Sa,T2 = np.linspace(0.02,3,299)):
    sa_avg = np.zeros(len(T2))
    # sa_avg2 = np.zeros(len(T2)) # en casi que queramos definir con media aritmetica
    for ind,tt in enumerate(T2):
        ta = 0.2*tt # entre 0.2T
        tb = 2.5*tt # y 2.5T
        periods = np.linspace(ta,tb,int(np.ceil((tb-ta)/0.01))) # intervalos de 0.1
        Sas2 = np.interp(periods,T,Sa)
        sa_avg[ind] = gmean(Sas2)
        # sa_avg2[ind] = np.mean(Sas2) # definir con media aritmetica
    return T2,sa_avg


def EAF(t,a):
    '''
    Parameters
    ----------
    t : numpy array
        array with the time of the ground motion.
    a : numpy array
        array with the ground motion acceleration.

    Returns
    -------
    T : numpy array
        array with the periods.
    A:  numpy array
        array with the fourier amplitud.

    '''
    
    
    N = len(a)
    td = t[-1]
    # dt = td/N
    dw = 2*np.pi/td
    Nf = int(N/2 + 1)
    NT = np.linspace(0,Nf-1,Nf)
    W = NT*dw
    
    TF = fft(a)/N*td
    TF1 = TF[0:Nf]
    A = np.zeros(Nf)
    for i in range(Nf):
        A[i] = np.linalg.norm(TF1[i])
    
    T = 2*np.pi*np.reciprocal(W[1::])
    # F = W/(2*np.pi)
    return T,A[1::]


def cumAI(tiempo,sismo1,plot=1,cum=[0.05,0.95]):
    '''  
    Parameters
    ----------
    tiempo : time of the seismic record
    sismo1 : history of accelerations
    plot : optional
        Select 1 if you want the plot. The default is 1.
    cum : optional
        Range of the . The default is [0.05,0.95].

    Returns
    -------
    a2 : cumulative of the arias intensity.
    t1 : start and end time where the record is between the percentages specified in cum.

    '''
  
    IA = np.pi/(2)*np.trapz(sismo1**2,tiempo)
    a = np.pi/(2)*cumulative_trapezoid(sismo1**2,tiempo)/IA
    a2 = np.hstack((0,a))
    t1 = np.interp(cum, a2, tiempo)
    if plot == 1:
        plt.plot(tiempo,a2)
        plt.plot(t1,[0.05,0.95],'ro')
        plt.show()
    return a2,t1

def e20Lobatto(Gfc,Lel,npint,fc,E,e0):
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

def col_materials(fcn=28,fy=420,detailing='DES',tension = 'tension',steeltag = int(100), unctag = int(102), conftag = int(101)):
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
    e20_c = e20Lobatto2(fcn, 3000, 5, fcn, Ec/1000, ec)
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
        e20_cc = e20Lobatto2(2*fcn*k, 3000, 5, fcn*k, Ec/1000, ecc)
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
        e20_cc = e20Lobatto2(2*fcn*k, 3000, 5, fcn*k, Ec/1000, ecc)
        eucc=e20_cc
        if tension == 'tension':
            uniaxialMaterial('Concrete02', conftag, fcc, ecc, fucc, eucc)
        else:
            uniaxialMaterial('Concrete01', conftag, fcc, ecc, fucc, eucc)
    return [unctag,conftag,steeltag]


def create_elements(coordx,coordy,coltag,beamtag,dia = 1):
    '''
    Function to create columns and beam elements. By default it uses one column and one beam.
    
    Parameters
    ----------
    coordx : list
        coordinates in X direction. Must have been used before in the creategrid() command.
    coordy : list
        coordinates in Y direction. Must have been used before in the creategrid() command.
    coltag : integer
        tag of the columns.
    beamtag : integer
        tag of the beams.
    dia : integer, optional
        Use 1 if you want a rigid diaphragm, any number if you don't. The default is 1.

    Returns
    -------
    TagColumns : list
        list with the tags of the columns.
    TagVigas : list
        list with the tags of the beams.

    '''
    lineal = 1
    geomTransf('Linear',lineal)
    pdelta = 2
    geomTransf('PDelta',pdelta)
    nx = len(coordx)
    ny = len(coordy)
    TagColumns = []
    for i in range(ny-1):
        for j in range(nx):
            nodeI = 1000*(j+1)+i
            nodeJ = 1000*(j+1)+(i+1)
            eltag = 100*(j+1) + i
            TagColumns.append(eltag)
            element('forceBeamColumn',eltag,nodeI,nodeJ,pdelta,coltag)
    TagVigas = []
    for i in range(1,ny):
        for j in range(nx-1):
            nodeI = 1000*(j+1)+i
            nodeJ = 1000*(j+2)+i
            eltag = 10000*(j+1) + i
            TagVigas.append(eltag)
            element('forceBeamColumn',eltag,nodeI,nodeJ,lineal,beamtag)
    if dia == 1:
        for j in range(1,ny):
            for i in range(1,nx):
                masternode = 1000 + j
                slavenode = 1000*(i+1) + j
                equalDOF(masternode,slavenode,1)
                
    
    return TagColumns, TagVigas



def create_elements3D(coordx,coordy,coordz,coltag,beamtagX,beamtagY,dia = 1):
    '''
    Function to create columns and beam elements. By default it uses one column and one beam.
    
    Parameters
    ----------
    coordx : list
        coordinates in X direction. Must have been used before in the creategrid() command.
    coordy : list
        coordinates in Y direction. Must have been used before in the creategrid() command.
    coordZ : list
        coordinates in Z direction. Must have been used before in the creategrid() command
    coltag : integer or list
        tag of the columns. IF the user inputs a list, it must be one column section per floor
    beamtagX : integer
        tag of the beams in X direction.
    beamtagY : integer
        tag of the beams in X direction.


    Returns
    -------
    TagColumns : list
        list with the tags of the columns.
    TagVigasX : list
        list with the tags of the beams in the X direction.
    TagVigasY : list
        list with the tags of the beams in the Y direction.

    '''
           
    coltrans = 1000000
    vigXtrans = 2000000
    vigYtrans = 3000000
    geomTransf('PDelta',coltrans,*[0, -1, 0])
    geomTransf('Linear',vigXtrans,*[0, -1, 0])
    geomTransf('Linear',vigYtrans,*[1, 0, 0])
    nx = len(coordx)
    ny = len(coordy)
    nz = len(coordz)
    if len(coltag) != nz-1:
        print('ERROR: Number of column tags does not match number of floors')
        
    TagColumns = []
    sectag = []
    for i in range(nx):
        for j in range(ny):
            for z in range(nz-1):
                nodeI = 10000*(i+1)+100*j+z
                nodeJ = 10000*(i+1)+100*j+(z+1)
                eltag = 10000*(z+1) + 100*i + j
                TagColumns.append(eltag)
                if type(coltag) == list:
                    element('forceBeamColumn',eltag,nodeI,nodeJ,coltrans,coltag[z])
                    sectag.append(coltag[z])    
                else:    
                    element('forceBeamColumn',eltag,nodeI,nodeJ,coltrans,coltag)
    # print(TagColumns)
    TagVigasX = []
    for z in range(1,nz):
        for i in range(nx-1):
            for j in range(ny): 
                nodeI = 10000*(i+1)+100*j+z
                nodeJ = 10000*(i+2)+100*j+z
                eltag = 1000000*z + 10000*(j+1) + 100*i
                TagVigasX.append(eltag)
                element('forceBeamColumn',eltag,*[nodeI,nodeJ],vigXtrans,beamtagX)
                # print(nodeI,nodeJ,':',eltag)
                # print('viga ok')
    TagVigasY = []
    for z in range(1,nz):
        for i in range(nx):
            for j in range(ny-1): 
                nodeI = 10000*(i+1)+100*j+z
                nodeJ = 10000*(i+1)+100*(j+1)+z
                eltag = 10000000*z + 10000*(i+1) + 100*j
                TagVigasY.append(eltag)
                element('forceBeamColumn',eltag,*[nodeI,nodeJ],vigYtrans,beamtagY)
                # print(nodeI,nodeJ,':',eltag)
    
    return TagColumns, TagVigasX, TagVigasY, sectag

def load_beams(floor_load,roof_load,tagbeams,tag = 1):
    '''
    Loads the beams of the model

    Parameters
    ----------
    floor_load : float. 
        Load for the beam at floor levels. Negative values are in gravity direction
    roof_load : float
        Load for the beam at floor levels. Negative values are in gravity direction
    tagbeams : list
        tags of the beams. If you want it to work properly use the same tags generated by the create_elements function.
    tag : integer, optional
        Integer to assign to the tag of the load pattern. The default is 1.

    Returns
    -------
    None.

    '''
    spans = int(str(tagbeams[-1])[0])
    roof_tags = tagbeams[-spans:]
    floor_tags = tagbeams[:-spans]
    timeSeries('Linear', tag)
    pattern('Plain',tag,tag)
    eleLoad('-ele',*floor_tags,'-type','beamUniform',floor_load)
    eleLoad('-ele',*roof_tags,'-type','beamUniform',roof_load)
    
def load_beams3D(floor_load_x, roof_load_x, floor_load_y, roof_load_y, tagbeamsX, tagbeamsY, coordx, coordy, tag = 1):
    '''
    Function to add loads to the beams. Inputs for tags must be the ones from the create_elements3D command

    Parameters
    ----------
    floor_load_x : float
        Value of the load for floor beams in the X direction. USE NEGATIVE for gravity direction.
    roof_load_x : float
        Value of the load for roof beams in the X direction. USE NEGATIVE for gravity direction.
    floor_load_y : float
        Value of the load for floor beams in the Y direction. USE NEGATIVE for gravity direction.
    roof_load_y : float
        Value of the load for roof beams in the Y direction. USE NEGATIVE for gravity direction.
    tagbeamsX : list
        tag of beams in the X direction as returned by the create_elements3D command
    tagbeamsY : TYPE
        DESCRIPTION.
    coordx : TYPE
        DESCRIPTION.
    coordy : TYPE
        DESCRIPTION.
    tag : TYPE, optional
        DESCRIPTION. The default is 1.

    Returns
    -------
    None.

    '''
    
    roofx = tagbeamsX[-(len(coordx)-1)*len(coordy):] # gets the tags of the roof beams in X direction
    floorx = tagbeamsX[:-(len(coordx)-1)*len(coordy)] # gets the tags of the floor beams in X direction
    roofy = tagbeamsY[-(len(coordy)-1)*len(coordx):] # gets the tags of the roof beams in Y direction
    floory = tagbeamsY[:-(len(coordy)-1)*len(coordx)] # gets the tags of the floor beams in Y direction
    timeSeries('Linear',tag)
    pattern('Plain',tag,tag)
    eleLoad('-ele',*floorx,'-type','beamUniform',floor_load_x,0.0)
    eleLoad('-ele',*floory,'-type','beamUniform',floor_load_y,0.0)
    eleLoad('-ele',*roofx,'-type','beamUniform',roof_load_x,0.0)
    eleLoad('-ele',*roofy,'-type','beamUniform',roof_load_y,0.0)
    
    
def pushover_loads(coordy, tag_pattern = 1001, nodes = 0):
    '''
    Generates a pushover pattern proportional to the each floor height. Works in combination with the creategrid command. 

    Parameters
    ----------
    coordy : List
        List with the y coordinates of the model including the base coordinate
    tag_pattern : int, optional
        Integer with the pattern tag for the pushover. The default is 1001.
    nodes : list, optional
        List with the node tags where to create the pushover. The default is 0. If you create the nodes using the creategrid command, you shouldn't change it.
        
            
    Returns
    -------
    None. It creates the pattern.

    '''
    puntos = len(coordy)-1
    suma = np.sum(coordy)
    timeSeries('Linear', tag_pattern)
    pattern('Plain',tag_pattern,tag_pattern)
    if nodes == 0:
        for i in range(puntos):
            load(int(1001+i),coordy[i+1]/suma,0,0)
    else:
        for i in range(puntos):
            load(nodes[i],coordy[i+1]/suma,0,0)
            
def pushover_loads3D(coordz, pushdir = 'x', tag_pattern = 1001, nodes = 0):
    '''
    Generates a pushover pattern proportional to the each floor height. Works in combination with the creategrid3D command. 

    Parameters
    ----------
    coordz : List
        List with the z coordinates of the model including the base coordinate
    pushdir : string, optional
        Enter 'x' for the X direction (default), 'y' for the Y direction
    tag_pattern : int, optional
        Integer with the pattern tag for the pushover. The default is 1001.
    nodes : list, optional
        List with the node tags where to create the pushover. The default is 0. If you create the nodes using the creategrid command, you shouldn't change it.
        
            
    Returns
    -------
    None. It creates the pattern.

    '''
    puntos = len(coordz)-1
    suma = np.sum(coordz)
    timeSeries('Linear', tag_pattern)
    pattern('Plain',tag_pattern,tag_pattern)
    if nodes == 0:
        for i in range(puntos):
            if pushdir == 'x':
                load(int(1+i),coordz[i+1]/suma,0,0,0,0,0)
            elif pushdir == 'y':
                load(int(1+i),0,coordz[i+1]/suma,0,0,0,0)
    else:
        for i in range(puntos):
            load(nodes[i],coordz[i+1]/suma,0,0,0,0,0)

def create_rect_RC_section(ID,HSec,BSec,cover,coreID,coverID,steelID,numBarsTop,barAreaTop,numBarsBot,barAreaBot,numBarsIntTot=2,barAreaInt=1e-10):
    BuildRCSection(ID,HSec,BSec,cover,cover,coreID,coverID,steelID,numBarsTop,barAreaTop,numBarsBot,barAreaBot,numBarsIntTot,barAreaInt,10,10,8,8)
    beamIntegration('Lobatto',ID,ID,5)
    
def create_slabs(coordx, coordy, coordz, hslab, Eslab, pois, seclosa = 12345, dens = 0.0, starttag = 0):
    '''
    create_slabs create a solid slab in the area of the model specified by the coordinates. It uses an ElasticMembratePlateSection formulation and the Shell DKGQ. 

    Parameters
    ----------
    coordx : list
        DESCRIPTION.
    coordy : list
        DESCRIPTION.
    coordz : list
        DESCRIPTION.
    hslab : float
        slab height.
    Eslab : float
        modulus of elasticity of the slab material.
    pois : float
        poisson ratio.
    seclosa : integer, optional
        tag for the slab section. The default is 12345 tro avoid conflicts with other tags.
    dens : float, optional
        Density of the slab. The default is 0.0.
    starttag: integer, optional
        Integer to use to start in another numbering scheme. Useful when you want to call the function several times. You enter the last integer of the previous call and it works.
    Returns
    -------
    slabtags : list
        List with the tags of the slabs

    '''
    section('ElasticMembranePlateSection', seclosa, 1.0*Eslab, pois, hslab, dens)
    nx = len(coordx)
    ny = len(coordy)
    nz = len(coordz)
    slabtags = []
    for z in range(1,nz):
        tag = 100*z + starttag
        for i in range(nx - 1):
            for j in range(ny - 1):
                n1 = 10000*(i+1) + 100*j + z
                n2 = 10000*(i+2) + 100*j + z
                n3 = 10000*(i+2) + 100*(j+1) + z
                n4 = 10000*(i+1) + 100*(j+1) + z
                nodoslosa = [n1,n2,n3,n4]
                element('ShellDKGQ', tag, *nodoslosa, seclosa)
                slabtags.append(tag)
                tag = tag + 1
    return slabtags           
    