# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 17:59:44 2022

@author: Orlando
"""
from openseespy.opensees import *
import matplotlib.pyplot as plt
import numpy as np

# ANALISIS DE GRAVEDAD
# =============================
def gravedad():
    
# Create the system of equation, a sparse solver with partial pivoting
    system('BandGeneral')

# Create the constraint handler, the transformation method
    constraints('Transformation')

# Create the DOF numberer, the reverse Cuthill-McKee algorithm
    numberer('RCM')

# Create the convergence test, the norm of the residual with a tolerance of
# 1e-12 and a max number of iterations of 10
    test('NormDispIncr', 1.0e-12, 10, 3)

# Create the solution algorithm, a Newton-Raphson algorithm
    algorithm('Newton')

# Create the integration scheme, the LoadControl scheme using steps of 0.1
    integrator('LoadControl', 0.1)

# Create the analysis object
    analysis('Static')

    ok = analyze(10)
    
    if ok != 0:
        print('Análisis de gravedad fallido')
        sys.exit()
    else:
        print('Análisis de gravedad completado')
        


# ANALISIS PUSHOVER
# =============================

def pushover(Dmax,Dincr,IDctrlNode,IDctrlDOF):
    
    recorder('Node','-file','techo.out','-time','-node',IDctrlNode,'-dof',IDctrlDOF,'disp')
    maxNumIter = 6
    Tol = 1e-8
      
    
    wipeAnalysis()
    constraints('Transformation')
    numberer('Plain')
    system('BandGeneral')
    test('EnergyIncr', Tol, maxNumIter)
    algorithm('Newton')
    
    integrator('DisplacementControl', IDctrlNode, IDctrlDOF, Dincr)
    analysis('Static')
    
    
    Nsteps =  int(Dmax/ Dincr)
    
    ok = analyze(Nsteps)
    print(ok)
    print('Pushover completado sin problemas')
    
    
    tests = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    algoritmo = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}
    
    
    for i in tests:
        for j in algoritmo:
    
            if ok != 0:
                if j < 4:
                    algorithm(algoritmo[j], '-initial')
                    
                else:
                    algorithm(algoritmo[j])
                    
                test(tests[i], Tol, 1000)
                ok = analyze(Nsteps)                            
                print(tests[i], algoritmo[j], ok)             
                if ok == 0:
                    break
            else:
                continue
            
def pushover2(Dmax,Dincr,IDctrlNode,IDctrlDOF,norm=[-1,1],Tol=1e-8):
    '''
    Function to calculate the pushover

    Parameters
    ----------
    Dmax : float
        Maximum displacement of the pushover.
    Dincr : float
        Increment in the displacement.
    IDctrlNode : int
        control node for the displacements.
    IDctrlDOF : int
        DOF for the displacement.
    norm : list, optional
        List that includes the roof displacement and the building weight to normalize the pushover and display the roof drift vs V/W plot. The default is [-1,1].
    Tol : float, optional
        Norm tolerance. The default is 1e-8.

    Returns
    -------
    techo : numpy array
        Numpy array with the roof displacement recorded during the Pushover.
    V : numpy array
        Numpy array with the base shear (when using an unitary patter) recorded during the Pushover. If pattern if not unitary it returns the multiplier

    '''
    # creación del recorder de techo y definición de la tolerancia
    # recorder('Node','-file','techo.out','-time','-node',IDctrlNode,'-dof',IDctrlDOF,'disp')
    maxNumIter = 10
    
      
    # configuración básica del análisis
    wipeAnalysis()
    constraints('Transformation')
    numberer('RCM')
    system('BandGeneral')
    # system('SparseSYM')
    # system('BandSPD')
    # system('ProfileSPD')
    test('EnergyIncr', Tol, maxNumIter)
    algorithm('Newton')    
    integrator('DisplacementControl', IDctrlNode, IDctrlDOF, Dincr)
    analysis('Static')
    
    # Otras opciones de análisis    
    tests = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    algoritmo = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}

    # rutina del análisis
    
    Nsteps =  int(Dmax/ Dincr) 
    dtecho = [nodeDisp(IDctrlNode,IDctrlDOF)]
    Vbasal = [getTime()]
    
    for k in range(Nsteps):
        ok = analyze(1)
        # ok2 = ok;
        # En caso de no converger en un paso entra al condicional que sigue
        if ok != 0:
            print('configuración por defecto no converge en desplazamiento: ',nodeDisp(IDctrlNode,IDctrlDOF))
            for j in algoritmo:
                if j < 4:
                    algorithm(algoritmo[j], '-initial')
    
                else:
                    algorithm(algoritmo[j])
                
                # el test se hace 50 veces más
                test('EnergyIncr', Tol, maxNumIter*50)
                ok = analyze(1)
                if ok == 0:
                    # si converge vuelve a las opciones iniciales de análisi
                    test('EnergyIncr', Tol, maxNumIter)
                    algorithm('Newton')
                    break
                    
        if ok != 0:
            print('Pushover analisis fallido')
            print('Desplazamiento alcanzado: ',nodeDisp(IDctrlNode,IDctrlDOF),'m')
            break
    
        
        dtecho.append(nodeDisp(IDctrlNode,IDctrlDOF))
        Vbasal.append(getTime())
        
    plt.figure()
    plt.plot(dtecho,Vbasal)
    plt.xlabel('desplazamiento de techo (m)')
    plt.ylabel('corte basal (kN)')
    
    techo = np.array(dtecho)
    V = np.array(Vbasal)
    
    
    if norm[0] != -1:
        deriva = techo/norm[0]*100
        VW = V/norm[1]
        plt.figure()
        plt.plot(deriva,VW)
        plt.xlabel('Deriva de techo (%)')
        plt.ylabel('V/W')
    
    return techo, V

def pushover2BD(Dmax,Dincr,IDctrlNode,IDctrlDOF,norm=[-1,1],Tol=1e-8):
    '''
    Runs a bidirectional pushover analysis

    Parameters
    ----------
    Dmax : float
        Maximum displacement of the pushover in the IDctrlDOF direction
    Dincr : float
        Increment in the displacement.
    IDctrlNode : int
        control node for the displacements.
    IDctrlDOF : int
        DOF to control the displacement.
    norm : list, optional
        List that includes the roof displacement and the building weight to normalize the pushover and display the roof drift vs V/W plot. The default is [-1,1].
    Tol : float, optional
        Norm tolerance. The default is 1e-8.

    Returns
    -------
    techo : numpy array
        Numpy array with the roof displacement recorded during the Pushover in the IDctrlDOF direction
    techo2 : numpy array
        Numpy array with the roof displacement recorded during the Pushover in the IDctrlDOF perpendicular direction
    techoT : numpy array
        Numpy array with the resultant roof displacement recorded during the Pushover.
    V : numpy array
        Numpy array with the base shear (when using an unitary patter) recorded during the Pushover. If pattern if not unitary it returns the multiplier

    '''
    
    
    # Para correr análisis Pushover Bidireccionales
    # creación del recorder de techo y definición de la tolerancia
    # recorder('Node','-file','techo.out','-time','-node',IDctrlNode,'-dof',IDctrlDOF,'disp')
    maxNumIter = 10
    
      
    # configuración básica del análisis
    wipeAnalysis()
    constraints('Transformation')
    numberer('RCM')
    system('BandGeneral')
    test('NormUnbalance', Tol, maxNumIter)
    algorithm('Newton')    
    integrator('DisplacementControl', IDctrlNode, IDctrlDOF, Dincr)
    analysis('Static')
    
    # Otras opciones de análisis    
    tests = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    algoritmo = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}

    # rutina del análisis
    
    Nsteps =  int(Dmax/ Dincr) 
    dtecho1 = [nodeDisp(IDctrlNode,IDctrlDOF)]
    if IDctrlDOF == 1:
        dir2 = 2
    else:
        dir2 = 1
    
    dtecho2 = [nodeDisp(IDctrlNode,dir2)]
    
    Vbasal = [getTime()]
    
    for k in range(Nsteps):
        ok = analyze(1)
        # ok2 = ok;
        # En caso de no converger en un paso entra al condicional que sigue
        if ok != 0:
            print('configuración por defecto no converge en desplazamiento: ',nodeDisp(IDctrlNode,IDctrlDOF))
            for j in algoritmo:
                if j < 4:
                    algorithm(algoritmo[j], '-initial')
    
                else:
                    algorithm(algoritmo[j])
                
                # el test se hace 50 veces más
                test('NormUnbalance', Tol, maxNumIter*50)
                ok = analyze(1)
                if ok == 0:
                    # si converge vuelve a las opciones iniciales de análisi
                    test('NormUnbalance', Tol, maxNumIter)
                    algorithm('Newton')
                    break
                    
        if ok != 0:
            print('Pushover analisis fallido')
            print('Desplazamiento alcanzado: ',nodeDisp(IDctrlNode,IDctrlDOF),'m')
            break
    
        
        dtecho1.append(nodeDisp(IDctrlNode,IDctrlDOF)) # direccion de control
        dtecho2.append(nodeDisp(IDctrlNode,dir2)) # direccion perpendicular a la de control
        Vbasal.append(getTime()) 
        
    plt.figure()
    plt.plot(dtecho1,Vbasal)
    plt.xlabel('desplazamiento de techo (m)')
    plt.ylabel('corte basal (kN)')
    
    techo = np.array(dtecho1)
    techo2 = np.array(dtecho2)
    techoT = np.sqrt(techo**2 + techo2**2)
    V = np.array(Vbasal)
    
    
    if norm[0] != -1:
        deriva = techo/norm[0]*100
        VW = V/norm[1]
        plt.figure()
        plt.plot(deriva,VW)
        plt.xlabel('Deriva de techo (%)')
        plt.ylabel('V/W')
    
    return techo, techo2, techoT, V

def pushover2MP(Dmax,Dincr,IDctrlNode,IDctrlDOF,norm=[-1,1],Tol=1e-8):
    
    # creación del recorder de techo y definición de la tolerancia
    recorder('Node','-file','techo.out','-time','-node',IDctrlNode,'-dof',IDctrlDOF,'disp')
    maxNumIter = 20
    
      
    # configuración básica del análisis
    wipeAnalysis()
    constraints('Transformation')
    numberer('RCM')
    system('BandGeneral')
    test('EnergyIncr', Tol, maxNumIter)
    algorithm('Newton')    
    integrator('DisplacementControl', IDctrlNode, IDctrlDOF, Dincr)
    analysis('Static')
    
    # Otras opciones de análisis    
    tests = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    algoritmo = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}

    # rutina del análisis
    
    Nsteps =  int(Dmax/ Dincr) 
    dtecho = [nodeDisp(IDctrlNode,IDctrlDOF)]
    Vbasal = [getTime()]
    
    for k in range(Nsteps):
        ok = analyze(1)
        # ok2 = ok;
        # En caso de no converger en un paso entra al condicional que sigue
        if ok != 0:
            print('configuración por defecto no converge en desplazamiento: ',nodeDisp(IDctrlNode,IDctrlDOF))
            for j in algoritmo:
                if j < 4:
                    algorithm(algoritmo[j], '-initial')
    
                else:
                    algorithm(algoritmo[j])
                
                # el test se hace 50 veces más
                test('EnergyIncr', Tol, maxNumIter*25)
                ok = analyze(1)
                if ok == 0:
                    # si converge vuelve a las opciones iniciales de análisi
                    test('EnergyIncr', Tol, maxNumIter)
                    algorithm('Newton')
                    break
                    
        if ok != 0:
            print('Pushover analisis fallido')
            print('Desplazamiento alcanzado: ',nodeDisp(IDctrlNode,IDctrlDOF),'m')
            break
    
        
        dtecho.append(nodeDisp(IDctrlNode,IDctrlDOF))
        Vbasal.append(getTime())
        

    techo = np.array(dtecho)
    V = np.array(Vbasal)
    
    
    if norm[0] != -1:
        deriva = techo/norm[0]*100
        VW = V/norm[1]
        plt.figure()
        plt.plot(deriva,VW)
        plt.xlabel('Deriva de techo (%)')
        plt.ylabel('V/W')
    
    return techo, V


def pushover2T(Dmax,Dincr,IDctrlNode,IDctrlDOF,norm=[-1,1],Tol=1e-8):
    '''
    Function to calculate the pushover and the building period during this one.

    Parameters
    ----------
    Dmax : float
        Maximum displacement of the pushover.
    Dincr : float
        Increment in the displacement.
    IDctrlNode : int
        control node for the displacements.
    IDctrlDOF : int
        DOF for the displacement.
    norm : list, optional
        List that includes the roof displacement and the building weight to normalize the pushover and display the roof drift vs V/W plot. The default is [-1,1].
    Tol : float, optional
        Norm tolerance. The default is 1e-8.

    Returns
    -------
    techo : numpy array
        Numpy array with the roof displacement recorded during the Pushover.
    V : numpy array
        Numpy array with the base shear (when using an unitary pattern) recorded during the Pushover. If pattern is not unitary it returns the multiplier
    T : numpy array
        Numpy array with the building period recorded during the Pushover.
    '''
    # creación del recorder de techo y definición de la tolerancia
    recorder('Node','-file','techo.out','-time','-node',IDctrlNode,'-dof',IDctrlDOF,'disp')
    maxNumIter = 10
    
      
    # configuración básica del análisis
    wipeAnalysis()
    constraints('Transformation')
    numberer('RCM')
    system('BandGeneral')
    test('EnergyIncr', Tol, maxNumIter)
    algorithm('Newton')    
    integrator('DisplacementControl', IDctrlNode, IDctrlDOF, Dincr)
    analysis('Static')
    
    # Otras opciones de análisis    
    tests = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    algoritmo = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}

    # rutina del análisis
    eig = eigen(1)
    TT = 2*3.1416/np.sqrt(eig[0])
    Nsteps =  int(Dmax/ Dincr) 
    dtecho = [nodeDisp(IDctrlNode,IDctrlDOF)]
    Vbasal = [getTime()]
    periods = [TT]
    fibras1 = [0]*8
    for k in range(Nsteps):
        ok = analyze(1)
        # ok2 = ok;
        # En caso de no converger en un paso entra al condicional que sigue
        if ok != 0:
            print('configuración por defecto no converge en desplazamiento: ',nodeDisp(IDctrlNode,IDctrlDOF))
            for j in algoritmo:
                if j < 4:
                    algorithm(algoritmo[j], '-initial')
    
                else:
                    algorithm(algoritmo[j])
                
                # el test se hace 50 veces más
                test('EnergyIncr', Tol, maxNumIter*50)
                ok = analyze(1)
                if ok == 0:
                    # si converge vuelve a las opciones iniciales de análisi
                    test('EnergyIncr', Tol, maxNumIter)
                    algorithm('Newton')
                    break
                    
        if ok != 0:
            print('Pushover analisis fallido')
            print('Desplazamiento alcanzado: ',nodeDisp(IDctrlNode,IDctrlDOF),'m')
            break
        
       
        eig = eigen(1)
        TT = 2*3.1416/np.sqrt(eig[0])
        periods.append(TT)
         
        dtecho.append(nodeDisp(IDctrlNode,IDctrlDOF))
        Vbasal.append(getTime())
        
    plt.figure()
    plt.plot(dtecho,Vbasal)
    plt.xlabel('desplazamiento de techo (m)')
    plt.ylabel('corte basal (kN)')
    
    techo = np.array(dtecho)
    V = np.array(Vbasal)
    PER = np.array(periods)
    
    
    if norm[0] != -1:
        deriva = techo/norm[0]*100
        VW = V/norm[1]
        plt.figure()
        plt.plot(deriva,VW)
        plt.xlabel('Deriva de techo (%)')
        plt.ylabel('V/W')
    
    plt.figure()
    plt.plot(dtecho,periods)
    plt.xlabel('desplazamiento de techo (m)')
    plt.ylabel('Periodo (s)')
    
    return techo, V, PER

def pushover3T(Dmax,Dincr,IDctrlNode,IDctrlDOF,elements,norm=[-1,1],Tol=1e-8):
    
    # creación del recorder de techo y definición de la tolerancia
    recorder('Node','-file','techo.out','-time','-node',IDctrlNode,'-dof',IDctrlDOF,'disp')
    maxNumIter = 10
    
      
    # configuración básica del análisis
    wipeAnalysis()
    constraints('Transformation')
    numberer('RCM')
    system('BandGeneral')
    test('EnergyIncr', Tol, maxNumIter)
    algorithm('Newton')    
    integrator('DisplacementControl', IDctrlNode, IDctrlDOF, Dincr)
    analysis('Static')
    
    # Otras opciones de análisis    
    tests = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    algoritmo = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}

    # rutina del análisis
    eig = eigen(1)
    TT = 2*3.1416/np.sqrt(eig[0])
    Nsteps =  int(Dmax/ Dincr) 
    dtecho = [nodeDisp(IDctrlNode,IDctrlDOF)]
    Vbasal = [getTime()]
    periods = [TT]
    
    
    nels = len(elements)
    Eds = np.zeros((nels, Nsteps+1, 6)) # para grabar las fuerzas de los elementos
    Curv = np.zeros((nels,Nsteps+1)) # para grabar la curvatura de los elementos
    # Strains = np.zeros((Nsteps+1, 8, nels)) # # para grabar las deformaciones de los muros en las 8 fibras que tienen los elementos
    Strains = np.zeros((nels, Nsteps+1, 8))
    cStress = np.zeros((nels, Nsteps+1, 8)) # # para grabar los esfuerzos del concreto de los muros en las 8 fibras que tienen los elementos
    sStress = np.zeros((nels, Nsteps+1, 8)) # # para grabar los esfuerzos del acero de los muros en las 8 fibras que tienen los elementos
    
    for k in range(Nsteps):
        ok = analyze(1)
        # ok2 = ok;
        # En caso de no converger en un paso entra al condicional que sigue
        if ok != 0:
            print('configuración por defecto no converge en desplazamiento: ',nodeDisp(IDctrlNode,IDctrlDOF))
            for j in algoritmo:
                if j < 4:
                    algorithm(algoritmo[j], '-initial')
    
                else:
                    algorithm(algoritmo[j])
                
                # el test se hace 50 veces más
                test('EnergyIncr', Tol, maxNumIter*50)
                ok = analyze(1)
                if ok == 0:
                    # si converge vuelve a las opciones iniciales de análisi
                    test('EnergyIncr', Tol, maxNumIter)
                    algorithm('Newton')
                    break
                    
        if ok != 0:
            print('Pushover analisis fallido')
            print('Desplazamiento alcanzado: ',nodeDisp(IDctrlNode,IDctrlDOF),'m')
            break
        
        
        for el_i, ele_tag in enumerate(elements):
            
            # Curv[k+1, el_i] = [eleResponse(ele_tag,'Curvature')]
            Eds[el_i , k+1, :] = [eleResponse(ele_tag,'globalForce')[0],
                                 eleResponse(ele_tag,'globalForce')[1],
                                 eleResponse(ele_tag,'globalForce')[2],
                                 eleResponse(ele_tag,'globalForce')[3],
                                 eleResponse(ele_tag,'globalForce')[4],
                                 eleResponse(ele_tag,'globalForce')[5]]
            
            Strains[el_i , k+1, :] = [eleResponse(ele_tag,'Fiber_Strain')[0],
                                 eleResponse(ele_tag,'Fiber_Strain')[1],
                                 eleResponse(ele_tag,'Fiber_Strain')[2],
                                 eleResponse(ele_tag,'Fiber_Strain')[3],
                                 eleResponse(ele_tag,'Fiber_Strain')[4],
                                 eleResponse(ele_tag,'Fiber_Strain')[5],
                                 eleResponse(ele_tag,'Fiber_Strain')[6],
                                 eleResponse(ele_tag,'Fiber_Strain')[7]]
            
            cStress[el_i , k+1, :] = [eleResponse(ele_tag,'Fiber_Stress_Concrete')[0],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[1],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[2],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[3],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[4],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[5],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[6],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[7]]
            
            sStress[el_i , k+1, :] = [eleResponse(ele_tag,'Fiber_Stress_Steel')[0],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[1],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[2],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[3],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[4],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[5],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[6],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[7]]
        
        eig = eigen(1)
        TT = 2*3.1416/np.sqrt(eig[0])
        periods.append(TT)
         
        dtecho.append(nodeDisp(IDctrlNode,IDctrlDOF))
        Vbasal.append(getTime())
        
    plt.figure()
    plt.plot(dtecho,Vbasal)
    plt.xlabel('desplazamiento de techo (m)')
    plt.ylabel('corte basal (kN)')
    
    techo = np.array(dtecho)
    V = np.array(Vbasal)
    PER = np.array(periods)
    
    
    if norm[0] != -1:
        deriva = techo/norm[0]*100
        VW = V/norm[1]
        plt.figure()
        plt.plot(deriva,VW)
        plt.xlabel('Deriva de techo (%)')
        plt.ylabel('V/W')
    
    plt.figure()
    plt.plot(dtecho,periods)
    plt.xlabel('desplazamiento de techo (m)')
    plt.ylabel('Periodo (s)')
    
    return techo, V, PER, Eds, Strains, cStress, sStress



def pushover3Tn(Dmax,Dincr,IDctrlNode,IDctrlDOF,elements,norm=[-1,1],Tol=1e-8):
    
    # creación del recorder de techo y definición de la tolerancia
    recorder('Node','-file','techo.out','-time','-node',IDctrlNode,'-dof',IDctrlDOF,'disp')
    maxNumIter = 10
    
      
    # configuración básica del análisis
    wipeAnalysis()
    constraints('Transformation')
    numberer('RCM')
    system('BandGeneral')
    test('EnergyIncr', Tol, maxNumIter)
    algorithm('Newton')    
    integrator('DisplacementControl', IDctrlNode, IDctrlDOF, Dincr)
    analysis('Static')
    
    # Otras opciones de análisis    
    tests = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    algoritmo = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}

    # rutina del análisis
    eig = eigen(1)
    TT = 2*3.1416/np.sqrt(eig[0])
    Nsteps =  int(Dmax/ Dincr) 
    dtecho = [nodeDisp(IDctrlNode,IDctrlDOF)]
    Vbasal = [getTime()]
    periods = [TT]
    
    
    nels = len(elements)
    Eds = np.zeros((nels, Nsteps+1, 6)) # para grabar las fuerzas de los elementos
    Curv = np.zeros((nels,Nsteps+1)) # para grabar la curvatura de los elementos
    # Strains = np.zeros((Nsteps+1, 8, nels)) # # para grabar las deformaciones de los muros en las 8 fibras que tienen los elementos
    Strains = np.zeros((nels, Nsteps+1, 14))
    cStress = np.zeros((nels, Nsteps+1, 14)) # # para grabar los esfuerzos del concreto de los muros en las 8 fibras que tienen los elementos
    sStress = np.zeros((nels, Nsteps+1, 14)) # # para grabar los esfuerzos del acero de los muros en las 8 fibras que tienen los elementos
    
    for k in range(Nsteps):
        ok = analyze(1)
        # ok2 = ok;
        # En caso de no converger en un paso entra al condicional que sigue
        if ok != 0:
            print('configuración por defecto no converge en desplazamiento: ',nodeDisp(IDctrlNode,IDctrlDOF))
            for j in algoritmo:
                if j < 4:
                    algorithm(algoritmo[j], '-initial')
    
                else:
                    algorithm(algoritmo[j])
                
                # el test se hace 50 veces más
                test('EnergyIncr', Tol, maxNumIter*50)
                ok = analyze(1)
                if ok == 0:
                    # si converge vuelve a las opciones iniciales de análisi
                    test('EnergyIncr', Tol, maxNumIter)
                    algorithm('Newton')
                    break
                    
        if ok != 0:
            print('Pushover analisis fallido')
            print('Desplazamiento alcanzado: ',nodeDisp(IDctrlNode,IDctrlDOF),'m')
            break
        
        
        for el_i, ele_tag in enumerate(elements):
            
            # Curv[k+1, el_i] = [eleResponse(ele_tag,'Curvature')]
            Eds[el_i , k+1, :] = [eleResponse(ele_tag,'globalForce')[0],
                                 eleResponse(ele_tag,'globalForce')[1],
                                 eleResponse(ele_tag,'globalForce')[2],
                                 eleResponse(ele_tag,'globalForce')[3],
                                 eleResponse(ele_tag,'globalForce')[4],
                                 eleResponse(ele_tag,'globalForce')[5]]
            
            Strains[el_i , k+1, :] = [eleResponse(ele_tag,'Fiber_Strain')[0],
                                 eleResponse(ele_tag,'Fiber_Strain')[1],
                                 eleResponse(ele_tag,'Fiber_Strain')[2],
                                 eleResponse(ele_tag,'Fiber_Strain')[3],
                                 eleResponse(ele_tag,'Fiber_Strain')[4],
                                 eleResponse(ele_tag,'Fiber_Strain')[5],
                                 eleResponse(ele_tag,'Fiber_Strain')[6],
                                 eleResponse(ele_tag,'Fiber_Strain')[7],
                                 eleResponse(ele_tag,'Fiber_Strain')[8],
                                 eleResponse(ele_tag,'Fiber_Strain')[9],
                                 eleResponse(ele_tag,'Fiber_Strain')[10],
                                 eleResponse(ele_tag,'Fiber_Strain')[11],
                                 eleResponse(ele_tag,'Fiber_Strain')[12],
                                 eleResponse(ele_tag,'Fiber_Strain')[13]]
                                 
            
            cStress[el_i , k+1, :] = [eleResponse(ele_tag,'Fiber_Stress_Concrete')[0],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[1],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[2],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[3],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[4],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[5],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[6],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[7],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[8],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[9],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[10],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[11],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[12],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[13]]
                                 
            
            sStress[el_i , k+1, :] = [eleResponse(ele_tag,'Fiber_Stress_Steel')[0],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[1],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[2],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[3],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[4],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[5],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[6],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[7],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[8],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[9],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[10],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[11],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[12],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[13]]


        eig = eigen(1)
        TT = 2*3.1416/np.sqrt(eig[0])
        periods.append(TT)
         
        dtecho.append(nodeDisp(IDctrlNode,IDctrlDOF))
        Vbasal.append(getTime())
        
    plt.figure()
    plt.plot(dtecho,Vbasal)
    plt.xlabel('desplazamiento de techo (m)')
    plt.ylabel('corte basal (kN)')
    
    techo = np.array(dtecho)
    V = np.array(Vbasal)
    PER = np.array(periods)
    
    
    if norm[0] != -1:
        deriva = techo/norm[0]*100
        VW = V/norm[1]
        plt.figure()
        plt.plot(deriva,VW)
        plt.xlabel('Deriva de techo (%)')
        plt.ylabel('V/W')
    
    plt.figure()
    plt.plot(dtecho,periods)
    plt.xlabel('desplazamiento de techo (m)')
    plt.ylabel('Periodo (s)')
    
    return techo, V, PER, Eds, Strains, cStress, sStress


# ANALISIS DINAMICO
# =============================   

# dinamico es el más sencillo de todos, corre un terremoto creando un recorder para el techo.
# dinamicoBD corre una pareja de sismos aplicadas al modelo en tres dimensiones.
# dinamicoIDA crea el recorder en función del factor escalar.
# dinamicoAnim es dinamico pero guarda la información para animar el registro
# dinamicoIDA2 está modificado para ser utilizado cuando se desee correr en paralelo los cálculos. Devuelve solo el desplazamiento de techo
# dinamicoIDA3 PARA SER UTILIZADO PARA CORRER EN PARALELO LOS SISMOS Y EXTRAYENDO LAS FUERZAS DE LOS ELEMENTOS INDICADOS EN ELEMENTS. SOLO PUEDEN SER LOS MUROS DE MOMENTO. También extrae desplazamientos de los nodos
# dinamicoIDA4 PARA SER UTILIZADO PARA CORRER EN PARALELO LOS SISMOS Y EXTRAYENDO LAS FUERZAS DE LOS ELEMENTOS INDICADOS EN ELEMENTS. SOLO PUEDEN SER LOS MUROS DE MOMENTO. También extrae desplazamientos de los nodos, aceleraciones, derivas, velocidades, esfuerzos en concreto, acero y deformaciones unitarias de cada muro indicado en elements
# dinamicoIDA4P PARA SER UTILIZADO PARA CORRER EN PARALELO LOS SISMOS Y EXTRAYENDO LAS FUERZAS DE LOS ELEMENTOS INDICADOS EN ELEMENTS. ADMITE VIGAS Y COLUMNAS. También extrae desplazamientos de los nodos, aceleraciones, derivas, velocidades, esfuerzos en concreto, acero y deformaciones unitarias de cada muro indicado en elements
# dinamicoIDA5 es lo mismo que IDA4, pero en lugar del nombre del registro, recibe una lista con las aceleraciones.
# dinamicoIDA6 es lo mismo que IDA4P pero recibe las dos componentes del sismo

def dinamico(recordName,dtrec,nPts,dtan,fact,damp,IDctrlNode,IDctrlDOF,modes = [0,2],Kswitch = 1,Tol=1e-4):
    
    # record es el nombre del registro, incluyendo extensión. P.ej. GM01.txt
    # dtrec es el dt del registro
    # nPts es el número de puntos del análisis
    # dtan es el dt del análisis
    # fact es el factor escalar del registro
    # damp es el porcentaje de amortiguamiento (EN DECIMAL. p.ej: 0.03 para 3%)
    # Kswitch recibe: 1: matriz inicial, 2: matriz actual
    # IDctrlNode,IDctrlDOF son respectivamente el nodo y desplazamiento de control deseados
    
    # creación del recorder de techo y definición de la tolerancia
    recorder('Node','-file','techo.out','-time','-node',IDctrlNode,'-dof',IDctrlDOF,'disp')
    maxNumIter = 25
    
    # creación del pattern
    
    timeSeries('Path',1000,'-filePath',recordName,'-dt',dtrec,'-factor',fact)
    pattern('UniformExcitation',  1000,   IDctrlDOF,  '-accel', 1000)
    
    # damping
    nmodes = max(modes)+1
    eigval = eigen(nmodes)
    
    eig1 = eigval[modes[0]]
    eig2 = eigval[modes[1]]
    
    w1 = eig1**0.5
    w2 = eig2**0.5
    
    beta = 2.0*damp/(w1 + w2)
    alfa = 2.0*damp*w1*w2/(w1 + w2)
    
    if Kswitch == 1:
        rayleigh(alfa, 0.0, beta, 0.0)
    else:
        rayleigh(alfa, beta, 0.0, 0.0)
    
    # configuración básica del análisis
    wipeAnalysis()
    constraints('Transformation')
    numberer('RCM')
    system('BandGeneral')
    test('NormDispIncr', Tol, maxNumIter)
    algorithm('Newton')    
    integrator('Newmark', 0.5, 0.25)
    analysis('Transient')
    
    # Otras opciones de análisis    
    tests = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    algoritmo = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}

    # rutina del análisis
    
    Nsteps =  int(dtrec*nPts/dtan)
    dtecho = [nodeDisp(IDctrlNode,IDctrlDOF)]
    t = [getTime()]
    
    for k in range(Nsteps):
        ok = analyze(1,dtan)
        # ok2 = ok;
        # En caso de no converger en un paso entra al condicional que sigue
        if ok != 0:
            print('configuración por defecto no converge en desplazamiento: ',nodeDisp(IDctrlNode,IDctrlDOF))
            for j in algoritmo:
                if j < 4:
                    algorithm(algoritmo[j], '-initial')
    
                else:
                    algorithm(algoritmo[j])
                
                # el test se hace 50 veces más
                test('NormDispIncr', Tol, maxNumIter*50)
                ok = analyze(1,dtan)
                if ok == 0:
                    # si converge vuelve a las opciones iniciales de análisi
                    test('NormDispIncr', Tol, maxNumIter)
                    algorithm('Newton')
                    break
                    
        if ok != 0:
            print('Análisis dinámico fallido')
            print('Desplazamiento alcanzado: ',nodeDisp(IDctrlNode,IDctrlDOF),'m')
            break
    
        
        dtecho.append(nodeDisp(IDctrlNode,IDctrlDOF))
        t.append(getTime())
        
    # plt.figure()
    # plt.plot(t,dtecho)
    # plt.xlabel('tiempo (s)')
    # plt.ylabel('desplazamiento (m)')  
    
    techo = np.array(dtecho)
    tiempo = np.array(t)
    
    
    
    return tiempo,techo

def dinamicoBD(recordName,dtrec,nPts,dtan,fact,damp,IDctrlNode,IDctrlDOF,modes = [0,2],Kswitch = 1,Tol=1e-4):
    '''
    Performs a dynamic analysis applying both components of a ground motion and recording the displacement of a user selected node.
    
    Parameters
    ----------
    recordName : string
        Name of the record including file extension (i.e., 'GM01.txt'). It must have one record instant per line. 
    dtrec : float
        time increment of the record.
    nPts : integer
        number of points of the record.
    dtan : float
        time increment to be used in the analysis. If smaller than dtrec, OpenSeesPy interpolates.
    fact : float
        scale factor to apply to the record.
    damp : float
        Damping percentage in decimal (i.e., use 0.03 for 3%).
    IDctrlNode : int
        control node for the displacements.
    IDctrlDOF : int
        DOF for the displacement.
    modes : list, optional
        Modes of the structure to apply the Rayleigh damping. The default is [0,2] which uses the first and third mode.
    Kswitch : int, optional
        Use it to define which stiffness matrix should be used for the ramping. The default is 1 that uses initial stiffness. Input 2 for current stifness.
    Tol : float, optional
        Tolerance for the analysis. The default is 1e-4 because it uses the NormUnbalance test.

    Returns
    -------
    tiempo : numpy array
        Numpy array with analysis time.
    techo : numpy array
        Displacement of the control node.
    techo2 : numpy array
            Numpy array with the roof displacement recorded during the analysis in the IDctrlDOF perpendicular direction
    techoT : numpy array
            Numpy array with the resultant roof displacement recorded during the analysis.

    '''
    # record es el nombre de los registros en una lista, incluyendo extensión. P.ej. GM01.txt
    # dtrec es el dt del registro. Debe ser el mismo para ambos
    # nPts es el número de puntos del análisis. Debe ser el mismo para ambos
    # dtan es el dt del análisis
    # fact es el factor escalar del registro. Debe ser el mismo para ambos
    # damp es el porcentaje de amortiguamiento (EN DECIMAL. p.ej: 0.03 para 3%)
    # Kswitch recibe: 1: matriz inicial, 2: matriz actual
    # IDctrlNode,IDctrlDOF son respectivamente el nodo y desplazamiento de control deseados
    
    # creación del recorder de techo y definición de la tolerancia
    # recorder('Node','-file','techo.out','-time','-node',IDctrlNode,'-dof',IDctrlDOF,'disp')
    maxNumIter = 25
    
    if IDctrlDOF == 1:
        dir2 = 2
    else:
        dir2 = 1
        
    # creación del pattern
    timeSeries('Path',1000,'-filePath',recordName[0],'-dt',dtrec,'-factor',fact)
    pattern('UniformExcitation',  1000,   IDctrlDOF,  '-accel', 1000)
    timeSeries('Path',1001,'-filePath',recordName[1],'-dt',dtrec,'-factor',fact)
    pattern('UniformExcitation',  1001,   dir2,  '-accel', 1001)
    
    # damping
    nmodes = max(modes)+1
    eigval = eigen(nmodes)
    
    eig1 = eigval[modes[0]]
    eig2 = eigval[modes[1]]
    
    w1 = eig1**0.5
    w2 = eig2**0.5
    
    beta = 2.0*damp/(w1 + w2)
    alfa = 2.0*damp*w1*w2/(w1 + w2)
    
    if Kswitch == 1:
        rayleigh(alfa, 0.0, beta, 0.0)
    else:
        rayleigh(alfa, beta, 0.0, 0.0)
    
    # configuración básica del análisis
    wipeAnalysis()
    constraints('Transformation')
    numberer('RCM')
    system('BandGeneral')
    test('NormDispIncr', Tol, maxNumIter)
    algorithm('Newton')    
    integrator('Newmark', 0.5, 0.25)
    analysis('Transient')
    
    # Otras opciones de análisis    
    tests = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    algoritmo = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}

    # rutina del análisis
    
    Nsteps =  int(dtrec*nPts/dtan)
    
    dtecho1 = [nodeDisp(IDctrlNode,IDctrlDOF)]
    dtecho2 = [nodeDisp(IDctrlNode,dir2)]

    t = [getTime()]
    
    for k in range(Nsteps):
        ok = analyze(1,dtan)
        # ok2 = ok;
        # En caso de no converger en un paso entra al condicional que sigue
        if ok != 0:
            print('configuración por defecto no converge en desplazamiento: ',nodeDisp(IDctrlNode,IDctrlDOF))
            for j in algoritmo:
                if j < 4:
                    algorithm(algoritmo[j], '-initial')
    
                else:
                    algorithm(algoritmo[j])
                
                # el test se hace 50 veces más
                test('NormDispIncr', Tol, maxNumIter*50)
                ok = analyze(1,dtan)
                if ok == 0:
                    # si converge vuelve a las opciones iniciales de análisi
                    test('NormDispIncr', Tol, maxNumIter)
                    algorithm('Newton')
                    break
                    
        if ok != 0:
            print('Análisis dinámico fallido')
            print('Desplazamiento alcanzado: ',nodeDisp(IDctrlNode,IDctrlDOF),'m')
            break
    
        
        dtecho1.append(nodeDisp(IDctrlNode,IDctrlDOF))
        dtecho2.append(nodeDisp(IDctrlNode,dir2))
        t.append(getTime())
    
    techo1 = np.array(dtecho1)
    techo2 = np.array(dtecho2)
    techoT = np.sqrt(techo1**2+techo2**2)
    tiempo = np.array(t)
    
    return tiempo,techo1,techo2,techoT

def dinamicoBD2(recordName,dtrec,nPts,dtan,fact,damp,IDctrlNode,IDctrlDOF,nodes_control,elements,modes = [0,2],Kswitch = 1,Tol=1e-4):
    '''
    Performs a dynamic analysis applying both components of a ground motion and recording the displacement of a user selected node. To be used ONLY with quadrilateral elements with 24DOF.
    
    ----------
    recordName : list
        list with the names of the record pair including file extension (i.e., 'GM01.txt'). It must have one record instant per line and each record. the records must be pairs, so the function expects that they are of the same length. 
    dtrec : float
        time increment of the record.
    nPts : integer
        number of points of the record.
    dtan : float
        time increment to be used in the analysis. If smaller than dtrec, OpenSeesPy interpolates.
    fact : float
        scale factor to apply to the record.
    damp : float
        Damping percentage in decimal (i.e., use 0.03 for 3%).
    IDctrlNode : int
        control node for the displacements.
    IDctrlDOF : int
        DOF to apply the record. The 
    nodes_control : list
        nodes to compute displacements and inter-story drift. You must input one per floor, otherwise you'll get an error..
    elements : list
        list of elements to record forces.
    modes : list, optional
        Modes of the structure to apply the Rayleigh damping. The default is [0,2] which uses the first and third mode.
    Kswitch : int, optional
        Use it to define which stiffness matrix should be used for the ramping. The default is 1 that uses initial stiffness. Input 2 for current stifness.
    Tol : float, optional
        Tolerance for the analysis. The default is 1e-4 because it uses the NormUnbalance test.

    Returns
    -------
    tiempo : numpy array
        Numpy array with analysis time.
    techo : numpy array
        Numpy array with displacement of the control node in the IDctrlDOF direction
    techo2 : numpy array
        Numpy array with the roof displacement recorded during the analysis in the IDctrlDOF perpendicular direction
    techoT : numpy array
        Numpy array with the resultant roof displacement recorded during the analysis
    node_disp:
        Numpy array with displacement of the control nodes in the IDctrlDOF direction
    node_vel:
        Numpy array with velocoties of the control nodes in the IDctrlDOF direction
    node_acel:
        Numpy array with relative accelerations of the control nodes in the IDctrlDOF direction
    node_disp2:
        Numpy array with displacement of the control nodes in the IDctrlDOF perpendicular direction
    node_acel2:
        Numpy array with relative accelerations of the control nodes in the IDctrlDOF perpendicular direction
    Eds:
        Element forces recorded for the element with tags defined in the input variable elements.
    '''
    
    # Realiza un análisis dinámico aplicando dos componentes ortogonales del terremoto
    # Graba información de elementos y nodos
    
    # record es el nombre de los registros en una lista, incluyendo extensión. P.ej. GM01.txt
    # dtrec es el dt del registro. Debe ser el mismo para ambos
    # nPts es el número de puntos del análisis. Debe ser el mismo para ambos
    # dtan es el dt del análisis
    # fact es el factor escalar del registro. Debe ser el mismo para ambos
    # damp es el porcentaje de amortiguamiento (EN DECIMAL. p.ej: 0.03 para 3%)
    # Kswitch recibe: 1: matriz inicial, 2: matriz actual
    # IDctrlNode,IDctrlDOF son respectivamente el nodo y desplazamiento de control deseados
    # nodes_control son los nodos de los diafragmas de piso
    # elements son los elementos de los que va a guardar la información
    
    # creación del recorder de techo y definición de la tolerancia
    # recorder('Node','-file','techo.out','-time','-node',IDctrlNode,'-dof',IDctrlDOF,'disp')
    maxNumIter = 25
    
    if IDctrlDOF == 1:
        dir2 = 2
    else:
        dir2 = 1
        
    # creación del pattern
    timeSeries('Path',1000,'-filePath',recordName[0],'-dt',dtrec,'-factor',fact)
    pattern('UniformExcitation',  1000,   IDctrlDOF,  '-accel', 1000)
    timeSeries('Path',1001,'-filePath',recordName[1],'-dt',dtrec,'-factor',fact)
    pattern('UniformExcitation',  1001,   dir2,  '-accel', 1001)
    
    # damping
    nmodes = max(modes)+1
    eigval = eigen(nmodes)
    
    eig1 = eigval[modes[0]]
    eig2 = eigval[modes[1]]
    
    w1 = eig1**0.5
    w2 = eig2**0.5
    
    beta = 2.0*damp/(w1 + w2)
    alfa = 2.0*damp*w1*w2/(w1 + w2)
    
    if Kswitch == 1:
        rayleigh(alfa, 0.0, beta, 0.0)
    else:
        rayleigh(alfa, beta, 0.0, 0.0)
    
    # configuración básica del análisis
    wipeAnalysis()
    constraints('Transformation')
    numberer('RCM')
    system('BandGeneral')
    test('NormDispIncr', Tol, maxNumIter)
    algorithm('Newton')    
    integrator('Newmark', 0.5, 0.25)
    analysis('Transient')
    
    # Otras opciones de análisis    
    tests = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    algoritmo = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}

    # rutina del análisis
    
    Nsteps =  int(dtrec*nPts/dtan)
    
    dtecho1 = [nodeDisp(IDctrlNode,IDctrlDOF)]
    dtecho2 = [nodeDisp(IDctrlNode,dir2)]

    t = [getTime()]
    nels = len(elements)
    nnodos = len(nodes_control)
    Eds = np.zeros((nels, Nsteps+1, 24))
    node_disp = np.zeros((Nsteps + 1, nnodos)) # para grabar los desplazamientos de los nodos
    node_vel = np.zeros((Nsteps + 1, nnodos)) # para grabar los desplazamientos de los nodos
    node_acel = np.zeros((Nsteps + 1, nnodos)) # para grabar los desplazamientos de los nodos
    node_disp2 = np.zeros((Nsteps + 1, nnodos))
    node_acel2 = np.zeros((Nsteps + 1, nnodos))
    
    for k in range(Nsteps):
        ok = analyze(1,dtan)
        # ok2 = ok;
        # En caso de no converger en un paso entra al condicional que sigue
        if ok != 0:
            print('configuración por defecto no converge en desplazamiento: ',nodeDisp(IDctrlNode,IDctrlDOF))
            for j in algoritmo:
                if j < 4:
                    algorithm(algoritmo[j], '-initial')
    
                else:
                    algorithm(algoritmo[j])
                
                # el test se hace 50 veces más
                test('NormDispIncr', Tol, maxNumIter*50)
                ok = analyze(1,dtan)
                if ok == 0:
                    # si converge vuelve a las opciones iniciales de análisi
                    test('NormDispIncr', Tol, maxNumIter)
                    algorithm('Newton')
                    break
                    
        if ok != 0:
            print('Análisis dinámico fallido')
            print('Desplazamiento alcanzado: ',nodeDisp(IDctrlNode,IDctrlDOF),'m')
            break
    
        
        dtecho1.append(nodeDisp(IDctrlNode,IDctrlDOF))
        dtecho2.append(nodeDisp(IDctrlNode,dir2))
        t.append(getTime())
        
        for node_i, node_tag in enumerate(nodes_control):           
            node_disp[k+1,node_i] = nodeDisp(node_tag,IDctrlDOF)
            node_disp2[k+1,node_i] = nodeDisp(node_tag,dir2)
            node_vel[k+1,node_i] = nodeVel(node_tag,IDctrlDOF)
            node_acel[k+1,node_i] = nodeAccel(node_tag,IDctrlDOF)
            node_acel2[k+1,node_i] = nodeAccel(node_tag,dir2)
            
        for el_i, ele_tag in enumerate(elements):
            Eds[el_i , k+1, :] = eleResponse(ele_tag,'globalForce')
        
    
    techo1 = np.array(dtecho1)
    techo2 = np.array(dtecho2)
    techoT = np.sqrt(techo1**2+techo2**2)
    tiempo = np.array(t)
    
    return tiempo,techo1,techo2,techoT,node_disp,node_vel,node_acel,node_disp2,node_acel2,Eds
     
def dinamicoIDA(recordName,dtrec,nPts,dtan,fact,damp,IDctrlNode,IDctrlDOF,modes = [0,2],Kswitch = 1,Tol=1e-4):
    
    # modelName
    # record es el nombre del registro, incluyendo extensión. P.ej. GM01.txt
    # dtrec es el dt del registro
    # nPts es el número de puntos del análisis
    # dtan es el dt del análisis
    # fact es el factor escalar del registro
    # damp es el porcentaje de amortiguamiento (EN DECIMAL. p.ej: 0.03 para 3%)
    # Kswitch recibe: 1: matriz inicial, 2: matriz actual
    # IDctrlNode,IDctrlDOF son respectivamente el nodo y desplazamiento de control deseados
     # creación del recorder de techo y definición de la tolerancia
    # nombre = str(int(fact/9.81*100))
    # recorder('Node','-file','techo'+nombre+'.out','-time','-node',IDctrlNode,'-dof',IDctrlDOF,'disp')
    maxNumIter = 25
    
    # creación del pattern
    
    timeSeries('Path',1000,'-filePath',recordName,'-dt',dtrec,'-factor',fact)
    pattern('UniformExcitation',  1000,   IDctrlDOF,  '-accel', 1000)
    
    # damping
    nmodes = max(modes)+1
    eigval = eigen(nmodes)
    
    eig1 = eigval[modes[0]]
    eig2 = eigval[modes[1]]
    
    w1 = eig1**0.5
    w2 = eig2**0.5
    
    beta = 2.0*damp/(w1 + w2)
    alfa = 2.0*damp*w1*w2/(w1 + w2)
    
    if Kswitch == 1:
        rayleigh(alfa, 0.0, beta, 0.0)
    else:
        rayleigh(alfa, beta, 0.0, 0.0)
    
    # configuración básica del análisis
    wipeAnalysis()
    constraints('Transformation')
    numberer('RCM')
    system('BandGeneral')
    test('NormDispIncr', Tol, maxNumIter)
    algorithm('Newton')    
    integrator('Newmark', 0.5, 0.25)
    analysis('Transient')
    
    # Otras opciones de análisis    
    tests = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    algoritmo = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}

    # rutina del análisis
    
    Nsteps =  int(dtrec*nPts/dtan)
    dtecho = [nodeDisp(IDctrlNode,IDctrlDOF)]
    t = [getTime()]
    
    for k in range(Nsteps):
        ok = analyze(1,dtan)
        # ok2 = ok;
        # En caso de no converger en un paso entra al condicional que sigue
        if ok != 0:
            print('configuración por defecto no converge en tiempo: ',getTime())
            for j in algoritmo:
                if j < 4:
                    algorithm(algoritmo[j], '-initial')
    
                else:
                    algorithm(algoritmo[j])
                
                # el test se hace 50 veces más
                test('NormDispIncr', Tol, maxNumIter*50)
                ok = analyze(1,dtan)
                if ok == 0:
                    # si converge vuelve a las opciones iniciales de análisi
                    test('NormDispIncr', Tol, maxNumIter)
                    algorithm('Newton')
                    break
                    
        if ok != 0:
            print('Análisis dinámico fallido')
            print('Desplazamiento alcanzado: ',nodeDisp(IDctrlNode,IDctrlDOF),'m')
            break
    
        
        dtecho.append(nodeDisp(IDctrlNode,IDctrlDOF))
        t.append(getTime())
        
    # plt.figure()
    # plt.plot(t,dtecho)
    # plt.xlabel('tiempo (s)')
    # plt.ylabel('desplazamiento (m)')
    
    # techo = np.array(dtecho)
    # tiempo = np.array(t)
    wipe()
   
def dinamicoAnim(recordName,dtrec,nPts,dtan,fact,damp,IDctrlNode,IDctrlDOF,modes = [0,2],Kswitch = 1,Tol=1e-8):
    
    # record es el nombre del registro, incluyendo extensión. P.ej. GM01.txt
    # dtrec es el dt del registro
    # nPts es el número de puntos del análisis
    # dtan es el dt del análisis
    # fact es el factor escalar del registro
    # damp es el porcentaje de amortiguamiento (EN DECIMAL. p.ej: 0.03 para 3%)
    # Kswitch recibe: 1: matriz inicial, 2: matriz actual
    # IDctrlNode,IDctrlDOF son respectivamente el nodo y desplazamiento de control deseados
    
    # creación del recorder de techo y definición de la tolerancia
    recorder('Node','-file','techo.out','-time','-node',IDctrlNode,'-dof',IDctrlDOF,'disp')
    maxNumIter = 10
    
    # creación del pattern
    
    timeSeries('Path',1000,'-filePath',recordName,'-dt',dtrec,'-factor',fact)
    pattern('UniformExcitation',  1000,   IDctrlDOF,  '-accel', 1000)
    
    # damping
    nmodes = max(modes)+1
    eigval = eigen(nmodes)
    
    eig1 = eigval[modes[0]]
    eig2 = eigval[modes[1]]
    
    w1 = eig1**0.5
    w2 = eig2**0.5
    
    beta = 2.0*damp/(w1 + w2)
    alfa = 2.0*damp*w1*w2/(w1 + w2)
    
    if Kswitch == 1:
        rayleigh(alfa, 0.0, beta, 0.0)
    else:
        rayleigh(alfa, beta, 0.0, 0.0)
    
    # configuración básica del análisis
    wipeAnalysis()
    constraints('Transformation')
    numberer('RCM')
    system('BandGeneral')
    test('EnergyIncr', Tol, maxNumIter)
    algorithm('Newton')    
    integrator('Newmark', 0.5, 0.25)
    analysis('Transient')
    
    # Otras opciones de análisis    
    tests = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    algoritmo = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}

    # rutina del análisis
    
    Nsteps =  int(dtrec*nPts/dtan)
    dtecho = [nodeDisp(IDctrlNode,IDctrlDOF)]
    t = [getTime()]
    el_tags = getEleTags()
    nels = len(el_tags)
    Eds = np.zeros((Nsteps+1, nels, 6))
    
    
    for k in range(Nsteps):
        ok = analyze(1,dtan)
        # ok2 = ok;
        # En caso de no converger en un paso entra al condicional que sigue
        if ok != 0:
            print('configuración por defecto no converge en desplazamiento: ',nodeDisp(IDctrlNode,IDctrlDOF))
            for j in algoritmo:
                if j < 4:
                    algorithm(algoritmo[j], '-initial')
    
                else:
                    algorithm(algoritmo[j])
                
                # el test se hace 50 veces más
                test('EnergyIncr', Tol, maxNumIter*50)
                ok = analyze(1,dtan)
                if ok == 0:
                    # si converge vuelve a las opciones iniciales de análisi
                    test('EnergyIncr', Tol, maxNumIter)
                    algorithm('Newton')
                    break
                    
        if ok != 0:
            print('Análisis dinámico fallido')
            print('Desplazamiento alcanzado: ',nodeDisp(IDctrlNode,IDctrlDOF),'m')
            break
    
        for el_i, ele_tag in enumerate(el_tags):
            nd1, nd2 = eleNodes(ele_tag)
            Eds[k+1, el_i, :] = [nodeDisp(nd1)[0],
                                  nodeDisp(nd1)[1],
                                  nodeDisp(nd1)[2],
                                  nodeDisp(nd2)[0],
                                  nodeDisp(nd2)[1],
                                  nodeDisp(nd2)[2]]
        
        
        dtecho.append(nodeDisp(IDctrlNode,IDctrlDOF))
        t.append(getTime())
        
    plt.figure()
    plt.plot(t,dtecho)
    plt.xlabel('tiempo (s)')
    plt.ylabel('desplazamiento (m)')  
    
    techo = np.array(dtecho)
    tiempo = np.array(t)
    
    return tiempo,techo,Eds
    
def dinamicoIDA2(recordName,dtrec,nPts,dtan,fact,damp,IDctrlNode,IDctrlDOF,modes = [0,2],Kswitch = 1,Tol=1e-3):
    '''  
    Performs a dynamic analysis recording the displacement of a user selected node.
    Parameters
    ----------
    recordName : string
        Name of the record including file extension (i.e., 'GM01.txt'). It must have one record instant per line. 
    dtrec : float
        time increment of the record.
    nPts : integer
        number of points of the record.
    dtan : float
        time increment to be used in the analysis. If smaller than dtrec, OpenSeesPy interpolates.
    fact : float
        scale factor to apply to the record.
    damp : float
        Damping percentage in decimal (i.e., use 0.03 for 3%).
    IDctrlNode : int
        control node for the displacements.
    IDctrlDOF : int
        DOF for the displacement.
    modes : list, optional
        Modes of the structure to apply the Rayleigh damping. The default is [0,2] which uses the first and third mode.
    Kswitch : int, optional
        Use it to define which stiffness matrix should be used for the ramping. The default is 1 that uses initial stiffness. Input 2 for current stifness.
    Tol : float, optional
        Tolerance for the analysis. The default is 1e-4 because it uses the NormUnbalance test.

    Returns
    -------
    tiempo : numpy array
        Numpy array with analysis time.
    techo : numpy array
        Displacement of the control node.

    '''
    # PARA SER UTILIZADO PARA CORRER EN PARALELO LOS SISMOS
    
    # record es el nombre del registro, incluyendo extensión. P.ej. GM01.txt
    # dtrec es el dt del registro
    # nPts es el número de puntos del análisis
    # dtan es el dt del análisis
    # fact es el factor escalar del registro
    # damp es el porcentaje de amortiguamiento (EN DECIMAL. p.ej: 0.03 para 3%)
    # Kswitch recibe: 1: matriz inicial, 2: matriz actual
    
    maxNumIter = 10
    
    # creación del pattern
    
    timeSeries('Path',1000,'-filePath',recordName,'-dt',dtrec,'-factor',fact)
    pattern('UniformExcitation',  1000,   IDctrlDOF,  '-accel', 1000)
    
    # damping
    nmodes = max(modes)+1
    eigval = eigen(nmodes)
    
    eig1 = eigval[modes[0]]
    eig2 = eigval[modes[1]]
    
    w1 = eig1**0.5
    w2 = eig2**0.5
    
    beta = 2.0*damp/(w1 + w2)
    alfa = 2.0*damp*w1*w2/(w1 + w2)
    
    if Kswitch == 1:
        rayleigh(alfa, 0.0, beta, 0.0)
    else:
        rayleigh(alfa, beta, 0.0, 0.0)
    
    # configuración básica del análisis
    wipeAnalysis()
    constraints('Transformation')
    numberer('RCM')
    system('BandGeneral')
    test('NormUnbalance', Tol, maxNumIter)
    algorithm('Newton')    
    integrator('Newmark', 0.5, 0.25)
    analysis('Transient')
    
    # Otras opciones de análisis    
    tests = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    algoritmo = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}

    # rutina del análisis
    
    Nsteps =  int(dtrec*nPts/dtan)
    dtecho = [nodeDisp(IDctrlNode,IDctrlDOF)]
    t = [getTime()]
    
    for k in range(Nsteps):
        ok = analyze(1,dtan)
        # ok2 = ok;
        # En caso de no converger en un paso entra al condicional que sigue
        if ok != 0:
            print('configuración por defecto no converge en tiempo: ',getTime())
            for j in algoritmo:
                if j < 4:
                    algorithm(algoritmo[j], '-initial')
    
                else:
                    algorithm(algoritmo[j])
                
                # el test se hace 50 veces más
                test('NormUnbalance', Tol, maxNumIter*10)
                ok = analyze(1,dtan)
                if ok == 0:
                    # si converge vuelve a las opciones iniciales de análisi
                    test('NormUnbalance', Tol, maxNumIter)
                    algorithm('Newton')
                    break
                    
        if ok != 0:
            print('Análisis dinámico fallido')
            print('Desplazamiento alcanzado: ',nodeDisp(IDctrlNode,IDctrlDOF),'m')
            break
    
        
        dtecho.append(nodeDisp(IDctrlNode,IDctrlDOF))
        t.append(getTime())
        
    # plt.figure()
    # plt.plot(t,dtecho)
    # plt.xlabel('tiempo (s)')
    # plt.ylabel('desplazamiento (m)')
    
    techo = np.array(dtecho)
    tiempo = np.array(t)
    wipe()
    return tiempo,techo


def dinamicoIDA3(recordName,dtrec,nPts,dtan,fact,damp,IDctrlNode,IDctrlDOF,elements,modes = [0,2],Kswitch = 1,Tol=1e-8):
    
    # PARA SER UTILIZADO PARA CORRER EN PARALELO LOS SISMOS Y EXTRAYENDO LAS FUERZAS DE LOS ELEMENTOS INDICADOS EN ELEMENTS
    
    # record es el nombre del registro, incluyendo extensión. P.ej. GM01.txt
    # dtrec es el dt del registro
    # nPts es el número de puntos del análisis
    # dtan es el dt del análisis
    # fact es el factor escalar del registro
    # damp es el porcentaje de amortiguamiento (EN DECIMAL. p.ej: 0.03 para 3%)
    # IDcrtlNode es el nodo de control para grabar desplazamientos
    # IDctrlDOF es el grado de libertad de control
    # elements son los elementos de los que se va a grabar información
    # nodes_control son los nodos donde se va a grabar las respuestas
    # Kswitch recibe: 1: matriz inicial, 2: matriz actual
    
    maxNumIter = 10
    
    # creación del pattern
    
    timeSeries('Path',1000,'-filePath',recordName,'-dt',dtrec,'-factor',fact)
    pattern('UniformExcitation',  1000,   IDctrlDOF,  '-accel', 1000)
    
    # damping
    nmodes = max(modes)+1
    eigval = eigen(nmodes)
    
    eig1 = eigval[modes[0]]
    eig2 = eigval[modes[1]]
    
    w1 = eig1**0.5
    w2 = eig2**0.5
    
    beta = 2.0*damp/(w1 + w2)
    alfa = 2.0*damp*w1*w2/(w1 + w2)
    
    if Kswitch == 1:
        rayleigh(alfa, 0.0, beta, 0.0)
    else:
        rayleigh(alfa, beta, 0.0, 0.0)
    
    # configuración básica del análisis
    wipeAnalysis()
    constraints('Transformation')
    numberer('RCM')
    system('BandGeneral')
    test('EnergyIncr', Tol, maxNumIter)
    algorithm('Newton')    
    integrator('Newmark', 0.5, 0.25)
    analysis('Transient')
    
    # Otras opciones de análisis    
    tests = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    algoritmo = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}

    # rutina del análisis
    
    Nsteps =  int(dtrec*nPts/dtan)
    dtecho = [nodeDisp(IDctrlNode,IDctrlDOF)]
    t = [getTime()]
    nels = len(elements)
    # nnodos = len(nodes_control)
    Eds = np.zeros((nels, Nsteps+1, 6)) # para grabar las fuerzas de los elementos
    Curv = np.zeros((nels,Nsteps+1)) # para grabar la curvatura de los elementos
    # Strains = np.zeros((Nsteps+1, 8, nels)) # # para grabar las deformaciones de los muros en las 8 fibras que tienen los elementos
    Strains = np.zeros((nels, Nsteps+1, 8))
    cStress = np.zeros((nels, Nsteps+1, 8)) # # para grabar los esfuerzos del concreto de los muros en las 8 fibras que tienen los elementos
    sStress = np.zeros((nels, Nsteps+1, 8)) # # para grabar los esfuerzos del acero de los muros en las 8 fibras que tienen los elementos
    # node_disp = np.zeros((Nsteps + 1, nnodos)) # para grabar los desplazamientos de los nodos
    # node_vel = np.zeros((Nsteps + 1, nnodos)) # para grabar los desplazamientos de los nodos
    # node_acel = np.zeros((Nsteps + 1, nnodos)) # para grabar los desplazamientos de los nodos
    
    for k in range(Nsteps):
        ok = analyze(1,dtan)
        # ok2 = ok;
        # En caso de no converger en un paso entra al condicional que sigue
        if ok != 0:
            print('configuración por defecto no converge en tiempo: ',getTime())
            for j in algoritmo:
                if j < 4:
                    algorithm(algoritmo[j], '-initial')
    
                else:
                    algorithm(algoritmo[j])
                
                # el test se hace 50 veces más
                test('EnergyIncr', Tol, maxNumIter*50)
                ok = analyze(1,dtan)
                if ok == 0:
                    # si converge vuelve a las opciones iniciales de análisi
                    test('EnergyIncr', Tol, maxNumIter)
                    algorithm('Newton')
                    break
                    
        if ok != 0:
            print('Análisis dinámico fallido')
            print('Desplazamiento alcanzado: ',nodeDisp(IDctrlNode,IDctrlDOF),'m')
            break
        
        # for node_i, node_tag in enumerate(nodes_control):
            
        #     node_disp[k+1,node_i] = nodeDisp(node_tag,1)
        #     node_disp[k+1,node_i] = nodeDisp(node_tag,1)
                           
        
        
        
        for el_i, ele_tag in enumerate(elements):
            
            # Curv[k+1, el_i] = [eleResponse(ele_tag,'Curvature')]
            Eds[el_i , k+1, :] = [eleResponse(ele_tag,'globalForce')[0],
                                 eleResponse(ele_tag,'globalForce')[1],
                                 eleResponse(ele_tag,'globalForce')[2],
                                 eleResponse(ele_tag,'globalForce')[3],
                                 eleResponse(ele_tag,'globalForce')[4],
                                 eleResponse(ele_tag,'globalForce')[5]]
            
            Strains[el_i , k+1, :] = [eleResponse(ele_tag,'Fiber_Strain')[0],
                                 eleResponse(ele_tag,'Fiber_Strain')[1],
                                 eleResponse(ele_tag,'Fiber_Strain')[2],
                                 eleResponse(ele_tag,'Fiber_Strain')[3],
                                 eleResponse(ele_tag,'Fiber_Strain')[4],
                                 eleResponse(ele_tag,'Fiber_Strain')[5],
                                 eleResponse(ele_tag,'Fiber_Strain')[6],
                                 eleResponse(ele_tag,'Fiber_Strain')[7]]
            
            cStress[el_i , k+1, :] = [eleResponse(ele_tag,'Fiber_Stress_Concrete')[0],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[1],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[2],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[3],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[4],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[5],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[6],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[7]]
            
            sStress[el_i , k+1, :] = [eleResponse(ele_tag,'Fiber_Stress_Steel')[0],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[1],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[2],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[3],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[4],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[5],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[6],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[7]]
            
            
            
        dtecho.append(nodeDisp(IDctrlNode,IDctrlDOF))
        t.append(getTime())
        
    # plt.figure()
    # plt.plot(t,dtecho)
    # plt.xlabel('tiempo (s)')
    # plt.ylabel('desplazamiento (m)')
    
    techo = np.array(dtecho)
    tiempo = np.array(t)
    wipe()
    return tiempo,techo,Eds,Strains,cStress,sStress


def dinamicoIDA4(recordName,dtrec,nPts,dtan,fact,damp,IDctrlNode,IDctrlDOF,elements,nodes_control,modes = [0,2],Kswitch = 1,Tol=1e-8):
    
    # PARA SER UTILIZADO PARA CORRER EN PARALELO LOS SISMOS Y EXTRAYENDO LAS FUERZAS DE LOS ELEMENTOS INDICADOS EN ELEMENTS
    
    # record es el nombre del registro, incluyendo extensión. P.ej. GM01.txt
    # dtrec es el dt del registro
    # nPts es el número de puntos del análisis
    # dtan es el dt del análisis
    # fact es el factor escalar del registro
    # damp es el porcentaje de amortiguamiento (EN DECIMAL. p.ej: 0.03 para 3%)
    # IDcrtlNode es el nodo de control para grabar desplazamientos
    # IDctrlDOF es el grado de libertad de control
    # elements son los elementos de los que se va a grabar información
    # nodes_control son los nodos donde se va a grabar las respuestas
    # Kswitch recibe: 1: matriz inicial, 2: matriz actual
    
    maxNumIter = 10
    
    # creación del pattern
    
    timeSeries('Path',1000,'-filePath',recordName,'-dt',dtrec,'-factor',fact)
    pattern('UniformExcitation',  1000,   IDctrlDOF,  '-accel', 1000)
    
    # damping
    nmodes = max(modes)+1
    eigval = eigen(nmodes)
    
    eig1 = eigval[modes[0]]
    eig2 = eigval[modes[1]]
    
    w1 = eig1**0.5
    w2 = eig2**0.5
    
    beta = 2.0*damp/(w1 + w2)
    alfa = 2.0*damp*w1*w2/(w1 + w2)
    
    if Kswitch == 1:
        rayleigh(alfa, 0.0, beta, 0.0)
    else:
        rayleigh(alfa, beta, 0.0, 0.0)
    
    # configuración básica del análisis
    wipeAnalysis()
    constraints('Transformation')
    numberer('RCM')
    system('BandGeneral')
    test('EnergyIncr', Tol, maxNumIter)
    algorithm('Newton')    
    integrator('Newmark', 0.5, 0.25)
    analysis('Transient')
    
    # Otras opciones de análisis    
    tests = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    algoritmo = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}

    # rutina del análisis
    
    Nsteps =  int(dtrec*nPts/dtan)
    dtecho = [nodeDisp(IDctrlNode,IDctrlDOF)]
    t = [getTime()]
    nels = len(elements)
    nnodos = len(nodes_control)
    Eds = np.zeros((nels, Nsteps+1, 24)) # para grabar las fuerzas de los elementos
    Curv = np.zeros((nels,Nsteps+1)) # para grabar la curvatura de los elementos
    
    Strains = np.zeros((nels, Nsteps+1, 8))
    cStress = np.zeros((nels, Nsteps+1, 8)) # # para grabar los esfuerzos del concreto de los muros en las 8 fibras que tienen los elementos
    sStress = np.zeros((nels, Nsteps+1, 8)) # # para grabar los esfuerzos del acero de los muros en las 8 fibras que tienen los elementos
    node_disp = np.zeros((Nsteps + 1, nnodos)) # para grabar los desplazamientos de los nodos
    node_vel = np.zeros((Nsteps + 1, nnodos)) # para grabar los desplazamientos de los nodos
    node_acel = np.zeros((Nsteps + 1, nnodos)) # para grabar los desplazamientos de los nodos
    drift = np.zeros((Nsteps + 1, nnodos - 1)) # para grabar la deriva de entrepiso
    
    for k in range(Nsteps):
        ok = analyze(1,dtan)
        # ok2 = ok;
        # En caso de no converger en un paso entra al condicional que sigue
        if ok != 0:
            print('configuración por defecto no converge en tiempo: ',getTime())
            for j in algoritmo:
                if j < 4:
                    algorithm(algoritmo[j], '-initial')
    
                else:
                    algorithm(algoritmo[j])
                
                # el test se hace 50 veces más
                test('EnergyIncr', Tol, maxNumIter*50)
                ok = analyze(1,dtan)
                if ok == 0:
                    # si converge vuelve a las opciones iniciales de análisi
                    test('EnergyIncr', Tol, maxNumIter)
                    algorithm('Newton')
                    break
                    
        if ok != 0:
            print('Análisis dinámico fallido')
            print('Desplazamiento alcanzado: ',nodeDisp(IDctrlNode,IDctrlDOF),'m')
            break
        
        for node_i, node_tag in enumerate(nodes_control):
            
            node_disp[k+1,node_i] = nodeDisp(node_tag,1)
            node_vel[k+1,node_i] = nodeVel(node_tag,1)
            node_acel[k+1,node_i] = nodeAccel(node_tag,1)
            if node_i != 0:
                drift[k+1,node_i-1] = (nodeDisp(node_tag,1) - nodeDisp(nodes_control[node_i-1],1))/(nodeCoord(node_tag,2) - nodeCoord(nodes_control[node_i-1],2))
                       

        for el_i, ele_tag in enumerate(elements):
            
            # Curv[k+1, el_i] = [eleResponse(ele_tag,'Curvature')]
            Eds[el_i , k+1, :] = [eleResponse(ele_tag,'globalForce')[0],
                                 eleResponse(ele_tag,'globalForce')[1],
                                 eleResponse(ele_tag,'globalForce')[2],
                                 eleResponse(ele_tag,'globalForce')[3],
                                 eleResponse(ele_tag,'globalForce')[4],
                                 eleResponse(ele_tag,'globalForce')[5]]
            
            Strains[el_i , k+1, :] = [eleResponse(ele_tag,'Fiber_Strain')[0],
                                 eleResponse(ele_tag,'Fiber_Strain')[1],
                                 eleResponse(ele_tag,'Fiber_Strain')[2],
                                 eleResponse(ele_tag,'Fiber_Strain')[3],
                                 eleResponse(ele_tag,'Fiber_Strain')[4],
                                 eleResponse(ele_tag,'Fiber_Strain')[5],
                                 eleResponse(ele_tag,'Fiber_Strain')[6],
                                 eleResponse(ele_tag,'Fiber_Strain')[7]]
            
            cStress[el_i , k+1, :] = [eleResponse(ele_tag,'Fiber_Stress_Concrete')[0],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[1],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[2],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[3],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[4],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[5],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[6],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[7]]
            
            sStress[el_i , k+1, :] = [eleResponse(ele_tag,'Fiber_Stress_Steel')[0],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[1],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[2],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[3],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[4],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[5],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[6],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[7]]
            
            
            
        dtecho.append(nodeDisp(IDctrlNode,IDctrlDOF))
        t.append(getTime())
        
    # plt.figure()
    # plt.plot(t,dtecho)
    # plt.xlabel('tiempo (s)')
    # plt.ylabel('desplazamiento (m)')
    
    techo = np.array(dtecho)
    tiempo = np.array(t)
    wipe()
    return tiempo,techo,Eds,Strains,cStress,sStress,node_disp,node_vel,node_acel,drift


def dinamicoIDA4P(recordName,dtrec,nPts,dtan,fact,damp,IDctrlNode,IDctrlDOF,elements,nodes_control,modes = [0,2],Kswitch = 1,Tol=1e-4):
    '''
    Performs a dynamic analysis for a ground motion, recording information about displacements, velocity, accelerations, forces. Only allows elements with six DOF per node.

    Parameters
    ----------
    recordName : string
        Name of the record including file extension (i.e., 'GM01.txt'). It must have one record instant per line. 
    dtrec : float
        time increment of the record.
    nPts : integer
        number of points of the record.
    dtan : float
        time increment to be used in the analysis. If smaller than dtrec, OpenSeesPy interpolates.
    fact : float
        scale factor to apply to the record.
    damp : float
        Damping percentage in decimal (i.e., use 0.03 for 3%).
    IDctrlNode : int
        control node for the displacements.
    IDctrlDOF : int
        DOF for the displacement.
    elements : list
        elements to record forces and stresses.
    nodes_control : list
        nodes to compute displacements and inter-story drift. You must input one per floor, otherwise you'll get an error.
    modes : list, optional
        Modes of the structure to apply the Rayleigh damping. The default is [0,2] which uses the first and third mode.
    Kswitch : int, optional
        Use it to define which stiffness matrix should be used for the ramping. The default is 1 that uses initial stiffness. Input 2 for current stifness.
    Tol : float, optional
        Tolerance for the analysis. The default is 1e-4 because it uses the NormUnbalance test.

    Returns
    -------
    tiempo : numpy array
        Numpy array with analysis time.
    techo : numpy array
        Displacement of the control node.
    Eds :
        Numpy array with the forces in the elements (columns and beams). The order is determined by the order used in the input variable elements. The array has three dimensions. The first one is the element, the second one the pushover instant and the third one is the DOF.
    node_disp : numpy array
        Displacement at each node in nodes_control. Each column correspond to a node and each row to an analysis instant.
    node_vel : numpy array
        Velocity at each node in nodes_control. Each column correspond to a node and each row to an analysis instant.
    node_acel : numpy array
        Relative displacement at each node in nodes_control. Each column correspond to a node and each row to an analysis instant.
    drift : numpy array
        Drift at story of the building. Each column correspond to a node and each row to an analysis instant.

    '''
    # PARA SER UTILIZADO PARA CORRER EN PARALELO LOS SISMOS Y EXTRAYENDO LAS FUERZAS DE LOS ELEMENTOS INDICADOS EN ELEMENTS
    
    # record es el nombre del registro, incluyendo extensión. P.ej. GM01.txt
    # dtrec es el dt del registro
    # nPts es el número de puntos del análisis
    # dtan es el dt del análisis
    # fact es el factor escalar del registro
    # damp es el porcentaje de amortiguamiento (EN DECIMAL. p.ej: 0.03 para 3%)
    # IDcrtlNode es el nodo de control para grabar desplazamientos
    # IDctrlDOF es el grado de libertad de control
    # elements son los elementos de los que se va a grabar información
    # nodes_control son los nodos donde se va a grabar las respuestas
    # Kswitch recibe: 1: matriz inicial, 2: matriz actual
    
    maxNumIter = 10
    
    # creación del pattern
    
    timeSeries('Path',1000,'-filePath',recordName,'-dt',dtrec,'-factor',fact)
    pattern('UniformExcitation',  1000,   IDctrlDOF,  '-accel', 1000)
    
    # damping
    nmodes = max(modes)+1
    eigval = eigen(nmodes)
    
    eig1 = eigval[modes[0]]
    eig2 = eigval[modes[1]]
    
    w1 = eig1**0.5
    w2 = eig2**0.5
    
    beta = 2.0*damp/(w1 + w2)
    alfa = 2.0*damp*w1*w2/(w1 + w2)
    
    if Kswitch == 1:
        rayleigh(alfa, 0.0, beta, 0.0)
    else:
        rayleigh(alfa, beta, 0.0, 0.0)
    
    # configuración básica del análisis
    wipeAnalysis()
    constraints('Transformation')
    numberer('RCM')
    system('BandGeneral')
    test('NormUnbalance', Tol, maxNumIter)
    algorithm('Newton')    
    integrator('Newmark', 0.5, 0.25)
    analysis('Transient')
    
    # Otras opciones de análisis    
    tests = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    algoritmo = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}

    # rutina del análisis
    
    Nsteps =  int(dtrec*nPts/dtan)
    dtecho = [nodeDisp(IDctrlNode,IDctrlDOF)]
    t = [getTime()]
    nels = len(elements)
    nnodos = len(nodes_control)
    Eds = np.zeros((nels, Nsteps+1, 12)) # para grabar las fuerzas de los elementos
    
    
    
    node_disp = np.zeros((Nsteps + 1, nnodos)) # para grabar los desplazamientos de los nodos
    node_vel = np.zeros((Nsteps + 1, nnodos)) # para grabar los desplazamientos de los nodos
    node_acel = np.zeros((Nsteps + 1, nnodos)) # para grabar los desplazamientos de los nodos
    drift = np.zeros((Nsteps + 1, nnodos - 1)) # para grabar la deriva de entrepiso
    
    for k in range(Nsteps):
        ok = analyze(1,dtan)
        # ok2 = ok;
        # En caso de no converger en un paso entra al condicional que sigue
        if ok != 0:
            print('configuración por defecto no converge en tiempo: ',getTime())
            for j in algoritmo:
                if j < 4:
                    algorithm(algoritmo[j], '-initial')
    
                else:
                    algorithm(algoritmo[j])
                
                # el test se hace 50 veces más
                test('NormUnbalance', Tol, maxNumIter*50)
                ok = analyze(1,dtan)
                if ok == 0:
                    # si converge vuelve a las opciones iniciales de análisi
                    test('NormUnbalance', Tol, maxNumIter)
                    algorithm('Newton')
                    break
                    
        if ok != 0:
            print('Análisis dinámico fallido')
            print('Desplazamiento alcanzado: ',nodeDisp(IDctrlNode,IDctrlDOF),'m')
            break
        
        for node_i, node_tag in enumerate(nodes_control):
            
            node_disp[k+1,node_i] = nodeDisp(node_tag,IDctrlDOF)
            node_vel[k+1,node_i] = nodeVel(node_tag,IDctrlDOF)
            node_acel[k+1,node_i] = nodeAccel(node_tag,IDctrlDOF)
            if node_i != 0:
                drift[k+1,node_i-1] = (nodeDisp(node_tag,IDctrlDOF) - nodeDisp(nodes_control[node_i-1],IDctrlDOF))/(nodeCoord(node_tag,3) - nodeCoord(nodes_control[node_i-1],3))
                       

        for el_i, ele_tag in enumerate(elements):
                      
            Eds[el_i , k+1, :] = eleResponse(ele_tag,'globalForce')
            
        dtecho.append(nodeDisp(IDctrlNode,IDctrlDOF))
        t.append(getTime())
        
    # plt.figure()
    # plt.plot(t,dtecho)
    # plt.xlabel('tiempo (s)')
    # plt.ylabel('desplazamiento (m)')
    
    techo = np.array(dtecho)
    tiempo = np.array(t)
    wipe()
    return tiempo,techo,Eds,node_disp,node_vel,node_acel,drift


def dinamicoIDA4G(recordName,dtrec,nPts,dtan,fact,damp,IDctrlNode,IDctrlDOF,nodes_control,Tol=1e-3,modes = [0,2],Kswitch = 1):
    
    # PARA SER UTILIZADO PARA CORRER EN PARALELO LOS SISMOS Y EXTRAYENDO INFORMACIÓN GLOBAL
    
    # record es el nombre del registro, incluyendo extensión. P.ej. GM01.txt
    # dtrec es el dt del registro
    # nPts es el número de puntos del análisis
    # dtan es el dt del análisis
    # fact es el factor escalar del registro
    # damp es el porcentaje de amortiguamiento (EN DECIMAL. p.ej: 0.03 para 3%)
    # IDcrtlNode es el nodo de control para grabar desplazamientos
    # IDctrlDOF es el grado de libertad de control
    # nodes_control son los nodos donde se va a grabar las respuestas
    # Kswitch recibe: 1: matriz inicial, 2: matriz actual
    
    maxNumIter = 10
    
    # creación del pattern
    
    timeSeries('Path',1000,'-filePath',recordName,'-dt',dtrec,'-factor',fact)
    pattern('UniformExcitation',  1000,   IDctrlDOF,  '-accel', 1000)
    
    # damping
    nmodes = max(modes)+1
    eigval = eigen('-genBandArpack',nmodes)
    
    eig1 = eigval[modes[0]]
    eig2 = eigval[modes[1]]
    
    w1 = eig1**0.5
    w2 = eig2**0.5
    
    beta = 2.0*damp/(w1 + w2)
    alfa = 2.0*damp*w1*w2/(w1 + w2)
    
    if Kswitch == 1:
        rayleigh(alfa, 0.0, beta, 0.0)
    else:
        rayleigh(alfa, beta, 0.0, 0.0)
    
    # configuración básica del análisis
    wipeAnalysis()
    constraints('Transformation')
    numberer('RCM')
    system('BandGeneral')
    test('NormDispIncr', Tol, maxNumIter)
    algorithm('Newton')    
    integrator('Newmark', 0.5, 0.25)
    analysis('Transient')
    
    # Otras opciones de análisis    
    tests = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    algoritmo = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}

    # rutina del análisis
    
    Nsteps =  int(dtrec*nPts/dtan)
    dtecho = [nodeDisp(IDctrlNode,IDctrlDOF)]
    t = [getTime()]
    nnodos = len(nodes_control)
    
    
    
    node_disp = np.zeros((Nsteps + 1, nnodos)) # para grabar los desplazamientos de los nodos
    node_vel = np.zeros((Nsteps + 1, nnodos)) # para grabar los desplazamientos de los nodos
    node_acel = np.zeros((Nsteps + 1, nnodos)) # para grabar los desplazamientos de los nodos
    drift = np.zeros((Nsteps + 1, nnodos - 1)) # para grabar la deriva de entrepiso
    
    for k in range(Nsteps):
        ok = analyze(1,dtan)
        # ok2 = ok;
        # En caso de no converger en un paso entra al condicional que sigue
        if ok != 0:
            print('configuración por defecto no converge en tiempo: ',getTime())
            for j in algoritmo:
                if j < 4:
                    algorithm(algoritmo[j], '-initial')
    
                else:
                    algorithm(algoritmo[j])
                
                # el test se hace 50 veces más
                test('NormDispIncr', Tol, maxNumIter*10)
                ok = analyze(1,dtan)
                if ok == 0:
                    # si converge vuelve a las opciones iniciales de análisi
                    test('NormDispIncr', Tol, maxNumIter)
                    algorithm('Newton')
                    break
                    
        if ok != 0:
            print('Análisis dinámico fallido')
            print('Desplazamiento alcanzado: ',nodeDisp(IDctrlNode,IDctrlDOF),'m')
            break
        
        for node_i, node_tag in enumerate(nodes_control):
            
            node_disp[k+1,node_i] = nodeDisp(node_tag,IDctrlDOF)
            node_vel[k+1,node_i] = nodeVel(node_tag,IDctrlDOF)
            node_acel[k+1,node_i] = nodeAccel(node_tag,IDctrlDOF)
            if node_i != 0:
                drift[k+1,node_i-1] = (nodeDisp(node_tag,IDctrlDOF) - nodeDisp(nodes_control[node_i-1],IDctrlDOF))/(nodeCoord(node_tag,3) - nodeCoord(nodes_control[node_i-1],3))
                                
        dtecho.append(nodeDisp(IDctrlNode,IDctrlDOF))
        t.append(getTime())
        
    # plt.figure()
    # plt.plot(t,dtecho)
    # plt.xlabel('tiempo (s)')
    # plt.ylabel('desplazamiento (m)')
    
    techo = np.array(dtecho)
    tiempo = np.array(t)
    wipe()
    return tiempo,techo,node_disp,node_vel,node_acel,drift

def dinamicoIDA5(acceleration,dtrec,nPts,dtan,fact,damp,IDctrlNode,IDctrlDOF,elements,nodes_control,modes = [0,2],Kswitch = 1,Tol=1e-8):
    
    # PARA SER UTILIZADO PARA CORRER EN PARALELO LOS SISMOS Y EXTRAYENDO LAS FUERZAS DE LOS ELEMENTOS INDICADOS EN ELEMENTS
    
    # acceleration es la lista de aceleraciones del registro. Se termina multiplicando por fact
    # dtrec es el dt del registro
    # nPts es el número de puntos del análisis
    # dtan es el dt del análisis
    # fact es el factor escalar del registro
    # damp es el porcentaje de amortiguamiento (EN DECIMAL. p.ej: 0.03 para 3%)
    # IDcrtlNode es el nodo de control para grabar desplazamientos
    # IDctrlDOF es el grado de libertad de control
    # elements son los elementos de los que se va a grabar información
    # nodes_control son los nodos donde se va a grabar las respuestas
    # Kswitch recibe: 1: matriz inicial, 2: matriz actual
    
    maxNumIter = 10
    
    # creación del pattern
    
    timeSeries('Path',1000,'-values',*acceleration,'-dt',dtrec,'-factor',fact)
    pattern('UniformExcitation',  1000,   IDctrlDOF,  '-accel', 1000)
    
    # damping
    nmodes = max(modes)+1
    eigval = eigen(nmodes)
    
    eig1 = eigval[modes[0]]
    eig2 = eigval[modes[1]]
    
    w1 = eig1**0.5
    w2 = eig2**0.5
    
    beta = 2.0*damp/(w1 + w2)
    alfa = 2.0*damp*w1*w2/(w1 + w2)
    
    if Kswitch == 1:
        rayleigh(alfa, 0.0, beta, 0.0)
    else:
        rayleigh(alfa, beta, 0.0, 0.0)
    
    # configuración básica del análisis
    wipeAnalysis()
    constraints('Transformation')
    numberer('RCM')
    system('BandGeneral')
    test('EnergyIncr', Tol, maxNumIter)
    algorithm('Newton')    
    integrator('Newmark', 0.5, 0.25)
    analysis('Transient')
    
    # Otras opciones de análisis    
    tests = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    algoritmo = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}

    # rutina del análisis
    
    Nsteps =  int(dtrec*nPts/dtan)
    dtecho = [nodeDisp(IDctrlNode,IDctrlDOF)]
    t = [getTime()]
    nels = len(elements)
    nnodos = len(nodes_control)
    Eds = np.zeros((nels, Nsteps+1, 6)) # para grabar las fuerzas de los elementos
    Curv = np.zeros((nels,Nsteps+1)) # para grabar la curvatura de los elementos
    
    Strains = np.zeros((nels, Nsteps+1, 8))
    cStress = np.zeros((nels, Nsteps+1, 8)) # # para grabar los esfuerzos del concreto de los muros en las 8 fibras que tienen los elementos
    sStress = np.zeros((nels, Nsteps+1, 8)) # # para grabar los esfuerzos del acero de los muros en las 8 fibras que tienen los elementos
    node_disp = np.zeros((Nsteps + 1, nnodos)) # para grabar los desplazamientos de los nodos
    node_vel = np.zeros((Nsteps + 1, nnodos)) # para grabar los desplazamientos de los nodos
    node_acel = np.zeros((Nsteps + 1, nnodos)) # para grabar los desplazamientos de los nodos
    drift = np.zeros((Nsteps + 1, nnodos - 1)) # para grabar la deriva de entrepiso
    
    for k in range(Nsteps):
        ok = analyze(1,dtan)
        # ok2 = ok;
        # En caso de no converger en un paso entra al condicional que sigue
        if ok != 0:
            print('configuración por defecto no converge en tiempo: ',getTime())
            for j in algoritmo:
                if j < 4:
                    algorithm(algoritmo[j], '-initial')
    
                else:
                    algorithm(algoritmo[j])
                
                # el test se hace 50 veces más
                test('EnergyIncr', Tol, maxNumIter*50)
                ok = analyze(1,dtan)
                if ok == 0:
                    # si converge vuelve a las opciones iniciales de análisi
                    test('EnergyIncr', Tol, maxNumIter)
                    algorithm('Newton')
                    break
                    
        if ok != 0:
            print('Análisis dinámico fallido')
            print('Desplazamiento alcanzado: ',nodeDisp(IDctrlNode,IDctrlDOF),'m')
            break
        
        for node_i, node_tag in enumerate(nodes_control):
            
            node_disp[k+1,node_i] = nodeDisp(node_tag,1)
            node_vel[k+1,node_i] = nodeVel(node_tag,1)
            node_acel[k+1,node_i] = nodeAccel(node_tag,1)
            if node_i != 0:
                drift[k+1,node_i-1] = (nodeDisp(node_tag,1) - nodeDisp(nodes_control[node_i-1],1))/(nodeCoord(node_tag,2) - nodeCoord(nodes_control[node_i-1],2))
                       

        for el_i, ele_tag in enumerate(elements):
            
            # Curv[k+1, el_i] = [eleResponse(ele_tag,'Curvature')]
            Eds[el_i , k+1, :] = [eleResponse(ele_tag,'globalForce')[0],
                                 eleResponse(ele_tag,'globalForce')[1],
                                 eleResponse(ele_tag,'globalForce')[2],
                                 eleResponse(ele_tag,'globalForce')[3],
                                 eleResponse(ele_tag,'globalForce')[4],
                                 eleResponse(ele_tag,'globalForce')[5]]
            
            Strains[el_i , k+1, :] = [eleResponse(ele_tag,'Fiber_Strain')[0],
                                 eleResponse(ele_tag,'Fiber_Strain')[1],
                                 eleResponse(ele_tag,'Fiber_Strain')[2],
                                 eleResponse(ele_tag,'Fiber_Strain')[3],
                                 eleResponse(ele_tag,'Fiber_Strain')[4],
                                 eleResponse(ele_tag,'Fiber_Strain')[5],
                                 eleResponse(ele_tag,'Fiber_Strain')[6],
                                 eleResponse(ele_tag,'Fiber_Strain')[7]]
            
            cStress[el_i , k+1, :] = [eleResponse(ele_tag,'Fiber_Stress_Concrete')[0],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[1],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[2],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[3],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[4],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[5],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[6],
                                 eleResponse(ele_tag,'Fiber_Stress_Concrete')[7]]
            
            sStress[el_i , k+1, :] = [eleResponse(ele_tag,'Fiber_Stress_Steel')[0],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[1],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[2],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[3],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[4],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[5],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[6],
                                 eleResponse(ele_tag,'Fiber_Stress_Steel')[7]]
            
            
            
        dtecho.append(nodeDisp(IDctrlNode,IDctrlDOF))
        t.append(getTime())
        
    # plt.figure()
    # plt.plot(t,dtecho)
    # plt.xlabel('tiempo (s)')
    # plt.ylabel('desplazamiento (m)')
    
    techo = np.array(dtecho)
    tiempo = np.array(t)
    wipe()
    return tiempo,techo,Eds,Strains,cStress,sStress,node_disp,node_vel,node_acel,drift

def dinamicoIDA6(recordName,dtrec,nPts,dtan,fact,damp,IDctrlNode,IDctrlDOF,elements,nodes_control,modes = [0,2],Kswitch = 1,Tol=1e-8):
    
    # PARA SER UTILIZADO PARA CORRER EN PARALELO LOS SISMOS Y EXTRAYENDO LAS FUERZAS DE LOS ELEMENTOS INDICADOS EN ELEMENTS
    
    # recordName recibe ima lista con la pareja de sismos
    # dtrec es el dt del registro
    # nPts es el número de puntos del análisis
    # dtan es el dt del análisis
    # fact es el factor escalar del registro
    # damp es el porcentaje de amortiguamiento (EN DECIMAL. p.ej: 0.03 para 3%)
    # IDcrtlNode es el nodo de control para grabar desplazamientos
    # IDctrlDOF es el grado de libertad de control
    # elements son los elementos de los que se va a grabar información
    # nodes_control son los nodos donde se va a grabar las respuestas
    # Kswitch recibe: 1: matriz inicial, 2: matriz actual
    
    maxNumIter = 10
    
    # creación del pattern
    
    timeSeries('Path',1000,'-filePath',recordName[0],'-dt',dtrec,'-factor',fact)
    pattern('UniformExcitation',  1000,   1,  '-accel', 1000)
    timeSeries('Path',1001,'-filePath',recordName[1],'-dt',dtrec,'-factor',fact)
    pattern('UniformExcitation',  1001,   2,  '-accel', 1001)
    
    # damping
    nmodes = max(modes)+1
    eigval = eigen(nmodes)
    
    eig1 = eigval[modes[0]]
    eig2 = eigval[modes[1]]
    
    w1 = eig1**0.5
    w2 = eig2**0.5
    
    beta = 2.0*damp/(w1 + w2)
    alfa = 2.0*damp*w1*w2/(w1 + w2)
    
    if Kswitch == 1:
        rayleigh(alfa, 0.0, beta, 0.0)
    else:
        rayleigh(alfa, beta, 0.0, 0.0)
    
    # configuración básica del análisis
    wipeAnalysis()
    constraints('Transformation')
    numberer('RCM')
    system('BandGeneral')
    test('EnergyIncr', Tol, maxNumIter)
    algorithm('Newton')    
    integrator('Newmark', 0.5, 0.25)
    analysis('Transient')
    
    # Otras opciones de análisis    
    tests = {1:'NormDispIncr', 2: 'RelativeEnergyIncr', 4: 'RelativeNormUnbalance',5: 'RelativeNormDispIncr', 6: 'NormUnbalance'}
    algoritmo = {1:'KrylovNewton', 2: 'SecantNewton' , 4: 'RaphsonNewton',5: 'PeriodicNewton', 6: 'BFGS', 7: 'Broyden', 8: 'NewtonLineSearch'}

    # rutina del análisis
    
    Nsteps =  int(dtrec*nPts/dtan)
    dtecho = [nodeDisp(IDctrlNode,IDctrlDOF)]
    t = [getTime()]
    nels = len(elements)
    nnodos = len(nodes_control)
    Eds = np.zeros((nels, Nsteps+1, 6)) # para grabar las fuerzas de los elementos
    
    
    
    node_disp = np.zeros((Nsteps + 1, nnodos)) # para grabar los desplazamientos de los nodos
    node_vel = np.zeros((Nsteps + 1, nnodos)) # para grabar los desplazamientos de los nodos
    node_acel = np.zeros((Nsteps + 1, nnodos)) # para grabar los desplazamientos de los nodos
    drift = np.zeros((Nsteps + 1, nnodos - 1)) # para grabar la deriva de entrepiso
    
    for k in range(Nsteps):
        ok = analyze(1,dtan)
        # ok2 = ok;
        # En caso de no converger en un paso entra al condicional que sigue
        if ok != 0:
            print('configuración por defecto no converge en tiempo: ',getTime())
            for j in algoritmo:
                if j < 4:
                    algorithm(algoritmo[j], '-initial')
    
                else:
                    algorithm(algoritmo[j])
                
                # el test se hace 50 veces más
                test('EnergyIncr', Tol, maxNumIter*50)
                ok = analyze(1,dtan)
                if ok == 0:
                    # si converge vuelve a las opciones iniciales de análisi
                    test('EnergyIncr', Tol, maxNumIter)
                    algorithm('Newton')
                    break
                    
        if ok != 0:
            print('Análisis dinámico fallido')
            print('Desplazamiento alcanzado: ',nodeDisp(IDctrlNode,IDctrlDOF),'m')
            break
        
        for node_i, node_tag in enumerate(nodes_control):
            
            node_disp[k+1,node_i] = nodeDisp(node_tag,1)
            node_vel[k+1,node_i] = nodeVel(node_tag,1)
            node_acel[k+1,node_i] = nodeAccel(node_tag,1)
            if node_i != 0:
                drift[k+1,node_i-1] = (nodeDisp(node_tag,1) - nodeDisp(nodes_control[node_i-1],1))/(nodeCoord(node_tag,2) - nodeCoord(nodes_control[node_i-1],2))
                       

        for el_i, ele_tag in enumerate(elements):
                      
            Eds[el_i , k+1, :] = [eleResponse(ele_tag,'globalForce')[0],
                                 eleResponse(ele_tag,'globalForce')[1],
                                 eleResponse(ele_tag,'globalForce')[2],
                                 eleResponse(ele_tag,'globalForce')[3],
                                 eleResponse(ele_tag,'globalForce')[4],
                                 eleResponse(ele_tag,'globalForce')[5]]
            
        dtecho.append(nodeDisp(IDctrlNode,IDctrlDOF))
        t.append(getTime())
        
    # plt.figure()
    # plt.plot(t,dtecho)
    # plt.xlabel('tiempo (s)')
    # plt.ylabel('desplazamiento (m)')
    
    techo = np.array(dtecho)
    tiempo = np.array(t)
    wipe()
    return tiempo,techo,Eds,node_disp,node_vel,node_acel,drift