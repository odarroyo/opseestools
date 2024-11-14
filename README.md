# opseestools

This is a library to be used with OpenSeesPy.

It includes a collection of analysis functions that are useful if you are analyzing frames or wall buildings both in two and three dimensions, storing the results in Python variables. 
It also includes another useful tools for processing results and creating scripts

## Installation

opseestools can be installed from the Python Package Index as:

#### pip install opseestools

## Description

opseestools comprises a set of functions in three modules:

1) 2D Analysis functions. This library is called analysis. You can import it as:

import opstools.analisis as an

2) 3D Analysis functions. This library is called analysis. You can import it as:

import opstools.analisis3D as an3D (or any name you want)

3) A set of utilities to support model building. You can import it as:

import opstools.utilidades as ut (or any other name)

4) A library to calculate and plot fragility functions. You can import it as:

import opstools.Lib_frag as ut (or any other name)

## Examples of use

#### You can find Jupyter Notebooks with examples of using the library on the examples folder in this GitHub repository. The list of examples are the following:

Example 1: how to install opseestools

Example 2: using functions of the utility module to calculate the spectrum, fourier spectrum, cummulative arias intensity, average pseudoacceleration spectrum of a ground motion

Example 3: using the analysis module to perform a pushover analysis of a RC frame

Example 4: using the analysis module to perform a dynamic analysis of a RC frame

Example 5: using the analysis module to perform an IDA of a RC frame

Example 6: using the analysis module and the lib_frag to perform an multiple record IDA and compute the fragility curves for a RC wall building

Example 7: demonstrating the functions added to the utility module in v0.4 that streamline creating reinforced concrete and steel frame models.

Example 8: demonstrating the functions added to the utility module in v0.45 that further streamline creating reinforced concrete and steel frame models.

Example 9: moment curvature of a T-shaped wall

Example 10: example demonstrating the use of new commands in v0.6. These allow generating 3D RC Frame buildings regular in plan and height. The example creates a 3 story building and runs a pushover in under 80 lines of code.

Example 11: example demonstrating the use of new commands in v0.6. These allow generating 2D RC Frame buildings. The example creates a 2 story building and runs a pushover in under 40 lines of code.

## Function help

Each function has a docstring with help for the function, detailing the inputs and outputs. 

You can summon this help by using the integrated Python help(), for example: help(an.gravedad) will return help for the gravedad function in the analysis module.

## Overview of some functions

You may find useful two functions that perform two of the most common analysis: pushover2 and dinamicoIDA4P:

### pushover2(Dmax,Dincr,IDctrlNode,IDctrlDOF)

This function performs a pushover analysis.
Dmax = maximum displacement
Dincr = increment of displacement
IDctrlNode = control node for the pushover
IDctrlDOF = DOF for pushover

The function returns two variables as numpy arrays: roof displacement and base shear. In that order. An example call  of this function is:

roof,shear = an.pushover2(1.5,0.01,9,1), which performs a pushover analysis using node 9 as a control node in the DOF = 1 until a displacement of 1.5 in 0.01 increments

### dinamicoIDA4P(recordName,dtrec,nPts,dtan,fact,damp,IDctrlNode,IDctrlDOF,elements,nodes_control,modes = [0,2],Kswitch = 1,Tol=1e-4)

This is a function to be used with any 2D model where each element has six degrees of freedom, e.g., any 2D frame (steel or concrete), with or without infills (as long as they're modeled as 6DOF elements) should work.

Inputs are:
- recordName: name of the txt file with the values of acceleration. e.g., 'GM01.txt'. Only one value per line
- dtrec: dt of the record, i.e, time interval in seconds
- nPts: number of points of the file
- dtan: dt to be used by OpenSeesPy in the analysis
- fact: scale factor to be applied to the record
- damp: damping ratio to be applied to the model
- IDctrlNode = control node for the pushover
- IDctrlDOF = DOF for pushover
- elements = list with elements to record forces. e.g, [3,5,6] records forces for elements 3, 5 and 6
- nodes_control = nodes to record displacement in the DOF direction. This is meant to calculate interstory drifts, thereforce you should input one node per building floor
- modes = list to include the number of modes for computing the damping in a Rayleight scheme. Default is [0,2], which uses the first and third mode. Adjust according to your needs.
- Kswitch = 1 to use initial stiffness, 2 for current stiffness. Default is 1
- Tol = tolerance for equilibrium. The algorithm uses the NormUnbalance test, which is appropriate for most circumstances and a 1e-4 value. Adjust the tolerance according to your problem size.

Outputs are:
- time: time of the analysis
- roof: roof displacement
- element_forces: numpy 3D array with the element forces. Think of this as a stack of 2D matrices, where each 2D matrix contains the response of one element: in the rows the analysis time time and in the columns the element DOF.
- node_disp: 2D matrix where each row is time and columns have the node displacement in the specified DOF
- node_vel: same as node_disp but for velocity
- node_acel: same as node_disp but for acceleration. Returns the floor relative acceleration.
- drift: same as node_disp but for interstory drifts.
