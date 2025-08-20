# opseestools

This is a library to be used with OpenSeesPy.

It includes a collection of analysis functions that are useful if you are analyzing frames or wall buildings both in two and three dimensions, storing the results in Python variables. 
It also includes another useful tools for processing results and creating scripts

## How to cite

The opseestools paper is available at SoftwareX, so you can cite it if you use it in your research: https://www.sciencedirect.com/science/article/pii/S2352711024002036

Arroyo, O., Feliciano, D., Novoa, D., & Valcárcel, J. (2024). opseestools: A Python library to streamline OpenSeesPy workflows. SoftwareX, 27, 101832.

## Installation

opseestools can be installed from the Python Package Index as:

#### pip install opseestools

## Description

opseestools comprises a set of functions in four modules:

1) 2D Analysis functions. This library is called analysis. You can import it as:

import opstools.analisis as an

2) 3D Analysis functions. This library is called analysis. You can import it as:

import opstools.analisis3D as an3D (or any name you want)

3) A set of utilities to support model building. You can import it as:

import opstools.utilidades as ut (or any other name)

4) A library to calculate and plot fragility functions. You can import it as:

import opstools.Lib_frag as ut (or any other name)

opseestools documentation is available at: https://opseestools.readthedocs.io/en/latest/index.html

## Examples of use

#### You can find Jupyter Notebooks and python scripts with examples of using the library on the examples folder in this GitHub repository. The list of examples are the following:

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

Example 12: using functions in v0.6 to generate a RC frame building and apply a cyclic pushover analysis using a function available in the analysis module.

Example 13: example of a 3D building subjected to a dynamic analysis using a function available in the analysis module.

Example 14: example of a new function implemented in v0.65 that calculates the confinement for an RC section according to Mander model.

Example 15: Showing how to use a new function where we integrate opseestools with opstool, another great library. It also demonstrated how to create a video for the dynamic analysis using that library.

Example 16: same as example 15, but for a pushover analysis.

Example 17: same as example 15 but for a 2D model created using opseestools functions.

Example 18: same as example 17 but for a dynamic analysis.

## Function help

Each function has a docstring with help for the function, detailing the inputs and outputs. 

You can summon this help by using the integrated Python help(), for example: help(an.gravedad) will return help for the gravedad function in the analysis module.
