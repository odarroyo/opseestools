# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 15:34:48 2023

@author: Dirsa
"""
#%%

import numpy as np                                                              # Importa librería estándar para operaciones matematicas.
                                                         # Importa librería estándar de análisis de datos.
import matplotlib.pyplot as plt                                                 # Importa librería estándar para plotear otras cosasy poder crear figuras.
                                       # Importa librería estándar que tiene un conjunto de herramientas para proporcionar una canalización ligera en Python
# Ajuste de distribuciones
# ==============================================================================
from scipy import stats, optimize
import pandas as pd

#%% CURVAS DE FRAGILIDAD AJUSTE FUNCIÓN LOGNORMAL
#-----Example calculations to demonstrate the use of fragility fitting
#-----functions. These calculations are based on the following paper:
#-----Baker, J. W. (2015). “Efficient analytical fragility function fitting 
#-----using dynamic structural analysis.” Earthquake Spectra, 31(1), 579-599.

#-----Created by Jack Baker
#-----2/4/2013
#-----Modified by Jack Baker, 1/25/2017, to update citation information

#-----objective function to be optimized
def mlefit(theta, num_gmrs, num_collapse, IM):
    if theta[0]<0:                                                               # don't let median of fragility function go below zero
        theta[0]=0
    #-----estimated probabilities of collapse, given the current fragility functionparameter estimates
    p = stats.norm.cdf(np.log(IM), np.log(theta[0]), theta[1])
    
    #-----likelihood of observing num_collapse(i) collapses, given num_gms
    #-----observations, using the current fragility function parameter estimates
    likelihood = stats.binom.pmf(np.transpose(num_collapse), np.transpose(num_gmrs), np.transpose(p))
    
    #-----sum negative log likelihood (we take the nevative value because we want
    #-----the maximum log likelihood, and the function is searching for a minimum)
    loglik = (-1)*sum(np.log(likelihood));
    return loglik
#---- example data: IM levels, number of analyses, and number of collapses
def fn_mle_pc(IM, num_gmrs, num_collapse):
    #-----by Jack Baker
    #-----10/9/2012
    #-----Modified by Gemma Cremen, 1/25/2017, to avoid estimating negative median
    #-----values for the fragility function
    #-----Modified by Jack Baker, 1/25/2017, to update citation information
    
    #-----This function fits a lognormal CDF to observed probability of collapse 
    #-----data using optimization on the likelihood function for the data. 
    #-----These calculations are based on equation 11 of the following paper:
    
    #-----Baker, J. W. (2015). “Efficient analytical fragility function fitting 
    #-----using dynamic structural analysis.” Earthquake Spectra, 31(1), 579-599.
    
    #-----INPUTS:
    #-----IM            1xn           IM levels of interest
    #-----num_gms       1x1 or 1xn    number of ground motions used at each IM level
    #-----num_collapse 	1xn           number of collapses observed at each IM level
    
    #-----OUTPUTS:
    #-----theta         1x1           median of fragility function
    #-----beta          1x1           lognormal standard deviation of fragility function
    
    
    
    #-----Initial guess for the fragility function parameters theta and beta. 
    #-----These initial choices should not need revision in most cases, but they 
    #-----could be altered if needed.
    #x0 = np.array([0.8, 0.4], dtype=float)
    x0 = [0.8, 0.4]
    args=(num_gmrs, num_collapse, IM)
    x = optimize.minimize(mlefit,x0,args=args,method='Nelder-Mead',options={'maxiter': 1000}) #Minimization of scalar function of one or more variables con el metodo de Lagarias, J.C., J. A. Reeds, M. H. Wright, and P. E. Wright, "Convergence Properties of the Nelder-Mead Simplex Method in Low Dimensions," SIAM Journal of Optimization, Vol. 9 Number 1, pp. 112-147, 1998.
    theta = x['x'][0]
    beta = x['x'][1]
    return theta, beta

def plotfrag(theta,beta,x = np.linspace(0,3,100)):
    y = stats.lognorm.cdf(x,s=beta,scale=theta)
    plt.plot(x,y)
    
def values_in_bins(data, bins='fd'):
    """
    Divide the values of a NumPy array into histogram bins and return the values in each bin.
    
    Parameters:
    - data: NumPy array of data points.
    - bins: Number of bins or a sequence of bin edges.
    
    Returns:
    - bins_values: A list of NumPy arrays, each containing the values in each bin.
    - bin_edges: The edges of the bins.
    - bin_midpoint: midpoint of the bin
    """
    # Determine the bins
    bin_counts, bin_edges = np.histogram(data, bins=bins)
    
    # Initialize a list to hold the values for each bin
    bins_values = [[] for _ in range(len(bin_edges)-1)]
    bins_indices = [[] for _ in range(len(bin_edges)-1)]
    bin_midpoints = (bin_edges[:-1] + bin_edges[1:]) / 2
    # Group the values into the corresponding bins
    for index, value in enumerate(data):
        for i in range(len(bin_edges)-1):
            if bin_edges[i] <= value < bin_edges[i+1]:
                bins_values[i].append(value)
                bins_indices[i].append(index)
                break
    
    # Convert lists to NumPy arrays for numerical operations
    bins_values = [np.array(bin) for bin in bins_values]
    bins_indices = [np.array(bin) for bin in bins_indices]
    
    return bins_values, bins_indices, bin_midpoints, bin_counts, bin_edges

def calculate_fragility(df,limit_name,limits,IM_column,EDP_column,plot=True):
    ''' Function to calculate fragility based on a dataframe.
    The dataframe must have a column named bin which containst the values for
    the intensity measure. It also must have a column that you 
    
    Inputs:
        df: dataframe
        limit_name: how you want to name your limits
        limits: values to define the limits
        column_limit: dataframe column to check the exceedance of the values in limits
        plot: by default is True, meaning that it returns the plots. Set to false if you do not want them
    
    Outputs:
        thetas: median of the lognormal distribution
        betas: deviation of the lognormal distribution
    
    '''
    ngm = df[IM_column].value_counts().to_numpy()
    
    for i in range(len(limits)):
        df[limit_name[i]] = df[EDP_column] > limits[i]

    true_counts = df.groupby(IM_column)[limit_name].apply(lambda x: x.sum()).reset_index()
    thetas, betas = [],[]
    
    for i,lim in enumerate(limit_name):
        theta,beta = fn_mle_pc(true_counts[IM_column], ngm, true_counts[lim])
        thetas.append(theta)
        betas.append(beta)
        if plot==True:
            plt.plot(true_counts[IM_column],true_counts[lim]/ngm[i],'x')
            plotfrag(theta, beta)
        
    if plot==True:
        plt.xlabel('Sa(g)')
        plt.ylabel('Probability')
        plt.show()
    
    return thetas,betas