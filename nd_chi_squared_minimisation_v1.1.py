# -*- coding: utf-8 -*-
"""
Created on Sun Apr 14 11:20:12 2024

@author: Jakub Dombrowski

Chi squared minimisation algorithm using scipy.optimize.fmin extended to any 
number of dimensions with uncertainty calculations made from meshgrids of 
the variables.

To use, use get_best_chi_fit(arguments) and get_best_chi_fit_error(arguemnts)
functions 
"""

import numpy as np
import scipy.optimize as op

###############################################################################
######################### CHI MINIMISATION VALUES #############################

def get_chi_squared(fun_inputs, x_observed, y_observed, dy_observed, function,\
                    known_constants = []):
    """
    x_observed: np array float
    y_observed: np array float
    dy_observed: np array float
    function: function type with inputs (x_observed, *fun_inputs, 
                                         *known_constants)
    ARG
    known_constants = python 1d list with known constants you want to input 
    into your function
    
    return chi squared: float
    standard chi squared formula SUM((yobs-yexp)**2/dyobs**2)
    returns the chi squared value (float)
    """
    return np.sum((y_observed-function(x_observed,\
                            *fun_inputs, *known_constants))**2/dy_observed**2)


def get_best_chi_fit(x_observed, y_observed, dy_observed, starting_values,\
                     function, known_constants = []):
    """
    USE THIS
    Performs the minimisation fit via scipy fmin
    x_observed: NUMPY ARRAY 1D
    y_observed: NUMPY ARRAY 1D
    dy_observed: NUMPY ARRAY 1D
    starting_values: 1D LIST/ARRAY THAT INDEXES WITH fitted_values[INDEX]
    function: the function that you want to model, has to have parameters in
    form function(X, *starting_values, *known_constants)

    ARGS #############################
    known_constants: python 1d list of constants that are known and necessary
                     for function
    returns python list with best return[0] = fit_values, return[1] = chi
    squared
    """
    try:
        vals = op.fmin(get_chi_squared, starting_values, \
        args = (x_observed, y_observed, dy_observed, function, \
                   known_constants), full_output = True, disp = True)
    except ValueError:
        vals = [np.array([None]*len(starting_values)), None]
        print("Check the dimensions and values of your input arrays")
    return vals


###############################################################################
########################### CHI MINIMISATION ERROR ############################
def get_best_chi_fit_error(x_observed, y_observed, dy_observed, \
                     function, fitted_vals, chi_val, \
                         percentage_of_variable=0.01, resolution = 100,\
                             known_constants = []):
    """
    USE THIS
    performs meshgrid operations to find the uncertainties on fitted values
    with chi_val + 1, i.e. 1 standard deviation
    x_observed: NUMPY ARRAY 1D
    y_observed: NUMPY ARRAY 1D
    dy_observed: NUMPY ARRAY 1D
    fitted_values: 1D LIST/ARRAY THAT INDEXES WITH fitted_values[INDEX]
    function: the function that you want to model, has to have parameters in
    form function(X, *starting_values)
    ARGS ########################
    resolution:     int, size/rank of mesh grid
    percentage_of_variable:     float < 1, selects range over which the mesh is
    taken.
    known_constants:    python 1d list of constants that are used in function

    returns python 1d length len(fitted_values) with uncertainties in same 
    order as fitted_vals
    """
    try:
        variables = []
        current_unc = [0] * len(fitted_vals)
        #generates initial meshes of values
        for variable in fitted_vals:
            variables.append(np.linspace(variable*(1-percentage_of_variable), \
            variable*(1+percentage_of_variable), resolution))
        meshes = np.meshgrid(*variables)
        chi_mesh = get_chi_squared_mesh(x_observed, y_observed, dy_observed, \
                                    meshes.copy(), function, known_constants)
        target = chi_val + 1
        dt = target*0.25
        indicies = d_contour(chi_mesh, target, dt)
        for index, fitted_value in enumerate(fitted_vals):
            for j in indicies:

                working_unc = abs(read_var_multidimension_index(meshes[index], j) \
                              - fitted_value)
                if working_unc > current_unc[index]:
                    current_unc[index] = working_unc
        if current_unc == [0]*len(fitted_vals):
            print("failed to find uncertainties check resolution and range")
    except ValueError:
        print("Check the dimensions and values of your input arrays")
        current_unc = np.array([None]*len(fitted_vals))
    return current_unc


def get_chi_squared_mesh(x_observed, y_observed, dy_observed, meshes,\
                         function, known_constants = []):

    """
    converts numpy meshes to python list, creates a chi squared mesh (outputs
    as dictionary for convenience) with indexed values corresponding to each
    mesh

    x_observed: NUMPY ARRAY 1D
    y_observed: NUMPY ARRAY 1D
    dy_observed: NUMPY ARRAY 1D
    meshes: python list of nd numpy arrays
    function: the function that you want to model, has to have parameters in
    form function(X, *starting_values)

    """
    chi_mesh = {}
    indicies = []
    dimension = len(meshes)
    coordinates = [0]*dimension
    rank = len(meshes[0])

    for i in range(0,dimension):
        meshes[i] = meshes[i].tolist()

    for i in range(0, dimension):
        work_indicies = [0 for i in range(0, dimension)]

        indicies.append(work_indicies)

    for j in range(0, int(rank)):
        for i in range(0, dimension):
            working_value = meshes.copy()[:]

            indicies[i].insert(i,j)
            indicies[i].pop(i+1)
            coordinates[i] = j
            indicies.reverse()
            for k in range(0, dimension):
                for l in range(0, dimension):

                    item = working_value[l].copy().pop(indicies[l][k])
                    working_value[l] = item

            indicies.reverse()
            final_value = working_value
            chi_mesh[str(coordinates)] = get_chi_squared(final_value,\
                                x_observed, y_observed, dy_observed, function,\
                                    known_constants)
    return chi_mesh



def d_contour(tensor_dict, target, dt):
    """
    finds the contour of a mesh (in dictionary format with key = "[x,y,...,w]" 
                                 for a target value with an uncertainty
    tensor_dict: dictionary of form dict["[x,y,..,w]"] = value
    target: float
    dt: float
    returns python list 2d (1d list of indicies which are themselves a 1d 
                            list each] of the contour
    """
    indicies = []
    for i in tensor_dict:
        index = i.split(",")
        index[0] = index[0].split("[")[1]
        index[-1] = index[-1].split("]")[0]
        for for_index, item in enumerate(index):
            index[for_index] = int(item)
        if target-dt <= tensor_dict[i] <= target+dt:
            indicies.append(index)
    return indicies


def read_var_multidimension_index(tensor, index):
    """
    tensor: any nd array
    index: 1d list in form [x,y,...,w] that contains the indicies that you want
    to access in tensor, len(index) < rank
           of tensor
    outputs: whatever the value inside of the tensor is indexed.
    """
    if len(index) == 1:
        return tensor[index[0]]
    else:
        values = read_var_multidimension_index(tensor[index[0]], index[1:])
        return values
