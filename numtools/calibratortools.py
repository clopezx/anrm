# Fits ANRM 1.0 (Irvin et. al 2013) against single-cell measurements
# of caspase reporters.
import numpy as np
import scipy as sp 
from scipy.interpolate import *
""" TODO
    insure that there is protein dynamics data and initial concentration data
    for each condition in experiment. For k in ydata.keys if k not in init_conc report Error:
    Exerimental Conditions must match. And the dictionaries have to be the same length.
    
    Should I require the same observable for each condition? Mathematically, it is possible to
    do the calibration without each condition having the same variable. I will allow them to be dipserate.
    """
def normalize(data, option = None):
    ymin = data.min(0)
    ymax = data.max(0)
    if option == None or option == 0:
        return (data - ymin) / (ymax - ymin) # max and min set to 1 and 0, respectively
    elif option == 1:
        return (data/ymax) # max  set to 1.
    else:
        print "ERROR: Normalize options must be 0 or 1"

def normalize_var(values, variance, option = None):
    ymin = values.min(0)
    ymax = values.max(0)
    if option == None or option == 0:
        return variance / (ymax - ymin)**2 # max and min set to 1 and 0, respectively
    elif option == 1:
        return variance / ymax**2 # max  set to 1.
    else:
        print "ERROR: Normalize options must be 0 or 1"

def normalize_array(data_array, option = None, normalize_variance = True):
    timepts = data_array[:,0]
    values  = normalize(data_array[:,1], option)
    if normalize_variance:
        variance = normalize_var(data_array[:,1], data_array[:, 2], option)
        return np.array(zip(*np.vstack((timepts,values,variance))))
    else:
        return np.array(zip(*np.vstack((timepts,values))))

def initial_conditions(names, values, model_icparams):
    """Return a list of initial condition values, eaching having the same index as in model.parameters_initial_conditions."""
    prot_names  = [p.name for p in model_icparams]
    prot_values = [p.value for p in model_icparams]

    for i,n in enumerate(prot_names):
        if n in names:
            prot_values[i] = values[names.index(n)]
    return prot_values

def extract_records(recarray, names):
    """Convert a record-type array and list of names into a float array"""
    return np.vstack([recarray[name] for name in names]).T

def cubic_spline(xsim, ysim, xdata, degree = None):
    """Approximates the simulated data at time = xdata using cubic spline"""
    if degree is not None:
        tck = sp.interpolate.splrep(xsim, ysim, k = degree)
        yspline = sp.interpolate.splev(xdata, tck)
    else:
        tck = sp.interpolate.splrep(xsim, ysim)
        yspline = sp.interpolate.splev(xdata, tck)
    return yspline

def calculate_time_delay(signal, tspan):
    norm_signal = normalize(signal, option = 0)
    
    if np.isnan(np.sum(norm_signal)):
        return None
    else:
        norm_signal = norm_signal.tolist()
        idx         = norm_signal.index(min(norm_signal, key = lambda x: abs(x-0.5)))
        return cubic_spline(norm_signal[idx-3:idx+3], tspan[idx-3:idx+3], [0.5], degree = 3)
