import csv
import numpy as np
import imageio
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import functools
from functools import partial
import scipy.stats as st

def read_CV(file, cycle,ref):
    """
    Objective
    ----------
    Read a CV DTA file  and output indicies of cycles, time of each data measurement, Voltage, Current,
    and number of cycles

    Parameters
    ----------
    file : :obj:`numpy.ndarray`
        CV file
    cycle : :obj:`numpy.ndarray`
        desired cycle(s) of data to be extracted
    ref : :obj:`numpy.float64`
        reference electrode voltage

    """
    #all voltage data is wrt Ag/AgCl ref electrode
    cycle=cycle-1
    file_as_text = []
    with open(file, newline='') as csvfile:
         spamreader = csv.reader(csvfile, delimiter='\t', quotechar='|')
         for row in spamreader:
            file_as_text.extend([row])


    header = []
    for row in file_as_text:
        if 'CURVE' not in row[0]:
            header.extend([row]) #tells us the indices of rows before CURVE
        else:
            break

    curve_data_raw = []
    curve_data_IV = []
    indicies = []
    time=[]
    V=[]
    I=[]
    numcyc=0

    for row in file_as_text[len(header)::]: #:: means start:stop:step minus the stop
    #start at line with CURVE and stop whenever we are done
    # string[::2] reads “default start index, default stop index, step size is two—take every second element”.
        i = -1
        if 'CURVE' in row[0]:
            j = 0 #row
            i+=1 #i=0 for CURVE row, i means column index #this is just a counting variable for rows in the empty array
            curve_data_raw.extend([[]]) #make empty arrays
            curve_data_IV.extend([[]])
            numcyc+=1
        else:
            j += 1 #j=1 for headers line
            curve_data_raw[i].extend([row]) #in each column slot in a row, add index  of stuff in matrix
            if j > 2:#this is the data
                curve_data_IV[i].extend([[float(val) for val in row[1:5]]])
            #    [cycle][row][index]
        numcyc=numcyc
    curve_data = [np.array(data).transpose() for data in curve_data_IV]
    #print(len(curve_data[0::][0]))

    #want to only look at 220
    if cycle < 0:
        for q in range(len(curve_data[0::])-1):
            indicies.append(curve_data[q][0])
            time.append(curve_data[q][1])
            V.append(curve_data[q][2])
            I.append(curve_data[q][3])
    else:
        indicies= curve_data[cycle][0]
        time= curve_data[cycle][1]
        V= curve_data[cycle][2]
        I= curve_data[cycle][3]

    return indicies, time, V, I,numcyc

def read_CCD(file):
    file_as_text = []
    with open(file, newline='') as csvfile:
         spamreader = csv.reader(csvfile, delimiter='\t', quotechar='|')
         for row in spamreader:
            file_as_text.extend([row])
    file_as_text.pop(-1)
    header = []
    for row in file_as_text:
        if 'CURVE' not in row[0]:
            header.extend([row]) #tells us the indices of rows before CURVE
        else:
            break

    curve_data_raw = []
    curve_data_IV = []
    indicies = []
    for row in file_as_text[len(header)::]: #:: means start:stop:step minus the stop
    #start at line with CURVE and stop whenever we are done
    # string[::2] reads “default start index, default stop index, step size is two—take every second element”.
        i = -1
        if 'CURVE' in row[0]:
            j = 0 #row
            i+=1 #i=0 for CURVE row, i means column index
            curve_data_raw.extend([[]]) #make empty arrays
            curve_data_IV.extend([[]])
        else:
            j += 1 #j=1 for headers line
            curve_data_raw[i].extend([row]) #in each column slot in a row, add index  of stuff in matrix
            if j > 2:#this is the data
                #if 'QUANT' in row[j]:
                    #break
                curve_data_IV[i].extend([[float(val) for val in row[1:5]]])
            #    [cycle][row][index]


    curve_data = [np.array(data).transpose() for data in curve_data_IV]
    #for row in curve_data_IV:
        #indicies= curve_data_IV[cycle][row][0]

    indicies= curve_data[cycle][0]
    time= curve_data[cycle][1]
    V= curve_data[cycle][2]
    I= curve_data[cycle][3]
    return indicies, time, V, I

def FindSOC(V,I,SRP,T,M):
    """
    Objective
    ----------
    To find the SOC of the electrolyte given its current and voltage from CV

    Parameters
    ----------
    V : :obj:`numpy.ndarray`
        An array of voltage values [V]
    I : :obj:`numpy.ndarray`
        An array of current values [A]
    SRP : :obj:`numpy.float64`
        Standard reduction Potential experimentally measured of theoretically found in literature [V]
    T : :obj: `numpy.float64`
        temperature [K]
    M : :obj: `numpy.float64`
        Total molarity of active material in battery system [M]

    """
    T=298#K
    n=1#electron
    F=96485 #Faraday's constant [coulomb/mol]
    R=8.314 #J/mol-K
    FRP_range = np.where(np.abs(I)<1e-9) #V is in V, I is in A
 
    
    FRP=np.mean(V[FRP_range[0]]) #FRP is in V

    #use the nernst equation to find the FindSOC
    z=(-FRP+SRP)*(F*n/(R*T)) #unitless
    #z=(FRP-SRP)*(F*n/(R*T)) #unitless
    p=np.exp(z) #unitless
    ox=p*M/(1+p) #solve for concentration of oxidized species#M
    SOC= ox*100/M #solve for FindSOC#%

    return SOC, FRP, ox

def convert_images_to_gif(filenames, filename_final, fps=10):
    #take the png images generated and make it into a gif
    images = []
    for filename in filenames:
        images.append(imageio.imread(filename))
    imageio.mimsave(filename_final, images, fps=fps)

def fit_ko(eta,ko,oxidation):
    """
    Objective
    ----------
    To fit current and overpotential data to the BV equation and return ko

    Parameters
    ----------
    eta : :obj:`numpy.ndarray`
        An array of overpotential values [V]
    oxidation : :obj:`numpy.float64`
        A concentration value of oxidized species [M]
    Formal : :obj: `numpy.float64`
        The formal reduction potential from this set of CV data [V]
    ko : :obj: to be fitted
        concentration of oxidized species resulting from CV Data [M]

    """
    
    A = 7.85e-7 # area in cm^2 of CF microelectrode
    Tsys = 298 #K
    Rgas = 8.314 #J/mol-K
    F = 96485 #C/mol
    alpha = 0.5
    C0 = oxidation #mol/L concentration of oxidized species
    CR = 1 - C0 #mol/L concentration of reduced species

    return -ko/1000*F*A*(C0*np.exp(-alpha*F/Rgas/Tsys*eta)-CR*np.exp((1-alpha)*F/Rgas/Tsys*eta))

  
  
