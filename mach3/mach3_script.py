import numpy as np
from array import array
import matplotlib.pyplot as plt
import math
from scipy.optimize import curve_fit

# global variables

filename = 'mb62lu.log'

# Functions
def read_scan(start_line,finish_line):
    """This function should read the lines in the
    given range and create an array from them."""
    
    f = open(filename)

    length = finish_line - start_line
    
    scan_array = np.empty([length+1,48], dtype = int)
    output_array = np.zeros(48,dtype=int)

    lines = f.readlines()
    
    for n in range(length+1):
        for i in range (48):
            if lines[start_line+n-1][i:i+1] in '0123456789' and lines[start_line+n-1][i:i+1] not in '':
                scan_array[n,i] = int(lines[start_line+n-1][i:i+1])
            else: scan_array[n,i] = 0


    for i in range(48):
        for n in range(length+1):
            output_array[i] += scan_array[n,i]*(10**(length-n))
    output_array = np.transpose(output_array)
    f.close()

    # Find out the omega range and the stepsize
    
    f = open(filename)
    line = f.readlines()
    omega_range = float(line[finish_line][51:56])
    f.close()
    step_size = omega_range/47.0

    # Create array of degrees 

    degrees_array = np.zeros(48,dtype = float)
    for i in range(48):
        degrees_array[i] = round(-omega_range/2 + step_size*i,4)

    # Combine angles and intensities
    #final = np.concatenate((degrees_array, output_array),axis = 1)

    return degrees_array, output_array

def scan_search(filename):
    """This function search for the scan in the log file and returns
    arrays of start and end numbers of scans"""
    f = open(filename)
    begin_array = array('i',[])
    end_array = array('i',[])
    num = 0
    for line in f:
        num += 1
        if 'Total' in line[:5]:
            begin = num + 1
            begin_array.append(begin)
        if 'ADAS' in line[:4]:
            end = num - 1
            end_array.append(end)
    f.close()
    return begin_array, end_array

def gauss(x,*p):
    """Gauss function"""
    A, mu, sigma, bkg = p
    return A*np.exp((-(x-mu)**2)/(2.*(sigma**2)))+bkg

def gauss_fit(x_data,y_data):
    """This function makes Gauss fit"""
    global coeff
    p0 = [100.0,0.0,0.15,0.0]
    try:
        coeff, var_matrix = curve_fit(gauss,x_data,y_data,p0=p0,maxfev = 10000)
    except RuntimeError:
        coeff = [0.,0.,0.,0.]
    fit = gauss(x_data,*coeff)
    return fit

def represents_flt(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def represents_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

def hkl_extract(start_line):
    f = open(filename)
    lines = f.readlines()
    integers = array('i',[0,0,0,0,0,0])
    counter1 = 0
    counter2 = 0
    while counter1 <= 5:
        broken_line = str.split(lines[start_line-counter2-3])
        broken_line = reversed(broken_line)
        for word in broken_line:
            if represents_int(word) and word != '+1':
                integers[counter1] = int(word)
                counter1 +=1
            elif represents_flt(word) and word != '+1':
                counter1 += 1
        counter2 +=1

    hkl = str(integers[5])+ ' '  + str(integers[4]) + ' ' + str(integers[3])
    f.close()
    return hkl

# Program body

grid_size = math.ceil(math.sqrt(len(scan_search(filename)[0])))

if grid_size*(grid_size-1) >= len(scan_search(filename)[0]):
    plt.figure(figsize=(4*grid_size, 5*(grid_size-1)))
else:
    plt.figure(figsize=(4*grid_size, 3*grid_size))

for refl_number in range(len(scan_search(filename)[0])):

       # determine the grid size based on the number of reflections found in the .log file
    if grid_size*(grid_size-1) >= len(scan_search(filename)[0]):
        ax = plt.subplot(grid_size,grid_size-1,refl_number+1)               # create subplot
    else:
        ax = plt.subplot(grid_size,grid_size,refl_number+1)

    x_data = read_scan(scan_search(filename)[0][refl_number],scan_search(filename)[1][refl_number])[0]
    y_data = read_scan(scan_search(filename)[0][refl_number],scan_search(filename)[1][refl_number])[1]

    plt.plot(x_data,gauss_fit(x_data,y_data))
    plt.plot(x_data,y_data)

    FWHM = math.fabs(round(2.35482*coeff[2],3))
    ax.text(0.95,0.9,'FWHM='+str(FWHM), ha ='right', va = 'center', size = 12,transform=ax.transAxes)
    ax.text(0.20,0.9,hkl_extract(scan_search(filename)[0][refl_number]), ha ='right', va = 'center', size = 12,transform=ax.transAxes)

plt.show()