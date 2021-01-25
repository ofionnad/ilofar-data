"""
This file is just for functions that will filter out possible electrostatic discharges from dynspecs.

There are a couple methods for doing this, and it is a work in progress. 
"""

import numpy as np 


def get_aeff(wavelength, d):
    """
    The the aeff of lofar lba at different frequencies
    """
    aeff1 = wavelength**2 /3.
    aeff2 = np.pi*d / 4. 
    return min([aeff1, aeff2])

def get_std_dev(bandwidth, aeff, int_time, wavelength, d):
    """
    Get the standard deviation of a radio antenna/interferometer
    """
    aeff = get_aeff(wavelength, d)


def get_bg(data):
    """
    Get the background of data (the mean value)

    data = 1d array representing a single frequency channel
    """
    return np.nanmean(data), np.nanstd(data)

def basic_filter(dynspec, trigger_threshold):
    """
    dynspec - 2d array of (freqs, time)

    Temporal resolution should be on the order of the maximum electrostatic discharge duration
    default = 300 ms
    """

    #initialize an empty boolean array
    mask = np.zeros_like(dynspec)

    for ii, i in enumerate(dynspec):
        bg_mean, std_deviation = get_bg(i)
        dyn_threshold = trigger_threshold*std_deviation
        #selects for each frequency band
        pre = bg_mean #initialize the first pre comparison to the mean
        post = i[1] #initialize this as the second datapoint in the array 
        for jj, j in enumerate(i):
            diff_1 = np.abs(j - pre)
            diff_2 = np.abs(j - post)
            if diff_1 > dyn_threshold and diff_2 > dyn_threshold and j > (bg_mean + dyn_threshold):
                mask[ii,jj] = 1
            else:
                pass
    
    return mask 
            