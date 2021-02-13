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

def basic_filter(dynspec, trigger_threshold, stokesV=False):
    """
    dynspec - 2d array of (freqs, time)

    Temporal resolution should be on the order of the maximum electrostatic discharge duration
    default = 300 ms
    """

    #initialize an empty boolean array
    mask = np.zeros_like(dynspec)

    for ii, i in enumerate(dynspec):
        if stokesV:
            """
            Need to split up the data if it is stokes V to get proper rcp and lcp comparison.
            Absolute pain
            """
            ipos = i.copy()
            ipos[ipos < 0] = 0
            ineg = i.copy()
            ineg[ineg > 0] = 0
            
            mean_pos, std_pos = get_bg(ipos)
            pos_thresh = trigger_threshold*std_pos
            mean_neg, std_neg = get_bg(ineg)
            neg_thresh = trigger_threshold*std_neg
            pre_pos, pre_neg = mean_pos, mean_neg
            post_pos, post_neg = ipos[1], ineg[1]
            for jj,j in enumerate(ipos):
                diff_1 = np.abs(j-pre_pos)
                diff_2 = np.abs(j-post_pos)
                if diff_1 > pos_thresh and diff_2 > pos_thresh and j > (mean_pos + pos_thresh):
                    mask[ii,jj] = 1
                else:
                    pass
                pre_pos = j
                try:
                    post_pos = i[jj+2]
                except:
                    post = mean_pos    
            for jj,j in enumerate(ineg):
                absj = np.abs(j)
                diff_1 = absj-np.abs(pre_neg)
                diff_2 = absj-np.abs(post_neg)
                if diff_1 > neg_thresh and diff_2 > neg_thresh and (absj > (np.abs(np.mean(j)+np.std(j))) or absj < (np.abs(np.mean(j) - np.std(j)))):
                    mask[ii,jj] = 1
                else:
                    pass
                pre_neg = j
                try:
                    post_neg = i[jj+2]
                except:
                    post_neg = mean_neg
        else:
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
                    with open('trigger_indices.dat', 'a+') as f:
                        f.write("{}\n".format(ii))
                else:
                    pass
                #increment the comparison pixels
                pre = j
                try:
                    post = i[jj+2]
                except:
                    post = bg_mean
            #print("Subchannel {} done".format(77+ii))
    
    return mask.T

def get_hist_dynspec(data, axis=1):
    """
    Quick function to the sum along one axis of an array

    Not really a histogram but for the 0,1 mask for flagging it works
    """

    axis_sum = np.sum(data, axis=axis)
    return axis_sum

