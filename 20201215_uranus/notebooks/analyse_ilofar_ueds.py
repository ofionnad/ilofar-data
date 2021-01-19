#!/usr/bin/env python

import numpy as np
import sys

import datetime as dt

#plotting modules
import matplotlib.pyplot as plt 
import matplotlib.dates as mdates
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm
from matplotlib import gridspec,transforms

#custom module
sys.path.insert(1, '../../20201013_jupiter/')
from sk import LofarRaw

def get_data_lims(sbs_number, trange):
    """
    Return the limits of the data in usable form
    sbs_number = number of subbands used
    trange = the length of time. 

    """
    flimits = np.array([0,sbs_number])
    freqlimits = raw.sb_to_f(flimits+sbs[0], obs_mode)
    xlims = list(map(dt.datetime.fromisoformat, [trange.start.value, trange.end.value]))
    xlims = mdates.date2num(xlims)
    return freqlimits, xlims

def plot_data(data, trange, flims, xlabel, ylabel, plotname):
    f = plt.figure(figsize=(12,6))
    f.set_facecolor('white')

    spec = gridspec.GridSpec(ncols=2, nrows=2, width_ratios=[6,1], height_ratios=[1,4])
    ax0 = f.add_subplot(spec[0])
    ax2 = f.add_subplot(spec[2])
    ax3 = f.add_subplot(spec[3])

    ax0.plot(xs, lw=0.5)
    ax0.set(xticklabels=[], xlim=[0,None])
    ax0.tick_params(bottom=False)

    rot = transforms.Affine2D().rotate_deg(270)
    ax3.plot(ys[::-1], lw=0.5, transform=rot+plt.gca().transData)
    ax3.set(yticklabels=[])
    ax3.tick_params(left=False)

    xlims = list(map(dt.datetime.fromisoformat, [trange.start.value, trange.end.value]))
    xlims = mdates.date2num(xlims)

    im = ax2.imshow(data.T, aspect='auto', origin='lower',
            vmin=np.nanpercentile(data.T, 5), 
            vmax=np.nanpercentile(data.T, 90),
            extent=[xlims[0], xlims[1], 
                    flims[0].value, flims[1].value])

    ax2.xaxis_date()
    date_format = mdates.DateFormatter('%H:%M:%S')
    ax2.xaxis.set_major_formatter(date_format)

    ax2.set_title("Uranus observation - Stokes I")
    ax2.set_ylabel(ylabel, fontsize=14)
    ax2.set_xlabel(xlabel, fontsize=14)

    divider = make_axes_locatable(ax2)
    cax = divider.append_axes("right", size="3%", pad=0.05)
    plt.colorbar(im, cax=cax)

    plt.tight_layout()
    plt.savefig(plotname, dpi=900)

def remove_rfi_wtih_std(data, stdy=2, stdx=2):
    """
    Function to mask out suspected rfi automatically

    data - data to be masked
    stdy - mask threshold for channels in number of standard deviations
    stdx - mask threshold for time in number of standard deviations

    Returns the masked data and the result of the sum along each axis
    """

    #sum along the time axis
    ys = data.T.sum(axis=1)
    sigma_y = stdy*np.std(ys) #2sigma
    ymean = np.mean(ys)
    print(ymean+sigma_y)

    for i,j in enumerate(ys[::-1]):
            if j > (ymean+sigma_y):
                    print("RFI trigger channel : {}".format(i))
                    ys[i] = np.NaN 
                    data[:, i] = np.NaN

    #ys = data.T.sum(axis=1)
    xs = np.nansum(data.T, axis=0)
    sigma_x = stdx.*np.nanstd(xs) #2sigma
    xmean = np.nanmean(xs)
    print(xmean + sigma_x)
    for i,j in enumerate(xs):
        if j> (xmean+sigma_x):
            print("RFI trigger time sample : {}".format(i))
            xs[i] = np.NaN
            data[i, :] = np.NaN

    return xs, ys, data 

def resample_dataset(data, f=12207):
    """
    Function takes the dataset, find the remainder after dividing by the resample factor (f) and pads the array. 
    Then takes the mean of the chunks of data of those resampled size. 

    f - resampling factor 
    """

    pad = np.ceil(float(data.shape[0])/f)*f - data.shape[0]
    pad_arr = np.zeros((int(pad), data.shape[1]))*np.NaN
    data_padded = np.append(data, pad_arr).reshape(-1, data.shape[1])
    r_data = np.nanmean(data_padded.T.reshape(data.shape[1], -1, dec), axis=2).T

    return r_data


if __name__=="__main__":

    #observation specs
    filename = '../udpoutput/uranus-stokesVectors_0_2020-12-15T20:04:00_19629670898060' #0 is I and 1 is V I assume
    off_fname = '../udpoutput/offsource-stokesVectors_0_2020-12-15T20:04:00_19629670898060'
    cal_fname = '../udpoutput/cygA-stokesVectors_0_2020-12-15T20:00:00_19629667968374' #calibrator has different obs length!

    frange = [15,30]
    sbs = np.arange(76,198)
    obs_mode = 3
    time_len_mins = 176.
    trange = TimeRange(filename.split('_')[-2], time_len_mins*u.min)
    xlabel = "Time from {} {}".format(filename.split('_')[-2].split('T')[0], filename.split('_')[-2].split('T')[1])
    ylabel = "Frequency (MHz)"
    title = filename.split('/')[2]

    #on-beam
    data = LofarRaw(fname=filename, sbs=sbs, obs_mode=obs_mode, frange=frange)
    #off-beam
    offdata = LofarRaw(fname=off_fname, sbs=sbs, obs_mode=obs_mode, frange=frange)

    ylims, xlims = get_data_lims(sbs.shape[0], trange)
    