import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 
import astropy.units as u
from sunpy.time import TimeRange
import time
import sys
import os

import datetime as dt
import matplotlib.dates as mdates
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm
from matplotlib import gridspec,transforms

sys.path.insert(1, '../../20201013_jupiter/')
from sk import LofarRaw
import ed_flagger as edf
import ilofar_ueds_rfi as ued


filename = '../udpoutput/uranus-stokesVectors_0_2020-12-15T20:04:00_19629670898060' #0 is I and 1 is V I assume
sV = '../udpoutput/uranus-stokesVectors_1_2020-12-15T20:04:00_19629670898060'
off_fname = '../udpoutput/offsource-stokesVectors_0_2020-12-15T20:04:00_19629670898060'
frange = [15,30]
sbs = np.arange(76,198)
obs_mode = 3
time_len_mins = 176.
trange = TimeRange(filename.split('_')[-2], time_len_mins*u.min)
xlabel = "Time on {} (UTC)".format(filename.split('_')[-2].split('T')[0])
ylabel = "Frequency (MHz)"
title = filename.split('/')[2]
no_sbs = 78
nsplit = 10


#on-beam
rawdata = LofarRaw(fname=filename, sbs=sbs, obs_mode=obs_mode, frange=frange)
rawdata.data = rawdata.data[:,:no_sbs] #need to do this because of the way subbands were set up for uranus observations! (only use 78 subbands!)
#off-beam
rawoffdata = LofarRaw(fname=off_fname, sbs=sbs, obs_mode=obs_mode, frange=frange)
rawoffdata.data = rawoffdata.data[:,:no_sbs]
#stokes V
sV_data = LofarRaw(fname=sV, sbs=sbs, obs_mode=obs_mode, frange=frange)
sV_data.data = sV_data.data[:,:no_sbs]


df_chunk = ued.data_chunker(rawdata.data, nsplit)
off_chunk = ued.data_chunker(rawoffdata.data, nsplit)
sV_chunks = ued.data_chunker(sV_data.data, nsplit)
tchunks = trange.split(nsplit)


strings = ["StokesI", "StokesI_OFF", "StokesV"]
#strings = ["StokesV"]
for i,df in enumerate([df_chunk, off_chunk, sV_chunks]):
    #for i,df in enumerate([sV_chunks]):
    total_f_sum = []
    print("Analysing {}".format(strings[i]))
    for n,df_split in enumerate(df):
        print("Analysing Chunk #{}".format(n+1))
        ylims, xlims = ued.get_data_lims(sbs, obs_mode, no_sbs, tchunks[n])
        print("Removing RFI")
        df_norfi = ued.resample_dataset(df_split, f=1220) #resampled to 100 ms resolution to mask rfi shorter than this
        rfi_mask = edf.basic_filter(df_norfi, 4.)
        rfi_begone = np.ma.MaskedArray(df_norfi, rfi_mask)
        df_rs = ued.resample_dataset(rfi_begone, f=4) #this is the 2nd resample so be wary of numbers here. 
        df_filtered = edf.basic_filter(df_rs.T, 4.) #4 sigma
        ued.plot_data(df_filtered, None, None, xlims, ylims, xlabel, ylabel, "ed_flags/{}".format(strings[i])+str(n+1), "Electrostatic Filtering {}".format(strings[i]), gs=True)
        f_sum = edf.get_hist_dynspec(df_filtered, axis=0)
        total_f_sum.append(f_sum)
        #t_sum = edf.get_hist_dynspec(df_filtered, axis=0)
    plt.close()
    f = plt.figure()
    f.set_facecolor('white')
    plt.bar(np.arange(f_sum.shape[0]), np.sum(total_f_sum, axis=0))
    plt.savefig("ed_flags/histogram_frequency_{}".format(strings[i]))
    os.rename('trigger_indices.dat', 'trigger_indices_{}.dat'.format(strings[i]))