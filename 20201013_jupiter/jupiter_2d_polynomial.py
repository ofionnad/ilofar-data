#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd 
import astropy.units as u
from sunpy.time import TimeRange
import time
import sys

import datetime as dt
import matplotlib.dates as mdates
from mpl_toolkits.axes_grid1 import make_axes_locatable

from sk import LofarRaw 

start = time.time()
filename = sys.argv[1]

#observation parameters
frange = [15, 60]
#fname = 'udpoutput/jupiter-stokesI_0_2020-10-13T17:47:00_19563125244140'
sbs = np.arange(76, 320)
beam1_sbs = np.arange(76, 198)
obs_mode = 3
trange = TimeRange("2020-10-13T17:45:00", 15.*u.min)
"""
fname = 'udpoutput/jupiter-stokesI_0_2020-10-13T17:47:00_19563125244140'
sbs = np.arange(76, 319)
obs_mode = 3
trange = TimeRange("2020-10-13T17:47:00", 10.*u.minute)
"""
#just plotting parameters
xlabel = "Time from 2020/10/13 17:45:00.0000"
ylabel = "Frequency (MHz)"
title = sys.argv[1].split('/')[1]

# read data
raw = LofarRaw(fname=filename, sbs=sbs, obs_mode=obs_mode, frange=frange)
# time resolution was scrubbed x16 in udpPM

flimits = np.array([40,94])
freqlimits = raw.sb_to_f(flimits+sbs[0], obs_mode)

xlims = list(map(dt.datetime.fromisoformat, [trange.start.value, trange.end.value]))
xlims = mdates.date2num(xlims)
print(freqlimits[0].value)

raw.data = raw.data[:, flimits[0]:flimits[1]]

#defining some of the rfi dominated bands and masking those
rfi_bands = [0, 1, 25, -1, -2, -3] #dominant rfi bands 0->20 + 31, 32, 38
rfi_freqs = raw.sb_to_f(np.array(rfi_bands) + flimits[0] + sbs[0], obs_mode)
print("\nMasked subbands: {}\nMasked frequencies: {}\n".format(np.array(rfi_bands)+flimits[0], rfi_freqs))

data_mask = np.zeros_like(raw.data) #make an empty mask
data_mask[:,rfi_bands] = 1 #mask the columns with rfi in them
raw.masked_data = np.ma.masked_array(raw.data, data_mask, copy=True) #apply mask to the data
print("RFI dominated channels masked")

raw.data = raw.masked_data

#normalise data
#raw.data_normalise(dark_obs=False, use_xarray=False, read_norms=False, title=title)


arrs = np.array_split(raw.data, 6000)

for arr in arrs:
    norms = np.quantile(arr, 0.1, axis=0)
    arr = arr/norms

print("Normalisation done...\n")

plt.figure(figsize=(12,6))
im = plt.imshow(raw.data.T, aspect='auto', origin='lower',
            vmin=np.nanpercentile(raw.data.T, 5), 
            vmax=np.nanpercentile(raw.data.T, 95),
            extent=[xlims[0], xlims[1], 
                    22.65625, 34.375])

plt.gca().xaxis_date()
date_format = mdates.DateFormatter('%H:%M:%S')
plt.gca().xaxis.set_major_formatter(date_format)

divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)

plt.tight_layout()
plt.savefig('datamask_test.png')
plt.close()

print("take single-handed pol")
rh = True #only take right handed polarisation
lh = False
if rh:
    raw.data[raw.data < 0] = 0
    print("Just right-handed pols")
elif lh:
    raw.data[raw.data > 0] = 0
    print("Just left-handed pols")
print("Done")

print(raw.data.max(), raw.data.min())

"""
for arr in arrs:
    
    #q95 = np.quantile(arr, 0.95, axis=0)
    #for i,limit in enumerate(q95):
    #    arr[np.argwhere( arr[:,i] > limit), i] = limit
    #arr[np.argwhere(arr > q95)] = q95
    
    for i in np.arange(arr.shape[1]):
        if arr[:, i].any() == None:
            pass
        else:
            coeffs = np.polyfit(np.linspace(0, len(arr[:, i]), len(arr[:, i])), arr[:, i], 2)
            ynew = np.poly1d(coeffs)
            arr[:, i] = ynew(np.linspace(0, len(arr[:, i]), len(arr[:, i])))

print("polyfits done")
"""

data_mask = np.zeros_like(raw.data) #make an empty mask
data_mask[:,rfi_bands] = 1 #mask the columns with rfi in them
raw.masked_data = np.ma.masked_array(raw.data, data_mask, copy=True) #apply mask to the data

raw.data = raw.masked_data

print("masked data")



arrs = np.array_split(raw.data, 6000)

for arr in arrs:
    norms = np.quantile(arr, 0.1, axis=0)
    arr = arr/norms

#raw.masked_data = raw.data
#raw.do_stft(None, ci=0.95, wind_div=300., skips=rfi_bands) #returns raw.clean_data


#normalise data again along frequency direction
#raw.data_normalise(dark_obs=False, use_xarray=False, read_norms=False, title=title)

plt.figure(figsize=(12,6))
im = plt.imshow(raw.data.T, aspect='auto', origin='lower',
            vmin=np.nanpercentile(raw.data.T, 5), 
            vmax=np.nanpercentile(raw.data.T, 95),
            extent=[xlims[0], xlims[1], 
                    freqlimits[0].value, freqlimits[1].value])

plt.gca().xaxis_date()
date_format = mdates.DateFormatter('%H:%M:%S')
plt.gca().xaxis.set_major_formatter(date_format)

plt.title("Jupiter auroral emission - Stokes V - RCP")
plt.ylabel(ylabel)
plt.xlabel('Time on 2020-10-13')

divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)

plt.tight_layout()
plt.savefig('jupiter-stokesV-righthanded.png')