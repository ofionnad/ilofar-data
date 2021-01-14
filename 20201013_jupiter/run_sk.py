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

### let's try remove that rfi with sk again
### do it first thing, no normalisation until afterwards


### now mask dominant channels and plot
rfi_bands = [0, 1, 2] #dominant rfi bands 0->20 + 31, 32, 38
rfi_freqs = raw.sb_to_f(np.array(rfi_bands) + flimits[0] + sbs[0], obs_mode)
print("\nMasked subbands: {}\nMasked frequencies: {}\n".format(np.array(rfi_bands)+flimits[0], rfi_freqs))

data_mask = np.zeros_like(raw.data) #make an empty mask
data_mask[:,rfi_bands] = 1 #mask the columns with rfi in them
raw.masked_data = np.ma.masked_array(raw.data, data_mask, copy=True) #apply mask to the data
print("RFI dominated channels masked")

#raw.data = raw.masked_data

raw.do_stft(None, ci=0.9999, wind_div=300., skips=rfi_bands) #returns raw.clean_data
#raw.data = raw.clean_data


###plot

plt.figure(figsize=(12,6))
im = plt.imshow(raw.clean_data.T, aspect='auto', origin='lower',
            vmin=np.nanpercentile(raw.clean_data.T, 5), 
            vmax=np.nanpercentile(raw.clean_data.T, 95),
            extent=[xlims[0], xlims[1], 
                    22.65625, 34.375])

plt.gca().xaxis_date()
date_format = mdates.DateFormatter('%H:%M:%S')
plt.gca().xaxis.set_major_formatter(date_format)

divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(im, cax=cax)

plt.tight_layout()
plt.savefig('sk_test.png')
plt.close()
