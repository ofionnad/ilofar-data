import sigpyproc as spp
import argparse
from astropy.time import Time
import numpy as np
from tqdm import trange
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as ticker
import matplotlib.dates as mdates

from mpl_toolkits import axes_grid1

import datetime
import scipy.fft
from scipy.signal import stft, hann


def rollingAverage(data, step = 8):
	rollingSum = np.cumsum(data)
	return rollingSum[step:] - rollingSum[:-step]

def decimate(data, step = 64):
	rollingSum = np.cumsum(data)
	return rollingSum[step::step] - rollingSum[:-step:step]


filename = '/mnt/ucc4_data1/data/Uranus/2021_11_30/Uranus_IE613_2021-11-30T23:00:00.000_stokesI.fil'

filReader = spp.FilReader(filename)

time = Time(filReader.header.tstart, format = "mjd")
time_filReader = Time(filReader.header.tstart, format = "mjd")
samplesPerBlock = filReader.header.nsamples
readTimestamp = 0
dataBlock_all = filReader.readBlock(readTimestamp, samplesPerBlock)

# P_AA = np.zeros([filReader.header.nchans,int(samplesPerBlock/4)])
# P_BB = np.zeros([filReader.header.nchans,int(samplesPerBlock/4)])
# P_AB = np.zeros([filReader.header.nchans,int(samplesPerBlock/4)])
# P_BA = np.zeros([filReader.header.nchans,int(samplesPerBlock/4)])

# for j in range(int(samplesPerBlock/4)):
#     P_AA[:,j] = dataBlock_all[:,j*4+0]
#     P_BB[:,j] = dataBlock_all[:,j*4+1]
#     P_AB[:,j] = dataBlock_all[:,j*4+2]
#     P_BA[:,j] = dataBlock_all[:,j*4+3]

stokesI = spp.Filterbank.FilterbankBlock(dataBlock_all, dataBlock_all.header)
# dataBlock_I=spp.Filterbank.FilterbankBlock(P_AA**2+P_BB**2,dataBlock_all.header)
# dataBlock_V=spp.Filterbank.FilterbankBlock(-2*P_BA,dataBlock_all.header)

#stokesI = spp.Filterbank.Filterbank(stokesI)
c = filReader.getChan(140)

# N = 500
# w = hann(N)
# yf, xt, ft = stft(stokesI[:,58600000:117200000], fs=10000, window='hann', nperseg=2048)

#xf = scipy.fft.fftfreq(N)
#xf = scipy.fft.fftshift(xf)

# print(ft.shape)
# print(yf.shape)
# print(xt.shape)

#plt.pcolormesh(xt, yf, np.abs(ft), shading='gouraud')
#plt.imshow(ft.real, aspect='auto', extent=[xf[0], xf[-1], None, None], origin='lower')

# plt.imshow(stokesI[58600000:117200000, :], origin='lower', aspect='auto')

plt.plot(c)
plt.tight_layout()
plt.savefig('FTs/sk_test.png')