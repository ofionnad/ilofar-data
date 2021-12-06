import sigpyproc as spp
import numpy as np

# import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.ticker as ticker
import matplotlib.dates as mdates
from matplotlib import gridspec,transforms
from mpl_toolkits.axes_grid1 import make_axes_locatable

from astropy.time import Time
from sunpy.time import TimeRange
import astropy.units as u
import datetime as dt


def sb_to_f(sbs, obs_mode):
    """
    Convert subband number to a frequency in MHz
    """
    nyq_dict = {3:1, 5:2, 7:3}
    nyq_zone = nyq_dict[obs_mode]
    clock_dict = {3:200, 4:160, 5:200, 6:160, 7:200} #MHz
    clock = clock_dict[obs_mode]
    nu = (nyq_zone-1. + sbs/512.) * (clock/2.)
    return nu * u.MHz


def get_data_lims(sbs, obs_mode, sbs_number, trange):
    """
    Return the limits of the data in usable form
    sbs_number = number of subbands used
    trange = the length of time. 

    """

    flimits = np.array([0,sbs_number])
    freqlimits = sb_to_f(flimits+sbs[0], obs_mode)
    xlims = list(map(dt.datetime.fromisoformat, [trange.start.value, trange.end.value]))
    xlims = mdates.date2num(xlims)

    return freqlimits, xlims


def remove_rfi_wtih_sd(data, stdy=100, stdx=100):
    """
    Function to mask out suspected rfi automatically

    data - data to be masked
    stdy - mask threshold for channels in number of standard deviations
    stdx - mask threshold for time in number of standard deviations

    Returns the masked data and the result of the sum along each axis
    """
    print('\nRemoving RFI...\n')

    #sum along the time axis
    ys = data.T.sum(axis=1)
    sigma_y = stdy*np.std(ys) #2sigma
    ymean = np.mean(ys)
    #print(ymean+sigma_y)

    for i,j in enumerate(ys[::-1]):
        rfi_removed_list = []
        if j > (ymean+sigma_y):
            rfi_removed_list.append(i)
            ys[i] = np.NaN 
            data[:, i] = np.NaN
    if rfi_removed_list is not None:
        print("RFI trigger channel : {}".format(rfi_removed_list))

    #ys = data.T.sum(axis=1)
    xs = np.nansum(data.T, axis=0)
    sigma_x = stdx*np.nanstd(xs) #2sigma
    xmean = np.nanmean(xs)
    #print(xmean + sigma_x)
    for i,j in enumerate(xs):
        rfi_removed_list = []
        if j > (xmean+sigma_x):
            xs[i] = np.NaN
            data[i, :] = np.NaN
    if rfi_removed_list is not None:
        print("RFI trigger time sample: {}".format(rfi_removed_list))

    return xs, ys, data 

def plot_data(data, xs, ys, xlims, flims, xlabel, ylabel, plotname, plot_title, minutes=True, gs=False):
    print("\nPlotting...\n")
    if gs:
        cmap = plt.cm.Greys
    else:
        cmap = plt.cm.viridis

    data = np.flip(data, axis=1)

    if xs is None and ys is None:
        print("Not using summed axes!")
        f, ax = plt.subplots(figsize=(12,6))
        f.set_facecolor('white')
        im = ax.imshow(data.T, aspect='auto', origin='lower', cmap=cmap,
            vmin=np.nanpercentile(data.T, 5), 
            vmax=np.nanpercentile(data.T, 90),
            extent=[xlims[0], xlims[1], 
                    flims[0].value, flims[1].value])

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.05)
        plt.colorbar(im, cax=cax)

        ax.xaxis_date()
        date_format = mdates.DateFormatter('%H:%M:%S.%f')
        ax.xaxis.set_major_formatter(date_format)
        if minutes:
            ax.xaxis.set_minor_locator(mdates.MinuteLocator())
        f.autofmt_xdate()

        ax.set_title(plot_title)
        ax.set_ylabel(ylabel, fontsize=14)
        ax.set_xlabel(xlabel, fontsize=14)

        plt.tight_layout()
        plt.savefig(plotname, dpi=900)

    else:

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

        im = ax2.imshow(data.T, aspect='auto', origin='lower',
                vmin=np.nanpercentile(data.T, 5), 
                vmax=np.nanpercentile(data.T, 90),
                extent=[xlims[0], xlims[1], 
                        flims[0].value, flims[1].value])

        ax2.xaxis_date()
        date_format = mdates.DateFormatter('%H:%M:%S.%f')
        ax2.xaxis.set_major_formatter(date_format)
        if minutes:
            ax2.xaxis.set_minor_locator(mdates.MinuteLocator())
        f.autofmt_xdate()
        
        ax2.set_title("Uranus observation - Stokes I")
        ax2.set_ylabel(ylabel, fontsize=14)
        ax2.set_xlabel(xlabel, fontsize=14)

        divider = make_axes_locatable(ax2)
        cax = divider.append_axes("right", size="3%", pad=0.05)
        plt.colorbar(im, cax=cax)

        plt.tight_layout()
        plt.savefig(plotname, dpi=900)

def data_chunker(df, n_split):
    """
    Function to split up data
    df - dataset to be split
    n_split - number of splits
    """
    chunked = np.array_split(df, n_split, axis=1)

    return chunked


if __name__=='__main__':

    filename = '/mnt/ucc4_data1/data/Uranus/2021_11_30/Uranus_IE613_2021-11-30T23:00:00.000_stokesI.fil'
    tstart = '/mnt/ucc4_data1/data/Uranus/2021_11_30/tstart'
    with open(tstart, 'r') as f:
        tstart = f.read()
    tstart = tstart.split()[2]

    xlabel = "Time on {} (UTC)".format(filename.split('_')[-1].strip('.fil').split('T')[0])
    ylabel = "Frequency (MHz)"
    title = filename.split('/')[-1]
    plot_title = "Uranus observation - Stokes I - IE613 - 2021-11-30"


    filReader = spp.FilReader(filename)

    time = Time(filReader.header.tstart, format = "mjd")
    time_filReader = Time(filReader.header.tstart, format = "mjd")
    samplesPerBlock = filReader.header.nsamples
    readTimestamp = 0
    dataBlock_all = filReader.readBlock(readTimestamp, samplesPerBlock)

    stokesI = spp.Filterbank.FilterbankBlock(dataBlock_all, dataBlock_all.header)
    stokesI = stokesI.normalise()

    time_len = stokesI.shape[1]
    time_res = 5.12e-6*512
    trange = TimeRange(tstart, time_len*time_res*u.second)

    sbs = np.arange(51,461)
    obs_mode=3
    no_sbs=len(sbs)
    ylims, xlims = get_data_lims(sbs, obs_mode, no_sbs, trange)

    nsplit = 10
    df_chunk = data_chunker(stokesI.data, nsplit)
    dts = trange.split(nsplit)
    plot_names = 'all_data_plots/Uranus_StokesI_normalised_'

    for i,j in enumerate(df_chunk):
        print(plot_names+str(i))
        ys = j.sum(axis=1)[::-1]
        xs = j.sum(axis=0)
        _, xlims = get_data_lims(sbs, obs_mode, no_sbs, dts[i])
        plot_data(j.T, xs, ys, xlims, ylims, xlabel, ylabel, plot_names+str(i), plot_title)
        plt.close()
        print('...next')
