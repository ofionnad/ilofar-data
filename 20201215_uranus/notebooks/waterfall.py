import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import ilofar_ueds_rfi as ued
import sys

import astropy.units as u
import matplotlib.dates as mdates
import matplotlib.patheffects as PathEffects
import matplotlib as mpl
sys.path.insert(1, '../../20201013_jupiter/')
from sk import LofarRaw

def waterfall_plot(data, xlims, sbs, obs_mode, sbs_interval=10, interval=10000):
    num_plots = int(data.shape[1]/sbs_interval)
    #print(num_plots)
    f, axs = plt.subplots(num_plots, 1, figsize=(10,10), sharex=True)
    x = np.linspace(xlims[0], xlims[1], data.shape[0])
    c = np.arange(1, int(data.shape[1]/sbs_interval)+1) #+1 so the last few plots are not completely white
    norm = mpl.colors.Normalize(vmin=c.min(), vmax=c.max())
    cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.cool)
    cmap.set_array([])

    leg = LofarRaw.sb_to_f(sbs[0]+np.linspace(0, data.shape[1], data.shape[1]+1)[::sbs_interval], obs_mode)
    leg = [round(x.value, 2) for x in leg]*u.MHz

    data = data[:,::sbs_interval]

    for i,ax in enumerate(axs):
        ax.plot(x, data[:,num_plots-i], c='k')
        ax.fill_between(x, data[:,num_plots-i], facecolor=cmap.to_rgba(i+1))
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.axhline(y=data.min(), lw=2, clip_on=False, color=cmap.to_rgba(i+1))
        ax.set_yticks([])
        ax.set_ylim([data.min(),data.max()])
        ax.patch.set_alpha(0)
        label(ax, cmap.to_rgba(i+1), leg[num_plots-i])
        date_format = mdates.DateFormatter('%H:%M:%S.%f')
    axs[-1].spines['bottom'].set_visible(True)
    ax.xaxis.set_major_formatter(date_format)
    ax.xaxis.set_minor_locator(mdates.MicrosecondLocator(interval=interval))
    ax.set_xlabel('Time', fontsize=14)
    ax.tick_params(labelsize=14)
    f.subplots_adjust(hspace=-.5)
    f.autofmt_xdate()


def waterfall_defunct(data):
    """
    This works for a 2d dynamic spectrum where data is [time, frequency]
    """

    flatten = np.ndarray.flatten(data)
    f = np.tile(list(np.arange(data.shape[1])), data.shape[0])
    df = pd.DataFrame(dict(x=flatten, y=f))

    pal = sns.cubehelix_palette(80, rot=-.25, light=.7)
    g = sns.FacetGrid(df, row="y", hue="y", aspect=15, height=.5, palette=pal)

    # Draw the densities in a few steps
    g.map(sns.kdeplot, "x",
        bw_adjust=.5, clip_on=False,
        fill=True, alpha=1, linewidth=1.5)
    g.map(sns.kdeplot, "x", clip_on=False, color="w", lw=2, bw_adjust=.5)
    g.map(plt.axhline, y=0, lw=2, clip_on=False)
    g.map(label, "x")
    g.fig.subplots_adjust(hspace=-0.8)
    g.set_titles("")
    #axes = g.axes
    #for i in g.axes:
    #    for ax in i:
    #        ax.set_xlim([0,70000])
    g.set(yticks=[])
    g.despine(bottom=True, left=True)

def label(ax, color, label):
    txt = ax.text(-0.1, 0.2, label, fontweight="normal", color=color, fontsize=14,
        ha="left", va="center", transform=ax.transAxes)
    txt.set_path_effects([PathEffects.withStroke(linewidth=2, foreground='k')])

def zoom_window(raw_obj, data, timechunk, n_splits, split_return, no_sbs):
    """
    I use this to zoom in on individual events to check if they look like EDs

    raw_obj - just the LofarRaw class from sk.py, let's me use all the observation variables
    data - dynspec [time, freq]
    timechunk - sunpy timerange object covering datarange (probably should have changed this in all my codes to be astropy timerange oh well)
    n_splits - how many splits of the data to make along time
    split_return - the number split to return from the function (from 0 to n_splits-1)
    """
    tsplit = timechunk.split(n_splits)[split_return]
    ylims, xlims = ued.get_data_lims(raw_obj.sbs, raw_obj.obs_mode, no_sbs, tsplit)
    newdata = np.array_split(data, n_splits)[split_return]
    return newdata, tsplit, ylims, xlims


def plot_event_1d(data, xlims, sbs, obs_mode, sbs_interval=10, interval=10000):
    """
    This just shows a couple 1d plots of an event up close acrdeductionoss sbs_interval frequencies (total_sbs/10)
    """
    f, ax = plt.subplots(figsize=(10,6))
    f.set_facecolor('white')
    plt.plot(np.linspace(xlims[0], xlims[1], data.shape[0]), data[:,::sbs_interval])
    plt.xlim(xlims[0], xlims[1])
    plt.gca().xaxis_date()
    date_format = mdates.DateFormatter('%H:%M:%S.%f')
    ax.xaxis.set_major_formatter(date_format)
    ax.xaxis.set_minor_locator(mdates.MicrosecondLocator(interval=interval))
    f.autofmt_xdate()
    plt.xlabel('Time')
    plt.ylabel('Arbitrary voltages')
    leg = LofarRaw.sb_to_f(sbs[0]+np.linspace(0, data.shape[1], data.shape[1]+1)[::sbs_interval], obs_mode)
    leg = [round(x.value, 2) for x in leg]*u.MHz
    plt.legend(leg)
    plt.grid(ls=':')
    plt.tight_layout()

def plot_event_waterfall(data, xlims, sbs, obs_mode, sbs_interval=10, interval=10000):

    f, axs = plt.subplots(int(data.shape[1]/sbs_interval), 1, figsize=(10,10), sharex=True)
    x = np.linspace(xlims[0], xlims[1], data.shape[0])
    c = np.arange(1, int(data.shape[1]/sbs_interval)+1) #+1 so the last few plots are not completely white
    norm = mpl.colors.Normalize(vmin=c.min(), vmax=c.max())
    cmap = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.YlGnBu)
    cmap.set_array([])

    leg = LofarRaw.sb_to_f(sbs[0]+np.linspace(0, data.shape[1], data.shape[1]+1)[::sbs_interval], obs_mode)
    leg = [round(x.value, 2) for x in leg]*u.MHz

    for i,ax in enumerate(axs):
        ax.plot(x, data[:,i], c='k')#c=cmap.to_rgba(i + 1), linecolour='w')
        ax.fill_between(x, data[:,i], facecolor=cmap.to_rgba(i+1))
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        #ax.axhline(y=0, lw=0.5, clip_on=False, color='k')
        ax.set_yticks([])
        ax.patch.set_alpha(0)
        label(ax, cmap.to_rgba(i+1), leg[i])
        date_format = mdates.DateFormatter('%H:%M:%S.%f')
    axs[-1].spines['bottom'].set_visible(True)
    ax.xaxis.set_major_formatter(date_format)
    ax.xaxis.set_minor_locator(mdates.MicrosecondLocator(interval=interval))
    f.subplots_adjust(hspace=-.7)
    f.autofmt_xdate()