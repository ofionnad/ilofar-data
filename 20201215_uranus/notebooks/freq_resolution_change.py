import analyse_ilofar_ueds as lofardynspec
import sys

sys.path.insert(1, '../../20201013_jupiter/')
from sk import LofarRaw

def fft_change_res(data):
    """
    Function to change the data frequency resolution
    """
    for subchannel in data:
        fft_channel = np.fft.fft()
    return new_data

if __name__=="__main__":

    #observation specs
    filename = '../udpoutput/uranus-stokesVectors_0_2020-12-15T20:04:00_19629670898060' #0 is I and 1 is V I assume
    off_fname = '../udpoutput/offsource-stokesVectors_0_2020-12-15T20:04:00_19629670898060'
    cal_fname = '../udpoutput/cygA-stokesVectors_0_2020-12-15T20:00:00_19629667968374' #calibrator has different obs length!

    plot_names = 'test_script/Uranus_StokesI_'
    plot_title = "Uranus observation - Stokes I"

    frange = [15,30]
    sbs = np.arange(76,198)
    obs_mode = 3
    time_len_mins = 176.
    trange = TimeRange(filename.split('_')[-2], time_len_mins*u.min)
    xlabel = "Time from {} {}".format(filename.split('_')[-2].split('T')[0], filename.split('_')[-2].split('T')[1])
    xlabel = "Time on {} (UTC)".format(filename.split('_')[-2].split('T')[0])
    ylabel = "Frequency (MHz)"
    title = filename.split('/')[2]
    no_sbs = 78 #number of usable subbands!

    #how much to split up data into
    nsplit = 10
    r_factor = 12207

    #on-beam
    rawdata = LofarRaw(fname=filename, sbs=sbs, obs_mode=obs_mode, frange=frange)
    rawdata.data = rawdata.data[:,:no_sbs] #need to do this because of the way subbands were set up for uranus observations! (only use 78 subbands!)

    ylims, xlims = get_data_lims(no_sbs, trange)

    df_chunk = data_chunker(rawdata.data, nsplit)

    for i,j in enumerate(df_chunk):
        data = fft_change_res()