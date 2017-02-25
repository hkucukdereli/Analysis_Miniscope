import numpy as np
import pandas as pd
from scipy.signal import butter, lfilter, freqz
from imagingIO import get_fs

def dffCalc(mice, dataList, lowest=False):
    """
    """
    dFFList = {}
    for mus, sess in mice:
        data = dataList[mus]

        if lowest:
            # F0 from the lowest <lowest>% of the trace
            t_half = data.index.max() * (lowest / 100) * 0.5
            minPoint = data.idxmin(axis= 0)
            dFF = pd.DataFrame()
            for i, m in enumerate(minPoint):
                F0 = data[data.columns[i]].loc[m-t_half : m+t_half].mean(axis= 0)
                dFF[data.columns[i]] = (data[data.columns[i]] - F0) / F0 * 100
        else:
            # Calculate the dFF
            F0 = data.mean(axis= 0)
            dFF = ((data - F0) / F0) * 100

        # Put each dataframe into the new dict
        dFFList[mus] = dFF

    return dFFList

def low_pass(data, cutoff, fs, order):
    """
    """
    # set the butterworth filter
    nyq = 0.5 * (1./fs)
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)

    # filter the data
    filtData = lfilter(b, a, data)

    return filtData

def filtData(mice, dataList, cutoff, order=4):
    """
    """
    filtList = {}
    for mus, sess in mice:
        data = dataList[mus]

        fs =  get_fs(data.index)

        for cell in data.columns:
            filtData = low_pass(data[cell], cutoff, fs, order)
            data[cell] = filtData
        filtList[mus] = data

    return filtList

def smoothData(mice, dataList, window=1):
    """
    """
    smoothList = {}
    for mus, sess in mice:
        data = dataList[mus]

        smoothData = data.rolling(window=window).mean()
        smoothList[mus] = smoothData

    return smoothList

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from imagingIO import loopMice, loadData, loadBeh, getBeh

    # Start gathering the data from FR1 exp
    # Mice
    mice = [(8404, 6), (8857, 3), (8864, 1)]
    # Parameters
    behType = 'FR1'
    fs = 0.05
    base = 5.0
    duration = 30.0
    trials=[1,12]

    # Get file paths
    fileList = loopMice(mice, behType)

    # Load the data
    dataList = loadData(mice, behType)
    trace_raw = dataList[8857][dataList[8857].columns[5]]


    # Plot the dFF data
    plt.figure(figsize=(20,10), facecolor="w", dpi= 150)

    # Only dFF
    dFFList = dffCalc(mice, dataList, lowest=.60)
    trace_dff = dFFList[8857][dFFList[8857].columns[5]]
    plt.plot(trace_dff.index, trace_dff.values, 'gray')

    # Low pass filter, then dFF
    filtList = filtData(mice, dFFList, cutoff=5.0, order=6)
    smoothList = smoothData(mice, filtList, window=3)
    trace_filt = smoothList[8857][smoothList[8857].columns[5]]
    plt.plot(trace_filt.index, trace_filt, 'r')
    #plt.xlim([400, 500])
    #plt.ylim([-1, 3])
    plt.show()
