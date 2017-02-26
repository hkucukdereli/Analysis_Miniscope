import pandas as pd
import numpy as np
from imagingIO import loopMice, loadData, getBeh
from events import find_nearest

def markTrials(mice, dataList, base, duration, eventType, behType, trials, baselining=False):
    """
    """
    # Load the events
    fileList = loopMice(mice, behType)
    eventList = getBeh(mice, fileList['Behaviour'], behType)

    eventsData = pd.DataFrame()
    for mus, sess in mice:

        # Find the events
        eventTimes = np.array([])
        df = eventList[mus][eventType].dropna().reset_index(drop=True)
        for i, event in enumerate(df):
            ind, nearest = find_nearest(dataList[mus].index.values, event)
            eventTimes = np.append(eventTimes, nearest)

        fs = dataList[mus].index[1] - dataList[mus].index[0]
        newTime = np.linspace(base, 10000*fs,10000)
        for col in dataList[mus].columns:
            if not col == 'Time (s)':
                data = dataList[mus][col]
                for i, event in enumerate(eventTimes[trials[0]-1:trials[1]]):
                    slicedData = pd.DataFrame()

                    if baselining:
                        basedFF = data.loc[event+base+0.0001:event+0.0001].mean()
                    elif not baselining:
                        basedFF = 0
                    slicedData['Fluoro'] = data.loc[event+base+0.0001:event+duration+0.0001].values - basedFF
                    slicedData['Cell'] = col
                    slicedData['Event'] = i+1
                    slicedData['New_Time'] = np.round(np.linspace(base,((len(slicedData['Fluoro'])-1)*fs)+base,len(slicedData['Fluoro'])), 2)
                    ##print np.linspace(base,((len(slicedData['Fluoro'])-1)*fs)+base,len(slicedData['Fluoro']))
                    ##print np.round(np.linspace(base,((len(slicedData['Fluoro'])-1)*fs)+base,len(slicedData['Fluoro'])), 2)
                    #np.arange(+base,duration,fs)

                    eventsData = eventsData.append(slicedData)

    print "\n", len(mice), " mice were loaded."
    if behType == 'FR1':
        for mus, sess in mice:
            print "Mouse number", mus, " had ", eventList[mus]['Right_Count'].max(), " total rewards."

    return eventsData

if __name__ == "__main__":
    import pandas as pd
    import numpy as np
    #import scipy.signal as sig
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from imagingIO import loopMice
    from trials import markTrials
    from analysis import dffCalc, filtData, smoothData, normData

    # Start running the analysis
    mice = [(8404, 6), (8857, 3), (8864, 1)]
    behType = 'FR1'
    fs = 0.05
    base = -5.0
    duration = 30
    trials=[1,12]

    # Locate the files
    fileList = loopMice(mice, behType)
    # Load the events
    eventList = getBeh(mice, fileList['Behaviour'], behType)
    # Load the data
    dataList = loadData(mice, behType)
    # dFF and plot again
    filtList = filtData(mice, dffCalc(mice, dataList, lowest=False), cutoff=5.0, order=6)
    dFFList = smoothData(mice, filtList, window=4)

    eventType = 'Eat_Start'
    eventsData = markTrials(mice, dataList, base=base, duration=duration, eventType=eventType, behType=behType, trials=trials, baselining=True)

    sem=eventsData.pivot_table(index=['Event', 'Cell'], columns='New_Time').sem()
    time = np.arange(base,duration,fs)

    print len(sem), len(time)
