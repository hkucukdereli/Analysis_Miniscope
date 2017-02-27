import pandas as pd
import numpy as np
from imagingIO import loopMice, loadData, getBeh
from events import find_nearest, getEvents
from statistics import getBoutMeans, getStats

def getBouts(mice, dataList, eventStart, eventEnd, baseEvent, base, behType, trials, baselining=False):
    """
    """
    eventsData = pd.DataFrame()
    for mus, sess in mice:

        for i in range(trials[0], trials[1]+1):
            ind1, nearest1 = find_nearest(dataList[mus].index.values, eventStart[mus].loc[i])
            ind2, nearest2 = find_nearest(dataList[mus].index.values, eventEnd[mus].loc[i])
            ind3, nearest3 = find_nearest(dataList[mus].index.values, baseEvent[mus].loc[i])
            eventStart[mus].loc[i] = nearest1
            eventEnd[mus].loc[i] = nearest2
            baseEvent[mus].loc[i] = nearest3

        fs = dataList[mus].index[1] - dataList[mus].index[0]
        newTime = np.arange(base, 10000, fs)

        for col in dataList[mus].columns:
            if not col == 'Time (s)':
                data = dataList[mus][col]
                for j in range(trials[0], trials[1]+1):
                    startT = eventStart[mus].loc[j]
                    endT = eventEnd[mus].loc[j]
                    baseT = baseEvent[mus].loc[j]

                    slicedData = pd.DataFrame()

                    if baselining:
                        if base < 0:
                            basedFF = data.loc[baseT+base:baseT+0.0001].mean()
                        elif base > 0:
                            basedFF = data.loc[baseT:baseT+base+0.0001].mean()
                    elif not baselining:
                        basedFF = 0.0

                    if base < 0:
                        slicePre = data.loc[baseT+base:baseT].values - basedFF
                        slicePost = data.loc[startT:endT].values - basedFF
                        dataSlice = np.concatenate((slicePre, slicePost), axis=0)
                        slicedData['Fluoro'] = dataSlice
                        slicedData['New_Time'] = newTime[0:len(dataSlice)]
                        #slicedData['New_Time'] = np.linspace(base,((len(dataSlice)-1)*fs)+base,len(dataSlice))
                    elif base > 0:
                        dataSlice = data.loc[startT+base:endT+0.0001].values - basedFF
                        slicedData['Fluoro'] = dataSlice
                        slicedData['New_Time'] = np.linspace(base,((len(dataSlice)-1)*fs)+base,len(dataSlice))

                    slicedData['Cell'] = col
                    slicedData['Event'] = j
                    eventsData = eventsData.append(slicedData)

    print "\n", len(mice), " mice were loaded."

    if baselining:
        print "\nBaseline was set", base, "sec before the event."
    elif not baselining:
        print "\nTraces were not baselined."

    return eventsData

def getBoutDur(mice, eventType, behType, trials):
    """
        DEPRICATED: Needs major revision!
    """
    print "DEPRICATED: Needs major revision!"
    # Load the data
    dataList = loadData(mice)

    # Load the events
    fileList = loopMice(mice, behType)
    eventList = getBeh(mice, fileList['Behaviour'], behType)

    durationData = pd.DataFrame()
    for mus, sess in mice:
        # Find the events
        boutDur = np.array([])
        startList = eventList[mus][eventType[0]].dropna().values
        endList = eventList[mus][eventType[1]].dropna().values
        for i in range(trials[0]-1, trials[1]):
            ind, start = find_nearest(dataList[mus].index.values, startList[i])
            ind, end = find_nearest(dataList[mus].index.values, endList[i])
            boutDur = np.append(boutDur, end-start)

        durationData[mus] = boutDur

    durationData.set_index(np.arange(trials[0], trials[1]+1), inplace=True)

    return durationData

def markBouts(mice, dataList, eventType, behType, trials, base, dff=True, baseline=False):
    """
        DEPRICATED: Use markBouts instead.
    """
    print "DEPRICATED: Use getBouts instead."
    fileList = loopMice(mice, behType)
    eventList = getBeh(mice, fileList['Behaviour'], behType)

    eventsData = pd.DataFrame()
    for mus, sess in mice:

        # Find the events
        eventTimes1 = np.array([])
        eventTimes2 = np.array([])
        startTimes = eventList[mus][eventType[0]].dropna().reset_index(drop=True)
        endTimes = eventList[mus][eventType[1]].dropna().reset_index(drop=True)
        for i, event in enumerate(startTimes):
            ind1, nearest1 = find_nearest(dataList[mus].index.values, startTimes.loc[i])
            ind2, nearest2 = find_nearest(dataList[mus].index.values, endTimes.loc[i])
            eventTimes1 = np.append(eventTimes1, nearest1)
            eventTimes2 = np.append(eventTimes2, nearest2)

        fs = dataList[mus].index[1] - dataList[mus].index[0]
        for col in dataList[mus].columns:
            if not col == 'Time (s)':
                data = dataList[mus][col]
                for i, event in enumerate(eventTimes1[trials[0]-1:trials[1]]):
                    start = eventTimes1[i]
                    end = eventTimes2[i]

                    slicedData = pd.DataFrame()

                    if baseline:
                        basedFF = data.loc[start+base:start+0.0001].mean()
                    elif not baseline:
                        basedFF = 0.0

                    slicedData['Fluoro'] = data.loc[start+base:end+0.0001].values - basedFF
                    slicedData['Cell'] = col
                    slicedData['Event'] = i+1
                    if len(data.loc[start+base:end+0.0001].values - basedFF) == len(np.arange(base,end-start+0.0001,fs)):
                        slicedData['New_Time'] = np.arange(base,end-start+0.0001,fs)
                    if len(data.loc[start+base:end+0.0001].values - basedFF) < len(np.arange(base,end-start+0.0001,fs)):
                        slicedData['New_Time'] = np.arange(base,end-start,fs)

                    eventsData = eventsData.append(slicedData)

    print "\n", len(mice), " mice were loaded."
    if behType == 'FR1':
        for mus, sess in mice:
            print "Mouse number", mus, " had ", eventList[mus]['Right_Count'].max(), " total rewards."

    if baseline:
        print "\nBaseline was set", baseline, "sec before the event."
    elif not baseline:
        print "\nTraces were not baselined."

    return eventsData

if __name__ == "__main__":
    import pandas as pd
    import numpy as np
    #import scipy.signal as sig
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from imagingIO import loopMice
    from trials import markTrials

    # Start running the analysis
    mice = [(8404, 6), (8857, 3), (8864, 1)]
    behType = 'FR1'
    fs = 0.05
    base = -5.0
    trials=[1,12]

    # Locate the files
    fileList = loopMice(mice, behType)
    # Load the events
    eventList = getBeh(mice, fileList['Behaviour'], behType)
    # Load the data
    dataList = loadData(mice, behType)

    eventStart = getEvents(mice, eventList, 'Eat_Start', trials)
    eventEnd = getEvents(mice, eventList, 'Eat_End', trials)
    baseEvent = getEvents(mice, eventList, 'Eat_Start', trials)

    #boutsData = markBoutsDouble(mice, dataList, eventStart, eventEnd, baseEvent, base=base, behType=behType, trials=trials, dff=True, baselining=True)
    #boutMeans = getBoutMeans(mice, boutsData, behType=behType, trials=trials, base=base, dff=True, baseline=False)

    # Run a bunch of stats
    #[stats, perc] = getStats(boutMeans, choice='Ranksum', p_val=0.05)

    #print boutsData
    #print boutsData.pivot_table(index=['Cell', 'Event'], columns='New_Time', values='Fluoro')
    #print perc

    #markBoutsDouble(mice, dataList, eventStart, eventEnd, baseEvent, base=base, behType=behType, trials=trials, dff=True, baselining=True)
    eventType = ['Eat_Start', 'Eat_End']
    #markBouts(mice, dataList, eventType=eventType, behType=behType, trials=trials, base=base, dff=True, baseline=False)

    boutsData =  getBouts(mice, dataList, eventStart, eventEnd, baseEvent, base, behType, trials, baselining=False)
