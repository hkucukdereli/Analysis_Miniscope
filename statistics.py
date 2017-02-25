import pandas as pd
import numpy as np
import scipy.stats as stat

def stdCalc(mice, dataList):
    """
        Calculate the standard deviation for each cell.
    """
    stdList = pd.DataFrame()
    for mus, sess in mice:
        stdList = stdList.append(pd.DataFrame({'Std':dataList[mus].std().values}, index=dataList[mus].columns))
    return stdList

def getBoutMeans(mice, boutsData, behType, trials, base, baselining=True):
    # Get the cell names
    cellNames = boutsData['Cell'].drop_duplicates().values

    # Calculate the bout and baseline means
    boutMeans = pd.DataFrame()
    for i, cell in enumerate(cellNames):
        boutDF = pd.DataFrame()
        boutDF['Baseline'] = boutsData.reset_index(drop=True).pivot_table(index=['Event'], columns=['Cell', 'New_Time'], values='Fluoro')[cell].T.loc[base:0.00].mean(axis=0).values
        boutDF['Bout_Mean'] = boutsData.reset_index(drop=True).pivot_table(index=['Event'], columns=['Cell', 'New_Time'], values='Fluoro')[cell].T.loc[0.00:].mean(axis=0).values
        boutDF['Cell'] = cell
        boutDF['Event'] = boutDF['Baseline'].index.values+1
        #boutDF['Event'] = np.arange(trials[0], trials[1]+1)
        boutMeans = boutMeans.append(boutDF)
    print "Dunnzo!"
    return boutMeans

def getBoutMax(mice, boutsData, behType, trials, base, baselining=True):
    # Get the cell names
    cellNames = boutsData['Cell'].drop_duplicates().values

    # Calculate the bout and baseline means
    boutMeans = pd.DataFrame()
    for i, cell in enumerate(cellNames):
        boutDF = pd.DataFrame()
        boutDF['Baseline'] = boutsData.reset_index(drop=True).pivot_table(index=['Event'], columns=['Cell', 'New_Time'], values='Fluoro')[cell].T.loc[base:0.00].max().values
        boutDF['Bout_Mean'] = boutsData.reset_index(drop=True).pivot_table(index=['Event'], columns=['Cell', 'New_Time'], values='Fluoro')[cell].T.loc[0.00:].max().values
        boutDF['Cell'] = cell
        boutDF['Event'] = boutDF['Baseline'].index.values+1
        #boutDF['Event'] = np.arange(trials[0], trials[1]+1)
        boutMeans = boutMeans.append(boutDF)
    print "Dunnzo!"
    return boutMeans

def getStats(boutMeans, choice, p_val=0.05):
    stats = pd.DataFrame()
    cellNames = boutMeans['Cell'].drop_duplicates().values
    for cell in cellNames:
        x = boutMeans.pivot_table(index='Event', columns='Cell', values='Baseline')[cell]
        y = boutMeans.pivot_table(index='Event', columns='Cell', values='Bout_Mean')[cell]
        T_wilc, p_wilc = stat.wilcoxon(x, y)
        T_rank, p_rank = stat.ranksums(x, y)
        if y.mean() > x.mean():
            PInd = 1 * np.absolute((y.mean() - x.mean()) / (y.mean() + x.mean()))
        elif y.mean() < x.mean():
            PInd = -1 * np.absolute((y.mean() - x.mean()) / (y.mean() + x.mean()))
        elif y.mean() == x.mean():
            PInd = 0

        temp = pd.DataFrame()
        temp['Cell'] = [cell]
        temp['Wilcoxon'] = [p_wilc]
        temp['Ranksum'] = [p_rank]
        temp['PInd'] = [PInd]
        temp['Baseline_Mean'] = [x.mean()]
        temp['Bout_Mean'] = [y.mean()]

        if p_wilc < p_val:
            temp['Wilcoxon Result'] = ['*']
        else:
            temp['Wilcoxon Result'] = ['ns']
        if p_rank < p_val:
            temp['Ranksum Result'] = ['*']
        else:
            temp['Ranksum Result'] = ['ns']

        if choice == 'Wilcoxon':
            p = p_wilc
        elif choice == 'Ranksum':
            p = p_rank

        if PInd > 0.0:
            temp['Preference'] = ['Positive']
            if p < p_val:
                temp['Class'] = ['Up']
            else:
                temp['Class'] = ['None']
        elif PInd < 0:
            temp['Preference'] = ['Negative']
            if p < p_val:
                temp['Class'] = ['Down']
            else:
                temp['Class'] = ['None']
        elif PInd == 0:
            if p < p_val:
                temp['Class'] = ['None']
            else:
                temp['Class'] = ['None']

        stats = stats.append(temp)
    stats = stats.reset_index(drop=True)

    upPercent = len(stats[stats['Class'] == 'Up']['Class'])*1.0 / len(stats['Class']) * 100.0
    downPercent = len(stats[stats['Class'] == 'Down']['Class'])*1.0 / len(stats['Class']) * 100.0
    nanPercent = len(stats[stats['Class'] == 'None']['Class'])*1.0 / len(stats['Class']) * 100.0

    percentages = pd.DataFrame({'Up' : [upPercent], 'Down' : [downPercent], 'None' : [nanPercent]}, index=['Percentage (%)'])

    stats.set_index('Cell', drop=True, inplace=True)

    return stats, percentages