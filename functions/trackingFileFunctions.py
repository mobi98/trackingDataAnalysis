# path = "vsFpvB_gm_RepBio1_RepT2_tl_CSV-Table2.csv"
# df = pd.read_csv(path)

# for every frame, count the number of cells that are WT and mutant
def count_strains(df):
    
    """
    Parameters
    ----------
    df: a dataframe of tracking results from Ilastik
       
    Returns
    -------
    time: time points to be plotted in x axis
    wt: counts of wt bacterial cells per time point
    mnt: counts of mutant bacterial cells per time point
       
    """
    
    time = []
    wt = []
    mnt = []
    
    for i in range(0,30):
        n_wt = len(df[(df['cell_type'] == 'wt') & (df['frame'] == i)])
        n_mnt = len(df[(df['cell_type'] == 'mnt') & (df['frame'] == i)])
        time.append(i*10)
        wt.append(n_wt)
        mnt.append(n_mnt)
        
        
    return time, wt, mnt
    
#time, wt, mnt = count_strains(df)

def cellCount(df):
    
    '''Given a dataframe of single strain, counts number of cells over time
    
    Parameters
    ----------
    
    df: pandas dataframe of tracking information
    
    Returns
    -------
    
    time: list of time stamps
    
    cellCount: list of cell counts at each time point
    '''

    cellCount = []
    for i in range(30):
        c = len(df[df['frame'] == i])
        cellCount.append(c)

    time = [i*10 for i in range(30)] 
    
    return time, cellCount


def calculateGrowthRate(cell_counts, time):
    """
    Parameters
    ----------
    cell_counts: a list of cell counts at each time point
    time: a list of time in minutes
    
    Returns
    -------
    r: growth rate of bacteria, r, as the constant of the exponent
    
    """
    from numpy import exp
    A = cell_counts[0]
    ft = cell_counts[15]
    t = time[15]
    
    r = (1/t)*np.log(ft/A)
    return round(r,3)


def pvd_conc(df):
    
    """
    Parameters 
    ----------
    
    df: a pandas dataframe of tracking results from Ilastik
    
    Returns
    -------
    
    time: list of time points
    pvd_wt: average pvd intensity of WT bacteria per time point in time list
    pvd_mnt: average pvd intensity of mutant bacteria per time point in time list
    OR
    time: list of time points
    pvd: average pvd intensity for each time point

    """
    
    time = []
    pvd_wt = []
    pvd_mnt = []
    pvd = []
    
    for i in range(0,30):
        
        time.append(i*10)
        
        if 'cell_type' in df.columns:
            wt = df[(df['frame'] == i) & (df['cell_type'] == 'wt')]
            wtPvd = sum(wt['pvd_intensity'])
            wtPvdAvg = round((wtPvd/len(wt['pvd_intensity'])), 2)
            pvd_wt.append(wtPvdAvg)

            mnt = df[(df['frame'] == i) & (df['cell_type'] == 'mnt')]
            mntPvd = sum(mnt['pvd_intensity'])
            mntPvdAvg = round((mntPvd/len(mnt['pvd_intensity'])), 2)
            pvd_mnt.append(mntPvdAvg)

        
        else:
            pvdSum = sum(df[df['frame'] == i]['pvd_intensity'])
            pvdAvg = pvdSum/len(df[df['frame'] == i])
            pvd.append(pvdAvg)
            
            
    if 'cell_type' in df.columns:
        return time,pvd_wt,pvd_mnt
    else:
        return time,pvd
    
    
    
    
def cellDivisionTime(df):
    
    """
    Parameters
    ----------
    
    df: a dataframe containing tracking information 
    
    Returns
    -------
    
    division_times: a list of division times for every cell present 
    totalCells: total number of cells across the entire
    
    """
    if 'cell_type' in df.columns:
        division_timesWt = []
        division_timesMnt = []
        trackIdsMnt = list(df[(df['trackId'] != -1) & (df['cell_type'] == 'mnt')]['trackId'])
        totalCellsMnt = len(trackIdsMnt)
        trackIdsWt = list(df[(df['trackId'] != -1) & (df['cell_type'] == 'wt')]['trackId'])
        totalCellsWt = len(trackIdsWt)
        
        for tid in trackIdsMnt:
            start_frame = min(df[df['trackId'] == tid]['frame'])
            rowIdxs = list(df.index[df['parentTrackId'] == tid])
            if rowIdxs:
                fod = df.iloc[rowIdxs[0],:]['frame']
                tod = (fod - start_frame)*10
                division_timesMnt.append(tod)
                
        for tid in trackIdsWt:
            start_frame = min(df[df['trackId'] == tid]['frame'])
            rowIdxs = list(df.index[df['parentTrackId'] == tid])
            if rowIdxs:
                fod = df.iloc[rowIdxs[0],:]['frame']
                tod = (fod - start_frame)*10
                division_timesWt.append(tod)

                
    else:
        division_times = []
    
        trackIds = list(df[df['trackId'] != -1]['trackId'])
        totalCells = len(trackIds)

        for tid in trackIds:
            start_frame = min(df[df['trackId'] == tid]['frame'])
            rowIdxs = list(df.index[df['parentTrackId'] == tid])
            if rowIdxs:
                fod = df.iloc[rowIdxs[0],:]['frame']
                tod = (fod - start_frame)*10
                division_times.append(tod)


    if 'cell_type' in df.columns:
        return division_timesWt, division_timesMnt
    else:
        return division_times
        
                       
            
def cellDivisionTime_pvd(df):
    
    """
    Parameters
    ----------
    
    df: a dataframe containing tracking information 
    
    Returns
    -------
    
    division_times: a list of division times for every cell present 
    pvd_intensity: pyoverdine intensity per cell that is recorded
    
    """
    if 'cell_type' in df.columns:
        
        division_timesWt = []
        pvd_intensityWt = []
        division_timesMnt = []
        pvd_intensityMnt = []
        
        trackIdsMnt = list(df[(df['trackId'] != -1) & (df['cell_type'] == 'mnt')]['trackId'])
        totalCellsMnt = len(trackIdsMnt)
        trackIdsWt = list(df[(df['trackId'] != -1) & (df['cell_type'] == 'wt')]['trackId'])
        totalCellsWt = len(trackIdsWt)
        
        for tid in trackIdsMnt:
            start_frame = min(df[df['trackId'] == tid]['frame'])
            rowIdxs = list(df.index[df['parentTrackId'] == tid])
            if rowIdxs:
                fod = df.iloc[rowIdxs[0],:]['frame']
                pvd_ad = df.loc[rowIdxs[0],'pvd_intensity']
                tod = (fod - start_frame)*10
                division_timesMnt.append(tod)
                pvd_intensityMnt.append(pvd_ad)
                
        for tid in trackIdsWt:
            start_frame = min(df[df['trackId'] == tid]['frame'])
            rowIdxs = list(df.index[df['parentTrackId'] == tid])
            if rowIdxs:
                fod = df.iloc[rowIdxs[0],:]['frame']
                pvd_ad = df.loc[rowIdxs[0],'pvd_intensity']
                tod = (fod - start_frame)*10
                division_timesWt.append(tod)
                pvd_intensityWt.append(pvd_ad)

                
    else:
        division_times = []
        pvd_intensity = []
    
        trackIds = list(df[df['trackId'] != -1]['trackId'])
        totalCells = len(trackIds)

        for tid in trackIds:
            start_frame = min(df[df['trackId'] == tid]['frame'])
            rowIdxs = list(df.index[df['parentTrackId'] == tid])
            if rowIdxs:
                fod = df.iloc[rowIdxs[0],:]['frame']
                pvd_ad = df.loc[rowIdxs[0], 'pvd_intensity']
                tod = (fod - start_frame)*10
                division_times.append(tod)
                pvd_intensity.append(pvd_ad)


    if 'cell_type' in df.columns:
        return division_timesWt, pvd_intensityWt, division_timesMnt, pvd_intensityMnt
    else:
        return division_times, pvd_intensity 
    
    
def errorRate(df):
    
    totalRows = len(df)
    totalErrors = 0
    for i in range(totalRows):
        lid = df.loc[i, 'lineageId']
        if lid == -1:
            totalErrors += 1
                
    return (totalErrors/totalRows)*100
        
