'''The functions listed below take csv files from Ilastik AFTER modification with GFP, RFP and Pvd intensity modifications (see fileManipulation)'''

import pandas as pd
import os 

class CalculateCompetitionMeasures:
    def __init__(self, comp_condition):
        self.comp = comp_condition

    @property
    def df(self):
        path=os.path.join('tracking_files_mod', 'ann_{}_tl-CSV.csv'.format(self.comp))
        try:
            data=pd.read_csv(path)
            return data 
        
        except Exception:
            raise 'Please check the annotated Ilatik files are correct for {}'.format(self.comp)
        

    def error_rate(self):
    
        totalRows = len(self.df)
        totalErrors = 0
        for i in range(totalRows):
            lid = df.loc[i, 'lineageId']
            if lid == -1:
                totalErrors += 1
                    
        return (totalErrors/totalRows)*100
        
    def count_strains(self):
        """
        
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
            n_wt = len(self.df[(self.df['cell_type'] == 'wt') & (self.df['frame'] == i)])
            n_mnt = len(self.df[(self.df['cell_type'] == 'mnt') & (self.df['frame'] == i)])
            time.append(i*10)
            wt.append(n_wt)
            mnt.append(n_mnt)
            
            
        return time, wt, mnt
        
    def cell_count(self):
    
        """
        Given a dataframe of single strain, counts number of cells over time
        
        Parameters
        ----------
        
        df: pandas dataframe of tracking information
        
        Returns
        -------
        
        time: list of time points
        
        cellCount: list of cell counts at each time point
        """

        cellCount = []
        for i in range(30):
            c = len(self.df[self.df['frame'] == i])
            cellCount.append(c)

        time = [i*10 for i in range(30)] 
        
        return time, cellCount
    
    def calculate_growth_rate(self):
        """
        Parameters
        ----------
        cell_counts: a list of cell counts at each time point
        time: a list of time in minutes
        
        Returns
        -------
        r: growth rate of bacteria, r, as the constant of the exponent

        """
        time, cell_counts = self.cell_count()
        from numpy import exp
        A = cell_counts[0]
        ft = cell_counts[15]
        t = time[15]
        
        r = (1/t)*np.log(ft/A)
        return round(r,3)
    
    
    def avg_pvd_intensity(self, frame, condition=None):

        if condition is None:
            pvd_sum = sum(self.df[self.df['frame'] == frame]['pvd_intensity'])
            pvd_avg = pvd_sum/len(self.df[self.df['frame'] == frame])

        else:
            pvd = self.df[(self.df['frame'] == frame) & (self.df['cell_type'] == condition)]
            pvd_sum = sum(pvd['pvd_intensity'])
            pvd_avg = round((pvd_sum/len(pvd['pvd_intensity'])), 2)
        
        return pvd_avg

    def pvd_conc_time(self):
    
        """
        Parameters 
        ----------
        
        df: a pandas dataframe of tracking results from Ilastik
        
        Returns
        -------
        
        time: list of time points
        pvd_wt: average pvd intensity of WT bacteria per time point 
        pvd_mnt: average pvd intensity of mutant bacteria per time point 
        
        OR
        
        time: list of time points
        pvd: average pvd intensity for each time point

        """
        time = []
        
        for i in range(0,30):
            
            time.append(i*10)
            
            if 'vs' in self.comp:

                pvd_wt = [self.avg_pvd_intensity(i, 'wt')]
                pvd_mnt = [self.avg_pvd_intensity(i, 'mnt')]

                return time, pvd_wt, pvd_mnt
            
            else:
                pvd = [self.avg_pvd_intensity(i)]

                return time, pvd
        

    def calculate_division_times_pvd(self, condition=None):

        division_times = []
        pvd_intensity =[]

        if condition is None:
            track_ids = list(self.df[self.df['trackId'] != -1]['trackId'])
        
        else:
            track_ids = list(self.df[(self.df['trackId'] != -1) & (self.df['cell_type'] == condition)]['trackId'])

        for tid in track_ids:
            start_frame = min(self.df[self.df['trackId'] == tid]['frame'])
            rowIdxs = list(self.df.index[self.df['parentTrackId'] == tid])
            if rowIdxs:
                pvd_ad = self.df.loc[rowIdxs[0],'pvd_intensity']
                fod = self.df.iloc[rowIdxs[0],:]['frame'] # frame of division
                tod = (fod - start_frame)*10 # time of division
                division_times.append(tod)
                pvd_intensity.append(pvd_ad)

        return division_times, pvd_intensity