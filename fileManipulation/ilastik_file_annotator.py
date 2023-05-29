import os
import re
import pandas as pd
from tifffile import imread
from matplotlib import pyplot as plt
import numpy as np
import argparse

def parser():
    parser = argparse.ArgumentParser(
                prog='ManipulateTrackingFiles',
                description='parses command line arguments for script')
    parser.add_argument('--working-dir', dest='working_dir', required=False,
                        help='directory to work in, all required tracking files and timelapse videos should be here')

class IlastikFileAnnotator:

    def __init__(self, working_dir, comp_condition, video_folder=None, tracking_file=None):
        self.working_dir = working_dir
        self.comp_condition = comp_condition
        self._video_folder=video_folder
        self._tracking_file=tracking_file


    @staticmethod
    def create_circular_mask(h, w, center, radius):

        Y, X = np.ogrid[:h, :w]
        dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

        mask = dist_from_center <= radius
        return mask
    
    def fluorescent_image_path(self, fluorescence_channel):
        
        """
        Parameters 
        --------------
        @condition: competition condition 
        @fp_color: fluorescence protein color, either 'GFP' or 'RFP'

        Returns 
        ---------------
        Path to first image in timelapse for given condition, in given fluorescence channel
        """
        fluo_folder = os.path.join(self.video_folder, fluorescence_channel)
        fluo_folder_contents = os.listdir(fluo_folder)

        for img in fluo_folder_contents:
            if re.search(r'T=0\b', img):
                return os.path.join(self.video_folder,
                                    fluorescence_channel,
                                    img)

    @property
    def video_folder(self):
        if self._video_folder is None:

            self._video_folder = os.path.join(
                self.working_dir,
                'timelapse',
                self.comp_condition
            )

        return self._video_folder

    @property
    def tracking_file(self):
        if self._tracking_file is None:
            self._tracking_file = os.path.join(self.working_dir, 
                                               'tracking_files', 
                                               '{}_tl_CSV-Table.csv'.format(self.comp_condition)
                                               )
        return self._tracking_file_dir
            

    def annotate_tracking_files(self):

        tracking_data = pd.read_csv(self.tracking_file, delimiter = ',') # reading the csv
        tracking_data['pvd_intensity'] = 0 # creating new pvd_intensity column
        pvd_folder = os.path.join(self.video_folder, 'Pvd_corrected_2') 
        
        for i in range(0,30): # for each frame in the timelapse of comp assay 
        
            for img in pvd_folder:
                if re.search(r'T={timepoint}\b'.format(timepoint=i), img): # identify the corresponding frame/time point (T = 0,1,2...30)
                    pvd_img_path = os.path.join(pvd_folder,img)
                    break
                else:
                    continue
                    
            imagePvd = imread(pvd_img_path) # open the pyoverdine image
            h, w = imagePvd.shape[0], imagePvd.shape[1]
            
            row_idxs = tracking_data.index[tracking_data['frame'] == i]
            for idx in row_idxs:
                center = tuple(tracking_data.loc[idx, ['Center_of_the_object_0', 'Center_of_the_object_1']]) # extract bacterial centre point from csv file
                mask = self.create_circular_mask(h,w,center = center, radius = 5)
                pvd_intensity = sum(imagePvd[mask]) # record fluorescence intensity in cyan channel
                tracking_data.loc[idx, 'pvd_intensity'] = pvd_intensity
            
        if self.comp_condition.startswith('vs'): # if condition is a competition assays

            tracking_data['cell_type'] = ''
            row_idxs2 = tracking_data.index[tracking_data['frame'] == 0]
            co_ch = {'wt':[], 'mnt':[]} # initialise dictionary of wt and mutant lineage ids

            imageGfp = imread(self.fluorescent_image_path(self.video_folder, 'GFP')) # open RFP and GFP images
            imageRfp = imread(self.fluorescent_image_path(self.video_folder, 'RFP'))

            co_ch = {}

            # identify if wt or mnt at t=0
            for idx2 in row_idxs2:
                lineageId  = tracking_data.loc[idx2, 'lineageId']
                if lineageId != -1:
                    
                    center = tuple(tracking_data.loc[idx2, ['Center_of_the_object_0', 'Center_of_the_object_1']])
                    mask = self.create_circular_mask(h,w,center = center, radius = 5)
                    gfp_intensity, rfp_intensity = sum(imageGfp[mask]), sum(imageRfp[mask])
                    if rfp_intensity > gfp_intensity:
                        co_ch[lineageId] = 'mnt'
                        tracking_data.loc[idx2,'cell_type'] = 'mnt'
                    if gfp_intensity > rfp_intensity:
                        co_ch[lineageId] = 'wt'
                        tracking_data.loc[idx2, 'cell_type'] = 'wt'

            # add above annotation to all all daughters of same lineageId
            for n in range(row_idxs2[-1]+1, len(tracking_data)):
                lineageId=tracking_data.loc[n, 'lineageId']

                if lineageId in co_ch.keys():
                    tracking_data.loc[n, 'cell_type'] = co_ch[lineageId]
                else:
                    continue
                    
            tracking_data.to_csv('output/ann_{}_tl-CSV.csv'.format(self.comp_condition), 
                                 index = False)

if __name__ == "__main__": 

args=parser.parse_args()

for comp in os.listdir('{}/timelapse'.format(args.working_dir)):
   ila = IlastikFileAnnotator(args.working_dir,
                              comp)
   ila.annotate_tracking_files()