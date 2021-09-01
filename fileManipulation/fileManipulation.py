'''This script manipulates tracking files produced by Ilastik, by adding fluorescence data for each cell detected in the tracking step'''


import os
import re
import pandas as pd
from tifffile import imread
from matplotlib import pyplot as plt
import numpy as np

def create_circular_mask(h, w, center, radius):
  '''Creates a circular area within a bacterium in which to measure fluorescence intensity
  Parameters
  ----------
  
  h = height of the input image
  w = width of the input image
  center = coordinates of bacterial centers within the image
  radius = radius of circular mask generated at the center
  
  Returns
  -------
  mask = a numpy array of dimension h*w 
  
  '''

    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    mask = dist_from_center <= radius
    return mask

os.chdir('~working_directory')
tracking_files = os.listdir('~path_to_tracking_files')


def get_folder_list(cwd, re_tag):
  
  '''Retreives all the folders that contain images of competition assays for each condition
  
  Parameters
  ----------
  
  cwd = current directory
  re_tag = strings contained within all folder names
 
  Returns
  -------
  
  vid_folders = list of all folders names that contain microscopy images of competition assays
  
  '''
    
    p = re.compile(r'(%s)'%re_tag)
    f = os.listdir(cwd)
    vid_folders = []

    for folder in f:
        if p.search(folder):
            vid_folders.append(folder)
    return vid_folders

vid_folders = get_folder_list(os.getcwd(), 'Fpv|WT')

for folder in vid_folders:
    
    for file in tracking_files:
        m = re.compile(r'^(%s)'%folder)
        if m.search(file): # finding the tracking CSV for the current folder 
            csv_path = 'tracking_results/'+file
            break
        else:
            continue
    tracking_data = pd.read_csv(csv_path, delimiter = ',') # opening the csv
    tracking_data['pvd_intensity'] = 0 # creating new pvd_intensity column
    pvd_folder_contents = os.listdir(folder + '/Pvd_corrected_2') # 
    
    for i in range(0,30):
        row_idxs = tracking_data.index[tracking_data['frame'] == i]
        for img in pvd_folder_contents:
            if re.search(r'T=(%s)\b'%i, img):
                pvd_img_path = folder+'/Pvd_corrected_2/'+img
                break
            else:
                continue
                
        imagePvd = imread(pvd_img_path)
        h, w = imagePvd.shape[0], imagePvd.shape[1]
        
        for idx in row_idxs:
            center = tuple(tracking_data.loc[idx, ['Center_of_the_object_0', 'Center_of_the_object_1']])
            mask = create_circular_mask(h,w,center = center, radius = 5)
            pvd_intensity = sum(imagePvd[mask])
            tracking_data.loc[idx, 'pvd_intensity'] = pvd_intensity
        
    if folder.startswith('vs'): # if folder contains competition assay

        tracking_data['cell_type'] = ''
        row_idxs2 = tracking_data.index[tracking_data['frame'] == 0]
        co_ch = {'wt':[], 'mnt':[]}

        gpf_folder = folder+"/GFP/"
        gfp_folder_contents = os.listdir(gpf_folder)
        for gfp in gfp_folder_contents:
            if re.search(r'T=0\b', gfp):
                gfp_img_path = folder+'/GFP/'+gfp
                break

        rfp_folder = folder+'/RFP/'
        rfp_folder_contents = os.listdir(rfp_folder)
        for rfp in rfp_folder_contents:
            if re.search(r'T=0\b', rfp):
                rfp_img_path = folder+'/RFP/'+rfp
                break

        imageGfp = imread(gfp_img_path)
        imageRfp = imread(rfp_img_path)

        for idx2 in row_idxs2:
            lineageId  = tracking_data.loc[idx2, 'lineageId']
            if lineageId != -1:

                center = tuple(tracking_data.loc[idx2, ['Center_of_the_object_0', 'Center_of_the_object_1']])
                mask = create_circular_mask(h,w,center = center, radius = 5)
                gfp_intensity, rfp_intensity = sum(imageGfp[mask]), sum(imageRfp[mask])
                if rfp_intensity > gfp_intensity:
                    co_ch['mnt'].append(lineageId)
                    tracking_data.loc[idx2,'cell_type'] = 'mnt'
                if gfp_intensity > rfp_intensity:
                    co_ch['wt'].append(lineageId)
                    tracking_data.loc[idx2, 'cell_type'] = 'wt'

        for n in range(row_idxs2[-1]+1, len(tracking_data)):

            if tracking_data.loc[n, 'lineageId'] in co_ch['wt']:
                tracking_data.loc[n, 'cell_type'] = 'wt'
            elif tracking_data.loc[n, 'lineageId'] in co_ch['mnt']:
                tracking_data.loc[n, 'cell_type'] = 'mnt'
            else:
                continue
    
                
    tracking_data.to_csv(csv_path, index = False)
            
