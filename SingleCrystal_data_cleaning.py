# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 12:26:39 2021

@author: Blackfang
"""

import pandas as pd
import numpy as np
import os
import tkinter as tk
from tkinter import filedialog as fd
from natsort import os_sorted
root = tk.Tk()

#first get the crystal directory
direct = fd.askdirectory(parent = root, title = 'Please select a crystal directory')
print('working with files in ' + direct)
root.destroy()
os.chdir(direct)
#%%
#grab sub-folder names

subdirs = []
for (dirpath, dirnames, filenames) in os.walk(direct):
    subdirs.extend(dirnames)
    break

#ask which nucleus you want to work with
nuc = input('Enter which nucleus you want to work with: ')

for fol in subdirs:
    compileMax(fol,nuc)





#helper functions
def compileMax(folder,nucleus):
    os.chdir(direct + '\\' + folder)
    txtFiles = os_sorted([x for x in os.listdir() if '.txt' in x and nucleus in x])
    df_cleaned = pd.DataFrame(columns = ["Angle", "Real", 'ppm'])

    for f in txtFiles:
        df = pd.read_csv(f, header = None).drop([0,1])
        df.columns = [df.iloc[0,0]]
        df = df.drop(2)
        if nuc == 'P31':
            orig_ref_freq = float(input("Please enter the reference frequency in MHz: "))
            obs_freq = float(input("Please enter transmitter frequency in MHz: "))
            
            #splits single column into three by tab delimiter of data
            df_split = df['Real\tImag\tHz'].apply(lambda x: pd.Series(x.split('\t')))
            df_split = df_split.apply(lambda x: pd.to_numeric(x, errors = 'coerce'))
            
            #converts Hz to ppm using ref freq    
            df_split.iloc[:,2] = (df_split.iloc[:,2] + 1e6*(obs_freq - orig_ref_freq))/orig_ref_freq 
            
        #otherwise, if nucleus is Li7, data is already in ppm, dont need to convert and can continue
        else:
            df_split = df['Real\tImag\tppm'].apply(lambda x: pd.Series(x.split('\t')))
            df_split = df_split.apply(lambda x: pd.to_numeric(x, errors = 'coerce'))
       
        df_split.columns = ['Real','Imag','ppm']
        df_split = df_split.drop(columns = ['Imag'])

        #find max intensity and peak position for selected spectra, and add it to existing table along with angle
        maxInt, peakPosition = df_split.loc[df_split.idxmax()["Real"],:]
        angles = np.linspace(0,180,13)
            
        currentAngle = [x for x in np.char.mod('%d',angles) if x in f]
        new_row = pd.Series([float(currentAngle[0]), maxInt, peakPosition], index = df_cleaned.columns)
        df_cleaned = df_cleaned.append(new_row, ignore_index = True)
      
        #insert work for P31 here
        
        