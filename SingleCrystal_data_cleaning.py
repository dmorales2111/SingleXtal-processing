# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 12:26:39 2021

@author: Blackfang
"""

import pandas as pd
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

    if nuc == "Li7":
        # frequency axis is already in ppm for Li7, time to clean data

        df_cleaned = pd.DataFrame(columns = ["Angle", "Real", 'ppm'])

        for f in txtFiles:
            df = pd.read_csv(f, header = None).drop([0,1])
            df.columns = [df.iloc[0,0]]
            df = df.drop(2)

            #splits single column into three by tab delimiter of data
            df_split = df['Real\tImag\tppm'].apply(lambda x: pd.Series(x.split('\t')))
            df_split.columns = ['Real','Imag','ppm']
            df_split = df_split.drop(columns = ['Imag'])
            df_split = df_split.apply(lambda x: pd.to_numeric(x, errors = 'coerce'))

            #find max intensity and peak position for selected spectra, and add it to existing table along with angle
            maxInt, peakPosition = df_split.loc[df2.idxmax()["Real"],:]
            angles = np.linspace(0,180,13)
            for a in np.char.mod('%d',angles):
                currentAngle = a if a in f
            df_cleaned.append([float(currentAngle),maxInt,peakPosition])
        else:
        #insert work for P31 here
