#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This file plots OCR for Pink1

"""


import pandas as pd
import matplotlib.pyplot as plt
import matplotlib


#Define the conversion units
sec2min = 60   # to convert from /s to /min

#to convert from /l of mito to total mito volume 
#Ward 2007 CGNs mitochondrial volume = 1.76 +/- 0.15 e-14 litres
litre2mito = 5e-14 

#to convert from 1 cell to # cells/well 
#300,000 = Pierre's seeding density
cell2wells = 300000  

#to convert # cells to ug protein 
#40-50 ug protein in 300000 cells / well (Pierre's measurements)
cells2protein = 45   

convert = (sec2min * litre2mito * cell2wells) / cells2protein;

mol2pmol = 1e12





    

#Read the file, Time OCR has timepoints for the OCR
Time_OCR = pd.read_csv('data/WT_NADH096_Hle07_OCR.csv', header=None)
OCR = pd.read_csv('data/C1_F1_Hle_defectsPink1_OCR.csv', header=None)

#Convert
OCR = OCR*convert*mol2pmol




Labels = ['WT','C1_0.7','F1_0.9','Hle_0.4', 'C1_0.7_F1_0.9','C1_0.7_F1_0.9_Hle_0.4']



matplotlib.font_manager.FontProperties(family='Calibri', weight='bold')

plt.figure(figsize=(10,7))
plt.rc('font', family='sans-serif')
plt.rc('xtick', labelsize='x-small')
plt.rc('ytick', labelsize='x-large')
#plt.grid(False)


for i in range(len(Labels)):
    
    plt.plot(Time_OCR.iloc[:, 0], OCR.iloc[:,i], label=Labels[i], linewidth=4.0)
    



plt.legend(loc='best', fontsize = 16)


plt.ylabel('OCR (pmol/min/ug protein)', fontsize=20, multialignment='center')
plt.xlabel('Time (min)', fontsize=20)

#plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.xlim(1, 70)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.tick_params(axis='both', which='minor', labelsize=18)
plt.tight_layout()
plt.savefig('OCR_Pink1_all_defects_Combinations.svg', dpi = 300)

