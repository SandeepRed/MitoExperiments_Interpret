#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This File Clusters the defects and given experimental data

"""

import numpy as np
import pandas as pd

import seaborn as sn
import matplotlib.pyplot as plt


def binning_fn(df,lb_10,lb_iqr, lb_out, ub_10, ub_iqr, ub_out):
    

    '''
    Inputs: a df of experimental data and baseline ranges(lower and upperbounds)
    Returns: Value, corresponding to the bin
    '''
    
    binned_vals = np.where((df >=  lb_10) & (df <= ub_10),  
                        0,
                        np.where((df >  ub_10) & (df <= ub_iqr),  
                        1,
                        np.where( (df > ub_iqr) & (df <= ub_out),
                        2,
                        np.where( (df > ub_out),
                        3,
                        np.where((df <  lb_10) & (df >= lb_iqr),  
                        -1,
                        np.where( (df <  lb_iqr) & (df >= lb_out),
                        -2,
                        np.where( (df <  lb_out),
                        -3,
                        'other')))))))
    
    return binned_vals


def transform_FC(col_fc):
    '''
    This function flips FC between 0 and 1 to 1-FC
    Inputs: columns with fold changes
    Outputs: Pandas series with tranformed fold changes
    '''
    
    FC = np.zeros(len(col_fc))

    
    i=0
    
    for x in col_fc:
        if x<=1:
            FC[i] = 1-x
            i =i+1
        else:
           FC[i] = x
           i=i+1
    return pd.Series(FC)
           



#Load all variables for baseline
Seahorse_Par_basal = pd.read_csv('data/WildType1000sims_exp1_2_3_50Sims_withROS_FC.csv', header=None)


#Select One variable for each OCR metric 2: Basal, 3: Leak, 7: Max 
Seahorse_Par_basal1 = Seahorse_Par_basal.iloc[:, np.r_[2, 4, 7, 14:75]]




#Calculate Spare Capacity
All_Par_basal = Seahorse_Par_basal1.assign(Spare_capacity = Seahorse_Par_basal1.iloc[:,2] - Seahorse_Par_basal1.iloc[:,0])



#Define all Column Names
Column_names = ['basal', 'Leak', 'maximal', 'delPsi', 'ATP_m', 'NADH','ATP_c', 'H2O2',
                'delPsi_Oligo', 'ATP_m_Oligo', 'NADH_m_Oligo',  'ATP_c_Oligo',  
                'delPsi_Oligo_FCCP', 'ATP_m_Oligo_FCCP', 'NADH_m_Oligo_FCCP',  'ATP_c_Oligo_FCCP',
                'delPsi_Rot', 'ATP_m_Rot', 'NADH_m_Rot', 'ATP_c_Rot',
                'delPsi_Rot_Oligo', 'ATP_m_Rot_Oligo', 'NADH_m_Rot_Oligo', 'ATP_c_Rot_Oligo', 
                'delPsi_AA', 'ATP_m_AA', 'NADH_m_AA', 'ATP_c_AA', 
                'delPsi_AA_Oligo', 'ATP_m_AA_Oligo', 'NADH_m_AA_Oligo', 'ATP_c_AA_Oligo', 
                'delPsi_FCCP', 'ATP_m_FCCP','NADH_FCCP', 'ATP_c_FCCP',
                'delPsi_FCCP_Rot', 'ATP_m_FCCP_Rot','NADH_FCCP_Rot', 'ATP_c_FCCP_Rot', 
                'delPsi_FC_Oligo',' delPsi_FC_Oligo_FCCP',' ATP_m_FC_Oligo',' ATP_m_FC_Oligo_FCCP',
                'NADH_m_FC_Oligo',' NADH_m_FC_Oligo_FCCP','delPsi_FC_Rot',' delPsi_FC_Rot_Oligo',' ATP_m_FC_Rot',' ATP_m_FC_Rot_Oligo',
                'NADH_m_FC_Rot',' NADH_m_FC_Rot_Oligo','delPsi_FC_AA',' delPsi_FC_AA_Oligo',' ATP_m_FC_AA',' ATP_m_FC_AA_Oligo',
                'NADH_m_FC_AA',' NADH_m_FC_AA_Oligo','delPsi_FC_FCCP',' delPsi_FC_FCCP_Rot',' ATP_m_FC_FCCP',' ATP_m_FC_FCCP_Rot',
                'NADH_m_FC_FCCP',' NADH_m_FC_FCCP_Rot','Spare_capacity']


All_Par_basal.columns = Column_names



###Change the datasets here

#Read appropriate file, with either delPsi or delPsiFC for the combined defects
All_data = pd.read_csv('data/C1_F1_Hle_defectsPink1_delPsi_Rot_AA.csv', header=None)



###Change for appropriate data
#In the all data file first 6 columns are Rotenone
All_Par = All_data.iloc[:,np.r_[0:6]]

#Second 6 columns are Antimycin
#All_Par = All_data.iloc[:,np.r_[6:12]]

#Name the defect combinations of those 6 columns
All_Par.columns = ['WT','C1_0.7','F1_0.9','Hle_0.4', 'C1_0.7_F1_0.9','C1_0.7_F1_0.9_Hle_0.4']





#Create a dataframe, same shape as All Par
binned_out = np.zeros((len(All_Par.index), len(All_Par.columns)))
binned_out = pd.DataFrame(binned_out)


###Change for appropriate data
###delPsi_Rot or delPsi_AA or delPsiFC

#Create a dataframe, same shape as All PPar
col_num_basal = 'delPsi_Rot'

All_Par_basal.loc[:,col_num_basal] = transform_FC(All_Par_basal.loc[:,col_num_basal])


#This loop bins the data by calculating the ranges and calling the binning function

for col_num in range(len(All_Par.columns)):
    
    #Transform the FC's if FC in column name
    
    if '_FC_' in All_Par.columns[col_num]:
        All_Par.loc[:,All_Par.columns[col_num]] = transform_FC(All_Par.loc[:,All_Par.columns[col_num]])
        
    
    
    #Calculate the range from the baseline and calls binning function

    df = All_Par.iloc[:, col_num]
    lb_10 = All_Par_basal.loc[:,col_num_basal].quantile(0.4)
    lb_iqr =All_Par_basal.loc[:,col_num_basal].quantile(0.25)
    lb_out =All_Par_basal.loc[:,col_num_basal].quantile(0.1) 


    
    ub_10 =All_Par_basal.loc[:,col_num_basal].quantile(0.6)   
    ub_iqr =All_Par_basal.loc[:,col_num_basal].quantile(0.75)   
    ub_out =All_Par_basal.loc[:,col_num_basal].quantile(0.9) 
    
    binned_out.iloc[:, col_num] = binning_fn(df,lb_10,lb_iqr, lb_out, ub_10, ub_iqr, ub_out)
 

 

#Name the defect combinations of those 6 columns

 
binned_out.columns = ['WT','C1_0.7','F1_0.9','Hle_0.4', 'C1_0.7_F1_0.9','C1_0.7_F1_0.9_Hle_0.4']


#Create figure with appropriate fonts and figsize
plt.figure(figsize=(10,7))
plt.rc('font', family='sans-serif')
plt.rc('xtick', labelsize='x-small')
plt.rc('ytick', labelsize='x-large')


#Define meanline properties
meanlineprops = dict(linestyle='--', linewidth=2, color='black')

#Create boxplot with meanline
sn.boxplot(data = All_Par, color = 'white', width = 0.5, linewidth=3, showmeans=True, meanline = True, meanprops=meanlineprops)



#Define pallate, vlag has 6 combinations blue to red, add white in between for 0
pal = sn.color_palette("vlag")
pal.insert(3,'white')

#Define all bins and associate appropriate colour
bins_all =np.arange(-3,4)
bins_pal_dict = {bins_all[i]: pal[i] for i in range(len(pal))}

#Get all unique bins and sort them
unique_bins= np.unique(binned_out)
bins = np.sort(np.asarray(unique_bins[0:-1], dtype=int))



#Create a list which will hold a list of values for each defect
list_test = []


for i in bins:
    
    #For the present bin loop, append the list of values for each defect
    #In All_Par
    for j in range(len(All_Par.columns)):
        list_test.append(All_Par.iloc[:,j].loc[binned_out.iloc[:,j] == str(i)])

    #Create a df of the list of lists and reset the list for next iteration
    binned_test_df = pd.DataFrame(list_test).T
    list_test = []

    binned_test_df.columns = ['WT','C1_0.7','F1_0.9','Hle_0.4', 'C1_0.7_F1_0.9','C1_0.7_F1_0.9_Hle_0.4']


    #Create the stripplot with the color from pallete
    sn.stripplot(data = binned_test_df, color = bins_pal_dict[i], palette=[bins_pal_dict[i]] * 3)
    



#Other plot labels and save the plot

plt.xticks(list(range(len(All_Par.columns))),
               labels=['WT','C1_0.7','F1_0.9','Hle_0.4', 
                        'C1_0.7\nF1_0.9','C1_0.7\nF1_0.9\nHle_0.4'])

plt.ylabel(' $\Delta\Psi $  (Rotenone)',  fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.tick_params(axis='both', which='minor', labelsize=20)

plt.tight_layout()

plt.savefig('Pink1_defectCombo_FCdelPsi_Rotenone.svg', dpi = 300)
















