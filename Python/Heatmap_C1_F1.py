#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file creates a clustermap for C1 with appropriate bioenregetic variable
Figure 3
"""

import numpy as np
import pandas as pd

import seaborn as sn


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
                     



#Load all variables for defects and baseline
Seahorse_Par = pd.read_csv('data/2PercentAllDefects_exp1_2_3_50Sims_withROS_FC.csv', header=None)
Seahorse_Par_basal = pd.read_csv('data/WildType1000sims_exp1_2_3_50Sims_withROS_FC.csv', header=None)

#Load all variables for both proton leak increase

Hle_Increase_All_Par = pd.read_csv('data/HleIncrease_exp1_2_3_50Sims_withROS_FC.csv', header=None)



#Concat all defects with proton leak increase

All_Par1 = pd.concat([Seahorse_Par, Hle_Increase_All_Par], axis = 0)

#Select One variable for each OCR metric 2: Basal, 3: Leak, 7: Max 
All_Par1 = All_Par1.iloc[:, np.r_[2, 4, 7, 14:75]]
All_Par_basal1 = Seahorse_Par_basal.iloc[:, np.r_[2, 4, 7, 14:75]]





#Calculate Spare Capacity

All_Par1 = All_Par1.assign(Spare_capacity = All_Par1.iloc[:,2] - All_Par1.iloc[:,0])
All_Par_basal1 = All_Par_basal1.assign(Spare_capacity = All_Par_basal1.iloc[:,2] - All_Par_basal1.iloc[:,0])


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
                'delPsi_FC_Oligo','delPsi_FC_Oligo_FCCP','ATP_m_FC_Oligo','ATP_m_FC_Oligo_FCCP',
                'NADH_m_FC_Oligo','NADH_m_FC_Oligo_FCCP','delPsi_FC_Rot','delPsi_FC_Rot_Oligo','ATP_m_FC_Rot','ATP_m_FC_Rot_Oligo',
                'NADH_m_FC_Rot','NADH_m_FC_Rot_Oligo','delPsi_FC_AA','delPsi_FC_AA_Oligo','ATP_m_FC_AA','ATP_m_FC_AA_Oligo',
                'NADH_m_FC_AA','NADH_m_FC_AA_Oligo','delPsi_FC_FCCP','delPsi_FC_FCCP_Rot','ATP_m_FC_FCCP','ATP_m_FC_FCCP_Rot',
                'NADH_m_FC_FCCP','NADH_m_FC_FCCP_Rot','Spare_capacity']






All_Par1.columns = Column_names
All_Par_basal1.columns = Column_names




###Change for appropriate data


# #OCR
# All_Par = All_Par1.loc[:,['basal', 'Leak', 'maximal', 'Spare_capacity']]
# All_Par_basal = All_Par_basal1.loc[:,['basal', 'Leak', 'maximal', 'Spare_capacity']]


#delPsi
# All_Par = All_Par1.iloc[:,[3]+list(np.arange(8,40,4)) +[40,41,46,47,52,53,58,59]]
# All_Par_basal = All_Par_basal1.iloc[:,[3]+list(np.arange(8,40,4))+[40,41,46,47,52,53,58,59]]


# #ATPm
# All_Par = All_Par1.iloc[:,[4]+list(np.arange(9,40,4))+[i+2 for i in [40,41,46,47,52,53,58,59]]]
# All_Par_basal = All_Par_basal1.iloc[:,[4]+list(np.arange(9,40,4))+[i+2 for i in [40,41,46,47,52,53,58,59]]]



#NADHm
# All_Par = All_Par1.iloc[:,[5]+list(np.arange(10,40,4))+[i+4 for i in [40,41,46,47,52,53,58,59]]]
# All_Par_basal = All_Par_basal1.iloc[:,[5]+list(np.arange(10,40,4))+[i+4 for i in [40,41,46,47,52,53,58,59]]]



# #H2O2
# All_Par = All_Par1.iloc[:,[7]]
# All_Par_basal = All_Par_basal1.iloc[:,[7]]


# #All_outputs
All_Par = All_Par1
All_Par_basal = All_Par_basal1



#Append baseline row
All_Par = All_Par.append(All_Par_basal.median(), ignore_index=True)




#Create a numpy array, same shape as All Par

binned_out = np.zeros((len(All_Par.index), len(All_Par.columns)))


#This loop bins the data by calculating the ranges and calling the binning function

for col_num in range(len(All_Par.columns)):
    
    ###Change, Comment out if plotting raw values
    
    #Transform the FC's if FC in column name

    # if '_FC_' in All_Par.columns[col_num]:
    #     All_Par.loc[:,All_Par.columns[col_num]] = transform_FC(All_Par.loc[:,All_Par.columns[col_num]])
    #     All_Par_basal.loc[:,All_Par_basal.columns[col_num]] = transform_FC(All_Par_basal.loc[:,All_Par_basal.columns[col_num]])
        
    
    #Calculate the range from the baseline and calls binning function
   
    
    df = All_Par.iloc[:, col_num]
    lb_10 = All_Par_basal.iloc[:,col_num].quantile(0.4)
    lb_iqr = All_Par_basal.iloc[:,col_num].quantile(0.25)
    lb_out = All_Par_basal.iloc[:,col_num].quantile(0.1) 


    
    ub_10 = All_Par_basal.iloc[:,col_num].quantile(0.6)   
    ub_iqr = All_Par_basal.iloc[:,col_num].quantile(0.75)   
    ub_out = All_Par_basal.iloc[:,col_num].quantile(0.9) 
    
    ###Change
    ###For binned
    #binned_out[:, col_num] = binning_fn(df,lb_10,lb_iqr, lb_out, ub_10, ub_iqr, ub_out)
    ###To Get raw values
    binned_out[:, col_num] = df

    

#Define all the names for defects
Percent = np.around(np.linspace(0.98, 0.02, 49),2)

Labels = ['DH','C1','C3','C4','F1','Hle','KDyn','KCons']

Row_name_all = [ label  + '_'+ str(x) for label in Labels for x in Percent]

Hle_increase_names = [ 'Hle'  + '_'+ str(x) for x in np.arange(2,11)]

All_Par_binned = pd.DataFrame(binned_out, columns= All_Par.columns)



Row_names = np.concatenate((Row_name_all,Hle_increase_names,['WT']), axis=0)

All_Par.index = Row_names


All_Par_binned = pd.DataFrame(binned_out, columns= All_Par.columns, index = Row_names)




###Change for appropriate data
##CI , 49:98
DH_data = All_Par_binned.iloc[np.r_[401,49:98], :]





#Define custom color with centre at 1
#https://stackoverflow.com/questions/52626103/custom-colormap

import matplotlib.colors

norm = matplotlib.colors.Normalize(0,2)
colors = [[norm(0.0), "darkblue"],
          [norm(1.0), "white"],
          [norm( 2.0), "red"]]

cmap_rwb = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)









# sns_clustermap = sn.clustermap(DH_data.iloc[:,2], center=1,row_cluster = False, col_cluster=False,  figsize=(2, 10), cmap =cmap_rwb)
# sns_clustermap.savefig("C1_delPsi_rawvals_FC_test.svg")


###Change for appropriate data

#OCR

# DH_OCR_basal_Leak = DH_data[['maximal', 'basal', 'Leak']]             
# sns_clustermap = sn.clustermap(DH_OCR_basal_Leak, xticklabels=True, col_cluster=False, row_cluster = False,  figsize=(4, 12))
# sns_clustermap.savefig("Heatmap_F1_OCR_basal_Leak.svg")


#ATP_m_FC_AA
#Centre color palette change depending on binning(0 centre, vlag) or raw values
#sns_clustermap = sn.clustermap(DH_data['ATP_m_FC_AA'], center=0,row_cluster = False, col_cluster=False,  figsize=(3, 10), cmap ='vlag')
#sns_clustermap = sn.clustermap(DH_data['ATP_m_FC_AA'], center=1,row_cluster = False, col_cluster=False,  figsize=(3, 10), cmap =cmap_rwb)

#sns_clustermap.savefig("Heatmap_C1_ATP_AAQuantitative_1.svg")



#sns_clustermap = sn.clustermap(DH_data['delPsi_FC_Oligo'], center=0,row_cluster = False, col_cluster=False,  figsize=(2, 10), cmap ='vlag')
sns_clustermap = sn.clustermap(DH_data['delPsi_FC_Oligo'], center=1,row_cluster = False, col_cluster=False,  figsize=(3, 10), cmap =cmap_rwb)

sns_clustermap.savefig("Heatmap_C1_delPsi_Oligo_quant.svg")


# sns_clustermap = sn.clustermap(DH_data['NADH_m_FC_Rot'], center=0,row_cluster = False, col_cluster=False,  figsize=(2, 10), cmap ='vlag')
# sns_clustermap = sn.clustermap(DH_data['NADH_m_FC_Rot'], center=1,row_cluster = False, col_cluster=False,  figsize=(2, 10), cmap =cmap_rwb)

# sns_clustermap.savefig("Heatmap_C1_NADH_m_Rot_binned.svg")



