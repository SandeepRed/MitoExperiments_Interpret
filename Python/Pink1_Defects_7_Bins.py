
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This File Clusters the defects and given experimental data

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

###Full clustermap with no experimental data, Figure 4
# All_Par = All_Par1.loc[:, ['basal', 'Leak','maximal', 'delPsi', 'H2O2', 'NADH', 'ATP_m']]
# All_Par_basal = All_Par_basal1.loc[:, ['basal', 'Leak','maximal','delPsi', 'H2O2', 'NADH', 'ATP_m']]


###Pink1
All_Par = All_Par1.loc[:, ['basal', 'Leak','maximal',  'delPsi_FC_Rot', 'delPsi_FC_AA']]
All_Par_basal = All_Par_basal1.loc[:, ['basal', 'Leak','maximal', 'delPsi_FC_Rot', 'delPsi_FC_AA']]


###Parkin
# All_Par = All_Par1.loc[:, ['basal', 'maximal', 'ATP_m']]
# All_Par_basal = All_Par_basal1.loc[:, ['basal', 'maximal', 'ATP_m',]]





#Create a numpy array, same shape as All Par

binned_out = np.zeros((len(All_Par.index), len(All_Par.columns)))


#This loop bins the data by calculating the ranges and calling the binning function

for col_num in range(len(All_Par.columns)):
    
    #Transform the FC's if FC in column name

    if '_FC_' in All_Par.columns[col_num]:
        All_Par.loc[:,All_Par.columns[col_num]] = transform_FC(All_Par.loc[:,All_Par.columns[col_num]])
        All_Par_basal.loc[:,All_Par_basal.columns[col_num]] = transform_FC(All_Par_basal.loc[:,All_Par_basal.columns[col_num]])
        
    
    #Calculate the range from the baseline and calls binning function
   
    
    df = All_Par.iloc[:, col_num]
    lb_10 = All_Par_basal.iloc[:,col_num].quantile(0.4)
    lb_iqr = All_Par_basal.iloc[:,col_num].quantile(0.25)
    lb_out = All_Par_basal.iloc[:,col_num].quantile(0.1) 


    
    ub_10 = All_Par_basal.iloc[:,col_num].quantile(0.6)   
    ub_iqr = All_Par_basal.iloc[:,col_num].quantile(0.75)   
    ub_out = All_Par_basal.iloc[:,col_num].quantile(0.9) 
    
    binned_out[:, col_num] = binning_fn(df,lb_10,lb_iqr, lb_out, ub_10, ub_iqr, ub_out)
    
    
    

#Define all the names for defects
Percent = np.around(np.linspace(0.98, 0.02, 49),2)

Labels = ['DH','C1','C3','C4','F1','Hle','KDyn','KCons']

Row_name_all = [ label  + '_'+ str(x) for label in Labels for x in Percent]

Hle_increase_names = [ 'Hle'  + '_'+ str(x) for x in np.arange(2,11)]

# Multiple_defect_names = ['C1_70_F1_90_Hle_50', 'C1_70_F1_90_Hle_40']



#Create a df from the defects numpy array
All_Par_binned = pd.DataFrame(binned_out, columns= All_Par.columns)


###Change for appropriate data

##Pink1
Test_Pink1 = {"basal": -3,
                  "Leak":-3,
          "maximal": -3,
            "delPsi_FC_Rot": 0,
            "delPsi_FC_AA": -3
            }

##Parkin

# Test_Pink1 = {"basal": 3,
#           "maximal": 0,
#             "ATP_m": -3
#             }


##Pink1
#Append the experimental data and name it
df_with_AlphaSyn_test = All_Par_binned.append(Test_Pink1, ignore_index=True)

###Change for appropriate data
test_array_Pink1 = np.array(['PD_Pink1'])
Row_names = np.concatenate((Row_name_all,Hle_increase_names, test_array_Pink1), axis=0)

df_with_AlphaSyn_test.index = Row_names





sn.set(font_scale = 1.5)

###Change for appropriate data
###Full Clustermap
#sns_clustermap = sn.clustermap(df_with_AlphaSyn_test, col_cluster=False, cmap ='vlag', figsize=(8,20))

###Partial Clustermap

sns_clustermap = sn.clustermap(df_with_AlphaSyn_test, col_cluster=False, yticklabels=True, cmap ='vlag', figsize=(5,120))


sns_clustermap.savefig("Pink1_Full_clustermap_FC.svg")



##Get all the labels from clustermap and save them in a csv file
labels = sns_clustermap.ax_heatmap.yaxis.get_majorticklabels()

#https://stackoverflow.com/questions/73404428/how-to-extract-the-labels-from-sns-clustermap

labels_list = []
number_list = []
for i in labels:
    i = str(i)
    name_start = i.find('\'')+1
    name_end = i.rfind('\'')
    name = i[name_start:name_end]
    number_start = name.rfind('-')+1
    number = name[number_start:]
    labels_list.append(name)
    number_list.append(number)
    
   
labels_dict = {'label': labels_list}
labels_df = pd.DataFrame(labels_dict)
labels_df.to_csv('Pink1_All_defects_Cluster_labels_FC_Full.csv')
