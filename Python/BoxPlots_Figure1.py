#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
This file creates boxplots with appropriate bins
Figure 1
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



C1 = pd.read_csv('data/C1_F1_Hle_50Percent_DelPsi_AA_FC.csv', header=None)




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




# df = All_Par_basal.iloc[:,7]
# df_C1 = C1.iloc[:,1]
           



df = All_Par_basal.loc[:, 'delPsi_AA']/All_Par_basal.loc[:,'delPsi']
#df_C1 = C1.iloc[:, 1]/C1.iloc[:,0]

df = transform_FC(df)

#df_C1 = All_Par_basal.loc[:, 'delPsi_Rot']/All_Par_basal.loc[:,'delPsi']
df_C1 = C1.iloc[:,1]

df_C1 = transform_FC(df_C1)



lb_10 = df.quantile(0.4)
lb_iqr = df.quantile(0.25)
lb_out = df.quantile(0.1) 


 
ub_10 = df.quantile(0.6)   
ub_iqr = df.quantile(0.75)   
ub_out = df.quantile(0.9) 


 
binned_out = binning_fn(df,lb_10,lb_iqr, lb_out, ub_10, ub_iqr, ub_out)




binned_out_C1 = binning_fn(df_C1,lb_10,lb_iqr, lb_out, ub_10, ub_iqr, ub_out)


# binned_out = binning_fn(C1.iloc[:,2],lb_10,lb_iqr, lb_out, ub_10, ub_iqr, ub_out)

 

 
 
test_df = pd.DataFrame([df,df_C1]).T
test_df.columns = ['Wild Type', 'F1 0.5']

plt.figure(figsize=(6,8))
plt.rc('font', family='sans-serif')
plt.rc('xtick', labelsize='x-small')
plt.rc('ytick', labelsize='x-large')
# delPsi_DH_30.boxplot()

# plt.savefig('delPsi_DH_60_Percent_cells.svg', dpi = 300)

meanlineprops = dict(linestyle='--', linewidth=2, color='black')
# sn.boxplot(data = C1.iloc[:,2]*convert*mol2pmol, color = 'white', width = 0.5, linewidth=3, showmeans=True, meanline = True, meanprops=meanlineprops)

sn.boxplot(data = df, color = 'white', width = 0.5, linewidth=3, showmeans=True, meanline = True, meanprops=meanlineprops)
sn.boxplot(data = test_df, color = 'white', width = 0.5, linewidth=3, showmeans=True, meanline = True, meanprops=meanlineprops)


# sn.violinplot(data = delPsi_DH_30, color = 'white')

pal = sn.color_palette("vlag")
pal.insert(3,'white')

bins_all =np.arange(-3,4)
bins_pal_dict = {bins_all[i]: pal[i] for i in range(len(pal))}


unique_bins= np.unique(binned_out)
bins = np.sort(np.asarray(unique_bins[0:-1], dtype=int))

for i in bins:
    
    
    #df = C1.loc[binned_out == str(i)]
    data_binned = df.loc[binned_out == str(i)]
    data_binned_C1 = df_C1.loc[binned_out_C1 == str(i)]

    
    binned_test_df = pd.DataFrame([data_binned,data_binned_C1]).T
    binned_test_df.columns = ['Wild Type', 'F1 0.5']

    
    #array_df = data_binned.T.to_numpy()
    # sn.stripplot(data = df.iloc[:,2]*convert*mol2pmol)

    sn.stripplot(data = binned_test_df, color = bins_pal_dict[i], palette=[bins_pal_dict[i]] * 3)
    



# sn.boxplot(data = C1.iloc[2,:]*convert*mol2pmol, color = 'white', width = 0.5, linewidth=3, showmeans=True, meanline = True, meanprops=meanlineprops)
# sn.stripplot(data = C1.iloc[2,:]*convert*mol2pmol, color = 'blue', alpha = 1)


#All_Par_basal1.loc[All_Par_basal1.iloc[:,2]*convert*mol2pmol >7, 2]
# All_Par_basal1 = All_Par_basal1.loc[All_Par_basal1.iloc[:,2]*convert*mol2pmol >7]
# sn.stripplot(data = All_Par_basal1.iloc[:,2]*convert*mol2pmol, color = 'green', alpha = 0.3)

# plt.axhline(y=(All_Par_basal1.iloc[:,2].mean()*convert*mol2pmol + 4.2), color='r', linestyle='--')
# plt.axhline(y=(All_Par_basal1.iloc[:,2].mean()*convert*mol2pmol - 4.2), color='r', linestyle='--')

# plt.axhline(y=0.00029681, color='r', linestyle='--')
# plt.axhline(y=0.0004061, color='r', linestyle='--')

# plt.xlabel('Fraction DH',  fontsize=20)
plt.ylabel(' $\Delta\Psi $  (FC)',  fontsize=18)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.tick_params(axis='both', which='minor', labelsize=20)
plt.tight_layout()
#plt.savefig('C1_50Percent_BoxPlot_bins_BeforeFCTransform.svg', dpi = 300)

#plt.savefig('Hle_50Percent_BoxPlot_bins_BeforeFCTransform.svg', dpi = 300)
#plt.savefig('F1_50Percent_BoxPlot_bins_AfterFCTransform.svg', dpi = 300)

#plt.savefig('delPsi_BoxPlot_bins_Oligo_delPsi_C1_80_FC_withROS.svg', dpi = 300)

