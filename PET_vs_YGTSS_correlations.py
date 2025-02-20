#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 14:59:11 2023

@author: gkotsoulias
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 12:04:38 2023

@author: gkotsoulias
"""

from scipy.io import savemat
import mat73
import scipy.io
import ants
import seaborn as sb
from matplotlib import pyplot as plt
import pandas as pd
import mat73
import scipy.stats as scs
import numpy as np


#Patients_7T_codes =  [4,6,8,12,11,14,7,15,16,18,20,21,17,19,22,25,23,gtsPp017(no 7T)]

"""
Create all needed YGTSS vectors------------------------------------------------
"""
YGTSS_Motor_3T =           [13,15,16,9,19,14,17,19,21,15,10,18,17,14,18,18,17,13]
YGTSS_Vocal_3T =           [14,12,8,7,10,10,0,18,11,9,6,13,13,16,22,14,5,9]
YGTSS_SubjImpairment_3T =  [10,0,20,0,20,20,10,40,10,20,30,30,10,10,30,30,20,10]

YGTSS_Full_3T = []
for i in range(len(YGTSS_Motor_3T)):
    YGTSS_Full_3T.append(YGTSS_Motor_3T[i] + YGTSS_Vocal_3T[i])

YGTSS_Total_3T = []
for i in range(len(YGTSS_Motor_3T)):
    YGTSS_Total_3T.append(YGTSS_Motor_3T[i] + YGTSS_Vocal_3T[i]+ YGTSS_SubjImpairment_3T[i])
    
GTS_BPnd_YGTSS  = np.zeros((9,18))
mat1 = scipy.io.loadmat('/---/---/PET_BPnd_Caudate_stats_FullCohorts.mat')
GTS_BPnd_YGTSS[0,:] = mat1['GTS_BPnd']
mat2 = scipy.io.loadmat('/---/---/PET_BPnd_Pallidum_stats_FullCohorts.mat')
GTS_BPnd_YGTSS[1,:] = mat2['GTS_BPnd']  
mat3 = scipy.io.loadmat('/---/---/PET_BPnd_Putamen_stats_FullCohorts.mat')
GTS_BPnd_YGTSS[2,:] = mat3['GTS_BPnd']    
GTS_BPnd_YGTSS[3,:] = (mat1['GTS_BPnd']+mat3['GTS_BPnd'])/2      #STR   
mat = scipy.io.loadmat('/---/---/PET_BPnd_Thalamus_stats_FullCohorts.mat')
GTS_BPnd_YGTSS[4,:] = mat['GTS_BPnd']          
GTS_BPnd_YGTSS[5,:] =    YGTSS_Motor_3T        
GTS_BPnd_YGTSS[6,:] =    YGTSS_Vocal_3T
GTS_BPnd_YGTSS[7,:] =    YGTSS_Full_3T        
GTS_BPnd_YGTSS[8,:] =    YGTSS_Total_3T   

final_data =  np.delete(GTS_BPnd_YGTSS, [], axis=1) #no exclusions

def calculate_pvalues(df):
    dfcols = pd.DataFrame(columns=df.columns)
    pvalues = dfcols.transpose().join(dfcols, how='outer')
    for r in df.columns:
        for c in df.columns:
            tmp = df[df[r].notnull() & df[c].notnull()]
            pvalues[r][c] = round(scs.pearsonr(tmp[r], tmp[c])[1], 4)
    return pvalues

# Print heatmap for pearsons coeff.
final_data_T = final_data.transpose()
frame_patients = pd.DataFrame(data=final_data_T,index=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18], columns=['Caudate','Pallidum','Putamen','STR','Thalamus','YGTSS_motor','YGTSS_vocal','YGTSS_full','YGTSS_total'])  
corr = frame_patients.corr(method = 'pearson')
corr
pvals = calculate_pvalues(frame_patients)
pvals
plt.figure(figsize=(20,20), dpi =300)
sb.heatmap(corr,annot=True, square = 'True',cmap="vlag")
imname = '/data/hu_gkotsoulias/Desktop/YGTSS_vs_BPnd.png'
plt.savefig(imname)      

plt.figure(figsize=(20,20), dpi =300)
sb.heatmap(corr,annot=pvals, square = 'True',cmap="rocket_r",vmin=0, vmax=1)
imname = '/data/hu_gkotsoulias/Desktop/YGTSS_vs_BPnd_pvals.png'
plt.savefig(imname)   

data = pd.DataFrame(data={'YGTSS Motor': final_data[0,:], 'BPnd': final_data[5,:] })  
sb.set_style("darkgrid")
g=sb.lmplot(x="BPnd", y="YGTSS Motor",palette='Set1',height=10, data=data)
sb.set_context("notebook", font_scale=1.5, rc={"lines.linewidth": 1.5})
g = (g.set_axis_labels("YGTSS Motor","CAUDATE BPnd").set(xlim=(8,22), ylim=(1.48, 2.73)))
imname = '/---/---/CAUDATE_BPnd_Vs_BPnd_correlationGraph.png'
plt.savefig(imname)         
    
data1 = pd.DataFrame(data={'x': final_data[3,:], 'y': final_data[5,:] })  #,'Weights': Weights
p = sb.regplot(x="x", y="y",data=data1)

#calculate slope and intercept of regression equation
slope, intercept, r, g, sterr = scipy.stats.linregress(x=p.get_lines()[0].get_xdata(),
                                                       y=p.get_lines()[0].get_ydata())
print(intercept, slope)

    
    
