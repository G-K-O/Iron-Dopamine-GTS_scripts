from scipy.io import savemat
import mat73
import scipy.io
import seaborn as sb
from matplotlib import pyplot as plt
import pandas as pd
import mat73
import scipy.stats as scs
import numpy as np

#dataset is constructed manually in matlab, here it is loaded and used for stats
mat = scipy.io.loadmat('/---/---/serumVsBPnd_fullCohorts_finalOrdered.mat')
alldata = mat['alldata']

def calculate_pvalues(df):
dfcols = pd.DataFrame(columns=df.columns)
pvalues = dfcols.transpose().join(dfcols, how='outer')
for r in df.columns:
for c in df.columns:
tmp = df[df[r].notnull() & df[c].notnull()]
pvalues[r][c] = round(scs.pearsonr(tmp[r], tmp[c])[1], 4)
return pvalues
final_data_T = alldata.transpose()
frame_patients = pd.DataFrame(data=final_data_T,index=
[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38],
columns=['Ferittin','Transferin','Iron','Transferin Sat','Sol Transf Rec','Putamen BPnd','Caudate
BPnd','Striatum BPnd','Pallidum BPnd','Thalamus BPnd','weights' ])
corr = frame_patients.corr(method = 'pearson')
pvals = calculate_pvalues(frame_patients)
plt.figure(figsize=(20,20), dpi =300)
sb.heatmap(corr,annot=True, square = 'True',cmap="vlag")
plt.figure(figsize=(20,20), dpi =300)
sb.heatmap(corr,annot=pvals, square = 'True',cmap="rocket_r",vmin=0, vmax=1)
Weights = alldata[10,:]
Weights_colour = []
for i in range(np.size(Weights)):
if (Weights[i] == 1):
Weights_colour.append("cornflowerblue")
else:
Weights_colour.append("maroon")
data = pd.DataFrame(data={'Putamen BPnd': final_data_T[:,5], 'FERRITIN': final_data_T[:,0]})
#,'Weights': Weights
data['color']= Weights_colour
plt.figure(figsize=(5,5), dpi =300)

sb.set_style("darkgrid")
sb.set_context("paper", font_scale=1.5, rc={"lines.linewidth": 2})
g = sb.regplot(x="FERRITIN", y="Putamen BPnd", data=data,scatter_kws={"s": 60,"color":
"black",'facecolors':data['color']})#,hue="Weights"
#g = (g.set(xlim=(-50, 480), ylim=(0.030, 0.09)))#-49, 270
data1 = pd.DataFrame(data={'x': final_data_T[:,5], 'y': final_data_T[:,0] }) #,'Weights': Weights
p = sb.regplot(x="y", y="x",data=data1)
#calculate slope and intercept of regression equation
slope, intercept, r, g, sterr = scipy.stats.linregress(x=p.get_lines()[0].get_xdata(),

y=p.get_lines()[0].get_ydata())

print(intercept, slope)
