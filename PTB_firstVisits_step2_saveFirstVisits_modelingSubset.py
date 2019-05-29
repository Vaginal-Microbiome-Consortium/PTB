#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 22:52:03 2017

@author: tarodz
"""
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt;
import seaborn as sns;

import pandas as pd;
import numpy as np;
import pickle
import pandas;
from pandas import DataFrame;
import scipy.stats;
import scipy;
from collections import defaultdict;
import statsmodels.stats.multitest;
import sys;
import os;

dataPrefix='ptb47';

np.random.seed(1234);
np.set_printoptions(precision=3,suppress=True,linewidth=200);

if (len(sys.argv)>1):
    binTypeIn=sys.argv[1];
else:
    binTypeIn='oneF24V'; #wide range of visit dates, all 16S modalities (MV1D, MRCD, MCKD, BRCD, BCKD)

binType=binTypeIn;


inputDumpFileName='data-sets-manuscript/'+dataPrefix+'_casecontrol_16StimeV8_%s_input.pickle'%binTypeIn;
thisSet=pickle.load( open( inputDumpFileName, "rb" ) );
(allData,cvcvTrains,cvcvTests,afam_all_indic,yAfAmClass_selIndices,y_all,y_all_int,n_all,all_days_org)=thisSet;


#df_X_Y_Z is a dataframe,
## all these dataframes should have the same rows - e.g. row 1 in all of them is the same subject (i.e., mother-baby pair), row 2 is another subject etc.
## columns describe the features, and differ
## name of the data frame has XYZ, they mean:
## X=1 or 2 (first or second sample collection range)
## Y=E or L (earliest or lates sample in the collection range)
## E or L should be in general different data; but often X=1 is the same as X=2 (for datasets starting with "one", for those with "two" they should be different)
## Z denotes what's actually there:
## T means "Tags" that is metadata including class, clinical, demographic info etc.
## R means some form of actual profiling data, 16S or WMTS or WMGS, from various bodysites, using various annotation
((df_1_E_T,df_2_E_T,df_1_L_T,df_2_L_T), # clinical info including class; all these dataframes should be the same except GA of sample collection
         (df_1_E_R,df_2_E_R,df_1_L_R,df_2_L_R), #mother vaginal; 'db99'
         (df_1_E_Rb,df_2_E_Rb,df_1_L_Rb,df_2_L_Rb), #mother rectal; 'db99'
         (df_1_E_Rc,df_2_E_Rc,df_1_L_Rc,df_2_L_Rc), #mother cheek; 'db99'
         (df_1_E_Rbb,df_2_E_Rbb,df_1_L_Rbb,df_2_L_Rbb), #baby rectal; 'db99'
         (df_1_E_Rbc,df_2_E_Rbc,df_1_L_Rbc,df_2_L_Rbc), #baby cheek; 'db99'
         (df_1_E_Ro,df_2_E_Ro,df_1_L_Ro,df_2_L_Ro), #mother vaingal; 'STIRRUPS v2'

         (df_1_E_RF,df_2_E_RF,df_1_L_RF,df_2_L_RF), #same in section above, except selected only top abundand MV1D species
         (df_1_E_RFb,df_2_E_RFb,df_1_L_RFb,df_2_L_RFb),#same in section above, except selected only top abundand species + top abundant in MV1D
         (df_1_E_RFc,df_2_E_RFc,df_1_L_RFc,df_2_L_RFc),
         (df_1_E_RFbb,df_2_E_RFbb,df_1_L_RFbb,df_2_L_RFbb),
         (df_1_E_RFbc,df_2_E_RFbc,df_1_L_RFbc,df_2_L_RFbc),
         (df_1_E_RFo,df_2_E_RFo,df_1_L_RFo,df_2_L_RFo),

         (df_1_E_Rt,df_2_E_Rt,df_1_L_Rt,df_2_L_Rt), #mother vaginal WMTS UniRef90 abundances
         (df_1_E_Rtp,df_2_E_Rtp,df_1_L_Rtp,df_2_L_Rtp), #mother vaginal WMTS HUMANN pathways
         (df_1_E_Rg,df_2_E_Rg,df_1_L_Rg,df_2_L_Rg), #mother vaginal WMGS UniRef90 abundances
         (df_1_E_Rgp,df_2_E_Rgp,df_1_L_Rgp,df_2_L_Rgp) #mother vaginal WMGS HUMAN pathways
         ) = allData;

#df_1_E_T.to_csv('modelingCohort-T.csv')
#df_1_E_R.to_csv('modelingCohort-R.csv')
#df_1_E_RF.to_csv('modelingCohort-RF.csv')



def clr(x,C=1.0,dj=0.001):
    x=x.copy();
    sCnt=x.shape[0];
    fCnt=x.shape[1];
    y=np.zeros((sCnt,fCnt));
    for i in range(0,sCnt):
        xx=x[i,:];
        xii=np.nonzero(xx);
        xi=xii[0];
        nzc=len(xi);
        x[i,:]=x[i,:]*(1.0-((nzc*dj)/C));
        zzz=dj*np.ones((fCnt,));
        zzz[xi]=0.0;
        x[i,:]=x[i,:]+zzz;
        lx=np.log(x[i,:]);
        delta=np.sum(lx)/fCnt;
        y[i,:]=lx-delta;
    return y;
    
def normalizeX(x,nt):
    x=x.copy();
    if (nt=='clr'):
        x[x<0.00001]=0;
        x=clr(x,1.0,0.00001);
        return x;    

x_for_clr=df_1_E_R.as_matrix().astype(dtype='float32'); #for each subject, earliest sample with first GA range
x_all_clr=normalizeX(x_for_clr,'clr');

sCnt=x_all_clr.shape[0];

clr_all_cols=list(df_1_E_R.columns);
clr_sel_cols=list(df_1_E_RF.columns);
fCnt=len(clr_sel_cols)

x_clr=np.zeros((sCnt,fCnt))
xi=0;
for ci in clr_sel_cols:
    idx=clr_all_cols.index(ci);
    x_clr[:,xi]=x_all_clr[:,idx];
    xi=xi+1;

#16S "~30 top abundant species" data
x_all_dict=dict();
x_all_dict['MV1D-1E']=df_1_E_RF.as_matrix().astype(dtype='float32'); #for each subject, earliest sample with first GA range
x_all_dict['MV1D-1L']=df_1_L_RF.as_matrix().astype(dtype='float32'); #for each subject, latest sample with first GA range

#column names (i.e. species)
c_all_dict=dict();
c_all_dict['MV1D-1E']=list(df_1_E_RF.columns);
c_all_dict['MV1D-1L']=list(df_1_L_RF.columns);

#gestational ages at sample collection
d_all_dict=dict();
d_all_dict['MV1D-1E']=df_1_E_T[['GA']].as_matrix().astype(dtype='float32');
d_all_dict['MV1D-1L']=df_1_L_T[['GA']].as_matrix().astype(dtype='float32');


#gestational age ranges for specific dataset
t_all_dict=dict();
t_all_dict['MV1D-1E']=(all_days_org[0],all_days_org[1]);
t_all_dict['MV1D-1L']=(all_days_org[0],all_days_org[1]);

x_data=x_all_dict['MV1D-1E'];
t_data=t_all_dict['MV1D-1E'];
x_colNames=c_all_dict['MV1D-1E'];
y_data=y_all.flatten();
x_clr_data=x_clr;

clinical_var_list=[
 'short_cervix',
 'cerclage',
 'vaginal_ph',
 'mdl_bmi',
 'mld_progesterone',
 'gravida',
 'parity',
 'gralesspar',
 'hhq_miscarriage_or_stillbirth_sign',
 'preterm_historyBinary',
 'mdl_antibiotics']
d_all_df=df_1_E_T.copy();
cl_data=d_all_df[clinical_var_list].as_matrix().astype(dtype='float32');

gaSample_all=df_1_E_T[['GA']].as_matrix().astype(dtype='float32');

gaDelivery_all=df_1_E_T[['GA_delivery']].as_matrix().astype(dtype='float32');

print((y_data.shape[0],np.count_nonzero(np.where(y_data==1))))
pickle.dump((x_data,x_clr_data,x_colNames,y_data,t_data,gaSample_all,gaDelivery_all,cl_data,clinical_var_list),open('data-sets-manuscript/firstVisits-modeling-oneF24V.pkl','wb'))
