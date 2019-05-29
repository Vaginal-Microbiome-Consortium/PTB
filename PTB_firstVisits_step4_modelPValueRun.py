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
from sklearn import preprocessing;
from scipy.stats.mstats import gmean
import sys;

from MLutils import *;
from model_training import *;

np.random.seed(1234);
np.set_printoptions(precision=3,suppress=True,linewidth=200);

dataPrefix='ptb47';
binTypeIn='oneF24V'; #wide range of visit dates, all 16S modalities (MV1D, MRCD, MCKD, BRCD, BCKD)
binType='oneF24V'; #wide range of visit dates, all 16S modalities (MV1D, MRCD, MCKD, BRCD, BCKD)

def normalizeX(x,nt):
    x=x.copy();
    if (nt=='log10-3'):
        x=x-0.001;
        x[x<0]=0;
        return np.log10((x+0.001)/0.001);
    
(x_data,x_clr_data,x_colNames,y_data,t_data,gaSample_all,gaDelivery_all,cl_data,clinical_var_list)=pickle.load(open('data-sets-manuscript/firstVisits-modeling-oneF24V.pkl','rb'))


expName='modelingCohortFirstVisits_16S';
normType='log10-3';


resultsDumpFileName='results/model_%s-skLR_PVALS__%s_%s.pickle'%(binType,normType,expName);

x_all_base=x_data;
x_all=normalizeX(x_all_base,normType);
y_all=y_data[:,np.newaxis];
n_all=x_all.shape[0];
colNames=x_colNames;    #try:

res=runSH_LR_LT_SHUFFLE(x_all,y_all,expName,10000)
pickle.dump(res,open('results/%s_PVAL.pkl'%expName,'wb'),protocol=4)
(p,roc,w,b,c,oTr,allWeights,allBs)=res;
print("Model:")
print("w="+str(w))
print("b=%f"%b);
print("p=%f"%p);

'''
Model:
w=[0.    0.    0.751 0.    0.    0.011 0.    0.    0.    0.775 0.    0.    0.    0.    0.    0.116 0.    0.    0.    0.    0.    0.    0.    0.    0.   ]
b=-0.163088
p=0.002400
'''