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
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys;

from model_training import *;
#from myHelpers import *;


np.random.seed(1234);

#from tho.yawnn.MLutils import *;


np.set_printoptions(precision=3,suppress=True,linewidth=200);

#from skbio.stats.composition import clr
#from scipy.stats.mstats import gmean


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

    #try:

resultsDumpFileName='results/model_%s-skLR_LOO__%s_%s.pickle'%(binType,normType,expName);

x_all_base=x_data;
x_all=normalizeX(x_all_base,normType);
y_all=y_data[:,np.newaxis];
n_all=x_all.shape[0];
colNames=x_colNames;


res=runLOO_LR_LT(x_all,y_all,expName)
pickle.dump(res,open('results/%s_LOO.pkl'%expName,'wb'),protocol=4)
(oTeAll,oTeWAll,y_all01,wAll,bAll,paramsAll,fRealCntAll, rocCurr,statsTe)=res;

'''
modelingCohortFirstVisits_16S 25 90 [4.289 0.723 0.723 0.767 0.763 0.774 0.632 0.865] (0.723346090759978, 0.76666665, 24, 45, 14, 7, 0.7627118644067796, 0.7741935483870968, 0.631578947368421, 0.8653846153846154, 0.7666666666666667)
ROC: 0.723
Sensitivity: 0.774
Specificity: 0.763

'''