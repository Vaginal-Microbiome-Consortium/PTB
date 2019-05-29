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

np.set_printoptions(precision=3,suppress=True,linewidth=200);

def span01(x):
    x=x.copy();
    sCnt=x.shape[0];
    fCnt=x.shape[1];
    y=np.zeros((sCnt,fCnt));
    for i in range(0,fCnt):
        xx=x[:,i];
        xx=xx-np.min(xx);
        mx=np.max(xx);
        if (mx>0.0):
            xx=xx/mx;
        y[:,i]=xx;
    return y;


dataPrefix='ptb47';
binTypeIn='oneF24V'; #wide range of visit dates, all 16S modalities (MV1D, MRCD, MCKD, BRCD, BCKD)
binType='oneF24V'; #wide range of visit dates, all 16S modalities (MV1D, MRCD, MCKD, BRCD, BCKD)


(x_data,x_clr_data,x_colNames,y_data,t_data,gaSample_all,gaDelivery_all,cl_data,clinical_var_list)=pickle.load(open('data-sets-manuscript/firstVisits-modeling-oneF24V.pkl','rb'))



expName='modelingCohortFirstVisits_Clinical';
normType='zeroone';

    #try:

resultsDumpFileName='results/model_%s-skLR_LOO__%s_%s.pickle'%(binType,normType,expName);

x_all=span01(cl_data);
y_all=y_data[:,np.newaxis];
n_all=x_all.shape[0];
colNames=clinical_var_list;


res=runLOO_LR_LT(x_all,y_all,expName)
pickle.dump(res,open('results/%s_LOO.pkl'%expName,'wb'),protocol=4)
(oTeAll,oTeWAll,y_all01,wAll,bAll,paramsAll,fRealCntAll, rocCurr,statsTe)=res;

'''
modelingCohortFirstVisits_Clinical 11 90 [4.111 0.764 0.764 0.711 0.695 0.742 0.561 0.837] (0.7638053581191908, 0.7111111, 23, 41, 18, 8, 0.6949152542372882, 0.7419354838709677, 0.5609756097560976, 0.8367346938775511, 0.7111111111111111)
ROC: 0.764
Sensitivity: 0.742
Specificity: 0.695
'''