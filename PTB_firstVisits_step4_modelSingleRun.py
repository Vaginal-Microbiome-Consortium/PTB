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
import sys;
from sklearn import preprocessing;
from scipy.stats.mstats import gmean

np.random.seed(1234);
np.set_printoptions(precision=3,suppress=True,linewidth=200);

from MLutils import *;
from model_training import *;

dataPrefix='ptb47';
binTypeIn='oneF24V'; #wide range of visit dates, all 16S modalities (MV1D, MRCD, MCKD, BRCD, BCKD)
binType='oneF24V'; #wide range of visit dates, all 16S modalities (MV1D, MRCD, MCKD, BRCD, BCKD)

def model_plot(w, n ):
    if (len(w)<1):
        return None;
    dictionary = dict(zip(w, n))
    w=sorted(w)
    params = {'back_color': 'w', 'autoscale': 1, 'orientation': 'h', 'format': 'png', 'all_feats': '', 'output_file': 'outi/jjj', 'rs': 0.1, 'class_legend_font_size': '10', 'height': 4.0, 'width': 7.0, 'feature_font_size': 7, 'title_font_size': '12', 'input_file': 'hmp_aerobiosis_small.in', 'title': '', 'max_feature_len': 60, 'n_scl': 1, 'dpi': 72, 'ls': 0.2}
    pos = np.arange(len(w))
    head = 0.75
    tail = 0.5
    ht = head + tail
    ints = max(len(pos)*0.2,1.5)
    #ints = 10


    fig = plt.figure(figsize=(params['width'], ints + ht), edgecolor=params['back_color'],facecolor=params['back_color'])
    ax = fig.add_subplot(111,frame_on=False)
    #ax.set_facecolor(params['back_color'])
    ls, rs = params['ls'], 1.0-params['rs']

    plt.subplots_adjust(left=ls,right=rs,top=1-head*(1.0-ints/(ints+ht)), bottom=tail*(1.0-ints/(ints+ht)))
    fig.canvas.set_window_title('results')


    l_align = {'horizontalalignment':'left', 'verticalalignment':'baseline'}
    r_align = {'horizontalalignment':'right', 'verticalalignment':'baseline'}
    added = []
    for i,v in enumerate(w):
        lab = n[i]
        if w[i] >0:
            col = 'g'
        else:
            col = 'r'
        ax.barh(pos[i],v, align='center', color=col, label=lab, height=0.8, edgecolor='k')

    mv = len(w)
    for i,r in enumerate(w):
        stringToPrint=dictionary[r].replace('_',' ');
        if (w[i]<0.2):
            ax.text(.11+mv/300.0, float(i)-0.3, stringToPrint, l_align, size=10,color='k')
        else:
            ax.text(.11+mv/300.0, float(i)-0.3, stringToPrint, l_align, size=10,color='w')

    ax.set_yticks([])
    ax.set_xlabel("WEIGHT")
    ax.xaxis.grid(True)
    xlim = ax.get_xlim()
    if params['autoscale']:
        ran = np.arange(0.0001,round(round((abs(xlim[0])+abs(xlim[1]))/10,4)*100,0)/100)
        if len(ran) > 1 and len(ran) < 100:
            ax.set_xticks(np.arange(xlim[0],xlim[1]+0.0001,min(xlim[1]+0.0001,round(round((abs(xlim[0])+abs(xlim[1]))/10,4)*100,0)/100)))
        ax.set_ylim((pos[0]-1,pos[-1]+1))
    #    leg = ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=5, mode="expand", borderaxespad=0., frameon=False,prop={'size':params['class_legend_font_size']})
#    plt.show()
    return fig;

def plotPs(expName,subName,x_all,y_all,colNames,roc,w,b,oTr,gaSample_all,gaDelivery_all):
    folder='paperFigs/';
    nzF=(w!=0);
    if (np.count_nonzero(nzF)>0):
        wSel=w[nzF];

        colSel = [i for indx,i in enumerate(colNames) if nzF[indx] == True]
        colSel= [i[5:50] for i in colSel];
        colSel2= [' ' for i in colSel];

        figL=model_plot(wSel,colSel);
        figL.savefig(folder+'Main-Fig2d.pdf',dpi=600);
        plt.close('all')

        #figL=model_plot(wSel,colSel2);
        #figL.savefig(folder+'Main-Fig2d_nonames.pdf',dpi=600);
        #plt.close('all')
        
        yVal=gaDelivery_all.flatten();
        yColName='gestational age at delivery [days]'
        yFileName='deliveryGA';
        
        figA, ax = plt.subplots()
        ax.scatter(oTr-b,yVal,s=16,c=y_all.flatten(),cmap='bwr',alpha=0.4,edgecolors='none')
        ax.set_ylabel(yColName);
        ax.set_xlabel('PTB predictive score')
        ax.set_xlim(xmin=-0.005,xmax=3.501)
        ax.set_ylim(ymin=140,ymax=300)
        figA.savefig(folder+'Extended-Fig3b.pdf',dpi=600);
        plt.close('all')

        yVal=gaDelivery_all.flatten();
        yColName=''
        yFileName='deliveryGAinset';
        figA, ax = plt.subplots()
        ax.scatter(oTr-b,yVal,s=16,c=y_all.flatten(),cmap='bwr',alpha=0.4,edgecolors='none')
        ax.set_ylabel(yColName);
        ax.set_xlabel('PTB predictive score')
        ax.set_xlim(xmin=-0.005,xmax=0.501)
        ax.set_ylim(ymin=270,ymax=290)
        figA.savefig(folder+'Extended-Fig3c.pdf',dpi=600);
        plt.close('all')


def normalizeX(x,nt):
    x=x.copy();
    if (nt=='log10-3'):
        x=x-0.001;
        x[x<0]=0;
        return np.log10((x+0.001)/0.001);
    
(x_data,x_clr_data,x_colNames,y_data,t_data,gaSample_all,gaDelivery_all,cl_data,clinical_var_list)=pickle.load(open('data-sets-manuscript/firstVisits-modeling-oneF24V.pkl','rb'))


expName='modelingCohortFirstVisits_16S';
normType='log10-3';


resultsDumpFileName='results/model_%s-skLR_SINGLE__%s_%s.pickle'%(binType,normType,expName);

x_all_base=x_data;
x_all=normalizeX(x_all_base,normType);
y_all=y_data[:,np.newaxis];
n_all=x_all.shape[0];
colNames=x_colNames;    #try:

res=runSH_LR_LT_ONCE(x_all,y_all,expName)
pickle.dump(res,open('results/%s_SINGLE.pkl'%expName,'wb'),protocol=4)
(roc,w,b,c,oTr)=res;
plotPs(expName,binType,x_all,y_all.flatten(),colNames,roc,w,b,oTr.flatten(),gaSample_all,gaDelivery_all);
print("Model:")
print("w="+str(w))
print("b=%f"%b);

'''
Model:
w=[0.    0.    0.751 0.    0.    0.011 0.    0.    0.    0.775 0.    0.    0.    0.    0.    0.116 0.    0.    0.    0.    0.    0.    0.    0.    0.   ]
b=-0.163088
'''