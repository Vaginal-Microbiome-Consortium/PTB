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
import pickle;
import scipy;
import numpy as np;
import pandas;
from pandas import DataFrame;
from collections import defaultdict;
import sys;
import statsmodels.stats.multitest;
from MLutils import *;


dataPrefix='ptb47';


np.random.seed(1234);
np.set_printoptions(precision=3,suppress=True,linewidth=200);


#Open file for reading to remap names

f = open ('Abbrev-Species.csv','r')

#Create dictionary for abbreviations

abbrev = {}

for line in iter(f):

    ln = line.split(",")

    abbrev[ln[0].strip()] = ln[1].strip()

f.close()

def boxPlotAll(x_all,isSignif,y,colName,pAdj,pRaw,log10='log10-3'):
    
    letters='abcdefghijklmnop';
    
    yp=y>0;
    yn=y<=0;

    fCnt=isSignif.shape[0];
    
    sigCnt=int(np.sum(isSignif));
    
    #!!!!!!!!!!!
    sigCnt=fCnt;

    f, ax = plt.subplots(nrows=1,ncols=sigCnt, sharex=False, sharey=False)
    
    #color of median
    medianprops = dict(color='white') 
    whiskerprops=dict(color='black',linestyle='-')
    boxprops=dict(color='black')
    
    sI=0;
    
    bplots=[];
    
    print("Blue: TB, Red: PTB")
    print("On the plot, from left to right:")
    for cI in range(0,fCnt):
        if (isSignif[cI]>=-10):
            
            print(colNames[cI][5:]+" FDR-adjusted p-value:"+str(pAdj[cI])+" raw p-value:"+str(pAdj[cI]))
        
            x=x_all[:,cI];

            xp=x[yp];
            xn=x[yn];

            nCnt=xn.shape[0];
            pCnt=xp.shape[0];
            
            axC=ax[sI];
            
            medianprops = dict(color='white') 
            whiskerprops=dict(color='black',linestyle='-')
            boxprops=dict(color='black');
            
            yl=xn;
            yr=xp;
            bplt=axC.boxplot((yl,yr), whis='range',labels=['-','+'], sym='',manage_xticks=True,positions=[-0.1,0.6],widths=[0.35,0.35],medianprops=medianprops,notch=False,patch_artist=True,whiskerprops=whiskerprops,boxprops=boxprops)  #whis=[25, 75]


            bplt['boxes'][0].set_facecolor('#4393c3')
            bplt['boxes'][1].set_facecolor('#d6604d')
            if (sI>0):
                axC.spines['left'].set_visible(False)
            axC.spines['right'].set_visible(False)
            axC.spines['top'].set_visible(False)
            axC.spines['bottom'].set_visible(False)
            axC.tick_params(labeltop='off', labelright='off',top='off', right='off')
            if (sI==0):
                if (log10=='log10-5'):
                    axC.set_ylim(ymin=-0.1,ymax=5.1)
                    axC.set_yticks([0,1,2,3,4,5])
                    axC.set_yticklabels(['0',"$10^{-4}$","$10^{-3}$", "$10^{-2}$", "$10^{-1}$", "1"])
                if (log10=='log10-3'):
                    axC.set_ylim(ymin=-0.1,ymax=3.1)
                    axC.set_yticks([0,1,2,3])
                    axC.set_yticklabels(['0', "$10^{-2}$", "$10^{-1}$", "1"])
            else:
                if (log10=='log10-5'):
                    axC.set_ylim(ymin=-0.1,ymax=5.1)
                    axC.set_yticks(())
                if (log10=='log10-3'):
                    axC.set_ylim(ymin=-0.1,ymax=3.1)
                    axC.set_yticks(())
        
              
#            axC.set_xlabel(letters[sI])
            axC.set_xlabel(abbrev[colName[cI][5:]], rotation=45) 
            sI=sI+1;

    ax[0].set_ylabel('vaginal 16S rRNA abundance')
    
    plt.show();
        
    return f;




def getMWUPs(x_all,y,log10=False):
    fCnt=x_all.shape[1];

    yp=y>0;
    yn=y<=0;

    xp=x_all[yp];
    xn=x_all[yn];

    p=np.ones((fCnt,))
    fStats=np.zeros((fCnt,100));

    rocs=np.zeros((fCnt,));

    for i in range(0,fCnt):
        try:
            j=0;
            xpf=xp[:,i].flatten();
            xnf=xn[:,i].flatten();
            xaf=x_all[:,i].flatten();

            _,p[i]=scipy.stats.mannwhitneyu(xpf,xnf,alternative='two-sided');
            
            rocs[i]=MLBinaryClassUtils.binaryAUCNP(x_all[:,i].flatten(),y.flatten());

            fStats[i,j+0]=np.mean(xaf);
            fStats[i,j+1]=np.mean(xpf);
            fStats[i,j+2]=np.mean(xnf);
            fStats[i,j+3]=fStats[i,j+1]-fStats[i,j+2];
            j=j+4;
            fStats[i,j+0]=np.std(xaf);
            fStats[i,j+1]=np.std(xpf);
            fStats[i,j+2]=np.std(xnf);
            fStats[i,j+3]=fStats[i,j+1]-fStats[i,j+2];
            j=j+4;
            fStats[i,j+0]=np.median(xaf);
            fStats[i,j+1]=np.median(xpf);
            fStats[i,j+2]=np.median(xnf);
            fStats[i,j+3]=fStats[i,j+1]-fStats[i,j+2];
            j=j+4;
            fStats[i,j+0]=np.percentile(xaf,75);
            fStats[i,j+1]=np.percentile(xpf,75);
            fStats[i,j+2]=np.percentile(xnf,75);
            fStats[i,j+3]=fStats[i,j+1]-fStats[i,j+2];
            j=j+4;
            fStats[i,j+0]=np.percentile(xaf,25);
            fStats[i,j+1]=np.percentile(xpf,25);
            fStats[i,j+2]=np.percentile(xnf,25);
            fStats[i,j+3]=fStats[i,j+1]-fStats[i,j+2];
            j=j+4;
            fStats[i,j+0]=np.percentile(xaf,75)-np.percentile(xaf,25);
            fStats[i,j+1]=np.percentile(xpf,75)-np.percentile(xpf,25);
            fStats[i,j+2]=np.percentile(xnf,75)-np.percentile(xnf,25);
            fStats[i,j+3]=fStats[i,j+1]-fStats[i,j+2];
            j=j+4;
        except ValueError:
            pass

    hAdj,pAdj, alphacSidak, alphacBonf=statsmodels.stats.multitest.multipletests(p,alpha=0.05,method='fdr_bh')
    return p,pAdj,fStats,rocs;
cstr='config';
cstr=cstr+';'+'featID';
cstr=cstr+';'+'p-adj';
cstr=cstr+';'+'p-raw';

cstr=cstr+';'+'mean_all';
cstr=cstr+';'+'mean_p';
cstr=cstr+';'+'mean_n';
cstr=cstr+';'+'mean_p-n';

cstr=cstr+';'+'std_all';
cstr=cstr+';'+'std_p';
cstr=cstr+';'+'std_n';
cstr=cstr+';'+'std_p-n';

cstr=cstr+';'+'median_all';
cstr=cstr+';'+'median_p';
cstr=cstr+';'+'median_n';
cstr=cstr+';'+'median_p-n';

cstr=cstr+';'+'pct75_all';
cstr=cstr+';'+'pct75_p';
cstr=cstr+';'+'pct75_n';
cstr=cstr+';'+'pct75_p-n';

cstr=cstr+';'+'pct25_all';
cstr=cstr+';'+'pct25_p';
cstr=cstr+';'+'pct25_n';
cstr=cstr+';'+'pct25_p-n';

cstr=cstr+';'+'IQR_all';
cstr=cstr+';'+'IQR_p';
cstr=cstr+';'+'IQR_n';
cstr=cstr+';'+'IQR_p-n';

cstr=cstr+';'+'taxa';
print(cstr);





nonFDRp=0.05;



def normalizeX(x,nt):
    x=x.copy();
    if (nt=='raw-3'):
        x[x<0.001]=0;
        return x;
    if (nt=='log10-3'):
        x=x-0.001;
        x[x<0]=0;
        return np.log10((x+0.001)/0.001);
    if (nt=='log10-5'):
        return np.log10((x+0.00001)/0.00001);

#        try:

            #resultsDumpFileName='data-sets-results/ptb47_casecontrol_16StimeV7_%s-skLR_PVAL__%s_%s.pickle'%(binType,normType,expName);
#pickle.dump((x_data,x_colNames,y_data),open('firstVisits.pkl','wb'))
if (len(sys.argv)>1):
    binTypeIn=sys.argv[1];
else:
    binTypeIn='oneF24V'; #wide range of visit dates, all 16S modalities (MV1D, MRCD, MCKD, BRCD, BCKD)

binType=binTypeIn;


(x_data,x_clr_data,x_colNames,y_data,t_data,df_data)=pickle.load(open('data-sets-manuscript/firstVisits-full-oneALLV.pkl','rb'))


df_xx=df_data[[ 'MV1D-Sneathia_amnii',
 'MV1D-Lachnospiraceae_BVAB1',
 'MV1D-TM7_OTU-H1',
 'MV1D-Prevotella_cluster2',
 'MV1D-Aerococcus_christensenii',
 'MV1D-Clostridiales_BVAB2',
 'MV1D-Coriobacteriaceae_OTU27',
 'MV1D-Dialister_cluster51',
 'MV1D-Dialister_micraerophilus',
 'MV1D-Megasphaera_OTU70_type1',
 'MV1D-Parvimonas_OTU142',
 'MV1D-Prevotella_amnii',
 'MV1D-Sneathia_sanguinegens',
 'MV1D-Lactobacillus_crispatus_cluster',
 'MV1D-Lactobacillus_iners',
 'MV1D-Gardnerella_vaginalis' ]];
       
        

x_all_base=df_xx.as_matrix().astype(dtype='float32');

x_all=normalizeX(x_all_base,'raw-3');
x_all2=normalizeX(x_all_base,'log10-3');

colNames=[ 'MV1D-Sneathia_amnii',
 'MV1D-Lachnospiraceae_BVAB1',
 'MV1D-TM7_OTU-H1',
 'MV1D-Prevotella_cluster2',
 'MV1D-Aerococcus_christensenii',
 'MV1D-Clostridiales_BVAB2',
 'MV1D-Coriobacteriaceae_OTU27',
 'MV1D-Dialister_cluster51',
 'MV1D-Dialister_micraerophilus',
 'MV1D-Megasphaera_OTU70_type1',
 'MV1D-Parvimonas_OTU142',
 'MV1D-Prevotella_amnii',
 'MV1D-Sneathia_sanguinegens',
 'MV1D-Lactobacillus_crispatus_cluster',
 'MV1D-Lactobacillus_iners',
 'MV1D-Gardnerella_vaginalis' ];

pRaw,pAdj,fStats,rocs=getMWUPs(x_all,y_data,log10=False);
fCnt=x_all.shape[1];

cntSignif=0;

isSignif=np.zeros((fCnt,))

for cI in range(0,fCnt):
    cn=colNames[cI][5:];
    cn=cn.replace('/','-')
    cstr='VCU';
    cstr=cstr+';'+str(cI);
    cstr=cstr+';'+str(pAdj[cI]);
    cstr=cstr+';ROC='+str(rocs[cI]);
    if (pAdj[cI]<=0.05):
        cntSignif=cntSignif+1;
        isSignif[cI]=1;
    cstr=cstr+';'+str(pRaw[cI]);
    for sfI in range(0,24):
        cstr=cstr+';'+str(fStats[cI,sfI]);
    cstr=cstr+';'+colNames[cI][5:];
    print(cstr)
    
print("SIGNIF:"+str(cntSignif))

f=boxPlotAll(x_all2,isSignif,y_data,colNames,pAdj,pRaw,log10='log10-3')
f.savefig('paperFigs/Extended-Fig10a.pdf',dpi=600,bbox_inches='tight')#,pad_inches=.05)

#x_all3=normalizeX(x_all_base,'log10-5');
#f=boxPlotAll(x_all3,isSignif,y_data,colNames,pAdj,pRaw,log10='log10-5')
#f.savefig('paperFigs/Extended-Fig10a-log5.pdf',dpi=300,bbox_inches='tight')#,pad_inches=.05)


