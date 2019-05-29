#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 13:10:55 2017

@author: tarodz
"""

import numpy as np;
np.set_printoptions(precision=3,suppress=True,linewidth=2000);

import gc;


import pickle;

from MLutils import *;

import shutil

#######

from sklearn.linear_model import LogisticRegression
from sklearn import metrics
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV

def runSK_LR_L_mannwhitney_grid(x_all,x_test,y_all01,clWdict):
    if (x_test is not None):
        x_test_org=x_test.copy();
    else:
        x_test_org=x_test;

    (sCnt,fCnt)=x_all.shape;
    yp=y_all01>0;
    yn=y_all01<=0;


    pv=np.zeros((fCnt,));

    for fI in range(0,fCnt):
        try:
            _,pv[fI]=scipy.stats.mannwhitneyu(x_all[yp,fI],x_all[yn,fI],alternative='two-sided');
        except ValueError:
            pv[fI]=1.0;


    goodF=np.zeros((fCnt,));
    goodFnoadj=np.zeros((fCnt,));
    goodFadj=np.zeros((fCnt,));

    usePabs=0.05;
    for fI in range(0,fCnt):
        if (pv[fI]<=usePabs):
            goodF[fI]=1;

    for fI in range(0,fCnt):
        if (pv[fI]<=usePabs):
            goodFnoadj[fI]=1;


    if (np.count_nonzero(goodF==1)>0):
        x_all=x_all[:,goodF==1];
        if (x_test is not None):
            x_test=x_test[:,goodF==1];
    else:
        x_all=x_all*0.0;
        if (x_test is not None):
            x_test=x_test*0.0;

    cL=np.power(2.0,np.arange(-10, 0,1));

    if (np.isscalar(cL)==False):
        cLrange=np.asarray(cL);
        cLrangeSK=(fCnt/sCnt)*(1.0/cLrange);
        logregC = LogisticRegression(penalty='l1',class_weight=clWdict);
        param_grid = {'C':cLrangeSK};
        skCLmx=np.zeros((20,));
        for gi in range(0,20):
            ordr=np.random.permutation(sCnt);
            cmod=GridSearchCV(logregC,param_grid,scoring='roc_auc');
            cmod.fit(x_all[ordr,:],y_all01[ordr]);
            grid_cv_results =  cmod.best_params_
            skCLmx[gi]=grid_cv_results['C'];
        skCL=np.mean(skCLmx);
        cL=(fCnt/sCnt)*(1.0/skCL);
        cL=cL*2.0;
#        print("Found C:"+str(cL)+" "+str(icL)+" from:"+str(cLrange))
        cLrange=tuple(cLrange);
    else:
        cLrange=[cL];

    params=(cL,cLrange);
    if (cL==0): cL=0.00001;

    if (cL==0):
        #sklearn does not support no regularization, have to simulate it by "huuge" C (C here is like SVM, i.e. 1/C)
        logreg = LogisticRegression(penalty='l2',C=1e10,class_weight=clWdict);
    else:
        logreg = LogisticRegression(penalty='l1',C=(fCnt/sCnt)*(1.0/cL),class_weight=clWdict);
    logreg.fit(x_all,y_all01);
    w2=np.squeeze(logreg.coef_)
    b=np.squeeze(logreg.intercept_)

    if (np.count_nonzero(goodF==1)>0):
        w=np.zeros((fCnt,))
        w[goodF==1]=w2;
    else:
        w=w2;

    oTr=np.squeeze(logreg.decision_function(x_all));

    eqthresh=MLBinaryClassUtils.findSensSpecEqThresh01sens(oTr,y_all01);

    orgB=b;
    ssB=b-eqthresh;

    # make sensitivity and specificity more even on the training set, but retain some of the original B
    newB=0.2*orgB+0.8*ssB;

    deltaB=newB-orgB;

    b=b+deltaB;
    oTr=oTr+deltaB;
    if (x_test is not None):
        oTe=np.squeeze(logreg.decision_function(x_test));
        oTe=oTe+deltaB;
        oTeW=np.matmul(x_test_org,w);
    else:
        oTe=None;
        oTeW=None;


    return w,b,cL,oTr,oTe,oTeW,params;


def runLOO_LR_LT(x_all,y_all,expName):

    x_all=x_all.astype(dtype='float32');

    fCnt=x_all.shape[1];
    sCnt=x_all.shape[0];

    y_all=np.squeeze(y_all)
    y_allNP=MLBinaryClassUtils.class01toNP(y_all);
    y_all01=MLBinaryClassUtils.classNPto01(y_allNP);

    oTeAll=np.zeros((sCnt,))
    oTeWAll=np.zeros((sCnt,))
    wAll=np.zeros((sCnt,fCnt))
    bAll=np.zeros((sCnt,))
    paramsAll=list(range(0,sCnt));
    fRealCntAll=np.zeros((sCnt,))

    for ti in range(0,sCnt):
        xTr=np.delete(x_all,ti,axis=0);
        xTe=x_all[[ti],:];
        yTr01=np.delete(y_all01,ti,axis=0);
        yTe01=y_all01[ti];
        yTrNP=np.delete(y_allNP,ti,axis=0);
        yTeNP=y_allNP[ti];

        (wp,wn)=MLBinaryClassUtils.classWeightsNP(yTrNP);

        clWdict=dict();
        clWdict[0]=wn;
        clWdict[1]=wp;
        w,b,c,oTr,oTe,oTeW,params=runSK_LR_L_mannwhitney_grid(xTr,xTe,yTr01,clWdict);
        wAll[ti,:]=w;
        bAll[ti]=b;
        oTeAll[ti]=oTe;
        oTeWAll[ti]=oTeW;
        paramsAll[ti]=params;
        fCntReal=np.count_nonzero(w);
        fRealCntAll[ti]=fCntReal;
    rocCurr=MLBinaryClassUtils.binaryAUC01(oTeAll,y_all01);
    statsTe=MLBinaryClassUtils.binaryClassifStats(oTeAll,y_allNP,True);
    (t_roc,t_acc,t_tp,t_tn,t_fp,t_fn,t_spec,t_sens,t_ppv,t_npv,t_accCMX)=statsTe;

    print(expName+" %d %d "%(fCnt,sCnt)+
          str(np.array([np.mean(fRealCntAll), rocCurr, t_roc,t_acc,t_spec,t_sens,t_ppv,t_npv ]))
          +" "+str(statsTe),flush=True)
    return (oTeAll,oTeWAll,y_all01,wAll,bAll,paramsAll,fRealCntAll, rocCurr,statsTe);

#####################################

def runSH_LR_LT_ONCE(x_all,y_all,expName):

    x_all=x_all.astype(dtype='float32');

    fCnt=x_all.shape[1];
    sCnt=x_all.shape[0];

    y_all=np.squeeze(y_all)
    y_allNP=MLBinaryClassUtils.class01toNP(y_all);
    y_all01=MLBinaryClassUtils.classNPto01(y_allNP);
    (wp,wn)=MLBinaryClassUtils.classWeightsNP(y_allNP);

    clWdict=dict();
    clWdict[0]=wn;
    clWdict[1]=wp;

    w,b,c,oTr,oTe,oTeW,params=runSK_LR_L_mannwhitney_grid(x_all.copy(),None,y_all01.copy(),clWdict);

    roc=MLBinaryClassUtils.binaryAUC01(oTr,y_all01);

    pvalue=0.0;
    isGood=0;
    if (pvalue<=0.05):
        isGood=1;
    fCntReal=np.count_nonzero(w);
    print(expName+" %d %d %d %f "%(fCnt,fCntReal,sCnt,c)+str(np.array([np.count_nonzero(w),c,roc ])),flush=True)
    return (roc,w,b,c,oTr);

def runSH_LR_LT_SHUFFLE(x_all,y_all,expName,shuffleCnt=10000):

    x_all=x_all.astype(dtype='float32');

    fCnt=x_all.shape[1];
    sCnt=x_all.shape[0];

    y_all=np.squeeze(y_all)
    y_allNP=MLBinaryClassUtils.class01toNP(y_all);
    y_all01=MLBinaryClassUtils.classNPto01(y_allNP);
    (wp,wn)=MLBinaryClassUtils.classWeightsNP(y_allNP);

    clWdict=dict();
    clWdict[0]=wn;
    clWdict[1]=wp;

    w,b,c,oTr,oTe,oTeW,params=runSK_LR_L_mannwhitney_grid(x_all.copy(),None,y_all01.copy(),clWdict);

    roc=MLBinaryClassUtils.binaryAUC01(oTr,y_all01);


    rocShAll=np.zeros((shuffleCnt,));
    allWeights=np.zeros((shuffleCnt, fCnt));
    allBs=np.zeros((shuffleCnt, ));
    allCs=np.zeros((shuffleCnt, ));

    for shI in range(0,shuffleCnt):
        np.random.shuffle(y_all01);
        wCurr,bCurr,cCurr,oTrCurr,oTeCurr,oTeWCurr,params=runSK_LR_L_mannwhitney_grid(x_all.copy(),None,y_all01.copy(),clWdict);
        rocCurr=MLBinaryClassUtils.binaryAUC01(oTrCurr,y_all01);
        rocShAll[shI]=rocCurr;
        wRow=np.transpose(wCurr);
        allWeights[shI,:]=wRow;
        allBs[shI]=bCurr;
        allCs[shI]=cCurr;
    pvalue=np.count_nonzero(rocShAll>=roc) / shuffleCnt;
    isGood=0;
    if (pvalue<=0.05):
        isGood=1;

    fCntReal=np.count_nonzero(w);
    print(expName+" %d %d %d %f "%(fCnt,fCntReal,sCnt,c)+str(np.array([np.count_nonzero(w),c,roc ])),flush=True)
    return (pvalue,roc,w,b,c,oTr,allWeights,allBs);
