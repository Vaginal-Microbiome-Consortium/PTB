#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
yawnn - yet another wrapper for nnets

Created on Tue Jun 20 21:47:18 2017

@author: tarodz@vcu.edu
"""

import sklearn.metrics;
import numpy as np;
import copy;
from sklearn.model_selection import StratifiedKFold;
import statsmodels.stats.multitest;
import scipy.stats;

class MLUtils:

    def normalize01(mx):
    # see sklearn.preprocessing.normalize for other normalizations
        sCnt=mx.shape[0];
        fCnt=mx.shape[1];
        for i in range(0,fCnt):
            minI=np.min(mx[:,i]);
            mx[:,i]=mx[:,i]-minI;
            maxI=np.max(mx[:,i]);
            if (maxI>0):
                mx[:,i]=mx[:,i]/maxI;
        return mx;

    def shapeWeights(mx):
        if (mx.ndim==1):
            fCnt=mx.shape[0];
            mx=np.reshape(mx,(1,fCnt))
        elif (mx.shape[1]==1):
            fCnt=mx.shape[0];
            mx=np.reshape(mx,(1,fCnt))
        rCnt=mx.shape[0];
        fCnt=mx.shape[1];
        return mx,rCnt,fCnt;

    def reduceWeights(mx,dx):
        mx,rCnt,fCnt=MLUtils.shapeWeights(mx);
        mxS=np.sign(mx);
        mxA=np.abs(mx);
        mxR=np.maximum(mxA-dx,0);
        mxRS=np.abs(np.sign(mxR));
        mxOrgSmallZeroed=mxRS*mx;
        mxReduced=mxR*mxS;
        amxReduced=np.sum(mxReduced,axis=0)/rCnt;
        amxOrgSmallZeroed=np.sum(mxOrgSmallZeroed,axis=0)/rCnt;
        return (amxOrgSmallZeroed,amxReduced,mxOrgSmallZeroed,mxReduced);

    def reduceWeightsByRange(mx,dxr):
        #use e.g. with dxr=np.linspace(2.0, 3.0, num=6) to get spacing by 0.2
        mx,rCnt,fCnt=MLUtils.shapeWeights(mx);
        dCnt=len(dxr);
        rMxOSZ=np.zeros((dCnt,fCnt));
        rMxR=np.zeros((dCnt,fCnt));
        for dxI in range(0,dCnt):
            dx=dxr[dxI];
            (amxOrgSmallZeroed,amxReduced,mxOrgSmallZeroed,mxReduced)=MLUtils.reduceWeights(mx,dx);
            rMxOSZ[dxI,:]=amxOrgSmallZeroed;
            rMxR[dxI,:]=amxReduced;
        return (rMxOSZ,rMxR);

    def cvcvGen(y_all_int,oN,iN):
        n_all=y_all_int.shape[0]
        x_all=y_all_int.astype(dtype='float32').reshape((n_all,1));
        cvcvTrains=[];
        cvcvTests=[];
        for cvi in range(0,oN):
            skf = StratifiedKFold(n_splits=iN,shuffle=True);
            cvTrain=[];
            cvTest=[];
            #for train_index, test_index in skf.split(X, y):
            for train_index, test_index in skf.split(x_all,y_all_int):
                #print("TRAIN:", train_index, "TEST:", test_index)
                cvTrain.append(train_index);
                cvTest.append(test_index);
            cvcvTrains.append(cvTrain);
            cvcvTests.append(cvTest);
        return cvcvTrains,cvcvTests;

    def net2Laplacian(mx):
        inDim=mx.shape[0];
        mx=0.5*(mx+np.transpose(mx));
        mx=-mx;
        for gi in range(0,inDim):
            ss=np.sum(mx[gi,:]);
            mx[gi,gi]=-ss;
        return mx;

    def makeGraphForRepeatedF(fCnt,rCnt):
        graph=np.zeros((fCnt*rCnt,fCnt*rCnt));
        for i in range(0,fCnt):
            for j in range(1,rCnt):
                graph[(j-1)*fCnt+i,j*fCnt+i]=1;
                graph[j*fCnt+i,(j-1)*fCnt+i]=1;
        return graph;

class MLMultiClassUtils:
    def toCategorical(mx,nClasses):
        y = np.array(mx, dtype='int').ravel()
        if not nClasses:
            nClasses = np.max(y) + 1
        n = y.shape[0]
        categorical = np.zeros((n, nClasses))
        categorical[np.arange(n), y] = 1
        return categorical

    def accuracyCategorical(pred,true):
        n = true.shape[0]
        nClasses=true.shape[1];
        py=MLMultiClassUtils.toCategorical(np.argmax(pred,axis=1),nClasses)
        return np.sum(np.multiply(py,true))/n;

class MLBinaryClassUtils:

    def findSensSpecEqThresh01(pred,true):
        sCnt=pred.shape[0];
        pred=np.reshape(pred,(sCnt,1));
        if (np.min(pred)==np.max(pred)):
            return np.min(pred);
        true01=np.reshape(MLBinaryClassUtils.classNPto01(true),(sCnt,1));
        sCnt=pred.shape[0];

        pt=np.zeros((sCnt+2,8));

        pt[1:-1,0:3]=np.concatenate((pred,true01,1-true01),axis=1)
        pts=pt[pt[:,0].argsort()]
        pts[:,3]=np.cumsum(pts[:,1])/np.sum(true01);

        xx=pts[1:,2]
        xx=xx[::-1];

        zz=np.cumsum(xx)/np.sum(1-true01);
        zzz=zz[::-1]
        pts[0:-1,4]=zzz;
        pts[0:-2,5]=(pts[0:-2,0]+pts[1:-1,0])/2
        pts[0:-2,6]=np.sign(pts[1:-1,0]-pts[0:-2,0]);
        pts[:,7]=np.abs(pts[:,3]-pts[:,4])+1000*(1-pts[:,6])

        ami=np.argmin(pts[:,7]);
        thresh=pts[ami,5];
        return np.mean(thresh);



    def findSensSpecEqThresh01sens(pred,true):
        sCnt=pred.shape[0];
        pred=np.reshape(pred,(sCnt,1));
        if (np.min(pred)==np.max(pred)):
            return np.min(pred);
        true01=np.reshape(MLBinaryClassUtils.classNPto01(true),(sCnt,1));
        sCnt=pred.shape[0];

        pt=np.zeros((sCnt+2,8));

        pt[1:-1,0:3]=np.concatenate((pred,true01,1-true01),axis=1)
        pts=pt[pt[:,0].argsort()]
        pts[:,3]=np.cumsum(pts[:,1])/np.sum(true01);

        xx=pts[1:,2]
        xx=xx[::-1];

        zz=np.cumsum(xx)/np.sum(1-true01);
        zzz=zz[::-1]
        pts[0:-1,4]=zzz;
        pts[0:-2,5]=(pts[0:-2,0]+pts[1:-1,0])/2
        pts[0:-2,6]=np.sign(pts[1:-1,0]-pts[0:-2,0]);
        pts[:,7]=np.abs(pts[:,3]-pts[:,4])+1000*(1-pts[:,6])

        ami=np.argmin(pts[:,7]);
        thresh=pts[ami-3,5];
        return np.mean(thresh);


    def class01toNP(y):
        y=copy.copy(y);
        y[y<=0.0001]=-1;
        y[y>-0.9999]=1;
        return y.astype(dtype='float32');

    def classNPto01(y):
        y=copy.copy(y);
        y[y<=0.0001]=0;
        y[y>0.0001]=1;
        return y.astype(dtype='float32');

    def countClasses(y):
        y=MLBinaryClassUtils.classNPto01(y);
        pe=np.extract(y>0,y);
        p=np.shape(pe)[0];
        a=y.shape[0];
        n=a-p;
        return (a,p,n);

    def binaryAccuracyNP(pred,true):
        #threshold on pred is at 0
        pred=np.sign(pred.astype(dtype='float32'));
        true=MLBinaryClassUtils.class01toNP(true);
        z=pred*true;
        zz=np.sign(1.0+z);
        zzm=np.mean(zz);
        return zzm;

    def binaryConfusionMxNP(pred,true):
        #threshold on pred is at 0
        pred=np.sign(pred.astype(dtype='float32'));
        true=MLBinaryClassUtils.class01toNP(true);
        pred=MLBinaryClassUtils.classNPto01(pred);
        true=MLBinaryClassUtils.classNPto01(true);
        z=2*pred+true;
        tn=np.shape(np.extract(z==0,z))[0];
        fn=np.shape(np.extract(z==1,z))[0];
        fp=np.shape(np.extract(z==2,z))[0];
        tp=np.shape(np.extract(z==3,z))[0];
        if (tp>0):
            sens=tp/(tp+fn); #recall
            ppv=tp/(tp+fp); #precision
        else:
            sens=0;
            ppv=0;
        if (tn>0):
            spec=tn/(tn+fp);
            npv=tn/(tn+fn);
        else:
            spec=0;
            npv=0;
        acc=(tp+tn)/(tp+tn+fn+fp); 
        return (tp,tn,fp,fn,spec,sens,ppv,npv,acc);

    def binaryAUCNP(pred,true):
        true=MLBinaryClassUtils.classNPto01(true); 
        (a,p,n)=MLBinaryClassUtils.countClasses(true);
        if (p*n==0):
            print("WARNING from MLBinaryClassUtils.binaryAUCNP: empty class (a,p,n)="+str((a,p,n)));
        sklearn.metrics.roc_curve(true, pred);
        auc = sklearn.metrics.roc_auc_score(true,pred);
        return auc;

    def binaryAccuracy01(pred,true):
        #threshold on pred is at 0.5
        true=MLBinaryClassUtils.class01toNP(true);
        pred=np.sign(pred.astype(dtype='float32')-0.5);
        return MLBinaryClassUtils.binaryAccuracyNP(pred,true);

    def binaryConfusionMx01(pred,true):
        #threshold on pred is at 0.5
        true=MLBinaryClassUtils.class01toNP(true);
        pred=np.sign(pred.astype(dtype='float32')-0.5);
        return MLBinaryClassUtils.binaryConfusionMxNP(pred,true);

    def binaryAUC01(pred,true):
        true=MLBinaryClassUtils.class01toNP(true);
        return MLBinaryClassUtils.binaryAUCNP(pred,true);

    def binaryAccuracy(pred,true,isNP=None):
        if (isNP is None):
            isNP=MLBinaryClassUtils.detectIsNP(true);
        if (isNP):
            return MLBinaryClassUtils.binaryAccuracyNP(pred,true);
        else:
            return MLBinaryClassUtils.binaryAccuracy01(pred,true);

    def binaryConfusionMx(pred,true,isNP=None):
        if (isNP is None):
            isNP=MLBinaryClassUtils.detectIsNP(true);
        if (isNP):
            return MLBinaryClassUtils.binaryConfusionMxNP(pred,true);
        else:
            return MLBinaryClassUtils.binaryConfusionMx01(pred,true);

    def binaryAUC(pred,true,isNP=None):
        if (isNP is None):
            isNP=MLBinaryClassUtils.detectIsNP(true);
        if (isNP):
            return MLBinaryClassUtils.binaryAUCNP(pred,true);
        else:
            return MLBinaryClassUtils.binaryAUC01(pred,true);

    def detectIsNP(true):
        mi=np.min(true);
        if (mi<0):
            isNP=True;
        else:
            isNP=False;
        return isNP;


    def binaryClassifStats(pred,true,isNP=None):
        # print(pred)
        if (isNP is None):
            isNP=MLBinaryClassUtils.detectIsNP(true);
        if (isNP):
            acc=MLBinaryClassUtils.binaryAccuracyNP(pred,true);
            roc=MLBinaryClassUtils.binaryAUCNP(pred,true);
            (tp,tn,fp,fn,spec,sens,ppv,npv,accCMX)=MLBinaryClassUtils.binaryConfusionMxNP(pred,true);
        else:
            acc=MLBinaryClassUtils.binaryAccuracy01(pred,true);
            roc=MLBinaryClassUtils.binaryAUC01(pred,true);
            (tp,tn,fp,fn,spec,sens,ppv,npv,accCMX)=MLBinaryClassUtils.binaryConfusionMx01(pred,true);
        return (roc,acc,tp,tn,fp,fn,spec,sens,ppv,npv,accCMX);

    def allFeatureStats(x,true):
        trueNP=MLBinaryClassUtils.class01toNP(true);
        true01=MLBinaryClassUtils.classNPto01(trueNP);
        fCnt=x.shape[1];
        sCnt=x.shape[0];

        r = np.zeros((fCnt,14));
        pv=np.zeros((fCnt,));

        for f in range(0,fCnt):
            pred=x[:,f];
            if (np.min(pred)<np.max(pred)):
                cc,cp=scipy.stats.spearmanr(pred,trueNP);
                if (cc<0):
                    trueNP=-trueNP;
                    true01=1-true01;
                t=MLBinaryClassUtils.findSensSpecEqThresh01(pred,true01)
                r[f,0]=cc;
                r[f,1]=t;
                r[f,2]=MLBinaryClassUtils.binaryAUCNP(pred,trueNP);
                r[f,3]=MLBinaryClassUtils.binaryAccuracyNP(pred-t,trueNP);
                (tp,tn,fp,fn,spec,sens,ppv,npv,accCMX)=MLBinaryClassUtils.binaryConfusionMxNP(pred-t,trueNP);
                r[f,4]=sens;
                r[f,5]=spec;
                r[f,6]=ppv;
                r[f,7]=npv;
                try:
                    [ss,p]=scipy.stats.mannwhitneyu(pred[true01==1],pred[true01==0]);
                except ValueError:
                    p=1.0;
                r[f,8]=p;
                if (p<=0.05):
                    r[f,9]=1;
                pv[f]=p;
            else:
                r[f,8]=1.0;
                pv[f]=1.0;
        hAdj,pAdj, alphacSidak, alphacBonf=statsmodels.stats.multitest.multipletests(pv,alpha=0.05,method='fdr_bh')
        hBonf,pBonf, alphacSidak, alphacBonf=statsmodels.stats.multitest.multipletests(pv,alpha=0.05,method='b')
        r[:,10]=pAdj;
        r[:,11]=hAdj;
        r[:,12]=pBonf;
        r[:,13]=hBonf;
        return r;

    def classWeightsNP(y):
        (a,p,n)=MLBinaryClassUtils.countClasses(y);
        wp=a/(2*p); # p*wp = a/2
        wn=a/(2*n); # n*wn = a/2
        return (wp,wn);

    def classToSampleWeightsNP(y):
        (wp,wn)=MLBinaryClassUtils.classWeightsNP(y);
        a=y.shape[0];
        w=wn*np.ones((a,1));
        np.place(w,y>0,wp);
        return w.astype(dtype='float32');

    def classUnitWeights(y):
        a=y.shape[0];
        w=np.ones((a,1));
        return w.astype(dtype='float32');

    def weights2ROC(mx,x,y):
        #mx is a NN layer, rCnt by fCnt, we get ROC for each of the rCnt neurons
        mx,rCnt,fCnt=MLUtils.shapeWeights(mx);
        rocAll=np.zeros((rCnt,4));
        for i in range(0,rCnt):
            wCurr=mx[i,:];
            nzW=np.sign(wCurr);
            wNZ=np.sum(np.abs(nzW));
            wP=0.5*np.sum(1+np.sign( nzW-0.5 ));
            wN=0.5*np.sum(1+np.sign( -1*nzW-0.5 ));
            predCurr=np.dot(x,wCurr.transpose());
            (roc,acc,tp,tn,fp,fn,spec,sens,ppv,npv,accCMX)=MLBinaryClassUtils.binaryClassifStats(predCurr,y);
            rocAll[i,0]=roc;
            rocAll[i,1]=wNZ;
            rocAll[i,2]=wP;
            rocAll[i,3]=wN;
        return rocAll;

class MLRegressionUtils:

    def regressionStats(pred,true):
        avgExpVar=sklearn.metrics.explained_variance_score(true,pred,multioutput='uniform_average');
        allExpVar=sklearn.metrics.explained_variance_score(true,pred,multioutput='raw_values');

        avgMAE=sklearn.metrics.mean_absolute_error(true,pred,multioutput='uniform_average');
        allMAE=sklearn.metrics.mean_absolute_error(true,pred,multioutput='raw_values');

        avgMSE=sklearn.metrics.mean_squared_error(true,pred,multioutput='uniform_average');
        allMSE=sklearn.metrics.mean_squared_error(true,pred,multioutput='raw_values');

        avgR2=sklearn.metrics.r2_score(true,pred,multioutput='uniform_average');
        allR2=sklearn.metrics.r2_score(true,pred,multioutput='raw_values');

        return (avgR2,avgMSE,avgMAE,avgExpVar,allR2,allMSE,allMAE,allExpVar);
