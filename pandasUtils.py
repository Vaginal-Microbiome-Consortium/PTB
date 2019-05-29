#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
yawnn - yet another wrapper for nnets

Created on Tue Jun 20 21:47:18 2017

@author: tarodz@vcu.edu
"""

import pandas as pd;
import numpy as np;

class PandasUtils:

    def normalizeRowsSum1(df):
        sdf=df.sum(axis=1);
        return df.div(sdf, axis=0),sdf;

    def normalizeColsMax1(df):
        mdf=df.max(axis=0).replace(0,1.0);
        return df.div(mdf,axis=1),mdf;

    def filterConstColumns(dfX,colN,mx=None):
        xMax=dfX.max(axis=0);
        xMin=dfX.min(axis=0);
        xDelta=xMax-xMin;
        xDeltaNZ=xDelta.nonzero()[0];
        xDeltaNZs=set(xDeltaNZ);
        allCols=list(dfX.columns)
        nzCols=[allCols[i] for i in xDeltaNZs];
        nzColN=[colN[i] for i in xDeltaNZs]
        dfXo=dfX[nzCols].copy();
        if (mx is not None):
            mxo=mx[np.ix_(xDeltaNZ,xDeltaNZ)]
            return dfXo,nzColN,mxo;
        else:
            return dfXo,nzColN;

    def df2idx(df):
        df.sort_index(inplace=True);
        c=list(df.columns);
        return df.drop(c,1);

    def innerNonJoin(dflist):
        odflist=[];
        dfAgg=PandasUtils.df2idx(dflist[0]);
        for df in dflist:
            dfAgg=dfAgg.join(PandasUtils.df2idx(df),how='inner');
        for df in dflist:
            dfO=df.join(dfAgg,how='inner');
            odflist.append(dfO);
        return odflist;

    def commonColsInner(dflist):
        cols=list(dflist[0].columns);
        for df in dflist:
            colsC=set(df.columns);
            cols=[x for x in cols if x in colsC];
        odflist=[];
        for df in dflist:
            dfO=df[cols];
            odflist.append(dfO);
        return odflist;

    def toNumpySetClassOnlyNonBinary(dfY,colName):
        selIndices=list(dfY.index);
        y_all=dfY[[colName]].as_matrix().astype(dtype='float32');
        n_all=y_all.shape[0]
        y_all_int=y_all.astype(dtype='int32').reshape((n_all,));
        return y_all,y_all_int,selIndices;

    def toNumpySetClassOnly(dfY,colName):
        dfYs=dfY.loc[dfY[colName].isin([0, 1])];
        selIndices=list(dfYs.index);
        y_all=np.sign(dfYs[[colName]].as_matrix().astype(dtype='float32'));
        n_all=y_all.shape[0]
        y_all_int=y_all.astype(dtype='int32').reshape((n_all,));
        return y_all,y_all_int,selIndices;

    def toNumpySetFromIndex(dfX,idx):
        dfXs=dfX.loc[idx];
        x_all=dfXs.as_matrix().astype(dtype='float32');
        return x_all;

    def toNumpySetInOrder(dfX):
        x_all=dfX.as_matrix().astype(dtype='float32');
        c_all=list(dfX.columns);
        return x_all,c_all;


    def toNumpySet(dfX,dfY,colName):
        dfXs=dfX.loc[dfY[colName].isin([0, 1])];
        dfYs=dfY.loc[dfY[colName].isin([0, 1])];
        x_all=dfXs.as_matrix().astype(dtype='float32');
        y_all=np.sign(dfYs[[colName]].as_matrix().astype(dtype='float32'));
        n_all=x_all.shape[0]
        y_all_int=y_all.astype(dtype='int32').reshape((n_all,));
        return x_all,y_all,n_all,y_all_int;

