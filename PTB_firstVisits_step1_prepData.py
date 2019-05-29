#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 14:56:45 2017

@author: tarodz
"""




import pickle;
import scipy;
import numpy as np;
import pandas;
from pandas import DataFrame;
from collections import defaultdict;
import sys;

from pandasUtils import *;
from MLutils import *;

import warnings;
warnings.simplefilter(action='ignore',category=FutureWarning)

import time;

dataPrefix='ptb47';

np.random.seed(1234);

#db16S='99'
db16S='S2'

useAAC=False;
#useAAC=True;

maxUsedDay=300;
goOne=0;
goThree=0;
minNumVisit=2;
justTry=False;

if (len(sys.argv)>1):
    binTypeIn=sys.argv[1];
    if (len(sys.argv)>2):
        otherCond=sys.argv[2];
    else:
        otherCond='';
    justTry=False;
else:
#    binTypeIn='oneALLV';
    binTypeIn='oneF24V';
    otherCond='';
    justTry=False;


def yesno(val):
    if (val=='Yes'):
        return 1;
    if (val=='YES'):
        return 1;
    if (val=='Y'):
        return 1;
    if (val=='No'):
        return -1;
    if (val=='NO'):
        return -1;
    if (val=='N'):
        return -1;
    return 0;
   

JN = [list of IDs];


if (binTypeIn[0:3]=='one'):
    minNumVisit=1;
    goOne=1;


useMCKD=False;
useMRCD=False;
useBRCD=False;
useWMTS=False;

if (binTypeIn[-1]=='V'):
    pass
if (binTypeIn[-1]=='B'):
    useMRCD=True;
if (binTypeIn[-1]=='C'):
    useMRCD=True;
    useMCKD=True;
if (binTypeIn[-1]=='D'):
    useBRCD=True;
    useMRCD=True;
    useMCKD=True;
if (binTypeIn[-1]=='P'):
    useBRCD=True;
if (binTypeIn[-1]=='Q'):
    useBRCD=True;
    useMRCD=True;
if (binTypeIn[-1]=='W'):
    useWMTS=True;
if (binTypeIn[-1]=='X'):
    useWMTS=True;
    useMRCD=True;
if (binTypeIn[-1]=='Y'):
    useWMTS=True;
    useBRCD=True;
if (binTypeIn[-1]=='Z'):
    useWMTS=True;
    useMRCD=True;
    useMCKD=True;
    useBRCD=True;

print("STARTING:"+binTypeIn+" "+otherCond);

binType=binTypeIn;

kidnum2pids=pickle.load( open( 'data/'+dataPrefix+'_ID2PID_fromSurvey.pickle', "rb" ) )


dfHHQ=pickle.load( open( 'data/'+dataPrefix+'_survey_HHQ.pickle', "rb" ) )
dfDND=pickle.load( open( 'data/'+dataPrefix+'_survey_DND.pickle', "rb" ) )
dfROSwAB=pickle.load( open( 'data/'+dataPrefix+'_survey_ROSwAB.pickle', "rb" ) )
dfROSHHQ=dfROSwAB.join(dfHHQ,how='left',on='PID')
dfSurvey=dfROSHHQ.join(dfDND,how='left',on='PID')

kitsDates4PIDS=pickle.load( open( 'data/'+dataPrefix+'_kitsDates4PID_fromSurvey.pickle', "rb" ) )
kitsDates4PIDC=pickle.load( open( 'data/'+dataPrefix+'_kitsDates4PID_fromClinical.pickle', "rb" ) );
nextTimeC=pickle.load( open( 'data/'+dataPrefix+'_clinical_nextTimePoint.pickle', "rb" ) )
nextTimeS=pickle.load( open( 'data/'+dataPrefix+'_survey_nextTimePoint.pickle', "rb" ) )

dfMRAbs=pickle.load( open( 'data/'+dataPrefix+'_MRAbs.pickle', "rb" ) );
dfTags=pickle.load( open( 'data/'+dataPrefix+'_tags.pickle', "rb" ) );

dfHHQ=pickle.load( open( 'data/'+dataPrefix+'_survey_HHQ.pickle', "rb" ) )
dfDND=pickle.load( open( 'data/'+dataPrefix+'_survey_DND.pickle', "rb" ) )
dfROSwAB=pickle.load( open( 'data/'+dataPrefix+'_survey_ROSwAB.pickle', "rb" ) )

dfHumAnn_wmts_pwy=pickle.load( open( 'data/'+dataPrefix+'_wmts_humann_pwy.pickle', "rb" ) ).fillna(0);
dfHumAnn_wmgs_pwy=pickle.load( open( 'data/'+dataPrefix+'_wmgs_humann_pwy.pickle', "rb" ) ).fillna(0);


dfBMI=pickle.load( open( 'data/'+dataPrefix+'_BMIetc.pickle', "rb" ) )



dfHumAnn_wmts_UR90l=pickle.load( open( 'data/'+dataPrefix+'_wmts_humann_UR90_log.pickle', "rb" ) )#.fillna(0);
dfHumAnn_wmgs_UR90l=pickle.load( open( 'data/'+dataPrefix+'_wmgs_humann_UR90_log.pickle', "rb" ) )#.fillna(0);

dfROSHHQ=dfROSwAB.join(dfHHQ,how='left',on='PID')
dfSurvey=dfROSHHQ.join(dfDND,how='left',on='PID')

dfClinical=pickle.load( open( 'data/'+dataPrefix+'_clinical.pickle', "rb" ) );

dfHHQsel=dfHHQ[['hhq_miscarriage_or_stillbirth','hhq_miscarriage_or_stillbirth_number','hhq_abortion_number','hhq_abortion']];

if (db16S=='99'):
    df16SpMV1D=pickle.load( open( 'data/'+dataPrefix+'_16S_MV1D_db99pct.pickle', "rb" ) ).applymap(float).fillna(0);
    dfXXX=df16SpMV1D;
    df16SpMV1D=dfXXX.div(dfXXX.sum(axis=1), axis=0).copy();
    colOrder=list(df16SpMV1D.mean(axis=0).sort_values(ascending=False).index);
    df16sort=df16SpMV1D[colOrder].copy();

    df16SpMV1Dalt=pickle.load( open( 'data/'+dataPrefix+'_16S_MV1D_STIRRUPSpctV2.pickle', "rb" ) ).applymap(float).fillna(0);
    dfXXX=df16SpMV1Dalt;
    df16SpMV1Dalt=dfXXX.div(dfXXX.sum(axis=1), axis=0).copy();
    colOrder=list(df16SpMV1Dalt.mean(axis=0).sort_values(ascending=False).index);
    df16SpMV1Dalt=df16SpMV1Dalt[colOrder].copy();

    df16SpMRCD=pickle.load( open( 'data/'+dataPrefix+'_16S_MRCD_db99pct.pickle', "rb" ) ).applymap(float).fillna(0);
    dfXXX=df16SpMRCD;
    df16SpMRCD=dfXXX.div(dfXXX.sum(axis=1), axis=0).copy();
    colOrder=list(df16SpMRCD.mean(axis=0).sort_values(ascending=False).index);
    df16SpMRCD=df16SpMRCD[colOrder].copy();

    df16SpMCKD=pickle.load( open( 'data/'+dataPrefix+'_16S_MCKD_db99pct.pickle', "rb" ) ).applymap(float).fillna(0);
    dfXXX=df16SpMCKD;
    df16SpMCKD=dfXXX.div(dfXXX.sum(axis=1), axis=0).copy();
    colOrder=list(df16SpMCKD.mean(axis=0).sort_values(ascending=False).index);
    df16SpMCKD=df16SpMCKD[colOrder].copy();

    df16SpBRCD=pickle.load( open( 'data/'+dataPrefix+'_16S_BRCD_db99pct.pickle', "rb" ) ).applymap(float).fillna(0);
    dfXXX=df16SpBRCD;
    df16SpBRCD=dfXXX.div(dfXXX.sum(axis=1), axis=0).copy();
    colOrder=list(df16SpBRCD.mean(axis=0).sort_values(ascending=False).index);
    df16SpBRCD=df16SpBRCD[colOrder].copy();

    df16SpBCKD=pickle.load( open( 'data/'+dataPrefix+'_16S_BCKD_db99pct.pickle', "rb" ) ).applymap(float).fillna(0);
    dfXXX=df16SpBCKD;
    df16SpBCKD=dfXXX.div(dfXXX.sum(axis=1), axis=0).copy();
    colOrder=list(df16SpBCKD.mean(axis=0).sort_values(ascending=False).index);
    df16SpBCKD=df16SpBCKD[colOrder].copy();

elif (db16S=='S2'):
    df16SpMV1D=pickle.load( open( 'data/'+dataPrefix+'_16S_MV1D_STIRRUPSpctV2.pickle', "rb" ) ).applymap(float).fillna(0);
    dfXXX=df16SpMV1D;
    df16SpMV1D=dfXXX.div(dfXXX.sum(axis=1), axis=0).copy();
    colOrder=list(df16SpMV1D.mean(axis=0).sort_values(ascending=False).index);
    df16sort=df16SpMV1D[colOrder].copy();

    df16SpMV1Dalt=pickle.load( open( 'data/'+dataPrefix+'_16S_MV1D_db99pct.pickle', "rb" ) ).applymap(float).fillna(0);
    dfXXX=df16SpMV1Dalt;
    df16SpMV1Dalt=dfXXX.div(dfXXX.sum(axis=1), axis=0).copy();
    colOrder=list(df16SpMV1Dalt.mean(axis=0).sort_values(ascending=False).index);
    df16SpMV1Dalt=df16SpMV1Dalt[colOrder].copy();

    df16SpMRCD=pickle.load( open( 'data/'+dataPrefix+'_16S_MRCD_db99pct.pickle', "rb" ) ).applymap(float).fillna(0);
    dfXXX=df16SpMRCD;
    df16SpMRCD=dfXXX.div(dfXXX.sum(axis=1), axis=0).copy();
    colOrder=list(df16SpMRCD.mean(axis=0).sort_values(ascending=False).index);
    df16SpMRCD=df16SpMRCD[colOrder].copy();

    df16SpMCKD=pickle.load( open( 'data/'+dataPrefix+'_16S_MCKD_db99pct.pickle', "rb" ) ).applymap(float).fillna(0);
    dfXXX=df16SpMCKD;
    df16SpMCKD=dfXXX.div(dfXXX.sum(axis=1), axis=0).copy();
    colOrder=list(df16SpMCKD.mean(axis=0).sort_values(ascending=False).index);
    df16SpMCKD=df16SpMCKD[colOrder].copy();

    df16SpBRCD=pickle.load( open( 'data/'+dataPrefix+'_16S_BRCD_db99pct.pickle', "rb" ) ).applymap(float).fillna(0);
    dfXXX=df16SpBRCD;
    df16SpBRCD=dfXXX.div(dfXXX.sum(axis=1), axis=0).copy();
    colOrder=list(df16SpBRCD.mean(axis=0).sort_values(ascending=False).index);
    df16SpBRCD=df16SpBRCD[colOrder].copy();

    df16SpBCKD=pickle.load( open( 'data/'+dataPrefix+'_16S_BCKD_db99pct.pickle', "rb" ) ).applymap(float).fillna(0);
    dfXXX=df16SpBCKD;
    df16SpBCKD=dfXXX.div(dfXXX.sum(axis=1), axis=0).copy();
    colOrder=list(df16SpBCKD.mean(axis=0).sort_values(ascending=False).index);
    df16SpBCKD=df16SpBCKD[colOrder].copy();

dfc=dfClinical[['days_rel2birth']].astype(dtype='float32').copy();

bc=df16SpBCKD.join(dfTags[['tags_PID']],how='inner');
bcc=bc.join(dfc,how='inner');
bcc0=bcc.loc[bcc['days_rel2birth']<=0].copy().set_index('tags_PID').copy();
bcc1=bcc.loc[bcc['days_rel2birth']>=1].copy().set_index('tags_PID').copy();

br=df16SpBRCD.join(dfTags[['tags_PID']],how='inner');
brr=br.join(dfc,how='inner');
brr0=brr.loc[brr['days_rel2birth']<=0].copy().set_index('tags_PID').copy();
brr1=brr.loc[brr['days_rel2birth']>=1].copy().set_index('tags_PID').copy();



def cp2num(mystr):
    if (mystr=='control'):
        return int(-1);
    if (mystr=='preterm'):
        return int(1);
    return int(0);

dfTags['this_KID']='';
dfTags['this_PID']='';
dfTags['tomID']=-1;

dfTags['class']=dfTags['class'].apply(cp2num);
dfTags=dfTags.loc[dfTags['class'] != 0].copy();
dfTags['class01']=(dfTags['class']+1)/2;
dfTags['KID_this']=np.nan;
dfTags['GA_delivery']=np.nan;
dfTags['GA']=np.nan;
dfTags['GA_clinical']=np.nan;
dfTags['days2birth']=np.nan;
dfTags['seqNum']=np.nan;
dfTags['seqNumRev']=np.nan;
dfTags['numVisits']=0;
dfTags['firstVisit']=np.nan;
dfTags['secondVisit']=np.nan;
dfTags['lastVisit']=np.nan;
dfTags['GA_next']=10000;
dfTags['KID_next']=np.nan;
dfTags['GA_prev']=-100;
dfTags['KID_prev']=np.nan;
dfTags['delta_next']=np.nan;
dfTags['delta_prev']=np.nan;
dfTags['has16S']=0;
dfTags['has16S_MRCD']=0;
dfTags['has16S_MCKD']=0;
dfTags['has16S_BRCD0']=0;
dfTags['has16S_BCKD0']=0;
dfTags['has16S_BRCD1']=0;
dfTags['has16S_BCKD1']=0;
dfTags['hasWMTS']=0;
dfTags['hasWMGS']=0;
dfTags['hasWMTSpwy']=0;
dfTags['hasWMGSpwy']=0;
dfTags['isAfAm']=0;
dfTags['isAfAmIsPTB']=0;

dfTags['isAsian']=0;
dfTags['isEuro']=0;
dfTags['isLatino']=0;
dfTags['isNative']=0;
dfTags['isEthnicOther']=0;
dfTags['isHawaii']=0;
dfTags['isAAC']=0;

dfTags['gravida']=np.nan;
dfTags['parity']=np.nan;
dfTags['preterm_history']=np.nan;
dfTags['miscarriage_history']=np.nan;
dfTags['no_of_live_births_history']=np.nan;
dfTags['AnyAntibiotic']=np.nan;
dfTags['AnyAntibioticToDate']=np.nan;
dfTags['hhq_miscarriage_or_stillbirth_binary']=np.nan;
dfTags['hhq_miscarriage_or_stillbirth_number']=np.nan;
dfTags['hhq_abortion_number']=np.nan;
dfTags['hhq_abortion_binary']=np.nan;
dfTags['hhq_miscarriage_or_stillbirth_sign']=0;
dfTags['hhq_abortion_sign']=0;
               
dfTags['short_cervix']=np.nan;
dfTags['vaginal_ph']=np.nan;
dfTags['cerclage']=np.nan;


dfTags['mdl_antibiotics']=np.nan;
dfTags['mdl_bmi']=np.nan;
dfTags['mld_progesterone']=np.nan;



dfTags['hhq_bv']='';
dfTags['hhq_bv_history']='';
dfTags['hhq_bv_lifetime_number']='';
dfTags['hhq_bv_recent']='';
dfTags['hhq_trichomoniasis_history']='';
dfTags['hhq_trichomoniasis_lifetime_number']='';
dfTags['hhq_trichomoniasis_recent']='';
dfTags['hhq_uti_history']='';
dfTags['hhq_uti_lifetime_number']='';
dfTags['hhq_uti_recent']='';

dfTags['baby_isBoy']=np.nan;
dfTags['baby_ICU']=np.nan;
dfTags['baby_weight']=np.nan;
dfTags['baby_apgar1']=np.nan;

dfTags['delivery_csection']=0;
dfTags['delivery_vaginal']=0;
dfTags['delivery_assisted']=0;





for PID,pList in kitsDates4PIDC.items():
    PIDn=int(PID[1:])
    tbI=0;
    seqNum=0;
    for lNode in pList:
        (toBirth,KID)=lNode;
        

        if ((KID in df16SpMV1D.index) and (KID in JN)):
            if (KID in df16SpMRCD.index):
                dfTags.set_value(KID,'has16S_MRCD',1);
            if (KID in df16SpMCKD.index):
                dfTags.set_value(KID,'has16S_MCKD',1);

            if (PID in bcc0.index):
                dfTags.set_value(KID,'has16S_BCKD0',1);
            if (PID in brr0.index):
                dfTags.set_value(KID,'has16S_BRCD0',1);
            if (PID in bcc1.index):
                dfTags.set_value(KID,'has16S_BCKD1',1);
            if (PID in brr1.index):
                dfTags.set_value(KID,'has16S_BRCD1',1);


            if (KID in dfHumAnn_wmts_pwy.index):
                dfTags.set_value(KID,'hasWMTSpwy',1);
            if (KID in dfHumAnn_wmgs_pwy.index):
                dfTags.set_value(KID,'hasWMGSpwy',1);
            if (KID in dfHumAnn_wmts_UR90l.index):
                dfTags.set_value(KID,'hasWMTS',1);
            if (KID in dfHumAnn_wmgs_UR90l.index):
                dfTags.set_value(KID,'hasWMGS',1);


            istb=0;
            gab=np.nan;
            dfTags.set_value(KID,'KID_this',KID);
            dfTags.set_value(KID,'this_KID',KID);
            dfTags.set_value(KID,'this_PID',PID);

                
            if (PID in dfDND.index):
                isAssisted=dfDND.at[PID,'dnd_delivery_assisted '];
                if (isAssisted[0]=='Y'):
                    dfTags.set_value(KID,'delivery_assisted',1);
                if (isAssisted[0]=='N'):
                    dfTags.set_value(KID,'delivery_assisted',-1);
                isCS=dfDND.at[PID,'dnd_delivery_csection'];
                if (isCS[0]=='Y'):
                    dfTags.set_value(KID,'delivery_csection',1);
                if (isCS[0]=='N'):
                    dfTags.set_value(KID,'delivery_csection',-1);
                isVaginal=dfDND.at[PID,'dnd_delivery_vaginal '];
                if (isVaginal[0]=='Y'):
                    dfTags.set_value(KID,'delivery_vaginal',1);
                if (isVaginal[0]=='N'):
                    dfTags.set_value(KID,'delivery_vaginal',-1);

                bbS=dfDND.at[PID,'dnd_baby1_sex '];
                if (bbS[0]=='F'):
                    dfTags.set_value(KID,'baby_isBoy',-1);
                if (bbS[0]=='M'):
                    dfTags.set_value(KID,'baby_isBoy',1);

                bbICU=dfDND.at[PID,'dnd_baby_icu'];
                if (bbICU[0]=='Y'):
                    dfTags.set_value(KID,'baby_ICU',1);
                if (bbICU[0:2]=='No'):
                    dfTags.set_value(KID,'baby_ICU',-1);

            if (PID in dfHHQ.index):
                isAAC=0;
                isAfAm=dfHHQ.at[PID,'hhq_african_american'];
                isEuro=dfHHQ.at[PID,'hhq_caucasian'];
                isEuroNum=0;
                isAfAmNum=0;
                if (isAfAm=='Yes'):
                    isAAC=1;
                    isAfAmNum=1;
                if (isEuro=='Yes'):
                    isEuroNum=1;
                    isAAC=1;

                if (dfHHQ.at[PID,'hhq_asian']=='Yes'): isAAC=1;
                if (dfHHQ.at[PID,'hhq_caucasian']=='Yes'): isAAC=1;

                if (dfHHQ.at[PID,'hhq_ethnic_other']=='Italian'): isAAC=10;isEuroNum=1;

                if (dfHHQ.at[PID,'hhq_ethnic_other']=='Moorish_American'): isAAC=1;isAfAmNum=1;

                dfTags.set_value(KID,'isAfAm',isAfAmNum);
                dfTags.set_value(KID,'isEuro',isEuroNum);
                dfTags.set_value(KID,'isAAC',isAAC);

                dfTags.set_value(KID,'isAsian', 1 if dfHHQ.at[PID,'hhq_asian']=='Yes' else 0);
                dfTags.set_value(KID,'isLatino', 1 if dfHHQ.at[PID,'hhq_hispanic_or_latino']=='Yes' else 0);
                dfTags.set_value(KID,'isNative', 1 if dfHHQ.at[PID,'hhq_american_indian_or_alaska_native']=='Yes' else 0);
                dfTags.set_value(KID,'isEthnicOther', 1 if dfHHQ.at[PID,'hhq_ethnic_other']=='Yes' else 0);
                dfTags.set_value(KID,'isHawaii', 1 if dfHHQ.at[PID,'hhq_native_hawaiian']=='Yes' else 0);


                dfTags.set_value(KID,'hhq_bv',yesno(dfHHQ.at[PID,'hhq_bv']));
                dfTags.set_value(KID,'hhq_bv_history',yesno(dfHHQ.at[PID,'hhq_bv_history']));
                dfTags.set_value(KID,'hhq_bv_lifetime_number',dfHHQ.at[PID,'hhq_bv_lifetime_number']);
                dfTags.set_value(KID,'hhq_bv_recent',dfHHQ.at[PID,'hhq_bv_recent']);
                dfTags.set_value(KID,'hhq_trichomoniasis_history',dfHHQ.at[PID,'hhq_trichomoniasis_history']);
                dfTags.set_value(KID,'hhq_trichomoniasis_lifetime_number',dfHHQ.at[PID,'hhq_trichomoniasis_lifetime_number']);
                dfTags.set_value(KID,'hhq_trichomoniasis_recent',dfHHQ.at[PID,'hhq_trichomoniasis_recent']);
                dfTags.set_value(KID,'hhq_uti_history',yesno(dfHHQ.at[PID,'hhq_uti_history']));
                dfTags.set_value(KID,'hhq_uti_lifetime_number',dfHHQ.at[PID,'hhq_uti_lifetime_number']);
                dfTags.set_value(KID,'hhq_uti_recent',yesno(dfHHQ.at[PID,'hhq_uti_recent']));


                hhq_mos=dfHHQ.at[PID,'hhq_miscarriage_or_stillbirth'];
                if (hhq_mos=='Yes'):
                    hhq_mos=1;
                    hhq_mos_sign=1;
                elif (hhq_mos=='No'):
                    hhq_mos=0;
                    hhq_mos_sign=-1;
                elif (hhq_mos=='Not_Sure'):
                    hhq_mos=0.1;
                    hhq_mos_sign=0;
                hhq_mosn=dfHHQ.at[PID,'hhq_miscarriage_or_stillbirth_number'];
                if (hhq_mosn=='NA'):
                    if (hhq_mos==0):
                        hhq_mosn=0;
                    else:
                        hhq_mosn=1.1;
                hhq_mosn=float(hhq_mosn);
                hhq_a=dfHHQ.at[PID,'hhq_abortion'];
                if (hhq_a=='Yes'):
                    hhq_a=1;
                    hhq_a_sign=1;
                elif (hhq_a=='No'):
                    hhq_a=0;
                    hhq_a_sign=-1;
                elif (hhq_a=='NA'):
                    hhq_a=0.5;
                    hhq_a_sign=0;
                hhq_an=dfHHQ.at[PID,'hhq_abortion_number'];
                if (hhq_an=='NA'):
                    if (hhq_a==0):
                        hhq_an=0;
                    else:
                        hhq_an=1.1;
                hhq_an=float(hhq_an);
                dfTags.set_value(KID,'hhq_miscarriage_or_stillbirth_binary',hhq_mos);
                dfTags.set_value(KID,'hhq_miscarriage_or_stillbirth_number',hhq_mosn);
                dfTags.set_value(KID,'hhq_abortion_number',hhq_an);
                dfTags.set_value(KID,'hhq_abortion_binary',hhq_a);
                dfTags.set_value(KID,'hhq_miscarriage_or_stillbirth_sign',hhq_mos_sign);
                dfTags.set_value(KID,'hhq_abortion_sign',hhq_a_sign);
            if (KID in dfSurvey.index):
                dfTags.set_value(KID,'AnyAntibiotic',dfSurvey.at[KID,'AnyAntibiotic']);
                dfTags.set_value(KID,'AnyAntibioticToDate',dfSurvey.at[KID,'AnyAntibioticToDate']);
            if (PID in dfMRAbs.index):
                shcvx=dfMRAbs.at[PID,'short_cervix'];
                print(shcvx);
                if (shcvx[0]=='Y'):
                    dfTags.set_value(KID,'short_cervix',1);
                if (shcvx[0]=='N'):
                    dfTags.set_value(KID,'short_cervix',-1);

                cclg=dfMRAbs.at[PID,'cerclage'];
                print(cclg);
                if (cclg[0]=='Y'):
                    dfTags.set_value(KID,'cerclage',1);
                if (cclg[0]=='N'):
                    dfTags.set_value(KID,'cerclage',-1);


                ph=dfClinical.at[KID,'vaginal_ph'];
                if (ph is not None):
                    ph=float(ph)
                    if (ph>0.1):
                        dfTags.set_value(KID,'vaginal_ph',ph);
                        



                bbw=dfMRAbs.at[PID,'birth_weight (gm)'];
                if (bbw<1):
                    bbw=1000*bbw;
                bba1=dfMRAbs.at[PID,'apgar_score_1min'];
                dfTags['baby_weight']=bbw;
                dfTags['baby_apgar1']=bba1;

                dfTags.set_value(KID,'baby_weight',gab);
                dfTags.set_value(KID,'baby_apgar1',gab);


                tgac=dfClinical.at[KID,'true_ga'];
                istb=dfMRAbs.at[PID,'isPreTerm'];
                gab=dfMRAbs.at[PID,'ga_at_delivery(days)'];
                ga=gab+toBirth;
                dfTags.set_value(KID,'GA_delivery',gab);
                dfTags.set_value(KID,'GA',ga);
                dfTags.set_value(KID,'GA_clinical',tgac);
                dfTags.set_value(KID,'days2birth',toBirth);
                dfTags.set_value(KID,'has16S',1);
                dfTags.set_value(KID,'gravida',float(dfMRAbs.at[PID,'gravida']));
                try:
                    dfTags.set_value(KID,'parity',float(dfMRAbs.at[PID,'parity']));
                except ValueError:
                    dfTags.set_value(KID,'parity',0.0);
                dfTags.set_value(KID,'preterm_history',float(dfMRAbs.at[PID,'preterm_history']));
                dfTags.set_value(KID,'miscarriage_history',float(dfMRAbs.at[PID,'miscarriage_history']));
                dfTags.set_value(KID,'no_of_live_births_history',float(dfMRAbs.at[PID,'no_of_live_births_history']));

            if (PID in dfBMI.index):
                dfTags.set_value(KID,'mdl_antibiotics',0.0);
                dfTags.set_value(KID,'mdl_bmi',30.587);
                dfTags.set_value(KID,'mld_progesterone',0.0);
                mab=dfBMI.at[PID,'model_antibiotics_anyclass'];
                bmi=dfBMI.at[PID,'BMI'];
                progestr=dfBMI.at[PID,'progesterone_pregnancy_selfreport'];
                print((mab,bmi,progestr));
                dfTags.set_value(KID,'mdl_antibiotics',mab);
                dfTags.set_value(KID,'mdl_bmi',bmi);
                dfTags.set_value(KID,'mld_progesterone',progestr);
            else:
                dfTags.set_value(KID,'mdl_antibiotics',0.0);
                dfTags.set_value(KID,'mdl_bmi',30.587);
                dfTags.set_value(KID,'mld_progesterone',0.0);


dfTags['isAfAmIsPTB']=dfTags['isAfAm']+2*dfTags['class01'];

dfRmedical=dfTags[['hhq_miscarriage_or_stillbirth_sign', 'hhq_abortion_sign','vaginal_ph','cerclage','short_cervix','mdl_antibiotics','mdl_bmi','mld_progesterone','gravida', 'parity', 'preterm_history', 'miscarriage_history', 'no_of_live_births_history',\
                     'hhq_miscarriage_or_stillbirth_binary','hhq_miscarriage_or_stillbirth_number','hhq_abortion_number','hhq_abortion_binary']].copy().fillna(0);
dfRmedical['gravidaBinary']=np.sign(dfRmedical['gravida']-1);
dfRmedical['parityBinary']=np.sign(dfRmedical['parity']);
dfRmedical['preterm_historyBinary']=np.sign(dfRmedical['preterm_history']);
dfRmedical['preterm_history3']=((dfRmedical['preterm_historyBinary']*2)-1)*dfRmedical['gravidaBinary'];
dfRmedical['miscarriage_historyBinary']=np.sign(dfRmedical['miscarriage_history']);
dfRmedical['miscarriage_historyRelative']=dfRmedical['miscarriage_history']/dfRmedical['gravida'];
dfRmedical['nonparityRelative']=1-(dfRmedical['parity'])/(dfRmedical['gravida']-1);
dfRmedical['gralesspar']=dfRmedical['gravida']-dfRmedical['parity'];
dfRmedical['no_of_live_births_historyBinary']=np.sign(dfRmedical['no_of_live_births_history']);
dfRmedical['abortionRelative']=dfRmedical['gravidaBinary']*(dfRmedical['hhq_abortion_number'])/(dfRmedical['gravida']-1);
dfRmedical['abortionRelative']=dfRmedical['abortionRelative'].apply( lambda x: 1 if (x>1) else x);
dfRmedical['miscarstillRelative']=dfRmedical['gravidaBinary']*(dfRmedical['hhq_miscarriage_or_stillbirth_number'])/(dfRmedical['gravida']-1);
dfRmedical['miscarstillRelative']=dfRmedical['miscarstillRelative'].apply( lambda x: 1 if (x>1) else x);
dfRmedical['abortion_history3']=((dfRmedical['hhq_abortion_binary']*2)-1)*dfRmedical['gravidaBinary'];
dfRmedical=dfRmedical[['hhq_miscarriage_or_stillbirth_sign', 'hhq_abortion_sign','vaginal_ph','cerclage','short_cervix','mdl_antibiotics','mdl_bmi','mld_progesterone','nonparityRelative','gralesspar','abortionRelative','miscarstillRelative','preterm_history3','abortion_history3','miscarriage_historyRelative','no_of_live_births_historyBinary','miscarriage_historyBinary','gravidaBinary','parityBinary','preterm_historyBinary','gravida','parity','preterm_history','miscarriage_history','no_of_live_births_history','hhq_miscarriage_or_stillbirth_binary','hhq_miscarriage_or_stillbirth_number','hhq_abortion_number','hhq_abortion_binary']].fillna(0).copy();

dfTags[['nonparityRelative','gralesspar','abortionRelative','miscarstillRelative','preterm_history3',
        'abortion_history3','miscarriage_historyRelative',
        'no_of_live_births_historyBinary','miscarriage_historyBinary',
        'gravidaBinary','parityBinary','preterm_historyBinary',
        'gravida','parity','preterm_history','miscarriage_history',
        'no_of_live_births_history','hhq_miscarriage_or_stillbirth_binary',
        'hhq_miscarriage_or_stillbirth_number','hhq_abortion_number','hhq_abortion_binary',\
        'hhq_miscarriage_or_stillbirth_sign', 'hhq_abortion_sign','vaginal_ph','cerclage',
        'short_cervix','mdl_antibiotics','mdl_bmi','mld_progesterone']]=dfRmedical[['nonparityRelative',
        'gralesspar','abortionRelative','miscarstillRelative','preterm_history3',
        'abortion_history3','miscarriage_historyRelative',
        'no_of_live_births_historyBinary','miscarriage_historyBinary',
        'gravidaBinary','parityBinary','preterm_historyBinary',
        'gravida','parity','preterm_history','miscarriage_history',
        'no_of_live_births_history','hhq_miscarriage_or_stillbirth_binary',
        'hhq_miscarriage_or_stillbirth_number','hhq_abortion_number','hhq_abortion_binary',\
         'hhq_miscarriage_or_stillbirth_sign', 'hhq_abortion_sign','vaginal_ph',   'cerclage',
         'short_cervix','mdl_antibiotics','mdl_bmi','mld_progesterone']]

dfTagsA=dfTags.loc[dfTags['has16S'] == 1].copy();
if (useMRCD):
    dfTagsA=dfTagsA.loc[dfTagsA['has16S_MRCD'] == 1].copy();
    dfTagsA=dfTagsA.loc[dfTagsA['has16S_MCKD'] == 1].copy();
if (useBRCD):
    dfTagsA=dfTagsA.loc[dfTagsA['has16S_BRCD1'] == 1].copy();
    dfTagsA=dfTagsA.loc[dfTagsA['has16S_BCKD1'] == 1].copy();
if (useWMTS):
    dfTagsA=dfTagsA.loc[dfTagsA['hasWMTSpwy'] == 1].copy();
    dfTagsA=dfTagsA.loc[dfTagsA['hasWMGSpwy'] == 1].copy();



dfTagsC=dfTagsA.loc[dfTagsA['GA']<maxUsedDay].copy();
if (useAAC is True):
    dfTagsC=dfTagsC.loc[dfTagsC['isAAC'] == 1].copy();

if (len(otherCond)>0):
    if (otherCond=='afam'):
        dfTagsC=dfTagsC.loc[dfTagsC['isAfAm']==1].copy();
    if (otherCond=='notafam'):
        dfTagsC=dfTagsC.loc[dfTagsC['isAfAm']==0].copy();

if (db16S=='99'):
    myList=list(df16sort.columns)[0:5];
    dfRNAW16b=0.0*df16sort[myList].copy();
    dfRNAW16c=0.0*df16sort[myList].copy();
    dfRNAW16o=0.0*df16sort[myList].copy();
    dfWMTSpwy=0.0*df16sort[myList].copy();
    dfWMGSpwy=0.0*df16sort[myList].copy();
    dfWMTS=0.0*df16sort[myList].copy();
    dfWMGS=0.0*df16sort[myList].copy();
if (db16S=='S2'):
    myList=list(df16SpMV1Dalt.columns)[0:5];
    dfRNAW16b=0.0*df16SpMV1Dalt[myList].copy();
    dfRNAW16c=0.0*df16SpMV1Dalt[myList].copy();
    dfRNAW16o=0.0*df16SpMV1Dalt[myList].copy();
    dfWMTSpwy=0.0*df16SpMV1Dalt[myList].copy();
    dfWMGSpwy=0.0*df16SpMV1Dalt[myList].copy();
    dfWMTS=0.0*df16SpMV1Dalt[myList].copy();
    dfWMGS=0.0*df16SpMV1Dalt[myList].copy();


if (useMRCD==False and useWMTS==False):
    (dfTagsIn,dfRNAW16,dfRNAW16o)=PandasUtils.innerNonJoin((dfTagsC.copy(),df16sort.copy(),df16SpMV1Dalt.copy()));
    
if (useMRCD==True and useMCKD==False and useWMTS==False):
    (dfTagsIn,dfRNAW16,dfRNAW16o,dfRNAW16b)=PandasUtils.innerNonJoin(
            (dfTagsC.copy(),df16sort.copy(),df16SpMV1Dalt.copy(),
            df16SpMRCD.copy()
            ));
if (useMRCD==True and useWMTS==False):
    (dfTagsIn,dfRNAW16,dfRNAW16o,dfRNAW16b,dfRNAW16c)=PandasUtils.innerNonJoin(
            (dfTagsC.copy(),df16sort.copy(),df16SpMV1Dalt.copy(),
            df16SpMRCD.copy(),df16SpMCKD.copy()
            ));
if (useMRCD==False and useWMTS==True):
    (dfTagsIn,dfRNAW16,dfRNAW16o,dfWMTSpwy,dfWMGSpwy,dfWMTS,dfWMGS)=PandasUtils.innerNonJoin(
            (dfTagsC.copy(),df16sort.copy(),df16SpMV1Dalt.copy(),
            dfHumAnn_wmts_pwy.copy(),dfHumAnn_wmgs_pwy.copy(),dfHumAnn_wmts_UR90l.copy(),dfHumAnn_wmgs_UR90l.copy()
            ));
if (useMRCD==True and useWMTS==True):
    (dfTagsIn,dfRNAW16,dfRNAW16o,dfRNAW16b,dfRNAW16c,dfWMTSpwy,dfWMGSpwy,dfWMTS,dfWMGS)=PandasUtils.innerNonJoin(
            (dfTagsC.copy(),df16sort.copy(),df16SpMV1Dalt.copy(),
            df16SpMRCD.copy(),df16SpMCKD.copy(),
            dfHumAnn_wmts_pwy.copy(),dfHumAnn_wmgs_pwy.copy(),dfHumAnn_wmts_UR90l.copy(),dfHumAnn_wmgs_UR90l.copy()
            ));


dfRNAW16['this_KID']=dfTagsIn['this_KID'];
dfRNAW16['this_PID']=dfTagsIn['this_PID'];

dfRNAW16o['this_KID']=dfTagsIn['this_KID'];
dfRNAW16b['this_KID']=dfTagsIn['this_KID'];
dfRNAW16c['this_KID']=dfTagsIn['this_KID'];
dfWMTSpwy['this_KID']=dfTagsIn['this_KID'];
dfWMGSpwy['this_KID']=dfTagsIn['this_KID'];
dfWMTS['this_KID']=dfTagsIn['this_KID'];
dfWMGS['this_KID']=dfTagsIn['this_KID'];

dfRNAW16o['this_PID']=dfTagsIn['this_PID'];
dfRNAW16b['this_PID']=dfTagsIn['this_PID'];
dfRNAW16c['this_PID']=dfTagsIn['this_PID'];
dfWMTSpwy['this_PID']=dfTagsIn['this_PID'];
dfWMGSpwy['this_PID']=dfTagsIn['this_PID'];
dfWMTS['this_PID']=dfTagsIn['this_PID'];
dfWMGS['this_PID']=dfTagsIn['this_PID'];

PIDnumVisit=defaultdict(int);
PIDfirstVisit=defaultdict(int);
PIDsecondVisit=defaultdict(int);
PIDlastVisit=defaultdict(int);
for PID,pList in kitsDates4PIDC.items():
    PIDnumVisit[PID]=0;
    PIDn=int(PID[1:])
    tbI=0;
    pGA=np.nan;
    pKID=np.nan;
    seqNum=0;
    for lNode in pList:
        (toBirth,KID)=lNode;
        if (KID in df16SpMV1D.index):
            istb=0;
            gab=np.nan;
            if (PID in dfMRAbs.index):
                tgac=float(dfClinical.at[KID,'true_ga']);
                istb=dfMRAbs.at[PID,'isPreTerm'];
                gab=dfMRAbs.at[PID,'ga_at_delivery(days)'];
                ga=gab+toBirth;
                if (tgac != ga):
                    print("MISMATCH:"+str(ga)+" "+str(tgac));
                if (ga<259):
                    PIDnumVisit[PID]=PIDnumVisit[PID]+1;
                    PIDlastVisit[PID]=ga;
                    if (seqNum==0):
                        PIDfirstVisit[PID]=ga;
                    if (seqNum==1):
                        PIDsecondVisit[PID]=ga;
                    seqNum=seqNum+1;

for PID,pList in kitsDates4PIDC.items():
    PIDn=int(PID[1:])
    tbI=0;
    pGA=np.nan;
    pKID=np.nan;
    seqNum=0;
    for lNode in pList:
        (toBirth,KID)=lNode;
        if (KID in dfTagsIn.index):
            istb=0;
            gab=np.nan;
            if (PID in dfMRAbs.index):
                tgac=dfClinical.at[KID,'true_ga'];
                istb=dfMRAbs.at[PID,'isPreTerm'];
                gab=dfMRAbs.at[PID,'ga_at_delivery(days)'];
                ga=gab+toBirth;
                nV=PIDnumVisit[PID];
                fV=PIDfirstVisit[PID];
                sV=PIDsecondVisit[PID];
                lV=PIDlastVisit[PID];
                dfTagsIn.set_value(KID,'seqNum',seqNum);
                dfTagsIn.set_value(KID,'seqNumRev',nV-seqNum-1);
                dfTagsIn.set_value(KID,'numVisits',nV);
                dfTagsIn.set_value(KID,'firstVisit',fV);
                dfTagsIn.set_value(KID,'secondVisit',sV);
                dfTagsIn.set_value(KID,'lastVisit',lV);
                seqNum=seqNum+1;
                if (np.isnan(pKID)==False):
                    dfTagsIn.set_value(pKID,'GA_next',ga);
                    dfTagsIn.set_value(pKID,'KID_next',int(KID));
                    dfTagsIn.set_value(pKID,'delta_next',ga-pGA);
                    dfTagsIn.set_value(KID,'GA_prev',pGA);
                    dfTagsIn.set_value(KID,'KID_prev',int(pKID));
                    dfTagsIn.set_value(KID,'delta_prev',ga-pGA);
                pKID=KID;
                pGA=ga;


if (binTypeIn[0:-1]=='oneF24'):
    t1s=42;t1e=168;
    t2s=168;t2e=245;
    
if (binTypeIn[0:-1]=='oneALL'):
    t1s=42;t1e=291;
    t2s=291;t2e=300;
    







binType=binTypeIn+otherCond;

tS=t1s;tE=t1e;
dfTemp=dfTagsIn.loc[dfTagsIn['GA']>=tS].copy();
dfTemp2=dfTemp.loc[dfTemp['GA_next']>=tE].copy();
dfTagsP1L=dfTemp2.loc[dfTemp2['GA']<tE].copy();
dfTagsP1sL=dfTagsP1L.sort_values(by='tags_PID').copy();
dfTagsP1siL=dfTagsP1sL.set_index('tags_PID').copy();

tS=t2s;tE=t2e;
dfTemp=dfTagsIn.loc[dfTagsIn['GA']>=tS].copy();
dfTemp2=dfTemp.loc[dfTemp['GA_next']>=tE].copy();
dfTagsP2L=dfTemp2.loc[dfTemp2['GA']<tE].copy();
dfTagsP2sL=dfTagsP2L.sort_values(by='tags_PID').copy();
dfTagsP2siL=dfTagsP2sL.set_index('tags_PID').copy()


tS=t1s;tE=t1e;
dfTemp=dfTagsIn.loc[dfTagsIn['GA']>=tS].copy();
dfTemp2=dfTemp.loc[dfTemp['GA_prev']<tS].copy();
dfTagsP1E=dfTemp2.loc[dfTemp2['GA']<tE].copy();
dfTagsP1sE=dfTagsP1E.sort_values(by='tags_PID').copy();
dfTagsP1siE=dfTagsP1sE.set_index('tags_PID').copy();

tS=t2s;tE=t2e;
dfTemp=dfTagsIn.loc[dfTagsIn['GA']>=tS].copy();
dfTemp2=dfTemp.loc[dfTemp['GA_prev']<tS].copy();
dfTagsP2E=dfTemp2.loc[dfTemp2['GA']<tE].copy();
dfTagsP2sE=dfTagsP2E.sort_values(by='tags_PID').copy();
dfTagsP2siE=dfTagsP2sE.set_index('tags_PID').copy()



if (goThree==1):
    tS=t3s;tE=t3e;
    dfTemp=dfTagsIn.loc[dfTagsIn['GA']>=tS].copy();
    dfTemp2=dfTemp.loc[dfTemp['GA_next']>=tE].copy();
    dfTagsP3L=dfTemp2.loc[dfTemp2['GA']<tE].copy();
    dfTagsP3sL=dfTagsP3L.sort_values(by='tags_PID').copy();
    dfTagsP3siL=dfTagsP3sL.set_index('tags_PID').copy()
    (dfTagsP1inj_E,dfTagsP2inj_E,dfTagsP1inj_L,dfTagsP2inj_L,dfTagsP3inj_L)=PandasUtils.innerNonJoin((dfTagsP1siE.copy(),dfTagsP2siE.copy(),dfTagsP1siL.copy(),dfTagsP2siL.copy(),dfTagsP3siL.copy()));
elif (goOne==1):
    (dfTagsP1inj_E,dfTagsP2inj_E,dfTagsP1inj_L,dfTagsP2inj_L)=PandasUtils.innerNonJoin((dfTagsP1siE.copy(),dfTagsP1siE.copy(),dfTagsP1siL.copy(),dfTagsP1siL.copy()));
else:
    (dfTagsP1inj_E,dfTagsP2inj_E,dfTagsP1inj_L,dfTagsP2inj_L)=PandasUtils.innerNonJoin((dfTagsP1siE.copy(),dfTagsP2siE.copy(),dfTagsP1siL.copy(),dfTagsP2siL.copy()));

pl=list(dfTagsP1inj_E.index);
pll=len(pl);
pli=np.random.permutation(pll);
dfTagsP1inj_E['tomID']=pli;
dfTagsP1inj_L['tomID']=pli;
dfTagsP2inj_E['tomID']=pli;
dfTagsP2inj_L['tomID']=pli;


dfTagsP1inj_E2=dfTagsP1inj_E.set_index('this_KID').copy();
dfTagsP1inj_L2=dfTagsP1inj_L.set_index('this_KID').copy();
dfTagsP2inj_E2=dfTagsP2inj_E.set_index('this_KID').copy();
dfTagsP2inj_L2=dfTagsP2inj_L.set_index('this_KID').copy();



(df_1_E_Tu,df_1_E_Ru,df_1_E_Rou,df_1_E_Rbu,df_1_E_Rcu,df_1_E_Rtpu,df_1_E_Rgpu,df_1_E_Rtu,df_1_E_Rgu)=PandasUtils.innerNonJoin((dfTagsP1inj_E2.copy(),
                                               dfRNAW16.copy(),dfRNAW16o.copy(),dfRNAW16b.copy(),dfRNAW16c.copy(),dfWMTSpwy.copy(),dfWMGSpwy.copy(),dfWMTS.copy(),dfWMGS.copy(),
                                               ));
(df_1_L_Tu,df_1_L_Ru,df_1_L_Rou,df_1_L_Rbu,df_1_L_Rcu,df_1_L_Rtpu,df_1_L_Rgpu,df_1_L_Rtu,df_1_L_Rgu)=PandasUtils.innerNonJoin((dfTagsP1inj_L2.copy(),
                                               dfRNAW16.copy(),dfRNAW16o.copy(),dfRNAW16b.copy(),dfRNAW16c.copy(),dfWMTSpwy.copy(),dfWMGSpwy.copy(),dfWMTS.copy(),dfWMGS.copy(),
                                               ));
(df_2_E_Tu,df_2_E_Ru,df_2_E_Rou,df_2_E_Rbu,df_2_E_Rcu,df_2_E_Rtpu,df_2_E_Rgpu,df_2_E_Rtu,df_2_E_Rgu)=PandasUtils.innerNonJoin((dfTagsP2inj_E2.copy(),
                                               dfRNAW16.copy(),dfRNAW16o.copy(),dfRNAW16b.copy(),dfRNAW16c.copy(),dfWMTSpwy.copy(),dfWMGSpwy.copy(),dfWMTS.copy(),dfWMGS.copy(),
                                               ));
(df_2_L_Tu,df_2_L_Ru,df_2_L_Rou,df_2_L_Rbu,df_2_L_Rcu,df_2_L_Rtpu,df_2_L_Rgpu,df_2_L_Rtu,df_2_L_Rgu)=PandasUtils.innerNonJoin((dfTagsP2inj_L2.copy(),
                                               dfRNAW16.copy(),dfRNAW16o.copy(),dfRNAW16b.copy(),dfRNAW16c.copy(),dfWMTSpwy.copy(),dfWMGSpwy.copy(),dfWMTS.copy(),dfWMGS.copy(),
                                               ));



df_1_E_R=df_1_E_Ru.sort_values(by='this_PID').set_index('this_PID').copy();
df_2_E_R=df_2_E_Ru.sort_values(by='this_PID').set_index('this_PID').copy();
df_1_L_R=df_1_L_Ru.sort_values(by='this_PID').set_index('this_PID').copy();
df_2_L_R=df_2_L_Ru.sort_values(by='this_PID').set_index('this_PID').copy();



df_1_E_Rb=df_1_E_Rbu.sort_values(by='this_PID').set_index('this_PID').copy();
df_2_E_Rb=df_2_E_Rbu.sort_values(by='this_PID').set_index('this_PID').copy();
df_1_L_Rb=df_1_L_Rbu.sort_values(by='this_PID').set_index('this_PID').copy();
df_2_L_Rb=df_2_L_Rbu.sort_values(by='this_PID').set_index('this_PID').copy();

df_1_E_Ro=df_1_E_Rou.sort_values(by='this_PID').set_index('this_PID').copy();
df_2_E_Ro=df_2_E_Rou.sort_values(by='this_PID').set_index('this_PID').copy();
df_1_L_Ro=df_1_L_Rou.sort_values(by='this_PID').set_index('this_PID').copy();
df_2_L_Ro=df_2_L_Rou.sort_values(by='this_PID').set_index('this_PID').copy();

df_1_E_Rc=df_1_E_Rcu.sort_values(by='this_PID').set_index('this_PID').copy();
df_2_E_Rc=df_2_E_Rcu.sort_values(by='this_PID').set_index('this_PID').copy();
df_1_L_Rc=df_1_L_Rcu.sort_values(by='this_PID').set_index('this_PID').copy();
df_2_L_Rc=df_2_L_Rcu.sort_values(by='this_PID').set_index('this_PID').copy();

df_1_E_Rtp=df_1_E_Rtpu.sort_values(by='this_PID').set_index('this_PID').copy();
df_2_E_Rtp=df_2_E_Rtpu.sort_values(by='this_PID').set_index('this_PID').copy();
df_1_L_Rtp=df_1_L_Rtpu.sort_values(by='this_PID').set_index('this_PID').copy();
df_2_L_Rtp=df_2_L_Rtpu.sort_values(by='this_PID').set_index('this_PID').copy();

df_1_E_Rt=df_1_E_Rtu.sort_values(by='this_PID').set_index('this_PID').copy();
df_2_E_Rt=df_2_E_Rtu.sort_values(by='this_PID').set_index('this_PID').copy();
df_1_L_Rt=df_1_L_Rtu.sort_values(by='this_PID').set_index('this_PID').copy();
df_2_L_Rt=df_2_L_Rtu.sort_values(by='this_PID').set_index('this_PID').copy();

df_1_E_Rgp=df_1_E_Rgpu.sort_values(by='this_PID').set_index('this_PID').copy();
df_2_E_Rgp=df_2_E_Rgpu.sort_values(by='this_PID').set_index('this_PID').copy();
df_1_L_Rgp=df_1_L_Rgpu.sort_values(by='this_PID').set_index('this_PID').copy();
df_2_L_Rgp=df_2_L_Rgpu.sort_values(by='this_PID').set_index('this_PID').copy();

df_1_E_Rg=df_1_E_Rgu.sort_values(by='this_PID').set_index('this_PID').copy();
df_2_E_Rg=df_2_E_Rgu.sort_values(by='this_PID').set_index('this_PID').copy();
df_1_L_Rg=df_1_L_Rgu.sort_values(by='this_PID').set_index('this_PID').copy();
df_2_L_Rg=df_2_L_Rgu.sort_values(by='this_PID').set_index('this_PID').copy();



df_1_E_T=df_1_E_Tu.sort_values(by='this_PID').set_index('this_PID').copy();
df_2_E_T=df_2_E_Tu.sort_values(by='this_PID').set_index('this_PID').copy();
df_1_L_T=df_1_L_Tu.sort_values(by='this_PID').set_index('this_PID').copy();
df_2_L_T=df_2_L_Tu.sort_values(by='this_PID').set_index('this_PID').copy();


if (useBRCD==False):
    brr1=0.0*df_1_E_Rb.copy();
    brr1.drop('this_KID', axis=1, inplace=True);
    bcc1=0.0*df_1_E_Rc.copy();
    bcc1.drop('this_KID', axis=1, inplace=True);
else:
    bcc1.sort_values('days_rel2birth',inplace=True);
    brr1.sort_values('days_rel2birth',inplace=True);
    bcc1 = bcc1[~bcc1.index.duplicated(keep='first')]
    brr1 = brr1[~brr1.index.duplicated(keep='first')]
    brr1.drop('days_rel2birth', axis=1, inplace=True);
    bcc1.drop('days_rel2birth', axis=1, inplace=True);

(df_1_E_T,df_1_E_R,df_1_E_Ro,df_1_E_Rb,df_1_E_Rc,df_1_E_Rbb,df_1_E_Rbc,df_1_E_Rtp,df_1_E_Rgp,df_1_E_Rt,df_1_E_Rg)=PandasUtils.innerNonJoin((
        df_1_E_T.copy(),df_1_E_R.copy(),df_1_E_Ro.copy(),df_1_E_Rb.copy(),df_1_E_Rc.copy(),brr1.copy(),bcc1.copy(),df_1_E_Rtp.copy(),df_1_E_Rgp.copy(),df_1_E_Rt.copy(),df_1_E_Rg.copy()
                                               ));
(df_1_L_T,df_1_L_R,df_1_L_Ro,df_1_L_Rb,df_1_L_Rc,df_1_L_Rbb,df_1_L_Rbc,df_1_L_Rtp,df_1_L_Rgp,df_1_L_Rt,df_1_L_Rg)=PandasUtils.innerNonJoin((
        df_1_L_T.copy(),df_1_L_R.copy(),df_1_L_Ro.copy(),df_1_L_Rb.copy(),df_1_L_Rc.copy(),brr1.copy(),bcc1.copy(),df_1_L_Rtp.copy(),df_1_L_Rgp.copy(),df_1_L_Rt.copy(),df_1_L_Rg.copy()
        ));
(df_2_E_T,df_2_E_R,df_2_E_Ro,df_2_E_Rb,df_2_E_Rc,df_2_E_Rbb,df_2_E_Rbc,df_2_E_Rtp,df_2_E_Rgp,df_2_E_Rt,df_2_E_Rg)=PandasUtils.innerNonJoin((
        df_2_E_T.copy(),df_2_E_R.copy(),df_2_E_Ro.copy(),df_2_E_Rb.copy(),df_2_E_Rc.copy(),brr1.copy(),bcc1.copy(),df_2_E_Rtp.copy(),df_2_E_Rgp.copy(),df_2_E_Rt.copy(),df_2_E_Rg.copy()
        ));
(df_2_L_T,df_2_L_R,df_2_L_Ro,df_2_L_Rb,df_2_L_Rc,df_2_L_Rbb,df_2_L_Rbc,df_2_L_Rtp,df_2_L_Rgp,df_2_L_Rt,df_2_L_Rg)=PandasUtils.innerNonJoin((
        df_2_L_T.copy(),df_2_L_R.copy(),df_2_L_Ro.copy(),df_2_L_Rb.copy(),df_2_L_Rc.copy(),brr1.copy(),bcc1.copy(),df_2_L_Rtp.copy(),df_2_L_Rgp.copy(),df_2_L_Rt.copy(),df_2_L_Rg.copy()
                                               ));



df_1_E_R.drop('this_KID', axis=1, inplace=True);
df_2_E_R.drop('this_KID', axis=1, inplace=True);
df_1_L_R.drop('this_KID', axis=1, inplace=True);
df_2_L_R.drop('this_KID', axis=1, inplace=True);

df_1_E_Rb.drop('this_KID', axis=1, inplace=True);
df_2_E_Rb.drop('this_KID', axis=1, inplace=True);
df_1_L_Rb.drop('this_KID', axis=1, inplace=True);
df_2_L_Rb.drop('this_KID', axis=1, inplace=True);

df_1_E_Ro.drop('this_KID', axis=1, inplace=True);
df_2_E_Ro.drop('this_KID', axis=1, inplace=True);
df_1_L_Ro.drop('this_KID', axis=1, inplace=True);
df_2_L_Ro.drop('this_KID', axis=1, inplace=True);

df_1_E_Rc.drop('this_KID', axis=1, inplace=True);
df_2_E_Rc.drop('this_KID', axis=1, inplace=True);
df_1_L_Rc.drop('this_KID', axis=1, inplace=True);
df_2_L_Rc.drop('this_KID', axis=1, inplace=True);

df_1_E_Rt.drop('this_KID', axis=1, inplace=True);
df_2_E_Rt.drop('this_KID', axis=1, inplace=True);
df_1_L_Rt.drop('this_KID', axis=1, inplace=True);
df_2_L_Rt.drop('this_KID', axis=1, inplace=True);

df_1_E_Rtp.drop('this_KID', axis=1, inplace=True);
df_2_E_Rtp.drop('this_KID', axis=1, inplace=True);
df_1_L_Rtp.drop('this_KID', axis=1, inplace=True);
df_2_L_Rtp.drop('this_KID', axis=1, inplace=True);

df_1_E_Rg.drop('this_KID', axis=1, inplace=True);
df_2_E_Rg.drop('this_KID', axis=1, inplace=True);
df_1_L_Rg.drop('this_KID', axis=1, inplace=True);
df_2_L_Rg.drop('this_KID', axis=1, inplace=True);

df_1_E_Rgp.drop('this_KID', axis=1, inplace=True);
df_2_E_Rgp.drop('this_KID', axis=1, inplace=True);
df_1_L_Rgp.drop('this_KID', axis=1, inplace=True);
df_2_L_Rgp.drop('this_KID', axis=1, inplace=True);


trueSCNT=df_1_E_R.shape[0];
if (trueSCNT>=10):

    sCnt=df_1_E_R.shape[0];
    fCnt=df_1_E_R.shape[1];

    allCols=list(df_1_E_R.columns);
    goodColsE=[];
    goodColsL=[];
    goodColsEL=[];
    idxI=0;

    colStats=np.zeros((fCnt,3));

    cI=0;
    for c in allCols:
        colStats[cI,2]=cI;
        nzc1e=(df_1_E_R[c] >= 0.001).sum();
        nzc2e=(df_2_E_R[c] >= 0.001).sum();
        nzc1l=(df_1_L_R[c] >= 0.001).sum();
        nzc2l=(df_2_L_R[c] >= 0.001).sum();

        nzc1eP=(df_1_E_R[c] >= 0.01).sum();
        nzc2eP=(df_2_E_R[c] >= 0.01).sum();
        nzc1lP=(df_1_L_R[c] >= 0.01).sum();
        nzc2lP=(df_2_L_R[c] >= 0.01).sum();

        colStats[cI,0]=(nzc1e)/sCnt;
        colStats[cI,1]=(nzc1eP)/sCnt;
        if (colStats[cI,0]>=0.15 or colStats[cI,1]>=0.05):
            goodColsEL.append(c);
        cI=cI+1;

    dfColStats=pandas.DataFrame(colStats,index=allCols)

    df_1_E_RF=df_1_E_R[goodColsEL].copy();
    df_2_E_RF=df_2_E_R[goodColsEL].copy();
    df_1_L_RF=df_1_L_R[goodColsEL].copy();
    df_2_L_RF=df_2_L_R[goodColsEL].copy();

    #if main 16S is 99, it should be used for reference
    if (db16S=='99'):
        goodColsV=goodColsEL.copy();

    ######

    ######
    sCnt=df_1_E_Ro.shape[0];
    fCnt=df_1_E_Ro.shape[1];
    allCols=list(df_1_E_Ro.columns);

    goodColsE=[];
    goodColsL=[];
    goodColsEL=[];
    idxI=0;

    colStats=np.zeros((fCnt,3));
    cI=0;
    for c in allCols:
        colStats[cI,2]=cI;
        nzc1e=(df_1_E_Ro[c] >= 0.001).sum();
        nzc2e=(df_2_E_Ro[c] >= 0.001).sum();
        nzc1l=(df_1_L_Ro[c] >= 0.001).sum();
        nzc2l=(df_2_L_Ro[c] >= 0.001).sum();

        nzc1eP=(df_1_E_Ro[c] >= 0.01).sum();
        nzc2eP=(df_2_E_Ro[c] >= 0.01).sum();
        nzc1lP=(df_1_L_Ro[c] >= 0.01).sum();
        nzc2lP=(df_2_L_Ro[c] >= 0.01).sum();

        colStats[cI,0]=(nzc1e)/sCnt;
        colStats[cI,1]=(nzc1eP)/sCnt;
        if (colStats[cI,0]>=0.15 or colStats[cI,1]>=0.05):
            goodColsEL.append(c);
        cI=cI+1;

    df_1_E_RFo=df_1_E_Ro[goodColsEL].copy();
    df_2_E_RFo=df_2_E_Ro[goodColsEL].copy();
    df_1_L_RFo=df_1_L_Ro[goodColsEL].copy();
    df_2_L_RFo=df_2_L_Ro[goodColsEL].copy();

    #if main 16S is S2, then the alt is 99, and should be used for reference
    if (db16S=='S2'):
        goodColsV=goodColsEL.copy();

    ######
    sCnt=2*df_1_E_Rb.shape[0];
    fCnt=df_1_E_Rb.shape[1];
    allCols=list(df_1_E_Rb.columns);

    idxI=0;
    goodColsAll=goodColsV.copy();

    for c in allCols:
        nzc1e=(df_1_E_Rb[c] >= 0.0001).sum();
        nzc2e=(df_2_E_Rb[c] >= 0.0001).sum();
        nzc1l=(df_1_L_Rb[c] >= 0.0001).sum();
        nzc2l=(df_2_L_Rb[c] >= 0.0001).sum();
        if (nzc1e+nzc2l>5):
            if (c not in goodColsV):
                goodColsAll.append(c);

    for c in goodColsV:
        if (c not in allCols):
            df_1_E_Rb[c]=0.0;
            df_1_L_Rb[c]=0.0;
            df_2_E_Rb[c]=0.0;
            df_2_L_Rb[c]=0.0;

    df_1_E_RFb=df_1_E_Rb[goodColsAll].copy();
    df_2_E_RFb=df_2_E_Rb[goodColsAll].copy();
    df_1_L_RFb=df_1_L_Rb[goodColsAll].copy();
    df_2_L_RFb=df_2_L_Rb[goodColsAll].copy();

    ######
    sCnt=2*df_1_E_Rc.shape[0];
    fCnt=df_1_E_Rc.shape[1];
    allCols=list(df_1_E_Rc.columns);

    idxI=0;
    goodColsAll=goodColsV.copy();

    for c in allCols:
        nzc1e=(df_1_E_Rc[c] >= 0.0001).sum();
        nzc2e=(df_2_E_Rc[c] >= 0.0001).sum();
        nzc1l=(df_1_L_Rc[c] >= 0.0001).sum();
        nzc2l=(df_2_L_Rc[c] >= 0.0001).sum();
        if (nzc1e+nzc2l>5):
            if (c not in goodColsV):
                goodColsAll.append(c);

    for c in goodColsV:
        if (c not in allCols):
            df_1_E_Rc[c]=0.0;
            df_1_L_Rc[c]=0.0;
            df_2_E_Rc[c]=0.0;
            df_2_L_Rc[c]=0.0;

    df_1_E_RFc=df_1_E_Rc[goodColsAll].copy();
    df_2_E_RFc=df_2_E_Rc[goodColsAll].copy();
    df_1_L_RFc=df_1_L_Rc[goodColsAll].copy();
    df_2_L_RFc=df_2_L_Rc[goodColsAll].copy();

    ######
    sCnt=2*df_1_E_Rbb.shape[0];
    fCnt=df_1_E_Rbb.shape[1];
    allCols=list(df_1_E_Rbb.columns);

    idxI=0;
    goodColsAll=goodColsV.copy();

    for c in allCols:
        nzc1e=(df_1_E_Rbb[c] >= 0.0001).sum();
        nzc2e=(df_2_E_Rbb[c] >= 0.0001).sum();
        nzc1l=(df_1_L_Rbb[c] >= 0.0001).sum();
        nzc2l=(df_2_L_Rbb[c] >= 0.0001).sum();
        if (nzc1e+nzc2l>5):
            if (c not in goodColsV):
                goodColsAll.append(c);

    for c in goodColsV:
        if (c not in allCols):
            df_1_E_Rbb[c]=0.0;
            df_1_L_Rbb[c]=0.0;
            df_2_E_Rbb[c]=0.0;
            df_2_L_Rbb[c]=0.0;

    df_1_E_RFbb=df_1_E_Rbb[goodColsAll].copy();
    df_2_E_RFbb=df_2_E_Rbb[goodColsAll].copy();
    df_1_L_RFbb=df_1_L_Rbb[goodColsAll].copy();
    df_2_L_RFbb=df_2_L_Rbb[goodColsAll].copy();

    ######
    sCnt=2*df_1_E_Rbc.shape[0];
    fCnt=df_1_E_Rbc.shape[1];
    allCols=list(df_1_E_Rbc.columns);

    idxI=0;
    goodColsAll=goodColsV.copy();

    for c in allCols:
        nzc1e=(df_1_E_Rbc[c] >= 0.0001).sum();
        nzc2e=(df_2_E_Rbc[c] >= 0.0001).sum();
        nzc1l=(df_1_L_Rbc[c] >= 0.0001).sum();
        nzc2l=(df_2_L_Rbc[c] >= 0.0001).sum();
        if (nzc1e+nzc2l>5):
            if (c not in goodColsV):
                goodColsAll.append(c);

    for c in goodColsV:
        if (c not in allCols):
            df_1_E_Rbc[c]=0.0;
            df_1_L_Rbc[c]=0.0;
            df_2_E_Rbc[c]=0.0;
            df_2_L_Rbc[c]=0.0;

    df_1_E_RFbc=df_1_E_Rbc[goodColsAll].copy();
    df_2_E_RFbc=df_2_E_Rbc[goodColsAll].copy();
    df_1_L_RFbc=df_1_L_Rbc[goodColsAll].copy();
    df_2_L_RFbc=df_2_L_Rbc[goodColsAll].copy();



    #
    #
    yAfAmClass_all,yAfAmClass_all_int,yAfAmClass_selIndices=PandasUtils.toNumpySetClassOnlyNonBinary(df_1_E_T,'isAfAmIsPTB');
    ya0=np.count_nonzero(yAfAmClass_all_int==0);
    ya1=np.count_nonzero(yAfAmClass_all_int==1);
    ya2=np.count_nonzero(yAfAmClass_all_int==2);
    ya3=np.count_nonzero(yAfAmClass_all_int==3);
    yaMin=np.min([ya0,ya1,ya2,ya3])


    dfCCC=df_1_E_T.loc[yAfAmClass_selIndices];
    y_all=dfCCC[['class']].as_matrix().astype(dtype='float32');


    www=np.round(np.sort(y_all.flatten())).astype(dtype='int32');

    print((np.count_nonzero(y_all==1),np.count_nonzero(y_all==0),np.count_nonzero(y_all==-1)))

    n_all=y_all.shape[0]
    y_all_int=y_all.astype(dtype='int32').reshape((n_all,));

    dfAFAM=df_1_E_T.loc[yAfAmClass_selIndices];
    afam_all=dfAFAM[['isAfAm']].as_matrix().astype(dtype='float32');
    afam_all_int=afam_all.astype(dtype='int32').reshape((n_all,));
    afam_all_indic=MLMultiClassUtils.toCategorical(afam_all_int,2);

    if yaMin>=5:
        afamPTB_cvcvTrains,afamPTB_cvcvTests=MLUtils.cvcvGen(yAfAmClass_all_int,1000,5);
    else:
        afamPTB_cvcvTrains,afamPTB_cvcvTests=MLUtils.cvcvGen(y_all_int,1000,5);

    yAfAm_all,yAfAm_all_int,yAfAm_selIndices=PandasUtils.toNumpySetClassOnlyNonBinary(df_1_E_T,'isAfAm');
    n=len(y_all_int)
    nPTB=np.count_nonzero(y_all_int+1);
    nAfAm=np.count_nonzero(yAfAm_all_int);

    fourway=2*(np.sign(y_all_int+1))+yAfAm_all_int

    pa=np.count_nonzero(fourway==3)
    pc=np.count_nonzero(fourway==2)
    ta=np.count_nonzero(fourway==1)
    tc=np.count_nonzero(fourway==0)

    dumpFileName='data-sets-manuscript/'+dataPrefix+'_casecontrol_16StimeV8_%s_input.pickle'%binType;

    cvcvFingerprint=afamPTB_cvcvTests[500][2][0];
    print('Set: %s Total: %d PTB(+1): %d TB(-1): %d AfAm: %d, NotAfAm: %d +AfAm:%d +NotAfAm:%d -AfAm:%d -NotAfAm:%d CVCV:%d'%(dumpFileName,n,nPTB,n-nPTB,nAfAm,n-nAfAm,pa,pc,ta,tc,cvcvFingerprint))

    if (justTry==False):


        fStat = open('data-sets/'+dataPrefix+'_casecontrol_16StimeV8_datasetsNumbers.txt', 'a');
        print('Set: %s Total: %d PTB(+1): %d TB(-1): %d AfAm: %d, NotAfAm: %d +AfAm:%d +NotAfAm:%d -AfAm:%d -NotAfAm:%d CVCV:%d'%(dumpFileName,n,nPTB,n-nPTB,nAfAm,n-nAfAm,pa,pc,ta,tc,cvcvFingerprint),file=fStat)
        fStat.close();

        allData=((df_1_E_T,df_2_E_T,df_1_L_T,df_2_L_T),
                 
                 (df_1_E_R,df_2_E_R,df_1_L_R,df_2_L_R),
                 (df_1_E_Rb,df_2_E_Rb,df_1_L_Rb,df_2_L_Rb),
                 (df_1_E_Rc,df_2_E_Rc,df_1_L_Rc,df_2_L_Rc),
                 (df_1_E_Rbb,df_2_E_Rbb,df_1_L_Rbb,df_2_L_Rbb),
                 (df_1_E_Rbc,df_2_E_Rbc,df_1_L_Rbc,df_2_L_Rbc),
                 (df_1_E_Ro,df_2_E_Ro,df_1_L_Ro,df_2_L_Ro),

                 (df_1_E_RF,df_2_E_RF,df_1_L_RF,df_2_L_RF),
                 (df_1_E_RFb,df_2_E_RFb,df_1_L_RFb,df_2_L_RFb),
                 (df_1_E_RFc,df_2_E_RFc,df_1_L_RFc,df_2_L_RFc),
                 (df_1_E_RFbb,df_2_E_RFbb,df_1_L_RFbb,df_2_L_RFbb),
                 (df_1_E_RFbc,df_2_E_RFbc,df_1_L_RFbc,df_2_L_RFbc),
                 (df_1_E_RFo,df_2_E_RFo,df_1_L_RFo,df_2_L_RFo),


                 (df_1_E_Rt,df_2_E_Rt,df_1_L_Rt,df_2_L_Rt),
                 (df_1_E_Rtp,df_2_E_Rtp,df_1_L_Rtp,df_2_L_Rtp),
                 (df_1_E_Rg,df_2_E_Rg,df_1_L_Rg,df_2_L_Rg),
                 (df_1_E_Rgp,df_2_E_Rgp,df_1_L_Rgp,df_2_L_Rgp)
                 );


        all_days_org=[t1s,t1e-1,t2s,t2e-1]

        pickle.dump( (allData,afamPTB_cvcvTrains,afamPTB_cvcvTests,afam_all_indic,yAfAmClass_selIndices,y_all,y_all_int,n_all,all_days_org), open( dumpFileName, "wb" ) );
