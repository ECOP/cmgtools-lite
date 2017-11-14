#!/bin/bash

S="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-18_ewktrees80X_M17_signals"
T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-25_ewkskims80X_M17_3l_bkg;/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-23_ewkskims80X_M17_3l_data;$S"
TT="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-25_ewkskims80X_M17_crwz_bkg;/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-23_ewkskims80X_M17_crwz_data;$S"
L=35.9
FL=35.9
QUEUE="-q all.q" #"-q all.q"
BLIND="-X blinding" #""

## DEFINE WHAT YOU WANNA DO
## =================================================================
## =================================================================

#DO="plots_all"
#DO="plots_jpt30to100"
#DO="plots_jpt100to200"
#DO="plots_jpt200toInf"
#DO="plots_ht0to100"
#DO="plots_ht100to200"
DO="plots_ht200toInf"
#DO="plots_eeemmm"
#DO="plots_eemmme"
#DO="plots_wz3"
#DO="plots_wz4"
#DO="plots_wz4p1"
#DO="plots_wz4p4"
#DO="plots_wz4p5"
#DO="plots_wz4p6"
#DO="wzcr"
#DO="pt30"
#DO="pt40"
#DO="pt50"
#DO="pt60"
#DO="wz3"
#DO="wz4"
#DO="wz4p1"
#DO="wz4p2"
#DO="wz4p3"
#DO="wz4p4"
#DO="wz4p5"
#DO="wz4p6"
#DO="wz4p6perBin"
#DO="wzx"
#DO="wz5"
#DO="wz5p4"



## PLOTS
## =================================================================
## =================================================================


## full selection
## -----------------------------------------------------------------

if [ $DO = "plots_all" ]; then
    O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-04-08_ewk80X_wzcorr" # Do NOT give a trailing /

    ## bkg versus sig
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots wzcorr -o bkgsig --flags "--perBin $BLIND" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE --sfriends leptonJetReCleanerSusyEWK3L_pt30 --srfriends leptonBuilderEWK
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots lep    -o bkgsig --flags "--perBin $BLIND" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE --sfriends leptonJetReCleanerSusyEWK3L_pt30 --srfriends leptonBuilderEWK
    
    ## wz versus sig
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots wzcorr -o wzsig  --flags "--perBin $BLIND" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE --sfriends leptonJetReCleanerSusyEWK3L_pt30 --srfriends leptonBuilderEWK
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots lep    -o wzsig  --flags "--perBin $BLIND" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE --sfriends leptonJetReCleanerSusyEWK3L_pt30 --srfriends leptonBuilderEWK
fi 


## full selection, jpt30to100
## -----------------------------------------------------------------

if [ $DO = "plots_jpt30to100" ]; then
    O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-04-08_ewk80X_wzcorr" # Do NOT give a trailing /

    ## bkg versus sig
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots wzcorr -o bkgsig_jpt30to100 --flags "--perBin $BLIND -A met isrjet 'nJetSel_Mini>0 && JetSel_Mini_pt[0]>=30 && JetSel_Mini_pt[0]<100'" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE --sfriends leptonJetReCleanerSusyEWK3L_pt30 --srfriends leptonBuilderEWK
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots lep    -o bkgsig_jpt30to100 --flags "--perBin $BLIND -A met isrjet 'nJetSel_Mini>0 && JetSel_Mini_pt[0]>=30 && JetSel_Mini_pt[0]<100'" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE --sfriends leptonJetReCleanerSusyEWK3L_pt30 --srfriends leptonBuilderEWK
    
    ## wz versus sig
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots wzcorr -o wzsig_jpt30to100  --flags "--perBin $BLIND -A met isrjet 'nJetSel_Mini>0 && JetSel_Mini_pt[0]>=30 && JetSel_Mini_pt[0]<100'" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE --sfriends leptonJetReCleanerSusyEWK3L_pt30 --srfriends leptonBuilderEWK
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots lep    -o wzsig_jpt30to100  --flags "--perBin $BLIND -A met isrjet 'nJetSel_Mini>0 && JetSel_Mini_pt[0]>=30 && JetSel_Mini_pt[0]<100'" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE --sfriends leptonJetReCleanerSusyEWK3L_pt30 --srfriends leptonBuilderEWK
fi 


## full selection, jpt100to200
## -----------------------------------------------------------------

if [ $DO = "plots_jpt100to200" ]; then
    O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-04-08_ewk80X_wzcorr" # Do NOT give a trailing /

    ## bkg versus sig
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots wzcorr -o bkgsig_jpt100to200 --flags "--perBin $BLIND -A met isrjet 'nJetSel_Mini>0 && JetSel_Mini_pt[0]>=100 && JetSel_Mini_pt[0]<200'" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE --sfriends leptonJetReCleanerSusyEWK3L_pt30 --srfriends leptonBuilderEWK
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots lep    -o bkgsig_jpt100to200 --flags "--perBin $BLIND -A met isrjet 'nJetSel_Mini>0 && JetSel_Mini_pt[0]>=100 && JetSel_Mini_pt[0]<200'" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE --sfriends leptonJetReCleanerSusyEWK3L_pt30 --srfriends leptonBuilderEWK
    
    ## wz versus sig
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots wzcorr -o wzsig_jpt100to200  --flags "--perBin $BLIND -A met isrjet 'nJetSel_Mini>0 && JetSel_Mini_pt[0]>=100 && JetSel_Mini_pt[0]<200'" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE --sfriends leptonJetReCleanerSusyEWK3L_pt30 --srfriends leptonBuilderEWK
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots lep    -o wzsig_jpt100to200  --flags "--perBin $BLIND -A met isrjet 'nJetSel_Mini>0 && JetSel_Mini_pt[0]>=100 && JetSel_Mini_pt[0]<200'" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE --sfriends leptonJetReCleanerSusyEWK3L_pt30 --srfriends leptonBuilderEWK
fi 


## full selection, pt200toInf
## -----------------------------------------------------------------

if [ $DO = "plots_jpt200toInf" ]; then
    O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-04-08_ewk80X_wzcorr" # Do NOT give a trailing /

    ## bkg versus sig
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots wzcorr -o bkgsig_jpt200toInf --flags "--perBin $BLIND -A met isrjet 'nJetSel_Mini>0 && JetSel_Mini_pt[0]>=200'" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE --sfriends leptonJetReCleanerSusyEWK3L_pt30 --srfriends leptonBuilderEWK
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots lep    -o bkgsig_jpt200toInf --flags "--perBin $BLIND -A met isrjet 'nJetSel_Mini>0 && JetSel_Mini_pt[0]>=200'" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE --sfriends leptonJetReCleanerSusyEWK3L_pt30 --srfriends leptonBuilderEWK
    
    ## wz versus sig
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots wzcorr -o wzsig_jpt200toInf  --flags "--perBin $BLIND -A met isrjet 'nJetSel_Mini>0 && JetSel_Mini_pt[0]>=200'" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE --sfriends leptonJetReCleanerSusyEWK3L_pt30 --srfriends leptonBuilderEWK
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots lep    -o wzsig_jpt200toInf  --flags "--perBin $BLIND -A met isrjet 'nJetSel_Mini>0 && JetSel_Mini_pt[0]>=200'" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE --sfriends leptonJetReCleanerSusyEWK3L_pt30 --srfriends leptonBuilderEWK
fi 


## full selection, ht0to100
## -----------------------------------------------------------------

if [ $DO = "plots_ht0to100" ]; then
    O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-04-08_ewk80X_wzcorr" # Do NOT give a trailing /

    ## bkg versus sig
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots evt    -o bkgsig_ht0to100 --flags "--perBin $BLIND -A met isrjet 'htJet30j_Mini>=0 && htJet30j_Mini<100'" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE 
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots lep    -o bkgsig_ht0to100 --flags "--perBin $BLIND -A met isrjet 'htJet30j_Mini>=0 && htJet30j_Mini<100'" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE 
    
    ## wz versus sig
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots evt    -o wzsig_ht0to100  --flags "--perBin $BLIND -A met isrjet 'htJet30j_Mini>=0 && htJet30j_Mini<100'" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE 
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots lep    -o wzsig_ht0to100  --flags "--perBin $BLIND -A met isrjet 'htJet30j_Mini>=0 && htJet30j_Mini<100'" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE 
fi 


## full selection, ht100to200
## -----------------------------------------------------------------

if [ $DO = "plots_ht100to200" ]; then
    O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-04-08_ewk80X_wzcorr" # Do NOT give a trailing /

    ## bkg versus sig
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots evt    -o bkgsig_ht100to200 --flags "--perBin $BLIND -A met isrjet 'htJet30j_Mini>=100 && htJet30j_Mini<200'" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE 
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots lep    -o bkgsig_ht100to200 --flags "--perBin $BLIND -A met isrjet 'htJet30j_Mini>=100 && htJet30j_Mini<200'" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE 
    
    ## wz versus sig
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots evt    -o wzsig_ht100to200  --flags "--perBin $BLIND -A met isrjet 'htJet30j_Mini>=100 && htJet30j_Mini<200'" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE 
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots lep    -o wzsig_ht100to200  --flags "--perBin $BLIND -A met isrjet 'htJet30j_Mini>=100 && htJet30j_Mini<200'" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE 
fi 


## full selection, ht200toInf
## -----------------------------------------------------------------

if [ $DO = "plots_ht200toInf" ]; then
    O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-04-08_ewk80X_wzcorr_paper" # Do NOT give a trailing /

    ## bkg versus sig
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots evt    -o bkgsig_ht200toInf --flags "--perBin $BLIND -A met isrjet 'htJet30j_Mini>=200'" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE --lspam "Simulation Preliminary" 
    #python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots lep    -o bkgsig_ht200toInf --flags "--perBin $BLIND -A met isrjet 'htJet30j_Mini>=200'" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE 
    
    ## wz versus sig
    #python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots evt    -o wzsig_ht200toInf  --flags "--perBin $BLIND -A met isrjet 'htJet30j_Mini>=200'" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE 
    #python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots lep    -o wzsig_ht200toInf  --flags "--perBin $BLIND -A met isrjet 'htJet30j_Mini>=200'" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE 
fi 


## full selection, eee, mmm
## -----------------------------------------------------------------

if [ $DO = "plots_eeemmm" ]; then
    O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-04-08_ewk80X_wzcorr" # Do NOT give a trailing /

    ## bkg versus sig
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots wzcorr -o bkgsig_eeemmm --flags "--perBin $BLIND -A met flav 'abs(LepSel_pdgId[0])+abs(LepSel_pdgId[1])+abs(LepSel_pdgId[2])==33 || abs(LepSel_pdgId[0])+abs(LepSel_pdgId[1])+abs(LepSel_pdgId[2])==39'" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE --sfriends leptonJetReCleanerSusyEWK3L_pt30 --srfriends leptonBuilderEWK
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots lep    -o bkgsig_eeemmm --flags "--perBin $BLIND -A met flav 'abs(LepSel_pdgId[0])+abs(LepSel_pdgId[1])+abs(LepSel_pdgId[2])==33 || abs(LepSel_pdgId[0])+abs(LepSel_pdgId[1])+abs(LepSel_pdgId[2])==39'" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE --sfriends leptonJetReCleanerSusyEWK3L_pt30 --srfriends leptonBuilderEWK
    
    ## wz versus sig
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots wzcorr -o wzsig_eeemmm  --flags "--perBin $BLIND -A met flav 'abs(LepSel_pdgId[0])+abs(LepSel_pdgId[1])+abs(LepSel_pdgId[2])==33 || abs(LepSel_pdgId[0])+abs(LepSel_pdgId[1])+abs(LepSel_pdgId[2])==39'" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE --sfriends leptonJetReCleanerSusyEWK3L_pt30 --srfriends leptonBuilderEWK
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots lep    -o wzsig_eeemmm  --flags "--perBin $BLIND -A met flav 'abs(LepSel_pdgId[0])+abs(LepSel_pdgId[1])+abs(LepSel_pdgId[2])==33 || abs(LepSel_pdgId[0])+abs(LepSel_pdgId[1])+abs(LepSel_pdgId[2])==39'" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE --sfriends leptonJetReCleanerSusyEWK3L_pt30 --srfriends leptonBuilderEWK
fi 


## full selection, eem, mme
## -----------------------------------------------------------------

if [ $DO = "plots_eemmme" ]; then
    O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-04-08_ewk80X_wzcorr" # Do NOT give a trailing /

    ## bkg versus sig
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots wzcorr -o bkgsig_eemmme --flags "--perBin $BLIND -A met flav 'abs(LepSel_pdgId[0])+abs(LepSel_pdgId[1])+abs(LepSel_pdgId[2])!=33 && abs(LepSel_pdgId[0])+abs(LepSel_pdgId[1])+abs(LepSel_pdgId[2])!=39'" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE --sfriends leptonJetReCleanerSusyEWK3L_pt30 --srfriends leptonBuilderEWK
    #python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots lep    -o bkgsig_eemmme --flags "--perBin $BLIND -A met flav 'abs(LepSel_pdgId[0])+abs(LepSel_pdgId[1])+abs(LepSel_pdgId[2])!=33 && abs(LepSel_pdgId[0])+abs(LepSel_pdgId[1])+abs(LepSel_pdgId[2])!=39'" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE --sfriends leptonJetReCleanerSusyEWK3L_pt30 --srfriends leptonBuilderEWK
    #
    ### wz versus sig
    #python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots wzcorr -o wzsig_eemmme  --flags "--perBin $BLIND -A met flav 'abs(LepSel_pdgId[0])+abs(LepSel_pdgId[1])+abs(LepSel_pdgId[2])!=33 && abs(LepSel_pdgId[0])+abs(LepSel_pdgId[1])+abs(LepSel_pdgId[2])!=39'" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE --sfriends leptonJetReCleanerSusyEWK3L_pt30 --srfriends leptonBuilderEWK
    #python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots lep    -o wzsig_eemmme  --flags "--perBin $BLIND -A met flav 'abs(LepSel_pdgId[0])+abs(LepSel_pdgId[1])+abs(LepSel_pdgId[2])!=33 && abs(LepSel_pdgId[0])+abs(LepSel_pdgId[1])+abs(LepSel_pdgId[2])!=39'" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE --sfriends leptonJetReCleanerSusyEWK3L_pt30 --srfriends leptonBuilderEWK
fi 


## full selection, wz3 binning
## -----------------------------------------------------------------

if [ $DO = "plots_wz3" ]; then
    O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-04-08_ewk80X_wzcorr" # Do NOT give a trailing /

    ## SR binning
    #python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --selPlots SR_wz3 -o bkgsig_wz3 --flags "--perBin $BLIND" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE
    #python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --selPlots SR_wz3 -o wzsig_wz3  --flags "--perBin $BLIND" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE 
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --selPlots SR_wz3 -o wzsig_wz3_norm  --flags "--plotmode norm --emptyStack --perBin $BLIND  --maxRatioRange 0.5 1.5 --fixRatioRange --ratioNums sig_TChiWZ_350_250 --ratioNums sig_TChiWZ_200_100 --ratioDen prompt_WZ" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE 


    #python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --selPlots SR_wz3 -o bkgsig_wz3 --flags "--perBin $BLIND" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE

    ### bkg versus sig
    #python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots wzcorr -o bkgsig_wz3 --flags "--perBin $BLIND" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE --sfriends leptonJetReCleanerSusyEWK3L_pt30 --srfriends leptonBuilderEWK
    #python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots lep    -o bkgsig_wz3 --flags "--perBin $BLIND" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE --sfriends leptonJetReCleanerSusyEWK3L_pt30 --srfriends leptonBuilderEWK
    #
    ### wz versus sig
    #python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots wzcorr -o wzsig_wz3  --flags "--perBin $BLIND" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE --sfriends leptonJetReCleanerSusyEWK3L_pt30 --srfriends leptonBuilderEWK
    #python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots lep    -o wzsig_wz3  --flags "--perBin $BLIND" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE --sfriends leptonJetReCleanerSusyEWK3L_pt30 --srfriends leptonBuilderEWK
fi 


## full selection, wz4 binning
## -----------------------------------------------------------------

if [ $DO = "plots_wz4" ]; then
    O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-04-08_ewk80X_wzcorr" # Do NOT give a trailing /

    ## SR binning
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --selPlots SR_wz4 -o bkgsig_wz4 --flags "--perBin $BLIND" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --selPlots SR_wz4 -o wzsig_wz4  --flags "--perBin $BLIND" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE 


    #python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --selPlots SR_wz3 -o bkgsig_wz3 --flags "--perBin $BLIND" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE

    ### bkg versus sig
    #python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots wzcorr -o bkgsig_wz3 --flags "--perBin $BLIND" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE --sfriends leptonJetReCleanerSusyEWK3L_pt30 --srfriends leptonBuilderEWK
    #python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots lep    -o bkgsig_wz3 --flags "--perBin $BLIND" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE --sfriends leptonJetReCleanerSusyEWK3L_pt30 --srfriends leptonBuilderEWK
    #
    ### wz versus sig
    #python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots wzcorr -o wzsig_wz3  --flags "--perBin $BLIND" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE --sfriends leptonJetReCleanerSusyEWK3L_pt30 --srfriends leptonBuilderEWK
    #python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots lep    -o wzsig_wz3  --flags "--perBin $BLIND" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE --sfriends leptonJetReCleanerSusyEWK3L_pt30 --srfriends leptonBuilderEWK
fi 


## full selection, wz4p1 binning
## -----------------------------------------------------------------

if [ $DO = "plots_wz4p1" ]; then
    O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-04-08_ewk80X_wzcorr" # Do NOT give a trailing /

    ## SR binning
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --selPlots SR_wz4p1 -o bkgsig_wz4 --flags "--perBin $BLIND" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --selPlots SR_wz4p1 -o wzsig_wz4  --flags "--perBin $BLIND" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE 

fi


## full selection, wz4p4 binning
## -----------------------------------------------------------------

if [ $DO = "plots_wz4p4" ]; then
    O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-04-08_ewk80X_wzcorr" # Do NOT give a trailing /

    ## SR binning
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --selPlots SR_wz4p4 -o bkgsig_wz4 --flags "--perBin $BLIND" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --selPlots SR_wz4p4 -o wzsig_wz4  --flags "--perBin $BLIND" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE 

fi


## full selection, wz4p5 binning
## -----------------------------------------------------------------

if [ $DO = "plots_wz4p5" ]; then
    O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-04-08_ewk80X_wzcorr" # Do NOT give a trailing /

    ## SR binning
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --selPlots SR_wz4p5 -o bkgsig_wz4 --flags "--perBin $BLIND" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE
    python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --selPlots SR_wz4p5 -o wzsig_wz4  --flags "--perBin $BLIND" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE 

fi


## full selection, wz4p6 binning
## -----------------------------------------------------------------

if [ $DO = "plots_wz4p6" ]; then
    O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-04-08_ewk80X_wzcorr" # Do NOT give a trailing /

    ## SR binning
    python susy-interface/plotmaker.py 3lcorr 3lA $T $O -l $L  --make mix --selPlots SR_wz4p6 -o bkgsig_wz4 --flags "--perBin $BLIND" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE
    python susy-interface/plotmaker.py 3lcorr 3lA $T $O -l $L  --make mix --selPlots SR_wz4p6 -o wzsig_wz4  --flags "--perBin $BLIND" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE 

fi


## LIMITS
## =================================================================
## =================================================================

## wz control region
## -----------------------------------------------------------------
if [ $DO = "wzcr" ]; then
    O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-04-08_ewk80X_wzCorrCards_wzcr"
    python susy-interface/scanmaker.py crwz crwz $TT $O -l $L -o SR --models TChiWZ --bkgOnly --redoBkg --flags "$BLIND" $QUEUE
    python susy-interface/scanmaker.py crwz crwz $TT $O -l $L -o SR --models TChiWZ --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpV_FL_data_comb --mpfr 2lss' --check --m1 350 --m2 250
    python susy-interface/scanmaker.py crwz crwz $TT $O -l $L -o SR --models TChiWZ --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpV_FL_data_comb --mpfr 2lss' --check --m1 300 --m2 200
    python susy-interface/scanmaker.py crwz crwz $TT $O -l $L -o SR --models TChiWZ --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpV_FL_data_comb --mpfr 2lss' --check --m1 250 --m2 150
    python susy-interface/scanmaker.py crwz crwz $TT $O -l $L -o SR --models TChiWZ --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpV_FL_data_comb --mpfr 2lss' --check --m1 200 --m2 100
fi


## full binning, 0/1 ISR jet, pt30
## -----------------------------------------------------------------

if [ $DO = "pt30" ]; then
    O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-04-08_ewk80X_wzCorrCards_1jet_pt30"
    #python susy-interface/scanmaker.py 3l "3lA1;3lA2" $T $O -l $L -o SR --models TChiWZ --bkgOnly --redoBkg --flags "$BLIND -A met isrjet 'nJet30_Mini<=1'" --mca susy-ewkino/3l/mca_forScan.txt $QUEUE --sfriends leptonJetReCleanerSusyEWK3L_pt30 --srfriends leptonBuilderEWK
    python susy-interface/scanmaker.py 3l "3lA1;3lA2" $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND -A met isjet 'nJet30_Mini<=1'" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 350 --m2 250 --sfriends leptonJetReCleanerSusyEWK3L_pt30 --srfriends leptonBuilderEWK
    python susy-interface/scanmaker.py 3l "3lA1;3lA2" $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND -A met isjet 'nJet30_Mini<=1'" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 300 --m2 200 --sfriends leptonJetReCleanerSusyEWK3L_pt30 --srfriends leptonBuilderEWK
    python susy-interface/scanmaker.py 3l "3lA1;3lA2" $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND -A met isjet 'nJet30_Mini<=1'" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 250 --m2 150 --sfriends leptonJetReCleanerSusyEWK3L_pt30 --srfriends leptonBuilderEWK
    python susy-interface/scanmaker.py 3l "3lA1;3lA2" $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND -A met isjet 'nJet30_Mini<=1'" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 200 --m2 100 --sfriends leptonJetReCleanerSusyEWK3L_pt30 --srfriends leptonBuilderEWK
fi


## full binning, 0/1 ISR jet, pt40
## -----------------------------------------------------------------

if [ $DO = "pt40" ]; then
    O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-04-08_ewk80X_wzCorrCards_1jet_pt40"
    #python susy-interface/scanmaker.py 3l "3lA1;3lA2" $T $O -l $L -o SR --models TChiWZ --bkgOnly --redoBkg --flags "$BLIND -A met isrjet 'nJet40_Mini<=1'" --mca susy-ewkino/3l/mca_forScan.txt $QUEUE --sfriends leptonJetReCleanerSusyEWK3L_pt40 --srfriends leptonBuilderEWK
    python susy-interface/scanmaker.py 3l "3lA1;3lA2" $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND -A met isjet 'nJet40_Mini<=1'" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 350 --m2 250 --sfriends leptonJetReCleanerSusyEWK3L_pt40 --srfriends leptonBuilderEWK
    python susy-interface/scanmaker.py 3l "3lA1;3lA2" $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND -A met isjet 'nJet40_Mini<=1'" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 300 --m2 200 --sfriends leptonJetReCleanerSusyEWK3L_pt40 --srfriends leptonBuilderEWK
    python susy-interface/scanmaker.py 3l "3lA1;3lA2" $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND -A met isjet 'nJet40_Mini<=1'" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 250 --m2 150 --sfriends leptonJetReCleanerSusyEWK3L_pt40 --srfriends leptonBuilderEWK
    python susy-interface/scanmaker.py 3l "3lA1;3lA2" $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND -A met isjet 'nJet40_Mini<=1'" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 200 --m2 100 --sfriends leptonJetReCleanerSusyEWK3L_pt40 --srfriends leptonBuilderEWK
fi


## full binning, 0/1 ISR jet, pt50
## -----------------------------------------------------------------

if [ $DO = "pt50" ]; then
    O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-04-08_ewk80X_wzCorrCards_1jet_pt50"
    #python susy-interface/scanmaker.py 3l "3lA1;3lA2" $T $O -l $L -o SR --models TChiWZ --bkgOnly --redoBkg --flags "$BLIND -A met isrjet 'nJet50_Mini<=1'" --mca susy-ewkino/3l/mca_forScan.txt $QUEUE --sfriends leptonJetReCleanerSusyEWK3L_pt50 --srfriends leptonBuilderEWK
    python susy-interface/scanmaker.py 3l "3lA1;3lA2" $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND -A met isjet 'nJet50_Mini<=1'" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 350 --m2 250 --sfriends leptonJetReCleanerSusyEWK3L_pt50 --srfriends leptonBuilderEWK
    python susy-interface/scanmaker.py 3l "3lA1;3lA2" $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND -A met isjet 'nJet50_Mini<=1'" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 300 --m2 200 --sfriends leptonJetReCleanerSusyEWK3L_pt50 --srfriends leptonBuilderEWK
    python susy-interface/scanmaker.py 3l "3lA1;3lA2" $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND -A met isjet 'nJet50_Mini<=1'" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 250 --m2 150 --sfriends leptonJetReCleanerSusyEWK3L_pt50 --srfriends leptonBuilderEWK
    python susy-interface/scanmaker.py 3l "3lA1;3lA2" $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND -A met isjet 'nJet50_Mini<=1'" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 200 --m2 100 --sfriends leptonJetReCleanerSusyEWK3L_pt50 --srfriends leptonBuilderEWK
fi


## full binning, 0/1 ISR jet, pt60
## -----------------------------------------------------------------

if [ $DO = "pt60" ]; then
    O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-04-08_ewk80X_wzCorrCards_1jet_pt60"
    #python susy-interface/scanmaker.py 3l "3lA1;3lA2" $T $O -l $L -o SR --models TChiWZ --bkgOnly --redoBkg --flags "$BLIND -A met isrjet 'nJet60_Mini<=1'" --mca susy-ewkino/3l/mca_forScan.txt $QUEUE --sfriends leptonJetReCleanerSusyEWK3L_pt60 --srfriends leptonBuilderEWK
    python susy-interface/scanmaker.py 3l "3lA1;3lA2" $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND -A met isjet 'nJet60_Mini<=1'" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 350 --m2 250 --sfriends leptonJetReCleanerSusyEWK3L_pt60 --srfriends leptonBuilderEWK
    python susy-interface/scanmaker.py 3l "3lA1;3lA2" $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND -A met isjet 'nJet60_Mini<=1'" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 300 --m2 200 --sfriends leptonJetReCleanerSusyEWK3L_pt60 --srfriends leptonBuilderEWK
    python susy-interface/scanmaker.py 3l "3lA1;3lA2" $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND -A met isjet 'nJet60_Mini<=1'" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 250 --m2 150 --sfriends leptonJetReCleanerSusyEWK3L_pt60 --srfriends leptonBuilderEWK
    python susy-interface/scanmaker.py 3l "3lA1;3lA2" $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND -A met isjet 'nJet60_Mini<=1'" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 200 --m2 100 --sfriends leptonJetReCleanerSusyEWK3L_pt60 --srfriends leptonBuilderEWK
fi


## full binning with HT and MT bin, wz3
## -----------------------------------------------------------------

if [ $DO = "wz3" ]; then
    O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-04-08_ewk80X_wzCorrCards_wz3binning"
    #python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --bkgOnly --redoBkg --flags "$BLIND" --mca susy-ewkino/3l/mca_forScan.txt $QUEUE
    python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 350 --m2 250
    python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 300 --m2 200
    python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 250 --m2 150
    python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 200 --m2 100
fi


## full binning with HT and MT bin, wz4
## -----------------------------------------------------------------

if [ $DO = "wz4" ]; then
    O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-04-08_ewk80X_wzCorrCards_wz4binning"
    #python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --bkgOnly --redoBkg --flags "$BLIND" --mca susy-ewkino/3l/mca_forScan.txt $QUEUE
    python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 350 --m2 250
    python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 300 --m2 200
    python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 250 --m2 150
    python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 200 --m2 100
fi


## full binning with HT and MT bin, wz4p1
## -----------------------------------------------------------------

if [ $DO = "wz4p1" ]; then
    O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-04-08_ewk80X_wzCorrCards_wz4p1binning"
    #python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --bkgOnly --redoBkg --flags "$BLIND" --mca susy-ewkino/3l/mca_forScan.txt $QUEUE
    python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 350 --m2 250
    python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 300 --m2 200
    python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 250 --m2 150
    python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 200 --m2 100
fi


## full binning with HT and MT bin, wz4p2
## -----------------------------------------------------------------

if [ $DO = "wz4p2" ]; then
    O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-04-08_ewk80X_wzCorrCards_wz4p2binning"
    #python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --bkgOnly --redoBkg --flags "$BLIND" --mca susy-ewkino/3l/mca_forScan.txt $QUEUE
    python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 350 --m2 250
    python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 300 --m2 200
    python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 250 --m2 150
    python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 200 --m2 100
fi


## full binning with HT and MT bin, wz4p3
## -----------------------------------------------------------------

if [ $DO = "wz4p3" ]; then
    O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-04-08_ewk80X_wzCorrCards_wz4p3binning"
    #python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --bkgOnly --redoBkg --flags "$BLIND" --mca susy-ewkino/3l/mca_forScan.txt $QUEUE
    python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 350 --m2 250
    python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 300 --m2 200
    python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 250 --m2 150
    python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 200 --m2 100
fi


## full binning with HT and MT bin, wz4p4
## -----------------------------------------------------------------

if [ $DO = "wz4p4" ]; then
    O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-04-08_ewk80X_wzCorrCards_wz4p4binning"
    #python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --bkgOnly --redoBkg --flags "$BLIND" --mca susy-ewkino/3l/mca_forScan.txt $QUEUE
    python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 350 --m2 250
    python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 300 --m2 200
    python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 250 --m2 150
    python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 200 --m2 100
fi


## full binning with HT and MT bin, wz4p5
## -----------------------------------------------------------------

if [ $DO = "wz4p5" ]; then
    O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-04-08_ewk80X_wzCorrCards_wz4p5binning_mcAtNlo"
    #python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --bkgOnly --redoBkg --flags "$BLIND" --mca susy-ewkino/3l/mca_forScan.txt $QUEUE
    python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly -q all.q --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check 
    #python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 350 --m2 250
    #python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 300 --m2 200
    #python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 250 --m2 150
    #python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 200 --m2 100
fi


## full binning with HT and MT bin, wz4p6
## -----------------------------------------------------------------

if [ $DO = "wz4p6" ]; then
    O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-04-08_ewk80X_wzCorrCards_wz4p6binning"
    #python susy-interface/scanmaker.py 3lcorr 3lA $T $O -l $L -o SR --models TChiWZ --bkgOnly --redoBkg --flags "$BLIND" --mca susy-ewkino/3l/mca_forScan.txt $QUEUE
    python susy-interface/scanmaker.py 3lcorr 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 350 --m2 250
    python susy-interface/scanmaker.py 3lcorr 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 300 --m2 200
    python susy-interface/scanmaker.py 3lcorr 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 250 --m2 150
    python susy-interface/scanmaker.py 3lcorr 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 200 --m2 100
fi


## full binning with HT and MT bin, wz4p6, perBin
## -----------------------------------------------------------------

if [ $DO = "wz4p6perBin" ]; then
    O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-04-08_ewk80X_wzCorrCards_wz4p6binningPB"
    python susy-interface/scanmaker.py 3lcorr 3lA $T $O -l $L -o SR --models TChiWZ --bkgOnly --redoBkg --flags "$BLIND" --mca susy-ewkino/3l/mca_forScan.txt $QUEUE --perBin
    #python susy-interface/scanmaker.py 3lcorr 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --perBin --m1 350 --m2 250
    #python susy-interface/scanmaker.py 3lcorr 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --perBin --m1 300 --m2 200
    #python susy-interface/scanmaker.py 3lcorr 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --perBin --m1 250 --m2 150
    #python susy-interface/scanmaker.py 3lcorr 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --perBin --m1 200 --m2 100
fi


## full binning with HT and MT bin, wz4x
## -----------------------------------------------------------------

if [ $DO = "wzx" ]; then
    O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-04-08_ewk80X_wzCorrCards_wzxbinning"
    #python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --bkgOnly --redoBkg --flags "$BLIND" --mca susy-ewkino/3l/mca_forScan.txt $QUEUE
    python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 350 --m2 250
    python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 300 --m2 200
    python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 250 --m2 150
    python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 200 --m2 100
fi


## full binning with HT and MT bin, wz5
## -----------------------------------------------------------------

if [ $DO = "wz5" ]; then
    O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-04-08_ewk80X_wzCorrCards_wz5binning"
    #python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --bkgOnly --redoBkg --flags "$BLIND" --mca susy-ewkino/3l/mca_forScan.txt $QUEUE
    python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 350 --m2 250
    python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 300 --m2 200
    python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 250 --m2 150
    python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 200 --m2 100
fi


## full binning with HT and MT bin, wz5p4
## -----------------------------------------------------------------

if [ $DO = "wz5p4" ]; then
    O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-04-08_ewk80X_wzCorrCards_wz5p4binning"
    #python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --bkgOnly --redoBkg --flags "$BLIND" --mca susy-ewkino/3l/mca_forScan.txt $QUEUE
    python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 350 --m2 250
    python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 300 --m2 200
    python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 250 --m2 150
    python susy-interface/scanmaker.py 3l 3lA $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 200 --m2 100
fi



