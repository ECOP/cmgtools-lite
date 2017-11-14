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
#DO="ssr1"
#DO="ssr2"
#DO="ssr3"
#DO="ssr4"
#DO="ssr5"
DO="ssr6"


## SSR HT>200 MT>100 MET>200
## -----------------------------------------------------------------

if [ $DO = "ssr1" ]; then
    O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-04-08_ewk80X_wzCorrCards_ssr1"
    #python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --bkgOnly --redoBkg --flags "$BLIND" --mca susy-ewkino/3lcorr/mca_forScan.txt $QUEUE --expr 'SuperSigcorr1(mllL_3l, mTL_3l, met_pt, htJet30j_Mini)' --bins '1,0.5,1.5'
    python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3lcorr/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --expr 'SuperSigcorr1(mllL_3l, mTL_3l, met_pt, htJet30j_Mini)' --bins '1,0.5,1.5' --m1 350 --m2 250
    python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3lcorr/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --expr 'SuperSigcorr1(mllL_3l, mTL_3l, met_pt, htJet30j_Mini)' --bins '1,0.5,1.5' --m1 300 --m2 200
    python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3lcorr/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --expr 'SuperSigcorr1(mllL_3l, mTL_3l, met_pt, htJet30j_Mini)' --bins '1,0.5,1.5' --m1 250 --m2 150
    python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3lcorr/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --expr 'SuperSigcorr1(mllL_3l, mTL_3l, met_pt, htJet30j_Mini)' --bins '1,0.5,1.5' --m1 200 --m2 100
fi


## SSR MT>160 MET>200
## -----------------------------------------------------------------

if [ $DO = "ssr2" ]; then
    O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-04-08_ewk80X_wzCorrCards_ssr2"
    #python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --bkgOnly --redoBkg --flags "$BLIND" --mca susy-ewkino/3lcorr/mca_forScan.txt $QUEUE --expr 'SuperSigcorr2(mllL_3l, mTL_3l, met_pt, htJet30j_Mini)' --bins '1,0.5,1.5'
    python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3lcorr/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --expr 'SuperSigcorr2(mllL_3l, mTL_3l, met_pt, htJet30j_Mini)' --bins '1,0.5,1.5' --m1 350 --m2 250
    python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3lcorr/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --expr 'SuperSigcorr2(mllL_3l, mTL_3l, met_pt, htJet30j_Mini)' --bins '1,0.5,1.5' --m1 300 --m2 200
    python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3lcorr/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --expr 'SuperSigcorr2(mllL_3l, mTL_3l, met_pt, htJet30j_Mini)' --bins '1,0.5,1.5' --m1 250 --m2 150
    python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3lcorr/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --expr 'SuperSigcorr2(mllL_3l, mTL_3l, met_pt, htJet30j_Mini)' --bins '1,0.5,1.5' --m1 200 --m2 100
fi


## SSR MET>200
## -----------------------------------------------------------------

if [ $DO = "ssr3" ]; then
    O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-04-08_ewk80X_wzCorrCards_ssr3"
    #python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --bkgOnly --redoBkg --flags "$BLIND" --mca susy-ewkino/3lcorr/mca_forScan.txt $QUEUE --expr 'SuperSigcorr3(mllL_3l, mTL_3l, met_pt, htJet30j_Mini)' --bins '1,0.5,1.5'
    python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3lcorr/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --expr 'SuperSigcorr3(mllL_3l, mTL_3l, met_pt, htJet30j_Mini)' --bins '1,0.5,1.5' --m1 350 --m2 250
    python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3lcorr/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --expr 'SuperSigcorr3(mllL_3l, mTL_3l, met_pt, htJet30j_Mini)' --bins '1,0.5,1.5' --m1 300 --m2 200
    python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3lcorr/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --expr 'SuperSigcorr3(mllL_3l, mTL_3l, met_pt, htJet30j_Mini)' --bins '1,0.5,1.5' --m1 250 --m2 150
    python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3lcorr/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --expr 'SuperSigcorr3(mllL_3l, mTL_3l, met_pt, htJet30j_Mini)' --bins '1,0.5,1.5' --m1 200 --m2 100
fi


## SSR onZ MET>200
## -----------------------------------------------------------------

if [ $DO = "ssr4" ]; then
    O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-04-08_ewk80X_wzCorrCards_ssr4"
    #python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --bkgOnly --redoBkg --flags "$BLIND" --mca susy-ewkino/3lcorr/mca_forScan.txt $QUEUE --expr 'SuperSigcorr4(mllL_3l, mTL_3l, met_pt, htJet30j_Mini)' --bins '1,0.5,1.5'
    python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3lcorr/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --expr 'SuperSigcorr4(mllL_3l, mTL_3l, met_pt, htJet30j_Mini)' --bins '1,0.5,1.5' --m1 350 --m2 250
    python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3lcorr/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --expr 'SuperSigcorr4(mllL_3l, mTL_3l, met_pt, htJet30j_Mini)' --bins '1,0.5,1.5' --m1 300 --m2 200
    python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3lcorr/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --expr 'SuperSigcorr4(mllL_3l, mTL_3l, met_pt, htJet30j_Mini)' --bins '1,0.5,1.5' --m1 250 --m2 150
    python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3lcorr/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --expr 'SuperSigcorr4(mllL_3l, mTL_3l, met_pt, htJet30j_Mini)' --bins '1,0.5,1.5' --m1 200 --m2 100
fi


## SSR onZ MET>200
## -----------------------------------------------------------------

if [ $DO = "ssr5" ]; then
    O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-04-08_ewk80X_wzCorrCards_ssr5"
    #python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --bkgOnly --redoBkg --flags "$BLIND" --mca susy-ewkino/3lcorr/mca_forScan.txt $QUEUE --expr 'SuperSigcorr5(mllL_3l, mTL_3l, met_pt, htJet30j_Mini)' --bins '1,0.5,1.5'
    python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3lcorr/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --expr 'SuperSigcorr5(mllL_3l, mTL_3l, met_pt, htJet30j_Mini)' --bins '1,0.5,1.5' --m1 350 --m2 250
    python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3lcorr/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --expr 'SuperSigcorr5(mllL_3l, mTL_3l, met_pt, htJet30j_Mini)' --bins '1,0.5,1.5' --m1 300 --m2 200
    python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3lcorr/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --expr 'SuperSigcorr5(mllL_3l, mTL_3l, met_pt, htJet30j_Mini)' --bins '1,0.5,1.5' --m1 250 --m2 150
    python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3lcorr/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --expr 'SuperSigcorr5(mllL_3l, mTL_3l, met_pt, htJet30j_Mini)' --bins '1,0.5,1.5' --m1 200 --m2 100
fi


## SSR onZ MET>200
## -----------------------------------------------------------------

if [ $DO = "ssr6" ]; then
    O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-04-08_ewk80X_wzCorrCards_ssr6"
    #python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --bkgOnly --redoBkg --flags "$BLIND" --mca susy-ewkino/3lcorr/mca_forScan.txt $QUEUE --expr 'SuperSigcorr6(mllL_3l, mTL_3l, met_pt, htJet30j_Mini)' --bins '1,0.5,1.5'
    python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3lcorr/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --expr 'SuperSigcorr6(mllL_3l, mTL_3l, met_pt, htJet30j_Mini)' --bins '1,0.5,1.5' --m1 350 --m2 250
    python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3lcorr/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --expr 'SuperSigcorr6(mllL_3l, mTL_3l, met_pt, htJet30j_Mini)' --bins '1,0.5,1.5' --m1 300 --m2 200
    python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3lcorr/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --expr 'SuperSigcorr6(mllL_3l, mTL_3l, met_pt, htJet30j_Mini)' --bins '1,0.5,1.5' --m1 250 --m2 150
    python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3lcorr/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --expr 'SuperSigcorr6(mllL_3l, mTL_3l, met_pt, htJet30j_Mini)' --bins '1,0.5,1.5' --m1 200 --m2 100
fi
