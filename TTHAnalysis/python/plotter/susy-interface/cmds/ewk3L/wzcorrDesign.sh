#!/bin/bash

T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-25_ewkskims80X_M17_3l_bkg;/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-23_ewkskims80X_M17_3l_data;/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-18_ewktrees80X_M17_signals"
O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-03-26_ewk80X_wzcorr" # Do NOT give a trailing /
L=35.9
FL=35.9
QUEUE="" #"-q all.q"
BLIND="-X blinding" #""

## PLOTS
## -----------------------------------------------------------------

#python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --selPlots SRwz2 -o bkgsig --flags "--perBin $BLIND" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE
#python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --selPlots SRwz2 -o wzsig  --flags "--perBin $BLIND" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE
#
#python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --selPlots SRwz2 -o bkgsig1j --flags "--perBin $BLIND -A met isrjet 'nJet30_Mini<=1'" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE
#python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --selPlots SRwz2 -o wzsig1j  --flags "--perBin $BLIND -A met isrjet 'nJet30_Mini<=1'" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE


#python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots evt -o bkgsig --flags "--perBin $BLIND" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE
#python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots evt -o wzsig  --flags "--perBin $BLIND --ratioNums sig_TChiWZ_200_100,sig_TChiWZ_350_250 --ratioDen prompt_WZ" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE
#python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots lep -o bkgsig --flags "--perBin $BLIND" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE
#python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots lep -o wzsig  --flags "--perBin $BLIND" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE
#python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots evt -o wz200  --flags "--perBin $BLIND --plotmode norm" -p "prompt_WZ;sig_TChiWZ_200_100" $QUEUE
#python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots perCateg -o bkgsig --flags "--perBin $BLIND" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE
#python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots perCateg -o wzsig  --flags "--perBin $BLIND" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE
#
#python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots evt -o bkgsig1j --flags "--perBin $BLIND -A met isrjet 'nJet30_Mini<=1'" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE
#python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots evt -o wzsig1j  --flags "--perBin $BLIND -A met isrjet 'nJet30_Mini<=1'" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE
#python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots lep -o bkgsig1j --flags "--perBin $BLIND -A met isrjet 'nJet30_Mini<=1'" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE
#python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots lep -o wzsig1j  --flags "--perBin $BLIND -A met isrjet 'nJet30_Mini<=1'" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE

#python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots evt -o bkgsig1jonz --flags "--perBin $BLIND -A met isrjet 'nJet30_Mini<=1' -A met onZ 'nOSSFL_3l>=1 && abs(mllL_3l-91)<15'" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE
#python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots evt -o wzsig1jonz  --flags "--perBin $BLIND -A met isrjet 'nJet30_Mini<=1' -A met onZ 'nOSSFL_3l>=1 && abs(mllL_3l-91)<15'" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE --pretend

#python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --plots perCateg -o bkgsig1jonz --flags "--perBin $BLIND -A met isrjet 'nJet30_Mini<=1' -A met onZ 'nOSSFL_3l>=1 && mllL_3l>=75 && mllL_3l<105'" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE




## LIMITS
## -----------------------------------------------------------------
S="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-18_ewktrees80X_M17_signals"
T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-25_ewkskims80X_M17_3l_bkg;/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-23_ewkskims80X_M17_3l_data;$S"
TT="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-25_ewkskims80X_M17_crwz_bkg;/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-23_ewkskims80X_M17_crwz_data;$S"
O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-03-13_ewk80X_wzCorrCards_onZ1jet_try2"
O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-03-13_ewk80X_wzCorrCards_min2jet"
O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-03-13_ewk80X_wzCorrCards_0jet"
#O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-03-13_ewk80X_wzCorrCards_1jet_upperPtCut"

## full dataset
LUMI=35.9
BLIND="-X blinding"

## queue
#QUEUE="--pretend --perBin"
QUEUE="-q all.q"
#QUEUE=""



## SRwz2
## =====================================================
#O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-03-13_ewk80X_wzCorrCards_SRwz2"
#python susy-interface/scanmaker.py 3l 3lA $T $O -l $LUMI -o SR --models TChiWZ --bkgOnly --redoBkg --flags "$BLIND" --mca susy-ewkino/3l/mca_forScan.txt $QUEUE
O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-03-13_ewk80X_wzCorrCards_SRwz2_1j"
#python susy-interface/scanmaker.py 3l 3lA $T $O -l $LUMI -o SR --models TChiWZ --bkgOnly --redoBkg --flags "$BLIND -A met isrjet 'nJet30_Mini<=1'" --mca susy-ewkino/3l/mca_forScan.txt $QUEUE
#python susy-interface/scanmaker.py 3l 3lA $T $O -l $LUMI -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND -A met isjet 'nJet30_Mini<=1'" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 350 --m2 250
#python susy-interface/scanmaker.py 3l 3lA $T $O -l $LUMI -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND -A met isjet 'nJet30_Mini<=1'" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 300 --m2 200
#python susy-interface/scanmaker.py 3l 3lA $T $O -l $LUMI -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND -A met isjet 'nJet30_Mini<=1'" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 250 --m2 150
#python susy-interface/scanmaker.py 3l 3lA $T $O -l $LUMI -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND -A met isjet 'nJet30_Mini<=1'" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 200 --m2 100



## BACKGROUND + DATA
## =====================================================

#python susy-interface/scanmaker.py crwz crwz $TT $O -l $LUMI -o SR --models TChiWZ --bkgOnly --redoBkg --flags "$BLIND" $QUEUE
#python susy-interface/scanmaker.py 3l wzcorr $T $O -l $LUMI -o SR --models TChiWZ --bkgOnly --redoBkg --flags "$BLIND" --mca susy-ewkino/3l/mca_forScan.txt $QUEUE
#python susy-interface/scanmaker.py 3l "3lA1;3lA2" $T $O -l $LUMI -o SR --models TChiWZ --bkgOnly --redoBkg --flags "$BLIND -A met isrjet 'nJet30_Mini<=1'" --mca susy-ewkino/3l/mca_forScan.txt $QUEUE --pretend
#python susy-interface/scanmaker.py 3l "3lA1;3lA2" $T $O -l $LUMI -o SR --models TChiWZ --bkgOnly --redoBkg --flags "$BLIND -A met isrjet 'nJet30_Mini<=1' -A met upperPtCut 'LepSel1_pt<200 && LepSel2_pt<150 && LepSel3_pt<150" --mca susy-ewkino/3l/mca_forScan.txt $QUEUE


## SIGNAL
## =====================================================

## CRWZ
#python susy-interface/scanmaker.py crwz crwz $TT $O -l $LUMI -o SR --models TChiWZ --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpV_FL_data_comb --mpfr 2lss' --check --m1 350 --m2 250 --pretend
#python susy-interface/scanmaker.py crwz crwz $TT $O -l $LUMI -o SR --models TChiWZ --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpV_FL_data_comb --mpfr 2lss' --check --m1 300 --m2 200
#python susy-interface/scanmaker.py crwz crwz $TT $O -l $LUMI -o SR --models TChiWZ --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpV_FL_data_comb --mpfr 2lss' --check --m1 250 --m2 150
#python susy-interface/scanmaker.py crwz crwz $TT $O -l $LUMI -o SR --models TChiWZ --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpV_FL_data_comb --mpfr 2lss' --check --m1 200 --m2 100
#
### WZCORR
#python susy-interface/scanmaker.py 3l wzcorr $T $O -l $LUMI -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 350 --m2 250
#python susy-interface/scanmaker.py 3l wzcorr $T $O -l $LUMI -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 200 --m2 100

### 3lA full binning
#python susy-interface/scanmaker.py 3l "3lA1;3lA2" $T $O -l $LUMI -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND -A met isjet 'nJet30_Mini<=0'" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 350 --m2 250 --pretend
#python susy-interface/scanmaker.py 3l "3lA1;3lA2" $T $O -l $LUMI -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND -A met isjet 'nJet30_Mini<=0'" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 300 --m2 200
#python susy-interface/scanmaker.py 3l "3lA1;3lA2" $T $O -l $LUMI -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND -A met isjet 'nJet30_Mini<=0'" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 250 --m2 150
#python susy-interface/scanmaker.py 3l "3lA1;3lA2" $T $O -l $LUMI -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND -A met isjet 'nJet30_Mini<=0'" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 200 --m2 100

### 3lA full binning
#python susy-interface/scanmaker.py 3l "3lA1;3lA2" $T $O -l $LUMI -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND -A met isjet 'nJet30_Mini<=1'" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 350 --m2 250
#python susy-interface/scanmaker.py 3l "3lA1;3lA2" $T $O -l $LUMI -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND -A met isjet 'nJet30_Mini<=1'" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 300 --m2 200
#python susy-interface/scanmaker.py 3l "3lA1;3lA2" $T $O -l $LUMI -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND -A met isjet 'nJet30_Mini<=1'" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 250 --m2 150
#python susy-interface/scanmaker.py 3l "3lA1;3lA2" $T $O -l $LUMI -o SR --models TChiWZ --mca susy-ewkino/3l/mca_forScan.txt --flags "$BLIND -A met isjet 'nJet30_Mini<=1'" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 200 --m2 100



