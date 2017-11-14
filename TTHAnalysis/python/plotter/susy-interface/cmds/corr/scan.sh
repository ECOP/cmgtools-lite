#!/bin/bash

## IMPORTANT REMARK: In principle, one could run all models per config/region all together,
## but this will generate a lot of jobs which are all accessing the same background ROOT file.
## This kills memory and performance. Keep in mind to split as good and effective as possible
## to optimize performance.

S="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-18_ewktrees80X_M17_signals"
T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-25_ewkskims80X_M17_3l_bkg;/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-23_ewkskims80X_M17_3l_data;$S"
TT="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-25_ewkskims80X_M17_crwz_bkg;/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-23_ewkskims80X_M17_crwz_data;$S"
TTT="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-25_ewkskims80X_M17_2lss_bkg;/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-23_ewkskims80X_M17_2lss_data;$S"
O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-05-16_ewk80X_combinationCorrISR_oldSF/"
#O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-05-16_ewk80X_combinationCorrISR2_perBin"

## full dataset
LUMI=35.9
BLIND="-X blinding"

## queue
QUEUE="-q all.q --perBin"
#QUEUE=""



## BACKGROUND + DATA
## =====================================================

#python susy-interface/scanmaker.py crwzcorr crwz $TT $O -l $LUMI -o SR --models TChiWZ --bkgOnly --redoBkg --flags "$BLIND" $QUEUE
#python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $LUMI -o SR --models TChiWZ --bkgOnly --redoBkg --flags "$BLIND" --mca susy-ewkino/3lcorr/mca_forScan.txt $QUEUE



## WINOS
## =====================================================

## TChiWH
#python susy-interface/scanmaker.py crwzcorr crwz $TT $O -l $LUMI -o SR --models TChiWH --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpV_FL_data_comb --mpfr 2lss' --check
#python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $LUMI -o SR --models TChiWH --mca susy-ewkino/3lcorr/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check

## TChiWZ
#python susy-interface/scanmaker.py crwzcorrWZ crwz $TT $O -l $LUMI -o SR --models TChiWZ --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpV_FL_data_comb --mpfr 2lss' --check 
python susy-interface/scanmaker.py 3lcorrWZ 3lW $T $O -l $LUMI -o SR --models TChiWZ --mca susy-ewkino/3lcorr/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 675 --m2 1




