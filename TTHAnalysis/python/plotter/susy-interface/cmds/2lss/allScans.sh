#!/bin/bash

## IMPORTANT REMARK: In principle, one could run all models per config/region all together,
## but this will generate a lot of jobs which are all accessing the same background ROOT file.
## This kills memory and performance. Keep in mind to split as good and effective as possible
## to optimize performance.

S="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-18_ewktrees80X_M17_signals"
T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-25_ewkskims80X_M17_2lss_bkg;/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-23_ewkskims80X_M17_2lss_data;$S"
#O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-03-06_ewk80X_approvalScans_TE"
O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-03-13_ewk80X_postApprovalScans"

## full dataset
#LUMI=35.9
#BLIND="-X blinding"

## ichep
LUMI=12.9
BLIND=""

## queue
QUEUE="-q all.q"

## background only first
#python susy-interface/scanmaker.py 2lss 2lss $T $O -l $LUMI -o SR --models TChiSlepSnu_FD_05 --flags "$BLIND" --bkgOnly --redoBkg $QUEUE --mca susy-ewkino/2lss_v1/mca_forScan.txt

## TChiWH
###python susy-interface/scanmaker.py 2lss 2lss $T $O -l $LUMI -o SR --models TChiWH            --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpV_FL_data_comb --mpfr 2lss' --mca susy-ewkino/2lss_v1/mca_forScan.txt --check

## TChiWZ
###python susy-interface/scanmaker.py 2lssWZ 2lss $T $O -l $LUMI -o SR --models TChiWZ            --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpV_FL_data_comb --mpfr 2lss' --mca susy-ewkino/2lss_v1/mca_forScan.txt --check

## TChiSlepSnu_FD
##python susy-interface/scanmaker.py 2lss 2lss $T $O -l $LUMI -o SR --models TChiSlepSnu_FD    --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpV_FL_data_comb --mpfr 2lss' --mca susy-ewkino/2lss_v1/mca_forScan.txt --check

## TChiSlepSnu_FD_05
python susy-interface/scanmaker.py 2lss 2lss $T $O -l $LUMI -o SR --models TChiSlepSnu_FD_05 --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpV_FL_data_comb --mpfr 2lss' --mca susy-ewkino/2lss_v1/mca_forScan.txt --check

## TChiSlepSnu_FD_95
python susy-interface/scanmaker.py 2lss 2lss $T $O -l $LUMI -o SR --models TChiSlepSnu_FD_95 --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpV_FL_data_comb --mpfr 2lss' --mca susy-ewkino/2lss_v1/mca_forScan.txt --check



## TChiSlepSnu_TE
#python susy-interface/scanmaker.py 2lss 2lss $T $O -l $LUMI -o SR --models TChiSlepSnu_TE    --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpV_FL_data_comb --mpfr 2lss' --mca susy-ewkino/2lss_v1/mca_forScan.txt --check

## TChiSlepSnu_TE_05
#python susy-interface/scanmaker.py 2lss 2lss $T $O -l $LUMI -o SR --models TChiSlepSnu_TE_05 --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpV_FL_data_comb --mpfr 2lss' --mca susy-ewkino/2lss_v1/mca_forScan.txt --check

## TChiSlepSnu_TE_95
#python susy-interface/scanmaker.py 2lss 2lss $T $O -l $LUMI -o SR --models TChiSlepSnu_TE_95 --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpV_FL_data_comb --mpfr 2lss' --mca susy-ewkino/2lss_v1/mca_forScan.txt --check
