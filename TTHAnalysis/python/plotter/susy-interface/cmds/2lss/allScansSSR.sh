#!/bin/bash

## IMPORTANT REMARK: In principle, one could run all models per config/region all together,
## but this will generate a lot of jobs which are all accessing the same background ROOT file.
## This kills memory and performance. Keep in mind to split as good and effective as possible
## to optimize performance.

S="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-18_ewktrees80X_M17_signals"
T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-25_ewkskims80X_M17_2lss_bkg;/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-23_ewkskims80X_M17_2lss_data;$S"
O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-03-06_ewk80X_approvalScans"
LUMI=35.9
QUEUE="-q all.q"

## background only first
#python susy-interface/scanmaker.py 2lss 2lss $T $O -l $LUMI -o SSR1 --models TChiWZ --bkgOnly --redoBkg --noRFlags --flags '-X blinding' --mca susy-ewkino/2lss_v1/mca_forScan.txt --bins "1,0.5,1.5" --expr "SSR_2lss(nJet40,LepGood1_conePt,LepGood1_phi, LepGood2_conePt,LepGood2_phi, metPt, metPhi, LepGood1_charge)" $QUEUE
#python susy-interface/scanmaker.py 2lss 2lss $T $O -l $LUMI -o SSR2 --models TChiWZ --bkgOnly --redoBkg --noRFlags --flags '-X blinding' --mca susy-ewkino/2lss_v1/mca_forScan.txt --bins "1,0.5,1.5" --expr "SSR_2lss(nJet40,LepGood1_conePt,LepGood1_phi, LepGood2_conePt,LepGood2_phi, metPt, metPhi, LepGood1_charge)" $QUEUE

## TChiSlepSnu_FD
python susy-interface/scanmaker.py 2lss 2lss $T $O -l $LUMI -o SSR1 --models TChiSlepSnu_FD --sigOnly --noRFlags --flags '-X blinding' --mca susy-ewkino/2lss_v1/mca_forScan.txt --bins "1,0.5,1.5" --expr "SSR_2lss(nJet40,LepGood1_conePt,LepGood1_phi, LepGood2_conePt,LepGood2_phi, metPt, metPhi, LepGood1_charge)" --mca susy-ewkino/2lss_v1/mca_forScan.txt $QUEUE
python susy-interface/scanmaker.py 2lss 2lss $T $O -l $LUMI -o SSR2 --models TChiSlepSnu_FD --sigOnly --noRFlags --flags '-X blinding' --mca susy-ewkino/2lss_v1/mca_forScan.txt --bins "1,0.5,1.5" --expr "SSR_2lss(nJet40,LepGood1_conePt,LepGood1_phi, LepGood2_conePt,LepGood2_phi, metPt, metPhi, LepGood1_charge)" --mca susy-ewkino/2lss_v1/mca_forScan.txt $QUEUE

## TChiSlepSnu_FD_05
python susy-interface/scanmaker.py 2lss 2lss $T $O -l $LUMI -o SSR1 --models TChiSlepSnu_FD_05 --sigOnly --noRFlags --flags '-X blinding' --mca susy-ewkino/2lss_v1/mca_forScan.txt --bins "1,0.5,1.5" --expr "SSR_2lss(nJet40,LepGood1_conePt,LepGood1_phi, LepGood2_conePt,LepGood2_phi, metPt, metPhi, LepGood1_charge)" --mca susy-ewkino/2lss_v1/mca_forScan.txt $QUEUE
python susy-interface/scanmaker.py 2lss 2lss $T $O -l $LUMI -o SSR2 --models TChiSlepSnu_FD_05 --sigOnly --noRFlags --flags '-X blinding' --mca susy-ewkino/2lss_v1/mca_forScan.txt --bins "1,0.5,1.5" --expr "SSR_2lss(nJet40,LepGood1_conePt,LepGood1_phi, LepGood2_conePt,LepGood2_phi, metPt, metPhi, LepGood1_charge)" --mca susy-ewkino/2lss_v1/mca_forScan.txt $QUEUE

## TChiSlepSnu_FD_95
python susy-interface/scanmaker.py 2lss 2lss $T $O -l $LUMI -o SSR1 --models TChiSlepSnu_FD_95 --sigOnly --noRFlags --flags '-X blinding' --mca susy-ewkino/2lss_v1/mca_forScan.txt --bins "1,0.5,1.5" --expr "SSR_2lss(nJet40,LepGood1_conePt,LepGood1_phi, LepGood2_conePt,LepGood2_phi, metPt, metPhi, LepGood1_charge)" --mca susy-ewkino/2lss_v1/mca_forScan.txt $QUEUE
python susy-interface/scanmaker.py 2lss 2lss $T $O -l $LUMI -o SSR2 --models TChiSlepSnu_FD_95 --sigOnly --noRFlags --flags '-X blinding' --mca susy-ewkino/2lss_v1/mca_forScan.txt --bins "1,0.5,1.5" --expr "SSR_2lss(nJet40,LepGood1_conePt,LepGood1_phi, LepGood2_conePt,LepGood2_phi, metPt, metPhi, LepGood1_charge)" --mca susy-ewkino/2lss_v1/mca_forScan.txt $QUEUE
