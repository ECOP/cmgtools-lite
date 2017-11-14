#!/bin/bash

T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-25_ewkskims80X_M17_3l_bkg;/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-23_ewkskims80X_M17_3l_data"
O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-03-02_ewk80X_arcFollowUp" # Do NOT give a trailing /
L=35.9
FL=35.9
QUEUE="" #"-q all.q"
BLIND="-X blinding" #""

## WITH LOGARITHMIC SCALE
## =================================================================

## signal regions and baseline regions plots
## -----------------------------------------------------------------

## data fakes
python susy-interface/plotmaker.py 3l "3lA" $T $O -l $L  --make data --selPlots lep_ptErrPt -o SR --flags "--perBin $BLIND -A bJetVeto onZ 'abs(mllL_3l-91)<15 && mTL_3l>100' -A met onlyMu 'abs(LepSel_pdgId[imTL_3l])==13'" $QUEUE
