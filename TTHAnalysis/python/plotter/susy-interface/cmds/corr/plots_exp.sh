#!/bin/bash

S="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-18_ewktrees80X_M17_signals"
T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-25_ewkskims80X_M17_3l_bkg;/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-23_ewkskims80X_M17_3l_data;$S"
TT="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-25_ewkskims80X_M17_crwz_bkg;/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-23_ewkskims80X_M17_crwz_data;$S"
O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-05-16_ewk80X_wzcorr" # Do NOT give a trailing /
L=35.9
FL=35.9
QUEUE="--lspam 'Simulation'" #"-q all.q"
BLIND="-X blinding" #""

## SR binning
python susy-interface/plotmaker.py 3lcorr 3lW $T $O -l $L  --make mix --plots perCateg -o SR --flags "--perBin $BLIND" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE
#python susy-interface/plotmaker.py 3lcorr 3lW $T $O -l $L  --make mix --plots perCateg -o WZ --flags "--perBin $BLIND" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE 

## event variables
python susy-interface/plotmaker.py 3lcorr 3lW $T $O -l $L  --make mix --plots evt      -o SR --flags "--perBin $BLIND" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE
python susy-interface/plotmaker.py 3lcorr 3lW $T $O -l $L  --make mix --plots evt      -o WZ --flags "--perBin $BLIND" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE 

## lepton variables
python susy-interface/plotmaker.py 3lcorr 3lW $T $O -l $L  --make mix --plots lep      -o SR --flags "--perBin $BLIND" -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE
python susy-interface/plotmaker.py 3lcorr 3lW $T $O -l $L  --make mix --plots lep      -o WZ --flags "--perBin $BLIND" -p "prompt_WZ;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE 

