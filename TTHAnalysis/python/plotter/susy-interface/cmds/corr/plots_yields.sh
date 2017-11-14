#!/bin/bash

S="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-18_ewktrees80X_M17_signals"
T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-25_ewkskims80X_M17_3l_bkg;/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-23_ewkskims80X_M17_3l_data;$S"
TT="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-25_ewkskims80X_M17_crwz_bkg;/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-23_ewkskims80X_M17_crwz_data;$S"
L=35.9
FL=35.9
QUEUE="" #"-q all.q"
BLIND="-X blinding" #""

## SR binning
O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-06-20_ewk80X_wzcorrbugfix2" # Do NOT give a trailing /
#python susy-interface/plotmaker.py 3lcorr 3lW $T $O -l $L  --make datasig --plots perCateg -o SR --flags "--perBin $BLIND" -p "fakes_appldata;promptsub;prompt_.*;rares_.*;convs;data;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE
#python susy-interface/plotmaker.py 3lcorr 3lW $T $O -l $L  --make datasig --selPlots htJet30j -o SR --flags "--perBin $BLIND" -p "prompt_WZ" $QUEUE
python susy-interface/plotmaker.py 3lcorr 3lW $T $O -l $L  --make datasig --plots evt      -o SR --flags "--perBin $BLIND" -p "fakes_appldata;promptsub;prompt_.*;rares_.*;convs;data;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE
python susy-interface/plotmaker.py 3lcorr 3lW $T $O -l $L  --make datasig --plots lep      -o SR --flags "--perBin $BLIND" -p "fakes_appldata;promptsub;prompt_.*;rares_.*;convs;data;sig_TChiWZ_200_100;sig_TChiWZ_350_250" $QUEUE

## CR WZ
WHAT="crwz_postFit"
O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-06-20_ewk80X_wzcorrbugfix2" # Do NOT give a trailing /
python susy-interface/plotmaker.py crwzcorr crwz $TT $O -l $L  --make data --plots evt   -o $WHAT --flags "--perBin $BLIND" -p "fakes_appldata;promptsub;prompt_.*;rares_.*;convs;data" $QUEUE
python susy-interface/plotmaker.py crwzcorr crwz $TT $O -l $L  --make data --plots lep   -o $WHAT --flags "--perBin $BLIND" -p "fakes_appldata;promptsub;prompt_.*;rares_.*;convs;data" $QUEUE


