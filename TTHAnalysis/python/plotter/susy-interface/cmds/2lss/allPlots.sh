#!/bin/bash

T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-25_ewkskims80X_M17_2lss_bkg;/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-23_ewkskims80X_M17_2lss_data;/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-18_ewktrees80X_M17_signals/"
O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-03-06_ewk80X_finalTest" # Do NOT give a trailing /
#O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-03-06_ewk80X_final" # Do NOT give a trailing /
O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-05-07_ewk80X_cwr4" # Do NOT give a trailing /
L=35.9
FL=35.9
QUEUE="" #"-q all.q"
BLIND="-X blinding --cmsprel ''" #""

## WITH LOGARITHMIC SCALE
## =================================================================

## signal regions and baseline regions plots
## -----------------------------------------------------------------

#python susy-interface/plotmaker.py 2lss 2lss $T $O -l $L  --make datasig --plots perCateg -o SRexerc --flags "--perBin $BLIND" -p "fakes_appldata;prompt_.*;rares_.*;convs;data;sig_TChiWH_150_24;sig_TChiWH_175_49;sig_TChiWH_175_1" $QUEUE
python susy-interface/plotmaker.py 2lss 2lss $T $O -l $L  --make datasig --plots perCateg -o SR --flags "--perBin $BLIND" -p "fakes_appldata;promptsub;flips_appldata;prompt_.*;rares_.*;convs;data;sig_TChiSlepSnu_05_500_350" $QUEUE
#python susy-interface/plotmaker.py 2lss 2lss $T $O -l $L --make datasig --selPlots METrebin -o SR0j -p "fakes_appldata;promptsub;flips_appldata;prompt_.*;rares_.*;convs;data;sig_TChiSlepSnu_05_500_350" --flags "--perBin $BLIND -A 0and1jet 0jet 'nJet40==0'" $QUEUE
#python susy-interface/plotmaker.py 2lss 2lss $T $O -l $L --make datasig --selPlots METrebin -o SR1j -p "fakes_appldata;promptsub;flips_appldata;prompt_.*;rares_.*;convs;data;sig_TChiSlepSnu_05_500_350" --flags "--perBin $BLIND -A 0and1jet 1jet 'nJet40==1'" $QUEUE
#python susy-interface/plotmaker.py 2lss 2lss $T $O -l $L --make datasig --plots evt         -o SR   -p "fakes_appldata;promptsub;flips_appldata;prompt_.*;rares_.*;convs;data;sig_TChiSlepSnu_05_500_350" --flags "--perBin $BLIND" $QUEUE

## data fakes
#python susy-interface/plotmaker.py 2lss 2lss $T $O -l $L --make data --plots perCateg -o SR --flags "--perBin $BLIND" $QUEUE
#python susy-interface/plotmaker.py 2lss 2lss $T $O -l $L --make data --plots evt      -o SR --flags "--perBin $BLIND" $QUEUE
#python susy-interface/plotmaker.py 2lss 2lss $T $O -l $L --make data --selPlots METrebin -o SR0j --flags "--perBin $BLIND -A 0and1jet 0jet 'nJet40==0'" $QUEUE
#python susy-interface/plotmaker.py 2lss 2lss $T $O -l $L --make data --selPlots METrebin -o SR1j --flags "--perBin $BLIND -A 0and1jet 1jet 'nJet40==1'" $QUEUE
