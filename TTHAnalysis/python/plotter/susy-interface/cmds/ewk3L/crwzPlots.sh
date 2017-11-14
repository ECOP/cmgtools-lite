#!/bin/bash

T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-25_ewkskims80X_M17_crwz_bkg;/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-23_ewkskims80X_M17_crwz_data;/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-18_ewktrees80X_M17_signals"
O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-05-12_ewk80X_crwzHT2" # Do NOT give a trailing /
#O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-05-09_ewk80X_crwz" # Do NOT give a trailing /
L=35.9
FL=35.9
QUEUE="" #"-q all.q"
BLIND="-X blinding -X 0or1jet" #""

#python susy-interface/plotmaker.py crwz crwz $T $O -l $L --make data --selPlots "htJet40j_new" -o amcatnlo --flags "$BLIND" $QUEUE
python susy-interface/plotmaker.py crwz crwz $T $O -l $L --make data --selPlots "htJet40j_new" -o powheg --flags "$BLIND" $QUEUE
#python susy-interface/plotmaker.py crwz crwz $T $O -l $L --make data --selPlots "nJet40_log;htJet40j_new;jet1Pt" -o powheg --flags "$BLIND" $QUEUE
#python susy-interface/plotmaker.py crwz crwz $T $O -l $L --make data --selPlots "nJet40_log;htJet40j_new;jet1Pt" -o amcatnlo --flags "$BLIND" $QUEUE

#python susy-interface/plotmaker.py crwz crwz $T $O -l $L --make data --selPlots "nJet40_log;htJet40j;jet1Pt" -o amcatnlo_ISR --flags "$BLIND" -W 'ISRweight*wBG' --sfriends ISRPtWeightModule $QUEUE

