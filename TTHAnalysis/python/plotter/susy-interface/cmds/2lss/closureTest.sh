#!/bin/bash

#T="/scratch/cheidegg/Prod23Jan/"
T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-25_ewktrees80X_M17_bkg/"
#T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-25_ewkskims80X_M17_2lss_bkg;/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-23_ewkskims80X_M17_2lss_data"
O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-03-03_ewk80X_closureTest/" # do NOT give a trailing /
QUEUE="" #"-q all.q"
LUMI=35.9


## ttbar
#python susy-interface/plotmaker.py 2lss 2lss $T $O -l $LUMI --make bkgs --plots perCateg --flags '-X blinding -X filters -X trigger --ratioYLabel "Obs./pred." --plotmode closure --numerator matched_ttbar --fixRatioRange --maxRatioRange 0.5 1.5 --neg --perBin --ratioOffset 0.03 --print C,png,pdf,txt -e' --bkgs "matched_ttbar;appl_ttbar" --mca susy-ewkino/2lss_v1/mca_closure.txt --plot susy-ewkino/2lss_v1/plots_closure.txt --lspam 'Simulation' --noWeight --noCFlags  -o light_ttbar $QUEUE
#python susy-interface/plotmaker.py 2lss 2lss $T $O -l $LUMI --make bkgs --plots evt      --flags '-X blinding -X filters -X trigger --ratioYLabel "Obs./pred." --plotmode closure --numerator matched_ttbar --fixRatioRange --maxRatioRange 0.5 1.5 --neg --perBin --ratioOffset 0.03 --print C,png,pdf,txt -e' --bkgs "matched_ttbar;appl_ttbar" --mca susy-ewkino/2lss_v1/mca_closure.txt --plot susy-ewkino/2lss_v1/plots_closure.txt --lspam 'Simulation' --noWeight --noCFlags  -o light_ttbar $QUEUE

## wjets
#python susy-interface/plotmaker.py 2lss 2lss $T $O -l $LUMI --make bkgs --plots perCateg --flags '-X blinding -X filters -X trigger --ratioYLabel "Obs./pred." --plotmode closure --numerator matched_wjets --fixRatioRange --maxRatioRange 0.5 1.5 --neg --perBin --ratioOffset 0.03 --print C,png,pdf,txt -e' --bkgs "matched_wjets;appl_wjets" --mca susy-ewkino/2lss_v1/mca_closure.txt --plot susy-ewkino/2lss_v1/plots_closure.txt --lspam 'Simulation' --noWeight --noCFlags  -o light_wjets $QUEUE
python susy-interface/plotmaker.py 2lss 2lss $T $O -l $LUMI --make bkgs --plots evt      --flags '-X blinding -X filters -X trigger --ratioYLabel "Obs./pred." --plotmode closure --numerator matched_wjets --fixRatioRange --maxRatioRange 0.5 1.5 --neg --perBin --ratioOffset 0.03 --print C,png,pdf,txt -e' --bkgs "matched_wjets;appl_wjets" --mca susy-ewkino/2lss_v1/mca_closure.txt --plot susy-ewkino/2lss_v1/plots_closure.txt --lspam 'Simulation' --noWeight --noCFlags  -o light_wjets $QUEUE
