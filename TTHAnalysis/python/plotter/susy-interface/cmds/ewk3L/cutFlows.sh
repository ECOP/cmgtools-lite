#!/bin/bash

T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-25_ewkskims80X_M17_3l_bkg;/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-23_ewkskims80X_M17_3l_data;/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-18_ewktrees80X_M17_signals"
O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-04-24_ewk80X_cutflows" # Do NOT give a trailing /
L=35.9
FL=35.9
QUEUE="" #"-q all.q"
BLIND="-X blinding" #""

### TChiWZ 200/100
#python susy-interface/accmaker.py 3l 3lA $T $O -l $L -o SR --flags '$BLIND -A met highMT "mTL_3l>100" -A met aboveZ "mllL_3l>75"' -p "sig_TChiWZ_200_100" --noRFlags --cuts susy-ewkino/3l/cuts_cutflow.txt
#mv /afs/cern.ch/user/c/cheidegg/www/heppy/2017-04-24_ewk80X_cutflows/acc/SR/3lA/35p9fb/accmap_SR/3lA_3lA.txt /afs/cern.ch/user/c/cheidegg/www/heppy/2017-04-24_ewk80X_cutflows/acc/SR/3lA/35p9fb/accmap_SR/3lA_sig_TChiWZ_200_100.txt
#
### TChiWZ 500/150
#python susy-interface/accmaker.py 3l 3lA $T $O -l $L -o SR --flags '$BLIND -A met highMT "mTL_3l>100" -A met aboveZ "mllL_3l>75"' -p "sig_TChiWZ_500_150" --noRFlags --cuts susy-ewkino/3l/cuts_cutflow.txt
#mv /afs/cern.ch/user/c/cheidegg/www/heppy/2017-04-24_ewk80X_cutflows/acc/SR/3lA/35p9fb/accmap_SR/3lA_3lA.txt /afs/cern.ch/user/c/cheidegg/www/heppy/2017-04-24_ewk80X_cutflows/acc/SR/3lA/35p9fb/accmap_SR/3lA_sig_TChiWZ_500_150.txt

### TChiSlepSnuTD 250/150
#python susy-interface/accmaker.py 3l 3lA $T $O -l $L -o SR --flags '$BLIND -A met lowMTT "mT2T_3l<100" -A met belowZ "mllL_3l<75"' -p "sig_TChiSlepSnu_TD_250_150" --noRFlags --cuts susy-ewkino/3l/cuts_cutflow.txt
#mv /afs/cern.ch/user/c/cheidegg/www/heppy/2017-04-24_ewk80X_cutflows/acc/SR/3lA/35p9fb/accmap_SR/3lA_3lA.txt /afs/cern.ch/user/c/cheidegg/www/heppy/2017-04-24_ewk80X_cutflows/acc/SR/3lA/35p9fb/accmap_SR/3lA_sig_TChiSlepSnu_TD_250_150.txt
#
### TChiSlepSnuTD 600/1
#python susy-interface/accmaker.py 3l 3lA $T $O -l $L -o SR --flags '$BLIND -A met lowMTT "mT2T_3l<100" -A met belowZ "mllL_3l<75"' -p "sig_TChiSlepSnu_TD_600_1" --noRFlags --cuts susy-ewkino/3l/cuts_cutflow.txt
#mv /afs/cern.ch/user/c/cheidegg/www/heppy/2017-04-24_ewk80X_cutflows/acc/SR/3lA/35p9fb/accmap_SR/3lA_3lA.txt /afs/cern.ch/user/c/cheidegg/www/heppy/2017-04-24_ewk80X_cutflows/acc/SR/3lA/35p9fb/accmap_SR/3lA_sig_TChiSlepSnu_TD_600_1.txt

### TChiSlepSnu05 500/350
#python susy-interface/accmaker.py 2lss 2lss $T $O -l $L -o SR --flags '$BLIND -A ISRjet lowMT "mtwmin(LepGood1_conePt,LepGood1_phi,LepGood2_conePt,LepGood2_phi,metPt,metPhi)<100" -A ISRjet highPt "pt_2_2lss(LepGood1_conePt,LepGood1_phi,LepGood2_conePt,LepGood2_phi)>100"' -p "sig_TChiSlepSnu_05_500_350" --noRFlags --cuts susy-ewkino/2lss_v1/cuts_cutflow.txt
#mv /afs/cern.ch/user/c/cheidegg/www/heppy/2017-04-24_ewk80X_cutflows/acc/SR/2lss/35p9fb/accmap_SR/2lss_2lss.txt /afs/cern.ch/user/c/cheidegg/www/heppy/2017-04-24_ewk80X_cutflows/acc/SR/2lss/35p9fb/accmap_SR/2lss_sig_TChiSlepSnu_05_500_350.txt

### TChiSlepSnuFD 500/350
#python susy-interface/accmaker.py 2lss 2lss $T $O -l $L -o SR --flags '$BLIND -A ISRjet lowMT "mtwmin(LepGood1_conePt,LepGood1_phi,LepGood2_conePt,LepGood2_phi,metPt,metPhi)<100" -A ISRjet highPt "pt_2_2lss(LepGood1_conePt,LepGood1_phi,LepGood2_conePt,LepGood2_phi)>100"' -p "sig_TChiSlepSnu_FD_500_350" --noRFlags --cuts susy-ewkino/2lss_v1/cuts_cutflow.txt
#mv /afs/cern.ch/user/c/cheidegg/www/heppy/2017-04-24_ewk80X_cutflows/acc/SR/2lss/35p9fb/accmap_SR/2lss_2lss.txt /afs/cern.ch/user/c/cheidegg/www/heppy/2017-04-24_ewk80X_cutflows/acc/SR/2lss/35p9fb/accmap_SR/2lss_sig_TChiSlepSnu_FD_500_350.txt

### TChiZZ 100/1
#python susy-interface/accmaker.py 3l 4lG $T $O -l $L -o SR --flags '$BLIND' -p "sig_TChiZZ4L_100_1" --noRFlags --cuts susy-ewkino/4l/cuts_cutflow.txt
#mv /afs/cern.ch/user/c/cheidegg/www/heppy/2017-04-24_ewk80X_cutflows/acc/SR/4lG/35p9fb/accmap_SR/4lG_4lG.txt /afs/cern.ch/user/c/cheidegg/www/heppy/2017-04-24_ewk80X_cutflows/acc/SR/4lG/35p9fb/accmap_SR/4lG_sig_TChiZZ4L_100_1.txt

### TChiZZ 800/1
python susy-interface/accmaker.py 3l 4lG $T $O -l $L -o SR --flags '$BLIND' -p "sig_TChiZZ4L_800_1" --noRFlags --cuts susy-ewkino/4l/cuts_cutflow.txt
mv /afs/cern.ch/user/c/cheidegg/www/heppy/2017-04-24_ewk80X_cutflows/acc/SR/4lG/35p9fb/accmap_SR/4lG_4lG.txt /afs/cern.ch/user/c/cheidegg/www/heppy/2017-04-24_ewk80X_cutflows/acc/SR/4lG/35p9fb/accmap_SR/4lG_sig_TChiZZ4L_800_1.txt






