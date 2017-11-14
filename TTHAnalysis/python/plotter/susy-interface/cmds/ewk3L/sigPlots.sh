#!/bin/bash

T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-25_ewkskims80X_M17_3l_bkg;/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-23_ewkskims80X_M17_3l_data;/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-18_ewktrees80X_M17_signals"
#O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-03-06_ewk80X_finalTest" # Do NOT give a trailing /
#O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-03-06_ewk80X_final" # Do NOT give a trailing /
O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-05-11_ewk80X_wzcorr" # Do NOT give a trailing /
O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-05-07_ewk80X_cwr4_poster" # Do NOT give a trailing /
L=35.9
FL=35.9
QUEUE="" #"-q all.q"
BLIND="-X blinding --cmsprel 'Preliminary'" #""

## SUSY exercise
## -----------------------------------------------------------------

#python susy-interface/plotmaker.py 3l 3lB $T $O -l $L  --make datasig --plots perCateg -o SRexerc --flags "--perBin $BLIND" -p "fakes_appldata;prompt_.*;rares_.*;convs;data;sig_TChiWH_150_24;sig_TChiWH_175_49;sig_TChiWH_175_1" $QUEUE
#python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make datasig --plots perCateg -o SRexerc --flags "--perBin $BLIND" -p "fakes_appldata;prompt_.*;rares_.*;convs;data;sig_TChiWZ_150_120;sig_TChiWZ_200_100;sig_TChiWZ_250_170" $QUEUE
#python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make datasig --plots perCateg -o SRexerc --flags "--perBin $BLIND" -p "fakes_appldata;prompt_.*;rares_.*;convs;data;sig_TChiWZ_400_150;sig_TChiWZ_400_175;sig_TChiWZ_450_50" $QUEUE
#python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make datasig --plots perCateg -o SRexerc --flags "--perBin $BLIND" -p "fakes_appldata;prompt_.*;rares_.*;convs;data;sig_TChiSlepSnu_500_400;sig_TChiSlepSnu_1000_1;sig_TChiSlepSnu_1050_700" $QUEUE

#python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make datasig --plots perCateg -o SRexerc --flags "--perBin $BLIND" -p "fakes_appldata;prompt_.*;rares_.*;convs;data;sig_TChiWZ_100_1;sig_TChiWZ_100_1_gen;sig_TChiWZ_100_1_avg" $QUEUE
#python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make datasig --selPlots "nJet30;.*_pt" -o SRexerc --flags "--perBin $BLIND" -p "fakes_appldata;prompt_.*;rares_.*;convs;data;sig_TChiWZ_100_1;sig_TChiWZ_100_1_gen" $QUEUE

#python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make datasig --plots perCateg -o SRexerc --flags "--perBin $BLIND" -p "fakes_appldata;prompt_.*;rares_.*;convs;data;sig_TChiWZ_175_115;sig_TChiWZ_200_100" $QUEUE
#python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make datasig --plots perCateg -o SRexerc --flags "--perBin $BLIND" -p "fakes_appldata;prompt_.*;rares_.*;convs;data;sig_TChiWZ_350_25;sig_TChiWZ_550_1" $QUEUE


## signal overlaid
## -----------------------------------------------------------------

#python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make datasig --plots perCateg -o SR   --flags "--perBin $BLIND" $QUEUE -p "fakes_appldata;promptsub;prompt_.*;rares_.*;convs;data;sig_TChiWZ_200_100;sig_TChiWZ_350_250"
#python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make datasig --plots perCateg -o SR   --flags "--perBin $BLIND" $QUEUE -p "fakes_appldata;promptsub;prompt_.*;rares_.*;convs;data;sig_TChiSlepSnu_FD_500_450;sig_TChiSlepSnu_FD_1100_1"
#python susy-interface/plotmaker.py 3l 3lB $T $O -l $L  --make datasig --plots perCateg -o SR   --flags "--perBin $BLIND" $QUEUE -p "fakes_appldata;promptsub;prompt_.*;rares_.*;convs;data;sig_TChiWH_150_24;sig_TChiWH_200_1"
#python susy-interface/plotmaker.py 3l 3lC $T $O -l $L  --make datasig --plots perCateg -o SR   --flags "--perBin $BLIND" $QUEUE -p "fakes_appldata;promptsub;prompt_.*;rares_.*;convs;data;sig_TChiSlepSnu_TE_400_350;sig_TChiSlepSnu_TE_1100_1"
#python susy-interface/plotmaker.py 3l 3lD $T $O -l $L  --make datasig --plots perCateg -o SR   --flags "--perBin $BLIND" $QUEUE -p "fakes_appldata;promptsub;prompt_.*;rares_.*;convs;data;sig_TChiSlepSnu_TD_250_150;sig_TChiSlepSnu_TD_600_1"
#python susy-interface/plotmaker.py 3l 3lE $T $O -l $L  --make datasig --plots perCateg -o SR   --flags "--perBin $BLIND" $QUEUE -p "fakes_appldata;promptsub;flips_appldata;prompt_.*;rares_.*;convs;data;sig_TChiSlepSnu_TD_250_150;sig_TChiSlepSnu_TD_600_1"
#python susy-interface/plotmaker.py 3l 3lF $T $O -l $L  --make datasig --plots perCateg -o SR   --flags "--perBin $BLIND" $QUEUE -p "fakes_appldata;promptsub;prompt_.*;rares_.*;convs;data;sig_TChiSlepSnu_TD_250_150;sig_TChiSlepSnu_TD_600_1"
#python susy-interface/plotmaker.py 3l 4lG $T $O -l $L  --make datasig --selPlots SR_4L -o SR   --flags "--perBin $BLIND --legendHeader 'SR G-K: 4l' -A alwaystrue BRcut 'BR >= 7'" -p "fakes_appldata;promptsub;prompt_.*;rares_.*;convs;data;sig_TChiZZ4L_100_1" --noRFlags $QUEUE

#python susy-interface/plotmaker.py 3l 4lG $T $O -l $L  --make datasig --plots perCateg -o SR   --flags "--perBin $BLIND" $QUEUE -p "fakes_appldata;promptsub;prompt_.*;rares_.*;convs;data;sig_TChiZZ4L_100_1;sig_TChiZZ4l_300_1"
#python susy-interface/plotmaker.py 3l 4lH $T $O -l $L  --make datasig --plots perCateg -o SR   --flags "--perBin $BLIND" $QUEUE -p "fakes_appldata;promptsub;prompt_.*;rares_.*;convs;data;sig_TChiZZ4L_100_1;sig_TChiZZ4l_300_1"
#python susy-interface/plotmaker.py 3l 4lI $T $O -l $L  --make datasig --plots perCateg -o SR   --flags "--perBin $BLIND" $QUEUE -p "fakes_appldata;promptsub;prompt_.*;rares_.*;convs;data;sig_TChiZZ4L_100_1;sig_TChiZZ4l_300_1"
#python susy-interface/plotmaker.py 3l 4lJ $T $O -l $L  --make datasig --plots perCateg -o SR   --flags "--perBin $BLIND" $QUEUE -p "fakes_appldata;promptsub;prompt_.*;rares_.*;convs;data;sig_TChiZZ4L_100_1;sig_TChiZZ4l_300_1"
#python susy-interface/plotmaker.py 3l 4lK $T $O -l $L  --make datasig --plots perCateg -o SR   --flags "--perBin $BLIND" $QUEUE -p "fakes_appldata;promptsub;prompt_.*;rares_.*;convs;data;sig_TChiZZ4L_100_1;sig_TChiZZ4l_300_1"
#python susy-interface/plotmaker.py 3l 4lJ $T $O -l $L  --make datasig --plots perCateg -o SR   --flags "--perBin $BLIND" $QUEUE -p "fakes_appldata;promptsub;prompt_.*;rares_.*;convs;data;sig_TChiZZ4L_800_1;sig_TChiZZ4L_1000_1"
#python susy-interface/plotmaker.py 3l 4lK $T $O -l $L  --make datasig --plots perCateg -o SR   --flags "--perBin $BLIND" $QUEUE -p "fakes_appldata;promptsub;prompt_.*;rares_.*;convs;data;sig_TChiZZ4L_400_1;sig_TChiZZ4L_425_1"


## kinematic variables
## -----------------------------------------------------------------
#python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make datasig --selPlots "met;mtW3l;mll3l" -o SR   --flags "--perBin $BLIND" $QUEUE -p "fakes_appldata;promptsub;prompt_.*;rares_.*;convs;data;sig_TChiSlepSnu_FD_500_450;sig_TChiSlepSnu_FD_1100_1"
python susy-interface/plotmaker.py 3l 3lC $T $O -l $L  --make datasig --selPlots mT2L3l -o SR   --flags "--perBin $BLIND" $QUEUE -p "fakes_appldata;promptsub;prompt_.*;rares_.*;convs;data;sig_TChiSlepSnu_TE_400_350;sig_TChiSlepSnu_TE_1100_1"
#python susy-interface/plotmaker.py 3l 3lF $T $O -l $L  --make datasig --selPlots met -o SR   --flags "--perBin $BLIND" $QUEUE -p "fakes_appldata;promptsub;prompt_.*;rares_.*;convs;data;sig_TChiSlepSnu_TD_250_150;sig_TChiSlepSnu_TD_600_1"
#python susy-interface/plotmaker.py 3l 4lG $T $O -l $L  --make datasig --selPlots met -o SR4lincl --flags "--perBin $BLIND --legendHeader 'SR G-K: 4l' -A alwaystrue BRcut 'BR >= 7'" --noRFlags $QUEUE -p "fakes_appldata;promptsub;prompt_.*;rares_.*;convs;data;sig_TChiZZ4L_100_1"

#python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make datasig --plots evt -o SR   --flags "--perBin $BLIND" $QUEUE -p "fakes_appldata;promptsub;prompt_.*;rares_.*;convs;data;sig_TChiSlepSnu_FD_500_450;sig_TChiSlepSnu_FD_1100_1"
#python susy-interface/plotmaker.py 3l 3lB $T $O -l $L  --make datasig --plots evt -o SR   --flags "--perBin $BLIND" $QUEUE -p "fakes_appldata;promptsub;prompt_.*;rares_.*;convs;data;sig_TChiWH_150_24;sig_TChiWH_200_1"
#python susy-interface/plotmaker.py 3l 3lC $T $O -l $L  --make datasig --plots evt -o SR   --flags "--perBin $BLIND" $QUEUE -p "fakes_appldata;promptsub;prompt_.*;rares_.*;convs;data;sig_TChiSlepSnu_TE_400_350;sig_TChiSlepSnu_TE_1100_1"
#python susy-interface/plotmaker.py 3l 3lD $T $O -l $L  --make datasig --plots evt -o SR   --flags "--perBin $BLIND" $QUEUE -p "fakes_appldata;promptsub;prompt_.*;rares_.*;convs;data;sig_TChiSlepSnu_TD_250_150;sig_TChiSlepSnu_TD_600_1"
#python susy-interface/plotmaker.py 3l 3lE $T $O -l $L  --make datasig --plots evt -o SR   --flags "--perBin $BLIND" $QUEUE -p "fakes_appldata;promptsub;flips_appldata;prompt_.*;rares_.*;convs;data;sig_TChiSlepSnu_TD_250_150;sig_TChiSlepSnu_TD_600_1"
#python susy-interface/plotmaker.py 3l 3lF $T $O -l $L  --make datasig --plots evt -o SR   --flags "--perBin $BLIND" $QUEUE -p "fakes_appldata;promptsub;prompt_.*;rares_.*;convs;data;sig_TChiSlepSnu_TD_250_150;sig_TChiSlepSnu_TD_600_1"
#python susy-interface/plotmaker.py 3l 4lG $T $O -l $L  --make datasig --plots evt -o SR4lincl --flags "--perBin $BLIND --legendHeader '4l: G-K' -A alwaystrue BRcut 'BR >= 7'" --noRFlags $QUEUE -p "fakes_appldata;promptsub;prompt_.*;rares_.*;convs;data;sig_TChiZZ4L_100_1"

#python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make datasig --plots lep -o SR   --flags "--perBin $BLIND" $QUEUE -p "fakes_appldata;promptsub;prompt_.*;rares_.*;convs;data;sig_TChiSlepSnu_FD_500_450;sig_TChiSlepSnu_FD_1100_1"
#python susy-interface/plotmaker.py 3l 3lB $T $O -l $L  --make datasig --plots lep -o SR   --flags "--perBin $BLIND" $QUEUE -p "fakes_appldata;promptsub;prompt_.*;rares_.*;convs;data;sig_TChiWH_150_24;sig_TChiWH_200_1"
#python susy-interface/plotmaker.py 3l 3lC $T $O -l $L  --make datasig --plots lep -o SR   --flags "--perBin $BLIND" $QUEUE -p "fakes_appldata;promptsub;prompt_.*;rares_.*;convs;data;sig_TChiSlepSnu_TE_400_350;sig_TChiSlepSnu_TE_1100_1"
#python susy-interface/plotmaker.py 3l 3lD $T $O -l $L  --make datasig --plots lep -o SR   --flags "--perBin $BLIND" $QUEUE -p "fakes_appldata;promptsub;prompt_.*;rares_.*;convs;data;sig_TChiSlepSnu_TD_250_150;sig_TChiSlepSnu_TD_600_1"
#python susy-interface/plotmaker.py 3l 3lE $T $O -l $L  --make datasig --plots lep -o SR   --flags "--perBin $BLIND" $QUEUE -p "fakes_appldata;promptsub;flips_appldata;prompt_.*;rares_.*;convs;data;sig_TChiSlepSnu_TD_250_150;sig_TChiSlepSnu_TD_600_1"
#python susy-interface/plotmaker.py 3l 3lF $T $O -l $L  --make datasig --plots lep -o SR   --flags "--perBin $BLIND" $QUEUE -p "fakes_appldata;promptsub;prompt_.*;rares_.*;convs;data;sig_TChiSlepSnu_TD_250_150;sig_TChiSlepSnu_TD_600_1"
#python susy-interface/plotmaker.py 3l 4lG $T $O -l $L  --make datasig --plots lep -o SR4lincl --flags "--perBin $BLIND --legendHeader '4l: G-K' -A alwaystrue BRcut 'BR >= 7'" --noRFlags $QUEUE -p "fakes_appldata;promptsub;prompt_.*;rares_.*;convs;data;sig_TChiZZ4L_100_1"


## new Score test
## -----------------------------------------------------------------

#O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-04-05_ewk80X_newScore" # Do NOT give a trailing /
#python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --selPlots SR_A -o SR   --flags "--perBin $BLIND" $QUEUE -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_500_150;sig_TChiWZ_200_100"
#python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --selPlots new_score -o SR   --flags "--perBin $BLIND --wide" $QUEUE -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_500_150"
#python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --selPlots new_score_200_100 -o SR   --flags "--perBin $BLIND --wide" $QUEUE -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100"
#python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --selPlots new_score50 -o SR   --flags "--perBin $BLIND --wide" $QUEUE -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_500_150"
#python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --selPlots new_score50_200_100 -o SR   --flags "--perBin $BLIND --wide" $QUEUE -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100"

#python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --selPlots newScore_200_100 -o SR   --flags "--perBin $BLIND --wide" $QUEUE -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_200_100"
#python susy-interface/plotmaker.py 3l 3lA $T $O -l $L  --make mix --selPlots newScore_500_150 -o SR   --flags "--perBin $BLIND --wide" $QUEUE -p "fakes_matched;prompt_.*;rares_.*;convs;sig_TChiWZ_500_150"


