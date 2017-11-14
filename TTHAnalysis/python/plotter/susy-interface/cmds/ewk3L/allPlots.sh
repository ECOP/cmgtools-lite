#!/bin/bash

#T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-01-27_ewkskims80X_M17_MERGED"
#O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-02-24_ewk80X_upload" # Do NOT give a trailing /
#O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-02-01_ewk80X_preapproval" # Do NOT give a trailing /
T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-25_ewkskims80X_M17_3l_bkg;/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-23_ewkskims80X_M17_3l_data;/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-18_ewktrees80X_M17_signals"
O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-03-06_ewk80X_finalTest" # Do NOT give a trailing /
#O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-03-06_ewk80X_final" # Do NOT give a trailing /
L=35.9
FL=35.9
QUEUE="" #"-q all.q"
BLIND="-X blinding" #""

## WITH LOGARITHMIC SCALE
## =================================================================

## signal regions and baseline regions plots
## -----------------------------------------------------------------

## data fakes
#python susy-interface/plotmaker.py 3l "3lA;3lB;3lC;3lD;3lF;4lG;4lH;4lI;4lJ;4lK" $T $O -l $L  --make data --plots perCateg -o SR   --flags "--perBin $BLIND" $QUEUE
#python susy-interface/plotmaker.py 3l "3lE"                                     $T $O -l $L  --make data --plots perCateg -o SR   --flags "--perBin $BLIND" $QUEUE -p "fakes_appldata;promptsub;flips_appldata;prompt.*;rares_.*;convs;data"
python susy-interface/plotmaker.py 3l "3lA;3lB;3lC;3lD;3lF;4lG;4lH;4lI;4lJ;4lK" $T $O -l $L  --make data --plots evt      -o SR   --flags "--perBin $BLIND" $QUEUE
python susy-interface/plotmaker.py 3l "3lE"                                     $T $O -l $L  --make data --plots evt      -o SR   --flags "--perBin $BLIND" $QUEUE -p "fakes_appldata;promptsub;flips_appldata;prompt.*;rares_.*;convs;data"
#python susy-interface/plotmaker.py 3l "3lA;3lB;3lC;3lD;3lF;4lG;4lH;4lI;4lJ;4lK" $T $O -l $L  --make data --plots lep      -o SR   --flags "--perBin $BLIND" $QUEUE
#python susy-interface/plotmaker.py 3l "3lE"                                     $T $O -l $L  --make data --plots lep      -o SR   --flags "--perBin $BLIND" $QUEUE -p "fakes_appldata;promptsub;flips_appldata;prompt.*;rares_.*;convs;data"
#
### data fakes inclusive 4l
#python susy-interface/plotmaker.py 3l "4lG" $T $O -l $L  --make data --selPlots SR_4L -o SR   --flags "--perBin $BLIND --legendHeader '4l: G-K' -A alwaystrue BRcut 'BR >= 7'" --noRFlags $QUEUE
python susy-interface/plotmaker.py 3l "4lG" $T $O -l $L  --make data --plots evt -o SR4lincl --flags "--perBin $BLIND --legendHeader '4l: G-K' -A alwaystrue BRcut 'BR >= 7'" --noRFlags $QUEUE
#python susy-interface/plotmaker.py 3l "4lG" $T $O -l $L  --make data --plots lep -o SR4lincl --flags "--perBin $BLIND --legendHeader '4l: G-K' -A alwaystrue BRcut 'BR >= 7'" --noRFlags $QUEUE
#
#
#
### MC fakes
#python susy-interface/plotmaker.py 3l "3lA;3lB;3lC;3lD;3lF;4lG;4lH;4lI;4lJ;4lK" $T $O -l $L  --make data --plots perCateg -o SRmc --flags "--perBin $BLIND" $QUEUE -p "fakes_matched;prompt_.*;rares_.*;convs;data" 
#python susy-interface/plotmaker.py 3l "3lA;3lB;3lC;3lD;3lF;4lG;4lH;4lI;4lJ;4lK" $T $O -l $L  --make data --plots evt      -o SRmc --flags "--perBin $BLIND" $QUEUE -p "fakes_matched;prompt_.*;rares_.*;convs;data" 
#python susy-interface/plotmaker.py 3l "3lA;3lB;3lC;3lD;3lF;4lG;4lH;4lI;4lJ;4lK" $T $O -l $L  --make data --plots lep      -o SRmc --flags "--perBin $BLIND" $QUEUE -p "fakes_matched;prompt_.*;rares_.*;convs;data" 
#python susy-interface/plotmaker.py 3l "3lE"                                     $T $O -l $L  --make data --plots perCateg -o SRmc --flags "--perBin $BLIND" $QUEUE -p "fakes_matched;flips_appldata;prompt_.*;rares_.*;convs;data" 
#python susy-interface/plotmaker.py 3l "3lE"                                     $T $O -l $L  --make data --plots evt      -o SRmc --flags "--perBin $BLIND" $QUEUE -p "fakes_matched;flips_appldata;prompt_.*;rares_.*;convs;data" 
#python susy-interface/plotmaker.py 3l "3lE"                                     $T $O -l $L  --make data --plots lep      -o SRmc --flags "--perBin $BLIND" $QUEUE -p "fakes_matched;flips_appldata;prompt_.*;rares_.*;convs;data" 
#
### MC fakes splitted
#python susy-interface/plotmaker.py 3l "3lA;3lB;3lC;3lD;3lF;4lG;4lH;4lI;4lJ;4lK" $T $O -l $L  --make data --plots perCateg -o SRmcproc --flags "--perBin $BLIND" $QUEUE -p "fakes_process_.*;prompt_.*;rares_.*;convs;data" 
#python susy-interface/plotmaker.py 3l "3lA;3lB;3lC;3lD;3lF;4lG;4lH;4lI;4lJ;4lK" $T $O -l $L  --make data --plots evt      -o SRmcproc --flags "--perBin $BLIND" $QUEUE -p "fakes_process_.*;prompt_.*;rares_.*;convs;data" 
#python susy-interface/plotmaker.py 3l "3lA;3lB;3lC;3lD;3lF;4lG;4lH;4lI;4lJ;4lK" $T $O -l $L  --make data --plots lep      -o SRmcproc --flags "--perBin $BLIND" $QUEUE -p "fakes_process_.*;prompt_.*;rares_.*;convs;data" 
#python susy-interface/plotmaker.py 3l "3lE"                                     $T $O -l $L  --make data --plots perCateg -o SRmcproc --flags "--perBin $BLIND" $QUEUE -p "fakes_process_.*;flips_appldata;prompt_.*;rares_.*;convs;data" 
#python susy-interface/plotmaker.py 3l "3lE"                                     $T $O -l $L  --make data --plots evt      -o SRmcproc --flags "--perBin $BLIND" $QUEUE -p "fakes_process_.*;flips_appldata;prompt_.*;rares_.*;convs;data" 
#python susy-interface/plotmaker.py 3l "3lE"                                     $T $O -l $L  --make data --plots lep      -o SRmcproc --flags "--perBin $BLIND" $QUEUE -p "fakes_process_.*;flips_appldata;prompt_.*;rares_.*;convs;data" 
#
### sideband
### -----------------------------------------------------------------
#
### MC fakes
#python susy-interface/plotmaker.py 3l "3lE" $T $O -l $FL --make data --plots perCateg -o SBfakes --flags '-X blinding -I SRevent --perBin' $QUEUE -p "fakes_bare;bare_.*;data" --plot susy-ewkino/3l/plots_ewkino_sideband.txt  
#python susy-interface/plotmaker.py 3l "3lE" $T $O -l $FL --make data --plots evt      -o SBfakes --flags '-X blinding -I SRevent --perBin' $QUEUE -p "fakes_bare;bare_.*;data" --plot susy-ewkino/3l/plots_ewkino_sideband.txt 
#python susy-interface/plotmaker.py 3l "3lE" $T $O -l $FL --make data --plots lep      -o SBfakes --flags '-X blinding -I SRevent --perBin' $QUEUE -p "fakes_bare;bare_.*;data" --plot susy-ewkino/3l/plots_ewkino_sideband.txt 
#python susy-interface/plotmaker.py 3l "3lA;3lB;3lC;3lD;3lE;3lF;4lG;4lH;4lI;4lJ;4lK" $T $O -l $FL --make data --plots perCateg -o SBfakes --flags '-X blinding -I SRevent --perBin' $QUEUE -p "fakes_bare;bare_.*;data" --plot susy-ewkino/3l/plots_ewkino_sideband.txt  
#python susy-interface/plotmaker.py 3l "3lA;3lB;3lC;3lD;3lE;3lF;4lG;4lH;4lI;4lJ;4lK" $T $O -l $FL --make data --plots evt      -o SBfakes --flags '-X blinding -I SRevent --perBin' $QUEUE -p "fakes_bare;bare_.*;data" --plot susy-ewkino/3l/plots_ewkino_sideband.txt 
#python susy-interface/plotmaker.py 3l "3lA;3lB;3lC;3lD;3lE;3lF;4lG;4lH;4lI;4lJ;4lK" $T $O -l $FL --make data --plots lep      -o SBfakes --flags '-X blinding -I SRevent --perBin' $QUEUE -p "fakes_bare;bare_.*;data" --plot susy-ewkino/3l/plots_ewkino_sideband.txt 
#
#
#
### NO LOGARITHMIC SCALE
### =================================================================
#O="${O}_noLog"
#
### signal regions and baseline regions plots
### -----------------------------------------------------------------
#
### data fakes
#python susy-interface/plotmaker.py 3l "3lA;3lB;3lC;3lD;3lF;4lG;4lH;4lI;4lJ;4lK" $T $O -l $L --make data --plots perCateg -o SR   --flags "--perBin $BLIND" $QUEUE --plot susy-ewkino/3l/plots_ewkino_nolog.txt
#python susy-interface/plotmaker.py 3l "3lA;3lB;3lC;3lD;3lF;4lG;4lH;4lI;4lJ;4lK" $T $O -l $L --make data --plots evt      -o SR   --flags "--perBin $BLIND" $QUEUE --plot susy-ewkino/3l/plots_ewkino_nolog.txt
#python susy-interface/plotmaker.py 3l "3lA;3lB;3lC;3lD;3lF;4lG;4lH;4lI;4lJ;4lK" $T $O -l $L --make data --plots lep      -o SR   --flags "--perBin $BLIND" $QUEUE --plot susy-ewkino/3l/plots_ewkino_nolog.txt
#python susy-interface/plotmaker.py 3l "3lE"                                     $T $O -l $L  --make data --plots perCateg -o SR   --flags "--perBin $BLIND" $QUEUE --plot susy-ewkino/3l/plots_ewkino_nolog.txt -p "fakes_appldata;promptsub;flips_appldata;prompt_.*;rares_.*;convs;data"
#python susy-interface/plotmaker.py 3l "3lE"                                     $T $O -l $L  --make data --plots evt      -o SR   --flags "--perBin $BLIND" $QUEUE --plot susy-ewkino/3l/plots_ewkino_nolog.txt -p "fakes_appldata;promptsub;flips_appldata;prompt_.*;rares_.*;convs;data"
#python susy-interface/plotmaker.py 3l "3lE"                                     $T $O -l $L  --make data --plots lep      -o SR   --flags "--perBin $BLIND" $QUEUE --plot susy-ewkino/3l/plots_ewkino_nolog.txt -p "fakes_appldata;promptsub;flips_appldata;prompt_.*;rares_.*;convs;data"
#
### MC fakes
#python susy-interface/plotmaker.py 3l "3lA;3lB;3lC;3lD;3lF;4lG;4lH;4lI;4lJ;4lK" $T $O -l $L --make data --plots perCateg -o SRmc --flags "--perBin $BLIND" $QUEUE -p "fakes_matched;prompt_.*;rares_.*;convs;data" --plot susy-ewkino/3l/plots_ewkino_nolog.txt
#python susy-interface/plotmaker.py 3l "3lA;3lB;3lC;3lD;3lF;4lG;4lH;4lI;4lJ;4lK" $T $O -l $L --make data --plots evt      -o SRmc --flags "--perBin $BLIND" $QUEUE -p "fakes_matched;prompt_.*;rares_.*;convs;data" --plot susy-ewkino/3l/plots_ewkino_nolog.txt
#python susy-interface/plotmaker.py 3l "3lA;3lB;3lC;3lD;3lF;4lG;4lH;4lI;4lJ;4lK" $T $O -l $L --make data --plots lep      -o SRmc --flags "--perBin $BLIND" $QUEUE -p "fakes_matched;prompt_.*;rares_.*;convs;data" --plot susy-ewkino/3l/plots_ewkino_nolog.txt
#python susy-interface/plotmaker.py 3l "3lE"                                     $T $O -l $L --make data --plots perCateg -o SRmc --flags "--perBin $BLIND" $QUEUE -p "fakes_matched;flips_appldata;prompt_.*;rares_.*;convs;data" --plot susy-ewkino/3l/plots_ewkino_nolog.txt
#python susy-interface/plotmaker.py 3l "3lE"                                     $T $O -l $L --make data --plots evt      -o SRmc --flags "--perBin $BLIND" $QUEUE -p "fakes_matched;flips_appldata;prompt_.*;rares_.*;convs;data" --plot susy-ewkino/3l/plots_ewkino_nolog.txt
#python susy-interface/plotmaker.py 3l "3lE"                                     $T $O -l $L --make data --plots lep      -o SRmc --flags "--perBin $BLIND" $QUEUE -p "fakes_matched;flips_appldata;prompt_.*;rares_.*;convs;data" --plot susy-ewkino/3l/plots_ewkino_nolog.txt
#
### MC fakes processsplit
#python susy-interface/plotmaker.py 3l "3lA;3lB;3lC;3lD;3lF;4lG;4lH;4lI;4lJ;4lK" $T $O -l $L --make data --plots perCateg -o SRmcproc --flags "--perBin $BLIND" $QUEUE -p "fakes_process_.*;prompt_.*;rares_.*;convs;data" --plot susy-ewkino/3l/plots_ewkino_nolog.txt
#python susy-interface/plotmaker.py 3l "3lA;3lB;3lC;3lD;3lF;4lG;4lH;4lI;4lJ;4lK" $T $O -l $L --make data --plots evt      -o SRmcproc --flags "--perBin $BLIND" $QUEUE -p "fakes_process_.*;prompt_.*;rares_.*;convs;data" --plot susy-ewkino/3l/plots_ewkino_nolog.txt
#python susy-interface/plotmaker.py 3l "3lA;3lB;3lC;3lD;3lF;4lG;4lH;4lI;4lJ;4lK" $T $O -l $L --make data --plots lep      -o SRmcproc --flags "--perBin $BLIND" $QUEUE -p "fakes_process_.*;prompt_.*;rares_.*;convs;data" --plot susy-ewkino/3l/plots_ewkino_nolog.txt
#python susy-interface/plotmaker.py 3l "3lE"                                     $T $O -l $L --make data --plots perCateg -o SRmcproc --flags "--perBin $BLIND" $QUEUE -p "fakes_process_.*;flips_appldata;prompt_.*;rares_.*;convs;data" --plot susy-ewkino/3l/plots_ewkino_nolog.txt
#python susy-interface/plotmaker.py 3l "3lE"                                     $T $O -l $L --make data --plots evt      -o SRmcproc --flags "--perBin $BLIND" $QUEUE -p "fakes_process_.*;flips_appldata;prompt_.*;rares_.*;convs;data" --plot susy-ewkino/3l/plots_ewkino_nolog.txt
#python susy-interface/plotmaker.py 3l "3lE"                                     $T $O -l $L --make data --plots lep      -o SRmcproc --flags "--perBin $BLIND" $QUEUE -p "fakes_process_.*;flips_appldata;prompt_.*;rares_.*;convs;data" --plot susy-ewkino/3l/plots_ewkino_nolog.txt
#
### sideband
### -----------------------------------------------------------------
#
### MC fakes
#python susy-interface/plotmaker.py 3l "3lE" $T $O -l $FL --make data --plots perCateg -o SBfakes --flags '-X blinding -I SRevent --perBin' $QUEUE -p "fakes_bare;bare_.*;data" --plot susy-ewkino/3l/plots_ewkino_nolog.txt 
#python susy-interface/plotmaker.py 3l "3lE" $T $O -l $FL --make data --plots evt      -o SBfakes --flags '-X blinding -I SRevent --perBin' $QUEUE -p "fakes_bare;bare_.*;data" --plot susy-ewkino/3l/plots_ewkino_nolog.txt
#python susy-interface/plotmaker.py 3l "3lE" $T $O -l $FL --make data --plots lep      -o SBfakes --flags '-X blinding -I SRevent --perBin' $QUEUE -p "fakes_bare;bare_.*;data" --plot susy-ewkino/3l/plots_ewkino_nolog.txt
#python susy-interface/plotmaker.py 3l "3lA;3lB;3lC;3lD;3lE;3lF;4lG;4lH;4lI;4lJ;4lK" $T $O -l $FL --make data --plots perCateg -o SBfakes --flags '-X blinding -I SRevent --perBin' $QUEUE -p "fakes_bare;bare_.*;data" --plot susy-ewkino/3l/plots_ewkino_nolog.txt 
#python susy-interface/plotmaker.py 3l "3lA;3lB;3lC;3lD;3lE;3lF;4lG;4lH;4lI;4lJ;4lK" $T $O -l $FL --make data --plots evt      -o SBfakes --flags '-X blinding -I SRevent --perBin' $QUEUE -p "fakes_bare;bare_.*;data" --plot susy-ewkino/3l/plots_ewkino_nolog.txt
#python susy-interface/plotmaker.py 3l "3lA;3lB;3lC;3lD;3lE;3lF;4lG;4lH;4lI;4lJ;4lK" $T $O -l $FL --make data --plots lep      -o SBfakes --flags '-X blinding -I SRevent --perBin' $QUEUE -p "fakes_bare;bare_.*;data" --plot susy-ewkino/3l/plots_ewkino_nolog.txt
#
#
#
### CONFERENCE MODE
### =================================================================
#
##O="/afs/cern.ch/user/c/cheidegg/www/heppy/2016-12-17_confTest" # Do NOT give a trailing /
#
### log scale
##python susy-interface/plotmaker.py 3l "3lA;3lB;3lC;3lD;3lE;3lF" $T $O -l $L  --make data --plots perCateg -o SR   --flags "--perBin $BLIND" -p "prompt_WZ;fakes_appldata;promptsub;conf;data" $QUEUE
##python susy-interface/plotmaker.py 3l "4lG;4lH;4lI;4lJ;4lK"     $T $O -l $L  --make data --plots perCateg -o SR   --flags "--perBin $BLIND" -p "conf_ZZH;fakes_appldata;promptsub;conf;data"  $QUEUE
#
### no log scale
##python susy-interface/plotmaker.py 3l "3lA;3lB;3lC;3lD;3lE;3lF" $T $O -l $L --make data --plots perCateg -o SR   --flags "--perBin $BLIND" -p "prompt_WZ;fakes_appldata;promptsub;conf;data" $QUEUE --plot susy-ewkino/3l/plots_ewkino_nolog.txt
##python susy-interface/plotmaker.py 3l "4lG;4lH;4lI;4lJ;4lK"     $T $O -l $L --make data --plots perCateg -o SR   --flags "--perBin $BLIND" -p "conf_ZZH;fakes_appldata;promptsub;conf;data"  $QUEUE --plot susy-ewkino/3l/plots_ewkino_nolog.txt
#
#
### DEBUGGING
### =================================================================
##python susy-interface/plotmaker.py 3l "3lA;3lB;3lC;3lD;3lF;4lG;4lH;4lI;4lJ;4lK" $T $O -l $L --make data --selPlots "mu_pt_badClone;mu_pt_badNotClone" -o SR --flags "--perBin $BLIND --emptyStack" $QUEUE -p "dataNew"
##python susy-interface/plotmaker.py 3l "4lG" $T $O -l $L  --make data --plots perCateg -o SR   --flags "--perBin $BLIND" $QUEUE
