#!/bin/bash

T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-25_ewkskims80X_M17_3l_bkg"
O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-03-01_ewk80X_sidebandMPFR" # Do NOT give a trailing /
FL=35.9
QUEUE="" #"-q all.q"

#python susy-interface/plotmaker.py 3l "3lA;3lB;3lC;3lD;3lF;4lG;4lH;4lI;4lJ;4lK" $T $O -l $FL --make bkgs --selPlots map_el_3l  -o el  --flags '-X blinding -I SRevent --perBin -A met lepprop "LepSel_isTight==0 && abs(LepSel_pdgId)==11"' $QUEUE -p fakes_bare --plot susy-ewkino/plots_mpfr.txt
#python susy-interface/plotmaker.py 3l "3lA;3lB;3lC;3lD;3lF;4lG;4lH;4lI;4lJ;4lK" $T $O -l $FL --make bkgs --selPlots map_mu_3l  -o mu  --flags '-X blinding -I SRevent --perBin -A met lepprop "LepSel_isTight==0 && abs(LepSel_pdgId)==13"' $QUEUE -p fakes_bare --plot susy-ewkino/plots_mpfr.txt
#python susy-interface/plotmaker.py 3l "3lA;3lB;3lC;3lD;3lF;4lG;4lH;4lI;4lJ;4lK" $T $O -l $FL --make bkgs --selPlots map_tau_3l -o tau --flags '-X blinding -I SRevent --perBin -A met lepprop "LepSel_isTight==0 && abs(LepSel_pdgId)==15"' $QUEUE -p fakes_bare --plot susy-ewkino/plots_mpfr.txt

python susy-interface/plotmaker.py 3l "3lE" $T $O -l $FL --make bkgs --selPlots map_el_3l  -o el  --flags '-X blinding -I SRevent --perBin -A met lepprop "LepSel_isTight==0 && abs(LepSel_pdgId)==11"' $QUEUE -p fakes_bare --plot susy-ewkino/plots_mpfr.txt
#python susy-interface/plotmaker.py 3l "3lE" $T $O -l $FL --make bkgs --selPlots map_mu_3l  -o mu  --flags '-X blinding -I SRevent --perBin -A met lepprop "LepSel_isTight==0 && abs(LepSel_pdgId)==13"' $QUEUE -p fakes_bare --plot susy-ewkino/plots_mpfr.txt
#python susy-interface/plotmaker.py 3l "3lE" $T $O -l $FL --make bkgs --selPlots map_tau_3l -o tau --flags '-X blinding -I SRevent --perBin -A met lepprop "LepSel_isTight==0 && abs(LepSel_pdgId)==15"' $QUEUE -p fakes_bare --plot susy-ewkino/plots_mpfr.txt
