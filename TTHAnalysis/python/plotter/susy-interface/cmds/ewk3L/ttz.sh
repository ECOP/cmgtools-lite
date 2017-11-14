#!/bin/bash

T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-25_ewkskims80X_M17_ttz_bkg;/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-25_ewkskims80X_M17_ttz_data"
O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-03-19_ewk80X_ttZ" # Do NOT give a trailing /
LGH=16
LBF=20
QUEUE="" #"-q all.q"
BLIND="-X blinding" #""


## runs B-F
#python susy-interface/plotmaker.py 3l "3lA" $T $O -l $LBF --make data --selPlots "jet1_csv;jet2_csv;jet3_csv;jet4_csv;jet5_csv;jet6_csv;jet7_csv" -o SRbf --noRFlags --cuts susy-ewkino/3l/cuts_ttz.txt --flags "--perBin $BLIND -X onlyGH" $QUEUE
python susy-interface/plotmaker.py 3l "3lA" $T $O -l $LBF --make data --selPlots "nJet30" -o SRbf0b --noRFlags --cuts susy-ewkino/3l/cuts_ttz.txt --flags "--perBin $BLIND -X onlyGH -A nJets nBJets 'nBJetMedium25_Mini==0'" $QUEUE
python susy-interface/plotmaker.py 3l "3lA" $T $O -l $LBF --make data --selPlots "nJet30" -o SRbf1b --noRFlags --cuts susy-ewkino/3l/cuts_ttz.txt --flags "--perBin $BLIND -X onlyGH -A nJets nBJets 'nBJetMedium25_Mini==1'" $QUEUE
python susy-interface/plotmaker.py 3l "3lA" $T $O -l $LBF --make data --selPlots "nJet30" -o SRbf2b --noRFlags --cuts susy-ewkino/3l/cuts_ttz.txt --flags "--perBin $BLIND -X onlyGH -A nJets nBJets 'nBJetMedium25_Mini>=2'" $QUEUE

## runs G-H
#python susy-interface/plotmaker.py 3l "3lA" $T $O -l $LGH --make data --selPlots "jet1_csv;jet2_csv;jet3_csv;jet4_csv;jet5_csv;jet6_csv;jet7_csv" -o SRgh --noRFlags --cuts susy-ewkino/3l/cuts_ttz.txt --flags "--perBin $BLIND -X onlyBF" $QUEUE
python susy-interface/plotmaker.py 3l "3lA" $T $O -l $LGH --make data --selPlots "nJet30" -o SRgh0b --noRFlags --cuts susy-ewkino/3l/cuts_ttz.txt --flags "--perBin $BLIND -X onlyBF -A nJets nBJets 'nBJetMedium25_Mini==0'" $QUEUE
python susy-interface/plotmaker.py 3l "3lA" $T $O -l $LGH --make data --selPlots "nJet30" -o SRgh1b --noRFlags --cuts susy-ewkino/3l/cuts_ttz.txt --flags "--perBin $BLIND -X onlyBF -A nJets nBJets 'nBJetMedium25_Mini==1'" $QUEUE
python susy-interface/plotmaker.py 3l "3lA" $T $O -l $LGH --make data --selPlots "nJet30" -o SRgh2b --noRFlags --cuts susy-ewkino/3l/cuts_ttz.txt --flags "--perBin $BLIND -X onlyBF -A nJets nBJets 'nBJetMedium25_Mini>=2'" $QUEUE
