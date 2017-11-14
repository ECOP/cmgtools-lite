#!/bin/bash

#T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-04_signalMerger/"
#O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-04_signalSplitter"
#T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-01-24_treeMerger/"
#O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-01-24_signalSplitter"
O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-03-22_signalSplitter"
T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-03-22_signalMerger"
QUEUE="-q long.q --email 'constantin.heidegger@cern.ch'"

#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --accept T6ttWW
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --accept T6ttHZ --flags '--branch1 GenSusyMStop --branch2 GenSusyMStop2'
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --accept T6bbllslepton_mSbottom_1000To1500_mLSP_120To1450
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --accept T6bbllslepton_mSbottom_400To950_mLSP_120To140
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --accept T6bbllslepton_mSbottom_400To575_mLSP_150To550
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --accept T6bbllslepton_mSbottom_600To775_mLSP_150To725
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --accept T6bbllslepton_mSbottom_800To950_mLSP_150To900 

#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --accept TChiSlepSnux0p5  --flags '--branch1 GenSusyMChargino --branch2 GenSusyMNeutralino'
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --accept TChiSlepSnux0p05 --flags '--branch1 GenSusyMChargino --branch2 GenSusyMNeutralino'
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --accept TChiSlepSnux0p95 --flags '--branch1 GenSusyMChargino --branch2 GenSusyMNeutralino'
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --accept TChiSlepSnux0p5ext  --flags '--branch1 GenSusyMChargino --branch2 GenSusyMNeutralino'
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --accept TChiSlepSnux0p05ext --flags '--branch1 GenSusyMChargino --branch2 GenSusyMNeutralino'
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --accept TChiSlepSnuTEx0p5  --flags '--branch1 GenSusyMChargino --branch2 GenSusyMNeutralino'
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --accept TChiSlepSnuTEx0p05 --flags '--branch1 GenSusyMChargino --branch2 GenSusyMNeutralino'
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --accept TChiSlepSnuTEx0p95 --flags '--branch1 GenSusyMChargino --branch2 GenSusyMNeutralino'
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --accept TChiWZ           --flags '--branch1 GenSusyMChargino --branch2 GenSusyMNeutralino'
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --accept TChiWH           --flags '--branch1 GenSusyMChargino --branch2 GenSusyMNeutralino'

#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --accept TChiStauStaux0p5 --exclude ext   --flags '--branch1 GenSusyMChargino --branch2 GenSusyMNeutralino'
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --accept TChiStauStaux0p5ext --flags '--branch1 GenSusyMChargino --branch2 GenSusyMNeutralino'

python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --accept TChiZZ4L            --flags '--branch1 GenSusyMNeutralino2 --branch2 GenSusyMNeutralino'
python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --accept TChiZZ2L            --flags '--branch1 GenSusyMNeutralino2 --branch2 GenSusyMNeutralino'
python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --accept TChiHZ              --flags '--branch1 GenSusyMNeutralino2 --branch2 GenSusyMNeutralino'
python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --accept TChiHH              --flags '--branch1 GenSusyMNeutralino2 --branch2 GenSusyMNeutralino'

#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --accept T1tttt  
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --accept T5qqqqVV
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --accept T6ttWW --flags '--branch1 GenSusyMSbottom --branch2 GenSusyMNeutralino'
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --accept T6ttHZ --flags '--branch1 GenSusyMStop2 --branch2 GenSusyMStop'





## Winos
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --hist --gen --minmass 120 --maxmass  700 --step 25 --accept TChiWH
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --hist --gen --minmass 100 --maxmass  700 --step 25 --accept TChiWZ
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --hist --gen --minmass 100 --maxmass 1200 --step 50 --accept TChiSlepSnux0p5 --exclude ext
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --hist --gen --minmass 100 --maxmass 1200 --step 50 --accept TChiSlepSnuTEx0p5
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --hist --gen --minmass 100 --maxmass  700 --step 25 --accept SMS_TChiStauStaux0p5ext

## Higgsinos
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --hist --gen --minmass 100 --maxmass 1000 --step 25 --accept TChiZZ4L
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --hist --gen --minmass 100 --maxmass 1000 --step 25 --accept TChiZZ2
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --hist --gen --minmass 125 --maxmass 1000 --step 25 --accept TChiHZ
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --hist --gen --minmass 125 --maxmass 1000 --step 25 --accept TChiHH

## Strong Signals
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --accept SMS_T5qqqqVV --flags '--hist'
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --accept SMS_T1tttt   --flags '--hist'
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --accept SMS_T6bbllslepton_mSbottom_400To575_mLSP_150To550 --flags '--hist'
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --accept SMS_T6bbllslepton_mSbottom_400To950_mLSP_120To140 --flags '--hist'
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --accept SMS_T6bbllslepton_mSbottom_600To775_mLSP_150To725 --flags '--hist'
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --accept SMS_T6bbllslepton_mSbottom_800To950_mLSP_150To900 --flags '--hist'
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --accept SMS_T6bbllslepton_mSbottom_1000To1500_mLSP_120To1450 --flags '--hist'

########################## Gen Info ###############################

## Winos
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE  --accept TChiWH --gen
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --flags '--hist' --gen --minmass 125 --maxmass  700 --step 25 --accept TChiWH
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --hist --gen --minmass 100 --maxmass  700 --step 25 --accept TChiWZ
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --hist --gen --minmass 100 --maxmass 1200 --step 50 --accept TChiSlepSnux0p5 --exclude ext
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --hist --gen --minmass 100 --maxmass 1200 --step 50 --accept TChiSlepSnuTEx0p5
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --gen --accept SMS_TChiStauStaux0p5ext

## Higgsinos
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --hist --gen --minmass 100 --maxmass 1000 --step 25 --accept TChiZZ4L
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --hist --gen --minmass 100 --maxmass 1000 --step 25 --accept TChiZZ2
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --hist --gen --minmass 125 --maxmass 1000 --step 25 --accept TChiHZ
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --hist --gen --minmass 125 --maxmass 1000 --step 25 --accept TChiHH

## Strong Signals
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --gen --accept SMS_T5qqqqVV --pdgId1 1000021 --pretend
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --gen --accept SMS_T1tttt --pdgId1 1000021 --pretend
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --gen --accept SMS_T6bbllslepton_mSbottom_400To950_mLSP_120To140 --pdgId1 1000005 --pdgId2 1000023
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --gen --accept SMS_T6bbllslepton_mSbottom_1000To1500_mLSP_120To1450 --pdgId1 1000005 --pdgId2 1000023
#python susy-interface/splitmaker.py 3l 3lA $T $O $QUEUE --gen --accept SMS_T6bbllslepton_mSbottom_400To575_mLSP_150To550 --pdgId1 1000005 --pdgId2 1000023


