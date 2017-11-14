#!/bin/bash

S="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-18_ewktrees80X_M17_signals"
T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-25_ewkskims80X_M17_3l_bkg;/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-23_ewkskims80X_M17_3l_data;$S"
TT="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-25_ewkskims80X_M17_crwz_bkg;/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-23_ewkskims80X_M17_crwz_data;$S"
L=35.9
FL=35.9
QUEUE="-q all.q" #"-q all.q"
BLIND="-X blinding" #""

## DEFINE WHAT YOU WANNA DO
## =================================================================
## =================================================================
#DO="wz4p7"
#DO="wz4p7perBin"
#DO="cr"
DO="crfit"


## full binning with HT and MT bin, wz4p7, perBin
## -----------------------------------------------------------------

if [ $DO = "wz4p7" ]; then
    O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-04-08_ewk80X_wzCorrCards_wz4p7"
    #python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --bkgOnly --redoBkg --flags "$BLIND" --mca susy-ewkino/3lcorr/mca_forScan.txt $QUEUE
    python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3lcorr/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 350 --m2 250
    python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3lcorr/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 300 --m2 200
    python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3lcorr/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 250 --m2 150
    python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3lcorr/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --m1 200 --m2 100
fi



## full binning with HT and MT bin, wz4p7, perBin
## -----------------------------------------------------------------

if [ $DO = "wz4p7perBin" ]; then
    O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-04-08_ewk80X_wzCorrCards_wz4p7perBin"
    #python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --bkgOnly --redoBkg --flags "$BLIND" --mca susy-ewkino/3lcorr/mca_forScan.txt --perBin $QUEUE
    python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3lcorr/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --perBin --m1 350 --m2 250
    python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3lcorr/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --perBin --m1 300 --m2 200
    python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3lcorr/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --perBin --m1 250 --m2 150
    python susy-interface/scanmaker.py 3lcorr 3lW $T $O -l $L -o SR --models TChiWZ --mca susy-ewkino/3lcorr/mca_forScan.txt --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpM_FL_data_comb --mpfr 3lA' --check --perBin --m1 200 --m2 100
fi


## WZ CR
## -----------------------------------------------------------------

if [ $DO = "cr" ]; then
    O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-04-08_ewk80X_wzCorrCards_wz4p7"
    #python susy-interface/scanmaker.py crwzcorr crwz $TT $O -l $L -o SR --models TChiWZ --bkgOnly --redoBkg --flags "$BLIND" $QUEUE
    python susy-interface/scanmaker.py crwzcorr crwz $TT $O -l $L -o SR --models TChiWZ --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpV_FL_data_comb --mpfr 2lss' --check --m1 350 --m2 250
    python susy-interface/scanmaker.py crwzcorr crwz $TT $O -l $L -o SR --models TChiWZ --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpV_FL_data_comb --mpfr 2lss' --check --m1 300 --m2 200
    python susy-interface/scanmaker.py crwzcorr crwz $TT $O -l $L -o SR --models TChiWZ --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpV_FL_data_comb --mpfr 2lss' --check --m1 250 --m2 150
    python susy-interface/scanmaker.py crwzcorr crwz $TT $O -l $L -o SR --models TChiWZ --flags "$BLIND" --sigOnly $QUEUE --postfix '--postfix-pred fakes_appldata*=fixFakePredictionForZeroEvts --frFile $CMSSW_BASE/src/CMGTools/TTHAnalysis/data/fakerate/fakerate_light_EWKino_Moriond2017_36p5fb.root --frMap FR_susy_wpV_FL_data_comb --mpfr 2lss' --check --m1 200 --m2 100
fi


## WZ CR ML FIT
## -----------------------------------------------------------------

if [ $DO = "crfit" ]; then

    T1="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-25_ewkskims80X_M17_crwz_bkg"
    T2="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-23_ewkskims80X_M17_crwz_data"
    O="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-06-16_ewk80X_wzcrfit_bugFix"
    python makeShapeCardsSusy.py susy-ewkino/crwzcorr/mca_includesM17.txt susy-ewkino/crwzcorr/cuts_crwz.txt 1 '1,0.5,1.5' susy-ewkino/crwzcorr/systs_crwz.txt -P $T1 -P $T2 --Fs {P}/leptonJetReCleanerSusyEWK2L_pt30 --Fs {RP}/leptonBuilderEWK --FMCs {RP}/bTagEventWeightFullSim2L -j 4 -l $L --s2v --tree treeProducerSusyMultilepton --mcc susy-ewkino/crwzcorr/mcc_crwz.txt --mcc susy-ewkino/mcc_triggerdefs.txt  -f   -p data -p prompt_.* -p promptsub.* -p convs.* -p fakes_.* -p rares_.*  -W 'puw_nInt_Moriond(nTrueInt) * leptonSF(getLepSF(LepGood_pt[iF_Recl[0] ], LepGood_eta[iF_Recl[0] ], LepGood_pdgId[iF_Recl[0] ], LepGood_isTight_Recl[iF_Recl[0] ], 1), getLepSF(LepGood_pt[iF_Recl[1] ], LepGood_eta[iF_Recl[1] ], LepGood_pdgId[iF_Recl[1] ], LepGood_isTight_Recl[iF_Recl[1] ], 1), getLepSF(LepGood_pt[iF_Recl[2] ], LepGood_eta[iF_Recl[2] ], LepGood_pdgId[iF_Recl[2] ], LepGood_isTight_Recl[iF_Recl[2] ], 1))*bTagWeight' --load-macro susy-ewkino/3l/functionsEWK.cc --load-macro susy-ewkino/3lcorr/functionsEWKcorr.cc --load-macro susy-ewkino/functionsSF.cc --load-macro susy-ewkino/functionsPUW.cc  --plotgroup fakes_appldata+=promptsub --plotgroup fakes_appldata_ewk_Up+=promptsub_ewk_Up --plotgroup fakes_appldata_ewk_Dn+=promptsub_ewk_Dn  --od $O -o binWZ -b binWZ --sp prompt_WZ --neglist promptsub --neglist promptsub_ewk_Up --neglist promptsub_ewk_Dn
fi
