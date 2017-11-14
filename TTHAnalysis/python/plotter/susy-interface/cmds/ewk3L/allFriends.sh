
F="--direct -q all.q --accept WZTo3LNu"
#F="--nosplit --direct -q all.q --accept TChiWZ --accept TChiWH --accept TChiHH --accept TChiHZ --accept TChiZZ4L"
T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-18_ewktrees80X_M17_bkg"
T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-18_ewktrees80X_M17_data"
T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-25_ewkskims80X_M17_ttz_data"
T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-25_ewkskims80X_M17_ttz_bkg"

T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-25_ewkskims80X_M17_3l_bkg"
#T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-25_ewkskims80X_M17_crwz_bkg"
#T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-23_ewkskims80X_M17_crwz_data"
#T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-23_ewkskims80X_M17_3l_data"
T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-18_ewktrees80X_M17_signals"
#T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-05-05_treeMerger"
#T="/mnt/t3nfs01/data01/shome/cheidegg/o/wznewSkimmed"

#python susy-interface/friendmaker.py 3l 3lA $T $T --bk --log $F --modules leptonJetReCleanerSusyEWK3L
#python susy-interface/friendmaker.py 3l 3lA $T $T --bk --log $F --modules leptonJetReCleanerSusyEWK2L
#python susy-interface/friendmaker.py 3l 3lA $T $T --bk --log $F --modules leptonBuilderEWK
#python susy-interface/friendmaker.py 3l 3lA $T $T --bk --log $F --modules bTagEventWeightFullSim3L --flags "-F sf/t {RP}/leptonJetReCleanerSusyEWK3L/evVarFriend_{cname}.root" 
#python susy-interface/friendmaker.py 3l 3lA $T $T --bk --log $F --modules bTagEventWeightFastSim3L --flags "-F sf/t {P}/leptonJetReCleanerSusyEWK3L/evVarFriend_{cname}.root" 
#python susy-interface/friendmaker.py 3l 3lA $T $T --bk --log $F --modules bTagEventWeightFastSim2L --flags "-F sf/t {P}/leptonJetReCleanerSusyEWK2L/evVarFriend_{cname}.root" 
#python susy-interface/friendmaker.py 3l 3lA $T $T --bk --log $F --modules bTagEventWeightFullSim2L --flags "-F sf/t {P}/leptonJetReCleanerSusyEWK2L/evVarFriend_{cname}.root" 
#python susy-interface/friendmaker.py 3l 3lA $T $T --bk --log $F --modules ISRPtWeightModule
#python susy-interface/friendmaker.py 3l 3lA $T $T --bk --log $F --modules PUWProducer2L --flags "-F sf/t {RP}/leptonJetReCleanerSusyEWK2L/evVarFriend_{cname}.root" --exclude "-99" --accept TChiWZ -q all.q
#python susy-interface/friendmaker.py 3l 3lA $T $T --bk --log $F --modules PUWProducer3L --flags "-F sf/t {RP}/leptonJetReCleanerSusyEWK3L/evVarFriend_{cname}.root -F sf/t {RP}/leptonBuilderEWK/evVarFriend_{cname}.root" --exclude "-99" --accept TChi
python susy-interface/friendmaker.py 3l 3lA $T $T --bk --log $F --modules PUWProducer4L --flags "-F sf/t {RP}/leptonJetReCleanerSusyEWK3L/evVarFriend_{cname}.root -F sf/t {RP}/leptonBuilderEWK/evVarFriend_{cname}.root" --exclude "-99" --accept TChiWH --nosplit

T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-04_signalMerger"
#python susy-interface/friendmaker.py 3l 3lA $T $T --bk --log $F --modules leptonJetReCleanerSusyEWK3L --accept TChiWZ --accept TChiSlepSnux0p05 --accept TChiSlepSnux0p95 --exclude ext
#python susy-interface/friendmaker.py 3l 3lA $T $T --bk --log $F --modules leptonBuilderEWK --accept TChiWZ --accept TChiSlepSnux0p05 --accept TChiSlepSnux0p95 --exclude ext


##S="/scratch/cheidegg"
##
###T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-25_ewktrees80X_M17_bkg"
###F="-q all.q --direct"
###R="/pnfs/psi.ch/cms/trivcat/store/user/cheidegg/sea/11/2017-02-25-10-00-00"
##
##T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-18_ewktrees80X_M17_signals"
###T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-23_ewkskims80X_M17_3l_data/"
###T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-25_ewkskims80X_M17_2lss_bkg/"
####T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-25_ewkskims80X_M17_3l_bkg/"
##F="--nosplit --direct --accept TChiStauStaux0p5_225_150 --accept TChiStauStaux0p5ext_225_150"
###F="--nosplit --direct"
###F="-F -q all.q --nosplit --direct --exclude T1 --exclude T6 --exclude T5 --exclude TChiSlepSnux0p"
###R="/pnfs/psi.ch/cms/trivcat/store/user/cheidegg/sea/11/2017-02-05-11-37-00"
###
###F="-F -q all.q --nosplit --direct --accept T1 --accept T5 --accept T6ttWW"
###R="/pnfs/psi.ch/cms/trivcat/store/user/cheidegg/sea/11/2017-02-18-17-00-00"
##
###T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-23_ewkskims80X_M17_3l_data"
###F="-q all.q --nosplit --direct"
##
###T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-25_ewkskims80X_M17_2lss_bkg"
###F="-F -q all.q --nosplit --direct --accept WZTo"
##
#### no remote
###python susy-interface/friendmaker.py 3l 3lA $T $T --bk --log $F --modules leptonJetReCleanerSusyEWK3L
###python susy-interface/friendmaker.py 3l 3lA $T $T --bk --log $F --modules leptonJetReCleanerSusyEWK2L
##python susy-interface/friendmaker.py 3l 3lA $T $T --bk --log $F --modules leptonBuilderEWK
###python susy-interface/friendmaker.py 3l 3lA $T $T --bk --log $F --modules bTagEventWeightFullSim3L --flags "-F sf/t {P}/leptonJetReCleanerSusyEWK3L/evVarFriend_{cname}.root" 
###python susy-interface/friendmaker.py 3l 3lA $T $T --bk --log $F --modules bTagEventWeightFullSim2L --flags "-F sf/t {P}/leptonJetReCleanerSusyEWK2L/evVarFriend_{cname}.root" 
##
#### remote
###python susy-interface/friendmaker.py 3l 3lA $T $S --bk --log $F --flags "--remote '$R'" --modules leptonJetReCleanerSusyEWK3L
###python susy-interface/friendmaker.py 3l 3lA $T $T --bk --log $F --modules leptonBuilderEWK
###python susy-interface/friendmaker.py 3l 3lA $T $S --bk --log $F --flags "--remote '$R'" --modules leptonBuilderEWK
###python susy-interface/friendmaker.py 3l 3lA $T $T --bk --log $F --modules leptonJetReCleanerSusyEWK2L
###python susy-interface/friendmaker.py 3l 3lA $T $S --bk --log $F --flags "--remote '$R'" --modules leptonJetReCleanerSusyEWK2L
###python susy-interface/friendmaker.py 3l 3lA $T $T --bk --log $F --modules TrivialPUWModule
###python susy-interface/friendmaker.py 3l 3lA $T $T --bk --log $F --modules ISRWeightModule
###python susy-interface/friendmaker.py 3l 3lA $T $T --bk --log $F --modules ISRPtWeightModule
###python susy-interface/friendmaker.py 3l 3lA $T $S --bk --log $F --modules bTagEventWeightFullSim3L --flags "-F sf/t {RP}/leptonJetReCleanerSusyEWK3L/evVarFriend_{cname}.root --remote '$R'" 
###python susy-interface/friendmaker.py 3l 3lA $T $S --bk --log $F --modules bTagEventWeightFastSim3L --flags "-F sf/t {RP}/leptonJetReCleanerSusyEWK3L/evVarFriend_{cname}.root --remote '$R'" 
###python susy-interface/friendmaker.py 3l 3lA $T $T --bk --log $F --modules bTagEventWeightFullSim2L --flags "-F sf/t {RP}/leptonJetReCleanerSusyEWK2L/evVarFriend_{cname}.root" 
###python susy-interface/friendmaker.py 3l 3lA $T $S --bk --log $F --modules bTagEventWeightFullSim2L --flags "-F sf/t {RP}/leptonJetReCleanerSusyEWK2L/evVarFriend_{cname}.root --remote '$R'" 
###python susy-interface/friendmaker.py 3l 3lA $T $S --bk --log $F --modules bTagEventWeightFastSim2L --flags "-F sf/t {RP}/leptonJetReCleanerSusyEWK2L/evVarFriend_{cname}.root --remote '$R'" 
##
##
###python susy-interface/friendmaker.py 3l 3lA $T $T --bk --log $F --modules PUWProducer2L --flags "-F sf/t {RP}/leptonJetReCleanerSusyEWK2L/evVarFriend_{cname}.root" --exclude "-99" --accept TChiSlepSnux0p -q all.q
###python susy-interface/friendmaker.py 3l 3lA $T $T --bk --log $F --modules PUWProducer3L --flags "-F sf/t {RP}/leptonJetReCleanerSusyEWK3L/evVarFriend_{cname}.root -F sf/t {RP}/leptonBuilderEWK/evVarFriend_{cname}.root" --exclude "-99" --accept TChi
###python susy-interface/friendmaker.py 3l 3lA $T $T --bk --log $F --modules PUWProducer4L --flags "-F sf/t {RP}/leptonJetReCleanerSusyEWK3L/evVarFriend_{cname}.root -F sf/t {RP}/leptonBuilderEWK/evVarFriend_{cname}.root" --exclude "-99" --accept TChiZZ --accept TChiH
##
##
###T="/mnt/t3nfs01/data01/shome/cheidegg/o/2016-10-04_leptonMVAtraining/trees_80X/"
###F="-q all.q --direct" #"-F" 
###python susy-interface/friendmaker.py 3l 3lA $T $T --bk --log $F --modules mvaSUSY_80X_TTZ_decoupleAndBTagCut
