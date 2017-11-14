#!/bin/bash

T="/mnt/t3nfs01/data01/shome/cheidegg/o/2017-02-23_ewkskims80X_M17_3l_data"
O="/afs/cern.ch/user/c/cheidegg/www/heppy/2017-05-03_ewk80X_eventlists"

python susy-interface/dumpmaker.py 2lss 2lss $T $O -l 35.9 --flags "-X blinding" -p data
#python susy-interface/dumpmaker.py 3l 3lA $T $O -l 35.9 --flags "-X blinding -A met notSR15 'SR!=15'" -p data
python susy-interface/dumpmaker.py 3l 3lB $T $O -l 35.9 --flags "-X blinding" -p data
python susy-interface/dumpmaker.py 3l 3lC $T $O -l 35.9 --flags "-X blinding" -p data
python susy-interface/dumpmaker.py 3l 3lD $T $O -l 35.9 --flags "-X blinding" -p data
python susy-interface/dumpmaker.py 3l 3lE $T $O -l 35.9 --flags "-X blinding" -p data
python susy-interface/dumpmaker.py 3l 3lF $T $O -l 35.9 --flags "-X blinding" -p data
python susy-interface/dumpmaker.py 3l 4lG $T $O -l 35.9 --flags "-X blinding" -p data
python susy-interface/dumpmaker.py 3l 4lH $T $O -l 35.9 --flags "-X blinding" -p data
python susy-interface/dumpmaker.py 3l 4lI $T $O -l 35.9 --flags "-X blinding" -p data
python susy-interface/dumpmaker.py 3l 4lJ $T $O -l 35.9 --flags "-X blinding" -p data
python susy-interface/dumpmaker.py 3l 4lK $T $O -l 35.9 --flags "-X blinding" -p data

