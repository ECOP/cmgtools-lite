#!/usr/bin/env python
import sys, os, pickle
import ROOT
import os.path as osp

def cacheLocally(infile, tmpDir='/tmp/'):
    tmpfile = osp.join(tmpDir, osp.basename(infile))

    # Copy locally if it's not there already
    if not osp.exists(tmpfile):
        xrcmd = "xrdcp %s %s" % (infile, tmpfile)
        print " transferring to %s" % tmpDir
        os.system(xrcmd)
        print "... copied successfully"

    infile = tmpfile
    return infile

_treepath = None
_allfiles = []

def get_file_or_copy_local(url, copy_local=True):
    if os.path.exists(url):
        return url
    if not os.path.exists('%s.url'%url):
        raise RuntimeError('File not found: %s'%url)

    url = open('%s.url'%url, 'r').read().strip()
    if "#dcap" in url:
        url=url.split("\n")[0]
    #print "---->> ",url,"<<<<"
    if not copy_local:
        return url

    else:
        #from CMGTools.TTHAnalysis.macros.leptons.LepMVAEfficiencies.runLepTnPFriendMaker import cacheLocally
        return cacheLocally(url, os.environ.get('TMPDIR', '/tmp'))


def load_dataset(name, trainclass, addw, path=None, friends=[]):
    if not path:
        path = _treepath
        if not path:
            raise RuntimeError('No default tree path set')

    fileloc = get_file_or_copy_local(
                os.path.join(path, name, 'treeProducerSusyMultilepton/tree.root'))
    print 'Using %s' % fileloc
    infile = ROOT.TFile.Open(fileloc)
    _allfiles.append(infile) # Dirty trick to keep the files and trees in memory

    tree = infile.Get('mt2')#tree

    # Check if tree was loaded
    try: tree.GetName()
    except ReferenceError:
        raise RuntimeError("'tree' not found in %s" % fileloc)

    for friend in friends:
        friendloc = get_file_or_copy_local(os.path.join(path, friend, 'evVarFriend_%s.root' % name))
        print "Adding friend from", friendloc
        tree.AddFriend('sf/t', friendloc)

    pckfile = os.path.join(path, name, "skimAnalyzerCount/SkimReport.pck")
    pckobj  = pickle.load(open(pckfile,'r'))
    counters = dict(pckobj)
    weight = 1.0*addw/(counters['Sum Weights'])
    print ('Added %s dataset, category %s, with weight %f/%f' %
             (name, trainclass, addw, counters['Sum Weights']))

    return tree, weight

def train_multiclass(fOutName, options):
    dsets = [
        ('TTHnobb_pow', 'ttH', 1),
        ('TTW_LO', 'ttV', 1),
        ('TTZ_LO', 'ttV', 1),
    ]

    if '_3l' in options.training:
        dsets += [
            ('TTJets_DiLepton',            'tt', 1),
#            ('TTJets_DiLepton_ext_skim3l', 'tt', 0.9),
        ]
    else:
        dsets += [
            ('TTJets_SingleLeptonFromT',        'tt', 0.1),
            ('TTJets_SingleLeptonFromTbar',     'tt', 0.1),
            ('TTJets_SingleLeptonFromT_ext',    'tt', 0.9),
            ('TTJets_SingleLeptonFromTbar_ext', 'tt', 0.9),
        ]

    datasets = []
    for name, trainclass, addw, cut in dsets:
        tree, weight = load_dataset(name, trainclass, addw,
                                    path=options.treepath,
                                    friends=options.friends)
        datasets.append((name, trainclass, tree, weight, cut))

    fOut = ROOT.TFile(fOutName,"recreate")
    fOut.cd()
    factory = ROOT.TMVA.Factory(options.training, fOut, "!V:!Color:Transformations=I:AnalysisType=Multiclass")
    allcuts = ROOT.TCut('1')
    for cut in options.addcuts:
        allcuts += cut

    allcuts += "nLepFO_Recl>=2"
    allcuts += "LepGood_conePt[iLepFO_Recl[0]]>25"
    allcuts += "LepGood_conePt[iLepFO_Recl[1]]>15"

    allcuts += "abs(mZ1_Recl-91.2) > 10"
    allcuts += "(met_pt*0.00397 + mhtJet25_Recl*0.00265 > 0.2)"
    allcuts += "(nBJetLoose25_Recl >= 2 || nBJetMedium25_Recl >= 1)"
    allcuts += "minMllAFAS_Recl>12"

    if '_3l' in options.training:
        allcuts += "nLepFO_Recl>=3"
        allcuts += "LepGood_conePt[iLepFO_Recl[2]]>15"
        allcuts += "nJet25_Recl>=2"
        # allcuts += "LepGood_isTight_Recl[iLepFO_Recl[0]]"
        # allcuts += "LepGood_isTight_Recl[iLepFO_Recl[1]]"
        # allcuts += "LepGood_isTight_Recl[iLepFO_Recl[2]]"
    else:
        allcuts += "nLepTight_Recl<=2"
        allcuts += "nJet25_Recl>=4"
        allcuts += "(LepGood_charge[iLepFO_Recl[0]]*LepGood_charge[iLepFO_Recl[1]] > 0)" #!
        # allcuts += "LepGood_isTight_Recl[iLepFO_Recl[0]]"
        # allcuts += "LepGood_isTight_Recl[iLepFO_Recl[1]]"

    factory.AddSpectator("iF0 := iLepFO_Recl[0]","F") # do not remove this!
    factory.AddSpectator("iF1 := iLepFO_Recl[1]","F") # do not remove this!
    factory.AddSpectator("iF2 := iLepFO_Recl[2]","F") # do not remove this!

    factory.AddVariable("higher_Lep_eta := max(abs(LepGood_eta[iLepFO_Recl[0]]),abs(LepGood_eta[iLepFO_Recl[1]]))", 'F')
    factory.AddVariable("MT_met_lep1 := MT_met_lep1", 'F')
    factory.AddVariable("numJets_float := nJet25_Recl", 'F')
    factory.AddVariable("mindr_lep1_jet := mindr_lep1_jet", 'F')
    factory.AddVariable("mindr_lep2_jet := mindr_lep2_jet", 'F')
    factory.AddVariable("LepGood_conePt[iLepFO_Recl[0]] := LepGood_conePt[iLepFO_Recl[0]]", 'F')
    factory.AddVariable("LepGood_conePt[iLepFO_Recl[1]] := LepGood_conePt[iLepFO_Recl[1]]", 'F')
    factory.AddVariable("avg_dr_jet : = avg_dr_jet", 'F')
    factory.AddVariable("met := min(met_pt, 400)", 'F')

    ## Add the datasets
    for name,trainclass,tree,weight,cut in datasets:
        factory.AddTree(tree, trainclass, weight,cut)

    fOut.cd()
    for trainclass in set([x[1] for x in dsets]):
        factory.SetWeightExpression("genWeight*xsec", trainclass)

    ## Start the training
    factory.PrepareTrainingAndTestTree(allcuts, "!V")
    factory.BookMethod(ROOT.TMVA.Types.kBDT, 'BDTG',
                            ':'.join([
                                '!H',
                                '!V',
                                'NTrees=200',
                                # 'NTrees=500',
                                'BoostType=Grad',
                                'Shrinkage=0.10',
                                '!UseBaggedGrad',
                                'nCuts=200',
                                # 'nCuts=2000',
                                'nEventsMin=100',
                                'MaxDepth=8',
                                'NegWeightTreatment=PairNegWeightsGlobal',
                                ]))
    factory.TrainAllMethods()
    factory.TestAllMethods()
    factory.EvaluateAllMethods()

    fOut.Close()

def train_single(allcuts, variables, dsets, fOutName, options, spectators=[]):
    datasets = []
    for name, trainclass, addw, cut in dsets:
        #print "--->>>> ",options.treepath,options.friends
        tree, weight = load_dataset(name, trainclass, addw,
                                    path=options.treepath,
                                    friends=options.friends)
        datasets.append((name, trainclass, tree, weight,ROOT.TCut(cut)))

    fOut = ROOT.TFile(fOutName,"recreate")
    fOut.cd()
    factory = ROOT.TMVA.Factory(options.training, fOut, "!V:!Color:Transformations=I")

    for cut in options.addcuts:
        allcuts += cut

    #factory.AddSpectator("iF0 := iLepFO_Recl[0]","F") # do not remove this!
    #factory.AddSpectator("iF1 := iLepFO_Recl[1]","F") # do not remove this!
    #factory.AddSpectator("iF2 := iLepFO_Recl[2]","F") # do not remove this!

    for spec in spectators:
        factory.AddSpectator(spec)

    #factory.AddSpectator("iJ0 := ","F")
    #factory.AddSpectator("iJ1 := ","F")

    ## Add the variables
    for var in variables:
        factory.AddVariable(var, 'F')

    ## Add the datasets
    for name,trainclass,tree,weight,cut in datasets:
        factory.AddTree(tree, trainclass, weight,cut)

    fOut.cd()
    #for trainclass in set([x[1] for x in dsets]):
    #    factory.SetWeightExpression("genWeight*xsec", trainclass)

    ## Start the training
    print "======= prepare datasets ========"
    factory.PrepareTrainingAndTestTree(allcuts, "V")#:nTrain_Signal=1000:nTest_Signal=1000:nTrain_Background=1000:nTest_Background=1000")
    print "======= datasets ready for training and testing ========"
    factory.BookMethod(ROOT.TMVA.Types.kBDT, 'BDTG',
                            ':'.join([
                                '!H',
                                'V',
                                'NTrees=200',
                                'BoostType=Grad',
                                'Shrinkage=0.10',
                                '!UseBaggedGrad',
                                'nCuts=200',
                                'nEventsMin=100',
                                'NNodesMax=5',
                                'MaxDepth=8',
                                'NegWeightTreatment=PairNegWeightsGlobal',
                                'CreateMVAPdfs',
                                ]))
    print "=============== method booked, now training ========================================="
    factory.TrainAllMethods()
    print "===============  training done, now testing ========================================="
    factory.TestAllMethods()
    print "===============  testing done, now evaluation ========================================="
    factory.EvaluateAllMethods()
    print "=============================== successful! ==========================================="
    fOut.Close()

def train_2d(fOutName, training, options):
    allcuts = ROOT.TCut('1')
    if '2lss' in training:
        allcuts += "nLepFO_Recl>=2"
        allcuts += "LepGood_conePt[iLepFO_Recl[0]]>25"
        allcuts += "LepGood_conePt[iLepFO_Recl[1]]>15"
        allcuts += "LepGood_charge[iLepFO_Recl[0]] == LepGood_charge[iLepFO_Recl[1]]"
        allcuts += "(nBJetLoose25_Recl >= 2 || nBJetMedium25_Recl >= 1)"
        allcuts += "nJet25_Recl >= 4"
    elif '3l' in training:
        allcuts += "nLepFO_Recl>=3"
        allcuts += "abs(mZ1_Recl-91.2)>10"
        allcuts += "LepGood_conePt[iLepFO_Recl[0]]>25"
        allcuts += "LepGood_conePt[iLepFO_Recl[1]]>15"
        allcuts += "LepGood_conePt[iLepFO_Recl[2]]>15"
        allcuts += "(nJet25_Recl >= 4 || (met_pt*0.00397 + mhtJet25_Recl*0.00265 - 0.184 > 0.0 + 0.1*(mZ1_Recl > 0)))"
        allcuts += "nBJetLoose25_Recl >= 2"

    variables = [ # Common variables
        "max_Lep_eta := max(abs(LepGood_eta[iLepFO_Recl[0]]),abs(LepGood_eta[iLepFO_Recl[1]]))",
        "numJets_float := nJet25_Recl",
        "mindr_lep1_jet := mindr_lep1_jet",
        "mindr_lep2_jet := mindr_lep2_jet",
        "MT_met_lep1 := MT_met_lep1",
    ]
    
    dsets = []
    if '3l' in training and 'ttv' in training and 'mem' in training:
        if not 'skim3l2j2b1B' in options.treepath: raise RuntimeError
        dsets += [('TTHnobb_pow', 'Signal', 6.)]
    else:
        dsets += [('TTHnobb_pow', 'Signal', 1)]

    if '2lss' in training and 'ttw' in training:
        variables += ["LepGood_conePt[iLepFO_Recl[1]] := LepGood_conePt[iLepFO_Recl[1]]"]
        dsets += [('TTW_LO', 'Background', 1)]

    if '2lss' in training and 'ttv' in training:
        variables += [
            "LepGood_conePt[iLepFO_Recl[1]] := LepGood_conePt[iLepFO_Recl[1]]",
            "LepGood_conePt[iLepFO_Recl[0]] := LepGood_conePt[iLepFO_Recl[0]]"
        ]
        dsets += [
            ('TTW_LO', 'Background', 1),
            ('TTZ_LO', 'Background', 1),
        ]
    if '2lss' in training and 'ttbar' in training:
        variables += [
#            "met := min(met_pt, 400)",
#            "avg_dr_jet : = avg_dr_jet",
        ]
        dsets += [
            ('TTJets_SingleLeptonFromT',        'Background', 0.2),
            ('TTJets_SingleLeptonFromTbar',     'Background', 0.2),
            ('TTJets_SingleLeptonFromT_ext',    'Background', 0.8),
            ('TTJets_SingleLeptonFromTbar_ext', 'Background', 0.8),
        ]


    if '3l' in training and 'ttw' in training:
        variables += ["LepGood_conePt[iLepFO_Recl[2]] := LepGood_conePt[iLepFO_Recl[2]]"]
        dsets += [('TTW_LO', 'Background', 1)]
    if '3l' in training and 'ttv' in training:
        variables += [
            "LepGood_conePt[iLepFO_Recl[2]] := LepGood_conePt[iLepFO_Recl[2]]",
            "LepGood_conePt[iLepFO_Recl[0]] := LepGood_conePt[iLepFO_Recl[0]]"
        ]
        if "mem" in training:
            if not 'skim3l2j2b1B' in options.treepath: raise RuntimeError
            dsets += [
                ('TTW_LO', 'Background', 4.),
                ('TTZ_LO', 'Background', 26.),
                ]
        else:
            dsets += [
                ('TTW_LO', 'Background', 1),
                ('TTZ_LO', 'Background', 1),
                ]
    if '3l' in training and 'ttbar' in training:
        variables += [
#            "mhtJet25 := mhtJet25_Recl",
#            "avg_dr_jet : = avg_dr_jet",
        ]
        dsets += [
            ('TTJets_DiLepton',            'Background', 1.0/6),
            ('TTJets_DiLepton_ext_part1',            'Background', 2.0/6),
            ('TTJets_DiLepton_ext_part2',            'Background', 2.0/6),
            ('TTJets_DiLepton_ext_part3',            'Background', 1.0/6),
            ('TTJets_SingleLeptonFromT',        'Background', 0.2),
            ('TTJets_SingleLeptonFromTbar',     'Background', 0.2),
            ('TTJets_SingleLeptonFromT_ext',    'Background', 0.8),
            ('TTJets_SingleLeptonFromTbar_ext', 'Background', 0.8),
        ]

    if 'bdtv8_bestchoice' in training:
        variables += [
            'BDTv8_eventReco_mvaValue := max(-0.2,BDTv8_eventReco_mvaValue)',
            "BDTv8_eventReco_bJet_fromHadTop_CSV := max(-0.2,BDTv8_eventReco_bJet_fromHadTop_CSV)",
            "BDTv8_eventReco_HadTop_pT := max(-10,BDTv8_eventReco_HadTop_pT)",
            "BDTv8_eventReco_HadTop_mass := max(-10,BDTv8_eventReco_HadTop_mass)",
            ]
    if 'bdtv8_onlymass' in training:
        variables += [
            'BDTv8_eventReco_mvaValue := max(-0.2,BDTv8_eventReco_mvaValue)',
            "BDTv8_eventReco_HadTop_mass := max(-10,BDTv8_eventReco_HadTop_mass)",
            ]
    if 'bdtv8_value' in training:
        variables += [
            'BDTv8_eventReco_mvaValue := max(-1.1,BDTv8_eventReco_mvaValue)',
            ]
    if 'bdtv8_reco' in training:
        variables += [
            "BDTv8_eventReco_bJet_fromLepTop_CSV := max(-1.1,BDTv8_eventReco_bJet_fromLepTop_CSV)",
            "BDTv8_eventReco_bJet_fromHadTop_CSV := max(-1.1,BDTv8_eventReco_bJet_fromHadTop_CSV)",
            "BDTv8_eventReco_qJet1_fromW_fromHadTop_CSV := max(-1.1,BDTv8_eventReco_qJet1_fromW_fromHadTop_CSV)",
            "BDTv8_eventReco_HadTop_pT := BDTv8_eventReco_HadTop_pT",
            "BDTv8_eventReco_W_fromHadTop_mass := BDTv8_eventReco_W_fromHadTop_mass",
            "BDTv8_eventReco_HadTop_mass := BDTv8_eventReco_HadTop_mass",
            "BDTv8_eventReco_W_fromHiggs_mass := BDTv8_eventReco_W_fromHiggs_mass",
            "BDTv8_eventReco_LepTop_HadTop_dR := BDTv8_eventReco_LepTop_HadTop_dR",
            ]
    if 'hadtopsimple' in training:
        variables += [
            "HadTopSimple_bJet_fromHadTop_CSV := max(-1.1,bJet_fromHadTop_CSV)",
            "HadTopSimple_lJet_fromHadTop_CSV2 := max(-1.1,lJet_fromHadTop_CSV2)",
            "HadTopSimple_HadTop_Mass := HadTop_Mass",
            "HadTopSimple_HadTop_Pt := HadTop_Pt",
            "HadTopSimple_W_fromHadTop_Mass := W_fromHadTop_Mass",
            "HadTopSimple_bJet_notFromHadTop_CSV := max(-1.1,bJet_notFromHadTop_CSV)"
            ]
    if 'memlogs' in training:
        variables += [
            "MEM_TTHvsTTW := max(-10.,min(50.,log(MEM_TTH)-log(MEM_TTW)))",
            "MEM_TTHvsTTZ := max(-10.,min(50.,log(MEM_TTH)-log(MEM_TTLL)))",
            ]
    if 'memvars' in training:
        variables += [
            "MEM_TTH := max(-100.,min(10.,log(MEM_TTH)))",
            "MEM_TTW := max(-100.,min(10.,log(MEM_TTW)))",
            "MEM_TTZ := max(-100.,min(10.,log(MEM_TTLL)))",
            ]
    if 'memfixvars' in training:
        variables += [
            "MEM_TTH := min(0,log(max(3.72e-44,MEM_TTH)))",
            "MEM_TTW := min(0,log(max(3.72e-44,MEM_TTW)))",
            "MEM_TTZ := min(0,log(max(3.72e-44,MEM_TTLL)))",
            ]
    if 'memlr' in training:
        variables += [
            "MEM_LR := -log((0.00389464*MEM_TTLL*(MEM_TTLL<1) + 3.12221e-14*MEM_TTW*(MEM_TTW<1)) / (0.00389464*MEM_TTLL*(MEM_TTLL<1) + 3.12221e-14*MEM_TTW*(MEM_TTW<1)+9.99571e-05*(MEM_TTHfl*(MEM_TTHfl<1)+MEM_TTHsl*(MEM_TTHsl<1))/2))"
            ]
    if 'hj_value' in training:
        variables += [
            'BDTv8_eventReco_Hj_score := max(-1.1,BDTv8_eventReco_Hj_score)',
#            'BDTv8_eventReco_Hjj_score := max(-1.1,BDTv8_eventReco_Hjj_score)',
            ]

    outname = fOutName+'_'+training+'.root'
    train_single(allcuts, variables, dsets, outname, options)



def train_top(fOutName, training, options):
    allcuts = ROOT.TCut('1')
    allcuts += "(njet >= 3)"
    
    variables = [ # Common variables
        "numJets_float := njet",

        "lJetPt  := jet_pt[0]",
        "slJetPt := jet_pt[1]",
        #"JetPt3 := jet_pt[2]",
        #"JetPt4 := jet_pt[3]",
        #"JetPt5 := jet_pt[4]",
        #"JetPt6 := jet_pt[5]",

        "lBTag  := Max$(jet_btagCSV)",
        "slBTag := MaxIf$(jet_btagCSV, jet_btagCSV!=Max$(jet_btagCSV) )",
    ]

    dsets = []
    dsets += [('top', 'Signal',1., '(Sum$(jet_mcMatchId==6)==3 || Sum$(jet_mcMatchId==-6)==3)')]
    #if '3l' in training and 'wz' in training:
    dsets += [('top', 'Background',1., '!(Sum$(jet_mcMatchId==6)==3 || Sum$(jet_mcMatchId==-6)==3)')]

    outname = fOutName+'_'+training+'.root'
    train_single(allcuts, variables, dsets, outname, options)



def train_ttZ(fOutName, training, options):
    allcuts = ROOT.TCut('1')
    allcuts += "nLepTight_Recl>=3"
    allcuts += "abs(mZ_Extra-91.2)>25"
    #allcuts += "LepGood_Pt[ iLepFO_Recl[0] ]>25"
    #allcuts += "LepGood_Pt[ iLepFO_Recl[1] ]>15"
    #allcuts += "LepGood_Pt[ iLepFO_Recl[2] ]>15"
    allcuts += "(nJetSel_Recl >= 2)"
    #allcuts += "Entry$<10000"
    
    variables = [ # Common variables
        #"max_Lep_eta := max(abs(LepGood_eta[iLepFO_Recl[0]]),abs(LepGood_eta[iLepFO_Recl[1]]))",
        "numJets_float := nJetSel_Recl",
        #"mindr_lep1_jet := mindr_lep1_jet",
        #"mindr_lep2_jet := mindr_lep2_jet",
        #"MT_met_lep1 := MT_met_lep1",

        "lJetPt  := JetSel_Recl_pt[0]",
        "slJetPt := JetSel_Recl_pt[1]",
        
        "lBTag  := Max$(JetSel_Recl_btagCSV)",
        "slBTag := MaxIf$(JetSel_Recl_btagCSV, JetSel_Recl_btagCSV!=Max$(JetSel_Recl_btagCSV) )",

    ]

    spectators=[
        "iJPt0 := JetSel_Recl_pt[0]",
        "iJPt1 := JetSel_Recl_pt[1]",
        ]


    #if '3l' in training and 'WZ' in training:
    #    variables += []
    if '3l' in training and 'tth' in training:
        variables += [ 
            "mZ1TT := mZ1TT_Recl"
            "mindr_lep1_jet := mindr_lep1_jet",
            "mindr_lep2_jet := mindr_lep2_jet",
            #"MT_met_lep1 := MT_met_lep1",
            "max_Lep_eta := max(abs(LepGood_eta[iT_Recl[0]]),abs(LepGood_eta[iT_Recl[1]]))"
            ]
        spectators += [
            "iL0 := iT_Recl[0]",
            "iL1 := iT_Recl[1]"
            ]
    if '3l' in training and 'ttx' in training:
        variables += [
            #"mZ1TT := mZ1TT_Recl",
            #"max_Lep_eta := max(abs(LepGood_eta[iT_Recl[0]]),abs(LepGood_eta[iT_Recl[1]]))",
            "mindr_lep1_jet := mindr_lep1_jet_Recl",
            "mindr_lep2_jet := mindr_lep2_jet_Recl",
            "max_Lep_eta := max(abs(LepGood_eta[iT_Recl[0]]),abs(LepGood_eta[iT_Recl[1]]))",
            "mZ := mZ_Extra",
            "pTZ := pTZ_Extra",
            "mtWLep := mTWLep_Extra",
            "pTWLep := pTWLep_Extra",
            "minDRWLepBJet := minDRWLepBJet_Extra",
            "mWHad1 := mWHad1_Extra",
            "pTWHad1 := pTWHad1_Extra",
            "minDrWHadBJet1 := minDrWHadBJet1_Extra",
            "mWHad2 := mWHad2_Extra",
            "pTWHad2 := pTWHad2_Extra",
            "minDrWHadBJet2 := minDrWHadBJet2_Extra",
            "mTHad := mTHad_Extra",
            "pTTHad := pTTHad_Extra",
            "mTSLep := mTSLep_Extra",
            "pTSLep := pTSLep_Extra",
            ]
        spectators += [
            "iL0 := iT_Recl[0]",
            "iL1 := iT_Recl[1]"
            ]

    dsets = []
    dsets += [('TTZToLLNuNu_ext2', 'Signal', 1.)]
    if '3l' in training and 'wz' in training:
        dsets += [('WZTo3LNu', 'Background', 1.)]
    if '3l' in training and 'tth' in training:
        dsets += [('TTHnobb_pow', 'Background', 1.)]
    if '3l' in training and 'ttx' in training:
        dsets += [('tZq_ll_tiny', 'Background', 1.)]
        


    outname = fOutName+'_'+training+'.root'
    train_single(allcuts, variables, dsets, outname, options)

def main(args, options):
    global _treepath
    _treepath = options.treepath
    if 'MultiClassICHEP16' in options.training:
        train_multiclass(args[0]+'.root', options)
        return
    #print "-->>",options.training,options.training.lower()
    if len(options.training):
        #train_2d(args[0], options.training.lower(), options)
        train_ttZ(args[0], options.training.lower(), options)
    else:
        #train_2d(args[0], '2lss_ttv',   options)
        #train_2d(args[0], '2lss_ttbar', options)
        #train_2d(args[0], '3l_ttv',     options)
        #train_2d(args[0], '3l_ttbar',   options)
        #train_ttZ(args[0], '3l_ttvVsWZ', options)
        train_top(args[0], 'topBDTvsTop', options)


if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] OUTFILE")
    parser.add_option("-T","--training", dest="training",
                      type="string", default="",
                      help=('Either "MultiClassICHEP16" or "2lss_ttw",'
                            '"3l_ttbar", etc. Default will run all 2D trainings'))
    parser.add_option("-P","--treepath", dest="treepath",
                      type="string", default=None)
    parser.add_option("-F","--friend", dest="friends",
                      type="string", default=[], action="append")
    parser.add_option("-c","--cut", dest="addcuts",
                      type="string", default=[], action="append")
    (options, args) = parser.parse_args()

    if not len(args) or not options.treepath:
        parser.print_help()
        sys.exit(-1)

    main(args, options)
