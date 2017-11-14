import os, sys
from optparse import OptionParser
from lib import maker
from lib import functions as func

def deformat(line):
	line = line.replace("+", "_")
	line = line.replace("#", "_")
	line = line.replace(":", "_")
	line = line.replace("*", "_")
	return line

def getXS(xs, mass, factor):
	themass = [p[0] for p in xs]
	thexs   = [p[1] for p in xs]
	return str(thexs[themass.index(float(mass))]) + "*" + str(factor)

def prepareJob(mm, name, mp, baseSig, baseSys, binning, bkgpath, outpath, xslist, options):

	mp = [m.strip() for m in mp]
	#name = "mp_"+mm.model.name+"_"+mp[0]+"_"+mp[1]
	sig  = mm.model.name +"_"+ mp[0] +"_"+ mp[1]
	xs = getXS(xslist, float(mp[0]), mm.model.brcorr)

	cardpath = (outpath +"/mps/"+mp[0]+"_"+mp[1], "sig_"+sig+"/SR.card.txt")
	source   = open("susy-interface/scripts/job_scanmaker.py", "r").readlines()
	tmpfile  = open(mm.srcpath +"/submitJob_"+ name +".py", "w")

	thedummy = mm.tmppath +"/cuts_dummy.txt"
	fdummy   = open(thedummy, "w")
	fdummy.write("alwaystrue : 1")
	fdummy.close()

	for line in source:
		line = line.replace("THENAME"      , name                                                    )
		line = line.replace("THESIGNAL"    , sig                                                     )
		line = line.replace("THEMASS1"     , mp[0]                                                   )
		line = line.replace("THEMASS2"     , mp[1]                                                   )
		line = line.replace("THEOFFSET"    , func.getOffset(mm.getVariable("expr",""), binning)      )
		line = line.replace("THETREEDIR"   , mm.treedir                                              )
		line = line.replace("THETREENAME"  , mm.getVariable("treename","treeProducerSusyMultilepton"))
		line = line.replace("THEFILE"      , mp[2]                                                   )
		line = line.replace("THEXS"        , xs                                                      ) 
		#line = line.replace("THEQ2FILE"    , mm.getVariable("q2accfile","")                          )
		#line = line.replace("THEQ2SYNTAX"  , mm.getVariable("q2accformat","")                        )
		line = line.replace("THEWEIGHTSTR" , mm.getVariable("mcaWeightFS","1.0")                     )
		line = line.replace("THEFRFILES"   , "["+",".join(["\""+f+"\"" for f in mm.getVariable("frFilesFS","").split(";")])+"]")
		line = line.replace("THEJEC"       , mm.getVariable("jec","")                                )
		line = line.replace("THEMET"       , mm.getVariable("met","")                                )
		line = line.replace("THEQ2ACC"     , mm.getVariable("q2acc","")                              )
		#line = line.replace("THEPUW"       , mm.getVariable("puw","")                                )
		line = line.replace("THEFRJEC"     , "["+",".join(["\""+f+"\"" for f in mm.getVariable("frFilesFSJec","").split(";")])+"]")
		line = line.replace("THEWVJEC"     , "["+",".join(["\""+f+"\"" for f in mm.getVariable("wVarsFSJec"  ,"").split(";")])+"]")
		line = line.replace("THEFRMET"     , "["+",".join(["\""+f+"\"" for f in mm.getVariable("frFilesFSMet","").split(";")])+"]")
		line = line.replace("THEWVMET"     , "["+",".join(["\""+f+"\"" for f in mm.getVariable("wVarsFSMet"  ,"").split(";")])+"]")
		line = line.replace("THEWEIGHTVARS", "{"+",".join(["\""+k+"\":[" + ",".join("\""+v+"\"" for v in vals) + "]" for k,vals in mm.getVariable("wVarsFS",{}).iteritems()])+"}")
		line = line.replace("THEMCADIR"    , mm.tmppath                                        )
		line = line.replace("THEBKGDIR"    , bkgpath                                           )
		line = line.replace("THEOUTDIR"    , outpath                                           )
		line = line.replace("THEMCA"       , mm.getVariable("mcafile","")                      )
		line = line.replace("THECUTS"      , mm.getVariable("cutfile","")                      )
		line = line.replace("THEDUMMY"     , thedummy                                          )
		line = line.replace("THESYST"      , mm.getVariable("sysfile","")                      )
		line = line.replace("THEBASESIG"   , baseSig                                           )
		line = line.replace("THEBASESYS"   , baseSys                                           )
		tmpfile.write(line)
	tmpfile.close()
	return "python "+ mm.srcpath +"/submitJob_"+ name +".py", cardpath

	
parser = OptionParser(usage="%prog cfg regions treedir outdir [options]")
parser = maker.addMakerOptions(parser)
parser.add_option("--WFS"         , dest="weightFS"   , type="string"      , default=None , help="Overwrite the weightFS expression");
#parser.add_option("--q2accfile"   , dest="q2accfile"  , type="string"      , default=None , help="Overwrite the q2accfile expression");
#parser.add_option("--q2accformat" , dest="q2accformat", type="string"      , default=None , help="Overwrite the q2accfile expression");
parser.add_option("--bkgOnly"     , dest="bkgOnly"    , action="store_true", default=False, help="Only run the bkg only");
parser.add_option("--sigOnly"     , dest="sigOnly"    , action="store_true", default=False, help="Only run the signal (if bkg already is present)");
parser.add_option("--perBin"      , dest="perBin"     , action="store_true", default=False, help="Make datacards for every bin in 'expr' separately.");
parser.add_option("-m", "--models", dest="models"     , action="append"    , default=[]   , help="Fastsim signal models to loop upon.");
parser.add_option("--redoBkg"     , dest="redoBkg"    , action="store_true", default=False, help="Redo bkg if it already exists.");
parser.add_option("--postfix"     , dest="postfix"    , type="string"      , default=""   , help="Flags to be run for postfix (i.e. for signal only)");
parser.add_option("--check"       , dest="doCheck"    , action="store_true", default=False, help="Automatically run some (systematics) checks of the datacards that were produced.");
parser.add_option("--m1"          , dest="mass1"      , type="int"         , default=None , help="Only run signal points with this value as first mass");
parser.add_option("--m2"          , dest="mass2"      , type="int"         , default=None , help="Only run signal points with this value as second mass");
parser.add_option("--allowoverflow", dest="allowoverflow", action="store_true", default=False, help="Allow last bin to be overflow bin");

baseBkg = "python makeShapeCardsSusy.py {MCA} {CUTS} \"{EXPR}\" \"{BINS}\" -o SR --bin {TAG} {T} --tree {TREENAME} {MCCS} {MACROS} --s2v -f -l {LUMI} --od {O} {FRIENDS} {FLAGS} {OVERFLOWCUTS}"
baseSig = "python makeShapeCardsSusy.py [[[MCA]]] [[[CUTS]]] \\\"{EXPR}\\\" \\\"{BINS}\\\" [[[SYS]]] -o SR --bin {TAG} {T} --tree {TREENAME} {MCCS} {MACROS} --s2v -f -l {LUMI} --od [[[O]]] {FRIENDS} {FLAGS} {OVERFLOWCUTS} {POSTFIX}"
(options, args) = parser.parse_args()
options         = maker.splitLists(options)
options.models  = func.splitList(options.models)
mm              = maker.Maker("scanmaker", baseBkg, args, options, parser.defaults)
mm.loadModels()

sl      = mm.getVariable("lumi","12.9").replace(".","p")

combinePath = mm.getVariable("combineTool")
if options.doCheck and (not combinePath or not os.path.isdir(combinePath) or combinePath[-3:] != "src"):
	combinePath = combinePath.rstrip("/")
	print "WARNING: You want to run the automatic systematics check for your datacards"
	print "but you have not given all necessary information in your configuration!"
	print "Please give the path to the 'src' directory (e.g. /usr/cheidegg/CMSSW_7_1_5/src)"
	print "of the CMSSW release where HiggsCombine is installed. You should store this"
	print "path in your configuration with the variable 'combineTool'"
	cont = raw_input("Do you want to continue now without the systematics check? [y/n]\n>> ")
	if not cont in ["y"]: sys.exit()

## first do bkg
if not options.sigOnly:

	friends = mm.collectFriends(False)

	mm.reloadBase(baseBkg)
	mm.resetRegion()

	## looping over regions
	for r in range(len(mm.regions)):
		mm.iterateRegion()
		
		scenario = mm.getScenario  (True)
		mccs     = mm.collectMCCs  ()
		macros   = mm.collectMacros()	
		flags    = mm.collectFlags (["flagsScans"], True, True, False, True)

		## looping over binnings
		binnings = [mm.getVariable("bins","")] if not options.perBin else func.getAllBins(mm.getVariable("bins",""))
		for ib,b in enumerate(binnings):
			theSc  = scenario
			theBin = mm.getVariable("bins","")
			
			## change scenario if looping over all bins
			if options.perBin: 
				min, max = func.getMinMax(b)
				theSc    = scenario + "_" + min.replace(".","p")
				theBin   = b		

			## background first
			output = mm.outdir +"/scan/"+ theSc +"/"+ sl +"fb"
			bkgDir = output +"/bkg"
			if not options.redoBkg and os.path.exists(bkgDir+"/common/SR.input.root"): continue
			func.mkdir(bkgDir)
		
			mm.submit([mm.getVariable("mcafile",""), mm.getVariable("cutfile",""), mm.getVariable("expr",""), theBin, theSc.replace("/","_"), mm.treedirs, mm.getVariable("treename","treeProducerSusyMultilepton"), mccs, macros, mm.getVariable("lumi","12.9"), bkgDir, friends, flags, func.getCut(options, mm.getVariable("firstCut","alwaystrue"), mm.getVariable("expr",""), theBin)], theSc.replace("/", "_")+"_bkg",False)
	mm.runJobs()
	mm.clearJobs()
		


## second do models
if not options.bkgOnly:

	cards = []
	
	friends = mm.collectFriends(True)

	mm.reloadBase(baseSig)
	mm.resetRegion()

	## looping over regions
	for r in range(len(mm.regions)):
		mm.iterateRegion()
	
		scenario = mm.getScenario  (True)
		mccs     = mm.collectMCCs  ()
		macros   = mm.collectMacros()	
		flagsSig = mm.collectFlags (["flagsScans"              ], True, True, True, True)
		flagsSys = mm.collectFlags (["flagsScans", "flagsSysts"], True, True, True, True)
	
		## looping over binnings
		binnings = [mm.getVariable("bins","")] if not options.perBin else func.getAllBins(mm.getVariable("bins",""))
		for ib,b in enumerate(binnings):
			theSc = scenario
			theBin = mm.getVariable("bins","")
			
			## change scenario if looping over all bins
			if options.perBin: 
				min, max = func.getMinMax(b)
				theSc    = scenario + "_" + min.replace(".","p")
				theBin   = b		
		
			## background first
			output = mm.outdir +"/scan/"+ theSc +"/"+ sl +"fb" 
	
			## looping over models
			mm.resetModel()

			for m in range(len(mm.models)):
				mm.iterateModel()
	
				myDir  = output +"/"+ mm.model.name 
				bkgDir = output +"/bkg"
				func.mkdir(myDir +"/acc")
				func.mkdir(myDir +"/mps")
	
				## prepare jobs for masspoints
				xslist = [l.rstrip("\n").split(":") for l in open(mm.model.xsecfile  , "r").readlines()]
				xslist = [[float(m.strip()),float(xs.strip())] for [m,xs,err] in xslist ]
				mps    = [l.rstrip("\n").split(":") for l in open(mm.model.masspoints, "r").readlines()]
				mps    = [[m[0].strip(), m[1].strip(), m[2].strip()] for m in mps]
				if options.mass1:
					mps = [mp for mp in mps if int(mp[0])==options.mass1]
				if options.mass2:
					mps = [mp for mp in mps if int(mp[1])==options.mass2]

				## looping over masspoints
				for iiii,mp in enumerate(mps):
					thebasesig = mm.makeCmd([mm.getVariable("expr",""), b, theSc.replace("/","_"), mm.treedirs, mm.getVariable("treename","treeProducerSusyMultilepton"), mccs, macros, mm.getVariable("lumi","12.9"), friends, flagsSig, func.getCut(options, mm.getVariable("firstCut","alwaystrue"), mm.getVariable("expr",""), theBin), ""])
					thebasesys = mm.makeCmd([mm.getVariable("expr",""), b, theSc.replace("/","_"), mm.treedirs, mm.getVariable("treename","treeProducerSusyMultilepton"), mccs, macros, mm.getVariable("lumi","12.9"), friends, flagsSys, func.getCut(options, mm.getVariable("firstCut","alwaystrue"), mm.getVariable("expr",""), theBin), options.postfix])
					thecmd, cardpath = prepareJob(mm, theSc.replace("/", "_")+"_mp_"+deformat(mp[2]), mp, thebasesig, thebasesys, b, bkgDir, myDir, xslist, options)
					mm.registerCmd(thecmd, theSc.replace("/", "_")+"_mp_"+deformat(mp[2]),False,5)
					cards.append(cardpath)
					#break
	mm.runJobs()
	mm.clearJobs()


## do the HTML systematics check page
#if options.doCheck and os.path.isdir(combinePath) and os.path.exists(combinePath+"/HiggsAnalysis/CombinedLimit/test/systematicsAnalyzer.py"):
#	for card in cards:
#		func.mkdir(card[0]+"/html")
#		mm.registerCmd("python "+combinePath+"/HiggsAnalysis/CombinedLimit/test/systematicsAnalyzer.py "+card[0]+"/"+card[1]+" --all -f html > "+card[0]+"/html/SR.card.html", "sysCheck", False, 9999999999, work = combinePath, src = combinePath)
#	mm.runJobs()
#	mm.clearJobs()



