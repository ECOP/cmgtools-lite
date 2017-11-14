import datetime, sys, os
from optparse import OptionParser
from lib import maker
from lib import functions as func

def makeDummyMCAs(mm):
	if mm.options.allIn: samples = collectInputSamples(mm)
	else               : samples = collectMcaSamples  (mm)
	## noSplit
	if mm.options.noSplit:
		thePath = mm.tmppath+"/mca.txt"
		makeMCA(thePath, samples)
		return [thePath]
	## split
	mcas = []
	for i,s in enumerate(samples):
		thePath = mm.tmppath+"/mca"+str(i)+".txt"
		makeMCA(thePath, [s])
		mcas.append(thePath)
	return mcas

def makeMCA(path, samples):
	if len(samples)==0: return
	f = open(path,"w")
	for p in samples:
		f.write(p +" : "+ p +" : 1.0\n")
		#f.write(p +" : "+ p +" : 1.0 ; FakeRate=\"susy-ewkino/3l/fakerate_standard_central.txt\,susy-ewkino/3l/fakerate_load_isprompt.txt\", isFastSim\n")
	f.close()

def collectMcaSamples(mm):
	s = mm.options.samples if len(mm.options.samples)>0 and not mm.options.allSamples else mm.getAllSamples(m.options.allSamples)
	return s

def collectInputSamples(mm):
	samples = []	
	for d in os.listdir(mm.treedir):
		## only consider real samples
		if not os.path.isdir(mm.treedir +"/"+ d): continue
		if not os.path.exists(mm.treedir +"/"+ d +"/"+mm.getVariable("treename","treeProducerSusyMultilepton")+"/tree.root") and \
           not os.path.exists(mm.treedir +"/"+ d +"/"+mm.getVariable("treename","treeProducerSusyMultilepton")+"/tree.root.url"): continue
		## exclude or accept
		if mm.options.accept  != [] and all([d.find(a) == -1 for a in mm.options.accept ]): continue
		if mm.options.exclude != [] and any([d.find(e) >  -1 for e in mm.options.exclude]): continue
		## if exists
		if not mm.options.force and os.path.isdir(mm.outdir +"/"+ d): continue
		samples.append(d)
	return samples

parser = OptionParser(usage="%prog cfg regions treedir outdir [options]")
parser = maker.addMakerOptions(parser)
parser.add_option("--json"      , dest="json"      , type="string", default=None, help="Regexp of the JSON file you want to use to skim your tree")
parser.add_option("--samples"   , dest="samples"   , type="string", action="append", default=[], help="Only use these samples in the skimming and ignore procs, bkgs and sigs")
parser.add_option("--allSamples", dest="allSamples", action="store_true", default=False, help="Skim all samples in the MCA")
parser.add_option("--allIn"     , dest="allIn"     , action="store_true", default=False, help="Run all samples from input tree directory and not from MCA")
parser.add_option("--exclude"   , dest="exclude"   , type="string", action="append", default=[], help="Semicolon-separated list of samples to exclude (regexp)")
parser.add_option("--accept"    , dest="accept"    , type="string", action="append", default=[], help="Semicolon-separated list of samples to accept (regexp)")
parser.add_option("-F"          , dest="force"     , action="store_true", default=False, help="Force reproduction of skims that already were done.")
parser.add_option("--noSplit"   , dest="noSplit"   , action="store_true", default=False, help="Do not make jobs to run in parallel, only make one big job")

base = "python skimTrees.py {MCA} {CUTS} {O} -P {T} --tree {TREENAME} {MCCS} {MACROS} {FRIENDS} {JSON} {FLAGS}"
baseFriends = "python skimFTrees.py {T} {F} {O}" 

(options, args) = parser.parse_args()
options         = maker.splitLists(options)
options.samples = func.splitList(options.samples)
options.accept  = func.splitList(options.accept )
options.exclude = func.splitList(options.exclude)
mm              = maker.Maker("skimmaker", base, args, options, parser.defaults)

## skim main tree
friends = mm.collectFriends()	
flags   = " ".join(mm.options.flags) ## we do not want to have all flags, only the additional ones given here

for r in range(len(mm.regions)):
	mm.iterateRegion()
	mm.reloadBase(base)

	mccs   = mm.collectMCCs  ()
	macros = mm.collectMacros()

	output = mm.outdir
	func.mkdir(output, False)
	json   = options.json if options.json else ""
	mcas   = makeDummyMCAs(mm)
	
	for mca in mcas:
		mm.submit([mca, mm.getVariable("cutfile",""), output, mm.treedir, mm.getVariable("treename", "treeProducerSusyMultilepton"), mccs, macros, friends, json, flags],"main_"+mca[mca.rfind("/")+1:].rstrip(".txt"),False)

	mm.runJobs()
	mm.clearJobs()


	#### skim friend trees
	mm.reloadBase(baseFriends)
	fs, fm = mm.getFriendLocations()
	
	for i,f in enumerate(fs):
		#output = mm.outdir+"/"+fm[i]
		func.mkdir(output)	
		mm.submit([mm.outdir, f, mm.outdir],"friend_"+f[f.rfind("/")+1:],False)

	mm.runJobs()
	mm.clearJobs()


