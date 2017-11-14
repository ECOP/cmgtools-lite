import os, copy, ROOT

def deformatFile(rawFile):
	if "#" in rawFile:
		theList = []
		for line in rawFile.split("#"):
			sl = line.split("*")
			theList.append((sl[1],float(sl[0])))
		return theList
	return [(rawFile,1.0)]

## DO NOT TOUCH: THIS IS GOING TO BE CHANGED BY THE BATCH SUBMISSION
name     = "THENAME"
sig      = "THESIGNAL"
mass1    = "THEMASS1"
mass2    = "THEMASS2"
offset   = int("THEOFFSET")
treedir  = "THETREEDIR"
treename = "THETREENAME"
filelist = deformatFile("THEFILE")
xs       = THEXS
#q2file   = "THEQ2FILE"
#q2syntax = "THEQ2SYNTAX"
wstr     = "THEWEIGHTSTR"
frfiles  = THEFRFILES
thejec   = "THEJEC" # name of the jec in the systs file
themet   = "THEMET" # name of the met in the systs file
q2acc    = "THEQ2ACC" # name of the q2acc in the systs file
#puw      = ""   # name of the puw in the systs file
frjec    = THEFRJEC # central, jecUp, jecDn
wvjec    = THEWVJEC # central, jecUp, jecDn
frmet    = THEFRMET # pfMET, genMET
wvmet    = THEWVMET # pfMET, genMET
wVars    = THEWEIGHTVARS
bkgdir   = "THEBKGDIR"
mcadir   = "THEMCADIR"
outdir   = "THEOUTDIR"
themca   = "THEMCA"
thecuts  = "THECUTS"
thedummy = "THEDUMMY"
thesyst  = "THESYST"
#first    = "THECMDFIRST"
#second   = "THECMDSECOND"
sigbase  = "THEBASESIG" # without certain flags that affect the systematics
sysbase  = "THEBASESYS" # the final command including all systematic variations



## ---------

def cmd(cmd):
	print cmd
	os.system(cmd)

def cp(location, destination):
	cmd("cp " + location + " " + destination)

def mkdir(path):
	if os.path.isdir(path): return
	cmd("mkdir -p " + path)

def makeFakeRate(thefrfiles, more = [], index = 0):
	myfrfiles = copy.deepcopy(thefrfiles)
	if len(more) > 0 and index > 0:
		if more[0]:
			myfrfiles = [x if x != more[0] else more[index] for x in myfrfiles]
		else:
			myfrfiles.append(more[index])
	myfrfiles = filter(None, myfrfiles)
	if len(myfrfiles) == 0: return ""
	return ", FakeRate=\"" + "\,".join(myfrfiles) + "\""

def div(num, den):
	if den==0: return 1
	return float(num)/den

def cutOff(value, up, down):
	if value > up  : return up
	if value < down: return down
	return value

def openTree(theFile):
	global treedir, treename
	fpath = None
	for base in treedir.split(";"):
		if   os.path.exists(base+"/"+theFile+"/"+treename+"/tree.root"):
			fpath = base+"/"+theFile+"/"+treename+"/tree.root"
			break
		elif os.path.exists(base+"/"+theFile+"/"+treename+"/tree.root.url"):
			fpath = open(base+"/"+theFile+"/"+treename+"/tree.root.url","r").readlines()[0].rstrip("\n")
			break
		else:
			continue
	if not fpath: return None,None
	f = ROOT.TFile.Open(fpath, "read")
	if not f: return None,None
	t = f.Get("tree")
	if not t: return None,None
	if t.GetEntries() == 0: return None,None
	return f,t

def openSkimReports():
	global filelist, treedir, treename
	tote=0; evts = []
	for file in filelist:
		for theFile in file[0].split("+"):
			for base in treedir.split(";"):
				if not os.path.exists(base+"/"+theFile+"/skimAnalyzerCount/SkimReport.txt"): continue
				f = open(base+"/"+theFile+"/skimAnalyzerCount/SkimReport.txt","r")
				for line in f.readlines():
					sl = line.split()
					if "All Events" in line:
						evts.append(file[1]*float(sl[2]))
						tote+=evts[-1]
						break
				f.close()
				break
	return [e/tote for e in evts]	

def getWsum(cutvalue=20):
	global filelist
	print "retrieving wsums"
	evts  = openSkimReports()
	wsums = [0,0,0,0,0]
	for file in filelist:
		for iff,ff in enumerate(file[0].split("+")):
			f,t = openTree(ff)
			if not t: continue
			wsums[0] += file[1]*evts[iff]*t.GetEntries()
			wsums[1] += file[1]*evts[iff]*t.GetEntries("nVert>="+str(cutvalue))
			wsums[2] += file[1]*evts[iff]*t.GetEntries("nVert<" +str(cutvalue))
			for evt in t:
				wsums[3] += file[1]*evts[iff]*evt.LHEweight_wgt[4]/evt.LHEweight_wgt[0]
				wsums[4] += file[1]*evts[iff]*evt.LHEweight_wgt[8]/evt.LHEweight_wgt[0]
			f.Close()
	return wsums
	#return [total, upP, dnP, upW, dnW]


def doReweightingVariations(infile, outfile, procName, total = 0, up = 0, dn = 0):
	global sig
	f   = ROOT.TFile.Open(infile, "read")
	hcn = f.Get("x_sig_{s}_central".format(s=sig))
	hup = f.Get("x_sig_{s}_up"     .format(s=sig))
	hdn = f.Get("x_sig_{s}_down"   .format(s=sig))
	ff = ROOT.TFile.Open(outfile, "recreate")
	ff.cd()
	hOutUp = ROOT.TH1F(procName+"_Up",procName+"_Up",hcn.GetNbinsX(),hcn.GetXaxis().GetXmin(),hcn.GetXaxis().GetXmax())
	hOutDn = ROOT.TH1F(procName+"_Dn",procName+"_Dn",hcn.GetNbinsX(),hcn.GetXaxis().GetXmin(),hcn.GetXaxis().GetXmax())
	fUp = total/up if up>0 else 1.0
	fDn = total/dn if dn>0 else 1.0
	## loop over SR
	for sr in range(1,hcn.GetNbinsX()+1):
		valZero = float(hcn.GetBinContent(sr)    )
		errZero = float(hcn.GetBinError  (sr)    )
		valUp   = float(hup.GetBinContent(sr)*fUp)
		errUp   = float(hup.GetBinError  (sr)*fUp)
		valDn   = float(hdn.GetBinContent(sr)*fDn)
		errDn   = float(hdn.GetBinError  (sr)*fDn)
		if valUp==0 or valDn==0 or valZero==0: continue
		sysUp  = valUp/valZero - 1.0;
		sysDn  = 1.0 - valZero/valDn
		#RMS    = sysUp*sysUp + sysDn*sysDn
		hOutUp.SetBinContent(sr, cutOff(1+sysUp if sysUp>0 else -1/(sysUp-1), 2, 0.5))
		hOutDn.SetBinContent(sr, cutOff(1+sysDn if sysDn>0 else -1/(sysDn-1), 2, 0.5))	
	hOutUp.Write()	
	hOutDn.Write()	
	hOutUp.SetDirectory(0)
	hOutDn.SetDirectory(0)
	ff.Close()
	return [hOutUp, hOutDn]

##def doPuwVariation(infile, outfile, wsum):
##	## extracts relative correction for pileup according to susy recommendation
##	global sig, puw
##	print "running pileup variation"
##	return doReweightingVariations(infile, outfile, "x_sig_{s}_{j}_Up".format(s=sig, j=puw), wsum[0], wsum[1], wsum[2])


def doQ2Variation(infile, outfile, wsum):
	## extracts relative correction for Q2 Acceptance
	global sig, q2acc
	print "running q2 acceptance variation"
	return doReweightingVariations(infile, outfile, "x_sig_{s}_{j}_Up".format(s=sig, j=q2acc), wsum[3], wsum[4])


def doMetVariation(infile, outfile, sig, jec, met, wVars, q2vars = []):
	global mass1, mass2, q2acc
	print "running met variation"
	## retrieve all histograms
	f = ROOT.TFile.Open(infile, "read")
	pfMET  = f.Get("x_sig_{s}_pfMET" .format(s=sig))
	genMET = f.Get("x_sig_{s}_genMET".format(s=sig))
	jecUp  = f.Get("x_sig_{s}_{j}_Up".format(s=sig, j=jec))
	jecDn  = f.Get("x_sig_{s}_{j}_Dn".format(s=sig, j=jec))
	wvHist = {}
	for key, vals in wVars.iteritems():
		wvHist[key + "Up"] = f.Get("x_sig_{s}_{j}_Up".format(s=sig, j=key))
		wvHist[key + "Dn"] = f.Get("x_sig_{s}_{j}_Dn".format(s=sig, j=key))
	final = pfMET.Clone("x_sig_{s}"       .format(s=sig       )); final.Reset()
	metUp = pfMET.Clone("x_sig_{s}_{m}_Up".format(s=sig, m=met)); metUp.Reset()
	metDn = pfMET.Clone("x_sig_{s}_{m}_Dn".format(s=sig, m=met)); metDn.Reset()
	#if len(puwvars)>0:
	#	puwUp = pfMET.Clone("x_sig_{s}_{q}_Up".format(s=sig, q=puw)); puwUp.Reset()
	#	puwDn = pfMET.Clone("x_sig_{s}_{q}_Dn".format(s=sig, q=puw)); puwDn.Reset()
	if len(q2vars)>0:
		q2Up = pfMET.Clone("x_sig_{s}_{q}_Up".format(s=sig, q=q2acc)); q2Up.Reset()
		q2Dn = pfMET.Clone("x_sig_{s}_{q}_Dn".format(s=sig, q=q2acc)); q2Dn.Reset()
	## do the correction
	for bin in range(1,final.GetNbinsX()+1):
		diff = abs(pfMET.GetBinContent(bin)-genMET.GetBinContent(bin))/2
		avg  = float(pfMET.GetBinContent(bin)+genMET.GetBinContent(bin))/2
		sf   = div(avg, pfMET.GetBinContent(bin))
		final.SetBinContent(bin, avg     )
		final.SetBinError  (bin, diff    )
		metUp.SetBinContent(bin, avg+diff)
		metDn.SetBinContent(bin, avg-diff)
		jecUp.SetBinContent(bin, jecUp   .GetBinContent(bin)*sf)
		jecDn.SetBinContent(bin, jecDn   .GetBinContent(bin)*sf)
		#if len(puwvars)>0:
		#	puwUp.SetBinContent(bin, avg * puwvars[0].GetBinContent(bin))
		#	puwDn.SetBinContent(bin, avg * puwvars[1].GetBinContent(bin))
		if len(q2vars)>0:
			q2Up.SetBinContent(bin, avg * q2vars[0].GetBinContent(bin))
			q2Dn.SetBinContent(bin, avg * q2vars[1].GetBinContent(bin))
			#q2Up.SetBinContent(bin, avg * q2vars[0].GetBinContent(bin)/100.*sf)
			#q2Dn.SetBinContent(bin, avg * q2vars[1].GetBinContent(bin)/100.*sf)
		for key, vals in wVars.iteritems():
			wvHist[key + "Up"].SetBinContent(bin, wvHist[key + "Up"].GetBinContent(bin)*sf)
			wvHist[key + "Dn"].SetBinContent(bin, wvHist[key + "Dn"].GetBinContent(bin)*sf)
	## write all to output file
	ff = ROOT.TFile.Open(outfile, "recreate")
	ff.cd()
	final.Write()
	metUp.Write()
	metDn.Write()
	jecUp.Write()
	jecDn.Write()
	#if len(puwvars)>0:
	#	puwUp.Write()
	#	puwDn.Write()
	if len(q2vars)>0:
		q2Up.Write()
		q2Dn.Write()
	for key, vals in wVars.iteritems():
		wvHist[key + "Up"].Write()
		wvHist[key + "Dn"].Write()
	ff.Close()

def plugFiles(dirs):
	return " ".join(["--infile {d}".format(d=thedir) for thedir in dirs])

def makeWeight(wstr, wvar):
	if not wvar: return wstr
	if wvar == "1": return wstr
	return "("+wstr+")*"+wvar

cmdbase = sigbase.replace("{","[[").replace("}","]]")
cmdbase = cmdbase.replace("[[[","{").replace("]]]","}")
sysbase = sysbase.replace("{","[[").replace("}","]]")
sysbase = sysbase.replace("[[[","{").replace("]]]","}")
#cmdbase = "python {sc} {{MCA}} {FIRST} {{SYSTS}} --od {{OUTDIR}} ".format(sc=script, FIRST=first)
mcabase = "\n".join("sig_{{name}} : {file} : {xs} : {{ws}} ; Label=\"{{name}}\", isFastSim{{FRfiles}}".format(file=theFile[0], xs="%s*%f"%(xs,theFile[1])) for theFile in filelist)

short = mass1 + "_" + mass2
#puwdir = outdir + "/puw/" + short
q2dir  = outdir + "/q2/"  + short
accdir = outdir + "/acc/" + short
mpsdir = outdir + "/mps/" + short
mkdir(mcadir)
mkdir(outdir)
#mkdir(puwdir)
mkdir(q2dir )
mkdir(accdir)
mkdir(mpsdir)

wsum   = getWsum()

## first do pileup variation (simple)
##if puw:
##	f = open(mcadir + "/mca_puw_"+name+".txt", "w")
##	f.write(mcabase.format(name=sig+"_central", ws=wstr                        , FRfiles=makeFakeRate(frfiles)) + "\n")
##	f.write(mcabase.format(name=sig+"_up"     , ws=makeWeight(wstr, "nVert>20"), FRfiles=makeFakeRate(frfiles)) + "\n")
##	f.write(mcabase.format(name=sig+"_down"   , ws=makeWeight(wstr, "nVert<20"), FRfiles=makeFakeRate(frfiles)) + "\n")
##	f.close()
##	mybase = cmdbase.format(MCA=mcadir + "/mca_puw_"+name+".txt", SYS="", O=outdir + "/puw/"+short)
##	cmd(mybase.replace("[[","{").replace("]]","}") + " --asimov")
##
##	## get central value of acceptance 
##	puwvars = doPuwVariation(puwdir+"/common/SR.input.root", puwdir+"/puw_SR.input.root", wsum)


## then do Q2 acceptance variation
q2vars = []
if q2acc:
	f = open(mcadir + "/mca_q2_"+name+".txt", "w")
	f.write(mcabase.format(name=sig+"_central", ws=makeWeight(wstr, "LHEweight_wgt[0]/LHEweight_wgt[0]"), FRfiles=makeFakeRate(frfiles)) + "\n")
	f.write(mcabase.format(name=sig+"_up"     , ws=makeWeight(wstr, "LHEweight_wgt[4]/LHEweight_wgt[0]"), FRfiles=makeFakeRate(frfiles)) + "\n")
	f.write(mcabase.format(name=sig+"_down"   , ws=makeWeight(wstr, "LHEweight_wgt[8]/LHEweight_wgt[0]"), FRfiles=makeFakeRate(frfiles)) + "\n")
	f.close()
	mybase = cmdbase.format(MCA=mcadir + "/mca_q2_"+name+".txt", CUTS=thedummy, SYS="", O=outdir + "/q2/"+short)
	cmd(mybase.replace("[[","{").replace("]]","}") + " --asimov")

	## get central value of acceptance 
	q2vars = doQ2Variation(q2dir+"/common/SR.input.root", q2dir+"/q2_SR.input.root", wsum)
	


## then loop on pfMET and genMET, and the JEC
if len(frmet)==2:
	f = open(mcadir + "/mca_acc_"+name+".txt", "w")
	f.write(mcabase.format(name=sig+"_pfMET"           , ws=wstr                      , FRfiles=makeFakeRate(frfiles,frmet, 0)) + "\n")
	f.write(mcabase.format(name=sig+"_genMET"          , ws=makeWeight(wstr, wvmet[1]), FRfiles=makeFakeRate(frfiles,frmet, 1)) + "\n")
	if len(frjec)==3:
		f.write(mcabase.format(name=sig+"_"+thejec+"_Up"   , ws=makeWeight(wstr, wvjec[1]), FRfiles=makeFakeRate(frfiles,frjec, 1)) + "\n")
		f.write(mcabase.format(name=sig+"_"+thejec+"_Dn"   , ws=makeWeight(wstr, wvjec[2]), FRfiles=makeFakeRate(frfiles,frjec, 2)) + "\n")
	for k,vals in wVars.iteritems():
		first  = vals[0].replace("Sqrt","TMath::Sqrt").replace("Power","TMath::Power")
		second = vals[1].replace("Sqrt","TMath::Sqrt").replace("Power","TMath::Power")
		f.write(mcabase.format(name=sig+"_"+k+"_Up", ws=makeWeight(wstr,first), FRfiles=makeFakeRate(frfiles)) + "\n")
		f.write(mcabase.format(name=sig+"_"+k+"_Dn", ws=makeWeight(wstr,second), FRfiles=makeFakeRate(frfiles)) + "\n")
	f.close()
	mybase = cmdbase.format(MCA=mcadir + "/mca_acc_"+name+".txt", CUTS=thecuts, SYS="", O=outdir + "/acc/"+short)
	cmd(mybase.replace("[[","{").replace("]]","}") + " --asimov")
	
	
	## get central value of acceptance 
	doMetVariation(accdir + "/common/SR.input.root", accdir + "/acc_SR.input.root", sig, thejec, themet, wVars, q2vars)


## prepare the proper job
cp(themca, mcadir + "/mca_full_"+name+".txt")
f = open(mcadir + "/mca_full_"+name+".txt", "a")
f.write(mcabase.format(name=sig+"+"              , ws=wstr, FRfiles=makeFakeRate(frfiles, frjec, 0)) + "\n")
if len(frjec)==3:
	f.write(mcabase.format(name=sig+"_"+thejec+"_Up+", ws=wstr, FRfiles=makeFakeRate(frfiles, frjec, 1)).replace("\n", ",SkipMe=True\n") + ", SkipMe=True\n")
	f.write(mcabase.format(name=sig+"_"+thejec+"_Dn+", ws=wstr, FRfiles=makeFakeRate(frfiles, frjec, 2)).replace("\n", ",SkipMe=True\n") + ", SkipMe=True\n")
if len(frmet)==2:
	f.write(mcabase.format(name=sig+"_"+themet+"_Up+", ws=wstr, FRfiles=makeFakeRate(frfiles, frmet, 0)).replace("\n", ",SkipMe=True\n") + ", SkipMe=True\n")
	f.write(mcabase.format(name=sig+"_"+themet+"_Dn+", ws=wstr, FRfiles=makeFakeRate(frfiles, frmet, 0)).replace("\n", ",SkipMe=True\n") + ", SkipMe=True\n")
for k,vals in wVars.iteritems():
	f.write(mcabase.format(name=sig+"_"+k+"_Up+", ws=makeWeight(wstr,vals[0]), FRfiles=makeFakeRate(frfiles)).replace("\n", ",SkipMe=True\n") + ", SkipMe=True\n")
	f.write(mcabase.format(name=sig+"_"+k+"_Dn+", ws=makeWeight(wstr,vals[1]), FRfiles=makeFakeRate(frfiles)).replace("\n", ",SkipMe=True\n") + ", SkipMe=True\n")
if len(q2vars)>0:
	f.write(mcabase.format(name=sig+"_"+q2acc+"_Up+", ws=wstr, FRfiles=makeFakeRate(frfiles)).replace("\n", ",SkipMe=True\n") + ", SkipMe=True\n")
	f.write(mcabase.format(name=sig+"_"+q2acc+"_Dn+", ws=wstr, FRfiles=makeFakeRate(frfiles)).replace("\n", ",SkipMe=True\n") + ", SkipMe=True\n")
#if len(puwvars)>0:
#	f.write(mcabase.format(name=sig+"_"+puw  +"_Up+", ws=wstr, FRfiles=makeFakeRate(frfiles)) + ",SkipMe=True\n")
#	f.write(mcabase.format(name=sig+"_"+puw  +"_Dn+", ws=wstr, FRfiles=makeFakeRate(frfiles)) + ",SkipMe=True\n")
f.close()

## run the proper job, which is actually just making the cards
mybase = sysbase.format(MCA=mcadir + "/mca_full_"+name+".txt", CUTS=thecuts, SYS=thesyst, O=outdir+"/mps/"+short)
cmd(mybase.replace("[[","{").replace("]]","}") + " --ip x " + plugFiles([bkgdir+"/common/SR.input.root", accdir+"/acc_SR.input.root"]))


