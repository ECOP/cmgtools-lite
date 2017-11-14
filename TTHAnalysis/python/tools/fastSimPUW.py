from CMGTools.TTHAnalysis.treeReAnalyzer import *
import ROOT, os, math, copy
from array import *

if "functions_cc.so" not in ROOT.gSystem.GetLibraries():
    ROOT.gROOT.LoadMacro(os.environ["CMSSW_BASE"]+"/src/CMGTools/TTHAnalysis/python/plotter/functions.cc")

def poissonErrorDn(N):
	alpha = 1-0.6827 # 1 sigma
	L     = 0 if N==0 else ROOT.Math.gamma_quantile(alpha/2, N, 1)
	return N - L

def poissonErrorUp(N):
	alpha = 1-0.6827 # 1 sigma
	U     = ROOT.Math.gamma_quantile_c(alpha/2, N+1, 1)
	return U - N

def eff(p1, f1, p2, f2):
	return (p1/(p1+f1))*((p1+f1+p2+f2)/(p1+p2))

def sigEff(p1, f1, p2, f2):
	down = math.sqrt(math.pow(dEffP1(p1, f1, p2, f2),2)*math.pow(poissonErrorDn(p1),2) + math.pow(dEffF1(p1, f1, p2, f2),2)*math.pow(poissonErrorDn(f1),2) + math.pow(dEffP2(p1, f1, p2, f2),2)*math.pow(poissonErrorDn(p2),2) + math.pow(dEffF2(p1, f1, p2, f2),2)*math.pow(poissonErrorDn(f2),2))
	up   = math.sqrt(math.pow(dEffP1(p1, f1, p2, f2),2)*math.pow(poissonErrorUp(p1),2) + math.pow(dEffF1(p1, f1, p2, f2),2)*math.pow(poissonErrorUp(f1),2) + math.pow(dEffP2(p1, f1, p2, f2),2)*math.pow(poissonErrorUp(p2),2) + math.pow(dEffF2(p1, f1, p2, f2),2)*math.pow(poissonErrorUp(f2),2))
	return down, up

def dEffP1(p1, f1, p2, f2):
	return (f1*p2*(f1+2*p1+p2)+f2*(f1*p2-p1*p1))/(math.pow(p1+f1,2)*math.pow(p1+p2,2))

def dEffF1(p1, f1, p2, f2):
	return -1*p1*(p2+f2)/(math.pow(p1+f1,2)*(p1+p2))

def dEffP2(p1, f1, p2, f2):
	return -1*p1*(f1+f2)/((p1+f1)*math.pow(p1+p2,2))

def dEffF2(p1, f1, p2, f2):
	return p1/((p1+f1)*(p1+p2))


## CutsFile
## ___________________________________________________________________
class CutsFile():
    ## __init__
    ## _______________________________________________________________
	def __init__(self, path):
		f = open(path,"r")
		lines = [line.strip().rstrip("\n").strip() for line in f.readlines()]
		lines = filter(lambda x: x[0:1]!="#" and len(x)>0, lines)
		f.close()
		cuts = ["("+line.split(":")[1]+")" for line in lines]
		self.cutsstring = "&&".join(cuts)
			

## FastSimPUWProducer
## ___________________________________________________________________
class FastSimPUWProducer():


	## __init__
	## _______________________________________________________________
	def __init__(self, dataFile, dataHist, cutsFile):
		if not os.path.exists(cutsFile):
			raise RuntimeError, "Cannot load cuts file"
		self.cuts = CutsFile(cutsFile)
		if not os.path.exists(dataFile) or not dataHist:
			raise RuntimeError, "Cannot find TFile with PU profile in data" 
		f = ROOT.TFile.Open(dataFile,"read")
		self.dataProfile = copy.deepcopy(f.Get(dataHist))
		f.Close()
		if not self.dataProfile:
			raise RuntimeError, "Cannot load PU profile in data"
		int = self.dataProfile.Integral(0, self.dataProfile.GetNbinsX()+1)
		self.dataProfile.Scale(1.0/int)


	## init
	## _______________________________________________________________
	def init(self, tree):
		## preparational stuff before the event loop
		self.getCutValue(tree)
		self.getPUProfile(tree)
		self.getAcceptance(tree)
		self.doTheFit()


	## doTheFit
	## _______________________________________________________________
	def doTheFit(self):
		if self.pA==0 or self.fA==0 or self.pB==0 or self.fB==0:
			self.unc = {"PUuncFS_Up": 1.0, "PUuncFS_Dn": 1.0}
			return
		relEffAB           = eff(self.pA, self.fA, self.pB, self.fB)
		relEffBA           = eff(self.pB, self.fB, self.pA, self.fA)
		sig2ABDn, sig2ABUp = sigEff(self.pA, self.fA, self.pB, self.fB)
		sig2BADn, sig2BAUp = sigEff(self.pB, self.fB, self.pA, self.fA)
		xc                 = [self.meanlow, self.meanhigh]
		xe                 = [self.errlow , self.errhigh ]
		yc                 = [relEffAB    , relEffBA     ]
		yed                = [sig2ABDn    , sig2BADn     ]
		yeu                = [sig2ABUp    , sig2BAUp     ]
		tacc               = ROOT.TGraphAsymmErrors(2, array('d', xc), array('d', yc), array('d', xe), array('d', xe), array('d', yed), array('d', yeu))
		tacc.SetName("puwAcc")
		gfit               = ROOT.TF1("puwFit","pol1",0,50)
		ffit               = tacc.Fit(gfit, "NRQS")
		grap               = ROOT.TGraphErrors(50)
		grap.SetName("puwUnc")
		for i in range(0,50):
			grap.SetPoint(i, i, gfit.Eval(i))
		ffit.GetConfidenceIntervals(grap.GetN(), 1, 1, grap.GetX(), grap.GetEY(), 0.68, False)
		central = 0
		up      = 0
		dn      = 0
		for b in range(1,self.dataProfile.GetNbinsX()+1):
			nvtx     = self.dataProfile.GetBinContent(b)
			gpt      = grap.GetY()[b-1]
			gerr     = grap.GetEY()[b-1]
			central += nvtx*gpt
			up      += nvtx*(gpt+gerr)
			dn      += nvtx*(gpt-gerr)
		self.unc = {"PUuncFS_Up": float(up)/central, "PUuncFS_Dn": float(dn)/central}
		#### debugging
		##f = ROOT.TFile("/mnt/t3nfs01/data01/shome/cheidegg/d/CMSSW_8_0_20/src/CMGTools/TTHAnalysis/python/tools/debug.root","recreate")
		##f.cd()
		##tacc.Write()
		##gfit.Write()
		##grap.Write()
		##f.Close()


	## getAcceptance
	## _______________________________________________________________
	def getAcceptance(self, t):
		self.pA = float(t.GetEntries("nVert<" +str(self.cutvalue)+"&&"+self.cuts.cutsstring)          )
		self.fA = float(t.GetEntries("nVert<" +str(self.cutvalue)                          ) - self.pA)
		self.pB = float(t.GetEntries("nVert>="+str(self.cutvalue)+"&&"+self.cuts.cutsstring)          )
		self.fB = float(t.GetEntries("nVert>="+str(self.cutvalue)                          ) - self.pB)


	## getCutValue
	## _______________________________________________________________
	def getCutValue(self, t):
		self.cuvalue = 20
		h = ROOT.TH1F("cutval","cutval",50,0,50)
		t.Draw("nVert>>cutval")
		theValue = int(round(h.GetMean()))
		self.cutvalue = theValue if 0<theValue<50 else self.cutvalue
		del h


	## getPUProfile
	## _______________________________________________________________
	def getPUProfile(self, t):
		if not hasattr(self, "cutvalue"): self.cutvalue = 20
		h = ROOT.TH1F("cutval","cutval",50,0,50)
		t.Draw("nVert>>cutval","nVert<"+str(self.cutvalue))
		self.meanlow  = h.GetMean()
		self.errlow   = h.GetMeanError()
		h.Reset()
		t.Draw("nVert>>cutval","nVert>="+str(self.cutvalue))
		self.meanhigh = h.GetMean()
		self.errhigh  = h.GetMeanError()
		del h


	## __call__
	## _______________________________________________________________
	def __call__(self, event):
		return self.unc;


	## listBranches
	## _______________________________________________________________
	def listBranches(self):
		biglist = [
		    ("PUuncFS_Up", "F"),
		    ("PUuncFS_Dn", "F")]
		return biglist





if __name__ == '__main__':
    from sys import argv
    file = ROOT.TFile(argv[1])
    tree = file.Get("tree")
    tree.vectorTree = True
    class Tester(Module):
        def __init__(self, name):
            Module.__init__(self,name,None)
            self.sf1 = LeptonChoiceEWK("Old", 
                lambda lep : lep.relIso03 < 0.5, 
                lambda lep : lep.relIso03 < 0.1 and lep.sip3d < 4 and _susy2lss_lepId_CB(lep),
                cleanJet = lambda lep,jet,dr : (lep.pt > 10 and dr < 0.4))
            self.sf2 = LeptonChoiceEWK("PtRel", 
                lambda lep : lep.relIso03 < 0.4 or lep.jetPtRel > 5, 
                lambda lep : (lep.relIso03 < 0.1 or lep.jetPtRel > 14) and lep.sip3d < 4 and _susy2lss_lepId_CB(lep),
                cleanJet = lambda lep,jet,dr : (lep.pt > 10 and dr < 0.4))
            self.sf3 = LeptonChoiceEWK("MiniIso", 
                lambda lep : lep.miniRelIso < 0.4, 
                lambda lep : lep.miniRelIso < 0.05 and lep.sip3d < 4 and _susy2lss_lepId_CB(lep),
                cleanJet = lambda lep,jet,dr : (lep.pt > 10 and dr < 0.4))
            self.sf4 = LeptonChoiceEWK("PtRelJC", 
                lambda lep : lep.relIso03 < 0.4 or lep.jetPtRel > 5, 
                lambda lep : (lep.relIso03 < 0.1 or lep.jetPtRel > 14) and lep.sip3d < 4 and _susy2lss_lepId_CB(lep),
                cleanJet = lambda lep,jet,dr : (lep.pt > 10 and dr < 0.4 and not (lep.jetPtRel > 5 and lep.pt*(1/lep.jetPtRatio-1) > 25)))
            self.sf5 = LeptonChoiceEWK("MiniIsoJC", 
                lambda lep : lep.miniRelIso < 0.4, 
                lambda lep : lep.miniRelIso < 0.05 and lep.sip3d < 4 and _susy2lss_lepId_CB(lep),
                cleanJet = lambda lep,jet,dr : (lep.pt > 10 and dr < 0.4 and not (lep.jetDR > 0.5*10/min(50,max(lep.pt,200)) and lep.pt*(1/lep.jetPtRatio-1) > 25)))
        def analyze(self,ev):
            print "\nrun %6d lumi %4d event %d: leps %d" % (ev.run, ev.lumi, ev.evt, ev.nLepGood)
            print self.sf1(ev)
            print self.sf2(ev)
            print self.sf3(ev)
            print self.sf4(ev)
            print self.sf5(ev)
    el = EventLoop([ Tester("tester") ])
    el.loop([tree], maxEvents = 50)


