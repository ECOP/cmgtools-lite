
from CMGTools.TTHAnalysis.tools.leptonJetReCleaner import LeptonJetReCleaner
from CMGTools.TTHAnalysis.tools.conept import conept_TTH

MODULES=[]

from CMGTools.TTHAnalysis.tools.combinedObjectTaggerForCleaning import *
from CMGTools.TTHAnalysis.tools.fastCombinedObjectRecleaner import *

from CMGTools.TTHAnalysis.tools.functionsTTH import _ttH_idEmu_cuts_E2, clean_and_FO_selection_TTH
from CMGTools.TTHAnalysis.tools.TTVVariables import TTVVariables

def HLTEmul(lep):
    if (abs(lep.pdgId)!=11): return True
    if (lep.hadronicOverEm>=(0.10-0.03*(abs(lep.etaSc)>1.479))): return False
    if (abs(lep.dEtaScTrkIn)>=(0.01-0.002*(abs(lep.etaSc)>1.479))): return False
    if (abs(lep.dPhiScTrkIn)>=(0.04+0.03*(abs(lep.etaSc)>1.479))): return False
    if (lep.eInvMinusPInv<=-0.05): return False
    if (lep.eInvMinusPInv>=(0.01-0.005*(abs(lep.etaSc)>1.479))): return False
    if (lep.sigmaIEtaIEta>=(0.011+0.019*(abs(lep.etaSc)>1.479))): return False
    return True
    
#def tightHLTEmul(lep):


def looseTTVLepton(lep):
    if abs(lep.pdgId)==11:
        if abs(lep.eta)>2.5: return False
        if lep.pt<=10: return False
        if abs(lep.dz)>=0.1: return False
        if abs(lep.dxy)>=0.05: return False
        if lep.sip3d>=4: return False
        if lep.relIso03>=1: return False
        if lep.mvaIdSpring16GP <(-0.76+(0.76-0.52)*(abs(lep.eta)>0.8 and abs(lep.eta)<1.479)+(0.76-0.23)*(abs(lep.eta)>1.479)  ): return False
        if not HLTEmul(lep): return False
    if abs(lep.pdgId)==13:
        if abs(lep.eta)>2.4: return False
        if lep.pt<=10: return False
        if abs(lep.dz)>=0.1: return False
        if abs(lep.dxy)>=0.05: return False
        if lep.sip3d>=4: return False
        if lep.relIso03>=0.7: return False
        if lep.mediumMuonId==0: return False
    return True

def tightTTVLepton(lep):
    if abs(lep.pdgId)==11:
        if abs(lep.eta)>2.5: return False
        if lep.pt<=10: return False
        if abs(lep.dz)>=0.1: return False
        if abs(lep.dxy)>=0.05: return False
        if lep.sip3d>=4: return False
        if lep.relIso03>=(0.0994+(0.0013*(abs(lep.eta)>1.479))): return False
        if lep.mvaIdSpring16GP < (0.837+(-0.837+0.715)*(abs(lep.eta)>0.8 and abs(lep.eta)<1.479)+(-0.837+0.357)*(abs(lep.eta)>1.479)  ): return False
        if not HLTEmul(lep): return False
    if abs(lep.pdgId)==13:
        if abs(lep.eta)>2.4: return False
        if lep.pt<=10: return False
        if abs(lep.dz)>=0.1: return False
        if abs(lep.dxy)>=0.05: return False
        if lep.sip3d>=4: return False
        if lep.relIso03>=0.25: return False
        if lep.mediumMuonId==0: return False
    return True


MODULES.append( ('leptonJetReCleanerTTVCB', lambda : LeptonJetReCleaner("ReclCB",
                                                                      looseLeptonSel = lambda lep : looseTTVLepton(lep),#lep.miniRelIso < 0.4 and lep.sip3d < 8,
                                                                      cleaningLeptonSel = lambda lep : looseTTVLepton(lep),#clean_and_FO_selection_TTH(lep),
                                                                      FOLeptonSel = lambda lep : looseTTVLepton(lep),#clean_and_FO_selection_TTH(lep),
                                                                      tightLeptonSel = lambda lep : tightTTVLepton(lep), #clean_and_FO_selection_TTH(lep) and (abs(lep.pdgId)!=13 or lep.mediumMuonId>0) and lep.mvaTTH > 0.90,
                                                                      cleanJet = lambda lep,jet,dr : dr<0.4,
                                                                      selectJet = lambda jet: abs(jet.eta)<2.4 and jet.pt>30,
# and (abs(jet.eta)<3 or jet.puId)
                                                                      cleanTau = lambda lep,tau,dr: True,#tau.pt > 20 and abs(tau.eta)<2.3 and abs(tau.dxy) < 1000 and abs(tau.dz) < 0.2 and tau.idMVAdR03 >=2  and tau.idDecayMode,
                                                                      looseTau = lambda tau: True,#tau.pt > 20 and abs(tau.eta)<2.3 and abs(tau.dxy) < 1000 and abs(tau.dz) < 0.2 and tau.idMVAdR03 >=2  and tau.idDecayMode,
                                                                      tightTau = lambda tau: True,#tau.idMVAdR03 >= 3,
                                                                      cleanJetsWithTaus = False,
                                                                      cleanTausWithLoose = False,
                                                                      doVetoZ = False,
                                                                      doVetoLMf = False,
                                                                      doVetoLMt = False,
                                                                      jetPt = 30,
                                                                      bJetPt = 30,
                                                                      coneptdef = lambda lep: 0,#conept_TTH(lep),
                                                                      storeJetVariables = True) ))



MODULES.append( ('leptonJetReCleanerTTVMVA', lambda : LeptonJetReCleaner("ReclMVA",
                                                                      looseLeptonSel = lambda lep : lep.miniRelIso < 0.4 and lep.sip3d < 8,
                                                                      cleaningLeptonSel = lambda lep : clean_and_FO_selection_TTH(lep),
                                                                      FOLeptonSel = lambda lep : clean_and_FO_selection_TTH(lep),
                                                                      tightLeptonSel = lambda lep : clean_and_FO_selection_TTH(lep) and (abs(lep.pdgId)!=13 or lep.mediumMuonId>0) and lep.mvaTTH > 0.90,
                                                                      cleanJet = lambda lep,jet,dr : dr<0.4,
                                                                      selectJet = lambda jet: abs(jet.eta)<2.4 and jet.pt>30,
# and (abs(jet.eta)<3 or jet.puId)
                                                                      cleanTau = lambda lep,tau,dr: True,#tau.pt > 20 and abs(tau.eta)<2.3 and abs(tau.dxy) < 1000 and abs(tau.dz) < 0.2 and tau.idMVAdR03 >=2  and tau.idDecayMode,
                                                                      looseTau = lambda tau: True,#tau.pt > 20 and abs(tau.eta)<2.3 and abs(tau.dxy) < 1000 and abs(tau.dz) < 0.2 and tau.idMVAdR03 >=2  and tau.idDecayMode,
                                                                      tightTau = lambda tau: True,#tau.idMVAdR03 >= 3,
                                                                      cleanJetsWithTaus = False,
                                                                      cleanTausWithLoose = False,
                                                                      doVetoZ = False,
                                                                      doVetoLMf = False,
                                                                      doVetoLMt = False,
                                                                      jetPt = 30,
                                                                      bJetPt = 30,
                                                                      coneptdef = lambda lep: conept_TTH(lep),
                                                                      storeJetVariables = True) ))






MODULES.append( ('ttZVariables', lambda: TTVVariables(label="Extra",
                                                      leptonSel=lambda lep: lep.isTight_Recl,
                                                      jetSel=lambda jet:True
                                                      ) ))
