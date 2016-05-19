from CMGTools.TTHAnalysis.treeReAnalyzer import *

def coneptRA5(ptlep,minireliso,ptratio,ptrel,pdgid,wp):
    assert (wp==2)
    if (abs(pdgid)!=11 and abs(pdgid)!=13):
        return ptlep
    A = 0.12 if (abs(pdgid)==11) else 0.16
    B = 0.80 if (abs(pdgid)==11) else 0.76
    C = 7.2 if (abs(pdgid)==11) else 7.2
    if (ptrel>C):
        return ptlep*(1+max(minireliso-A,0))
    else:
        return max(ptlep,ptlep/ptratio*B)


def coneptRA7(ptlep,minireliso,ptratio,ptrel,pdgid,wp):
    assert (wp==2)
    if (abs(pdgid)!=11 and abs(pdgid)!=13):
        return ptlep
    A = 0.16 if (abs(pdgid)==11) else 0.20
    B = 0.76 if (abs(pdgid)==11) else 0.69
    C = 7.2 if (abs(pdgid)==11) else 6.0
    if (ptrel>C):
        return ptlep*(1+max(minireliso-A,0))
    else:
        return max(ptlep,ptlep/ptratio*B)


