import ROOT
import os, sys

# this script can be used to verify that a friend tree has the same number of entries as the main tree
# (useful to check that all chunks in the friend tree production have been done successfully)
#
# usage: python verifyFTree BIGTREE_DIR FTREE_DIR DATASET_NAME ...

dsets = sys.argv[3:]
if len(sys.argv)<4:
    dsets = [d.replace('evVarFriend_','').replace('.root','') for d in os.listdir(sys.argv[2]) if ('evVarFriend' in d and 'chunk' not in d)]

def openRootOrUrl(myfile):
    _f_t = None
    if os.path.exists(myfile):
        _f_t = ROOT.TFile.Open("dcap://t3se01.psi.ch:22125/"+myfile if "/pnfs/psi.ch" in myfile else myfile)
    elif os.path.exists(myfile+'.url'):
        with open(myfile+'.url','r') as urlf:
            myfile = urlf.readline().replace('\n','')
            if myfile.startswith("root://"):
                _f_t = ROOT.TXNetFile(myfile)
            else:
                _f_t = ROOT.TFile.Open(myfile)
    return _f_t

for dset in dsets:
    print "running " + dset
    f_t = openRootOrUrl(sys.argv[1]+'/'+dset+'/treeProducerSusyMultilepton/tree.root')
    if not f_t: continue
    t_t = f_t.Get("tree")
    if not t_t: continue
    n_t = t_t.GetEntries()
    f_t.Close()
    f_f = openRootOrUrl(sys.argv[2]+'/evVarFriend_'+dset+'.root')
    t_f = f_f.Get("sf/t")
    if not t_f: continue
    n_f = t_f.GetEntries()
    f_f.Close()
    print '%s: %d - %d : %s'%(dset,n_t,n_f,'OK' if n_t==n_f else 'ERROR '*15+' !!!')


