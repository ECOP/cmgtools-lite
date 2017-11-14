#include "TH2F.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TSystem.h"
#include <iostream>

using namespace std;

TString CMSSW_BASE_SF = gSystem->ExpandPathName("${CMSSW_BASE}");
TString DATA_SF = CMSSW_BASE_SF+"/src/CMGTools/TTHAnalysis/data";

float getSF(TH2F* hist, float pt, float eta){
    int xbin = std::max(1, std::min(hist->GetNbinsX(), hist->GetXaxis()->FindBin(pt)));
    int ybin = std::max(1, std::min(hist->GetNbinsY(), hist->GetYaxis()->FindBin(eta)));
    return hist->GetBinContent(xbin,ybin);
}

float getUnc(TH2F* hist, float pt, float eta){
    int xbin = std::max(1, std::min(hist->GetNbinsX(), hist->GetXaxis()->FindBin(pt)));
    int ybin = std::max(1, std::min(hist->GetNbinsY(), hist->GetYaxis()->FindBin(eta)));
    return hist->GetBinError(xbin,ybin);
}


// TRIGGER SCALE FACTORS FULLSIM
// -------------------------------------------------------------

// files
TFile* f_trigEff_0tau = new TFile(DATA_SF+"/triggerSF/triggerSF_0tau_EWKino_fastsim_M17_36p5fb.root", "read");
TFile* f_trigEff_1tau = new TFile(DATA_SF+"/triggerSF/triggerSF_1tau_EWKino_fastsim_M17_36p5fb.root", "read");
TFile* f_trigSF_1lep  = new TFile(DATA_SF+"/triggerSF/triggerSF_1lep_EWKino_M17_36p5fb.root"        , "read");

// fullsim
TH2F* h_trigSF_1el   = (TH2F*) f_trigSF_1lep -> Get("ele27_sf");
TH2F* h_trigSF_1mu   = (TH2F*) f_trigSF_1lep -> Get("isomu24_sf");

float triggerSFfullsim(float pt1  , float eta1  , int pdg1  ,
                       float pt2  , float eta2  , int pdg2  ,
                       float pt3=0, float eta3=0, int pdg3=0) {

    int nTaus = (abs(pdg1)==15)+(abs(pdg2)==15)+(abs(pdg3)==15);

    // 2lss or 3l 1tau region
    if(pdg3==0 || nTaus==1) {
        float lpt1 = pt1; float lpt2 = pt2; int lpdg1 = pdg1; int lpdg2 = pdg2;
        if(abs(pdg1)==15) {lpt1 = pt2; lpdg1 = pdg2; lpt2 = pt3; lpdg2 = pdg3; } 
        if(abs(pdg2)==15) {lpt1 = pt1; lpdg1 = pdg1; lpt2 = pt3; lpdg2 = pdg3; } 
        if(abs(lpdg1)+abs(lpdg2)!=22) return 1.0;
        if(lpt1>30 || lpt2>20       ) return 1.0;
        return 0.95;
    }

    // 3l 2tau
    if(nTaus==2){ 
        float lpt = pt1; float leta = eta1; int lpdg = pdg1;
        if(abs(pdg2)!=15) {lpt = pt2; leta = eta2; lpdg = pdg2; } 
        if(abs(pdg3)!=15) {lpt = pt3; leta = eta3; lpdg = pdg3; } 
        TH2F* hist = (abs(lpdg)==13)?h_trigSF_1mu:h_trigSF_1el;
        return getSF(hist, lpt, abs(leta));
    }

    return 1.0;
}

// fastsim
TH2F* h_trigEff_0tau = (TH2F*) f_trigEff_0tau -> Get("h_0tau_num");
TH2F* h_trigEff_1tau = (TH2F*) f_trigEff_1tau -> Get("h_1tau_num");
TH2F* h_trigEff_1el  = (TH2F*) f_trigSF_1lep  -> Get("hist2dnum_Ele27_WPTight_Gsf_fromemu__HLT_Ele27_WPTight_Gsfpt");
TH2F* h_trigEff_1mu  = (TH2F*) f_trigSF_1lep  -> Get("hist2dnum_IsoMu24orIsoTkMu24_fromemu__HLT_IsoMu24orIsoTkMu24pt");


float triggerSFfastsim(float pt1  , float eta1  , int pdg1  ,
                       float pt2  , float eta2  , int pdg2  ,
                       float pt3=0, float eta3=0, int pdg3=0) {

    int nTaus = (abs(pdg1)==15)+(abs(pdg2)==15)+(abs(pdg3)==15);

    // 2lss or 3l 1tau region
    if(pdg3==0 || nTaus==1) {
        float lpt1 = pt1; float lpt2 = pt2;
        if(abs(pdg1)==15) {lpt1 = pt2; lpt2 = pt3; } 
        if(abs(pdg2)==15) {lpt1 = pt1; lpt2 = pt3; } 
        return getSF(h_trigEff_1tau, lpt1, lpt2);
    }

    // 3l 2tau
    if(nTaus==2){ 
        float lpt = pt1; float leta = eta1; int lpdg = pdg1;
        if(abs(pdg2)!=15) {lpt = pt2; leta = eta2; lpdg = pdg2; } 
        if(abs(pdg3)!=15) {lpt = pt3; leta = eta3; lpdg = pdg3; } 
        TH2F* hist = (abs(lpdg)==13)?h_trigEff_1mu:h_trigEff_1el;
        return getSF(hist, lpt, abs(leta));
    }

    // 3l 0tau
    if(nTaus==0){
        return getSF(h_trigEff_0tau, pt2, pt3);
    } 
    return 1.0;
}


//TFile* f_trigSF       = new TFile(DATA_SF+"/triggerSF/triggerSF_EWKino_fullsim_ICHEP2016_9p2fb.root"       , "read");
//TFile* f_trigSF_ele27 = new TFile(DATA_SF+"/triggerSF/triggerSF_Ele27_EWKino_fullsim_ICHEP2016_12p9fb.root", "read");
//
//TH2F* h_trigSF_3l_mu = (TH2F*) f_trigSF      ->Get("eff_3l_mu" );
//TH2F* h_trigSF_3l_el = (TH2F*) f_trigSF      ->Get("eff_3l_ele");
//TH2F* h_trigSF_2l_mu = (TH2F*) f_trigSF      ->Get("eff_2l_mu" );
//TH2F* h_trigSF_2l_el = (TH2F*) f_trigSF      ->Get("eff_2l_ele");
//TH2F* h_trigSF_ele27 = (TH2F*) f_trigSF_ele27->Get("hist2dnum_Ele27_WPLoose_Gsf__HLT_Ele27_WPLoose_Gsf");
//
//float triggerSFBR6(float pt1, float eta1, int pdg1,
//                   float pt2, float eta2, int pdg2,
//                   float pt3, float eta3, int pdg3) {
//
//        if(abs(pdg1)+abs(pdg2)+abs(pdg3)==43) return 0.86;
//
//        float pt = pt1; float eta = eta1;
//        if(abs(pdg2)==11) {pt = pt2; eta=eta2; }
//        if(abs(pdg3)==11) {pt = pt3; eta=eta3; }
//        TH2F* hist = h_trigSF_ele27;
//        int xbin = std::max(1, std::min(hist->GetNbinsX(), hist->GetXaxis()->FindBin(pt)));
//        int ybin = std::max(1, std::min(hist->GetNbinsY(), hist->GetYaxis()->FindBin(abs(eta))));
//        return hist->GetBinContent(xbin,ybin);
//}
//
//float triggerSF(int BR, float pt1, int pdg1, 
//                        float pt2, int pdg2, 
//                        float pt3 = 0, int pdg3 = 0, 
//                        float pt4 = 0, int pdg4 = 0){
//    // Lesya's mail:
//    // - split for trailing ele or trailing mu
//    // - 3l: subleading vs trailing lepton pt (1l + 2l triggers)
//    // - 2l: leading light lepton vs subleading light lepton ==> good for both 2l+tau and 2lSS cases (1l + 2l triggers)
//    // - l+tautau: use flat 86% everywhere; pt_e > 35 GeV; pt_mu > 25 GeV (1l + l/tau triggers)
//
//    // 3l: 2tau (flat 86% in dedicated function)
//    if(BR == 6) return 1.0;
//
//    // 2lss 
//    if(BR == -1){
//        TH2F* hist = (pdg2 == 13)?h_trigSF_2l_mu:h_trigSF_2l_el;
//        int xbin = std::max(1, std::min(hist->GetNbinsX(), hist->GetXaxis()->FindBin(pt1)));
//        int ybin = std::max(1, std::min(hist->GetNbinsY(), hist->GetYaxis()->FindBin(pt2)));
//        return hist->GetBinContent(xbin,ybin);
//    }
//
//    // 3l: 3light
//    if(BR <= 2) {
//        TH2F* hist = (abs(pdg3) == 13)?h_trigSF_3l_mu:h_trigSF_3l_el;
//        int xbin = std::max(1, std::min(hist->GetNbinsX(), hist->GetXaxis()->FindBin(pt2)));
//        int ybin = std::max(1, std::min(hist->GetNbinsY(), hist->GetYaxis()->FindBin(pt3)));
//        return hist->GetBinContent(xbin,ybin);
//    } 
//
//    // 3l: 2light + 1tau
//    if(BR >= 3 && BR <= 5){
//        vector<int> pdgs; vector<float> pts;
//        if(abs(pdg1)!=15) { pdgs.push_back(abs(pdg1)); pts.push_back(pt1); }
//        if(abs(pdg2)!=15) { pdgs.push_back(abs(pdg2)); pts.push_back(pt2); }
//        if(abs(pdg3)!=15) { pdgs.push_back(abs(pdg3)); pts.push_back(pt3); }
//        TH2F* hist = (pdgs[1] == 13)?h_trigSF_2l_mu:h_trigSF_2l_el;
//        int xbin = std::max(1, std::min(hist->GetNbinsX(), hist->GetXaxis()->FindBin(pts[0])));
//        int ybin = std::max(1, std::min(hist->GetNbinsY(), hist->GetYaxis()->FindBin(pts[1])));
//        return hist->GetBinContent(xbin,ybin);
//    }
//
//    // others: (4l, crwz) 
//    return 1;
//}



// LEPTON SCALE FACTORS FULLSIM
// -------------------------------------------------------------


// electrons
TFile* f_elSF_id   = new TFile(DATA_SF+"/leptonSF/electronSF_id_EWKino_fullsim_M17_36p5fb.root"    , "read");
TFile* f_elSF_eff  = new TFile(DATA_SF+"/leptonSF/electronSF_trkEff_EWKino_fullsim_M17_36p5fb.root", "read");
TH2F* h_elSF_mvaVT = (TH2F*) f_elSF_id ->Get("GsfElectronToLeptonMvaVTIDEmuTightIP2DSIP3D8mini04");
TH2F* h_elSF_mvaM  = (TH2F*) f_elSF_id ->Get("GsfElectronToLeptonMvaMIDEmuTightIP2DSIP3D8mini04");
TH2F* h_elSF_id    = (TH2F*) f_elSF_id ->Get("GsfElectronToMVAVLooseTightIP2D");
TH2F* h_elSF_trk   = (TH2F*) f_elSF_eff->Get("EGamma_SF2D");

// muons
TFile* f_muSF_mvaVT = new TFile(DATA_SF+"/leptonSF/muonSF_mvaVT_EWKino_fullsim_M17_36p5fb.root", "read");
TFile* f_muSF_mvaM  = new TFile(DATA_SF+"/leptonSF/muonSF_mvaM_EWKino_fullsim_M17_36p5fb.root" , "read");
TFile* f_muSF_id    = new TFile(DATA_SF+"/leptonSF/muonSF_id_EWKino_fullsim_M17_36p5fb.root"   , "read");
TFile* f_muSF_eff   = new TFile(DATA_SF+"/leptonSF/muonSF_trk_EWKino_fullsim_M17_36p5fb.root"  , "read"); 
TH2F* h_muSF_mvaVT = (TH2F*) f_muSF_mvaVT->Get("SF" );
TH2F* h_muSF_mvaM  = (TH2F*) f_muSF_mvaM ->Get("SF" );
TH2F* h_muSF_id    = (TH2F*) f_muSF_id   ->Get("SF" );
TGraphAsymmErrors* h_muSF_trk = (TGraphAsymmErrors*) f_muSF_eff->Get("ratio_eff_eta3_dr030e030_corr");

float getElectronSF(float pt, float eta, int wp = 0){
    TH2F* hist = (wp == 1)?h_elSF_mvaVT:h_elSF_mvaM;
    return getSF(hist, pt, abs(eta))*getSF(h_elSF_id, pt, abs(eta))*getSF(h_elSF_trk, eta, pt);
}

float getElectronUnc(float pt, float eta, int wp = 0, int var = 0){
    TH2F* hist = (wp == 1)?h_elSF_mvaVT:h_elSF_mvaM;
    float error1 = getUnc(hist      , pt , abs(eta));
    float error2 = getUnc(h_elSF_id , pt , abs(eta));
    float error3 = getUnc(h_elSF_trk, eta, pt);
    return TMath::Sqrt(error1*error1 + error2*error2 + error3*error3);
}

float getMuonSF(float pt, float eta, int wp = 0){
    TH2F* hist = (wp == 1)?h_muSF_mvaVT:h_muSF_mvaM;
    return h_muSF_trk->Eval(eta)*getSF(hist, pt, abs(eta))*getSF(h_muSF_id, pt, abs(eta)); 
}

float getMuonUnc(float pt, int var = 0) {
    if (pt<20)  //FIXME: check uncertainty on tracking efficiency once it is available
         return TMath::Sqrt(0.03*0.03+0.01*0.01+0.01*0.01);
    return TMath::Sqrt(0.02*0.02+0.01*0.01);  
}

float getLepSF(float pt, float eta, int pdgId, int isTight, int wp = 0, int var = 0){
    if(!isTight) return 1.0;
    float sf  = 1.0; 
    float err = 0.0;
    if(abs(pdgId) == 11) { sf = getElectronSF(pt, eta, wp); err = getElectronUnc(pt, eta, wp, var); }
    if(abs(pdgId) == 13) { sf = getMuonSF    (pt, eta, wp); err = sf*getMuonUnc (pt, var);          } // only relative error
    if(abs(pdgId) == 15) { sf = 0.95                      ; err = 0.05;                             }
    return (var==0)?sf:(sf+var*err)/sf;

}

float leptonSF(float lepSF1, float lepSF2, float lepSF3 = 1, float lepSF4 = 1){
    return lepSF1*lepSF2*lepSF3*lepSF4;
}



// LEPTON SCALE FACTORS FASTSIM SPRING16
// -------------------------------------------------------------

// electrons
TFile* f_elSF_FS_mvaVT = new TFile(DATA_SF+"/leptonSF/electronSF_mvaVT_EWKino_fastsim_M17_36p5fb.root", "read");
TFile* f_elSF_FS_mvaM  = new TFile(DATA_SF+"/leptonSF/electronSF_mvaM_EWKino_fastsim_M17_36p5fb.root" , "read");
TFile* f_elSF_FS_id    = new TFile(DATA_SF+"/leptonSF/electronSF_id_EWKino_fastsim_M17_36p5fb.root"   , "read");
TH2F* h_elSF_FS_mvaVT  = (TH2F*) f_elSF_FS_mvaVT->Get("histo2D");
TH2F* h_elSF_FS_mvaM   = (TH2F*) f_elSF_FS_mvaM ->Get("histo2D");
TH2F* h_elSF_FS_id     = (TH2F*) f_elSF_FS_id   ->Get("histo2D");

// muons
TFile* f_muSF_FS_mvaVT = new TFile(DATA_SF+"/leptonSF/muonSF_mvaVT_EWKino_fastsim_M17_36p5fb.root", "read");
TFile* f_muSF_FS_mvaM  = new TFile(DATA_SF+"/leptonSF/muonSF_mvaM_EWKino_fastsim_M17_36p5fb.root" , "read");
TFile* f_muSF_FS_id    = new TFile(DATA_SF+"/leptonSF/muonSF_id_EWKino_fastsim_M17_36p5fb.root"   , "read");
TH2F* h_muSF_FS_mvaVT  = (TH2F*) f_muSF_FS_mvaVT->Get("histo2D");
TH2F* h_muSF_FS_mvaM   = (TH2F*) f_muSF_FS_mvaM ->Get("histo2D");
TH2F* h_muSF_FS_id     = (TH2F*) f_muSF_FS_id   ->Get("histo2D");

// taus
TFile* f_tauSF_FS_id = new TFile(DATA_SF+"/leptonSF/tauSF_id_EWKino_fastsim_M17_36p5fb.root", "read");
TH2F* h_tauSF_FS_id  = (TH2F*) f_tauSF_FS_id->Get("histo2D" );

float getElectronSFFS(float pt, float eta, int wp = 0){
    TH2F* hist = (wp == 1)?h_elSF_FS_mvaVT:h_elSF_FS_mvaM;
    return getSF(hist, pt, abs(eta))*getSF(h_elSF_FS_id, pt, abs(eta));
}

float getElectronUncFS(int var = 0){
	return 0.02;
}

float getMuonSFFS(float pt, float eta, int wp = 0){
    TH2F* hist = (wp == 1)?h_muSF_FS_mvaVT:h_muSF_FS_mvaM;
    return getSF(hist, pt, abs(eta))*getSF(h_muSF_FS_id, pt, abs(eta)); 
}

float getMuonUncFS(float pt, int var = 0) {
	return 0.02;
}

float getTauSFFS(float pt, float eta){
    return getSF(h_tauSF_FS_id, pt, abs(eta));
}

float getTauUncFS(float pt, float eta, int var = 0) {
	return getUnc(h_tauSF_FS_id, pt, abs(eta));
}

float getLepSFFS(float pt, float eta, int pdgId, int isTight, int wp = 0, int var = 0){
    if(!isTight) return 1.0;
    float sf  = 1.0; 
    float err = 0.0;
    if(abs(pdgId) == 11) { sf = getElectronSFFS(pt, eta, wp); err = sf*getElectronUncFS(var); } // relative uncertainty
    if(abs(pdgId) == 13) { sf = getMuonSFFS    (pt, eta, wp); err = sf*getMuonUncFS    (var); } // relative uncertainty
    if(abs(pdgId) == 15) { sf = getTauSFFS     (pt, eta    ); err = getTauUncFS(pt, eta, var);}
    return (var==0)?sf:(sf+var*err)/sf;
}

float leptonSFFS(float lepSF1, float lepSF2, float lepSF3 = 1.0, float lepSF4 = 1.0){
    return lepSF1*lepSF2*lepSF3*lepSF4;
}



// LEPTON SCALE FACTORS FASTSIM SUMMER16
// -------------------------------------------------------------

// electrons
TFile* f_elSF_FS16_mvaVT = new TFile(DATA_SF+"/leptonSF/electronSF_mvaVT_EWKino_fastsim_Summer16_M17_36p5fb.root", "read");
TFile* f_elSF_FS16_mvaM  = new TFile(DATA_SF+"/leptonSF/electronSF_mvaM_EWKino_fastsim_Summer16_M17_36p5fb.root" , "read");
TH2F* h_elSF_FS16_mvaVT  = (TH2F*) f_elSF_FS16_mvaVT->Get("histo2D");
TH2F* h_elSF_FS16_mvaM   = (TH2F*) f_elSF_FS16_mvaM ->Get("histo2D");

// muons
TFile* f_muSF_FS16_mvaVT = new TFile(DATA_SF+"/leptonSF/muonSF_mvaVT_EWKino_fastsim_Summer16_M17_36p5fb.root", "read");
TFile* f_muSF_FS16_mvaM  = new TFile(DATA_SF+"/leptonSF/muonSF_mvaM_EWKino_fastsim_Summer16_M17_36p5fb.root" , "read");
TH2F* h_muSF_FS16_mvaVT  = (TH2F*) f_muSF_FS16_mvaVT->Get("histo2D");
TH2F* h_muSF_FS16_mvaM   = (TH2F*) f_muSF_FS16_mvaM ->Get("histo2D");

float getElectronSFFS16(float pt, float eta, int wp = 0){
    TH2F* hist = (wp == 1)?h_elSF_FS16_mvaVT:h_elSF_FS16_mvaM;
    return getSF(hist, pt, abs(eta));
}

float getElectronUncFS16(int var = 0){
	return 0.02;
}

float getMuonSFFS16(float pt, float eta, int wp = 0){
    TH2F* hist = (wp == 1)?h_muSF_FS16_mvaVT:h_muSF_FS16_mvaM;
    return getSF(hist, pt, abs(eta));
}

float getMuonUncFS16(float pt, int var = 0) {
	return 0.02;
}

float getTauSFFS16(float pt, float eta){
	return 1.0;
}

float getTauUncFS16(float pt, float eta, int var = 0) {
	return 0.02;
}

float getLepSFFS16(float pt, float eta, int pdgId, int isTight, int wp = 0, int var = 0){
    if(!isTight) return 1.0;
    float sf  = 1.0; 
    float err = 0.0;
    if(abs(pdgId) == 11) { sf = getElectronSFFS16(pt, eta, wp); err = sf*getElectronUncFS16(var); } // relative uncertainty
    if(abs(pdgId) == 13) { sf = getMuonSFFS16    (pt, eta, wp); err = sf*getMuonUncFS16    (var); } // relative uncertainty
    if(abs(pdgId) == 15) { sf = getTauSFFS16     (pt, eta    ); err = sf*getTauUncFS16(pt, eta, var);} // relative uncertainty
    return (var==0)?sf:(sf+var*err)/sf;
}

float leptonSFFS16(float lepSF1, float lepSF2, float lepSF3 = 1.0, float lepSF4 = 1.0){
    return lepSF1*lepSF2*lepSF3*lepSF4;
}



void functionsSF() {}
