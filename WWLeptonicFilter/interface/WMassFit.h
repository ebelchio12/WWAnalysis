#ifndef WWFilter_WMassFit_h
#define WWFilter_WMassFit_h

// default include files
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// user include files
#include "FWCore/Utilities/interface/Exception.h"
#include "WWAnalysis/WWUtils/interface/FitResults2.h"
#include "WWAnalysis/WWUtils/interface/CTPPS.h"
#include "WWAnalysis/WWUtils/interface/WWLeptonicCandidate.h"
#include "DataFormats/METReco/interface/MET.h"
#include "PhysicsTools/Utilities/interface/RootMinuit.h"
#include "TROOT.h"
#include "TMinuit.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TString.h"
#include "TTree.h"
#include "TFile.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cctype>

enum {MASS_W1,MASS_W2,MET_X,MET_Y,MET_E,MX,Y,ACOPLANARITY,MET_PT,MET_PHI,WW_PXPY,NConstraints};
double Data[NConstraints];
double ErrData[NConstraints];
int NPAR = 6;
int CF=1.0; // Confidence Level on data weights

TLorentzVector muon1;
TLorentzVector muon2;

double ChiSQ(double (*ConstrainFunction)(const TLorentzVector&,const TLorentzVector&), const TLorentzVector& v1,const TLorentzVector& v2, int idx);

double Wmass(const TLorentzVector& mu, const TLorentzVector& nu) {
   return (mu+nu).M();
};
double METX(const TLorentzVector& nu1,const TLorentzVector& nu2) {
   return (nu1+nu2).Px();
};
double METY(const TLorentzVector& nu1,const TLorentzVector& nu2) {
   return (nu1+nu2).Py();
};
double METE(const TLorentzVector& nu1,const TLorentzVector& nu2) {
   return (nu1+nu2).E();
};
double Rapidity(const TLorentzVector& w1, const TLorentzVector& w2){
   return (w1+w2).Rapidity();
}
double CentralMass(const TLorentzVector& w1, const TLorentzVector& w2) {
   return (w1+w2).M();
}
double Acoplanarity(const TLorentzVector& w1, const TLorentzVector& w2) {
   return 1.-fabs((w1.DeltaPhi(w2))/TMath::Pi());
}
double WWSumPxPy(const TLorentzVector& w1, const TLorentzVector& w2) {
   return w1.Px()+w1.Py()+w2.Px()+w2.Py();
}
double METPhi(const TLorentzVector& nu1, const TLorentzVector& nu2) {
   return (nu1+nu2).Phi();
}
double METPt(const TLorentzVector& nu1, const TLorentzVector& nu2) {
   return (nu1+nu2).Pt();
}

void FitFunction(int &npar, double *gin, double &f, double *par, int iflag);

bool WmassFit(TLorentzVector& nu1,TLorentzVector& nu2,double& chi2);

bool WPairMatching(edm::Handle<edm::View<reco::GenMET> > genmet,
                   edm::Handle<edm::View<reco::PFMET> > met,
                   edm::Handle<Recons::WWLeptonicCandidate> ww,
                   edm::Handle<Recons::CTPPS> ctpps,
                   Recons::FitResults2* fitResults
                  );


// class declaration
class WMassFit : public edm::EDFilter {

  public:

    explicit WMassFit(const edm::ParameterSet&);
    ~WMassFit(){};


  private:

    virtual bool filter(edm::Event&, const edm::EventSetup&);

    edm::EDGetTokenT<edm::View<reco::GenMET> > genmet_tk;
    edm::EDGetTokenT<edm::View<reco::PFMET> > met_tk;
    edm::EDGetTokenT<Recons::WWLeptonicCandidate> ww_tk;
    edm::EDGetTokenT<Recons::CTPPS> ctpps_tk;

    edm::Handle<edm::View<reco::GenMET> > genmet; 
    edm::Handle<edm::View<reco::PFMET> > met;
    edm::Handle<Recons::WWLeptonicCandidate> ww;
    edm::Handle<Recons::CTPPS> ctpps;

    //----------
    Recons::FitResults2* fitResults;

};

#endif 
