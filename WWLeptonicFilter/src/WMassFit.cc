#include "WWAnalysis/WWLeptonicFilter/interface/WMassFit.h"

WMassFit::WMassFit(const edm::ParameterSet& cfg) :
  genmet_tk(consumes<edm::View<reco::GenMET> >(cfg.getParameter<edm::InputTag>("genmet")))
, met_tk(consumes<edm::View<reco::PFMET> >(cfg.getParameter<edm::InputTag>("met")))
, ww_tk(consumes<Recons::WWLeptonicCandidate>(edm::InputTag("wwLeptonicBuilder","ww")))
, ctpps_tk(consumes<Recons::CTPPS>(cfg.getParameter<edm::InputTag>("mixLabel")))
{
  produces<Recons::FitResults2>("fitResults");
  fitResults = new Recons::FitResults2();

}

double ChiSQ(double (*ConstrainFunction)(const TLorentzVector&,const TLorentzVector&), const TLorentzVector& v1,const TLorentzVector& v2, int idx) {
   return pow((Data[idx]-ConstrainFunction(v1,v2))/(CF*ErrData[idx]),2);
}

void FitFunction(int &npar, double *gin, double &f, double *par, int iflag){
/*
     switch(iflag) {
           case 1: //read input data - done at the steering routine
                   break;
           case 2: //calculate the first derivative in gin - not used (optional)
                   break;
           case 3: // end of fit
                   break;
     }
*/
     TLorentzVector nu1,nu2;
     nu1.SetPxPyPzE(par[0],par[1],par[2],sqrt(par[0]*par[0]+par[1]*par[1]+par[2]*par[2]));
     nu2.SetPxPyPzE(par[3],par[4],par[5],sqrt(par[3]*par[3]+par[4]*par[4]+par[5]*par[5]));
     TLorentzVector w1=muon1+nu1;
     TLorentzVector w2=muon2+nu2;
     double chisq=0.;
//
     chisq+=ChiSQ(Wmass,muon1,nu1,MASS_W1);
     chisq+=ChiSQ(Wmass,muon2,nu2,MASS_W2);
     chisq+=ChiSQ(METX,nu1,nu2,MET_X);
     chisq+=ChiSQ(METY,nu1,nu2,MET_Y);
     chisq+=ChiSQ(METPt,nu1,nu2,MET_PT);
     chisq+=ChiSQ(METPhi,nu1,nu2,MET_PHI);
     chisq+=ChiSQ(METE,nu1,nu2,MET_E);
     chisq+=ChiSQ(CentralMass,w1,w2,MX);
     chisq+=ChiSQ(Rapidity,w1,w2,Y);
     chisq+=ChiSQ(Acoplanarity,w1,w2,ACOPLANARITY);
     chisq+=ChiSQ(WWSumPxPy,w1,w2,WW_PXPY);
     f=chisq;
}
bool WmassFit(TLorentzVector& nu1,TLorentzVector& nu2,double& chi2) {
     TMinuit *gMinuit = new TMinuit(NPAR);  //initialize TMinuit with a maximum of 6 params
     gMinuit->SetFCN(FitFunction);
     gMinuit->Command("SET PRINT -1");

     Double_t arglist[10];
     Int_t ierflg = 0;

     arglist[0] = 1;
     gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

     // Set starting values and step sizes for parameters - starts with the muon momentum
     static Double_t vstart[7] = {nu1.Px(), nu1.Py() , nu1.Pz(), nu2.Px(),nu2.Py(),nu2.Pz()};
     static Double_t step[7] =   {0.1     , 0.1      , 1.       , 0.1     ,0.1     ,1.};
     gMinuit->mnparm(0, "n1_px", vstart[0], step[0], 0,0,ierflg);
     gMinuit->mnparm(1, "n1_py", vstart[1], step[1], 0,0,ierflg);
     gMinuit->mnparm(2, "n1_pz", vstart[2], step[2], 0,0,ierflg);
     gMinuit->mnparm(3, "n2_px", vstart[3], step[3], 0,0,ierflg);
     gMinuit->mnparm(4, "n2_py", vstart[4], step[4], 0,0,ierflg);
     gMinuit->mnparm(5, "n2_pz", vstart[5], step[5], 0,0,ierflg);

     // Now ready for minimization step
     arglist[0] = 2000;
     arglist[1] = 1.;
     gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
     gMinuit->mnexcm("HESSE", arglist ,2,ierflg);


     // Print results
     Double_t amin,edm,errdef;
     Int_t nvpar,nparx,icstat;
     gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
     std::cout << "Fit status : " << icstat << std::endl;
     //if (icstat==3) {
        double px,py,pz;
        TString name;
        double pxer,pyer,pzer;
        double low,high;
        int ipar;
        gMinuit->mnpout(0,name,px,pxer,low,high,ipar);
        gMinuit->mnpout(1,name,py,pyer,low,high,ipar);
        gMinuit->mnpout(2,name,pz,pzer,low,high,ipar);
        nu1.SetPxPyPzE(px,py,pz,sqrt(px*px+py*py+pz*pz));
        gMinuit->mnpout(3,name,px,pxer,low,high,ipar);
        gMinuit->mnpout(4,name,py,pyer,low,high,ipar);
        gMinuit->mnpout(5,name,pz,pzer,low,high,ipar);
        nu2.SetPxPyPzE(px,py,pz,sqrt(px*px+py*py+pz*pz));
        chi2=amin;
     //}

  return icstat==3;
}

bool WPairMatching(edm::Handle<edm::View<reco::GenMET> > genmet,
                   edm::Handle<edm::View<reco::PFMET> > met,
                   edm::Handle<Recons::WWLeptonicCandidate> ww,
                   edm::Handle<Recons::CTPPS> ctpps,
                   Recons::FitResults2* fitResults
                  ) {
     TLorentzVector nu1gen,nu2gen;
     TLorentzVector mu1gen,mu2gen;
     TLorentzVector w1gen,w2gen;
     TLorentzVector metgen;
     TLorentzVector nu1,nu2;
     TLorentzVector w1,w2;
     double chi2=0.;
     const double ecms = 13000.0;
     double muonmass = 0.10565837; // mass in GeV
     
//
     Data[MASS_W1] = 80.385; // W 1 mass
     Data[MASS_W2] = 80.385; // W 2 mass
     Data[ACOPLANARITY] = 0.; // Acoplanarity
     Data[WW_PXPY]=0.;//-ecms*sin(180e-6); // require transverse momentum conservation
     ErrData[MASS_W1] = 0.015;//2.085; // take the W decay width as weight
     ErrData[MASS_W2] = 0.015;//2.085;
     ErrData[ACOPLANARITY] = 0.002;// assume the dominant uncertainty the one in MET_phi ;//fabs(5.*TMath::Pi()/180.);// allow a 10 degree around zero

      double pxnu1 = ww->genNeutrino1().Px();
      double pynu1 = ww->genNeutrino1().Py();
      double pznu1 = ww->genNeutrino1().Pz();
      double enu1 = ww->genNeutrino1().E();
      double pxnu2 = ww->genNeutrino2().Px();
      double pynu2 = ww->genNeutrino2().Py();
      double pznu2 = ww->genNeutrino2().Pz();
      double enu2 = ww->genNeutrino2().E();
      double pxmu1g = ww->genMuon1().Px();
      double pymu1g = ww->genMuon1().Py();
      double pzmu1g = ww->genMuon1().Pz();
      double emu1g = ww->genMuon1().E();
      double pxmu2g = ww->genMuon2().Px();
      double pymu2g = ww->genMuon2().Py();
      double pzmu2g = ww->genMuon2().Pz();
      double emu2g = ww->genMuon2().E();
      double pxmetg = genmet.product()->front().px();
      double pymetg = genmet.product()->front().py();
      double pzmetg = genmet.product()->front().pz();
      double emetg = genmet.product()->front().energy();
//      pxmu1 = ww->recoMuon1().Px();
//      pymu1 = ww->recoMuon1().Py();
//      pzmu1 = ww->recoMuon1().Pz();
      double ptmu1 = ww->recoMuon1().Pt();
      double etamu1 = ww->recoMuon1().Eta();
      double phimu1 = ww->recoMuon1().Phi();
//      emu1 = ww->recoMuon1().E();
//      pxmu2 = ww->recoMuon2().Px();
//      pymu2 = ww->recoMuon2().Py();
//      pzmu2 = ww->recoMuon2().Pz();
      double ptmu2 = ww->recoMuon2().Pt();
      double etamu2 = ww->recoMuon2().Eta();
      double phimu2 = ww->recoMuon2().Phi();
//      emu2 = ww->recoMuon2().E();
      double metx = met.product()->front().px();
      double mety = met.product()->front().py();
      double metpt = met.product()->front().pt();
      double metphi = met.product()->front().phi();
      double mete = met.product()->front().energy();
      double xi45near = (ctpps->xi45near()!=Dummy::double_ && ctpps->xi45near()!=0) ? ctpps->xi45near() : 0;
      double xi45far = (ctpps->xi45far()!=Dummy::double_ && ctpps->xi45far()!=0) ? ctpps->xi45far() : 0;
      double xi56near = (ctpps->xi56near()!=Dummy::double_ && ctpps->xi56near()!=0) ? ctpps->xi56near() : 0;
      double xi56far = (ctpps->xi56far()!=Dummy::double_ && ctpps->xi56far()!=0) ? ctpps->xi56far() : 0;
      double xi45D = 0;
      double xi56D = 0;
//
// chose the best xi reconstruction
//
        double xi1 = 0;
        double xi2 = 0;

        xi45D=0.;xi56D=0.; // use the single RP by now
        if (xi45D) {xi1=xi45D;}
        else if (xi45near) {xi1=xi45near;}
        else {xi1=xi45far;};
        if (xi56D) {xi1=xi56D;}
        else if (xi56near) {xi1=xi56near;}
        else {xi1=xi56far;};

        xi1 = xi45near; xi2=xi56near;

        if (xi1*xi2==0||xi1*xi2>1.) 
          {std::cout << "NOT CEP event." << std::endl; return false;}
//          throw cms::Exception("WMassFit") << "\n +" << __FILE__ << ":" << __LINE__ << ":"
//          << "\n filter(): NOT CEP event. Aborting.";
        //{std::cout << "NOT CEP event." << std::endl;continue;}

//        std::cout << "Processing event "<<nevt<<std::endl;
      
        nu1gen.SetPxPyPzE(pxnu1,pynu1,pznu1,enu1);
        nu2gen.SetPxPyPzE(pxnu2,pynu2,pznu2,enu2);
        mu1gen.SetPxPyPzE(pxmu1g,pymu1g,pzmu1g,emu1g);
        mu2gen.SetPxPyPzE(pxmu2g,pymu2g,pzmu2g,emu2g);
        metgen.SetPxPyPzE(pxmetg,pymetg,pzmetg,emetg);
        w1gen=mu1gen+nu1gen;
        w2gen=mu2gen+nu2gen;
//
        muon1.SetPtEtaPhiM(ptmu1,etamu1,phimu1,muonmass);
        muon2.SetPtEtaPhiM(ptmu2,etamu2,phimu2,muonmass);
//        std::cout << "Muon 1 " << ptmu1 << " " << etamu1 << " " << phimu1 << std::endl;
//        std::cout << "Muon 2 " << ptmu2 << " " << etamu2 << " " << phimu2 << std::endl;
//
	Data[MET_X] = metx; // MET_X - read from file
	Data[MET_Y] = mety; // MET_Y
	Data[MET_E] = mete; // MET_SumET
	Data[MX] = sqrt(xi1*xi2)*ecms; // M_X - read from file
	Data[Y] = 0.5*log(xi1/xi2); // Rapidy
        Data[MET_PT] = metpt;
        Data[MET_PHI] = metphi;
        double sig_xi=0.055;
        double sig_xi1 = sig_xi*xi1;
        double sig_xi2 = sig_xi*xi2;
        ErrData[MET_X] = fabs(0.13*metx);
        ErrData[MET_Y] = fabs(0.15*mety);
        ErrData[WW_PXPY] =ecms*sin(180e-6);
        ErrData[MET_E] = fabs(0.08*mete);
        ErrData[MX] = 0.5*sqrt((pow(sig_xi1*xi2,2)+pow(sig_xi2*xi1,2))/(xi1*xi2))*ecms;
        ErrData[Y] =  0.5*sqrt(pow(sig_xi1*xi2,2)+pow(sig_xi2*xi1,2))/(xi1*xi2);
        ErrData[MET_PT] = fabs(0.08*Data[MET_PT]);
        ErrData[MET_PHI] = 0.007;//fabs(2*TMath::Pi()/180.);//0.0065;//fabs(10*TMath::Pi()/180.); // 10 degree 
        double pz=mete-(muon1.E()+muon2.E()+sqrt(metx*metx+mety*mety));
        double ediff=mete-(muon1.E()+muon2.E());
        nu1.SetPxPyPzE(metx/2.,mety/2.,pz,ediff);
        nu2.SetPxPyPzE(metx/2.,mety/2.,-pz,ediff);
/*
        std::cout << "\n [before] nu1.Px()=" << nu1.Px();
        std::cout << "\n [before] nu1.Py()=" << nu1.Py();
        std::cout << "\n [before] nu1.Pz()=" << nu1.Pz();
        std::cout << "\n [before] nu2.Px()=" << nu2.Px();
        std::cout << "\n [before] nu2.Py()=" << nu2.Py();
        std::cout << "\n [before] nu2.Pz()=" << nu2.Pz();
*/

        bool converged = WmassFit(nu1,nu2,chi2);

        w1=(muon1+nu1);
        w2=(muon2+nu2);
//
        double wmass_1 = Wmass(muon1,nu1);
        double wmass_2 = Wmass(muon2,nu2);
        double ChiW1 = ChiSQ(Wmass,muon1,nu1,MASS_W1)/NPAR;
        double ChiW2 = ChiSQ(Wmass,muon2,nu2,MASS_W2)/NPAR;
        double ChiMETX=ChiSQ(METX,nu1,nu2,MET_X)/NPAR;
        double ChiMETY=ChiSQ(METY,nu1,nu2,MET_Y)/NPAR;
        double ChiMETPhi=ChiSQ(METPhi,nu1,nu2,MET_PHI)/NPAR;
        double ChiMETPt=ChiSQ(METPt,nu1,nu2,MET_PT)/NPAR;
        double ChiMETE=ChiSQ(METE,nu1,nu1,MET_E)/NPAR;
        double ChiCM=ChiSQ(CentralMass,w1,w2,MX)/NPAR;
        double ChiY=ChiSQ(Rapidity,w1,w2,Y)/NPAR;
        double ChiAcop=ChiSQ(Acoplanarity,w1,w2,ACOPLANARITY)/NPAR;
        double ChiPxPy=ChiSQ(WWSumPxPy,w1,w2,WW_PXPY)/NPAR;
//

        double chi2ndf = chi2/NPAR;

        fitResults->clear();

        if(fitResults->nu1ResXNorm()==-1) std::cout << "\n fitResults->nu1ResXNorm()=" << fitResults->nu1ResXNorm();
        if(fitResults->nu1ResYNorm()==-1) std::cout << "\n fitResults->nu1ResYNorm()=" << fitResults->nu1ResYNorm();
        if(fitResults->nu1ResZNorm()==-1) std::cout << "\n fitResults->nu1ResZNorm()=" << fitResults->nu1ResZNorm();

        if(fitResults->nu2ResXNorm()==-1) std::cout << "\n fitResults->nu2ResXNorm()=" << fitResults->nu2ResXNorm();
        if(fitResults->nu2ResYNorm()==-1) std::cout << "\n fitResults->nu2ResYNorm()=" << fitResults->nu2ResYNorm();
        if(fitResults->nu2ResZNorm()==-1) std::cout << "\n fitResults->nu2ResZNorm()=" << fitResults->nu2ResZNorm();


        fitResults->setChi2(chi2);
        fitResults->setChi2ndf(chi2ndf);
        fitResults->setChiW1(ChiW1);
        fitResults->setChiW2(ChiW2);
        fitResults->setChiMETX(ChiMETX);
        fitResults->setChiMETY(ChiMETY);
        fitResults->setChiMETPhi(ChiMETPhi);
        fitResults->setChiMETPt(ChiMETPt);
        fitResults->setChiMET(ChiMETE);
        fitResults->setChiCM(ChiCM);
        fitResults->setChiY(ChiY);
        fitResults->setChiAcop(ChiAcop);
        fitResults->setChiPxPy(ChiPxPy);
        fitResults->setNu1Gen(nu1gen);
        fitResults->setNu2Gen(nu2gen);
        fitResults->setNu1Fitted(nu1);
        fitResults->setNu2Fitted(nu2);
        fitResults->setRes(); // resolutions

       std::cout << "\n nu1gen(px,py,pz,eta,phi)=(" << nu1gen.Px() << "," << nu1gen.Py() << "," << nu1gen.Pz() << "," << nu1gen.Eta() << "," << nu1gen.Phi() << ")";
       std::cout << "\n nu2gen(px,py,pz,eta,phi)=(" << nu2gen.Px() << "," << nu2gen.Py() << "," << nu2gen.Pz() << "," << nu2gen.Eta() << "," << nu2gen.Phi() << ")";
       std::cout << "\n nu1(px,py,pz,eta,phi)=(" << nu1.Px() << "," << nu1.Py() << "," << nu1.Pz() << "," << nu1.Eta() << "," << nu1.Phi() << ")";
       std::cout << "\n nu2(px,py,pz,eta,phi)=(" << nu2.Px() << "," << nu2.Py() << "," << nu2.Pz() << "," << nu2.Eta() << "," << nu2.Phi() << ")\n\n";
/*
        std::cout << "\n [after] nu1.Px()=" << nu1.Px();
        std::cout << "\n [after] nu1.Py()=" << nu1.Py();
        std::cout << "\n [after] nu1.Pz()=" << nu1.Pz();
        std::cout << "\n [after] nu2.Px()=" << nu2.Px();
        std::cout << "\n [after] nu2.Py()=" << nu2.Py();
        std::cout << "\n [after] nu2.Pz()=" << nu2.Pz();
  std::cout << "\n [after] muon1.Px()=" << muon1.Px();
  std::cout << "\n [after] muon1.Py()=" << muon1.Py();
  std::cout << "\n [after] muon1.Pz()=" << muon1.Pz();
  std::cout << "\n [after] muon1.E()=" << muon1.E();
  std::cout << "\n [after] muon2.Px()=" << muon2.Px();
  std::cout << "\n [after] muon2.Py()=" << muon2.Py();
  std::cout << "\n [after] muon2.Pz()=" << muon2.Pz();
  std::cout << "\n [after] muon2.E()=" << muon2.E();
        std::cout << "\n fitResults->nu1ResX()=" << fitResults->nu1ResX();
        std::cout << "\n fitResults->nu1ResY()=" << fitResults->nu1ResY();
        std::cout << "\n fitResults->nu1ResZ()=" << fitResults->nu1ResZ();
        std::cout << "\n fitResults->nu2ResX()=" << fitResults->nu2ResX();
        std::cout << "\n fitResults->nu2ResY()=" << fitResults->nu2ResY();
        std::cout << "\n fitResults->nu2ResZ()=" << fitResults->nu2ResZ();
*/


/*
        std::cout << "Rapidity : " << 0.5*log((muon1.E()+nu1.E()+muon1.Pz()+nu1.Pz())/(muon1.E()+nu1.E()-muon1.Pz()-nu1.Pz())) << " " << w1.Rapidity() << std::endl;
        std::cout << "Rapidity : " << 0.5*log((muon2.E()+nu2.E()+muon2.Pz()+nu2.Pz())/(muon2.E()+nu2.E()-muon2.Pz()-nu2.Pz())) << " " << w2.Rapidity() << std::endl;
        std::cout <<"Fit chi2 : " << chi2 << std::endl;
        std::cout <<"Chi2/NDF : " << chi2/NPAR << std::endl;
        std::cout <<"Chi2  MX : " << ChiSQ(CentralMass,(muon1+nu1),(muon2+nu2),MX) << std::endl;
        std::cout <<"Neutrino 1 : " << nu1.Px() << " " << nu1.Py() << " " << nu1.Pz() << std::endl;
        std::cout <<"Neutrino 2 : " << nu2.Px() << " " << nu2.Py() << " " << nu2.Pz() << std::endl;
        std::cout <<"Fitted parameters "<<std::endl;
        std::cout << "W_1 Mass     : "<< wmass_1 << std::endl;
        std::cout << "W_2 Mass     : "<< wmass_2 << std::endl;
        std::cout << "Acoplanarity : " << Acoplanarity((muon1+nu1),(muon2+nu2)) << std::endl;
        std::cout << "Central Mass : "<< Data[MX] << " --- " << CentralMass((muon1+nu1),(muon2+nu2))<<std::endl;
        std::cout << "Central Rap. : "<< Data[Y] << " --- " << Rapidity((muon1+nu1),(muon2+nu2))<<std::endl;
        std::cout << "MET X        : "<< Data[MET_X] << " --- " << (nu1+nu2).Px() << std::endl;
        std::cout << "MET Y        : "<< Data[MET_Y] << " --- " << (nu1+nu2).Py() << std::endl;
        std::cout << "MET Phi      : "<< Data[MET_PHI] << " --- " << (nu1+nu2).Phi() << std::endl;
        std::cout << "MET Pt       : "<< Data[MET_PT] << " --- " << (nu1+nu2).Pt()  << std::endl;
        std::cout << "WW Px+Py     : "<< Data[WW_PXPY] << " --- " << WWSumPxPy((muon1+nu1),(muon2+nu2))  << std::endl;
        std::cout << "==============================================================================\n" << std::endl;
*/
  return converged;
}


bool WMassFit::filter(edm::Event& event, const edm::EventSetup& setup)
{
  event.getByToken(genmet_tk,genmet);
  event.getByToken(ww_tk,ww);
  event.getByToken(met_tk,met);
  event.getByToken(ctpps_tk,ctpps);

  // only if collections are available
  if( !genmet.isValid() 
   || !ww.isValid() 
   || !met.isValid() 
   || !ctpps.isValid() 
    ) {
    throw cms::Exception("WMassFit") << "\n +" << __FILE__ << ":" << __LINE__ << ":"
    << "\n filter(): One of the following collections could not be retrieved" 
    << "\n filter(): from root file:"
    << "\n filter(): genmet.isValid()=" << genmet.isValid()
    << "\n filter(): ww.isValid()=" << ww.isValid()
    << "\n filter(): met.isValid()=" << met.isValid()
    << "\n filter(): ctpps.isValid()=" << ctpps.isValid()
    << "\n filter(): Aborting.\n\n";
  }


  muon1.SetPxPyPzE(0,0,0,0);
  muon2.SetPxPyPzE(0,0,0,0);
/*
  std::cout << "\n [before] muon1.Px()=" << muon1.Px();
  std::cout << "\n [before] muon1.Py()=" << muon1.Py();
  std::cout << "\n [before] muon1.Pz()=" << muon1.Pz();
  std::cout << "\n [before] muon1.E()=" << muon1.E();
  std::cout << "\n [before] muon2.Px()=" << muon2.Px();
  std::cout << "\n [before] muon2.Py()=" << muon2.Py();
  std::cout << "\n [before] muon2.Pz()=" << muon2.Pz();
  std::cout << "\n [before] muon2.E()=" << muon2.E();
*/
  bool result = WPairMatching(genmet,met,ww,ctpps,fitResults);

  // put in the event
  std::unique_ptr<Recons::FitResults2> pfitResults(new Recons::FitResults2); 
  *pfitResults = *fitResults;

  event.put(std::move(pfitResults),"fitResults");

  return result;
}



#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(WMassFit);
