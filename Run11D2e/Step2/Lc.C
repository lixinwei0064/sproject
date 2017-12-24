#include "TLorentzVector.h"
#include "TRandom.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include <iostream>
#include "charm_decay.h"
using namespace std;

TLorentzVector myBoost(TLorentzVector parent,TLorentzVector daughter){
       double mass               = parent.M();
    TVector3 eta = (-1./mass)*parent.Vect();            // gamma*beta
    double gamma              = fabs(parent.E())/mass;
    //    cout<<daughter.vect()<<" "<<eta<<endl;
    TVector3 pl  = ((daughter.Vect()*eta)/(eta*eta))*eta;  // longitudinal momentum
    return TLorentzVector(daughter.Vect() + (gamma-1.)*pl - daughter.E()*eta,
			  gamma*daughter.E() - daughter.Vect()*eta);
}

TLorentzVector decay_kinematics(TLorentzVector Lcp,TH1D *formfactor){
   double e2,px,py,pz,ptot2,phi,costheta2;
   TLorentzVector leptonp3;
   ptot2 = formfactor->GetRandom();
   costheta2 = myRandom->Uniform(-1.0,1.0);
   phi=myRandom->Uniform(0.,6.283);
   pz = ptot2*costheta2;  // theta lepton relative to thetav see p482 PRD 54 (1996)
   px = ptot2*sqrt(1.0-costheta2*costheta2)*cos(phi);
   py = ptot2*sqrt(1.0-costheta2*costheta2)*sin(phi);
   e2 = sqrt(ptot2*ptot2+Me*Me);
   TVector3 l3v(px,py,pz);
   TLorentzVector leptonp(l3v,e2);
   leptonp3 = myBoost(Lcp,leptonp);
   return leptonp3;
}

void Lc(const int nEvents = 10000000, const int cen = 0, const int seed = 32553, TLegend *legend = new TLegend(0.1,0.7,0.48,0.9)){

  const double c2eBR = 0.045;
  const double BRLc  = 1;
  const double wt = 20./2./3.14159265*c2eBR/BRLc/2.;

  myRandom->SetSeed(seed);

  TH1F *HistoLc;

  TFile *fin = new TFile("cHadron.root");
  if(cen == 0) HistoLc = (TH1F *)fin->Get("hLc_spectra_10_60_model1");
  if(cen == 1) HistoLc = (TH1F *)fin->Get("hLc_spectra_10_60_model2");
  if(cen == 2) HistoLc = (TH1F *)fin->Get("hLc_spectra_10_60_model3");
  if(cen == 3) HistoLc = (TH1F *)fin->Get("hLc_spectra_10_60_mean");
  if(!HistoLc) exit(0);
  HistoLc->SetName("HistoLc");
  
  float ptwtnorm = 20./HistoLc->Integral(0.,20.);
  TH1D *einvyLc = new TH1D("einvyLc","",200,0.,20.0);
  TH1D *Lcinvy = new TH1D("Lcinvy","",200,0.,20.0);
  double p1,p2,p3,ptRandom,phiRandom,yRandom;
  TLorentzVector leptonp(0.,0.,0.,0.);


  TFile *fformfactor = new TFile("charm_edecayforms_atrest.root");
  TH1D *myatrest = (TH1D*)fformfactor->Get("a15");

  for (int Event = 0;Event<nEvents;Event++){
    if(Event%1000000==0)cout<<"Events "<<Event<<endl;
    ptRandom = myRandom->Uniform(0.,20.);
    double value;
    value = HistoLc->Interpolate(ptRandom);
    float ptwt = value*ptwtnorm;
    phiRandom = myRandom->Uniform(0,6.283);
    p1 = ptRandom*cos(phiRandom);
    p2 = ptRandom*sin(phiRandom);
    yRandom = myRandom->Gaus(0.,1.9);
    p3 = sqrt(ML*ML+ptRandom*ptRandom)*(exp(yRandom)-exp(-yRandom))/2.0;
    double e = sqrt(ML*ML+ptRandom*ptRandom+p3*p3);

    TLorentzVector Lcp(TVector3(p1,p2,p3),e);
    leptonp=decay_kinematics(Lcp,myatrest);
    
    if(fabs(Lcp.Rapidity())<1.) 
      Lcinvy->Fill(ptRandom,ptwt*wt/ptRandom);  
    if(fabs(leptonp.PseudoRapidity())<1.)
      einvyLc->Fill(leptonp.Perp(),ptwt*wt/leptonp.Perp());
  }

  Lcinvy->Scale(85);
  einvyLc->Scale(85);

  einvyLc->SetLineStyle(7);
  einvyLc->SetLineColor(kRed);
  einvyLc->Draw("SAME");
  Lcinvy->SetLineStyle(1);
  Lcinvy->SetLineColor(kRed);
  Lcinvy->Draw("SAME");

  legend->AddEntry(einvyLc,"e inv. yield #Lambda_{c}","l");
  legend->AddEntry(Lcinvy,"#Lambda_{c} inv. yield","l");
  legend->Draw();
}
