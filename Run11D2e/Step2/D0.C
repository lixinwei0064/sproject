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
const char name[2][250] = {"0_80","10_40"};

TLorentzVector myBoost(TLorentzVector parent,TLorentzVector daughter){
       double mass               = parent.M();
    TVector3 eta = (-1./mass)*parent.Vect();            // gamma*beta
    double gamma = fabs(parent.E())/mass;
    TVector3 pl  = ((daughter.Vect()*eta)/(eta*eta))*eta;  // longitudinal momentum
    return TLorentzVector(daughter.Vect() + (gamma-1.)*pl - daughter.E()*eta,
			  gamma*daughter.E() - daughter.Vect()*eta);
}

TLorentzVector decay_kinematics(TLorentzVector D0p,TH1D *formfactor){
   double e2,px,py,pz,ptot2,phi,costheta2;
   TLorentzVector leptonp3;
   ptot2 = formfactor->GetRandom();
   costheta2 = myRandom->Uniform(-1.0,1.0);
   phi=myRandom->Uniform(0.,6.283);
   pz = ptot2*costheta2;  
   px = ptot2*sqrt(1.0-costheta2*costheta2)*cos(phi);
   py = ptot2*sqrt(1.0-costheta2*costheta2)*sin(phi);
   e2 = sqrt(ptot2*ptot2+Me*Me);
   TVector3 l3v(px,py,pz);
   TLorentzVector leptonp(l3v,e2);
   leptonp3 = myBoost(D0p,leptonp);
   return leptonp3;
}
void D0(const int nEvents = 10000000, const int cen = 0, const int seed = 32553, TLegend *legend = new TLegend(0.1,0.7,0.48,0.9)){

  const double c2eBR = 0.065;
  const double BRD0  = 1;
  const double wt = 20./2./3.14159265*c2eBR/BRD0/2.;
  char Name[250];
  myRandom->SetSeed(seed);

  TF1 *Levy;  
  TFile *fin = new TFile("cHadron.root");
  sprintf(Name, "D0_%s_levy",name[cen]);
  Levy = (TF1 *)fin->Get(Name);
  if(!Levy) exit(0);  
  Levy->SetName("Levy");

  float ptwtnorm = 20./Levy->Integral(0.,20.);

  TH1D *einvy = new TH1D("einvyD0","",200,0.,20.0);
  TH1D *Dinvy = new TH1D("Dinvy","",200,0.,20.0);
  double p1,p2,p3,ptRandom,phiRandom,yRandom;
  TLorentzVector leptonp(0.,0.,0.,0.);

  TFile *fformfactor = new TFile("charm_edecayforms_atrest.root");
  TH1D *myatrest = (TH1D*)fformfactor->Get("a15");

  for (int Event = 0;Event<nEvents;Event++){
    if(Event%1000000==0)cout<<"Events "<<Event<<endl;
    ptRandom = myRandom->Uniform(0.,20.);
    float ptwt = Levy->Eval(ptRandom)*ptwtnorm;
    ptwt *= ptRandom;
    phiRandom = myRandom->Uniform(0,6.283);
    p1 = ptRandom*cos(phiRandom);
    p2 = ptRandom*sin(phiRandom);
    yRandom = myRandom->Gaus(0.,1.9);
    p3 = sqrt(MD*MD+ptRandom*ptRandom)*(exp(yRandom)-exp(-yRandom))/2.0;
    double e = sqrt(MD*MD+ptRandom*ptRandom+p3*p3);

    TLorentzVector D0p(TVector3(p1,p2,p3),e);
    leptonp=decay_kinematics(D0p,myatrest);

    if(fabs(D0p.Rapidity())<1.) 
      Dinvy->Fill(ptRandom,ptwt*wt/ptRandom);
    if(fabs(leptonp.PseudoRapidity())<1.) 
      einvy->Fill(leptonp.Perp(),ptwt*wt/leptonp.Perp());
  }
  
  einvy->SetStats(kFALSE);
  einvy->SetLineStyle(7);
  einvy->SetLineColor(kBlack);
  einvy->Draw(); 
  Dinvy->SetLineStyle(1);
  Dinvy->SetLineColor(kBlack);
  Dinvy->Draw("SAME");

  legend->AddEntry(einvy,"e inv. yield D^{0}","l");
  legend->AddEntry(Dinvy,"D^{0} inv. yield","l");
}
