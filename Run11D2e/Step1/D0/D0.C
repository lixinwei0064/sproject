#include "TLorentzVector.h"
#include "TRandom.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TCanvas.h"
#include <iostream>
#include "charm_decay.h"
using namespace std;

Double_t LevyFcnPt(Double_t *x, Double_t *par)
{
  Double_t A  = par[0];
  Double_t n  = par[1];
  Double_t T  = par[2];
  Double_t m0 = par[3];
  Double_t mT = sqrt(x[0]*x[0]+m0*m0);

  Double_t a1 = A*(n-1)*(n-2);
  Double_t a2 = n*T*(n*T+m0*(n-2));
  Double_t a3 = pow(1+(mT-m0)/n/T,-n);

  return  a1/a2*a3 * x[0];
}

TLorentzVector myBoost(TLorentzVector parent,TLorentzVector daughter){
       double mass               = parent.M();
    TVector3 eta = (-1./mass)*parent.Vect();            // gamma*beta
    double gamma              = fabs(parent.E())/mass;
    //    cout<<daughter.vect()<<" "<<eta<<endl;
    TVector3 pl  = ((daughter.Vect()*eta)/(eta*eta))*eta;  // longitudinal momentum
    return TLorentzVector(daughter.Vect() + (gamma-1.)*pl - daughter.E()*eta,
			  gamma*daughter.E() - daughter.Vect()*eta);
}

}
TLorentzVector decay_kinematics(TLorentzVector D0p,TH1D *formfactor){
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
   leptonp3 = myBoost(D0p,leptonp);
   return leptonp3;

}
TLorentzVector decay_kinematics(TLorentzVector D0p,TF1 *formfactor){
   double e2,px,py,pz,ptot2,phi,costheta2;
   TLorentzVector leptonp3;
   ptot2 = formfactor->GetRandom(0.,1.);
   costheta2 = myRandom->Uniform(-1.0,1.0);
   phi=myRandom->Uniform(0.,6.283);
   pz = ptot2*costheta2;  // theta lepton relative to thetav see p482 PRD 54 (1996)
   px = ptot2*sqrt(1.0-costheta2*costheta2)*cos(phi);
   py = ptot2*sqrt(1.0-costheta2*costheta2)*sin(phi);
   e2 = sqrt(ptot2*ptot2+Me*Me);
   TVector3 l3v(px,py,pz);
   TLorentzVector leptonp(l3v,e2);
   leptonp3 = myBoost(D0p,leptonp);
   return leptonp3;

}
void D0(const int nEvents = 10000000, const int cen = 0, const int seed = 32553){

  const double c2eBR = 0.105;
  const double BRD0  = 0.565;
  const double wt = 20./2./3.14159265*c2eBR/BRD0/2.;
  const char nameCent[4][250] = {"0_80", "0_10", "10_40", "40_80"};
  myRandom->SetSeed(seed);

  TF1 *Levy;

  TCanvas* c1 = new TCanvas("c1", "A Canvas",10,10,1000,1000);

  
  TFile *fin = new TFile("cHadron.root");
  if(cen == 0) Levy = (TF1 *)fin->Get("D0_0_80_levy");
  if(cen == 1) Levy = (TF1 *)fin->Get("D0_0_10_levy");
  if(cen == 2) Levy = (TF1 *)fin->Get("D0_10_40_levy");
  if(cen == 3) Levy = (TF1 *)fin->Get("D0_40_80_levy");
  if(!Levy) exit(0);  
  Levy->SetName("Levy");

  float ptex = 3.5;

  TF1 *charmspectrum = new TF1("charmspectrum",LevyFcnPt,0.,20.,4);
  charmspectrum->SetParameters(Levy->GetParameters());

  float ptwtnorm = 20./Levy->Integral(0.,20.);

  TH1D *espectrum = new TH1D("espectrum","",800,0.,20.0);
  TH1D *einvy = new TH1D("einvy","",800,0.,20.0);
  //einvy->Sumw2();
  TH1D *charmrandom = new TH1D("charmrandom","charm spectra",800,0.,20.0);
  TH1D *Dinvy = new TH1D("Dinvy","",800,0.,20.0);
  //Dinvy->Sumw2();
  TH1D *charmy = new TH1D("charmy","charm rapidity spectra",200,-10.,10.);
  TH2D *Dey = new TH2D("Dey","electron and charm  rapidity spectra",200,-10.,10.,200,-10.,10.);
  double p1,p2,p3,ptRandom,phiRandom,yRandom;
  TLorentzVector leptonp(0.,0.,0.,0.);
  decayMode = pseudoscalar;

  TFile *fformfactor = new TFile("charm_edecayforms_atrest.root");
  TH1D *myatrest = (TH1D*)fformfactor->Get("a15");

  for (int Event = 0;Event<nEvents;Event++){
    if(Event%1000000==0)cout<<"Events "<<Event<<endl;
    //ptRandom = charmspectrum->GetRandom(0.0,20.0);
    ptRandom = myRandom->Uniform(0.,20.);
    //float ptwt = charmspectrum->Eval(ptRandom)*ptwtnorm;
    float ptwt = Levy->Eval(ptRandom)*ptwtnorm;
    //cout << exwt << endl;
    ptwt *= ptRandom;
    phiRandom = myRandom->Uniform(0,6.283);
    //    ptRandom = 0.001;//0.244; //
    p1 = ptRandom*cos(phiRandom);
    p2 = ptRandom*sin(phiRandom);
    yRandom = myRandom->Gaus(0.,1.9);
    //    yRandom= 0.0; //
    p3 = sqrt(MD*MD+ptRandom*ptRandom)*(exp(yRandom)-exp(-yRandom))/2.0;
    double e = sqrt(MD*MD+ptRandom*ptRandom+p3*p3);

    TLorentzVector D0p(TVector3(p1,p2,p3),e);
    //    double rapidity = 0.5*::log((e+p3)/(e-p3));
    charmy->Fill(D0p.Rapidity());
    if(fabs(D0p.Rapidity())<1.) {
      charmrandom->Fill(ptRandom,ptwt);
      Dinvy->Fill(ptRandom,ptwt*wt/ptRandom);
    }
    
    //TLorentzVector  D0p(TVector3(0.,0.,0.),MD);
    //    leptonp=decay_kinematics(D0p,MD,Mkstar,Me);
    leptonp=decay_kinematics(D0p,myatrest);
    Dey->Fill(leptonp.Rapidity(),D0p.Rapidity());

    if(fabs(leptonp.PseudoRapidity())<1.) {
      espectrum->Fill(leptonp.Perp(),ptwt); // dy*pT*dpT*Events
      einvy->Fill(leptonp.Perp(),ptwt*wt/leptonp.Perp());
    }

    //      espectrum->Fill(leptonp.P()); // dy*pT*dpT*Events
  }
  
  char outname[100];
  sprintf(outname,"D0_%s_levy.root",nameCent[cen]);
  TFile f1(outname,"recreate");
  espectrum->Draw();
  sprintf(outname,"espectrum_%s.png",nameCent[cen]);
  c1->SaveAs(outname);
  espectrum->Write();
  einvy->Draw();
  sprintf(outname,"einvy_%s.png",nameCent[cen]);
  c1->SaveAs(outname);
  einvy->Write();
  charmspectrum->Draw();
  sprintf(outname,"charmspectrum_%s.png",nameCent[cen]);
  c1->SaveAs(outname);
  charmspectrum->Write();
  Levy->Draw();
  sprintf(outname,"Levy_%s.png",nameCent[cen]);
  c1->SaveAs(outname);
  Levy->Write();
  Dinvy->Draw();
  sprintf(outname,"Dinvy_%s.png",nameCent[cen]);
  c1->SaveAs(outname);
  Dinvy->Write();
  charmrandom->Draw();
  sprintf(outname,"charmrandom_%s.png",nameCent[cen]);
  c1->SaveAs(outname);
  charmrandom->Write();
  charmy->Draw();
  sprintf(outname,"charmy_%s.png",nameCent[cen]);
  c1->SaveAs(outname);
  charmy->Write();
  Dey->Draw();
  sprintf(outname,"Dey_%s.png",nameCent[cen]);
  c1->SaveAs(outname);
  Dey->Write();
  f1.Close();
}
