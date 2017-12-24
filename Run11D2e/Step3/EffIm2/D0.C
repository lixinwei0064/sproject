#include "TRandom.h"
#include "TLorentzVector.h"
#include "TF1.h"
#include "TH2.h"
#include "TFile.h"
#include "TPad.h"
#include "TGraphErrors.h"
#include "TString.h"
#include <iostream>
using namespace std;
const char name[2][250] = {"0_80","10_40"};
const Double_t Me  = 0.005;
const Double_t MD = 1.8646; 

Double_t Func2(Double_t *x, Double_t *par){
  Double_t v2 = par[0];
  return 1+2*v2*TMath::Cos(2*x[0]);
}

TLorentzVector myBoost(TLorentzVector parent,TLorentzVector daughter){
  double mass  = parent.M();
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
  costheta2 = gRandom->Uniform(-1.0,1.0);
  phi= gRandom->Uniform(0.,6.283);
  pz = ptot2*costheta2;  
  px = ptot2*sqrt(1.0-costheta2*costheta2)*cos(phi);
  py = ptot2*sqrt(1.0-costheta2*costheta2)*sin(phi);
  e2 = sqrt(ptot2*ptot2+Me*Me);
  TVector3 l3v(px,py,pz);
  TLorentzVector leptonp(l3v,e2);
  leptonp3 = myBoost(D0p,leptonp);
  return leptonp3;
}

void D0(const int cen = 0, const int seed = 32553, const int nEvents = 10000000){

  TRandom *myRandom = new TRandom(seed);

  TH2D *phiD0 = new TH2D("phiD0","#phi distribution of D^{0};#phi;p_{T};count",18,0,6.283,80,0,8);
  TH2D *phie = new TH2D("phie","#phi distribution of electron;#phi;p_{T};count",18,0,6.283,80,0,8);
  TFile *fina = new TFile("my_D0_v2.root");
  TGraphErrors *f1 = (TGraphErrors *)fina->Get(Form("v2fun_%s_error",name[cen]));
  TFile *fformfactor = new TFile("charm_edecayforms_atrest.root");
  TH1D *myatrest = (TH1D*)fformfactor->Get("a15");
  TF1 *f2 = new TF1("f2",Func2,0,6.283,1);

  Double_t p1,p2,p3,ptRandom,phiRandom,yRandom,E,v2;
  TLorentzVector D0p,leptonp;

  for (int Event = 0;Event<nEvents;Event++){
    if(Event%1000000==0)cout<<"Events "<<Event<<endl;
    ptRandom = myRandom->Uniform(0.,8.);
    v2 = f1->Eval(ptRandom);
    if(v2<0) v2=0;
    f2->SetParameter(0,v2);
    phiRandom = f2->GetRandom(0.,6.283);
    phiD0->Fill(phiRandom,ptRandom);

    p1 = ptRandom*cos(phiRandom);
    p2 = ptRandom*sin(phiRandom);
    yRandom = myRandom->Gaus(0.,1.9);
    p3 = sqrt(MD*MD+ptRandom*ptRandom)*(exp(yRandom)-exp(-yRandom))/2.0;
    E = sqrt(MD*MD+ptRandom*ptRandom+p3*p3);  
    D0p.SetPxPyPzE(p1,p2,p3,E);

    leptonp = decay_kinematics(D0p,myatrest);
    if(fabs(leptonp.PseudoRapidity())<1.)
      phie->Fill(phiRandom,leptonp.Perp());
  }

  TFile *phifile = new TFile(Form("output/e_D0_%s.root",name[cen]),"RECREATE");
  phiD0->Write();
  phie->Write();
  phiD0->Draw();
  gPad->SaveAs(Form("D0_%s.png",name[cen]));
  phie->Draw();
  gPad->SaveAs(Form("e_from_D0_%s.png",name[cen]));

  TH1D *v2D = new TH1D("v2D","v_{2} vs p_{T} of electron from D^{0} decay;p_{T};v_{2}",50,0,5);
  TF1 *func = new TF1("func","[1]*(1+2*[0]*TMath::Cos(2*x))",0,6.283);
  TH1D *phicount;
  for(int i=1;i<=50;i++){
    phicount = phie->ProjectionX("",i,i);
    phicount->Fit(func);
    v2D->SetBinContent(i,func->GetParameter(0));
    v2D->SetBinError(i,func->GetParError(0));
  }
  v2D->Draw("E");
  v2D->Write();
  gPad->SaveAs(Form("v2_vs_pT_%s_e_D0_newer.png",name[cen]));
}

