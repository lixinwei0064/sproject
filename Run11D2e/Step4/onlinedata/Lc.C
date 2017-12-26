#include "TRandom.h"
#include "TLorentzVector.h"
#include "TF1.h"
#include "TH2.h"
#include "TFile.h"
#include "TPad.h"
#include "TGraphErrors.h"
#include "TString.h"
#include <iostream>
#include <fstream>
using namespace std;
const Double_t Me  = 0.005;
const Double_t ML = 2.286; 

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

void Lc(const int seed = 32553, const int nEvents = 100000){

  TRandom *myRandom = new TRandom(seed);

  TH2D *phiLc = new TH2D("phiLc","#phi distribution of #Lambda_{c};#phi;p_{T};count",18,0,6.283,60,0,6);
  TH2D *phie = new TH2D("phie","#phi distribution of electron;#phi;p_{T};count",18,0,6.283,60,0,6);
  /*TF1 *v2funBaryon = new TF1("v2funBaryon","[0]*[4]/(1+exp(-(x/[4]-[1])/[2]))-[3]*[4]",0.,5.);
  v2funBaryon->SetParameters(10.,0.35,0.2,3,3);*/
  ifstream input("inputdata.txt");
  double inputx[84] = {0};
  double inputy[84] = {0};
  for(int i=0;i<84;i++){
    input >> inputx[i];
    input >> inputy[i];
  }
  input.close();
  TGraph *v2funBaryon = new TGraph(84,inputx,inputy);
  TFile *fformfactor = new TFile("charm_edecayforms_atrest.root");
  TH1D *myatrest = (TH1D*)fformfactor->Get("a15");
  TF1 *f2 = new TF1("f2",Func2,0,6.283,1);

  Double_t p1,p2,p3,ptRandom,phiRandom,yRandom,E,v2;
  TLorentzVector Lcp,leptonp;

  for (int Event = 0;Event<nEvents;Event++){
    if(Event%1000000==0)cout<<"Events "<<Event<<endl;
    
    ptRandom = myRandom->Uniform(0,6);
    v2 = v2funBaryon->Eval(ptRandom);
    if(v2 < 0) v2 = 0;
    f2->SetParameter(0,v2);
    phiRandom = f2->GetRandom(0.,6.283);
    phiLc->Fill(phiRandom,ptRandom); 

    p1 = ptRandom*cos(phiRandom);
    p2 = ptRandom*sin(phiRandom);
    yRandom = myRandom->Gaus(0.,1.9);//?
    p3 = sqrt(ML*ML+ptRandom*ptRandom)*(exp(yRandom)-exp(-yRandom))/2.0;//?
    E = sqrt(ML*ML+ptRandom*ptRandom+p3*p3);  
    Lcp.SetPxPyPzE(p1,p2,p3,E);

    leptonp = decay_kinematics(Lcp,myatrest);
    if(fabs(leptonp.PseudoRapidity())<1.)  //?
      phie->Fill(phiRandom,leptonp.Perp());
  }

  TFile *phifile = new TFile(Form("output/e_Lc_1226.root"),"RECREATE");
  phiLc->Write();
  phie->Write();

  TF1 *func = new TF1("func","[1]*(1+2*[0]*TMath::Cos(2*x))",0,6.283);
  TH1D *v2D = new TH1D("v2D","distribution of v_{2} vs p_{T} of electron from #Lambda_{c} decay;p_{T};v_{2}",40,0,4.);
  TH1D *phicount;
  for(int i=1;i<=40;i++){
    phicount = phie->ProjectionX("",i,i);
    phicount->Fit(func);
    v2D->SetBinContent(i,func->GetParameter(0));
    v2D->SetBinError(i,func->GetParError(0));
  }
  //v2D->GetYaxis()->SetRangeUser(-0.2,0.2);
  v2D->Draw("E");
  v2D->Write();
  phifile->Close();
  gPad->SaveAs("v2_vs_pT_e_Lc_1226.png");
}

