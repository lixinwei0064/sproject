#include "TLorentzVector.h"
#include "TRandom.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "charm_decay.h"

TLorentzVector myBoost(TLorentzVector parent,TLorentzVector daughter){
       double mass               = parent.M();
    TVector3 eta = (-1./mass)*parent.Vect();            // gamma*beta
    double gamma              = fabs(parent.E())/mass;
    //    cout<<daughter.vect()<<" "<<eta<<endl;
    TVector3 pl  = ((daughter.Vect()*eta)/(eta*eta))*eta;  // longitudinal momentum
    return TLorentzVector(daughter.Vect() + (gamma-1.)*pl - daughter.E()*eta,
			  gamma*daughter.E() - daughter.Vect()*eta);
}

TLorentzVector decay_kinematics(TLorentzVector D0p,TH1D *formfactor){
   double e2,px,py,pz,ptot2,phi,costheta2;
   TLorentzVector leptonp3;
   ptot2 = formfactor->GetRandom();
   costheta2 = myRandom->Uniform(-1.0,1.0);
   phi=myRandom->Uniform(0.,2*TMath::Pi());
   pz = ptot2*costheta2;  // theta lepton relative to thetav see p482 PRD 54 (1996)
   px = ptot2*sqrt(1.0-costheta2*costheta2)*cos(phi);
   py = ptot2*sqrt(1.0-costheta2*costheta2)*sin(phi);
   e2 = sqrt(ptot2*ptot2+Me*Me);
   TVector3 l3v(px,py,pz);
   TLorentzVector leptonp(l3v,e2);
   leptonp3 = myBoost(D0p,leptonp);
   return leptonp3;

}

double get_value(TH1D* h, double x){
	if(x==0)
		return h->GetBinContent(0);
	if(x==9.95)
		return h->GetBinContent(200);
	int bin=(int)(x/0.05)+1;
	double b=x-0.05*(bin-1);
	if(b==0.025)
		return h->GetBinContent(bin);
	if(b<0.025){
		if(bin==1)
			return h->GetBinContent(1);
		return (h->GetBinContent(bin-1)*(0.025+b)+h->GetBinContent(bin)*(0.025-b))/0.05;
	}
	if(b>0.025){
		if(bin==199)
			return h->GetBinContent(199);
		return (h->GetBinContent(bin)*b+h->GetBinContent(bin+1)*(0.05-b))/0.05;
	}
	return 0;
}
	
void charm_Lc_Levy(const int nEvents = 1e9, const int cen = 0, const int seed = 32553){

  const double c2eBR = 0.045;
  const double BRLc  = 1;
  const double wt = 20./2./TMath::Pi()*c2eBR/BRLc/2.;
  char name[250]="10_60";
  char Name[250];

  myRandom->SetSeed(seed);

  cen==0?sprintf(Name,"hLc_spectra_10_60_mean"):sprintf(Name,"hLc_spectra_10_60_model%d",cen);
  TFile *fin = new TFile("cHadron.root");
  TH1D *Levy = (TH1D *)fin->Get(Name);

  double ptwtnorm = 9.95/Levy->Integral();

  TH1D *espectrum = new TH1D("espectrum","",800,0.,20.0);
  TH1D *einvy = new TH1D("einvy","",800,0.,20.0);
  TH1D *charmrandom = new TH1D("charmrandom","charm spectra",800,0.,20.0);
  TH1D *Linvy = new TH1D("Linvy","",800,0.,20.0);
  TH1D *charmy = new TH1D("charmy","charm rapidity spectra",200,-10.,10.);
  TH2D *Dey = new TH2D("Dey","electron and charm rapidity spectra",200,-10.,10.,200,-10.,10.);
  double p1,p2,p3,ptRandom,phiRandom,yRandom;
  TLorentzVector leptonp(0.,0.,0.,0.);
  decayMode = pseudoscalar;

  TFile *fformfactor = new TFile("charm_edecayforms_atrest.root");
  TH1D *myatrest = (TH1D*)fformfactor->Get("a15");

  for (int Event = 0;Event<nEvents;Event++){
    if(Event%20000000==0)cout<<"Events "<<Event<<endl;
    ptRandom = myRandom->Uniform(0.,9.95);
    double ptwt = get_value(Levy,ptRandom)*ptwtnorm;
    phiRandom = myRandom->Uniform(0,2*TMath::Pi());
    p1 = ptRandom*cos(phiRandom);
    p2 = ptRandom*sin(phiRandom);
    yRandom = myRandom->Gaus(0.,1.9);
    p3 = sqrt(MLc*MLc+ptRandom*ptRandom)*(exp(yRandom)-exp(-yRandom))/2.0;
    double e = sqrt(MLc*MLc+ptRandom*ptRandom+p3*p3);

    TLorentzVector D0p(TVector3(p1,p2,p3),e);
    charmy->Fill(D0p.Rapidity());
    if(fabs(D0p.Rapidity())<1.) {
      charmrandom->Fill(ptRandom,ptwt);
      Linvy->Fill(ptRandom,ptwt*wt/ptRandom);
    }
    leptonp=decay_kinematics(D0p,myatrest);
    Dey->Fill(leptonp.Rapidity(),D0p.Rapidity());

    if(fabs(leptonp.PseudoRapidity())<1.) {
      espectrum->Fill(leptonp.Perp(),ptwt); // dy*pT*dpT*Events
      einvy->Fill(leptonp.Perp(),ptwt*wt/leptonp.Perp());
    }
  }
  
  char outname[250];
  cen==0?sprintf(outname,"L2e_run11AuAu_%s_mean.root",name):sprintf(outname,"L2e_run11AuAu_%s_model%d.root",name,cen);
  TFile fout(outname,"recreate");
  espectrum->Write();
  einvy->Write();
  Linvy->Write();
  charmrandom->Write();
  charmy->Write();
  Dey->Write();
  fout.Close();
}
