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

Double_t powerlaw(Double_t *x, Double_t *par)
{
  double pT = x[0];
  double dNdy = par[0];
  double meanPt = par[1];
  double n = par[2];
  double p0 = meanPt * (n-3.)/2.;
  double A = dNdy*4.*(n-1)*(n-2)/(n-3)/(n-3)/meanPt/meanPt;
  return pT*A*TMath::Power(1+pT/p0, -n);
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

TLorentzVector decay_kinematics(TLorentzVector D0p,double cmass,double hmass,double lmass){
   double q2, e1,e2,px,py,pz,ptot1,ptot2,costheta1,phi,costheta2;
   double q2min = lmass*lmass,q2max = (cmass-hmass)*(cmass-hmass);
   bool q2random = true; int wasteloops=0;
   TLorentzVector leptonp2,leptonp3;
   double mpole2 = 1.89*1.89;
   double mpole2A = 2.5*2.5,mpole2V=2.1*2.1,Rv=1.62,R2=0.83;
   while(q2random){
     q2=myRandom->Uniform(q2min,q2max);
     e1 = (cmass*cmass-hmass*hmass+q2)/2./cmass; 
     double qmass = sqrt(q2);
     ptot1 = sqrt((cmass*cmass-(hmass+qmass)*(hmass+qmass))*(cmass*cmass-(hmass-qmass)*(hmass-qmass)))/2./cmass;
     costheta2 = myRandom->Uniform(-1.0,1.0);
     //     if (hmass<0.5){  //pseudoscalar 
     if (decayMode==pseudoscalar){  //pseudoscalar 
       double phasecorr = pow(ptot1,3)/pow((1.-q2/mpole2),2);
       wasteloops++;
       if(phasecorr>1.0)cout<<"phasecorr>1 something is wrong!!!"<<endl;
       if(phasecorr<myRandom->Uniform(0.,1.)) {
  //                if(wasteloops>10) cout<<"loops "<<wasteloops<<" "<<q2<<" < "<<phasecorr<<endl; 
	 continue;}
       wasteloops=0;
       q2random = false;
       phi=myRandom->Uniform(0.,6.283);
       costheta1 = myRandom->Uniform(-1.0,1.0);
       pz = ptot1*costheta1;
       px = ptot1*sqrt(1.0-costheta1*costheta1)*cos(phi);
       py = ptot1*sqrt(1.0-costheta1*costheta1)*sin(phi);
       TVector3 q3v(px,py,pz);
       TLorentzVector q2p(q3v,e1);
      //       cout<<sqrt(q2)<<"  "<<q2p.M()<<endl;
       e2 = (q2+lmass*lmass)/2./qmass; //neutrino mass = 0 
       ptot2 = (q2-lmass*lmass)/2./qmass;
       phi=myRandom->Uniform(0.,6.283);
       pz = ptot2*costheta2;  // theta lepton relative to thetav see p482 PRD 54 (1996)
       px = ptot2*sqrt(1.0-costheta2*costheta2)*cos(phi);
       py = ptot2*sqrt(1.0-costheta2*costheta2)*sin(phi);
       TVector3 l3v(px,py,pz);
       TLorentzVector leptonp(l3v,e2);
       leptonp2 = myBoost(q2p,leptonp);
     }
     else { //vector meson 
       double phasecorr = ptot1*q2*((
        pow((1-costheta2),2)*pow((cmass+hmass)/(1-q2/mpole2A)-2.*cmass*ptot1/(cmass+hmass)*Rv/(1-q2/mpole2V),2)+
        pow((1+costheta2),2)*pow((cmass+hmass)/(1-q2/mpole2A)+2.*cmass*ptot1/(cmass+hmass)*Rv/(1-q2/mpole2V),2),2)*4./3.+
        8./3.*(1-costheta2*costheta2)*pow(1./(2.*hmass*qmass)*((cmass*cmass-hmass*hmass-q2)*(cmass+hmass)/(1-q2/mpole2A)-(4.*cmass*cmass*ptot1*ptot1/(cmass+hmass)*R2/(1-q2/mpole2A))),2));
       wasteloops++;
       if(phasecorr>20.0)cout<<"phasecorr>1 something is wrong!!!"<<endl;
       if(phasecorr<myRandom->Uniform(0.,20.)) {
	 //       if(wasteloops>10) cout<<"loops "<<wasteloops<<" "<<q2<<" < "<<phasecorr<<endl; 
	 continue;}
       wasteloops=0;
       q2random = false;
       phi=myRandom->Uniform(0.,6.283);
       costheta1 = myRandom->Uniform(-1.0,1.0);
       pz = ptot1*costheta1;
       px = ptot1*sqrt(1.0-costheta1*costheta1)*cos(phi);
       py = ptot1*sqrt(1.0-costheta1*costheta1)*sin(phi);
       //     double phasecorr;
       TVector3 q3v(px,py,pz);
       TLorentzVector q2p(q3v,e1);
       e2 = (q2+lmass*lmass)/2./qmass; //neutrino mass = 0 
       ptot2 = (q2-lmass*lmass)/2./qmass;
       phi=myRandom->Uniform(0.,6.283);
       pz = -ptot2*costheta2;
       px = ptot2*sqrt(1.0-costheta2*costheta2)*cos(phi);
       py = ptot2*sqrt(1.0-costheta2*costheta2)*sin(phi);
       TVector3 l3v(px,py,pz);
       TVector3 qdirection = q3v.Unit();
       l3v.RotateUz(qdirection);
       TLorentzVector leptonp(l3v,e2);
       leptonp2 = myBoost(q2p,leptonp);
     }
     //         leptonp2 = leptonp.boost(q2p);
     leptonp3 = myBoost(D0p,leptonp2);
   }
   return leptonp3;

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
void Lc(const int nEvents = 10000000, const int cen = 0, const int seed = 32553){

  const double c2eBR = 0.105;
  const double BRD0  = 0.565;
  const double wt = 20./2./3.14159265*c2eBR/BRD0/2.;
  const char nameCent[4][250] = {"model1_diquark", "model2_threequark", "model3_Greco", "mean"};

  myRandom->SetSeed(seed);

  TH1F *Histo;

  TCanvas* c1 = new TCanvas("c1", "A Canvas",10,10,1000,1000);

  TFile *fin = new TFile("cHadron.root");
  if(cen == 0) Histo = (TH1F *)fin->Get("hLc_spectra_10_60_model1");
  if(cen == 1) Histo = (TH1F *)fin->Get("hLc_spectra_10_60_model2");
  if(cen == 2) Histo = (TH1F *)fin->Get("hLc_spectra_10_60_model3");
  if(cen == 3) Histo = (TH1F *)fin->Get("hLc_spectra_10_60_mean");
  if(!Histo) exit(0);
  Histo->SetName("Histo");

  float ptex = 3.5;
  TF1 *funPL = new TF1("funPL",powerlaw,3.,20.,3);
  if(cen == 0) { funPL->SetParameters(4.18514e+02,1.39882e-01,8.31462); }
  if(cen == 1) { funPL->SetParameters(2.48162e+03*20.,1.03483e-01,7.84250); ptex = 5.;}
  if(cen == 3) { funPL->SetParameters(1.09053e+02,9.64796e-02,7.34148); }

  float ptwtnorm = 20./Histo->Integral(0.,20.);

  TH1D *espectrum = new TH1D("espectrum","",800,0.,20.0);
  TH1D *einvy = new TH1D("einvy","",800,0.,20.0);
  //einvy->Sumw2();
  TH1D *charmrandom = new TH1D("charmrandom","charm spectra",800,0.,20.0);
  TH1D *Lcinvy = new TH1D("Lcinvy","",800,0.,20.0);
  //Lcinvy->Sumw2();
  TH1D *charmy = new TH1D("charmy","charm rapidity spectra",200,-10.,10.);
  TH2D *Lcey = new TH2D("Lcey","electron and charm  rapidity spectra",200,-10.,10.,200,-10.,10.);
  double p1,p2,p3,ptRandom,phiRandom,yRandom;
  TLorentzVector leptonp(0.,0.,0.,0.);
  decayMode = pseudoscalar;

  TFile *fformfactor = new TFile("charm_edecayforms_atrest.root");
  TH1D *pythia = (TH1D*)fformfactor->Get("a20");
  TH1D *myatrest = (TH1D*)fformfactor->Get("a15");
  TH1D *mypsi = (TH1D*)fformfactor->Get("a3");
  TF1 *vogt = (TF1*)fformfactor->Get("vogt");

  for (int Event = 0;Event<nEvents;Event++){
    if(Event%1000000==0)cout<<"Events "<<Event<<endl;
    //ptRandom = charmspectrum->GetRandom(0.0,20.0);
    ptRandom = myRandom->Uniform(0.,20.);
    int binNum;
    double binLowEdge, binUpEdge, binWidth, value;
    TAxis *Axis = Histo->GetXaxis();
    binNum = Axis->FindBin(ptRandom);
    binLowEdge = Axis->GetBinLowEdge(binNum);
    binUpEdge = Axis->GetBinUpEdge(binNum);
    binWidth = Axis->GetBinWidth(binNum);
    value = Histo->GetBinContent(binNum)+(ptRandom-binLowEdge)/binWidth*(Histo->GetBinContent(binNum+1)-Histo->GetBinContent(binNum));
    float ptwt = value*ptwtnorm;
    //cout << exwt << endl;
    phiRandom = myRandom->Uniform(0,6.283);
    //    ptRandom = 0.001;//0.244; //
    p1 = ptRandom*cos(phiRandom);
    p2 = ptRandom*sin(phiRandom);
    yRandom = myRandom->Gaus(0.,1.9);
    //    yRandom= 0.0; //
    p3 = sqrt(ML*ML+ptRandom*ptRandom)*(exp(yRandom)-exp(-yRandom))/2.0;
    double e = sqrt(ML*ML+ptRandom*ptRandom+p3*p3);

    TLorentzVector D0p(TVector3(p1,p2,p3),e);
    //    double rapidity = 0.5*::log((e+p3)/(e-p3));
    charmy->Fill(D0p.Rapidity());
    if(fabs(D0p.Rapidity())<1.) {
      charmrandom->Fill(ptRandom,ptwt);
      Lcinvy->Fill(ptRandom,ptwt*wt/ptRandom);
    }
    
    //TLorentzVector  D0p(TVector3(0.,0.,0.),ML);
    //    leptonp=decay_kinematics(D0p,ML,Mkstar,Me);
    leptonp=decay_kinematics(D0p,myatrest);
    Lcey->Fill(leptonp.Rapidity(),D0p.Rapidity());

    if(fabs(leptonp.PseudoRapidity())<1.) {
      espectrum->Fill(leptonp.Perp(),ptwt); // dy*pT*dpT*Events
      einvy->Fill(leptonp.Perp(),ptwt*wt/leptonp.Perp());
    }

    //      espectrum->Fill(leptonp.P()); // dy*pT*dpT*Events
  }
  
  char outname[100];
  sprintf(outname,"hLc_spectra_10_60_%s.root",nameCent[cen]);
  TFile f1(outname,"recreate");
  espectrum->Draw();
  sprintf(outname,"espectrum_%s.png",nameCent[cen]);
  c1->SaveAs(outname);
  espectrum->Write();
  einvy->Draw();
  sprintf(outname,"einvy_%s.png",nameCent[cen]);
  c1->SaveAs(outname);
  einvy->Write();
  Histo->Draw();
  sprintf(outname,"charmspectrum_%s.png",nameCent[cen]);
  c1->SaveAs(outname);
  Histo->Write();
  Lcinvy->Draw();
  sprintf(outname,"Lcinvy_%s.png",nameCent[cen]);
  c1->SaveAs(outname);
  Lcinvy->Write();
  charmrandom->Draw();
  sprintf(outname,"charmrandom_%s.png",nameCent[cen]);
  c1->SaveAs(outname);
  charmrandom->Write();
  charmy->Draw();
  sprintf(outname,"charmy_%s.png",nameCent[cen]);
  c1->SaveAs(outname);
  charmy->Write();
  Lcey->Draw();
  sprintf(outname,"Lcey_%s.png",nameCent[cen]);
  c1->SaveAs(outname);
  Lcey->Write();
  f1.Close();
}
