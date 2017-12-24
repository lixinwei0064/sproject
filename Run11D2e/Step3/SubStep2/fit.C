/*#include <TFile.h>
#include <TH2.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TPad.h>
#include <iostream>*/

void fit(){
    TFile *fi = new TFile("try_e_D0_0_80.root");
    TH2D *phie = (TH2D *)fi->Get("phie");
    TH1D *phicount;
    TF1 *func = new TF1("func","[1]*(1+2*[0]*TMath::Cos(2*x))",0,6.283);
    /*TGraphErrors *v2D = new TGraphErrors(60); 
    v2D->SetNameTitle("v2D","distribution of v_{2} vs p_{T} of electron from D^{0} decay;p_{T};v_{2}");*/
    TGraph *v2D = new TGraph(60); 
    v2D->SetNameTitle("v2D","distribution of v_{2} vs p_{T} of electron from D^{0} decay;p_{T};v_{2}");
    for(int i=1;i<=60;i++){
        phicount = phie->ProjectionX("",i,i);
        phicount->Fit(func);
        v2D->SetPoint(i-1,0.1*i,func->GetParameter(0));
        /*v2D->SetPointError(i-1,0,func->GetParError(0));*/
    }
    v2D->Draw("AC");
    gPad->SaveAs("try_v2_vs_pT_0_80_e_D0.png");
}