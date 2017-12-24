#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TPad.h"
#include "TString.h"
const char name[2][250] = {"0_80","10_40"};
const char nameCent[4][250] = {"model1_diquark", "model2_threequark", "model3_Greco", "mean"};

void combine(const int cen = 0, const int cent = 0){
    TFile *D0f = new TFile(Form("d0v2/e_D0_%s.root",name[cen]));
    TFile *Lcf = new TFile(Form("lcv2/e_Lc.root"));
    TFile *ratiof = new TFile(Form("ratio/%s_%s.root",name[cen],nameCent[cent]));
    //v2 distribution
    TH1D *v2eD0 = (TH1D *)D0f->Get("v2D");
    TH1D *v2eLc = (TH1D *)Lcf->Get("v2D");
    //yield
    TH1D *yeD0 = (TH1D *)ratiof->Get("einvyD0");
    TH1D *yeLc = (TH1D *)ratiof->Get("einvyLc");
    //total
    TH1D *distri = new TH1D("combined","combined distribution of e. from D^{0} and #Lambda_{c};p_{T};v_{2}",50,0,5);
    Double_t D0e,D0y,Lce,Lcy,D0r,Lcr;
    for(int i = 1; i<=50; i++){    
        D0y = v2eD0->GetBinContent(i);   
        Lcy = v2eLc->GetBinContent(i);
        D0e = v2eD0->GetBinError(i);
        Lce = v2eLc->GetBinError(i);
        D0r = yeD0->GetBinContent(i);
        Lcr = yeLc->GetBinContent(i);

        distri->SetBinContent(i,(D0y*D0r+Lcy*Lcr)/(D0r+Lcr));
        distri->SetBinError(i,(D0e*D0r+Lce*Lcr)/(D0r+Lcr));
    }
    distri->SetMarkerColor(kBlue);
    distri->SetMarkerStyle(kFullCircle);
    distri->Draw("P");
    gPad->SaveAs(Form("%s_%s.png",name[cen],nameCent[cent]));
}
