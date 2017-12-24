#include "myFunction.h"
void plot_Ds() {
    globalSetting();
    char buf[250];
    char dir[250];
    char name[250];
    char CMD[250];
    char title[250];
    TLegend* legend;
    TH1F* h0;
    TCanvas* c1 = new TCanvas("c1", "A Canvas",10,10,1000,1000);
    setPad(c1);
    
    sprintf(dir,"pic");
    sprintf(CMD,"[ -d %s ] || mkdir -p %s",dir,dir);
    gSystem->Exec(CMD);
    
    const int ncent = 4;
    const char nameCent[ncent][250] = {"0_10", "10_40", "40_80", "0_80"};
    const char nameCent1[ncent][250] = {"0-10%", "10-40%", "40-80%", "0-80%"};
    const float Nbin[ncent] = {938.80170, 386.08527, 56.99229, 301.05848};
    
    const char namePar[250] = "D_{s}^{+}";
    const char namePar1[250] = "Ds";
    const char namePrint[250] = Form("AuAu 200 GeV");
    
    const int npt = 14;
    const float ptbin[npt+1] = {0, 0.5, 1., 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7, 8, 10};
    const float xmin = 0;
    const float xmax = 10.;
    
    // read
    TGraphErrors* gSpectra[ncent];
    TGraphErrors* gBand[ncent];
    TF1* fSpectra[ncent];
    TFile* fin = new TFile("data/Ds_spectra_shifted.root");
    for(int icent=0; icent<ncent; icent++) {
        gSpectra[icent] = (TGraphErrors*)fin->Get(Form("Ds_spectra_shifted_%s",nameCent[icent]));
        gBand[icent] = (TGraphErrors*)fin->Get(Form("flevy_err_band_%s",nameCent[icent]));
        fSpectra[icent] = (TF1*)fin->Get(Form("Levy_shifted_%s",nameCent[icent]));
    }
    
    // calculte
    const float factor = 42.*1000; // ub
    const float twoPI = 2.*TMath::Pi();
    float xc[ncent], xc_up[ncent], xcErr[ncent];
    float xc_check[ncent];
    for(int icent=0; icent<ncent; icent++) {
        // init
        xc[icent] = 0; xc_up[icent] = 0; xcErr[icent] = 0; xc_check[icent] = 0;
        
        for(int i=0; i<gBand[icent]->GetN()-1; i++) {
            float ptw = gBand[icent]->GetX()[2] - gBand[icent]->GetX()[1];
            float pt = gBand[icent]->GetX()[i]+0.5*ptw;
            float ymean = 0.5*(gBand[icent]->GetY()[i]+gBand[icent]->GetY()[i+1]);
            float yup = ymean + 0.5*(gBand[icent]->GetEY()[i]+gBand[icent]->GetEY()[i+1]);
            //cout << pt << "\t" << ptw << "\t" << ymean << "\t" << yup << endl;
            //cout << gBand[icent]->GetN() << endl;
            xc[icent] += ymean*twoPI*pt*ptw*factor / Nbin[icent];
            xc_up[icent] += yup*twoPI*pt*ptw*factor / Nbin[icent];
            xcErr[icent] = fabs(xc_up[icent] - xc[icent]);
        }
        cout << "from fit: xc = " << xc[icent] << ",\t err = " << xc_up[icent]-xc[icent] << endl;
    }
    
    // plot
    for(int icent=0; icent<ncent; icent++) {
        float ymin = 0.5*fSpectra[icent]->Eval(xmax);
        float ymax = 10.*fSpectra[icent]->Eval(xmin);
        h0 = new TH1F("","",1,xmin,xmax);
        h0->GetYaxis()->SetRangeUser(ymin,ymax);
        setHisto(h0,"","p_{T} (GeV/c)", "d^{2}N/(N_{ev}2#pip_{T}dp_{T}dy) (GeV/c)^{-2}");
        h0->Draw();
        
        //gBand[icent]->SetFillStyle(3001);
        //gBand[icent]->SetFillColor(kGreen+2);
        gBand[icent]->SetFillColorAlpha(kGreen,0.4);
        gBand[icent]->Draw("e3same");
        
        fSpectra[icent]->SetLineColor(kRed);
        fSpectra[icent]->SetLineWidth(2);
        fSpectra[icent]->SetLineStyle(7);
        fSpectra[icent]->Draw("same");
        
        gSpectra[icent]->SetMarkerStyle(kFullCircle);
        gSpectra[icent]->SetMarkerSize(1.5);
        gSpectra[icent]->SetMarkerColor(kBlack);
        gSpectra[icent]->SetLineColor(kBlack);
        gSpectra[icent]->SetLineWidth(2);
        gSpectra[icent]->Draw("psame");
        
        sprintf(name,"%s, %s",namePrint,nameCent1[icent]);
        drawLatex(0.52,0.89,name,132,0.04,1);
        sprintf(name,"%s",namePar);
        drawLatex(0.52,0.8,name,132,0.04,1);
        sprintf(name,"cross section = %.1f #pm %.1f #mub",xc[icent],xcErr[icent]);
        drawLatex(0.52,0.75,name,132,0.04,1);
        
        c1->SetLogy();
        sprintf(name,"%s/%s_%s.gif",dir,namePar1,nameCent[icent]);
        c1->SaveAs(name);
    }
}
