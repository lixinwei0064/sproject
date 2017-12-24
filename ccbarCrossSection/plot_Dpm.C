#include "myFunction.h"
void plot_Dpm() {
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
    
    /*const int ncent = 4;
    const char nameCent[ncent][250] = {"0_10", "10_40", "40_80", "0_80"};
    const char nameCent1[ncent][250] = {"0-10%", "10-40%", "40-80%", "0-80%"};
    const float Nbin[ncent] = {938.80170, 386.08527, 56.99229, 301.05848};*/
    const int ncent = 2;
    const char nameCent[ncent][250] = {"0_10", "10_40"};
    const char nameCent1[ncent][250] = {"0-10%", "10-40%"};
    const float Nbin[ncent] = {938.80170, 386.08527};
    const float scale[ncent] = {1/20., 1./5.};
    
    const char namePar[250] = "D^{#pm}";
    const char namePar1[250] = "Dpm";
    const char namePrint[250] = Form("AuAu 200 GeV");
    
    const int npt = 4;
    const double ptbin[npt+1] = {1.5, 2.5, 3.5, 5.0, 8.0};
    const float xmin = 0;
    const float xmax = 10.;
    
    // read
    TGraphErrors* gSpectraUp[ncent];
    TGraphErrors* gSpectraSys[ncent];
    TGraphErrors* gSpectraErr[ncent];
    TGraphErrors* gBand[ncent];
    TF1* fSpectra[ncent];
    /*TFile* fin = new TFile("data/Run14_Ds_spectra_w_sys.root");
    for(int icent=0; icent<ncent; icent++) {
        gSpectraUp[icent] = (TGraphErrors*)fin->Get(Form("Ds_spectra_shifted_%s",nameCent[icent]));
        gSpectraErr[icent] = (TGraphErrors*)fin->Get(Form("Ds_spectra_shifted_%s",nameCent[icent]));
        gSpectraSys[icent] = (TGraphErrors*)fin->Get(Form("Ds_spectra_shifted_sys_%s",nameCent[icent]));
        //gBand[icent] = (TGraphErrors*)fin->Get(Form("flevy_err_band_%s",nameCent[icent]));
        //fSpectra[icent] = (TF1*)fin->Get(Form("flevy_%s",nameCent[icent]));
    }
    fin->Close();*/
    TFile* fin = new TFile("data/Run14_Dpm_MyRaa_pT1.0.xiaolong.root");
    for(int icent=0; icent<ncent; icent++) {
        gSpectraUp[icent] = (TGraphErrors*)fin->Get(Form("gMyAuAuSpectra_%s",nameCent[icent]));
        gSpectraErr[icent] = (TGraphErrors*)fin->Get(Form("gMyAuAuSpectra_%s",nameCent[icent]));
        gSpectraSys[icent] = (TGraphErrors*)fin->Get(Form("gMyAuAuSpectraSys_%s",nameCent[icent]));
        
        gSpectraUp[icent]->Sort();
        const int nDel = 2;
        //int ndel = 0;
        for(int idel=0; idel<nDel; idel++) {
            gSpectraUp[icent]->RemovePoint(0);
            gSpectraErr[icent]->RemovePoint(0);
            gSpectraSys[icent]->RemovePoint(0);
        }
        //gBand[icent] = (TGraphErrors*)fin->Get(Form("flevy_err_band_%s",nameCent[icent]));
        //fSpectra[icent] = (TF1*)fin->Get(Form("flevy_%s",nameCent[icent]));
    }
    fin->Close();
    
    // sclae
    for(int icent=0; icent<ncent; icent++) {
        for(int ipt=0; ipt<gSpectraUp[icent]->GetN(); ipt++) {
            gSpectraUp[icent]->GetY()[ipt] *= scale[icent];
            gSpectraUp[icent]->GetEY()[ipt] *= scale[icent];
            
            gSpectraErr[icent]->GetY()[ipt] *= scale[icent];
            gSpectraErr[icent]->GetEY()[ipt] *= scale[icent];
            
            gSpectraSys[icent]->GetY()[ipt] *= scale[icent];
            gSpectraSys[icent]->GetEY()[ipt] *= scale[icent];
        }
    }
    
    // shift sys. err
    for(int icent=0; icent<ncent; icent++) {
        for(int ipt=0; ipt<gSpectraUp[icent]->GetN(); ipt++) {
            //float err = gSpectraErr[icent]->GetEY()[ipt];
            float sys = gSpectraSys[icent]->GetEY()[ipt];
            float mean = gSpectraErr[icent]->GetY()[ipt];
            gSpectraUp[icent]->GetY()[ipt] = mean+sys;
            
            //gSpectraErr[icent]->GetEX()[ipt] = 0.5*(ptbin[ipt+1] - ptbin[ipt]);
            //gSpectraUp[icent]->GetEX()[ipt] = 0.5*(ptbin[ipt+1] - ptbin[ipt]);
        }
    }
    
    // define histograph not pt shifted
    /*TH1F* hSpectraErr[ncent];
    TH1F* hSpectraUp[ncent];
    for(int icent=0; icent<ncent; icent++) {
        hSpectraErr[icent] = new TH1F(Form("hSpectraErr%i",icent),"",npt,ptbin);
        hSpectraUp[icent] = new TH1F(Form("hSpectraUp%i",icent),"",npt,ptbin);
        for(int ipt=0; ipt<npt; ipt++) {
            float num = gSpectraErr[icent]->GetY()[ipt];
            float err = gSpectraErr[icent]->GetEY()[ipt];
            hSpectraErr[icent]->SetBinContent(ipt+1, num);
            hSpectraErr[icent]->SetBinError(ipt+1, err);
            
            num = gSpectraUp[icent]->GetY()[ipt];
            err = gSpectraUp[icent]->GetEY()[ipt];
            hSpectraUp[icent]->SetBinContent(ipt+1, num);
            hSpectraUp[icent]->SetBinError(ipt+1, err);
        }
    }*/
    
    // fit
    //define fit function, now use levy function
    char funcString[200];
    char funcString_time_pt[200];
    double m0 = 1.8693;//D0-1.8645, D+/- - 1.8693;
    sprintf(funcString,"(1/(2*TMath::Pi()))*[0]*([2]-1)*([2]-2)/([2]*[1]*([2]*[1]+%f*([2]-2)))*TMath::Power(([2]*[1]+TMath::Sqrt(x[0]*x[0]+%f*%f)-%f)/([2]*[1]),-[2])",m0,m0,m0,m0); // dN/pTdpTdy
    sprintf(funcString_time_pt,"(1/(2*TMath::Pi()))*[0]*([2]-1)*([2]-2)/([2]*[1]*([2]*[1]+%f*([2]-2)))*TMath::Power(([2]*[1]+TMath::Sqrt(x[0]*x[0]+%f*%f)-%f)/([2]*[1]),-[2])*x[0]",m0,m0,m0,m0); // dN/dpTdy
    float fitR_lw = 0;
    float fitR_up = 10;
    const char fitOpt[100] = "INOR";
    TF1* flevyUp[ncent];
    TF1* fpowLaw[ncent];
    for(int icent=0; icent<ncent; icent++) {
        fSpectra[icent] = new TF1("flevy",funcString,0,10);
        fSpectra[icent]->SetParameters(6.25342e-01,3.34550e-01,1.86042e+01);
        gSpectraErr[icent]->Fit(fSpectra[icent], fitOpt, "", fitR_lw, fitR_up);
        gSpectraErr[icent]->Fit(fSpectra[icent], fitOpt, "", fitR_lw, fitR_up);
        gSpectraErr[icent]->Fit(fSpectra[icent], fitOpt, "", fitR_lw, fitR_up);
        gSpectraErr[icent]->Fit(fSpectra[icent], fitOpt, "", fitR_lw, fitR_up);
        //hSpectraErr[icent]->Fit(fSpectra[icent], fitOpt, "", fitR_lw, fitR_up);
        //hSpectraErr[icent]->Fit(fSpectra[icent], fitOpt, "", fitR_lw, fitR_up);
        //hSpectraErr[icent]->Fit(fSpectra[icent], fitOpt, "", fitR_lw, fitR_up);
        //hSpectraErr[icent]->Fit(fSpectra[icent], fitOpt, "", fitR_lw, fitR_up);
        
        // ==== Create a TGraphErrors to hold the confidence intervals
        const Int_t NCL = 200;
        //TGraphErrors *grint;
        gBand[icent] = new TGraphErrors(NCL);
        gBand[icent]->SetTitle("Fitted line with .68 conf. band");
        for (int i=0; i<NCL; i++) gBand[icent]->SetPoint(i, 0 + (float)i/20., 0);
        
        (TVirtualFitter::GetFitter())->GetConfidenceIntervals(gBand[icent],0.68);
        gBand[icent]->SetName(Form("flevy_err_band_%s",nameCent1[0]));
        
        fpowLaw[icent] = new TF1("fpow-law",PowLaw,0,10,3);
        fpowLaw[icent]->SetParameters(1.80265e-01,3.90637e+01,1.16876e+00);
        gSpectraErr[icent]->Fit(fpowLaw[icent], fitOpt, "", fitR_lw, fitR_up);
        gSpectraErr[icent]->Fit(fpowLaw[icent], fitOpt, "", fitR_lw, fitR_up);
        gSpectraErr[icent]->Fit(fpowLaw[icent], fitOpt, "", fitR_lw, fitR_up);
        gSpectraErr[icent]->Fit(fpowLaw[icent], fitOpt, "", fitR_lw, fitR_up);
        
        flevyUp[icent] = new TF1("flevy",funcString,0,10);
        flevyUp[icent]->SetParameters(6.25342e-01,3.34550e-01,1.86042e+01);
        gSpectraUp[icent]->Fit(flevyUp[icent], fitOpt, "", fitR_lw, fitR_up);
        gSpectraUp[icent]->Fit(flevyUp[icent], fitOpt, "", fitR_lw, fitR_up);
        gSpectraUp[icent]->Fit(flevyUp[icent], fitOpt, "", fitR_lw, fitR_up);
        gSpectraUp[icent]->Fit(flevyUp[icent], fitOpt, "", fitR_lw, fitR_up);
        //hSpectraUp[icent]->Fit(flevyUp[icent], fitOpt, "", fitR_lw, fitR_up);
        //hSpectraUp[icent]->Fit(flevyUp[icent], fitOpt, "", fitR_lw, fitR_up);
        //hSpectraUp[icent]->Fit(flevyUp[icent], fitOpt, "", fitR_lw, fitR_up);
        //hSpectraUp[icent]->Fit(flevyUp[icent], fitOpt, "", fitR_lw, fitR_up);
        
        // calculate band
        for(int ipt=0; ipt<gBand[icent]->GetN(); ipt++) {
            float pt = gBand[icent]->GetX()[ipt];
            float sys = fabs(flevyUp[icent]->Eval(pt) - fSpectra[icent]->Eval(pt));
            float err = gBand[icent]->GetEY()[ipt];
            //cout << err << "\t" << sys << "\t" << pt << "\t" << gBand[icent]->GetY()[ipt] << "\t" << sqrt(pow(err,2)+pow(sys,2)) << endl;
            float mean = (fpowLaw[icent]->Eval(pt) + fSpectra[icent]->Eval(pt))/2.;
            float extErr = fabs(fpowLaw[icent]->Eval(pt) - fSpectra[icent]->Eval(pt))/2./mean;
            float old = fSpectra[icent]->Eval(pt);
            gBand[icent]->GetY()[ipt] = mean;
            gBand[icent]->GetEY()[ipt] = sqrt((pow(err,2)+pow(sys,2))/old/old+pow(extErr,2))*mean;
        }
    }
    
    
    // calculte
    const float factor = 42.*1000; // ub
    const float twoPI = 2.*TMath::Pi();
    float xc[ncent], xc_up[ncent], xcErr[ncent];
    float xc_check[ncent], xcErr_check[ncent];
    //float xc_pw[ncent];
    for(int icent=0; icent<ncent; icent++) {
        // init
        xc[icent] = 0; xc_up[icent] = 0; xcErr[icent] = 0; xc_check[icent] = 0; xcErr_check[icent] = 0;
        //xc_pw[icent] = 0;
        
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
        cout << nameCent1[icent] << endl;
        cout << "from fit: xc = " << xc[icent] << ",\t err = " << xc_up[icent]-xc[icent] << endl;
        
        for(int i=0; i<npt; i++) {
            float pt = (ptbin[i] + ptbin[i+1])*0.5;
            float ptw = ptbin[i+1] - ptbin[i];
            float ymean = gSpectraErr[icent]->GetY()[i];
            float yerr = gSpectraErr[icent]->GetEY()[i];
            xc_check[icent] += ymean*twoPI*pt*ptw*factor / Nbin[icent];
            xcErr_check[icent] += pow(yerr*twoPI*pt*ptw*factor / Nbin[icent],2);
        }
        xcErr_check[icent] = sqrt(xcErr_check[icent]);
        cout << "from data point: xc = " << xc_check[icent] << ",\t err = " << xcErr_check[icent] << endl;
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
        
        fpowLaw[icent]->SetLineColor(kBlue);
        fpowLaw[icent]->SetLineWidth(2);
        fpowLaw[icent]->SetLineStyle(7);
        fpowLaw[icent]->Draw("same");
        
        /*flevyUp[icent]->SetLineColor(kBlue);
        flevyUp[icent]->SetLineWidth(2);
        flevyUp[icent]->SetLineStyle(7);
        flevyUp[icent]->Draw("same");*/
        
        gSpectraErr[icent]->SetMarkerStyle(kFullCircle);
        gSpectraErr[icent]->SetMarkerSize(1.5);
        gSpectraErr[icent]->SetMarkerColor(kBlack);
        gSpectraErr[icent]->SetLineColor(kBlack);
        gSpectraErr[icent]->SetLineWidth(2);
        gSpectraErr[icent]->Draw("psame");
        
        //draw systematic error
        const float sysw = 0.15;
        for(int i=0; i<gSpectraSys[icent]->GetN(); i++) {
            const float sysl = gSpectraSys[icent]->GetY()[i] * 0.05;
            TLine *llw = new TLine(gSpectraSys[icent]->GetX()[i]-sysw,gSpectraSys[icent]->GetY()[i]-gSpectraSys[icent]->GetEY()[i],gSpectraSys[icent]->GetX()[i]+sysw,gSpectraSys[icent]->GetY()[i]-gSpectraSys[icent]->GetEY()[i]);
            llw->SetLineWidth(2);
            llw->SetLineColor(kBlack);
            llw->Draw("same");
            TLine *lhi = new TLine(gSpectraSys[icent]->GetX()[i]-sysw,gSpectraSys[icent]->GetY()[i]+gSpectraSys[icent]->GetEY()[i],gSpectraSys[icent]->GetX()[i]+sysw,gSpectraSys[icent]->GetY()[i]+gSpectraSys[icent]->GetEY()[i]);
            lhi->SetLineWidth(2);
            lhi->SetLineColor(kBlack);
            lhi->Draw("same");
            TLine *lv = new TLine(gSpectraSys[icent]->GetX()[i]-sysw,gSpectraSys[icent]->GetY()[i]-gSpectraSys[icent]->GetEY()[i],gSpectraSys[icent]->GetX()[i]-sysw,gSpectraSys[icent]->GetY()[i]-gSpectraSys[icent]->GetEY()[i]+sysl);
            lv->SetLineWidth(2);
            lv->SetLineColor(kBlack);
            lv->Draw("same");
            TLine *lv = new TLine(gSpectraSys[icent]->GetX()[i]+sysw,gSpectraSys[icent]->GetY()[i]-gSpectraSys[icent]->GetEY()[i],gSpectraSys[icent]->GetX()[i]+sysw,gSpectraSys[icent]->GetY()[i]-gSpectraSys[icent]->GetEY()[i]+sysl);
            lv->SetLineWidth(2);
            lv->SetLineColor(kBlack);
            lv->Draw("same");
            TLine *lv = new TLine(gSpectraSys[icent]->GetX()[i]-sysw,gSpectraSys[icent]->GetY()[i]+gSpectraSys[icent]->GetEY()[i],gSpectraSys[icent]->GetX()[i]-sysw,gSpectraSys[icent]->GetY()[i]+gSpectraSys[icent]->GetEY()[i]-sysl);
            lv->SetLineWidth(2);
            lv->SetLineColor(kBlack);
            lv->Draw("same");
            TLine *lv = new TLine(gSpectraSys[icent]->GetX()[i]+sysw,gSpectraSys[icent]->GetY()[i]+gSpectraSys[icent]->GetEY()[i],gSpectraSys[icent]->GetX()[i]+sysw,gSpectraSys[icent]->GetY()[i]+gSpectraSys[icent]->GetEY()[i]-sysl);
            lv->SetLineWidth(2);
            lv->SetLineColor(kBlack);
            lv->Draw("same");
        }
        
        sprintf(name,"%s, %s",namePrint,nameCent1[icent]);
        drawLatex(0.52,0.89,name,132,0.04,1);
        sprintf(name,"%s",namePar);
        drawLatex(0.52,0.8,name,132,0.04,1);
        sprintf(name,"cross section = %.1f #pm %.1f #mub",xc[icent],xcErr[icent]);
        drawLatex(0.52,0.75,name,132,0.04,1);
        
        //const float legy_lw = 0.9 - 0.05*nMod;
        legend = new TLegend(0.52,0.6,0.9,0.7);
        legend->SetFillStyle(0);
        legend->SetFillColor(10);
        legend->SetBorderSize(0);
        legend->SetTextSize(0.04);
        legend->SetTextFont(132);
        legend->AddEntry(fSpectra[icent],"Levy","l");
        legend->AddEntry(fpowLaw[icent],"Power-law","l");
        //legend->AddEntry(gModel[2],nameMod[2],"f");
        legend->Draw();
        
        c1->SetLogy();
        sprintf(name,"%s/%s_%s.gif",dir,namePar1,nameCent[icent]);
        c1->SaveAs(name);
    }
    
    // write
    const int N = gBand[0]->GetN() - 1;
    const float binWidth = gBand[0]->GetX()[2] - gBand[0]->GetX()[1];
    const float y_lw = 0;
    const float y_up = binWidth*N;
    TH1F* hSpectra[ncent];
    for(int icent=0; icent<ncent; icent++) {
        sprintf(name, "h%s_spectra_%s", namePar1, nameCent[icent]);
        sprintf(title, ";p_{T} (GeV/c);d^{2}N/(N_{ev}dp_{T}dy) (GeV/c)^{-1}");
        hSpectra[icent] = new TH1F(name,title,N,y_lw,y_up);
        for(int i=0; i<gBand[icent]->GetN()-1; i++) {
            float ptw = gBand[icent]->GetX()[2] - gBand[icent]->GetX()[1];
            float pt = gBand[icent]->GetX()[i]+0.5*ptw;
            float ymean = 0.5*(gBand[icent]->GetY()[i]+gBand[icent]->GetY()[i+1]);
            float yup = ymean + 0.5*(gBand[icent]->GetEY()[i]+gBand[icent]->GetEY()[i+1]);
            //cout << pt << "\t" << ptw << "\t" << ymean << "\t" << yup << endl;
            //cout << gBand[icent]->GetN() << endl;
            ymean = ymean*twoPI*pt;
            yup = yup*twoPI*pt;
            hSpectra[icent]->SetBinContent(i+1,ymean);
            hSpectra[icent]->SetBinError(i+1,fabs(yup-ymean));
        }
    }
    TFile* fout = new TFile(Form("out/%s.root",namePar1), "RECREATE");
    for(int icent=0; icent<ncent; icent++) {
        gSpectraErr[icent]->SetTitle(";p_{T} (GeV/c);d^{2}N/(N_{ev}2#pip_{T}dp_{T}dy) (GeV/c)^{-2}");
        sprintf(name, "%s_%s_err", namePar1, nameCent[icent]);
        gSpectraErr[icent]->Write(name);
        gSpectraSys[icent]->SetTitle(";p_{T} (GeV/c);d^{2}N/(N_{ev}2#pip_{T}dp_{T}dy) (GeV/c)^{-2}");
        sprintf(name, "%s_%s_sys", namePar1, nameCent[icent]);
        gSpectraSys[icent]->Write(name);
        sprintf(name, "%s_%s_levy", namePar1, nameCent[icent]);
        fSpectra[icent]->Write(name);
        hSpectra[icent]->Write();
    }
    fout->Close();
}
