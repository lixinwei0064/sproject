//===============================
//Version: 1.0 
//Time: Tue Apr  5 22:56:14 EDT 2016 
//Author: Long Zhou 
//Discribe: first version used to draw input D0 v2 spectra QA 

#include <iostream>
#include <iomanip>
#include "TCanvas.h"
#include "TLegend.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TPDF.h"
#include "TGraph.h"
#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TVirtualFitter.h"

TH1D* graph2hist(TGraphErrors *gr);

TGraphErrors *gr;
// Histogram fit 
Double_t fun_tamu(double* x, double* pars) {
  double xx = x[0]; // use x[1] to get 2nd dimension, x[2] for 3rd ...
  // fit parameters, i.e. the histogram weights
  double amp   = pars[0];
  //Get content of the histograms for this point
  double y = gr->Eval(xx);

  return amp*y;
}

void SaveToPDF(TCanvas *c1,TString filename,Int_t pageFlag = 1)
{
  filename += ".pdf";
  if(pageFlag==0) 
    {
      cout<<"save first page to pdf file : "<<filename.Data()<<endl;
      c1->Print(Form("%s(",filename.Data()));
    }
  else if(pageFlag==2)
    {
      cout<<"save last page to pdf file : "<<filename.Data()<<endl;
      c1->Print(Form("%s)",filename.Data()));
      c1->Close();
    }
  else 
    {
      c1->Print(Form("%s",filename.Data()));
    }
  c1->Clear();
}

Int_t addShade(TGraph *shade, TGraph *up, TGraph *down)
{
  Int_t N1 = up->GetN();
  Int_t N2 = down->GetN();
  Int_t N  = shade->GetN();

  if (N1 != N2 || N != (N1 + N2)) {
    cout << "Input two graph do not have equal point. please check it ! " << endl;
    return 0;
  }

  Int_t NPoint = N1;

  for (int i = 0; i < NPoint; i++) {
    float xu = up->GetX()[i];
    float yu = up->GetY()[i];
    float xd = down->GetX()[NPoint - i - 1];
    float yd = down->GetY()[NPoint - i - 1];
    shade->SetPoint(i, xu, yu);
    shade->SetPoint(NPoint + i, xd, yd);
  }

  shade->SetFillColor(16);
  return 1;
}

void inputV2()
{
  // gROOT->Reset();

  TCanvas *c1 = new TCanvas("c1", "c1",0,0,800,650);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetEndErrorSize(0);
  gStyle->SetOptDate(0);
  c1->SetFillColor(10);
  c1->SetFillStyle(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(0);
  c1->SetFrameFillColor(10);
  c1->SetFrameFillStyle(0);
  c1->SetFrameBorderMode(0);
  //c1->SetLogy();
  c1->SetGridx(0);
  c1->SetGridy(0);
  c1->SetLeftMargin(0.16);
  c1->SetBottomMargin(0.18);
  c1->SetTopMargin(0.03);
  c1->SetRightMargin(0.03);

  double x1 = 0.0;
  double x2 = 8.0;
  double y1 = -0.03;
  double y2 = 0.3;
  TH1 *h0 = new TH1D("h0","",100,x1, x2);
  h0->SetName("hframe");
  h0->SetMinimum(y1);
  h0->SetMaximum(y2);
  h0->GetXaxis()->SetNdivisions(208);
  h0->GetXaxis()->CenterTitle();
  h0->GetXaxis()->SetTitle("Transverse Momentum p_{T} (GeV/c)");
  h0->GetXaxis()->SetTitleOffset(1.0);
  h0->GetXaxis()->SetTitleSize(0.07);
  h0->GetXaxis()->SetLabelOffset(0.01);
  h0->GetXaxis()->SetLabelSize(0.05);
  h0->GetXaxis()->SetLabelFont(42);
  h0->GetXaxis()->SetTitleFont(42);
  h0->GetYaxis()->SetNdivisions(505);
  h0->GetYaxis()->SetTitle("Anisotropy Parameter v_{2}");
  h0->GetYaxis()->SetTitleOffset(0.8);
  h0->GetYaxis()->SetTitleSize(0.07);
  h0->GetYaxis()->SetLabelOffset(0.015);
  h0->GetYaxis()->SetLabelSize(0.05);
  h0->GetYaxis()->SetLabelFont(42);
  h0->GetYaxis()->SetTitleFont(42);

  TString cenname[4] = {"0_80", "0_10", "10_40", "40_80"};
  // vnStat_0_80; -- ststistic error 
  // vnSyst_0_80; -- systematic error 
  TFile *finput = new TFile("D0_v2_star_new/Systematics_D0vn_SL16d_2016-10-14.ME.root");
  TString graph_name = "";
  TString sysname = "vnSyst";
  TString statname = "vnStat";
  TString outname = "star_d0_v2_qa";
  TString centrality = "";
  TGraphErrors *gr_data_sys;	// sys error only 
  TGraphErrors *gr_data_sta;	// stat error only
  TGraphErrors *gr_data_all = NULL;	// combine sys + stat error 

  const Int_t n_TAMU = 30;
  TGraphErrors *gr_tamu = new TGraphErrors(n_TAMU);	
  TGraph *gr_tamu_u = new TGraph(n_TAMU);	
  TGraph *gr_tamu_d = new TGraph(n_TAMU);	
  TGraph *gr_tamu_shade = new TGraph(2*n_TAMU);	
  Double_t x_TAMU[n_TAMU], y_TAMU[n_TAMU], ye_TAMU[n_TAMU];

  gr_tamu->SetName("gr_tamu");
  gr_tamu_u->SetName("gr_tamu_u");
  gr_tamu_d->SetName("gr_tamu_d");
  gr_tamu_shade->SetName("gr_tamu_shade");
  ifstream inData("TAMU_v2_0_80_09222015.txt");
  cout << " ======== TAMU ========= " << endl;
  for(int i=0;i<n_TAMU;i++) {
    double a, b, c;
    inData >> a >> b >> c;
    x_TAMU[i] = a;
    y_TAMU[i] = (b+c)/2.;
    ye_TAMU[i] = fabs(b-c)/2.;
    gr_tamu->SetPoint(i,x_TAMU[i],y_TAMU[i]);
    gr_tamu->SetPointError(i,0,ye_TAMU[i]);
    gr_tamu_u->SetPoint(i,a,b);
    gr_tamu_d->SetPoint(i,a,c);
    cout <<  x_TAMU[i] << "\t" << y_TAMU[i] << endl;
  }
  inData.close();

  if (!addShade(gr_tamu_shade, gr_tamu_u, gr_tamu_d)) gr_tamu_shade = NULL; 
  gr_tamu_shade->SetFillColor(kCyan);

  TH1D *h_tamu = graph2hist(gr_tamu);
  h_tamu->SetName("h_tamu");

  TFile *fout = new TFile("my_D0_v2.root","RECREATE");
  for(int nc=0; nc<4; nc++) {
    h0->Draw("c");

    graph_name = sysname + "_" + cenname[nc];
    gr_data_sys = (TGraphErrors *)finput->Get(graph_name.Data());
    gr_data_sys->SetName(graph_name.Data());
    gr_data_sys->SetMarkerStyle(20);
    gr_data_sys->SetMarkerColor(1);
    gr_data_sys->SetMarkerSize(1.8);
    gr_data_sys->SetLineWidth(2);
    // gr_data_sys->Draw("[]same");
    gr_data_sys->RemovePoint(0);

    graph_name = statname + "_" + cenname[nc];
    gr_data_sta = (TGraphErrors *)finput->Get(graph_name.Data());
    gr_data_sta->SetName(graph_name.Data());
    gr_data_sta->SetMarkerStyle(20);
    gr_data_sta->SetMarkerColor(1);
    gr_data_sta->SetMarkerSize(1.8);
    gr_data_sta->SetLineWidth(2);
    // gr_data_sta->Draw("psame");
    gr_data_sta->RemovePoint(0);

    // Combine sys + stat error 
    Int_t Nv2point = gr_data_sta->GetN();
    graph_name = "vnAll_" + cenname[nc];
    gr_data_all = new TGraphErrors(Nv2point);
    gr_data_all->SetName(graph_name.Data());
    gr_data_all->SetMarkerStyle(20);
    gr_data_all->SetMarkerColor(1);
    gr_data_all->SetMarkerSize(1.8);
    gr_data_all->SetLineWidth(2);

    for(int j=0;j<Nv2point;j++)
      {
	float pt = gr_data_sta->GetX()[j];
	float v2 = gr_data_sta->GetY()[j];
	float pt_err  = gr_data_sta->GetEX()[j];
	float sys_err = gr_data_sys->GetEY()[j];
	float sta_err = gr_data_sta->GetEY()[j];
	float err = sqrt(pow(sys_err,2) + pow(sta_err,2));

	gr_data_all->SetPoint(j, pt, v2);
	gr_data_all->SetPointError(j, pt_err, err);
      }

    TF1 *v2funMeson = new TF1("v2funMeson","[0]*[4]/(1+exp(-(x/[4]-[1])/[2]))-[3]*[4]",0.35,8.0);
    // v2funMeson->SetParameters(10.,0.35,0.2,3,2); // Defaults parameters 
    v2funMeson->SetName(Form("v2fun_%s",cenname[nc].Data()));
    v2funMeson->SetParameters(1.22096e-01,6.64587e-01,2.03237e-01,5.55109e-02,2.00000e+00);
    v2funMeson->FixParameter(4,2);
    v2funMeson->SetLineStyle(9);
    v2funMeson->SetLineWidth(4);
    v2funMeson->SetLineColor(4);
    // gr_data_all->Fit(v2funMeson,"0QEM");
    gr_data_all->Fit(v2funMeson,"0EM");
    // gr_data_sta->Fit(v2funMeson,"");
    // v2funMeson->Draw("same");
    gr_data_all->GetXaxis()->SetRangeUser(0,20.0);
    // gr_data_sta->Draw("samep");

    // ==== Create a TGraphErrors to hold the confidence intervals
    const Int_t NCL = 200;
    TGraphErrors *grint = new TGraphErrors(NCL);
    grint->SetTitle("Fitted line with .95 conf. band");
    for (int i=0; i<NCL; i++) grint->SetPoint(i, 1.0 + (float)i/10., 0);

    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint,0.68);
    Int_t  NPoint      = grint->GetN();
    TGraph *grintshade = new TGraph(2 * NCL);
    TGraph *grint_u    = new TGraph(NCL);
    TGraph *grint_d    = new TGraph(NCL);

    for(int i=0; i<NCL; i++) 
      {
    	grint_u->SetPoint(i,grint->GetX()[i],(grint->GetY()[i] + grint->GetEY()[i]));
    	grint_d->SetPoint(i,grint->GetX()[i],(grint->GetY()[i] - grint->GetEY()[i]));
      }

    if (!addShade(grintshade, grint_u, grint_d)) return;
    grintshade->SetFillColor(kCyan);
    // grintshade->Draw("fsame");
    v2funMeson->SetRange(1.0,8.0);
    v2funMeson->Draw("same");
    gr_data_all->Draw("samep");

    TF1 *fun = NULL;
    TGraphErrors *fun_tamu_fit_shade = NULL;
    TGraph *fun_tamu_fit_u = NULL;
    TGraph *fun_tamu_fit_d = NULL;
    if(nc == 0 ) 
      {
	h0->Draw("c");
	gr_tamu_shade->Draw("fsame");
	// gr_tamu->SetFillColor(kCyan);
	// gr_tamu->Draw("E3same");

	gr = new TGraphErrors(*gr_tamu);
	fun = new TF1("fun",fun_tamu,0,6,1);
	fun->SetName("fun_tamu_fit");
	fun->SetParNames("amp");
	fun->SetParameters(0,0.5);
	fun->SetLineColor(3);
	gr_data_all->Fit(fun,"N");
	//d0LSeD->Fit(fun,"N");
	// fun->Print();


	// TGraphErrors *fun_tamu_fit_shade = new TGraphErrors(NCL);
	fun_tamu_fit_shade = new TGraphErrors(61);
	fun_tamu_fit_shade->SetName("fun_tamu_fit_shade");
	fun_tamu_fit_shade->SetFillColor(kOrange);
	fun_tamu_fit_shade->SetTitle("Fitted line with .68 conf. band");
	for (int i=0; i<61; i++) fun_tamu_fit_shade->SetPoint(i, 0.0 + (float)i/10., 0);

	(TVirtualFitter::GetFitter())->GetConfidenceIntervals(fun_tamu_fit_shade,0.68);
	Int_t  NPoint  = fun_tamu_fit_shade->GetN();
	fun_tamu_fit_u = new TGraph(NPoint);
	fun_tamu_fit_d = new TGraph(NPoint);
	fun_tamu_fit_u->SetName("fun_tamu_fit_u");
	fun_tamu_fit_d->SetName("fun_tamu_fit_d");

	for(int i=0; i<NPoint; i++) 
	  {
	    fun_tamu_fit_u->SetPoint(i,fun_tamu_fit_shade->GetX()[i],(fun_tamu_fit_shade->GetY()[i] + fun_tamu_fit_shade->GetEY()[i]));
	    fun_tamu_fit_d->SetPoint(i,fun_tamu_fit_shade->GetX()[i],(fun_tamu_fit_shade->GetY()[i] - fun_tamu_fit_shade->GetEY()[i]));
	  }

	fun_tamu_fit_shade->Draw("E3same");
	gr_data_all->Draw("psame");
	fun->Draw("same");
      }

    // TF1 *v2fun  = new TF1("v2fun", "[0]*[4]/(1+exp(-(x/[4]-[1])/[2]))-[3]*[4]", 0.1, 20.);
    // v2fun->SetParameters(0.1, 0.35, 0.2, 0.03, 2);
    // v2fun->SetLineColor(2);
    // v2fun->Draw("same");

    TLatex *tex = new TLatex(x2*0.65, 0.29*0.9, "Au+Au@200 GeV");
    tex->SetTextFont(22);
    tex->SetTextSize(0.05);
    tex->Draw("same");
  
    centrality = cenname[nc].ReplaceAll("_","-");
    TLatex *tex2 = new TLatex(x2*0.75, 0.29*0.8, Form("%s%%",centrality.Data()));
    tex2->SetTextFont(22);
    tex2->SetTextSize(0.05);
    tex2->Draw("same");

    // TLegend *leg = new TLegend(0.5, 0.5, 0.7, 0.7);
    // leg->SetFillColor(10);
    // leg->SetLineStyle(4000);
    // leg->SetLineColor(10);
    // leg->SetLineWidth(0.);
    // leg->SetTextSize(0.035);
    // leg->AddEntry(gr_data_all, " D^{0} (P16id)", "p");
    // leg->AddEntry(grintshade,"95% conf.band", "f");

    TLegend *leg = new TLegend(0.17, 0.82, 0.4, 0.95);
    leg->SetFillColor(10);
    leg->SetLineStyle(4000);
    leg->SetLineColor(10);
    leg->SetLineWidth(0.);
    leg->SetTextSize(0.04);
    leg->AddEntry(gr_data_all, " D^{0} (P16id)", "p");
    // leg->AddEntry(grintshade,"68% conf.band", "f");
    if(nc == 0)
      {
	leg->AddEntry(gr_tamu_shade,"TAMU c-quark diff.", "f");
	leg->AddEntry(fun,"TAMU template fit", "l");
	leg->AddEntry(fun_tamu_fit_shade,"68% conf.band", "f");
      }
    leg->Draw();
  
    // c1->Update();
    // gSystem->Exec("mkdir -p plots");
    // c1->SaveAs(Form("plots/D0v2_input_%s.png",cenname[nc].Data()));
    // c1->SaveAs(Form("plots/D0v2_input_%s.eps",cenname[nc].Data()));
    // c1->Clear();

    if(nc==0) {c1->SaveAs("D0_v2_0_80_tamu.eps"); SaveToPDF(c1,outname,0);}
    else if(nc==3) SaveToPDF(c1,outname,2);
    else SaveToPDF(c1,outname);

    fout->WriteTObject(gr_data_sta);
    fout->WriteTObject(gr_data_sys);
    fout->WriteTObject(gr_data_all);
    fout->WriteTObject(v2funMeson);
    fout->WriteTObject(fun);
    fout->WriteTObject(fun_tamu_fit_shade);
    fout->WriteTObject(fun_tamu_fit_u);
    fout->WriteTObject(fun_tamu_fit_d);
  }
  h0->Write();
  h_tamu->Write();
  gr_tamu->Write();
  gr_tamu_u->Write();
  gr_tamu_d->Write();
  gr_tamu_shade->Write();
  // fun->Write();
  fout->Close();    
}

TH1D* graph2hist(TGraphErrors *gr)
{
  // Comments: the tgraph should have a uniform bins.
  Int_t NPoint = gr->GetN();
  float xmin = gr->GetX()[0];
  // float ex   = gr->GetEX()[0];
  float xmax = gr->GetX()[NPoint-1];
  float ex  = (xmax - xmin)/(NPoint - 1)/2;
  TH1D* hgraph = new TH1D("hgraph","",NPoint,xmin-ex,xmax+ex);
  for (int i = 0; i < gr->GetN(); i++) {
    float x = gr->GetX()[i];
    float y = gr->GetY()[i];
    float ey = gr->GetEY()[i];
    int bin = hgraph->GetXaxis()->FindBin(x);
    hgraph->SetBinContent(bin,y);
    hgraph->SetBinError(bin,ey);
  }
  return hgraph;
}
