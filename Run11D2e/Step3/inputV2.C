#include <iostream>
#include <iomanip>
#include "TCanvas.h"
#include "TLegend.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TVirtualFitter.h"

void inputV2()
{
  TString cenname[2] = {"0_80","10_40"};
  // vnStat_0_80; -- statistic error 
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

  TFile *fout = new TFile("my_D0_v2.root","RECREATE");
  for(int nc=0; nc<2; nc++) {
    graph_name = sysname + "_" + cenname[nc];
    gr_data_sys = (TGraphErrors *)finput->Get(graph_name.Data());
    gr_data_sys->SetName(graph_name.Data());
    gr_data_sys->RemovePoint(0);

    graph_name = statname + "_" + cenname[nc];
    gr_data_sta = (TGraphErrors *)finput->Get(graph_name.Data());
    gr_data_sta->SetName(graph_name.Data());
    gr_data_sta->RemovePoint(0);

    // Combine sys + stat error 
    Int_t Nv2point = gr_data_sta->GetN();
    graph_name = "vnAll_" + cenname[nc];
    gr_data_all = new TGraphErrors(Nv2point);
    gr_data_all->SetName(graph_name.Data());

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
    // new V2 function, Add one item to let v2 decress at high pt. 
    // -------------------------------------------------------------
    TF1 *v2funMeson = new TF1("v2funMeson","[0]*[4]/(1+exp(-(x/[4]-[1])/[2]))-([3]+[5]*x)*[4]",0.35,8.0);
    // TF1 *v2funMeson = new TF1("v2funMeson","[0]*[4]/(1+exp(-(x/[4]-[1])/[2]))-[3]*[4]",0.35,8.0);
    v2funMeson->SetName(Form("v2fun_%s",cenname[nc].Data()));
    // v2funMeson->SetParameters(1.33246e-01,6.17872e-01,2.08257e-01,6.81131e-02,2.00000e+00);
    // v2funMeson->SetParameters(1.33246e-01,6.17872e-01,2.08257e-01,6.81131e-02,2.00000e+00,1);
    v2funMeson->SetParameters(-1.60104e+03,-4.10613e+00,-4.92650e-01,-1.22096e-01,2.00000e+00,1.24718e-02);
    v2funMeson->FixParameter(4,2);

    // gr_data_all->Fit(v2funMeson,"0QEM");
    gr_data_all->Fit(v2funMeson,"0EM");
    // gr_data_sta->Fit(v2funMeson,"");
    // v2funMeson->Draw("same");
    gr_data_all->GetXaxis()->SetRangeUser(0,20.0);
    // gr_data_sta->Draw("samep");

    // ==== Create a TGraphErrors to hold the confidence intervals
    const Int_t NCL = 200;
    TGraphErrors *grint = new TGraphErrors(NCL);
    grint->SetName(Form("v2fun_%s_error",cenname[nc].Data()));
    grint->SetTitle("Fitted line with .95 conf. band");
    for (int i=0; i<NCL; i++) grint->SetPoint(i, 1.0 + (float)i/10., 0);

    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(grint,0.68);
    Int_t  NPoint      = grint->GetN();
    TGraph *grintshade = new TGraph(2 * NCL);
    TGraph *grint_u    = new TGraph(NCL);
    TGraph *grint_d    = new TGraph(NCL);
    grint_u->SetName(Form("v2fun_%s_u",cenname[nc].Data()));
    grint_d->SetName(Form("v2fun_%s_d",cenname[nc].Data()));

    fout->WriteTObject(gr_data_sta);
    fout->WriteTObject(gr_data_sys);
    fout->WriteTObject(gr_data_all);
    fout->WriteTObject(v2funMeson);
    fout->WriteTObject(grint);
    fout->WriteTObject(grint_u);
    fout->WriteTObject(grint_d);
  }
  fout->Close();    
}
