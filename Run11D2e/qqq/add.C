#include "TH1.h"
#include "TFile.h"
void add(const int cen = 0, const int cent = 0){
	char name[4][250] = {"0_80","0_10","10_40","40_80"};
	char Name[250];
	
	sprintf(Name, "D2e_run11AuAu_%s.root", name[cen]);
	TFile* fin1 = new TFile(Name);
	TH1D* einvy_D0 = (TH1D *)fin1->Get("einvy");
	TH1D* Dinvy = (TH1D *)fin1->Get("Dinvy");
	cent==0?sprintf(Name,"L2e_run11AuAu_10_60_mean.root"):sprintf(Name,"L2e_run11AuAu_10_60_model%d.root",cent);
	TFile* fin2 = new TFile(Name);
	TH1D* einvy_Lc = (TH1D *)fin2->Get("einvy");
	TFile* fin3 = new TFile("cHadron.root");
	sprintf(Name, "D0_%s_levy", name[cen]);
	TF1* D0_spectra = (TF1 *)fin3->Get(Name);
	
	TH1::SetDefaultSumw2();
	TH1D* einvy = new TH1D(*einvy_D0);
	einvy->Add(einvy_D0,einvy_Lc,1,85);
	einvy->Scale((0.42*D0_spectra->Integral(0,20))/(Dinvy->Integral(0,20)*Dinvy->GetBinWidth(1)));
	
	cent==0?sprintf(Name,"run11AuAu_mean.root"):sprintf(Name,"run11AuAu_model%d.root",cent);
	TFile* fout = new TFile(Name, "RECREATE");
	einvy->Write("einvy");
    fout->Close();
	
	TCanvas* c1 = new TCanvas("c1"," ",10,10,1000,1000);
	einvy->Draw();
	c1->SetLogy();
	c1->SaveAs("einvy.gif");
}