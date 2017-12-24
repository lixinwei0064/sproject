#include "TH1.h"
#include "TFile.h"
void calc(const int cen = 0, const int cent = 0){
	TCanvas* c1 = new TCanvas("c1"," ",10,10,1000,1000);
	char name[4][250] = {"0_80","0_10","10_40","40_80"};
	char Name[250];
	
	sprintf(Name, "D2e_run11AuAu_%s.root", name[cen]);
	TFile* fin1 = new TFile(Name);
	//TH1D* einvy_D0 = (TH1D *)fin1->Get("einvy");
	TH1D* Dinvy = (TH1D *)fin1->Get("Dinvy");
	cent==0?sprintf(Name,"L2e_run11AuAu_10_60_mean.root"):sprintf(Name,"L2e_run11AuAu_10_60_model%d.root",cent);
	TFile* fin2 = new TFile(Name);
	//TH1D* einvy_Lc = (TH1D *)fin2->Get("einvy");
	TH1D* Linvy = (TH1D *)fin2->Get("Linvy");
	TFile* fin3 = new TFile("cHadron.root");
	sprintf(Name, "hD0_spectra_%s", name[cen]);
	TH1D* D0_spectra = (TH1D *)fin3->Get(Name);
	cent==0?sprintf(Name,"hLc_spectra_10_60_mean"):sprintf(Name,"hLc_spectra_10_60_model%d",cent);
	TH1D* Lc_spectra = (TH1D *)fin3->Get(Name);

	TH1D* L1 = new TH1D(*Lc_spectra);
	L1->Sumw2();
	L1->Divide(Lc_spectra,D0_spectra,1,1);
	TH1D* L2 = new TH1D(*Linvy);
	L2->Sumw2();
	L2->Divide(Linvy,Dinvy,85,1);
	
	L1->Draw();
	L2->Draw("same");
	c1->SetLogy();
	cent==0?sprintf(Name,"calc_mean.gif"):sprintf(Name,"calc_model%d.gif",cent);
	c1->SaveAs(Name);
}