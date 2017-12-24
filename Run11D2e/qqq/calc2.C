#include "TH1.h"
#include "TF1.h"
#include "TFile.h"
void calc2(const int cen = 0, const int cent = 0){
	TCanvas* c1 = new TCanvas("c1"," ",10,10,1000,1000);
	char name[4][250] = {"0_80","0_10","10_40","40_80"};
	char Name[250];
	
	sprintf(Name, "D2e_run11AuAu_%s.root", name[cen]);
	TFile* fin1 = new TFile(Name);
	TH1D* Dinvy = (TH1D *)fin1->Get("Dinvy");
	TFile* fin2 = new TFile("cHadron.root");
	sprintf(Name, "D0_%s_levy", name[cen]);
	TF1* D0_spectra = (TF1 *)fin2->Get(Name);

	TH1D* D=new TH1D(*Dinvy);
	cout<<D0_spectra->Integral(0,20)/Dinvy->Integral(0,20)<<endl;
	cout<<Dinvy->GetBinWidth(1)<<endl;
	D->Scale((0.42*D0_spectra->Integral(0,20))/(Dinvy->Integral(0,20)*Dinvy->GetBinWidth(1)));
	
	D0_spectra->SetLineColor(kRed);
	D0_spectra->Draw();
	D->SetLineColor(kBlue);
	D->Draw("same");
	c1->SetLogy();
	cent==0?sprintf(Name,"calc2_mean.gif"):sprintf(Name,"calc2_model%d.gif",cent);
	c1->SaveAs(Name);
}