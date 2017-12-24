#include "TH1.h"
#include "TFile.h"
TLatex* drawLatex(Double_t x, Double_t y, char* text, Int_t textFont, Double_t textSize, Int_t colorIndex)
{
    TLatex *latex = new TLatex(x,y,text);
    latex->SetNDC();
    latex->SetTextFont(textFont);
    latex->SetTextSize(textSize);
    latex->SetTextColor(colorIndex);
    latex->Draw("same");
    return latex;
}
void draw2(const int cen = 0, const int cent = 0){
	TCanvas* c1 = new TCanvas("c1"," ",10,10,1000,1000);
	char name[4][250] = {"0_80","0_10","10_40","40_80"};
	char Name[250];
	
	sprintf(Name, "D2e_run11AuAu_%s.root", name[cen]);
	TFile* fin1 = new TFile(Name);
	TH1D* einvy_D0 = (TH1D *)fin1->Get("einvy");
	TH1D* Dinvy = (TH1D *)fin1->Get("Dinvy");
	cent==0?sprintf(Name,"L2e_run11AuAu_10_60_mean.root"):sprintf(Name,"L2e_run11AuAu_10_60_model%d.root",cent);
	TFile* fin2 = new TFile(Name);
	TH1D* einvy_Lc = (TH1D *)fin2->Get("einvy");
	TH1D* Linvy = (TH1D *)fin2->Get("Linvy");
	
	Linvy->Scale(85);
	einvy_Lc->Scale(85);
	
	einvy_D0->SetLineColor(kBlue);
	einvy_D0->Draw("same");
	Dinvy->SetLineColor(kGreen);
	Dinvy->Draw("same");
	einvy_Lc->SetLineColor(kRed);
	einvy_Lc->Draw("same");
	Linvy->SetLineColor(kYellow);
	Linvy->Draw("same");
	c1->SetLogy();
	
	drawLatex(0.45,0.86,"--- e inv. yield of D_{0}",132,0.04,kBlue);
	drawLatex(0.45,0.81,"--- D_{0} inv. yield",132,0.04,kGreen);
	drawLatex(0.45,0.76,"--- e inv. yield of L_{c}",132,0.04,kRed);
	drawLatex(0.45,0.71,"--- L_{c} inv. yield",132,0.04,kYellow);
	
	cent==0?sprintf(Name,"invy2_mean.gif"):sprintf(Name,"invy2_model%d.gif",cent);
	c1->SaveAs(Name);
}