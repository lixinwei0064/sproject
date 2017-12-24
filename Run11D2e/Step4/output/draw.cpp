void draw(){
    TFile fi("e_Lc.root");
    TH2D *phiLc;
    fi.GetObject("phiLc",phiLc);
    TF1 *func = new TF1("func","[1]*(1+2*[0]*TMath::Cos(2*x))",0,6.283);
    TH1D *v2D = new TH1D("v2D","distribution of v_{2} vs p_{T} of #Lambda_{c};p_{T};v_{2}",50,0.,5.);
    TH1D *phicount;
    for(int i=1;i<=50;i++){
        phicount = phiLc->ProjectionX("",i,i);
        phicount->Fit(func);
        v2D->SetBinContent(i,func->GetParameter(0));
        v2D->SetBinError(i,func->GetParError(0));
    }
    v2D->Draw("E");
    gPad->SaveAs("v2_vs_pT_Lc_1208.png");
}