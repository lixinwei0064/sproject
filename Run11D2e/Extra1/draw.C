void draw(){
    //D0
    TFile *D0file = new TFile("e_D0_10_40.root");
    TH2D *phiD0 = (TH2D *)D0file->Get("phiD0");

    TF1 *func = new TF1("func","[1]*(1+2*[0]*TMath::Cos(2*x))",0,6.283);
    TGraph *v2D0 = new TGraph(80); 
    v2D0->SetNameTitle("v2D0","D^{0}");
    TH1D *phicount;
    for(int i=1;i<=80;i++){
        phicount = phiD0->ProjectionX("",i,i);
        phicount->Fit(func);
        v2D0->SetPoint(i-1,0.1*i,func->GetParameter(0));
    }
    v2D0->SetFillColor(0);
    v2D0->SetLineWidth(2);
    
    TGraph *v2eD0 = (TGraph *)D0file->Get("v2D");//Was generated in former steps in the same manner as above
    v2eD0->SetTitle("e. from D^{0} decay");
    v2eD0->SetFillColor(0);
    v2eD0->SetLineWidth(2);

    //Lc
    TFile *Lcfile = new TFile("e_Lc.root");
    TH2D *phiLc = (TH2D *)Lcfile->Get("phiLc");

    TGraph *v2Lc = new TGraph(80); 
    v2Lc->SetNameTitle("v2Lc","#Lambda_{c}");
    for(int i=1;i<=80;i++){
        phicount = phiLc->ProjectionX("",i,i);
        phicount->Fit(func);
        v2Lc->SetPoint(i-1,0.1*i,func->GetParameter(0));
    }
    v2Lc->SetFillColor(0);
    v2Lc->SetLineWidth(2);
 
    TGraph *v2eLc = (TGraph *)Lcfile->Get("v2D");
    v2eLc->SetTitle("e. from #Lambda_{c} decay");
    v2eLc->SetFillColor(0);
    v2eLc->SetLineWidth(2);
    //Draw
    TMultiGraph *mg = new TMultiGraph();
    mg->SetTitle(";p_{T};v_{2}");
    mg->Add(v2D0);
    mg->Add(v2eD0);
    mg->Add(v2Lc);
    mg->Add(v2eLc);
    mg->Draw("AC PLC");
    gPad->BuildLegend(0.1,0.1,0.4,0.3);
    gPad->SaveAs("v2_10_40.png");

}