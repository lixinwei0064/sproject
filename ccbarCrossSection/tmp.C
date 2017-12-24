//draw systematic error
const float sysw = 0.15;
for(int i=0; i<gRatioSys[icent]->GetN(); i++) {
    const float sysl = gRatioSys[icent]->GetY()[i] * 0.05;
    TLine *llw = new TLine(gRatioSys[icent]->GetX()[i]-sysw,gRatioSys[icent]->GetY()[i]-gRatioSys[icent]->GetEY()[i],gRatioSys[icent]->GetX()[i]+sysw,gRatioSys[icent]->GetY()[i]-gRatioSys[icent]->GetEY()[i]);
    llw->SetLineWidth(2);
    llw->SetLineColor(kBlack);
    llw->Draw("same");
    TLine *lhi = new TLine(gRatioSys[icent]->GetX()[i]-sysw,gRatioSys[icent]->GetY()[i]+gRatioSys[icent]->GetEY()[i],gRatioSys[icent]->GetX()[i]+sysw,gRatioSys[icent]->GetY()[i]+gRatioSys[icent]->GetEY()[i]);
    lhi->SetLineWidth(2);
    lhi->SetLineColor(kBlack);
    lhi->Draw("same");
    TLine *lv = new TLine(gRatioSys[icent]->GetX()[i]-sysw,gRatioSys[icent]->GetY()[i]-gRatioSys[icent]->GetEY()[i],gRatioSys[icent]->GetX()[i]-sysw,gRatioSys[icent]->GetY()[i]-gRatioSys[icent]->GetEY()[i]+sysl);
    lv->SetLineWidth(2);
    lv->SetLineColor(kBlack);
    lv->Draw("same");
    TLine *lv = new TLine(gRatioSys[icent]->GetX()[i]+sysw,gRatioSys[icent]->GetY()[i]-gRatioSys[icent]->GetEY()[i],gRatioSys[icent]->GetX()[i]+sysw,gRatioSys[icent]->GetY()[i]-gRatioSys[icent]->GetEY()[i]+sysl);
    lv->SetLineWidth(2);
    lv->SetLineColor(kBlack);
    lv->Draw("same");
    TLine *lv = new TLine(gRatioSys[icent]->GetX()[i]-sysw,gRatioSys[icent]->GetY()[i]+gRatioSys[icent]->GetEY()[i],gRatioSys[icent]->GetX()[i]-sysw,gRatioSys[icent]->GetY()[i]+gRatioSys[icent]->GetEY()[i]-sysl);
    lv->SetLineWidth(2);
    lv->SetLineColor(kBlack);
    lv->Draw("same");
    TLine *lv = new TLine(gRatioSys[icent]->GetX()[i]+sysw,gRatioSys[icent]->GetY()[i]+gRatioSys[icent]->GetEY()[i],gRatioSys[icent]->GetX()[i]+sysw,gRatioSys[icent]->GetY()[i]+gRatioSys[icent]->GetEY()[i]-sysl);
    lv->SetLineWidth(2);
    lv->SetLineColor(kBlack);
    lv->Draw("same");
}
