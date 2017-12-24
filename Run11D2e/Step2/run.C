R__LOAD_LIBRARY(D0.C)
R__LOAD_LIBRARY(Lc.C)
void run(const int cen = 0, const int cent = 0,const int seed = 32553, const int seed2 = 32553)
{
  const char nameCent[4][250] = {"model1_diquark", "model2_threequark", "model3_Greco", "mean"};
  const char name[2][250] = {"0_80","10_40"};
  gSystem->Load("Lc_C.so");
  gSystem->Load("D0_C.so");
  
  TCanvas* c1 = new TCanvas("c1", "A Canvas",10,10,1000,1000);
  auto legend = new TLegend(0.65,0.7,0.9,0.9);
  TRandom *myRandom = new TRandom();

  D0(10000000, cen, seed, legend);
  Lc(10000000, cent, seed2, legend);
  c1->SetLogy();  

  char outname[100];
  sprintf(outname,"Decay_%s_%s.png",name[cen],nameCent[cent]);
  c1->SaveAs(outname);
}

