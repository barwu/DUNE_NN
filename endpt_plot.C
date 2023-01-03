//Created on Tue Jul  6 01:01:41 2021

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
using namespace std;

int endpt_plot();

int endpt_plot()
{
TChain chain("muonEfficiency");
chain.Add("/home/barwu/repos/MuonEffNN/endpt_test/0.0005/FHC.100*.CAF_MuonEff.root");
TTree *muonEfficiency=&chain;
TCanvas* c1=new TCanvas("c1", "muon ending positions", 1600, 900);
c1->Divide(2,2);
c1->cd(1);
muonEfficiency->GetEntries();
muonEfficiency->Draw("muon_endpoint_x:muon_endpoint_y:muon_endpoint_z:2.5e>>hist", "");
c1->cd(2);
muonEfficiency->Draw("muon_endpoint_x:muon_endpoint_y>>hist_xy(100,-100,100,100,-100,100)");
TH2D* histogram_xy=(TH2D*) gDirectory->Get("hist_xy");
histogram_xy->Draw("COLHIST");
c1->cd(3);
muonEfficiency->Draw("muon_endpoint_y:muon_endpoint_z>>hist_yz(100,-100,100,100,-100,100)");
TH2D* histogram_yz=(TH2D*) gDirectory->Get("hist_yz");
histogram_yz->Draw("COLHIST");
c1->cd(4);
muonEfficiency->Draw("muon_endpoint_x:muon_endpoint_z>>hist_xz(100,-100,100,100,-100,100)");
TH2D* histogram_xz=(TH2D*) gDirectory->Get("hist_xz");
//histogram->Fill(muon_endpoint_x,muon_endpoint_z);
//muonEfficiency->GetYaxis()->SetRange(-9,9);
//muonEfficiency->GetXaxis()->SetRange(-9,9);
histogram_xz->Draw("COLHIST");
c1->Update();
return 0;
}