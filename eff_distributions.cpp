#include <iostream>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <algorithm>
#include <cmath>
#include <filesystem>
using namespace std;
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include <vector>

vector<vector<double>>* muon_cont_eff=nullptr;
vector<vector<double>>* muon_tra_eff=nullptr;
vector<vector<double>>* muon_sel_eff=nullptr;
vector<vector<double>>* hadr_eff=nullptr;
vector<vector<double>>* comb_eff=nullptr;

struct sel_type
{
  const char* sel_name;
  const char* eff_name;
  bool calced=false;
  vector<vector<double>>* eff_value=nullptr;
  sel_type() {}
  sel_type(const char* sn, const char* en, bool c)
  :sel_name(sn),eff_name(en),calced(c) {}
};

vector<sel_type> br=
{
  //sel_type("muon_contained", "muon_contained_eff", false, muon_cont_eff),
  //sel_type("muon_tracker", "muon_tracker_eff", false, muon_tra_eff),
  //sel_type("muon_selected", "muon_sel_eff", true, muon_sel_eff),
  sel_type("hadron_selected", "hadron_selected_eff", false),
  //sel_type("combined", "combined_eff", false, &comb_eff)
};

void populate_histograms(char* eff,TH1D* hists[9],int j)
{
  TFile eff_file(eff);
  TTree *event_data=(TTree*)eff_file.Get("event_data");
  for(auto sel:br) //efficiencies are in 3d array, but energy is in 1d array
  {
    if(sel.calced) continue;
    event_data->SetBranchAddress(sel.eff_name, &(sel.eff_value));
    hadr_eff=sel.eff_value;
  }

  Long64_t nentries=event_data->GetEntries();
  for (int i=0;i<nentries;i++) {
    event_data->GetEntry(i);
    for (int lar_pos=0;lar_pos<hadr_eff->size();lar_pos++) {
      for (int vtx_pos=0;vtx_pos<hadr_eff->at(j).size();vtx_pos++) {
        //calculation for the muon-selected cut
        //muon_sel_eff[lar_pos][vtx_pos]=muon_cont_eff[lar_pos][vtx_pos]+muon_tra_eff[lar_pos][vtx_pos];
        //if (muon_sel_eff<0.0001) {cout<<"event "<<i<<" of file #"<<j<<" has small selected efficiency"<<endl;}

	    int n=0;
        for (auto sel:br) {
          TH1D* hist=hists[n];
          n++;
          vector<vector<double>> eff_value2=*(hadr_eff);
          double geo_eff=eff_value2[lar_pos][vtx_pos];
          hist->Fill(geo_eff);
        }
      }
    }
  }
  eff_file.Close();
}

void eff_distributions()
{
  char eff[99];
  char caf[99];
  TH1D *histograms[9];
  int m=0;
  for(auto sel:br)
  {
    const char* dt=sel.sel_name;
    histograms[m]=new TH1D(Form("%s_hist", dt), Form("%s geometric efficiency", dt), 100, 0, 1.2);
    m++;
  }

  for (int j=0; j<=9; j++)
  {
    memset(eff, 0, 99); // clear array each time
    sprintf(eff,"/storage/shared/barwu/10thTry/FDEff/FDGeoEff_62877585_99%d_Eff.root",j);
    if(access(eff, 0)==0)
    {
      populate_histograms(eff,histograms,j);
    } else {
      cout<<"Warning: missing file:"<<eff<<endl;
      continue;
    }
  }

  //gStyle->SetOptStat(000000000);
  gStyle->SetOptStat(111111111);
  TCanvas *c=new TCanvas("c","Energy Distributions",1800,1000);
  //c->Divide(3,2);
  int n=0;
  for(auto sel:br)
  {
    const char *dt=sel.sel_name;
    //c->cd(i);
    TH1D *hist=histograms[n];
    hist->SetLineColor(kAzure);
    hist->Draw("histS");

    float max1=1.01*hist->GetMaximum();
    hist->SetAxisRange(0.,max1,"Y");
    hist->SetTitle(dt);
    hist->GetXaxis()->SetTitle("Geometric Efficiency");
    hist->GetYaxis()->SetTitle("# of events");

    TLegend *leg=new TLegend(0.1,0.8,0.35,0.9);
    leg->SetHeader("comparison"); 
    leg->AddEntry(hist, "geometric efficiency");
    leg->Draw();
    n++;
  }
  c->Update();
  //c->SaveAs("/home/barwu/repos/MuonEffNN/test_eff_hist.png");
}