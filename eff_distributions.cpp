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
double vertex_position[22]={-299.,-292.,-285.,-278.,-271.,-264.,-216.,-168.,-120.,-72.,-24.,24.,72.,120.,168.,
                216.,264.,271.,278.,285.,292.,299.};

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
  // sel_type("muon_contained", "muon_contained_eff", false),
  // sel_type("muon_tracker", "muon_tracker_eff", false),
  // sel_type("muon_selected", "muon_sel_eff", true),
  sel_type("hadron_selected", "hadron_selected_eff", false),
  // sel_type("combined", "combined_eff", false)
};

void populate_histograms(char* eff,vector<vector<TH1D*>>& hists,int j)
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
    for (unsigned long lar_pos=0;lar_pos<hadr_eff->size();lar_pos++) {
      for (unsigned long vtx_pos=0;vtx_pos<hadr_eff->at(j).size();vtx_pos++) {
        //calculation for the muon-selected cut
        //muon_sel_eff[lar_pos][vtx_pos]=muon_cont_eff[lar_pos][vtx_pos]+muon_tra_eff[lar_pos][vtx_pos];
        //if (muon_sel_eff<0.0001) {cout<<"event "<<i<<" of file #"<<j<<" has small selected efficiency"<<endl;}

	      int n=0;
        for (auto sel:br) {
          TH1D* hist=hists[n][lar_pos];
          n++;
          vector<vector<double>> eff_value2=*(hadr_eff);
          double geo_eff=eff_value2[lar_pos][vtx_pos];
          hist->Fill(vertex_position[vtx_pos], geo_eff);
        }
      }
    }
  }
  eff_file.Close();
}

void eff_distributions()
{
  char eff[99];
  vector<vector<TH1D*>> histograms;
  for(auto sel:br)
  {
    const char* dt=sel.sel_name;
    vector<TH1D*> histset;
    histograms.push_back(histset);
    for (int i=0;i<15;i++)
    {
      histograms.back().push_back(new TH1D(Form("%s_hist_%d",dt,i), Form("%s geometric efficiency %d", dt, i), 100, -300., 450.));
    }
  }

  for (int j=0; j<10; j++)
  {
    memset(eff, 0, 99); // clear array each time
    sprintf(eff, "/storage/shared/barwu/10thTry/FDEff/FDGeoEff_62877585_99%d_Eff.root", j);
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
  c->Divide(5,3);
  for (int i=1;i<16;i++)
  {
  int n=0;
  c->cd(i);
    for(auto sel:br)
    {
      const char *dt=sel.sel_name;
      TH1D *hist=histograms[n][i-1];
      //hist->SetLineColor(kAzure);
      hist->Draw("histL");

      float max1=1.15*hist->GetMaximum();
      hist->SetAxisRange(0.,max1,"Y");
      hist->SetTitle(Form("Geometric Efficiencies", dt));
      hist->GetXaxis()->SetTitle("Geometric Efficiency");
      hist->GetYaxis()->SetTitle("# of events");

      TLegend *leg=new TLegend(0.1,0.8,0.35,0.9);
      leg->SetHeader("comparison"); 
      leg->AddEntry(hist, "geometric efficiency");
      leg->Draw();
      n++;
    }
  c->Update();
  }
  c->SaveAs("/home/barwu/repos/MuonEffNN/10thTry/FD_geo_effs.png");
}