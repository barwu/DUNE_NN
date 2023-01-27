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
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TAxis.h"
#include <vector>

const int NUM_VTX=22, NUM_LAR_DTR=15;
double LAr_position[NUM_LAR_DTR]={-2800.,-2575.,-2400.,-2175.,-2000.,-1775.,-1600.,-1375.,-1200.,-975.,-800.,-575.,-400.,-175.,0.};
double vertex_position[NUM_VTX]={-299.,-292.,-285.,-278.,-271.,-264.,-216.,-168.,-120.,-72.,-24.,24.,72.,120.,168.,216.,264.,271.,
                      278.,285.,292.,299.};
double total_detected[5][NUM_LAR_DTR][NUM_VTX]={};
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
  sel_type("muon_contained", "muon_contained_eff", false),
  sel_type("muon_tracker", "muon_tracker_eff", false),
  //sel_type("muon_selected", "muon_sel_eff", true),
  sel_type("hadron_selected", "hadron_selected_eff", false),
  sel_type("combined", "combined_eff", false)
};

void accumulate_effs(char* eff,int j)
{
  TFile eff_file(eff);
  TTree *event_data=(TTree*)eff_file.Get("event_data");
  for(auto& sel:br) //efficiencies are in 3d array, but energy is in 1d array
  {
    if(sel.calced) continue;
    event_data->SetBranchAddress(sel.eff_name, &(sel.eff_value));
  }

  Long64_t nentries=event_data->GetEntries();
  for (int i_event=0;i_event<1;i_event++) {
    event_data->GetEntry(i_event);
	  int i_select=0;
    for (sel_type& sel:br) {
      vector<vector<double>>* eff_value=sel.eff_value;
      vector<vector<double>>& eff_value2= *eff_value;
      //vector<vector<double>>* eff_value2 = hadr_eff_ptr;
      for (unsigned long lar_pos=0;lar_pos<NUM_LAR_DTR;lar_pos++) {
        for (unsigned long vtx_pos=0;vtx_pos<NUM_VTX;vtx_pos++) {
        //calculation for the muon-selected cut
        //muon_sel_eff[lar_pos][vtx_pos]=muon_cont_eff[lar_pos][vtx_pos]+muon_tra_eff[lar_pos][vtx_pos];
        //if (muon_sel_eff[lar_pos][vtx_pos]<0.0001) {
          //cout<<"event "<<i<<" of file #"<<j<<" has small selected efficiency"<<endl;}
          // TGraph* hist=hists[n][lar_pos];
          double geo_eff=eff_value2[lar_pos][vtx_pos];
          // hist->Fill(vertex_position[vtx_pos].q, geo_eff);
          total_detected[i_select][lar_pos][vtx_pos]+=geo_eff;
        }
      }
      i_select++;
    }
  }
  eff_file.Close();
}

void eff_distributions()
{
  char eff[99];
  std::memset(total_detected, 0, 5 * NUM_LAR_DTR * NUM_VTX * sizeof(double));
  for (int j=1; j<9; j++)
  {
    memset(eff, 0, 99); // clear array each time
    sprintf(eff, "/storage/shared/barwu/10thTry/FDEff/FDGeoEff_62877585_99%d_Eff.root", j);
    if(access(eff, 0)==0)
    {
      accumulate_effs(eff,j);
    } else {
      cout<<"Warning: missing file:"<<eff<<endl;
      continue;
    }
  }
  for (int j=0;j<NUM_LAR_DTR;j++) {
    for (int k=0;k<NUM_VTX;k++) {
      for (int i=0;i<4;i++) {
        cout<<"LAr pos="<<LAr_position[j]<<", vtx pos="<<vertex_position[k]<<", eff="<<total_detected[i][j][k]<<endl;
      }
    }
  }

  vector<vector<TGraph*>> graphs;
  int i_select=0;
  for(auto& sel:br)
  {
    vector<TGraph*> vgraphs;
    graphs.push_back(vgraphs);
    for (int i=0;i<NUM_LAR_DTR;i++)
    {
      graphs.back().push_back(new TGraph(NUM_VTX,vertex_position,total_detected[i_select][i]));
    }
    i_select++;
  }
  for (int j=0;j<NUM_LAR_DTR;j++) {
    for (int k=0;k<NUM_VTX;k++) {
      for (int i=0;i<4;i++) {
        cout<<"LAr pos="<<LAr_position[j]<<", vtx pos="<<vertex_position[k]<<", eff="<<total_detected[i][j][k]<<endl;
      }
    }
  }

  gStyle->SetOptStat(111111111);
  TCanvas *c=new TCanvas("c","Energy Distributions",1800,1000);
  TMultiGraph *mgs[NUM_LAR_DTR];
  c->Divide(5,3);
  for (int i=0;i<NUM_LAR_DTR;i++) //pad
  {
    TMultiGraph *mg=new TMultiGraph();
    mgs[i]=mg;
    TGraph *linegraphs[5];
    TVirtualPad *pad=c->cd(i+1);
    int i_select=0;
    for (auto sel:br)
    {
      const char *selection=sel.sel_name;
      linegraphs[i_select]=graphs[i_select][i];
      linegraphs[i_select]->SetMarkerStyle(2*i_select+58);
      linegraphs[i_select]->SetMarkerColor(i_select+1);
      linegraphs[i_select]->SetLineColor(i_select+1);
      linegraphs[i_select]->SetTitle(selection);
      mg->Add(linegraphs[i_select]);
      i_select++;
    }
    int det_pos=LAr_position[i];
    mg->SetTitle(Form("LAr pos=%d", det_pos));
    mg->GetXaxis()->SetTitle("Vertex Position (cm)");
    mg->GetYaxis()->SetTitle("# of detected events");
    mg->Draw("ALP");
    pad->BuildLegend();
  }
  c->Update();
}