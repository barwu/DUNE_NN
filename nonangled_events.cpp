#include <iostream>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <vector>
using namespace std;
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveStats.h"

float max_value=0.;
vector<vector<vector<double>>>* xyz_mom=nullptr;
vector<vector<double>>* muon_tra=nullptr;
double numu_e=0;
double e_vis_true=0;
const int NUM_VTX=22, NUM_LAR_DTR=15;
double LAr_position[NUM_LAR_DTR]={-2800.,-2575.,-2400.,-2175.,-2000.,-1775.,-1600.,-1375.,-1200.,-975.,-800.,-575.,-400.,-175.,0.};
double vertex_position[NUM_VTX]={-299.,-292.,-285.,-278.,-271.,-264.,-216.,-168.,-120.,-72.,-24.,24.,72.,120.,168.,216.,264.,271.,
                      278.,285.,292.,299.};
double total_detected[5][NUM_LAR_DTR][NUM_VTX]={};
//float scale[30]={.19,.18,.18,.7,.7,1.05,.04,.034,.034,.07,.07,.07,.23,.21,.21,.73,.65,1.05,1.,1.,1.,1.05,1.05,1.05,.23,.2,.2,.7,.7,1.};
float scale[45]={.19,.18,.18,.7,.7,1.05,1.,1.05,1.,.04,.034,.034,.07,.07,.07,1.,.075,.075,.23,.21,.21,.73,.65,1.05,1.,1.05,1.05,1.,
                1.,1.,1.05,1.05,1.05,1.,1.05,1.05,.23,.2,.2,.7,.7,1.,1.,1.05,1.};

struct Para
{
  //static constexpr const char *const S;
  //constexpr const *char , VTX_X="vtx_x", *VTX_Y="vtx_y", *VTX_Z="vtx_z";
  //const char *LMX="LepMomX", *LMY="LepMomY", *LMZ="LepMomZ";
  char field[30];
  double l;
  double h;
  //vector<vector<vector<double>>>* field_value;
};
Para pr[]=
{
  {"LepMomZ", -0.3, 5.5},
  {"LepMomTot", 0., 5.5}
};
/*
struct sel_type
{
  const char* sel_name;
  const char* eff_name;
  vector<vector<double>>* eff_value=nullptr;
  sel_type() {}
  sel_type(const char* sn, const char* en)
  :sel_name(sn),eff_name(en) {}
};

vector<sel_type> br=
{
  sel_type("muon_contained", "muon_contained_eff"),
  sel_type("muon_tracker", "muon_tracker_eff"),
  sel_type("muon_selected", "muon_selected_eff"),
  sel_type("hadron_selected", "hadron_selected_eff"),
  sel_type("combined", "combined_eff")
};
*/
void populate_histograms(char* eff,char* caf,vector<TH1D*>& hists1,vector<TH1D*>& hists2,int j)
{
  TFile eff_file(eff);
  TFile caf_file(caf);
  TTree *event_data=(TTree*)eff_file.Get("event_data");
  TTree *caf_tree=(TTree*)caf_file.Get("effTreeND");
  gInterpreter->GenerateDictionary("vector<vector<double>>", "vector");
  gInterpreter->GenerateDictionary("vector<vector<vector<double>>>", "vector");
  event_data->SetBranchAddress("muon_tracker_eff", &muon_tra);
  caf_tree->SetBranchAddress("ND_OffAxis_Sim_mu_start_p_xyz_LAr", &xyz_mom);

  Long64_t nentries1=event_data->GetEntries();
  Long64_t nentries2=caf_tree->GetEntries();
  if (nentries1!=nentries2) {cout<<"The efficiency file"<<eff<<"has"<<nentries2
  <<" events, and the CAF file"<<caf<<"has"<<nentries1<<"events."<<endl;}
  for (int i=0;i<nentries2;i++) {
    event_data->GetEntry(i);
    caf_tree->GetEntry(i);
    unsigned long lar_pos=14;
    int k=0;
    for (Para item:pr) {
      double variable_type=0.0;
      TH1D* hist1=hists1[k];
      TH1D* hist2=hists2[k];
      for (unsigned long vtx_pos=11;vtx_pos>9;vtx_pos--) {
      double cross_sectional_angle=atan((*xyz_mom)[lar_pos][vtx_pos][1]/(*xyz_mom)[lar_pos][vtx_pos][2]);
      if (abs(cross_sectional_angle)>0.2) {continue;}
      cout<<cross_sectional_angle<<endl;
      if (k==0) {
      variable_type=(*xyz_mom)[lar_pos][vtx_pos][2];
      } else if (k==1) {
        variable_type=sqrt(pow((*xyz_mom)[lar_pos][vtx_pos][0],2)+pow((*xyz_mom)[lar_pos][vtx_pos][1],2)
        +pow((*xyz_mom)[lar_pos][vtx_pos][2],2));
      }/*
	  if ((k==1)&&(max_value<variable_type)) {
	    max_value=variable_type;
	    cout<<max_value<<", ";
	  }*/
      vector<vector<double>>* eff_value=muon_tra;
      vector<vector<double>>& eff_value2=*eff_value;
      double geo_eff=eff_value2[lar_pos][vtx_pos];
      if (geo_eff>1.) {
        cout<<"efficiency of event "<<i<<" at position "<<LAr_position[lar_pos]<<", "
        <<vertex_position[vtx_pos]<<" is "<<geo_eff<<endl;
      }
      hist1->Fill(variable_type);
      hist2->Fill(variable_type, geo_eff);
      }
      k++;
    }
  }
  eff_file.Close();
  caf_file.Close();
}

void nonangled_events()
{
  char eff[99];
  char caf[99];
  vector<TH1D*> histograms1;
  vector<TH1D*> histograms2;
  for (Para item:pr)
  {
    const char *fd=item.field;
    double lowerbound=item.l;
    double upperbound=item.h;
    histograms1.push_back(new TH1D(Form("%s_hist",fd), Form("raw %s", fd), 300, lowerbound, upperbound));
    histograms2.push_back(new TH1D(Form("%s_hist",fd), Form("selected %s", fd), 300, lowerbound, upperbound));
  }

  for (int j=0; j<10; j++)
  {
    memset(eff, 0, 99); // clear array each time
    memset(caf, 0, 99);
    sprintf(eff, "/storage/shared/barwu/10thTry/FDEff/FDGeoEff_62877585_99%d_Eff.root", j);
    sprintf(caf, "/storage/shared/fyguo/FDGeoEff_nnhome/FDGeoEff_62877585_99%d.root", j);
    if(access(eff, 0)==0)
    {
      populate_histograms(eff,caf,histograms1,histograms2,j);
    } else {
      cout<<"Warning: missing file:"<<eff<<endl;
      continue;
    }
  }

  //gStyle->SetOptStat(000000000);
  gStyle->SetOptStat(111111111);
  TCanvas *c=new TCanvas("c","histograms",2000,1000);
  c->Divide(2,2);
  int n=1;
  for (Para& item:pr)
  {
    const char *fd=item.field;
    double lowerbound=item.l;
    double upperbound=item.h;
    c->cd(2*n-1);
    TH1D *hist2=histograms2[n-1];
    hist2->SetLineColor(kTeal-3);
    hist2->Draw("hist");
    TH1D *hist1=histograms1[n-1];
    hist1->SetLineColor(kPink);
    hist1->Draw("samehist");

    float max1=1.16*hist1->GetMaximum();
    hist2->SetAxisRange(lowerbound,upperbound,"X");
    hist2->SetAxisRange(0.,max1,"Y");
    hist2->SetTitle(Form("%s Selection Cut", fd));
    hist2->GetXaxis()->SetTitle(fd);
    hist2->GetYaxis()->SetTitle("# of events");
    TLegend *leg=new TLegend(0.1,0.8,0.35,0.9);
    leg->SetHeader("comparison");
    leg->AddEntry(hist1, "raw distribution");
    leg->AddEntry(hist2, "selection-cut distribution");
    leg->Draw();

    c->cd(2*n);
    TH1D *rplot=(TH1D*)hist2->Clone();
    rplot->Divide(hist1);
    //rplot->SetAxisRange(0.,scale[(i_select-1)*44+n-1],"Y");
    rplot->SetAxisRange(0.,0.1,"Y");
    rplot->SetLineColor(kBlue);
    rplot->Draw("hist");
    n++;
  }
  c->Update();
  c->SaveAs("/home/barwu/repos/MuonEffNN/10thTry/FD_non-angled_events.png");
}