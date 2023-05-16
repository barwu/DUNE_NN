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

vector<vector<vector<double>>>* xyz_pos=nullptr;
vector<vector<vector<double>>>* xyz_mom=nullptr;
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
  {"vtx_x", -300., 300.},
  {"vtx_y", 5.385, 5.39},
  {"vtx_z",  659.995, 660.005},
  {"LepMomX", -2.1, 3.},
  {"LepMomY", -8., 2.},
  {"LepMomZ", -1., 14.},
  {"LepMomTot", 0., 16.},
  //{"LepNuAngle", 0., 1.}
  //{"LongMom", 0., 16.}
  {"ND_Gen_numu_E", 0., 10.},
  {"ND_E_vis_true", 0., 10.}
};

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

void populate_histograms(char* eff,char* caf,vector<vector<TH1D*>>& hists1,vector<vector<TH1D*>>& hists2,int j)
{
  TFile eff_file(eff);
  TFile caf_file(caf);
  TTree *event_data=(TTree*)eff_file.Get("event_data");
  TTree *thing=(TTree*)caf_file.Get("effTreeND");
  //gSystem->Exec("rm -f AutoDict*vector*vector*vector*double*");
  gInterpreter->GenerateDictionary("vector<vector<vector<double>>>", "vector");
  for(auto& sel:br) //efficiencies are in 3d array, but energy is in 1d array
  {
    event_data->SetBranchAddress(sel.eff_name, &(sel.eff_value));
  }
  thing->SetBranchAddress("ND_OffAxis_Sim_mu_start_v_xyz_LAr", &xyz_pos);
  thing->SetBranchAddress("ND_OffAxis_Sim_mu_start_p_xyz_LAr", &xyz_mom);
  thing->SetBranchAddress("ND_Gen_numu_E", &numu_e);
  thing->SetBranchAddress("ND_E_vis_true", &e_vis_true);

  Long64_t nentries1=event_data->GetEntries();
  Long64_t nentries2=thing->GetEntries();
  if (nentries1!=nentries2) {cout<<"The efficiency file"<<eff<<"has"<<nentries2
  <<" events, and the CAF file"<<caf<<"has"<<nentries1<<"events."<<endl;}
  for (int i=0;i<nentries2;i++) {
    event_data->GetEntry(i);
    thing->GetEntry(i);
    unsigned long lar_pos=14;
    int k=0;
    for (Para item:pr) {
      double var_type=0.0;
      if (k==7) {var_type=numu_e;}
      if (k==8) {var_type=e_vis_true;}

      for (unsigned long vtx_pos=0;vtx_pos<NUM_VTX;vtx_pos++) {
	      int n=0;
        for (auto& sel:br) {
          TH1D* hist1=hists1[n][k];
          TH1D* hist2=hists2[n][k];
          n++;
          if (k<3) {
            var_type=(*xyz_pos)[lar_pos][vtx_pos][k];
          } else if (k<6) {
            var_type=(*xyz_mom)[lar_pos][vtx_pos][k-3];
          } else if (k==6) {
            var_type=sqrt(pow((*xyz_mom)[lar_pos][vtx_pos][0],2)+pow((*xyz_mom)[lar_pos][vtx_pos][1],2)
              +pow((*xyz_mom)[lar_pos][vtx_pos][2],2));
          }
          //if (k==1) {cout<<var_type<<" ";}
          vector<vector<double>>* eff_value=sel.eff_value;
          vector<vector<double>>& eff_value2=*eff_value;
          double geo_eff=eff_value2[lar_pos][vtx_pos];
          if (geo_eff>1.) {
            cout<<"efficiency of event "<<i<<" at position "<<LAr_position[lar_pos]<<", "
            <<vertex_position[vtx_pos]<<" is "<<geo_eff<<endl;
          }
          hist1->Fill(var_type);
          hist2->Fill(var_type, geo_eff);
        }
      }
      k++;
    }
  }
  eff_file.Close();
  caf_file.Close();
}

void FD_selection_cuts()
{
  char eff[99];
  char caf[99];
  vector<vector<TH1D*>> histograms1;
  vector<vector<TH1D*>> histograms2;
  for(auto& sel:br)
  {
    const char* dt=sel.sel_name;
    vector<TH1D*> histset1, histset2;
    histograms1.push_back(histset1);
    histograms2.push_back(histset2);
    int i=0;
    for (Para item:pr)
    {
      double lowerbound=item.l;
      double upperbound=item.h;
      histograms1.back().push_back(new TH1D(Form("%s_hist_%d",dt,i), Form("raw %s  %d", dt, i), 300, lowerbound, upperbound));
      histograms2.back().push_back(new TH1D(Form("%s_hist_%d",dt,i), Form("selected %s %d", dt, i), 300, lowerbound, upperbound));
    i++;
    }
  }

  for (int j=1; j<=10; j++)
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
  TCanvas *c=new TCanvas("c","Energy Distributions",2000,1000);
  TCanvas *r=new TCanvas("r","Ratio Plots",2000,1000);
  c->Divide(9,5);
  r->Divide(9,5);
  int i=1;
  int i_select=0;
  for (auto& sel:br)
  {
    const char *dt=sel.sel_name;
    int n=0;
    for(Para item:pr)
    {
      const char *fd=item.field;
      double lowerbound=item.l;
      double upperbound=item.h;
      TVirtualPad *pad=c->cd(i);
      if (n%9==7) {pad->SetLogy();} //pad needs to be made logarithmic, not canvas
      TH1D *hist2=histograms2[i_select][n];
      hist2->SetLineColor(kTeal-3);
      hist2->Draw("hist");
      TH1D *hist1=histograms1[i_select][n];
      hist1->SetLineColor(kPink);
      hist1->Draw("samehist");

      float max1=hist1->GetMaximum();
      //hist2->SetAxisRange(lowerbound,upperbound,"X");
      if (n!=7) hist2->SetAxisRange(0.,1.16*max1,"Y");
      hist2->SetTitle(Form("%s %s Selection Cut", fd, dt));
      hist2->GetXaxis()->SetTitle(fd);
      hist2->GetYaxis()->SetTitle("# of events");
      TLegend *leg=new TLegend(0.1,0.8,0.35,0.9);
      leg->SetHeader("comparison");
      leg->AddEntry(hist1, "raw distribution");
      leg->AddEntry(hist2, "selection-cut distribution");
      leg->Draw();

      r->cd(i);
      TH1D *rplot=(TH1D*)hist2->Clone();
      rplot->Divide(hist1);
      rplot->SetAxisRange(0.,scale[i-1],"Y");
      rplot->SetLineColor(kBlue);
      rplot->Draw("hist");
      n++;
      i++;
    }
    c->Update();
    r->Update();
    i_select++;
  }
  c->SaveAs("/home/barwu/repos/MuonEffNN/10thTry/FD_hists_all.png");
  r->SaveAs("/home/barwu/repos/MuonEffNN/10thTry/FD_hists_ratios.png");

  TCanvas *cs[5];
  TCanvas *rs[5];
  i=1;
  for (auto& sel:br)
  {
    const char *dt=sel.sel_name;
    cs[i-1]=new TCanvas(Form("c%01d",i),dt,2000,1000);
    cs[i-1]->Divide(3,3);
    rs[i-1]=new TCanvas(Form("r%01d",i),dt,2000,1000);
    rs[i-1]->Divide(3,3);
    int n=0;
    for(Para& item:pr)
    {
      const char *fd=item.field;
      double lowerbound=item.l;
      double upperbound=item.h;
      TVirtualPad *pads=cs[i-1]->cd(n+1);
      if (n%9==7) {pads->SetLogy();}
      TH1D *hist2=histograms2[i-1][n];
      hist2->SetLineColor(kTeal-3);
      hist2->Draw("hist");
      TH1D *hist1=histograms1[i-1][n];
      hist1->SetLineColor(kPink);
      hist1->Draw("samehist");

      float max1=1.16*hist1->GetMaximum();
      if (n!=7) hist2->SetAxisRange(lowerbound,upperbound,"X");
      hist2->SetAxisRange(0.,max1,"Y");
      hist2->SetTitle(Form("%s %s Selection Cut", fd, dt));
      hist2->GetXaxis()->SetTitle(fd);
      hist2->GetYaxis()->SetTitle("# of events");
      TLegend *leg=new TLegend(0.1,0.8,0.35,0.9);
      leg->SetHeader("comparison");
      leg->AddEntry(hist1, "raw distribution");
      leg->AddEntry(hist2, "selection-cut distribution");
      leg->Draw();

      rs[i-1]->cd(n+1);
      TH1D *rplot=(TH1D*)hist2->Clone();
      rplot->Divide(hist1);
      rplot->SetAxisRange(0.,scale[(i-1)*9+n],"Y");
      rplot->SetLineColor(kBlue);
      rplot->Draw("hist");
      n++;
    }
    cs[i-1]->Update();
    rs[i-1]->Update();
    cs[i-1]->SaveAs(Form("/home/barwu/repos/MuonEffNN/10thTry/FD_%s_hists.png", dt));
    rs[i-1]->SaveAs(Form("/home/barwu/repos/MuonEffNN/10thTry/FD_%s_hists_ratios.png", dt));
    i++;
  }
}