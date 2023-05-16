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
#include "TBranch.h"

struct Para
{
  //static constexpr const char *const S;
  //constexpr const *char , VTX_X="vtx_x", *VTX_Y="vtx_y", *VTX_Z="vtx_z";
  //const char *LMX="LepMomX", *LMY="LepMomY", *LMZ="LepMomZ";
  char field[30];
  bool iscaf;
  double l;
  double h;
  double* field_value;
};

struct Sel_type
{
  const char* sel_name;
  const char* eff_name;
  bool calced=false;
  int* sel_value;
  double* eff_value;
  Sel_type() {}
  Sel_type(const char* sn, const char* en, bool c, int* sv, double* ev)
  :sel_name(sn),eff_name(en),calced(c),sel_value(sv),eff_value(ev) {}
};

int muon_cont, muon_tra, muon_sel, hadr, comb;
double muon_cont_eff, muon_tra_eff, muon_sel_eff, hadr_eff, comb_eff;
double x_pos, y_pos, z_pos, XLepMom, YLepMom, ZLepMom;
double TotalMom, cos_angle, LongitudinalMom;
const char* list_of_directories[40]={"0mgsimple","0m","1.75m","2m","4m","5.75m","8m","9.75m","12m","13.75m","16m","17.75m","20m","21.75m","24m","25.75m","26.75m","28m",
"28.25m","28.5m","0mgsimpleRHC","0mRHC","1.75mRHC","2mRHC","4mRHC","5.75mRHC","8mRHC","9.75mRHC","12mRHC","13.75mRHC","16mRHC","17.75mRHC","20mRHC","21.75mRHC","24mRHC",
"25.75mRHC","26.75mRHC","28mRHC","28.25mRHC","28.5mRHC"};
const int NUM_VTX=22, NUM_LAR_DTR=15;
double LAr_position[NUM_LAR_DTR]={-2800.,-2575.,-2400.,-2175.,-2000.,-1775.,-1600.,-1375.,-1200.,-975.,-800.,-575.,-400.,-175.,0.};
double vertex_position[NUM_VTX]={-299.,-292.,-285.,-278.,-271.,-264.,-216.,-168.,-120.,-72.,-24.,24.,72.,120.,168.,216.,264.,271.,278.,285.,292.,299.};

Para pr[]= //position is in units of cm, momentum is in units of GeV/c, angle is in units of rad, and energy is in  units of GeV
{
  {"vtx_x", true, -300., 300. , &x_pos},
  {"vtx_y", true, -100., 100., &y_pos},
  {"vtx_z", true, 50., 350., &z_pos},
  {"LepMomX", true, -2., 2., &XLepMom},
  {"LepMomY", true, -4., 2., &YLepMom},
  {"LepMomZ", true, -0.5, 4.5, &ZLepMom},
  {"TotMom", false, 0., 5.,&TotalMom},
  {"cos_LepNuAngle", false, 0., 1.,&cos_angle},
  {"LongMom", false, -1., 5.,&LongitudinalMom}
};

vector<Sel_type> br=
{
  Sel_type("muon_contained", "muon_contained_eff", false, &muon_cont, &muon_cont_eff),
  Sel_type("muon_tracker", "muon_tracker_eff", false, &muon_tra, &muon_tra_eff),
  Sel_type("muon_selected", "muon_sel_eff", true, &muon_sel, &muon_sel_eff),
  Sel_type("hadron_selected", "hadron_selected_eff", false, &hadr, &hadr_eff ),
  Sel_type("combined", "combined_eff", false, &comb, &comb_eff)
};

void ttree_histograms()
{
  gROOT->SetBatch(kTRUE);
  char eff[99];
  char caf[99];
  vector<TH1D*> histograms1, histograms2, histograms3;
  for(auto sel:br)
  {
    const char* dt=sel.sel_name;
    for(auto item:pr)
    {
      char *fd=item.field;
      double l=item.l;
      double h=item.h;
      histograms1.push_back(new TH1D(Form("h1_%s_%s", fd, dt), Form("raw %s %s",fd, dt), 100, l, h));
      histograms2.push_back(new TH1D(Form("h2_%s_%s", fd, dt), Form("selected %s %s",fd, dt), 100, l, h));
      histograms3.push_back(new TH1D(Form("h3_%s_%s", fd, dt), Form("geo corrected %s %s",fd, dt), 100, l, h));
    }
  }

  const float directory_number=0; // Sometimes the directory # is an integer, sometimes its a fraction. Remember to change the wildcard and variable type accordingly.
  cout<<directory_number<<endl;
  for (int j=0; j<30; j++)
  {
    memset(eff, 0, 99); // clear array each time
    memset(caf, 0, 99); 
    // sprintf(eff,"/storage/shared/barwu/10thTry/combined1/%02dm/%02d/FHC.%03d%04d.CAF_Eff.root",directory_number,j/1000,int((directory_number+1)*100+j/10000),j%10000);
    // sprintf(caf,"/storage/shared/barwu/10thTry/NDCAF/%02dm/%02d/FHC.%03d%04d.CAF.root",directory_number,j/1000,int((directory_number+1)*100+j/10000),j%10000);
    sprintf(eff,"/storage/shared/barwu/10thTry/combined1/0m/%02d/FHC.10%05d.CAF_Eff.root",j/1000,j);
    sprintf(caf,"/storage/shared/wshi/CAFs/NDFHC_PRISM/%02d/FHC.10%05d.CAF.root",j/1000,j);
    if(access(eff, 0)!=0) {continue;}
    TFile eff_file(eff);
    TFile caf_file(caf);
    TTree *event_data=(TTree*)eff_file.Get("event_data");
    TTree *caf_tree=(TTree*)caf_file.Get("caf");

    int isCC=0, inFV=0;
    for(auto sel:br) 
    {
      if(sel.calced) continue;
      event_data->SetBranchAddress(sel.sel_name, sel.sel_value);
      event_data->SetBranchAddress(sel.eff_name, sel.eff_value);
    }
    caf_tree->SetBranchAddress("isCC",&isCC);
    event_data->SetBranchAddress("inFV",&inFV);   
    for (auto item:pr) 
    {
      TTree *tree=item.iscaf?caf_tree:event_data;
      tree->SetBranchAddress(item.field, item.field_value);
    }

    //there are some non-CC events
    Long64_t nentries1=caf_tree->GetEntries();
    Long64_t nentries2=event_data->GetEntries();
    if (nentries1!=nentries2) {cout<<"The efficiency file"<<eff<<"has"<<nentries2<<" events, and the CAF file"<<caf<<"has"<<nentries1<<"events."<<endl;}
    for (int i=0;i<nentries2;i++) {
      caf_tree->GetEntry(i);
      event_data->GetEntry(i);
      int sel_cut=isCC*inFV;
      if (sel_cut==0) {continue;}
      if (sel_cut!=1) {
        cout<<"sel_cut="<<sel_cut<<endl;
        continue;
      }
      //calculation for the muon-selected cut
      muon_sel=muon_cont+muon_tra;
      muon_sel_eff=muon_cont_eff+muon_tra_eff;
      if (muon_sel!=0&&muon_sel!=1) {
        cout<<"bad val for muon-selected check! "<<muon_sel<<endl;
        cout<<"Event #"<<i<<" in file "<<eff<<endl;
        cout<<"contained-"<<muon_cont<<", tracker-matched-"<<muon_tra<<endl;
        continue;
      }
     
      int n=0;
      for (auto sel:br) {
        for (auto item:pr) {
          const char *fd=item.field;
          TH1D* hist1=histograms1.at(n);
          TH1D* hist2=histograms2.at(n);
          TH1D* hist3=histograms3.at(n);
          n++;
          hist1->Fill(*item.field_value);
          hist2->Fill(*item.field_value,*sel.sel_value);
          double geo_eff=*sel.eff_value;
          if (geo_eff<=0.0001) {
            continue;
          } else {
            hist3->Fill(*item.field_value,*sel.sel_value/geo_eff);
          }
        }
      }
    }
    eff_file.Close();
    caf_file.Close();
  }

  TFile *hist_file=new TFile("/storage/shared/barwu/10thTry/xxxx.root","recreate");
  TTree *raw_hists=new TTree("raw data histograms","raw_hists");
  TTree *sel_hists=new TTree("selection-cut histograms","sel_hists");
  TTree *geo_hists=new TTree("geometrically-corrected histograms","geo_hists");
  int index=0;
  for (Para item:pr)
  {
    const char *fd=item.field;
    TH1D* hist1=histograms1.at(index);
    TBranch *raw=raw_hists->Branch(Form("raw_%s",fd),"TH1D",&hist1);
    index++;
  }
  index=0;
  for (auto sel:br)
  {
    const char *dt=sel.sel_name;
    for (Para item:pr)
    {
      const char *fd=item.field;
      TH1D* hist2=histograms2.at(index);
      TH1D* hist3=histograms3.at(index);
      TBranch *selected=sel_hists->Branch(Form("selected_%s_%s",dt,fd),"TH1D",&hist2);
      TBranch *geo_corrected=geo_hists->Branch(Form("geo-corrected_%s_%s",dt,fd),"TH1D",&hist3);
      index++;
    }
  }
  hist_file->Write();
  delete hist_file;
}