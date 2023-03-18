#include <iostream>
#include <vector>
#include <filesystem>
using namespace std;
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveStats.h"

//float max_value=0.;
int inFV=0, isCC=0;
double x_pos;
double y_mom;
double z_mom;
const int NUM_VTX=22, NUM_LAR_DTR=15;
double LAr_position[NUM_LAR_DTR]={-2800.,-2575.,-2400.,-2175.,-2000.,-1775.,-1600.,-1375.,-1200.,-975.,-800.,-575.,
                   -400.,-175.,0.};
double vertex_position[NUM_VTX]={-299.,-292.,-285.,-278.,-271.,-264.,-216.,-168.,-120.,-72.,-24.,24.,72.,120.,168.,
                      216.,264.,271.,278.,285.,292.,299.};
double separate_positions[12]={300.,296.,288.5,281.5,274.5,267.5,240.,168.,144.,96.,48.,0.};

void populate_histograms(char* eff,char* caf,TH1D* hists[22],int j)
{
  TFile eff_file(eff);
  TFile caf_file(caf);
  TTree *eff_tree=(TTree*)eff_file.Get("event_data");
  TTree *caf_tree=(TTree*)caf_file.Get("caf");
  eff_tree->SetBranchAddress("isCC", &isCC);
  eff_tree->SetBranchAddress("inFV", &inFV);
  caf_tree->SetBranchAddress("vtx_x", &x_pos);
  caf_tree->SetBranchAddress("LepMomY", &y_mom);
  caf_tree->SetBranchAddress("LepMomZ", &z_mom);
  Long64_t nentries1=eff_tree->GetEntries();
  Long64_t nentries2=caf_tree->GetEntries();
  if (nentries1!=nentries2) {cout<<"The efficiency file"<<eff<<"has"<<nentries1
  <<" events, and the CAF file"<<caf<<"has"<<nentries2<<"events."<<endl;}
  for (int i=0;i<nentries2;i++) {
    eff_tree->GetEntry(i);
    caf_tree->GetEntry(i); 
    int sel_cut=isCC*inFV;
    if (sel_cut==0) {continue;}
    if (x_pos<-300.||x_pos>300.) {
      cout<<x_pos<<endl;
      continue;
    }
    int vtx_pos;
    if (x_pos<0.) {
      vtx_pos=0;
      while ((-x_pos)<separate_positions[vtx_pos+1]) {vtx_pos++;}
    } else {
      vtx_pos=21;
      while (x_pos<separate_positions[22-vtx_pos]) {vtx_pos--;}
    }
    TH1D* hist=hists[vtx_pos];
    double cross_sectional_angle=atan(y_mom/z_mom);
    hist->Fill(cross_sectional_angle);
  }
  eff_file.Close();
  caf_file.Close();
}

void ND_muon_angles()
{
  char eff[99];
  char caf[99];
  TH1D* pos_selected_hists[22];
  for (int vtx_pos=0;vtx_pos<NUM_VTX;vtx_pos++)
  {
    if (vtx_pos<11)
    {
      pos_selected_hists[vtx_pos]=new TH1D(Form("c%d", vtx_pos+1),Form("cross-sectional angle between %f and %f cm",
      separate_positions[vtx_pos],separate_positions[vtx_pos+1]),200,-3.14159265358979323846,3.14159265358979323846);
    } else {
      pos_selected_hists[vtx_pos]=new TH1D(Form("c%d", vtx_pos+1),Form("cross-sectional angle between %f and %f cm",
      separate_positions[vtx_pos-11],separate_positions[vtx_pos-10]),200,-3.14159265358979323846,3.14159265358979323846);
    }
  }

  for (int j=0; j<30000; j++)
  {
    memset(eff, 0, 99); // clear array each time
    memset(caf, 0, 99);
    sprintf(eff,"/storage/shared/barwu/10thTry/combined1/0m/FHC.10%05d.CAF_Eff.root",j);
    sprintf(caf,"/storage/shared/wshi/CAFs/NDFHC_PRISM/%02d/FHC.10%05d.CAF.root",j/1000,j);
    //sprintf(eff,"/storage/shared/barwu/9thTry/eff_trees/FHC.100%04d.CAF_MuonEff.root",j);
    //sprintf(caf,"/storage/shared/cvilela/CAF/ND_v7/0%01d/FHC.100%04d.CAF.root",j/1000,j);
    if(access(eff, 0)==0) {
    populate_histograms(eff,caf,pos_selected_hists,j);
    } else {
      //cout<<"Warning: missing file:"<<eff<<endl;
      continue;
    }
  }

  gStyle->SetOptStat(111111111);
  TCanvas *c=new TCanvas("Angle of FD muon momentum","c",2000,1000);;
  c->Divide(6,4);
  for (int vtx_pos=0;vtx_pos<NUM_VTX;vtx_pos++)
  {
    c->cd(vtx_pos+1);
    TH1D *hist=pos_selected_hists[vtx_pos];
    hist->SetLineColor(kBlue);
    hist->Draw("hist");
    float max1=hist->GetMaximum();
    hist->SetAxisRange(0.,1.16*max1,"Y");
    //hist->SetTitle(Form("%s %s Selection Cut", fd, dt));
    hist->GetXaxis()->SetTitle("Muon angle (radians)");
    hist->GetYaxis()->SetTitle("# of events");
  }
  c->Update();
  c->SaveAs("/home/barwu/repos/MuonEffNN/10thTry/ND_on-axis_muon_angle_hists_radians.png");
}