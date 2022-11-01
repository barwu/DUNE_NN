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
#include "TH2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"

void populate_histograms(char* eff,char* CAF,TH2D* heat_map,TH1D* diff,int file_num)
{
  TFile eff_file(eff);
  TFile caf_file(CAF);
  TTree *event_data=(TTree*)eff_file.Get("event_data");
  TTree *caf=(TTree*)caf_file.Get("caf");
  int isCC=0, inFV=0;
  double Elep_reco=0., eRecoP=0., eRecoN=0., eRecoPip=0., eRecoPim=0., eRecoPi0=0., eRecoOther=0.;
  double LepE=0., eP=0., eN=0., ePip=0., ePim=0., ePi0=0., eOther=0.;
  int nipi0=0; //why is this int?
  caf->SetBranchAddress("isCC",&isCC);
  event_data->SetBranchAddress("inFV",&inFV);
  caf->SetBranchAddress("Elep_reco", &Elep_reco);
  caf->SetBranchAddress("eRecoP", &eRecoP);
  caf->SetBranchAddress("eRecoN", &eRecoP);
  caf->SetBranchAddress("eRecoPip", &eRecoPip);
  caf->SetBranchAddress("eRecoPim", &eRecoPim);
  caf->SetBranchAddress("eRecoPi0", &eRecoPi0);
  caf->SetBranchAddress("eRecoOther", &eRecoOther); //reconstructed energy
  caf->SetBranchAddress("LepE", &LepE);
  caf->SetBranchAddress("eP", &eP);
  caf->SetBranchAddress("eN", &eP);
  caf->SetBranchAddress("ePip", &ePip);
  caf->SetBranchAddress("ePi0", &ePi0);
  caf->SetBranchAddress("eOther", &eOther);
  caf->SetBranchAddress("nipi0", &nipi0);
  const double pi0_mass=0.1349768; //GeV //true energy

  //there are some non-CC events
  Long64_t nentries1=caf->GetEntries();
  Long64_t nentries2=event_data->GetEntries();
  if (nentries1!=nentries2) {cout<<"The efficiency file "<<file_num<<"has"<<nentries2
  <<" events, and the CAF file #"<<file_num<<"has"<<nentries1<<"events."<<endl;}
  for (int i=0;i<nentries2;i++) {
    caf->GetEntry(i);
    event_data->GetEntry(i); 
    int sel_cut=isCC*inFV;
    if (sel_cut==0) {continue;}
    if (sel_cut!=1) {
      cout<<"sel_cut="<<sel_cut<<endl;
      continue;
    }
    double E_vis_reco=Elep_reco+eRecoP+eRecoN+eRecoPip+eRecoPim+eRecoPi0+eRecoOther;
    double E_vis_true=LepE+eP+eN+ePip+ePim+ePi0+eOther+nipi0*pi0_mass;
    double energy_diff=E_vis_true-E_vis_reco;
    //if (abs(energy_diff)>0.5) {cout<<file_num<<","<<i<<","<<E_vis_true<<","<<E_vis_reco<<","<<energy_diff<<","<<Elep_reco<<endl;}

    if (Elep_reco<0.001) {continue;}
    heat_map->Fill(E_vis_reco,E_vis_true);
    diff->Fill(energy_diff);
  }
  eff_file.Close();
  caf_file.Close();
}

void energy_compare()
{
  char eff[99];
  char caf[99];
  TH2D *heat_map=new TH2D("energy", "reconstructed vs true particle collision energies", 1000, 0., 100., 1000, 0., 100.);
  heat_map->GetXaxis()->SetTitle("reconstructed visible energy (GeV)");
  heat_map->GetYaxis()->SetTitle("true visible energy (GeV)");
  TH1D *ediff_hist=new TH1D("diff", "difference between reco and true energies", 1000, -5., 15.);
  ediff_hist->GetXaxis()->SetTitle("energy difference (GeV)");
  ediff_hist->GetYaxis()->SetTitle("# of events");

  for (int j=0; j<=9999; j++)
  {
    memset(eff, 0, 99); // clear array each time
    memset(caf, 0, 99);
    sprintf(eff,"/storage/shared/barwu/9thTry/eff_trees/FHC.100%04d.CAF_MuonEff.root",j);
    sprintf(caf,"/storage/shared/cvilela/CAF/ND_v7/0%01d/FHC.100%04d.CAF.root",j/1000,j);
    if(access(eff, 0)==0)
    {
      populate_histograms(eff,caf,heat_map,ediff_hist,j);
    } else {
      cout<<"Warning: missing file:"<<eff<<endl;
      continue;
    }
  }
  TCanvas *c1=new TCanvas("c1", "energy distribution", 1200, 1200);
  heat_map->Draw("COL");
  TCanvas *c2=new TCanvas("c2", "energy difference", 2000, 1000);
  ediff_hist->Draw("COL");
}