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

struct Para
{
  char field[30];
  double l;
  double h;
  double* field_value;
};
double ND_E_vis_true=0., ND_Gen_numu_E=0.;
vector<vector<double>>* muon_cont_eff=nullptr;
vector<vector<double>>* muon_tra_eff=nullptr;
vector<vector<double>>* muon_sel_eff=nullptr;
vector<vector<double>>* hadr_eff=nullptr;
vector<vector<double>>* comb_eff=nullptr;
Para pr[]=
{
  {"ND_E_vis_true",0,8,&ND_E_vis_true},
  {"ND_Gen_numu_E",0,8,&ND_Gen_numu_E}
};

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

void populate_histograms(char* eff,char* caf,TH1D* hists[2][10],int j)
{
  TFile eff_file(eff);
  TFile caf_file(caf);
  TTree *event_data=(TTree*)eff_file.Get("event_data");
  TTree *detections=(TTree*)caf_file.Get("effTreeND");
  for(auto item:br) //efficiencies are in 3d array, but energy is in 1d array
  {
    if(item.calced) continue;
    event_data->SetBranchAddress(item.eff_name, &(item.eff_value));
    hadr_eff=item.eff_value;
  }
  //energy is in gev
  detections->SetBranchAddress("ND_E_vis_true", &ND_E_vis_true); //true energy
  detections->SetBranchAddress("ND_Gen_numu_E", &ND_Gen_numu_E); //true neutrino energy

  Long64_t nentries1=detections->GetEntries(); //no fv filter for fd //cc events already filtered
  Long64_t nentries2=event_data->GetEntries();
  cout<<nentries2;
  if (nentries1!=nentries2) {cout<<"the efficiency file #"<<j<<" has "<<nentries2
  <<" events, and the caf file #"<<j<<" has "<<nentries1<<" events."<<endl;}
  for (int i=0;i<nentries2;i++) {
    detections->GetEntry(i);
    event_data->GetEntry(i);
    for (int lar_pos=0;lar_pos<hadr_eff->size();lar_pos++) {
      for (int vtx_pos=0;vtx_pos<hadr_eff->at(j).size();vtx_pos++) {
        //calculation for the muon-selected cut
        //muon_sel_eff[lar_pos][vtx_pos]=muon_cont_eff[lar_pos][vtx_pos]+muon_tra_eff[lar_pos][vtx_pos];
        //if (muon_sel_eff<0.0001) {cout<<"event "<<i<<" of file #"<<j<<" has small selected efficiency"<<endl;}
	
	      int n=0;
        for (auto item:pr) {
          for (auto sel:br) {
            TH1D* hist1=hists[0][n];
            TH1D* hist2=hists[1][n];
            n++;
            hist1->Fill(*item.field_value);
            vector<vector<double>> eff_value2=*(hadr_eff);
            double geo_eff=eff_value2[lar_pos][vtx_pos];
            //cout<<"geoeff="<<geo_eff<<endl;
            if (geo_eff<=0.001) {
              hist2->Fill(*item.field_value,0.);
            } else {
              hist2->Fill(*item.field_value,1/geo_eff);
	          }
          }
        }
      }
    }
  }
  eff_file.Close();
  caf_file.Close();
}

void FD_energy_histograms()
{
  char eff[99];
  char caf[99];
  TH1D *histograms[2][10];
  int m=0;
  for(auto sel:br)
  {
    const char* dt=sel.sel_name;
    for(auto item:pr)
    {
      char *fd=item.field;
      double l=item.l;
      double h=item.h;
      histograms[0][m]=new TH1D(Form("h1_%s_%s", fd, dt), Form("raw %s %s",fd, dt), 100, l, h);
      histograms[1][m]=new TH1D(Form("h3_%s_%s", fd, dt), Form("geo corrected %s %s",fd, dt), 100, l, h);
      m++;
    }
  }

  for (int j=0; j<=9; j++)
  {
    memset(eff, 0, 99); // clear array each time
    memset(caf, 0, 99);
    sprintf(eff,"/storage/shared/barwu/10thTry/FDEff/FDGeoEff_62877585_99%01d_Eff.root",j);
    sprintf(caf,"/storage/shared/fyguo/FDGeoEff_nnhome/FDGeoEff_62877585_99%01d.root",j);
    if(access(eff, 0)==0)
    {
      populate_histograms(eff,caf,histograms,j);
    } else {
      cout<<"Warning: missing file:"<<eff<<endl;
      continue;
    }
  }

  gStyle->SetOptStat(000000000);
  TCanvas *c=new TCanvas("c","Energy Distributions",1800,1000);
  c->Divide(2,2);
  int n=0;
  for(int i=0;i<2;i++)
  {
    Para item=pr[i];
    const char *fd=item.field;
    for(auto sel:br)
    {
      const char *dt=sel.sel_name;
      c->cd(2*i+1);
      TH1D *hist2=histograms[1][n];
      hist2->SetLineColor(kAzure);
      hist2->Draw("samehistS");
      TH1D *hist1=histograms[0][n];
      hist1->SetLineColor(kPink);
      hist1->Draw("samehistS");

      float max1=hist1->GetMaximum();
      float max2=hist2->GetMaximum();
      float upper_y_bound=max(max2, max1)*1.2;
      hist2->SetAxisRange(0.,upper_y_bound,"Y");
      hist2->SetTitle(Form("%s: %s",fd,dt));
      hist2->GetXaxis()->SetTitle(Form("%s (GeV)",fd));
      hist2->GetYaxis()->SetTitle("# of events");
      TLegend *leg=new TLegend(0.1,0.8,0.35,0.9);
      leg->SetHeader("comparison"); 
      leg->AddEntry(hist1, "raw distribution");
      leg->AddEntry(hist2, "geo corrected distribution");
      leg->Draw();

      c->cd(2*i+2);
      TH1D *rplot=(TH1D*)hist2->Clone();
      rplot->Divide(hist1);
      rplot->SetTitle(Form("%s: %s ratio plot",fd,dt));
      rplot->SetAxisRange(0.,3.,"Y");
      rplot->SetLineColor(kBlue);
      rplot->Draw("hist");
      n++;
    }
    c->Update();
    //c->SaveAs("10thTry/true_energy_distributions.png");
  }
}