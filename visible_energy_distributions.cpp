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

struct Para
{
  char field[30];
  double l;
  double h;
};
int muon_cont, muon_tra, muon_sel, hadr, comb;
double muon_cont_eff, muon_tra_eff, muon_sel_eff, hadr_eff, comb_eff;
Para pr[]=
{
  {"E_vis_true",0,6},
  {"Ev",0,6}
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

vector<Sel_type> br=
{
  Sel_type("muon_contained", "muon_contained_eff", false, &muon_cont, &muon_cont_eff),
  Sel_type("muon_tracker", "muon_tracker_eff", false, &muon_tra, &muon_tra_eff),
  Sel_type("muon_selected", "muon_sel_eff", true, &muon_sel, &muon_sel_eff),
  Sel_type("hadron_selected", "hadron_selected_eff", false, &hadr, &hadr_eff ),
  Sel_type("combined", "combined_eff", false, &comb, &comb_eff)
};

void populate_histograms(char* eff,char* CAF,TH1D* hists[3][10],int j)
{
  TFile eff_file(eff);
  TFile caf_file(CAF);
  TTree *event_data=(TTree*)eff_file.Get("event_data");
  TTree *caf=(TTree*)caf_file.Get("caf");
  int isCC=0, inFV=0;
  for(auto item:br) 
  {
    if(item.calced) continue;
    event_data->SetBranchAddress(item.sel_name, item.sel_value);
    event_data->SetBranchAddress(item.eff_name, item.eff_value);
  }
  //double Elep_reco=0., eRecoP=0., eRecoPip=0., eRecoPim=0., eRecoPi0=0., eRecoOther=0.;
  double LepE=0., eP=0., ePip=0., ePim=0., ePi0=0., eOther=0.;
  int nipi0=0; //why is this int?
  double Ev=0.;
  caf->SetBranchAddress("isCC",&isCC);
  event_data->SetBranchAddress("inFV",&inFV);
  // caf->SetBranchAddress("Elep_reco", &Elep_reco);
  // caf->SetBranchAddress("eRecoP", &eRecoP);
  // caf->SetBranchAddress("eRecoPip", &eRecoPip);
  // caf->SetBranchAddress("eRecoPim", &eRecoPim);
  // caf->SetBranchAddress("eRecoPi0", &eRecoPi0);
  // caf->SetBranchAddress("eRecoOther", &eRecoOther); //reconstructed energy
  caf->SetBranchAddress("LepE", &LepE);
  caf->SetBranchAddress("eP", &eP);
  caf->SetBranchAddress("ePip", &ePip);
  caf->SetBranchAddress("ePi0", &ePi0);
  caf->SetBranchAddress("eOther", &eOther);
  caf->SetBranchAddress("nipi0", &nipi0);
  const double pi0_mass=0.1349768; //GeV //true energy
  caf->SetBranchAddress("Ev", &Ev); //true neutrino energy

  //there are some non-CC events
  Long64_t nentries1=caf->GetEntries();
  Long64_t nentries2=event_data->GetEntries();
  if (nentries1!=nentries2) {cout<<"The efficiency file #"<<j<<"has"<<nentries2
  <<" events, and the CAF file #"<<j<<"has"<<nentries1<<"events."<<endl;}
  for (int i=0;i<nentries2;i++) {
    caf->GetEntry(i);
    event_data->GetEntry(i); 
    int sel_cut=isCC*inFV;
    if (sel_cut==0) {continue;}
    if (sel_cut!=1) {cout<<"sel_cut="<<sel_cut<<endl;
      continue;
    }
    //if (Elep_reco<0.001) {continue;} //energy should always be positive
    //calculation for the muon-selected cut
    muon_sel=muon_cont+muon_tra;
    muon_sel_eff=muon_cont_eff+muon_tra_eff;
    if (muon_sel!=0&&muon_sel!=1) {
      cout<<"bad val for muon-selected check!"<<muon_sel<<endl;
      continue;
    }
    //if (muon_sel_eff<0.0001) {cout<<"event "<<i<<" of file #"<<j<<" has small selected efficiency"<<endl;}
    //double E_vis_reco=Elep_reco+eRecoP+eRecoPip+eRecoPim+eRecoPi0+eRecoOther;
    double E_vis_true=LepE+eP+ePip+ePim+ePi0+eOther+nipi0*pi0_mass;

    //double energy[2]={E_vis_reco,Ev_reco};
    double energy[2]={E_vis_true,Ev};
    int n=0;
    for (int k=0; k<2; k++) {
      for (auto sel:br) {
        TH1D* hist1=hists[0][n];
        TH1D* hist2=hists[1][n];
        TH1D* hist3=hists[2][n];
        n++;
	      //if (Elep_reco<0.001&&k==0) {
        //  hist1->Fill(energy[k],0);
        //} else {
        hist1->Fill(energy[k]);
        //}
        double geo_eff=*sel.eff_value;
        if (geo_eff<=0.1) {
          hist2->Fill(energy[k],0.);
	        hist3->Fill(energy[k],0.);
	      } else {
          hist2->Fill(energy[k],*sel.sel_value);
	        hist3->Fill(energy[k],*sel.sel_value/geo_eff);
	      }
      }
    }
  }
  eff_file.Close();
  caf_file.Close();
}

void visible_energy_distributions()
{
  char eff[99];
  char caf[99];
  TH1D *histograms[3][10];
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
      histograms[1][m]=new TH1D(Form("h2_%s_%s", fd, dt), Form("selected %s %s",fd, dt), 100, l, h);
      histograms[2][m]=new TH1D(Form("h3_%s_%s", fd, dt), Form("geo corrected %s %s",fd, dt), 100, l, h);
      m++;
    }
  }

  for (int j=0; j<=9999; j++)
  {
    memset(eff, 0, 99); // clear array each time
    memset(caf, 0, 99);
    sprintf(eff,"/storage/shared/barwu/9thTry/eff_trees/FHC.100%04d.CAF_MuonEff.root",j);
    sprintf(caf,"/storage/shared/cvilela/CAF/ND_v7/0%01d/FHC.100%04d.CAF.root",j/1000,j);
    if(access(eff, 0)==0)
    {
      populate_histograms(eff,caf,histograms,j);
    } else {
      cout<<"Warning: missing file:"<<eff<<endl;
      continue;
    }
  }

  gStyle->SetOptStat(000000000);
  TCanvas *cs[2];
  TCanvas *rs[2];
  int n=0;
  for(int i=0;i<2;i++)
  {
    Para item=pr[i];
    const char *fd=item.field;
    cs[i]=new TCanvas(Form("c%01d",i+1),fd,1800,1000);
    cs[i]->Divide(3,2);
    rs[i]=new TCanvas(Form("r%01d",i+1),fd,1800,1000);
    rs[i]->Divide(3,2);
    int k=0;
    for(auto sel:br)
    {
      const char *dt=sel.sel_name;
      TVirtualPad *p=cs[i]->cd(k+1);
      TH1D *hist3=histograms[2][n];
      hist3->SetLineColor(kAzure);
      hist3->Draw("histS");
      TH1D *hist2=histograms[1][n];
      hist2->SetLineColor(kSpring);
      hist2->Draw("samehistS");
      TH1D *hist1=histograms[0][n];
      hist1->SetLineColor(kPink);
      hist1->Draw("samehistS");

      float max1=hist1->GetMaximum();
      float max2=hist2->GetMaximum();
      float max3=hist3->GetMaximum();
      float upper_y_bound=max(max(max2,max3), max1)*1.2;
      hist3->SetAxisRange(0.,upper_y_bound,"Y");
      hist3->SetTitle(Form("%s: %s",fd,dt));
      hist3->GetXaxis()->SetTitle(Form("%s (GeV)",fd));
      hist3->GetYaxis()->SetTitle("# of events");
      TLegend *leg=new TLegend(0.1,0.8,0.35,0.9);
      leg->SetHeader("comparison"); 
      leg->AddEntry(hist1, "raw distribution");
      leg->AddEntry(hist2, "selection-cut distribution");
      leg->AddEntry(hist3, "geo corrected distribution");
      leg->Draw();

      rs[i]->cd(k+1);
      TH1D *rplot=(TH1D*)hist3->Clone();
      rplot->Divide(hist1);
      rplot->SetAxisRange(0.,1.,"Y");
      rplot->SetLineColor(kViolet);
      rplot->Draw("hist");
      n++;
      k++;
    }
    cs[i]->Update();
    rs[i]->Update();
    cs[i]->SaveAs(Form("images/true_%s_distributions_with_veto_cut.png",fd));
    rs[i]->SaveAs(Form("images/true_%s_correction_ratios_with_veto_cut.png",fd));
  }
}