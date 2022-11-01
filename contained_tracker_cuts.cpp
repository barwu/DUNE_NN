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
  char event_cut[30];
  double l;
  double h;
  int key;
};
int muon_cont, muon_tra;//, muon_sel, hadr, comb;
double muon_cont_eff, muon_tra_eff;//, muon_sel_eff, hadr_eff, comb_eff;
Para pr[]=
{
  {"",0,8,1},
  {"veto_cut_passed",0,8,2},
  {"veto_cut_failed",0,8,3},
  {"eff_cut",0,8,4},
};

struct Sel_type
{
  const char* sel_name;
  const char* eff_name;
  //bool calced=false;
  int* sel_value;
  double* eff_value;
  Sel_type() {}
  Sel_type(const char* sn, const char* en,  int* sv, double* ev)
  :sel_name(sn),eff_name(en),sel_value(sv),eff_value(ev) {}
};
vector<Sel_type> br=
{
  Sel_type("muon_contained", "muon_contained_eff", &muon_cont, &muon_cont_eff),
  Sel_type("muon_tracker", "muon_tracker_eff", &muon_tra, &muon_tra_eff)
  //Sel_type("muon_selected", "muon_sel_eff", true, &muon_sel, &muon_sel_eff),
  //Sel_type("hadron_selected", "hadron_selected_eff", false, &hadr, &hadr_eff ),
  //Sel_type("combined", "combined_eff", false, &comb, &comb_eff)
};

void populate_histograms(char* eff,char* CAF,TH1D* hists[6][8],int j)
{
  TFile eff_file(eff);
  TFile caf_file(CAF);
  TTree *event_data=(TTree*)eff_file.Get("event_data");
  TTree *caf=(TTree*)caf_file.Get("caf");
  int isCC=0, inFV=0;
  for(auto item:br) 
  {
    event_data->SetBranchAddress(item.sel_name, item.sel_value);
    event_data->SetBranchAddress(item.eff_name, item.eff_value);
  }
  double LepE=0., eP=0., ePip=0., ePim=0., ePi0=0., eOther=0.;
  int nipi0=0; //why is this int?
  double Ev=0.;
  caf->SetBranchAddress("isCC",&isCC);
  event_data->SetBranchAddress("inFV",&inFV);
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
    if (sel_cut!=1) {
      cout<<"sel_cut="<<sel_cut<<endl;
      continue;
    }
    double E_vis_true=LepE+eP+ePip+ePim+ePi0+eOther+nipi0*pi0_mass;
    double energy[2]={E_vis_true,Ev};

    int n=0;
    for (auto sel:br) {
      double geo_eff=*sel.eff_value;
      for (auto item:pr) {
        TH1D* hist1=hists[0][n];
        TH1D* hist2=hists[1][n];
        TH1D* hist3=hists[2][n];
        TH1D* hist4=hists[3][n];
        TH1D* hist5=hists[4][n];
        TH1D* hist6=hists[5][n];
        n++;
        int key=item.key;
        hist1->Fill(energy[0]);
        hist4->Fill(energy[1]);

        if (key==1) {
          if (geo_eff<0.0001) {
            hist2->Fill(energy[0],0.);
	          hist3->Fill(energy[0],0.);
            hist5->Fill(energy[1],0.);
	          hist6->Fill(energy[1],0.);
	        } else {
            hist2->Fill(energy[0],*sel.sel_value);
	          hist3->Fill(energy[0],*sel.sel_value/geo_eff);
            hist5->Fill(energy[1],*sel.sel_value);
	          hist6->Fill(energy[1],*sel.sel_value/geo_eff);
	        }
        }

        if (key==2) {
          if (geo_eff<=0.1) {
            hist2->Fill(energy[0],0.);
	          hist3->Fill(energy[0],0.);
            hist5->Fill(energy[1],0.);
	          hist6->Fill(energy[1],0.);
          } else {
            hist2->Fill(energy[0],*sel.sel_value);
	          hist3->Fill(energy[0],*sel.sel_value/geo_eff);
            hist5->Fill(energy[1],*sel.sel_value);
	          hist6->Fill(energy[1],*sel.sel_value/geo_eff);
	        }
        }
        
        if (key==3) {
          if (geo_eff>0.1||geo_eff<0.0001) {
            hist2->Fill(energy[0],0.);
	          hist3->Fill(energy[0],0.);
            hist5->Fill(energy[1],0.);
	          hist6->Fill(energy[1],0.);
          } else {
            hist2->Fill(energy[0],*sel.sel_value);
	          hist3->Fill(energy[0],*sel.sel_value/geo_eff);
            hist5->Fill(energy[1],*sel.sel_value);
	          hist6->Fill(energy[1],*sel.sel_value/geo_eff);
	        }
        }

        if (key==4) {
	        if (n<5) {
	          if (muon_tra_eff>muon_cont_eff) {
              hist2->Fill(energy[0],0.);
              hist3->Fill(energy[0],0.);
              hist5->Fill(energy[1],0.);
              hist6->Fill(energy[1],0.);
            } else {
	            //cout<<"events in top graph= " << muon_tra_eff << "," << muon_cont_eff <<endl;
              hist2->Fill(energy[0],muon_cont);
              hist3->Fill(energy[0],muon_cont/muon_cont_eff);
              hist5->Fill(energy[1],muon_cont);
              hist6->Fill(energy[1],muon_cont/muon_cont_eff);
            }
	        } else {
	          if (muon_cont_eff>muon_tra_eff) {
              hist2->Fill(energy[0],0.);
              hist3->Fill(energy[0],0.);
              hist5->Fill(energy[1],0.);
              hist6->Fill(energy[1],0.);
            } else {
              //cout<<"events bt graph= " << muon_tra_eff << "," << muon_cont_eff <<endl;
              hist2->Fill(energy[0],muon_tra);
              hist3->Fill(energy[0],muon_tra/muon_tra_eff);
              hist5->Fill(energy[1],muon_tra);
              hist6->Fill(energy[1],muon_tra/muon_tra_eff);
          	}
	        }
        }
      }
    }
  }
  eff_file.Close();
  caf_file.Close();
}

void contained_tracker_cuts()
{
  char eff[99];
  char caf[99];
  TH1D *histograms[6][8];
  int m=0;
  for(auto sel:br)
  {
    const char* dt=sel.sel_name;
    for(auto item:pr)
    {
      char *fd=item.event_cut;
      double l=item.l;
      double h=item.h;
      histograms[0][m]=new TH1D(Form("h1_%s_E_vis_true %s", fd, dt), Form("raw %s %s",fd, dt), 100, l, h);
      histograms[1][m]=new TH1D(Form("h2_%s_E_vis_true %s", fd, dt), Form("selected %s %s",fd, dt), 100, l, h);
      histograms[2][m]=new TH1D(Form("h3_%s_E_vis_true %s", fd, dt), Form("geo corrected %s %s",fd, dt), 100, l, h);
      histograms[3][m]=new TH1D(Form("h1_%s_Ev %s", fd, dt), Form("raw %s %s",fd, dt), 100, l, h);
      histograms[4][m]=new TH1D(Form("h2_%s_Ev %s", fd, dt), Form("selected %s %s",fd, dt), 100, l, h);
      histograms[5][m]=new TH1D(Form("h3_%s_Ev %s", fd, dt), Form("geo corrected %s %s",fd, dt), 100, l, h);
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

  //gStyle->SetOptStat(000000000);
  TCanvas *cs[4];
  cs[0]=new TCanvas("c1","E_vis_true cuts",1800,1000);
  cs[1]=new TCanvas("c2","Ev cuts",1800,1000);
  cs[2]=new TCanvas("c3","E_vis_true ratio plots",1800,1000);
  cs[3]=new TCanvas("c4","Ev ratio plots",1800,1000);
  for (int i=0;i<4;i++) {cs[i]->Divide(4,2);}
  int n=0;
  for(int j=0;j<2;j++) {
    const char *dt=br[j].sel_name;
    for(int k=0;k<4;k++)
    {
      const char *fd=pr[k].event_cut;

      cs[0]->cd(n+1);
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
      float upper_y_bound1=max(max(max2,max3),max1)*1.2;
      hist3->SetAxisRange(0.,upper_y_bound1,"Y");
      hist3->SetTitle(Form("E_vis_true: %s %s",dt,fd));
      hist3->GetXaxis()->SetTitle("Energy (GeV)");
      hist3->GetYaxis()->SetTitle("# of events");
      TLegend *leg1=new TLegend(0.1,0.8,0.35,0.9);
      leg1->SetHeader("comparison"); 
      leg1->AddEntry(hist1, "raw distribution");
      leg1->AddEntry(hist2, "selection-cut distribution");
      leg1->AddEntry(hist3, "geo corrected distribution");
      leg1->Draw();

      cs[1]->cd(n+1);
      TH1D *hist6=histograms[5][n];
      hist6->SetLineColor(kAzure);
      hist6->Draw("histS");
      TH1D *hist5=histograms[4][n];
      hist5->SetLineColor(kSpring);
      hist5->Draw("samehistS");
      TH1D *hist4=histograms[3][n];
      hist4->SetLineColor(kPink);
      hist4->Draw("samehistS");

      float max4=hist4->GetMaximum();
      float max5=hist5->GetMaximum();
      float max6=hist6->GetMaximum();
      float upper_y_bound2=max(max(max5,max6),max4)*1.2;
      hist6->SetAxisRange(0.,upper_y_bound2,"Y");
      hist6->SetTitle(Form("Ev: %s %s",dt,fd));
      hist6->GetXaxis()->SetTitle("Energy (GeV)");
      hist6->GetYaxis()->SetTitle("# of events");
      TLegend *leg2=new TLegend(0.1,0.8,0.35,0.9);
      leg2->SetHeader("comparison");
      leg2->AddEntry(hist4, "raw distribution");
      leg2->AddEntry(hist5, "selection-cut distribution");
      leg2->AddEntry(hist6, "geo corrected distribution");
      leg2->Draw();

      cs[2]->cd(n+1);
      TH1D *rplot1=(TH1D*)hist3->Clone();
      rplot1->Divide(hist1);
      rplot1->SetAxisRange(0.,max3/max1*1.5,"Y");
      rplot1->SetLineColor(kBlack);
      rplot1->Draw("hist");
      cs[3]->cd(n+1);
      TH1D *rplot2=(TH1D*)hist6->Clone();
      rplot2->Divide(hist4);
      rplot2->SetAxisRange(0.,max6/max4*1.5,"Y");
      rplot2->SetLineColor(kBlack);
      rplot2->Draw("hist");
      n++;
    }
  }
  for (int i=0;i<4;i++) {cs[i]->Update();}
  cs[0]->SaveAs("energy_cuts/E_vis_true_cuts5.png");
  cs[1]->SaveAs("energy_cuts/Ev_cuts5.png");
  cs[2]->SaveAs("energy_cuts/E_vis_true_ratio_plots5.png");
  cs[3]->SaveAs("energy_cuts/Ev_ratio_plots5.png");
}