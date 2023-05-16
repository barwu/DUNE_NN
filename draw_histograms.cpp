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
#include "TPaveStats.h"

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
const char* list_of_directories[40]={"0mgsimple","0m","1.75m","2m","4m","5.75m","8m","9.75m","12m","13.75m","16m",
"17.75m","20m","21.75m","24m","25.75m","26.75m","28m","28.25m","28.5m","0mgsimpleRHC","0mRHC","1.75mRHC",
"2mRHC","4mRHC","5.75mRHC","8mRHC","9.75mRHC","12mRHC","13.75mRHC","16mRHC","17.75mRHC","20mRHC",
"21.75mRHC","24mRHC","25.75mRHC","26.75mRHC","28mRHC","28.25mRHC","28.5mRHC"};

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

void draw_histograms()
{
  TFile *raw_files[9];
  TFile *sel_files[45];
  TFile *geo_files[45];

  TH1D* raw_histograms[9];
  TH1D* sel_histograms[45];
  TH1D* geo_histograms[45];
  int index=0;
  int i_pr=0;
  for(Para item:pr)
  {
    TTree *sel_hists=(TTree*)raw_histograms.Get("sel_hist_sel");
    for (auto sel:br)
    {

      TTree *raw_hists=(TTree*)hist_file.Get("raw_hist_raw");
      TTree *geo_hists=(TTree*)hist_file.Get("hist_geo");
      raw_hists->SetBranchAddress(Form("raw_%s", item.field), raw_data[i_pr]);
                i_pr++;
      if(sel.calced) continue;
      for (Para item:pr)
      {
        sel_hists->SetBranchAddress(Form("selected_%s_%s", item.field,sel.sel_name), sel_data[index]);
        geo_hists->SetBranchAddress(Form("geo-corrected_%s_%s", item.field,sel.sel_name), geo_data[index]);
        index++;
      }
    }
  }

  TCanvas *cs[5];
  index=0;
  int i=0;
  for(auto sel:br)
  {
    cs[i]=new TCanvas(Form("c%01d",i+1),Form("c%01d",i+1),1800,1000);
    cs[i]->Divide(3,3);
    const char *dt=sel.sel_name;
    for(int k=0;k<9;k++)
    {
      Para item=pr[k];
      const char *fd=item.field;
      TVirtualPad *p=cs[i]->cd(k+1);
      if (k==7) {p->SetLogy();} //pad needs to be made logarithmic, not canvas
      TH1D *hist3=raw_histograms[index%9];
      hist3->SetLineColor(kBlue);
      hist3->Draw("histS");
      TH1D *hist2=sel_histograms[index];
      //hist2->SetLineColor(kGreen);
      hist2->SetLineColor(kTeal+10);
      hist2->Draw("samehistS");
      TH1D *hist1=geo_histograms[index];
      hist1->SetLineColor(kPink);
      hist1->Draw("samehistS");

      float max1=hist1->GetMaximum();
      float max2=hist2->GetMaximum();
      float max3=hist3->GetMaximum();
      float upper_y_bound=max(max(max2,max3), max1);
      if (k!=7) {hist3->SetAxisRange(0.,upper_y_bound,"Y");}
      //else {hist3->SetAxisRange(0.1,upper_y_bound,"Y");}
      hist3->SetTitle(Form("%s: %s",fd,dt));
      hist3->GetXaxis()->SetTitle(Form("%s",fd));
      hist3->GetYaxis()->SetTitle("# of events");
      TLegend *leg=new TLegend(0.1,0.75,0.33,0.9);
      leg->SetHeader("comparison"); 
      leg->AddEntry(hist1, "raw distribution");
      leg->AddEntry(hist2, "selection-cut distribution");
      leg->AddEntry(hist3, "geo corrected distribution");
      leg->Draw();
      gPad->Update();
      TPaveStats *ps;
      ps=(TPaveStats*)hist3->GetListOfFunctions()->FindObject("stats");
      ps->SetFillStyle(0);
      index++;
    }
    cs[i]->Update();
    //cs[i]->SaveAs(Form("/home/barwu/repos/MuonEffNN/images/0m_%s_PRISM_hists.png",directory_number,dt));
    i++;
  }
}