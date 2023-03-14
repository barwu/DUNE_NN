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

//float max_value=0.;
vector<vector<vector<double>>>* xyz_mom=nullptr;
const int NUM_VTX=22, NUM_LAR_DTR=15;
double LAr_position[NUM_LAR_DTR]={-2800.,-2575.,-2400.,-2175.,-2000.,-1775.,-1600.,-1375.,-1200.,-975.,-800.,-575.,-400.,-175.,0.};
double vertex_position[NUM_VTX]={-299.,-292.,-285.,-278.,-271.,-264.,-216.,-168.,-120.,-72.,-24.,24.,72.,120.,168.,216.,264.,271.,
                      278.,285.,292.,299.};

void populate_histograms(char* caf,TH1D* hists[22],int j)
{
  TFile caf_file(caf);
  TTree *caf_tree=(TTree*)caf_file.Get("effTreeND");
  gInterpreter->GenerateDictionary("vector<vector<vector<double>>>", "vector");
  caf_tree->SetBranchAddress("ND_OffAxis_Sim_mu_start_p_xyz_LAr", &xyz_mom);
  Long64_t nentries1=caf_tree->GetEntries();
  for (int i=0;i<nentries1;i++) {
    caf_tree->GetEntry(i);
    unsigned long lar_pos=14;
    for (unsigned long vtx_pos=0;vtx_pos<NUM_VTX;vtx_pos++) {
      TH1D* hist=hists[vtx_pos];
      double cross_sectional_angle=180/3.14159265358979323846*atan((*xyz_mom)[lar_pos][vtx_pos][1]/(*xyz_mom)[lar_pos][vtx_pos][2]);
      hist->Fill(cross_sectional_angle);
    }
  }
  caf_file.Close();
}

void FD_muon_angles()
{
  char caf[99];
  TH1D* pos_selected_hists[22];
  for (int vtx_pos=0;vtx_pos<NUM_VTX;vtx_pos++)
  {
    pos_selected_hists[vtx_pos]=new TH1D(Form("c%d", vtx_pos+1), Form("cross-sectional angle at %f cm",vertex_position[vtx_pos]),
    200, -100, 100);
  }

  for (int j=0; j<10; j++)
  { // clear array each time
    memset(caf, 0, 99);
    sprintf(caf, "/storage/shared/fyguo/FDGeoEff_nnhome/FDGeoEff_62877585_99%d.root", j);
    populate_histograms(caf,pos_selected_hists,j);
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
    hist->GetXaxis()->SetTitle("Muon angle (degrees)");
    hist->GetYaxis()->SetTitle("# of events");
  }
  c->Update();
  c->SaveAs("/home/barwu/repos/MuonEffNN/10thTry/FD_on-axis_muon_angle_hists_degrees.png");
}