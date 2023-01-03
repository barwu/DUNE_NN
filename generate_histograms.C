#include <iostream>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <algorithm>
using namespace std;
//#include "TString.h"
//#include "TH1D.h"
//#include "TCanvas.h"
int run(int first, int last);

struct Branch 
{
  char field[50];
  char detector[30];
  double l;
  double h;
  char expr[500];
};

void generate_histograms()
{
  run(0, 9999);
}

void h12ascii (TH1* h)
{
  FILE *f=fopen(Form("csv_files/%s.csv", h->GetName()), "w");
  Int_t n=h->GetNbinsX();

  fprintf(f, "bin,content\n");
  for (Int_t i=1; i<=n; i++) {
    fprintf(f, "%g,%g\n",
    h->GetBinLowEdge(i)+h->GetBinWidth(i)/2,
    h->GetBinContent(i));
  }
  fclose(f);
}

int run(int first, int last)
{
   char buff[100];
   char fnum[10];
   TChain chain("event_data");
   for (int i=first; i<=last; i++)
   {
          memset(buff, 0, 100); // clear array each time
   	  //sprintf(fnum, "%04d", i); //get home directory
          //strcpy(buff,getenv("HOME"));
   	  //strcat(buff, "/repos/MuonEffNN/8thTry/outTrees30test2/FHC.100");
   	  //strcat(buff, fnum);
   	  //strcat(buff, ".CAF_MuonEff.root");
    	  sprintf(buff,"%s/repos/MuonEffNN/8thTry/outTrees30test2/FHC.100%04d.CAF_MuonEff.root", getenv("HOME"), i);

        if(access(buff, 0)==0)
	{
   	   chain.Add(buff);
        } else {
           cout<<"Warning: missing file:"<<buff<<endl;
           continue;
        }
   }

Branch  br[]=
{  
{"vtx_x", "contained", -200., 200., ""},
{"vtx_y", "contained", -100., 100., ""},
{"vtx_z", "contained", 50, 450., ""},
{"LepMomX", "contained", -3, 3, ""},
{"LepMomY", "contained", -6, 2, ""},
{"LepMomZ", "contained", -1, 4, ""},
{"TotLepMom", "contained", 4, 5, "sqrt(LepMomX**2+LepMomY**2+LepMomZ**2)"},
{"vtx_x", "tracker", -200., 200., ""},
{"vtx_y", "tracker", -100., 100., ""},
{"vtx_z", "tracker", 50, 450., ""},
{"LepMomX", "tracker", -3, 3, ""},
{"LepMomY", "tracker", -6, 2, ""},
{"LepMomZ", "tracker", -1, 4, ""},
{"TotLepMom", "tracker", 4, 5, "sqrt(LepMomX**2+LepMomY**2+LepMomZ**2)"},
{"vtx_x", "undetected", -200., 200., ""},
{"vtx_y", "undetected", -100., 100., ""},
{"vtx_z", "undetected", 50, 450., ""},
{"LepMomX", "undetected", -3, 3, ""},
{"LepMomY", "undetected", -6, 2., ""},
{"LepMomZ", "undetected", -1, 4, ""},
{"TotLepMom", "undetected", 4, 5, "sqrt(LepMomX**2+LepMomY**2+LepMomZ**2)"},
};

//chain.Print();
// chain.Process("myselect.C", "something");
//TCanvas* c1=new TCanvas("c1", "test", 1600, 1000);
TTree *event_data=&chain;
//int entries=chain.GetEntries();
TCanvas* c1=new TCanvas("c1", "muon neuro-training - histogram", 1800, 1000);
TCanvas* c2=new TCanvas("c2", "box-matched selection", 1800, 1000);
TCanvas* c3=new TCanvas("c3", "tracker-matched selection", 1800, 1000);
TCanvas* c4=new TCanvas("c4", "undetected selection", 1800, 1000);
TCanvas* c5=new TCanvas("c5", "detected selection", 1800, 1000);
c1->Divide(7,3);
c2->Divide(3,3);
c3->Divide(3,3);
c4->Divide(3,3);
c5->Divide(3,3);
int j=0;
//gStyle->SetOptStat(000000000);
gStyle->SetOptStat(111111111);

for(auto item:br)
{
  char *fd=item.field;
  char *dt=item.detector;
  char *expr=item.expr; 
  double l=item.l;
  double h=item.h;
  c1->cd(++j);
  char hg1[100], hg2[100], hg3[100];
  if (strlen(expr)==0) expr=fd;

//I think all the data has isCC=1, but not sure
//hg1=detector selection with geometric-efficiency corrections, hg2=no filter, hg3=detector selection
  event_data->SetLineColor(kBlue);
  sprintf(hg1, "h1_%s_%s", fd, dt);
  //cout<<fd<<hg1<<endl;
  event_data->Draw(Form("%s>>%s(100,%f,%f)", expr, hg1, l, h), Form("(abs(%s_eff)<0.0001)?0:isCC*new_%s/%s_eff", dt, dt, dt), "HIST");
  //event_data->Draw(Form("%s>>%s(100,%f,%f)", fd, hg1, l, h), "isCC", "HIST");
  //event_data->Draw(Form("%s>>%s(100,%f,%f)", fd, hg1, l, h), Form("%s_eff", dt), "HIST");
  //event_data->Draw(Form("%s>>%s(100,%f,%f)", fd, hg1, l, h), Form("(abs(%c%s_eff)<0.0001)?0:isCC*muon_%c%s/%c%s_eff", tolower(dt[0])), "HIST");
  //event_data->Draw(Form("%s>>%s(100,%f,%f)", fd, hg1, l, h), "", "HIST");
  event_data->SetLineColor(kRed);
  //cout<<fd<<hg2<<endl;
  sprintf(hg2, "h2_%s_%s", fd, dt);
  //event_data->Draw(Form("%s>>%s", fd, hg2), Form("(abs(%s_eff)<0.0001)?0:isCC*new_%s/%s_eff", dt, dt, dt), "sameHIST"); 
  //event_data->Draw(Form("%s>>%s", fd, hg2), Form("isCC*new_%c%s*(abs(%c%s_eff)<0.0001)/%c%s_eff", tolower(dt[0]), dt+1, tolower(dt[0]), dt+1), "same HIST");
  event_data->Draw(Form("%s>>%s", expr, hg2), "isCC", "sameHIST");
  //event_data->Draw(Form("%s>>%s", fd, hg2), Form("isCC*new_%c%s*(abs(%c%s_eff)<0.0001)/%c%s_eff", tolower(dt[0]), dt+1, tolower(dt[0]), dt+1), "same HIST");
  //event_data->Draw(Form("%s>>%s", fd, hg2), "", "same HIST");
  event_data->SetLineColor(kGreen);
  //cout<<fd<<hg3<endl;
  sprintf(hg3, "h3_%s_%s", fd, dt);
  event_data->Draw(Form("%s>>%s", expr, hg3), Form("isCC*new_%s", dt), "sameHIST");
  TH1D* histogram1=(TH1D*) gDirectory->Get(hg1);
  TH1D* histogram2=(TH1D*) gDirectory->Get(hg2);
  TH1D* histogram3=(TH1D*) gDirectory->Get(hg3);
  histogram1->SetTitle(Form("%s-%s distribution", fd, dt));
  //h12ascii(histogram1);
  //h12ascii(histogram2);
  //h12ascii(histogram3);

  TLegend *leg=new TLegend(0.15,0.8,0.38,0.9);
  leg->SetHeader("comparison"); 
  leg->AddEntry(hg1, Form("nn_%s_corrected_Eff", dt)); 
  leg->AddEntry(hg2, Form("%s_unselected", dt));
  leg->AddEntry(hg3, Form("%s_selection", dt));
  leg->Draw();

  float max1=histogram1->GetMaximum();
  //float min1=histogram1->GetMinimum();
  float max2=histogram2->GetMaximum();
  //float min2=histogram2->GetMinimum();
  float max3=histogram3->GetMaximum();
  //float min3=histogram3->GetMinimum();
  double upper_y_bound1=std::max(std::max(max2,max3), max1)*1.2;
  histogram1->SetAxisRange(0.,upper_y_bound1,"Y");

  if (j<8)
    {c2->cd(j);}
  else if (j>14)
    {c4->cd(j-14);}
  else
    {c3->cd(j-7);}
//if statements distribute copies of the histograms to 3 new canvases, if you change the original histograms that will appear on these, too.
//Tried using switch, caused a bug

  TH1D* histogram_a=(TH1D*)histogram1->Clone();
  TH1D* histogram_b=(TH1D*)histogram2->Clone();
  TH1D* histogram_c=(TH1D*)histogram3->Clone();
  histogram_a->Draw("HIST");
  histogram_b->Draw("sameHIST");
  histogram_c->Draw("sameHIST");
  histogram_c->SetAxisRange(0.,upper_y_bound1,"Y");  
/*
  auto rp=new TRatioPlot(histogram2, histogram3);
  rp->SetH1DrawOpt("HIST");
  rp->SetH2DrawOpt("HIST");
  rp->SetGraphDrawOpt("AC");
  rp->SetSeparationMargin(0.0);
  rp->Draw();
  //rp->GetLowerRefGraph()->SetMinimum(0.5);
  //rp->GetLowerRefGraph()->SetMaximum(1.5);

  gPad->Modified();
  gPad->Update();
  TPad *p=rp->GetUpperPad();
  TLegend *leg=p->BuildLegend();
  leg->Clear();

  leg=new TLegend(0.15,0.8,0.38,0.9);
  leg->SetHeader("comparison"); 
  leg->AddEntry(hg1, Form("nn_%s_corrected_Eff", dt)); 
  leg->AddEntry(hg2, Form("%s_unselected", dt));
  leg->AddEntry(hg3, Form("%s_selection", dt));
  leg->Draw();
*/
//detected (contained_tracker), separate section of code
  if (j>7) continue;
  char hg4[100], hg5[100], hg6[100];
  c5->cd(j);
  event_data->SetLineColor(kBlue);
  sprintf(hg4, "h1_%s_detected", fd);
  //event_data->Draw(Form("%s>>%s(100,%f,%f)", fd, hg4, l, h), "((abs(contained_eff)<0.0001)?0:isCC*new_contained/contained_eff)+((abs(tracker_eff)<0.0001)?0:isCC*new_tracker/tracker_eff)", "HIST");
  event_data->Draw(Form("%s>>%s(100,%f,%f)", expr, hg4, l, h), "((abs(contained_eff)<0.0001)&&(abs(tracker_eff)<0.0001))?0:isCC*(new_contained+new_tracker)/(contained_eff+tracker_eff)", "HIST");
  event_data->SetLineColor(kRed);
  sprintf(hg5, "h2_%s_detected", fd);
  event_data->Draw(Form("%s>>%s", expr, hg5), "isCC", "sameHIST");
  event_data->SetLineColor(kGreen);
  sprintf(hg6, "h3_%s_detected", fd);
  event_data->Draw(Form("%s>>%s", expr, hg6), "isCC*(new_contained+new_tracker)", "sameHIST");
  TH1D* histogram4=(TH1D*) gDirectory->Get(hg4);
  TH1D* histogram5=(TH1D*) gDirectory->Get(hg5);
  TH1D* histogram6=(TH1D*) gDirectory->Get(hg6);
  histogram6->SetTitle(Form("%s-detector distribution", fd));

  leg=new TLegend(0.15,0.8,0.38,0.9);
  leg->SetHeader("comparison"); 
  leg->AddEntry(hg4, "nn_detected_corrected_Eff"); 
  leg->AddEntry(hg5, "detected_selected");
  leg->AddEntry(hg6, "detected_selection");
  leg->Draw();

  float max4=histogram4->GetMaximum();
  //float min1=histogram1->GetMinimum();
  float max5=histogram5->GetMaximum();
  //float min2=histogram2->GetMinimum();
  float max6=histogram6->GetMaximum();
  //float min3=histogram3->GetMinimum();
  double upper_y_bound2=std::max(std::max(max6,max5), max4)*1.2;
  histogram4->SetAxisRange(0.,upper_y_bound2,"Y");
}

c1->Update();
c2->Update();
c3->Update();
c4->Update();
c5->Update();
//c1->SaveAs("~/repos/MuonEffNN/8thTry/test.png");
//can't seem to auto-save if I have multiple canvases.
return 0;
}
