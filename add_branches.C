#include "iostream"
#include "algorithm"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"
using namespace std;
void upd(const char* new_file,const char* backup_file)
{
  TFile *read_file=new TFile(backup_file,"read"); 
  TTree *event_data=(TTree*)read_file->Get("event_data");
  double LepMomX, LepMomY, LepMomZ, LepNuAngle;
  double TotMom_dat, cos_LepNuAngle_dat, LongMom_dat;
  TFile *write_file=new TFile(new_file,"recreate");
  TTree *T=(TTree*)event_data->CloneTree();

  TBranch *TotMom=T->Branch("TotMom",&TotMom_dat,"TotMom_dat/D");
  TBranch *cos_LepNuAngle=T->Branch("cos_LepNuAngle",&cos_LepNuAngle_dat,"cos_LepNuAngle_dat/D");
  TBranch *LongMom=T->Branch("LongMom",&LongMom_dat,"LongMom_dat/D");
  T->SetBranchAddress("LepMomX",&LepMomX);
  T->SetBranchAddress("LepMomY",&LepMomY);
  T->SetBranchAddress("LepMomZ",&LepMomZ);
  T->SetBranchAddress("LepNuAngle",&LepNuAngle);
  Long64_t nentries=T->GetEntries();
  for (Long64_t i=0;i<nentries;i++)
  {
    T->GetEntry(i);
    cos_LepNuAngle_dat=cos(LepNuAngle);
    cos_LepNuAngle->Fill();
    TotMom_dat=sqrt(pow(LepMomX,2)+pow(LepMomY,2)+pow(LepMomZ,2));
    TotMom->Fill();
    LongMom_dat=cos_LepNuAngle_dat*TotMom_dat;
    LongMom->Fill();
  }

  //T->Print();
  T->Write();
  delete read_file;
  delete write_file;
}

int main()
{
   char buff[99];
   char backup[99];
   TChain chain("event_data");
   for (int i=0; i<=9999; i++)
   {
     memset(buff, 0, 100); // clear array each time
     memset(backup, 0, 100);
     sprintf(buff,"%s/repos/MuonEffNN/9thTry/combined1/FHC.100%04d.CAF_MuonEff.root", getenv("HOME"), i);
     sprintf(backup,"%s/repos/MuonEffNN/9thTry/combined2/FHC.100%04d.CAF_MuonEff.root", getenv("HOME"), i);

     if(access(backup, 0)==0)
     {
   	upd(buff,backup);
        //chain.Add(buff);
     } else {
        cout<<"Warning: missing file:"<<buff<<endl;
        continue;
     }
   }
return 0;
}