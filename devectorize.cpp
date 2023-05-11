#include <iostream>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <vector>
#include "algorithm"
#include "TFile.h"
#include "TTree.h"
#include "TInterpreter.h"
#include "TBranch.h"
using namespace std;

int main() {
    gInterpreter->GenerateDictionary("vector<vector<vector<double>>>", "vector");
    gInterpreter->GenerateDictionary("vector<vector<vector<uint64_t>>>", "vector");
    gInterpreter->GenerateDictionary("vector<vector<vector<vector<vector<uint64_t>>>>>", "vector");
    for (int i=0;i<10;i++)
    {
        const char *caf=Form("/storage/shared/barwu/FDGeoEffinND/FDGeoEff_62877585_99%01d.root",i);
        const char *newcaf=Form("/storage/shared/barwu/FDdevectorized/FDGeoEff_62877585_99%01d.root",i);
        TFile *read_file=new TFile(caf,"read");
        TTree *effTree=(TTree*)read_file->Get("effTreeND");
        TTree *throwsFD=(TTree*)read_file->Get("ThrowsFD");
        TFile *write_file=new TFile(newcaf,"recreate");
        TTree *FDCAF=new TTree("FDCAF", Form("caf_%d", i));
        TTree *T=(TTree*)throwsFD->CloneTree();
        vector<vector<vector<double>>>* xyz_pos=nullptr;
        vector<vector<vector<double>>>* xyz_mom=nullptr;
        vector<vector<vector<vector<vector<uint64_t>>>>>*  hadron_throw_result_LAr=nullptr;
        double vtx_vals[3]={};
        double lepmom_vals[3]={};
        int i_LAr_pos=0, i_vtx_x=0;
        vector<vector<vector<uint64_t>>> hadron_throws;
        effTree->SetBranchAddress("ND_OffAxis_Sim_mu_start_v_xyz_LAr", &xyz_pos);
        effTree->SetBranchAddress("ND_OffAxis_Sim_mu_start_p_xyz_LAr", &xyz_mom);
        effTree->SetBranchAddress("hadron_throw_result_LAr", &hadron_throw_result_LAr);
        TBranch *vtx_x=FDCAF->Branch("vtx_x", &(vtx_vals[0]), "vtx_vals/D"); //just put the whole corresponding
        TBranch *vtx_y=FDCAF->Branch("vtx_y", &(vtx_vals[1]), "vtx_vals/D"); //array in the "data type size"
        TBranch *vtx_z=FDCAF->Branch("vtx_z", &(vtx_vals[2]), "vtx_vals/D"); //section if you're using an array
        TBranch *lepmomx=FDCAF->Branch("LepMomX", &(lepmom_vals[0]), "lepmom_vals/D"); //element as the object
        TBranch *lepmomy=FDCAF->Branch("LepMomY", &(lepmom_vals[1]), "lepmom_vals/D");
        TBranch *lepmomz=FDCAF->Branch("LepMomZ", &(lepmom_vals[2]), "lepmom_vals/D");
        TBranch *hadron_throw_result=FDCAF->Branch("hadron_throw_result", &hadron_throws, "hadron_throws/L");
        TBranch *LAr_pos_index=FDCAF->Branch("i_LAr_pos", &i_LAr_pos, "i_LAr_pos/I");
        TBranch *vtx_pos_index=FDCAF->Branch("i_vtx_pos", &i_vtx_x, "i_vtx_x/I");
        Long64_t nentries=effTree->GetEntries();
        for (Long64_t i_event=0;i_event<nentries;i_event++)
        {
            effTree->GetEntry(i_event);
            for (i_LAr_pos=0;i_LAr_pos<15;i_LAr_pos++)
            {
                for (i_vtx_x=0;i_vtx_x<22;i_vtx_x++)
                {
                    for (int j=0;j<3;j++) {
                        vtx_vals[j]=(*xyz_pos)[i_LAr_pos][i_vtx_x][j];
                        lepmom_vals[j]=(*xyz_mom)[i_LAr_pos][i_vtx_x][j];
                    }
                    hadron_throws=(*hadron_throw_result_LAr)[i_LAr_pos][i_vtx_x];
                    // vtx_x->Fill();
                    // vtx_y->Fill();
                    // vtx_z->Fill();
                    // lepmomx->Fill();
                    // lepmomy->Fill();
                    // lepmomz->Fill();
                    // hadron_throw_result->Fill();
                    // LAr_pos_index->Fill();
                    // vtx_pos_index->Fill();
                    FDCAF->Fill();
                }
            }
        }
        write_file->Write();
        delete read_file;
        delete write_file;
    }
    return 0;
}