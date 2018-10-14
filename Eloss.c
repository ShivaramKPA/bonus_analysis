///////////////////////////////////////////////////////////////////
//
// Nate Dzbenski
// Open data from GEMC .root file, and compare simulation track
// to recustructed track.
//
// To run, type:
// root rec.c <filename> <number of events>
//
///////////////////////////////////////////////////////////////////


#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2F.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TChain.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <locale>
#include <sstream>
#include <math.h>


void Eloss()
{
    gROOT->SetStyle("Plain");
    
    //________________________________________________________________________________________________
    // ____________________________________________Variables__________________________________________
    //________________________________________________________________________________________________
    
    char filename[100];
    //const int NEve = argv[1];
    const int NEve = 100000;
    
    const double PI = TMath::ACos(-1);
    
    
    //________________________________________________________________________________________________
    //_____________________________________ Canvas and histograms ____________________________________
    //________________________________________________________________________________________________
    
    TCanvas *c1 = new TCanvas("c1","energy loss",800,600);
    
    TCanvas *c2 = new TCanvas("c2","energy loss v phi",800,600);
    
    TCanvas *c3 = new TCanvas("c3","energy loss v theta",800,600);
    
    TCanvas *c4 = new TCanvas("c4","delta theta",800,600);
    
    TH1D *h_eLoss = new TH1D("h_eLoss","Energy loss",100,0,0);
    
    TH1D *h_dTheta = new TH1D("h_dTheta","#it{#Delta#theta}",100,0,0);
    
    TH2F *h_Edep_v_phi = new TH2F("h_Edep_v_phi","Eloss vs. #it{#phi}",100,0,360,100,0,15);
    
    TH2F *h_Edep_v_theta = new TH2F("h_Edep_v_theta","Eloss vs. #it{#theta}",30,10,27,30,0,10);
    
    //________________________________________________________________________________________________
    // ___________________________________________ Openings __________________________________________
    //________________________________________________________________________________________________
    
    // Copy file name
    sprintf(filename,"/Users/nateMac/jlab_software/2.2/gemc/devel/data/gemc_Eloss.root");
    //TFile *myfile= new TFile(filename, "open");
    
    //sprintf(filename,argv[0]);
    
    TFile *myfile= new TFile(filename, "open");
    TTree *tree= (TTree*)myfile->Get("rtpc");
    TTree *gen_tree= (TTree*)myfile->Get("generated");
    
    std::vector<double>* CellID = 0;
    std::vector<double>* edep = 0;
    std::vector<double>* Xpos = 0;
    std::vector<double>* Ypos = 0;
    std::vector<double>* Zpos = 0;
    
    std::vector<double>* px_in = 0;
    std::vector<double>* py_in = 0;
    std::vector<double>* pz_in = 0;
    
    std::vector<double>* px_out = 0;
    std::vector<double>* py_out = 0;
    std::vector<double>* pz_out = 0;
    
    TBranch *bcellID =  0;
    tree->SetBranchAddress("CellID",&CellID,&bcellID);

    TBranch *bedep = 0;
    tree->SetBranchAddress("totEdep",&edep,&bedep);
    
    TBranch *bXpos = 0;
    tree->SetBranchAddress("avg_x",&Xpos,&bXpos);
    
    TBranch *bYpos = 0;
    tree->SetBranchAddress("avg_y",&Ypos,&bYpos);
    
    TBranch *bZpos = 0;
    tree->SetBranchAddress("avg_z",&Zpos,&bZpos);
    
    // generated info
    TBranch *bpx_in = 0;
    gen_tree->SetBranchAddress("px",&px_in,&bpx_in);
    
    TBranch *bpy_in = 0;
    gen_tree->SetBranchAddress("py",&py_in,&bpy_in);
    
    TBranch *bpz_in = 0;
    gen_tree->SetBranchAddress("pz",&pz_in,&bpz_in);
    
    // mom out info
    TBranch *bpx_out = 0;
    tree->SetBranchAddress("px",&px_out,&bpx_out);
    
    TBranch *bpy_out = 0;
    tree->SetBranchAddress("py",&py_out,&bpy_out);
    
    TBranch *bpz_out = 0;
    tree->SetBranchAddress("pz",&pz_out,&bpz_out);
    
    
    //________________________________________________________________________________________________
    // ___________________________________________ Readings __________________________________________
    //________________________________________________________________________________________________
    
    
    
    for(int i=0; i<tree->GetEntries(); i++){
        
        Long64_t tentry = tree->LoadTree(i);
        Long64_t gen_tentry = gen_tree->LoadTree(i);
        
        bcellID->GetEntry(tentry);
        bedep->GetEntry(tentry);
        
        // Get x,y,z positions
        bXpos->GetEntry(tentry);
        bYpos->GetEntry(tentry);
        bZpos->GetEntry(tentry);
        
        bpx_out->GetEntry(tentry);
        bpy_out->GetEntry(tentry);
        bpz_out->GetEntry(tentry);
        
        // Get generated mom info
        bpx_in->GetEntry(gen_tentry);
        bpy_in->GetEntry(gen_tentry);
        bpz_in->GetEntry(gen_tentry);
        
        double totEloss=0;
        
        double p_tot_in = 0;
        double p_tot_out = 0;
        
        double phi_track = 0;
        double theta_track = 0;
        
        double theta_out = 0;
        
        // ________________________________________ Reconstruction _______________________________________
        
        // find postition from Cell ID
        for (UInt_t s = 0; s < Xpos->size(); s++) {
            
            // generated r position
            double r_pos=TMath::Sqrt(((Xpos->at(s))*(Xpos->at(s)))+((Ypos->at(s))*(Ypos->at(s))));
            
            if (r_pos<80) continue;
            
            else totEloss += edep->at(s);
            
        }// end finding position from Cell ID
        
        
        
        p_tot_in = TMath::Sqrt((px_in->at(0)*px_in->at(0)) + (py_in->at(0)*py_in->at(0)) + (pz_in->at(0)*pz_in->at(0)));
        
        p_tot_out = TMath::Sqrt((px_out->back()*px_out->back()) + (py_out->back()*py_out->back()) + (pz_out->back()*pz_out->back()));
        
        phi_track = TMath::ATan2((py_in->at(0)),(px_in->at(0)))*(180/PI)+180.0;
        theta_track = TMath::ACos(pz_in->at(0)/p_tot_in)*(180/PI);

        theta_out = TMath::ACos(pz_in->back()/p_tot_out)*(180/PI);
        
        // if(totEloss != 0.0)
        h_eLoss->Fill(totEloss);
        h_dTheta->Fill(theta_track-theta_out);
        h_Edep_v_phi->Fill(phi_track, totEloss);
        h_Edep_v_theta->Fill(theta_track, totEloss);
        
    } // tree->GetEntries
    
    
    
    //________________________________________________________________________________________________
    // _____________________________________________ Displays ________________________________________
    //________________________________________________________________________________________________
    
    
    // --------------------------- Eloss histogram ---------------------------
    c1->cd();
    h_eLoss->Draw();
    //h_eLoss->SetLineColor(kBlue);
    h_eLoss->GetXaxis()->SetTitle("#it{E}_{dep}/track [MeV]");
    //h_r->GetYaxis()->SetTitle("Counts");
    c1->SaveAs("figs/Edep.png");
    // -------------------------------------------------------------------
    
    // --------------------------- Eloss vs Phi ---------------------------
    c2->cd();
    h_Edep_v_phi->Draw("CONT4");
    h_Edep_v_phi->SetTitle("Edep v #it{#phi}");
    //h_Edep_v_phi->SetLineColor(kBlue);
    h_Edep_v_phi->GetXaxis()->SetTitle("#it{#phi} [deg]");
    h_Edep_v_phi->GetYaxis()->SetTitle("Edep/track [MeV]");
    c2->SaveAs("figs/Edep_v_phi.png");
    // -------------------------------------------------------------------
    
    // --------------------------- Eloss vs Theta ---------------------------
    c3->cd();
    h_Edep_v_theta->Draw("CONTZ");
    h_Edep_v_theta->SetTitle("Edep v #it{#theta}");
    //h_Edep_v_theta->SetLineColor(kBlue);
    h_Edep_v_theta->GetXaxis()->SetTitle("#it{#theta} [deg]");
    h_Edep_v_theta->GetYaxis()->SetTitle("Edep/track [MeV]");
    c3->SaveAs("figs/Edep_v_theta.png");
    // -------------------------------------------------------------------
    
    // --------------------------- Change in theta histogram ---------------------------
    c4->cd();
    h_dTheta->Draw();
    h_dTheta->GetXaxis()->SetTitle("#it{#theta}_{in} - #it{#theta}_{out} [deg]");
    //h_r->GetYaxis()->SetTitle("Counts");
    c4->SaveAs("figs/dTheta.png");
    // -------------------------------------------------------------------
    
    
    gPad->Update();
    
    
    myfile->Close();
    return;
    
    
    
    
    
}
