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
    
    TH1D *h_eLoss = new TH1D("h_eLoss","Energy loss",100,0,0);
    
    TH2F *h_Edep_v_phi = new TH2F("h_Edep_v_phi","Eloss vs. Phi",100,0,360,100,0,15);
    
    TH2F *h_Edep_v_theta = new TH2F("h_Edep_v_theta","Eloss vs. Theta",100,0,0,100,0,0);
    
    int n=0;
    
    
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
    
    std::vector<double>* px = 0;
    std::vector<double>* py = 0;
    std::vector<double>* pz = 0;
    
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
    TBranch *bpx = 0;
    gen_tree->SetBranchAddress("px",&px,&bpx);
    
    TBranch *bpy = 0;
    gen_tree->SetBranchAddress("py",&py,&bpy);
    
    TBranch *bpz = 0;
    gen_tree->SetBranchAddress("pz",&pz,&bpz);
    
    
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
        
        // Get generated mom info
        bpx->GetEntry(gen_tentry);
        bpy->GetEntry(gen_tentry);
        bpz->GetEntry(gen_tentry);
        
        double totEloss=0;
        
        double p_tot = 0;
        double phi_track = 0;
        double theta_track = 0;
        
        // ________________________________________ Reconstruction _______________________________________
        
        // find postition from Cell ID
        for (UInt_t s = 0; s < CellID->size(); s++) {
            
            // generated r position
            double r_pos=TMath::Sqrt(((Xpos->at(s))*(Xpos->at(s)))+((Ypos->at(s))*(Ypos->at(s))));
            
            if (r_pos<80) continue;
            
            else totEloss += edep->at(s);
            
        }// end finding position from Cell ID
        
        
        
        p_tot = TMath::Sqrt((px->at(0)*px->at(0)) + (py->at(0)*py->at(0)));
        phi_track = TMath::ATan2((py->at(0)),(px->at(0)))*(180/PI)+180.0;
        theta_track = TMath::ATan2(p_tot,pz->at(0))*(180/PI);
        
        if(totEloss != 0.0) h_eLoss->Fill(totEloss);
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
    h_eLoss->GetXaxis()->SetTitle("Edep/track [MeV]");
    //h_r->GetYaxis()->SetTitle("Counts");
    c1->SaveAs("figs/Edep.png");
    // -------------------------------------------------------------------
    
    // --------------------------- Eloss vs Phi ---------------------------
    c2->cd();
    h_Edep_v_phi->Draw("CONT4");
    h_Edep_v_phi->SetTitle("Edep v #phi");
    //h_Edep_v_phi->SetLineColor(kBlue);
    h_Edep_v_phi->GetXaxis()->SetTitle("#phi [deg]");
    h_Edep_v_phi->GetYaxis()->SetTitle("Edep/track [MeV]");
    c2->SaveAs("figs/Edep_v_phi.png");
    // -------------------------------------------------------------------
    
    // --------------------------- Eloss vs Theta ---------------------------
    c3->cd();
    h_Edep_v_theta->Draw("CONT");
    h_Edep_v_theta->SetTitle("Edep v #theta");
    //h_Edep_v_theta->SetLineColor(kBlue);
    h_Edep_v_theta->GetXaxis()->SetTitle("#theta [deg]");
    h_Edep_v_theta->GetYaxis()->SetTitle("Edep/track [MeV]");
    c3->SaveAs("figs/Edep_v_theta.png");
    // -------------------------------------------------------------------
    
    gPad->Update();
    
    
    myfile->Close();
    return;
    
    
    
    
    
}
