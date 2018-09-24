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
    
    
    
    
    //________________________________________________________________________________________________
    //_____________________________________ Canvas and histograms ____________________________________
    //________________________________________________________________________________________________
    
    TCanvas *c1 = new TCanvas("c1","energy loss",800,600);
    
    TCanvas *c2 = new TCanvas("c2","energy loss v phi",800,600);
    
    TH1D *h_eLoss = new TH1D("h_eLoss","Energy loss",100,0,0);
    
    TH2F *h_EvPhi = new TH2F("h_EvPhi","Eloss vs. Phi",100,0,15,100,0,365);
    
    int n=0;
    
    
    //________________________________________________________________________________________________
    // ___________________________________________ Openings __________________________________________
    //________________________________________________________________________________________________
    
    // Copy file name
    sprintf(filename,"/Users/nateMac/jlab_software/2.2/gemc/devel/gemc_Eloss.root");
    //TFile *myfile= new TFile(filename, "open");
    
    //sprintf(filename,argv[0]);
    
    TFile *myfile= new TFile(filename, "open");
    TTree *tree= (TTree*)myfile->Get("rtpc");
    
    std::vector<double>* CellID = 0;
    std::vector<double>* edep = 0;
    std::vector<double>* Xpos = 0;
    std::vector<double>* Ypos = 0;
    std::vector<double>* Zpos = 0;
    
    TBranch *bcellID =  0;
    tree->SetBranchAddress("CellID",&CellID,&bcellID);
    
    // change after you get the name
    TBranch *bedep = 0;
    tree->SetBranchAddress("totEdep",&edep,&bedep);
    
    TBranch *bXpos = 0;
    tree->SetBranchAddress("avg_x",&Xpos,&bXpos);
    
    TBranch *bYpos = 0;
    tree->SetBranchAddress("avg_y",&Ypos,&bYpos);
    
    TBranch *bZpos = 0;
    tree->SetBranchAddress("avg_z",&Zpos,&bZpos);
    
    //TBranch *bpy = 0;
    //tree->SetBranchAddress("py",&py,&bpy);
    
    
    //________________________________________________________________________________________________
    // ___________________________________________ Readings __________________________________________
    //________________________________________________________________________________________________
    
    
    
    for(int i=0; i<tree->GetEntries(); i++){
        
        Long64_t tentry = tree->LoadTree(i);
        
        bcellID->GetEntry(tentry);
        bedep->GetEntry(tentry);
        
        // Get x,y,z positions
        bXpos->GetEntry(tentry);
        bYpos->GetEntry(tentry);
        bZpos->GetEntry(tentry);
        
        double totEloss=0;
        double phi_track=0;
        
        
        // ________________________________________ Reconstruction _______________________________________
        
        // find postition from Cell ID
        for (UInt_t s = 0; s < CellID->size(); s++) {
            
            
            // generated r position
            double r_pos=TMath::Sqrt(((Xpos->at(s))*(Xpos->at(s)))+((Ypos->at(s))*(Ypos->at(s))));
            
            if (r_pos<80) continue;
            
            else totEloss += edep->at(s);
            
            if (s == CellID->size()-1) {
                phi_track = TMath::ATan2((Ypos->at(s)),(Xpos->at(s)))*(180/acos(-1))+180;
                h_EvPhi->Fill(totEloss, phi_track);
            }
            
            
        }// end finding position from Cell ID
        
        if(totEloss != 0.0) h_eLoss->Fill(totEloss);

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
    
    // --------------------------- Eloss histogram ---------------------------
    c2->cd();
    h_EvPhi->Draw("CONT");
    h_EvPhi->SetTitle("#phi v Edep");
    //h_EvPhi->SetLineColor(kBlue);
    h_EvPhi->GetXaxis()->SetTitle("Edep/track [MeV]");
    h_EvPhi->GetYaxis()->SetTitle("#phi [deg]");
    c2->SaveAs("figs/phi_v_Edep.png");
    // -------------------------------------------------------------------
    
    gPad->Update();
    
    
    myfile->Close();
    return;
    
    
    
    
    
}
