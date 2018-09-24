///////////////////////////////////////////////////////////////////
//
// Nate Dzbenski
// Open data from GEMC .root file, and compare simulation track
// to recustructed track.
//
///////////////////////////////////////////////////////////////////


#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
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

void bonus_energy()
{
gROOT->SetStyle("Plain");

//________________________________________________________________________________________________
// ____________________________________________Variables__________________________________________
//________________________________________________________________________________________________
    
   char filename[100];
   const int NEve = 1000;
    
    double PAD_W = 2.79; // in mm
    double PAD_S = 80.0; //in mm
    double PAD_L = 4.0; // in mm
    double RTPC_L=400.0; // in mm
    
    int Num_of_Rows = (2.0*(TMath::Pi())*PAD_S)/PAD_W;
    int Num_of_Cols = RTPC_L/PAD_L;
    int TotChan = Num_of_Rows*Num_of_Cols;
    
    double PI=TMath::Pi();
    
    double z0 = -(RTPC_L/2.0); // front of RTPC in mm at the center of the pad

    double phi_per_pad = PAD_W/PAD_S; // in rad
    double Edep =0;

    // MagBoltz parameters
    float a=-0.177, b=0.0392,c=10.977; // for f(x)=time(radius)
    float d=0.027,e=-0.517,f=2.4615;   // for f(x)=dphi(radius)

    TRandom *gRandom = new TRandom();
    

//________________________________________________________________________________________________
//_____________________________________ Canvas and histograms ____________________________________
//________________________________________________________________________________________________

 TCanvas *c1 = new TCanvas("c1","Energy Deposition",800,600);
 TCanvas *c2 = new TCanvas("c2","Delta_r",800,600);
 //TH2D *h_Energy = new TH2D("h_Energy","E/mm",100,0,80,100,0,0);          // 2D s plot
  TH1D *h_Energy = new TH1D("h_Energy","E_{dep}/mm for He:DME at 90:10 with 2.5 kV",100,0,6000);
  TH1D *h_deltar = new TH1D("h_deltar","Step Size (#Deltar) [mm]",300,0,0);


//________________________________________________________________________________________________
// ___________________________________________ Openings __________________________________________
//________________________________________________________________________________________________

 // Copy file name
            sprintf(filename,"bonus_gemc_co2_out_1000.root");

            TFile *myfile= new TFile(filename, "open");
            TTree *tree= (TTree*)myfile->Get("xbonus");

            std::vector<double>* CellID = 0;
            std::vector<double>* ADC = 0;
            std::vector<double>* TDC = 0;
            std::vector<double>* step = 0;
            std::vector<double>* EDep = 0;
            std::vector<double>* Xpos = 0;
            std::vector<double>* Ypos = 0;
            std::vector<double>* Zpos = 0;
            std::vector<double>* hitn = 0;
            std::vector<double>* phi_rad = 0;
    
            TBranch *bcellID =  0;
            tree->SetBranchAddress("CellID",&CellID,&bcellID);
    
            TBranch *badc = 0;
            tree->SetBranchAddress("ADC",&ADC,&badc);
    
            TBranch *btdc = 0;
            tree->SetBranchAddress("TDC",&TDC,&btdc);
    
            TBranch *bstep = 0;
            tree->SetBranchAddress("step",&step,&bstep);
    
            TBranch *bEdep = 0;
            tree->SetBranchAddress("EDep",&EDep,&bEdep);
    
            TBranch *bXpos = 0;
            tree->SetBranchAddress("posX",&Xpos,&bXpos);
    
            TBranch *bYpos = 0;
            tree->SetBranchAddress("posY",&Ypos,&bYpos);
    
            TBranch *bZpos = 0;
            tree->SetBranchAddress("posZ",&Zpos,&bZpos);
    
            TBranch *bhitn = 0;
            tree->SetBranchAddress("hitn",&hitn,&bhitn);

            TBranch *bphi_rad = 0;
            tree->SetBranchAddress("phiRad",&phi_rad,&bphi_rad);
    


//________________________________________________________________________________________________
// ___________________________________________ Readings __________________________________________
//________________________________________________________________________________________________


    
   for(int i=0; i<tree->GetEntries(); i++){

      Long64_t tentry = tree->LoadTree(i);
       
      bcellID->GetEntry(tentry);
      btdc->GetEntry(tentry);
      badc->GetEntry(tentry);
      bstep->GetEntry(tentry);
      bEdep->GetEntry(tentry);
      bphi_rad->GetEntry(tentry);
       
      // Get x,y,z positions
      bXpos->GetEntry(tentry);
      bYpos->GetEntry(tentry);
      bZpos->GetEntry(tentry);
      bhitn->GetEntry(tentry);
 
           double r_pos=0;
       
       // ________________________________________ Reconstruction _______________________________________
       
       // find postition from Cell ID
       for (UInt_t s = 0; s < CellID->size(); s++) {

           //h_Energy->Fill((EDep->at(s)));
           Edep = (EDep->at(s)*(1E6));
           
           double delta_r =0;
           
           if(s==0){
               h_Energy->Fill(Edep);
           }
           else{

           // generated position of ionization in s
        
               delta_r =(TMath::Sqrt(((Xpos->at(s))*(Xpos->at(s)))+((Ypos->at(s))*(Ypos->at(s)))+((Zpos->at(s))*(Zpos->at(s)))))-r_pos;
           
           //--------------------------------------- Plotting ---------------------------------------------
           
               h_Energy->Fill(Edep/delta_r);
               h_deltar->Fill(delta_r);

           
           // --------------------------------------------------------------------------------------------
           
           r_pos=TMath::Sqrt(((Xpos->at(s))*(Xpos->at(s)))+((Ypos->at(s))*(Ypos->at(s)))+((Zpos->at(s))*(Zpos->at(s))));
           }
           
       }// end finding position from Cell ID
            


 
   } // tree->GetEntries




//________________________________________________________________________________________________
// _____________________________________________ Displays ________________________________________
//________________________________________________________________________________________________


    // --------------------------- Energy Deposition/mm histogram ---------------------------
    c1->cd();
        h_Energy->Draw();
        h_Energy->SetLineColor(kBlue);
        // h_temp->Draw("same");
        // h_temp->SetLineColor(kRed);
        //h_Energy->SetStats(0);
        h_Energy->GetXaxis()->SetTitle("E_{dep}/mm [eV/mm]");
        h_Energy->GetYaxis()->SetTitle("Counts");
    // -------------------------------------------------------------------

    c2->cd();
    h_deltar->Draw();
    h_deltar->GetXaxis()->SetTitle("Step size [mm]");
    h_deltar->GetYaxis()->SetTitle("Counts");
    
    gPad->Update();


    myfile->Close();
return;
    
    


    
}
