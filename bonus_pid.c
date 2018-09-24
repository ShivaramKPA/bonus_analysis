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

void bonus_pid()
{
gROOT->SetStyle("Plain");

//________________________________________________________________________________________________
// ____________________________________________Variables__________________________________________
//________________________________________________________________________________________________
    
   char filename[100];
   const int NEve = 10000;
    
   double Edep = 0.;
   double r_old = 0.;
   double r_new = 0.;
   double dr = 0.;
   double p = 0.;

    TRandom *gRandom = new TRandom();
    

//________________________________________________________________________________________________
//_____________________________________ Canvas and histograms ____________________________________
//________________________________________________________________________________________________

 TCanvas *c1 = new TCanvas("c1","Energy Deposition",800,600);
 TH2D *h_Energy = new TH2D("h_Energy","",100,1.,1000.,100,0.,5.);          // 2D s plot
 //TH1D *h_Energy = new TH1D("h_Energy","",100,0,0);



//________________________________________________________________________________________________
// ___________________________________________ Openings __________________________________________
//________________________________________________________________________________________________

 // Copy file name
    sprintf(filename,"bonus_gemc_out_proton.root");

    TFile *myfile= new TFile(filename, "open");
    TTree *tree= (TTree*)myfile->Get("RTPC");

    
    std::vector<double>* EDep = 0;
    std::vector<double>* Xpos = 0;
    std::vector<double>* Ypos = 0;
    std::vector<double>* Zpos = 0;
    std::vector<double>* Px = 0;
    std::vector<double>* Py = 0;
    std::vector<double>* Pz = 0;
    
    TBranch *bEdep = 0;
    tree->SetBranchAddress("totEdep",&EDep,&bEdep);
    
    TBranch *bXpos = 0;
    tree->SetBranchAddress("avg_x",&Xpos,&bXpos);
    
    TBranch *bYpos = 0;
    tree->SetBranchAddress("avg_y",&Ypos,&bYpos);
    
    TBranch *bZpos = 0;
    tree->SetBranchAddress("avg_z",&Zpos,&bZpos);
    
    TBranch *bPx = 0;
    tree->SetBranchAddress("px",&Px,&bPx);

    TBranch *bPy = 0;
    tree->SetBranchAddress("py",&Py,&bPy);
    
    TBranch *bPz = 0;
    tree->SetBranchAddress("pz",&Pz,&bPz);

//________________________________________________________________________________________________
// ___________________________________________ Readings __________________________________________
//________________________________________________________________________________________________


    
   for(int i=0; i<tree->GetEntries(); i++){

      Long64_t tentry = tree->LoadTree(i);
       
      bEdep->GetEntry(tentry);
      bXpos->GetEntry(tentry);
      bYpos->GetEntry(tentry);
      bZpos->GetEntry(tentry);
      bPx->GetEntry(tentry);
      bPy->GetEntry(tentry);
      bPz->GetEntry(tentry);
       
       
       // find postition from Cell ID
       for (UInt_t s = 0; s < EDep->size(); s++) {

           //h_Energy->Fill((EDep->at(s)));
           Edep = (EDep->at(s)*(1E3));

           // generated position of ionization in s
        
           r_new = (TMath::Sqrt(((Xpos->at(s))*(Xpos->at(s)))+((Ypos->at(s))*(Ypos->at(s)))+((Zpos->at(s))*(Zpos->at(s)))));

	   dr = TMath::Abs(r_new - r_old);
           
           p = (TMath::Sqrt(((Px->at(s))*(Px->at(s)))+((Py->at(s))*(Py->at(s)))+((Pz->at(s))*(Pz->at(s)))));
           
           //--------------------------------------- Plotting ---------------------------------------------
           
               h_Energy->Fill(p,Edep/dr); //,Edep/dr);

		r_old = r_new;

           
           // --------------------------------------------------------------------------------------------
           
           
           
       }// end finding position from Cell ID
            


 
   } // tree->GetEntries




//________________________________________________________________________________________________
// _____________________________________________ Displays ________________________________________
//________________________________________________________________________________________________


    // --------------------------- Energy Deposition/mm histogram ---------------------------
    c1->cd();
        h_Energy->Draw("COL");
        //h_Energy->SetLineColor(kBlue);
        //h_Energy->SetStats(0);
        h_Energy->GetXaxis()->SetTitle("p [MeV/c]");
        h_Energy->GetYaxis()->SetTitle("#frac{dE}{dx} [MeV/mm]");
    // -------------------------------------------------------------------

    gPad->Update();


    myfile->Close();
return;
    
    


    
}
