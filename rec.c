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


void rec()
{
    gROOT->SetStyle("Plain");
    
    //________________________________________________________________________________________________
    // ____________________________________________Variables__________________________________________
    //________________________________________________________________________________________________
    
    char filename[100];
    //const int NEve = argv[1];
    const int NEve = 10000;
    
    double PAD_W = 2.79; // in mm
    double PAD_S = 80.0; //in mm
    double PAD_L = 4.0; // in mm
    double RTPC_L= 384.0; // in mm
    
    int Num_of_Rows = (2.0*(TMath::Pi())*PAD_S)/PAD_W;
    int Num_of_Cols = RTPC_L/PAD_L;
    int TotChan = Num_of_Rows*Num_of_Cols;
    
    double PI=TMath::Pi();
    
    double z0 = -(RTPC_L/2.0); // front of RTPC in mm at the center of the pad
    
    double phi_per_pad = PAD_W/PAD_S; // in rad
    
    // MagBoltz parameters
    float a_t=1741.179712, b_t=-1.25E+02; // for f(x)=time(radius)
    float a_phi=0.161689123, b_phi=0.023505021;   // for f(x)=dphi(radius)
    
    float t_2GEM2 = 296.082;
    float t_2GEM3 = 296.131;
    float t_2PAD = 399.09;
    float t_gap = t_2GEM2 + t_2GEM3 + t_2PAD;
    
    float phi_2GEM2 = 0.0492538;
    float phi_2GEM3 = 0.0470817;
    float phi_2PAD = 0.0612122;
    float phi_gap = phi_2GEM2 + phi_2GEM3 + phi_2PAD;
    
    TRandom *gRandom = new TRandom();
    
    
    //________________________________________________________________________________________________
    //_____________________________________ Canvas and histograms ____________________________________
    //________________________________________________________________________________________________
    
    TCanvas *c1 = new TCanvas("c1","delta s",800,600);
    TCanvas *c2 = new TCanvas("c2","delta z",800,600);
    // TCanvas *c4 = new TCanvas("c4","delta y",800,600);
    TCanvas *c5 = new TCanvas("c5","delta phi, old solenoid map",800,600);
    //TCanvas *c6 = new TCanvas("c6","Generated 100 MeV/c proton Track",800,600);
    
    TH1D *h_r = new TH1D("h_r","#Deltas, no cut",400,-3.,3);
    //TH2D *h_r = new TH2D("h_r","#Deltas",100,0,80,100,0,0);          // 2D s plot
    TH1D *h_phi = new TH1D("h_phi","#Delta#phi, old solenoid map",400,-0.08,0.08);
    //TH1D *h_phi_t = new TH1D("h_phi_t","#Delta#phi",3000,0,0);
    //TH2F *h_phi = new TH2F("h_phi","#Delta#phi",300,0,0,300,0,0);     // 2D phi plot
    TH1D *h_x = new TH1D("h_x","#Deltax",200,-3.,3.);
    TH1D *h_y = new TH1D("h_y","#Deltay",200,-3.,3.);
    TH1D *h_z = new TH1D("h_z","#Deltaz, old solenoid map",400,-4,6 );
    //    TH2D *h_z = new TH2D("h_z","#Deltaz",4000,-200,200,1000,-1,101);            // 2D z plot
    
    // array of graphs that contains TDC versus hit number
    TGraph *g_TDCevol[NEve];
    for(int g=0;g<NEve;g++){
        g_TDCevol[g] = new TGraph();
    }
    
    // array of 2D graphs that contains simulation tracks
    TGraph *g_tracks[NEve];
    for(int g=0;g<NEve;g++){
        g_tracks[g] = new TGraph();
    }
    
    // array of 2D graphs that contains reconstructed tracks
    TGraph *g_tracks_rec[NEve];
    for(int g=0;g<NEve;g++){
        g_tracks_rec[g] = new TGraph();
    }
    
    // array of 3D graphs that contains pad hits
    TGraph *g_tracks_pad[NEve];
    for(int g=0;g<NEve;g++){
        g_tracks_pad[g] = new TGraph();
    }
    
    
    
    TMultiGraph *mgTracks = new TMultiGraph();
    
    int n=0;
    
    
    //________________________________________________________________________________________________
    // ___________________________________________ Openings __________________________________________
    //________________________________________________________________________________________________
    
    // Copy file name
    sprintf(filename,"/Users/nateMac/jlab_software/2.2/gemc/devel/gemc_resolution_old.root");
    //TFile *myfile= new TFile(filename, "open");
    
    //sprintf(filename,argv[0]);
    
    TFile *myfile= new TFile(filename, "open");
    TTree *tree= (TTree*)myfile->Get("rtpc");
    
    std::vector<double>* CellID = 0;
    std::vector<double>* TDC = 0;
    std::vector<double>* Xpos = 0;
    std::vector<double>* Ypos = 0;
    std::vector<double>* Zpos = 0;
    
    TBranch *bcellID =  0;
    tree->SetBranchAddress("CellID",&CellID,&bcellID);
    
    TBranch *btdc = 0;
    tree->SetBranchAddress("Time",&TDC,&btdc);
    
    TBranch *bXpos = 0;
    tree->SetBranchAddress("avg_x",&Xpos,&bXpos);
    
    TBranch *bYpos = 0;
    tree->SetBranchAddress("avg_y",&Ypos,&bYpos);
    
    TBranch *bZpos = 0;
    tree->SetBranchAddress("avg_z",&Zpos,&bZpos);
    
    
    //________________________________________________________________________________________________
    // ___________________________________________ Readings __________________________________________
    //________________________________________________________________________________________________
    
    
    
    for(int i=0; i<tree->GetEntries(); i++){
        
        Long64_t tentry = tree->LoadTree(i);
        
        bcellID->GetEntry(tentry);
        btdc->GetEntry(tentry);
        
        // Get x,y,z positions
        bXpos->GetEntry(tentry);
        bYpos->GetEntry(tentry);
        bZpos->GetEntry(tentry);
        
        double s_max=0;
        
        
        // ________________________________________ Reconstruction _______________________________________
        
        // find postition from Cell ID
        for (UInt_t s = 0; s < CellID->size(); s++) {
            
            int chan=0;
            double t_s2pad = 0;
            double dphi=0;
            
            double x_pad=0;
            double x_rec=0;
            double delta_x=0;
            
            double y_pad=0;
            double y_rec=0;
            double delta_y=0;
            
            double z_pad=0;
            double z_rec=0;
            double z_hit=0;  // position of the hit on a single pad in z
            double delta_z=0;
            double dz=0.;
            
            double r_pos=0;
            double r_rec=0;
            double r_temp=0;
            double delta_r=0;
            
            double phi_pos=0;
            double phi_rt=0;
            double phi_hit=0; // position of the hit on a single pad in phi
            double phi_pad=0;
            double phi_rec=0;
            double phi_rad_temp=0;
            double delta_phi=0;
            
            // generated position of ionization in phi
            phi_pos=(TMath::ATan2(Ypos->at(s),Xpos->at(s)));
            
            // generated position of ionization in s
            r_pos=TMath::Sqrt(((Xpos->at(s))*(Xpos->at(s)))+((Ypos->at(s))*(Ypos->at(s))));
            
            
            // ------------------ find z and phi of pad from CellID ------------------
            chan = CellID->at(s);
            
            int col = chan%Num_of_Cols;
            int row=(chan-col)/Num_of_Cols;
            float z_shift = row%4;
            
            phi_pad=(row*phi_per_pad)+(phi_per_pad/2.0);
            if(phi_pad>= 2.0*PI) phi_pad -= 2.0*PI;
            if(phi_pad<0.0) phi_pad += 2.0*PI;
            
            z_pad=z0+(col*PAD_L)+(PAD_L/2.0)+z_shift;
            // -----------------------------------------------------------------------
            
            
            // find reconstructed position of ionization from TDC info
            t_s2pad = TDC->at(s)-t_gap; // in ns for MagBoltz
            
            r_rec=((-(TMath::Sqrt(a_t*a_t+(4.*b_t*t_s2pad)))+a_t+(14.*b_t))/(2.*b_t))*10.0; //in mm
            dphi=a_phi*(7.-r_rec/10.)+b_phi*(7.-r_rec/10.)*(7.-r_rec/10.); // in rad
            
            dz=0.; // in mm
            
            phi_rec=phi_pad-dphi-phi_gap;
            if( phi_rec<0.0 )  phi_rec+=2.0*PI;
            if( phi_rec>2.0*PI )  phi_rec-=2.0*PI;
            
            // x,y,z pos of reconstructed track
            x_rec=r_rec*(TMath::Cos(phi_rec));
            y_rec=r_rec*(TMath::Sin(phi_rec));
            z_rec=z_pad-dz;
            
            // x,y,z pos of pad hit
            x_pad=(PAD_S)*(TMath::Cos(phi_pad));
            y_pad=(PAD_S)*(TMath::Sin(phi_pad));
            
            // actual position on pad of hits
            z_hit=Zpos->at(s)-z0-(col*PAD_L)-z_shift;
            
            // find differences (delta = generated-reconstructed)
            delta_x=Xpos->at(s)-x_rec;
            delta_y=Ypos->at(s)-y_rec;
            delta_z=Zpos->at(s)-z_rec;
            delta_r=r_pos-r_rec;
            delta_phi = phi_pos-phi_rec;
            
            
            //--------------------------------------- Plotting ---------------------------------------------
            
            
            //if (r_pos>69)
            h_r->Fill(delta_r);
            //if (r_pos>69)
            h_phi->Fill(delta_phi);
            //h_phi_t->Fill(delta_phi);
            //h_phi->Fill(phi_hit,z_hit);       // set up to display hits on pad
            //h_x->Fill(gaus_check);
            h_y->Fill(delta_y);
            
            h_z->Fill(delta_z);
            //h_z->Fill(delta_z,Zpos->at(s));
            
            
            g_tracks[i]->SetPoint(s,Xpos->at(s),Ypos->at(s));
            g_tracks_pad[i]->SetPoint(s,x_pad,y_pad);
            g_tracks_rec[i]->SetPoint(s,x_rec,y_rec);
            
            // --------------------------------------------------------------------------------------------
            
            
            
        }// end finding position from Cell ID
        
        
        //mgTDCevol->Add(g_TDCevol[i]);
        //g_TDCevol[i]->SetLineColor(i%5);
        
        mgTracks->Add(g_tracks[i]);
        mgTracks->Add(g_tracks_rec[i]);
        mgTracks->Add(g_tracks_pad[i]);
        
        g_tracks[i]->SetMarkerSize(0.5);
        g_tracks[i]->SetMarkerStyle(kFullCircle);
        g_tracks[i]->SetMarkerColorAlpha(kBlue, 0.4);
        
        g_tracks_rec[i]->SetMarkerSize(0.5);
        g_tracks_rec[i]->SetMarkerStyle(kFullCircle);
        g_tracks_rec[i]->SetMarkerColorAlpha(kRed, 0.4);
        
        g_tracks_pad[i]->SetMarkerSize(1.0);
        g_tracks_pad[i]->SetMarkerStyle(kFullCircle);
        g_tracks_pad[i]->SetMarkerColorAlpha(kBlack, 0.4);
        
    } // tree->GetEntries
    
    cout << n << endl;
    
    // mgTracks->Add(g_tracks[0]);
    // mgTracks->Add(g_tracks_rec[0]);
    //  mgTracks->Add(g_tracks_pad[0]);
    
    
    //________________________________________________________________________________________________
    // _____________________________________________ Displays ________________________________________
    //________________________________________________________________________________________________
    
    
    // --------------------------- s histogram ---------------------------
    c1->cd();
    h_r->Draw();
    h_r->SetLineColor(kBlue);
    //h_r->Fit("gaus");
    //h_r->SetStats(0);
    h_r->GetXaxis()->SetTitle("#Deltas [mm]");
    //h_r->GetYaxis()->SetTitle("Counts");
    // -------------------------------------------------------------------
    
    // --------------------------- z histogram ---------------------------
    c2->cd();
    // h_z->Fit("gaus");
    h_z->SetLineColor(kBlue);
    h_z->Draw("");
    h_z->GetXaxis()->SetTitle("#Deltaz [mm]");
    //h_z->GetYaxis()->SetTitle("Counts");
    // h_z->GetXaxis()->SetTitle("#Deltaz [mm]");
    // -------------------------------------------------------------------
    
    
    // --------------------------- x histogram ---------------------------
    
    // -------------------------------------------------------------------
    /*
     // --------------------------- y histogram ---------------------------
     c4->cd();
     //  h_y->Fit("gaus");
     h_y->SetLineColor(kBlue);
     h_y->Draw();
     h_y->GetXaxis()->SetTitle("#Deltay [mm]");
     // -------------------------------------------------------------------
     */
    
    // --------------------------- phi histogram -------------------------
    c5->cd();
    h_phi->Draw("cont");
    //h_phi->Fit("gaus");
    h_phi->SetLineColor(kBlue);
    h_phi->GetXaxis()->SetTitle("#Delta#phi [rad]");
    //h_phi->GetYaxis()->SetTitle("Count"); //z_{pad} mm");
    // -------------------------------------------------------------------
    
    
    
    // --------------------------- Track Plots -------------------------
    /*c6->cd();
    
    mgTracks->SetTitle("Gen Track [Blue], Rec Track [Red], Pad Hit [Black] at #theta=90#circ#pm30#circ; x [mm]; y [mm]");
    
    
    TEllipse *ground = new TEllipse(0,0,20,20);
    ground->SetFillColorAlpha(kCyan, 0.3);
    ground->Draw();
    
    TEllipse *cathode = new TEllipse(0,0,30,30);
    cathode->SetFillColorAlpha(kCyan, 0.3);
    cathode->Draw();
    
    TEllipse *target = new TEllipse(0,0,3,3);
    target->SetFillColorAlpha(kOrange, 1);
    target->Draw();
    
    TEllipse *gem1 = new TEllipse(0,0,70,70);
    gem1->SetFillColorAlpha(kWhite, 0.01);
    gem1->SetLineColor(8);
    gem1->Draw();
    
    TEllipse *gem2 = new TEllipse(0,0,73,73);
    gem2->SetFillColorAlpha(kWhite, 0.01);
    gem2->SetLineColor(8);
    gem2->Draw();
    
    TEllipse *gem3 = new TEllipse(0,0,76,76);
    gem3->SetFillColorAlpha(kWhite, 0.01);
    gem3->SetLineColor(8);
    gem3->Draw();
    
    TEllipse *pads = new TEllipse(0,0,80,80);
    pads->SetFillColorAlpha(kWhite, 0.01);
    pads->SetLineColor(4);
    pads->Draw();
    
    mgTracks->Draw("ap"); */
    // -------------------------------------------------------------------
    gPad->Update();
    
    
    myfile->Close();
    return;
    
    
    
    
    
}
