///////////////////////////////////////////////////////////////////
//
// Nate Dzbenski
// Prepare input data for Track Finder and Kalman Filter
// 31 Aug 2017
//
///////////////////////////////////////////////////////////////////



#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TH1F.h"
#include "TBenchmark.h"
#include "TRandom.h"
#include "TSystem.h"
#include <vector>
#include <map>


void PrepData_kf()
{

gROOT->SetStyle("Plain");

//______________________________________________________________________________________________
//  ___________________________________________Variables__________________________________________
//______________________________________________________________________________________________
    
    char filename[100];
    char outname[100];
    
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
    
    map<double, int> track_id;
    map<double, int> :: iterator itr;


//_____________________________________________________________________________________________
//  ___________________________________________ Openings ________________________________________
//_____________________________________________________________________________________________

    // Two trees must be read
    // One with generated vertex data
    // One with track data
    
    sprintf(filename,"shaping.root");
    sprintf(outname,"../KalmanFilter/KalRTPC/bonus_gemc_rec.root");
    
    TFile *myfile= TFile::Open(filename);
    if (!myfile) { return; }
    
    TTree *gen= (TTree*)myfile->Get("generated");
    TTree *rec= (TTree*)myfile->Get("RTPC");
    
    // Load generated data
    int evn_v = 0;
    std::vector<double>* z_v = 0;
    std::vector<double>* p_v = 0;
    std::vector<double>* pt_v = 0;
    std::vector<double>* th_v = 0;
    std::vector<double>* phi_v = 0;
    
    // Create branches for generated data
    TBranch *bevn_v = 0;
    gen->SetBranchAddress("evn_v",&evn_v,&bevn_v);
    TBranch *bz_v = 0;
    gen->SetBranchAddress("z_v",&z_v,&bz_v);
    TBranch *bp_v = 0;
    gen->SetBranchAddress("p_v",&p_v,&bp_v);
    TBranch *bpt_v = 0;
    gen->SetBranchAddress("pt_v",&pt_v,&bpt_v);
    TBranch *bth_v = 0;
    gen->SetBranchAddress("th_v",&th_v,&bth_v);
    TBranch *bphi_v = 0;
    gen->SetBranchAddress("phi_v",&phi_v,&bphi_v);
    
    // Load track data
    int evn = 0;
    std::vector<double> *CellID = 0;
    std::vector<double> *ADC = 0;
    std::vector<double> *Time = 0;
    std::vector<double> *time_shift = 0;
    std::vector<double> *vz_track = 0;
    std::vector<int> *tid = 0;
    
    TBranch *bevn = 0;
    rec->SetBranchAddress("evn",&evn,&bevn);
    TBranch *bCellID = 0;
    rec->SetBranchAddress("CellID",&CellID,&bCellID);
    TBranch *bADC = 0;
    rec->SetBranchAddress("ADC",&ADC,&bADC);
    TBranch *bTime = 0;
    rec->SetBranchAddress("Time",&Time,&bTime);
    TBranch *btime_shift = 0;
    rec->SetBranchAddress("TimeShift",&time_shift,&btime_shift);
    TBranch *btid = 0;
    rec->SetBranchAddress("TrackID",&tid,&btid);
    TBranch *bvz_track = 0;
    rec->SetBranchAddress("vz",&vz_track,&bvz_track);
    
    
    // Output root file for Jixie
    TFile *fout = TFile::Open(outname,"RECREATE");
    TTree *v_data = new TTree("Gen","Generated Vertex Data");
    TTree *rec_data = new TTree("Rec","Reconstructed Track Data");
    
    // Create variables for vertex data
    int Evn_v=0;
    std::vector<int>* Vtid = 0;
    std::vector<double>* Vz = 0;
    std::vector<double>* Vp = 0;
    std::vector<double>* Vpt = 0;
    std::vector<double>* Vth = 0;
    std::vector<double>* Vphi = 0;
    
    // Create branches of vertex data
    v_data->Branch("event_v", &Evn_v, "EventID/I");
    v_data->Branch("TrackID", &Vtid);
    v_data->Branch("z_v", &Vz);
    v_data->Branch("p_v", &Vp);
    v_data->Branch("pt_v", &Vpt);
    v_data->Branch("th_v", &Vth);
    v_data->Branch("phi_v", &Vphi);
    
    // Create branches of reconstructed data for new file
    int Evn_r;
    std::vector<int> *cell = new std::vector<int>;
    std::vector<int> *tid_r = new std::vector<int>;
    std::vector<double> *X_rec = new std::vector<double>;
    std::vector<double> *Y_rec = new std::vector<double>;
    std::vector<double> *Z_rec = new std::vector<double>;
    std::vector<double> *adc = new std::vector<double>;
    std::vector<double> *time_r = new std::vector<double>;
    std::vector<double> *shift_r = new std::vector<double>;
    std::vector<double> *z_track = new std::vector<double>;
    
    rec_data->Branch("Event", &Evn_r, "EventID/I");
    rec_data->Branch("ChanID", &cell);
    rec_data->Branch("TrackID", &tid_r);
    rec_data->Branch("X_rec", &X_rec);
    rec_data->Branch("Y_rec", &Y_rec);
    rec_data->Branch("Z_rec", &Z_rec);
    rec_data->Branch("ADC", &adc);
    rec_data->Branch("Time", &time_r);
    rec_data->Branch("TimeShift", &shift_r);
    rec_data->Branch("z_v", &z_track);
    
    
    //__________________________________________________________________________________
    //__________________________________ Readings ______________________________________
    //__________________________________________________________________________________
    
    int count = 0;
    for (Int_t eve = 0; eve < gen->GetEntries(); eve++){
        
        track_id.clear();
        
        cell->clear();
        X_rec->clear();
        Y_rec->clear();
        Z_rec->clear();
        adc->clear();
        time_r->clear();
        tid_r->clear();
        shift_r->clear();
        z_track->clear();
        
        Vz->clear();
        Vp->clear();
        Vpt->clear();
        Vth->clear();
        Vphi->clear();
        
        int chan=0;
        double t_s2pad = 0.;
        double dphi=0.;
        double phi_pad=0.;
        double z_pad=0.;
        double r_rec=0.;
        double phi_rec=0.;
        double dz=0.;
        int tid_i =0;
        
        Long64_t gentry = gen->LoadTree(eve);
        Long64_t rentry = rec->LoadTree(eve);
        
        // Load Generated Data
        bevn_v->GetEntry(gentry);
        bz_v->GetEntry(gentry);
        bp_v->GetEntry(gentry);
        bpt_v->GetEntry(gentry);
        bth_v->GetEntry(gentry);
        bphi_v->GetEntry(gentry);
        
        // Load RTPC data
        bevn->GetEntry(rentry);
        bCellID->GetEntry(rentry);
        bADC->GetEntry(rentry);
        bTime->GetEntry(rentry);
        btime_shift->GetEntry(rentry);
        btid->GetEntry(rentry);
        bvz_track->GetEntry(rentry);
        
        // Fills track data
        if(CellID->size()>8000){cout << "Too many hits, event skipped. Increase Max hit number of use vectors." << endl; continue;}
                          
        Evn_r = eve;
        
        for (int s = 0; s < CellID->size(); s++){
        
            if(Time->at(s)>7000.0) continue;
            
            if(s==0) track_id.insert(make_pair(vz_track->at(s), tid->at(s)));
            
            else if(tid->at(s) != tid->at(s-1)){
                track_id.insert(make_pair(vz_track->at(s), tid->at(s)));
                
                //cout << "Track ID: " << tid->at(s) << ", Z vertex: " << vz_track->at(s) << endl;
            }
            
            
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
            t_s2pad = Time->at(s)-t_gap; // in ns for MagBoltz
            
            r_rec=((-(TMath::Sqrt(a_t*a_t+(4.*b_t*t_s2pad)))+a_t+(14.*b_t))/(2.*b_t))*10.0; //in mm
            dphi=a_phi*(7.-r_rec/10.)+b_phi*(7.-r_rec/10.)*(7.-r_rec/10.); // in rad
            
            dz=0.; // in mm
            
            phi_rec=phi_pad-dphi-phi_gap;
            if( phi_rec<0.0 )  phi_rec+=2.0*PI;
            if( phi_rec>2.0*PI )  phi_rec-=2.0*PI;
            
            // Fill reconstructed data
            cell->push_back(CellID->at(s));
            X_rec->push_back(r_rec*(TMath::Cos(phi_rec)));
            Y_rec->push_back(r_rec*(TMath::Sin(phi_rec)));
            Z_rec->push_back(z_pad-dz);
            adc->push_back(ADC->at(s));
            time_r->push_back(Time->at(s));
            shift_r->push_back(time_shift->at(s));
            tid_r->push_back(tid->at(s));
            z_track->push_back(vz_track->at(s));
            
        }
        
        //cout << "\tVertex Z \tTrack ID" << endl;
        //for(int it = track_id.begin(); it != track_id.end(); it++){
        //    cout  <<  '\t' << itr->first
        //    <<  '\t' << itr->second << '\n';
        //}
        
        // Fill vertex data into variables
        for(int c=0; c < z_v->size(); c++){
            Evn_v = eve;
            Vz->push_back(z_v->at(c));
            Vp->push_back(p_v->at(c));
            Vpt->push_back(pt_v->at(c));
            Vth->push_back(th_v->at(c));
            Vphi->push_back(phi_v->at(c));
            
            
            tid_i = track_id.find(z_v->at(c))->second;
            
            
            Vtid->push_back(c+1);
            
            //cout << "For tid: " << Vtid->at(c) << ", vertex z: " << z_v->at(c) << endl;
        }
        
        // Fills the output trees
        v_data->Fill();
        rec_data->Fill();
        
    }
    
    
    //__________________________________________________________________________________________
    //  _______________________________________ Closing __________________________________________
    //__________________________________________________________________________________________
    
    fout->Write();
    delete fout;
    myfile->Close();
    
    // End all of the for statements for files
    
    
// End the program
}


