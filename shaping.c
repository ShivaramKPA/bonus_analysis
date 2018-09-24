///////////////////////////////////////////////////////////////////
//
// Gabriel Charles
// Shape then fit the signal from BONuS detector
// 10/05/2017
//
///////////////////////////////////////////////////////////////////

//Reminder to lunch Gemc on JLab farms from my own Gemc version:
//~/ODU/gemc/source.git/trunk/gemc ../gcards/bonus.gcard -USE_GUI=0 -N=10
//evio2root -B=rtpc -INPUTF=file.ev
//Note that the gcard contains the path to my own version of the BoNuS12 detector

//Data should be from Nathan Dzbenski

//This file is creating signal shape response from the time and energy obtained in Gemc.
//For the moment the shape is a gaussian, which is probably wrong.
//This part should be soon included in Gemc itself but I'm not sure the time window will be adjustable.
//The code creates one shape per hit and then the shapes add up if they are on the same pad
//It then fits the individual signals on each pad to determine the time
//The code creates an output root file containing a new tree with:
//  - cellID as before
//  - time newly defined
//  - ADC
// All hits are considered in the time window
// All units are ns

// First step creates the signal on each pad without integration from the electronics
// Then integrates the signal and stores 1 bin over NBinKept
// Then fits the signal produced and stores the data obtained from it


#include "TCanvas.h"
#include "TStyle.h"
#include "TH1.h"
#include "TGaxis.h"
#include "TRandom.h"
#include "TFile.h"
#include "TGraph.h"
#include "TF1.h"
#include "TMath.h"
#include "TTree.h"
#include "TROOT.h"
#include "TChain.h"
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <locale>
#include <sstream>
#include <algorithm>


using namespace std;

TRandom *delta = new TRandom();

TRandom *noise = new TRandom();
const double sigTelec = 5; // 5 ns uncertainty on the signal
const double sigTVdrift = 5; // 5 ns uncertainty on the

double EtoS(double, double, double); // convert energy deposited into signal on a pad
double noise_elec(double);
double drift_V(double);

int delta_t;

//void bonus_shaping()
void shaping()
{
    gROOT->SetStyle("Plain");
    
    //______________________________________________________________________________________________
    //  __________________________________________ Variables _________________________________________
    //______________________________________________________________________________________________
    
    char filename[100];
    
    const int StepSize = 10; // step size of the signal before integration (arbitrary value)
    const int BinSize = 40; // electronics integrates the signal over 40 ns
    const int NBinKept = 3; // only 1 bin over 3 is kept by the daq
    const int TrigWindSize = 10000; // Trigger window should be 10 micro
    const int NTrigSampl = TrigWindSize/BinSize; // number of time samples
    
    map<int, double[TrigWindSize]> R_adc; // Raw depositions for CellID, ADC
    
    vector<int> PadN;  // used to read only cell with signal, one entry for each hit
    vector<int> PadNum;// used to read only cell with signal, one entry for each cell
    vector<int> TidNum;
    vector<int> PidNum;
    vector<double> TimeNum;
    vector<double> ShiftNum;
    vector<double> VzNum;
    
    
    int count=0;
    int TIDprev = 0;
    
    double inte=0;
    
    double time_inte=0;
    double inte_tot; // integral of the signal in BinSize
    double max_inte=0; // maximum of the integral to help the fit
    double max_t=0;    // time at the maximum to help the fit

    
    bool flag_event = false;
    
    //______________________________________________________________________________________________
    //  ___________________________________ Canvas and histograms ____________________________________
    //______________________________________________________________________________________________
    
    TGraph *g_pad_inte = new TGraph();
    
    // Note used for fitting
    TF1 *doublegaus = new TF1("doublegaus","[0]+[2]*exp(-(x-[1])^2/(2*((x<[1])*[3]+(x>=[1])*[4])^2))/(0.5*([3]+[4])*sqrt(2*pi))",0,10000);
    doublegaus->SetParameters(0.0, 1, 5000.0, 100.0, 90.0);
    
    TF1 *mgaus = new TF1("mgaus","[0]*exp(-(x-[1])^2/(2*[2]^2))/([2]*sqrt(2*pi))",0,10000);
    
    
    //______________________________________________________________________________________________
    //  __________________________________________ Openings __________________________________________
    //______________________________________________________________________________________________
    
    sprintf(filename,"gemc.root");
    
    TFile *myfile= TFile::Open(filename);
    TTree *tree= (TTree*)myfile->Get("RTPC");
    
    std::vector<int>* CellID = 0;
    std::vector<double>* Time = 0;
    std::vector<double>* totEdep = 0;
    std::vector<int>* pid = 0;
    std::vector<int>* tid = 0;
    std::vector<double>* z_track = 0;
    std::vector<double>* time_shift = 0;
    
    TBranch *bcellID = 0;
    tree->SetBranchAddress("CellID",&CellID,&bcellID);
    
    TBranch *btime = 0;
    tree->SetBranchAddress("Time",&Time,&btime);
    
    TBranch *btime_shift = 0;
    tree->SetBranchAddress("TimeShift",&time_shift,&btime_shift);
    
    TBranch *btotEdep = 0;
    tree->SetBranchAddress("totEdep",&totEdep,&btotEdep);
    
    TBranch *bpid = 0;
    tree->SetBranchAddress("pid",&pid,&bpid);
    
    TBranch *btid = 0;
    tree->SetBranchAddress("TrackID",&tid,&btid);
    
    TBranch *bz_track = 0;
    tree->SetBranchAddress("vz",&z_track,&bz_track);
    
    
    
    TTree *gene= (TTree*)myfile->Get("generated");
    
    std::vector<int>* pid_g = 0;
    std::vector<double>* px = 0;
    std::vector<double>* py = 0;
    std::vector<double>* pz = 0;
    std::vector<double>* vz = 0;
    
    TBranch *bpid_g = 0;
    gene->SetBranchAddress("pid",&pid_g,&bpid_g);
    
    TBranch *bpx = 0;
    gene->SetBranchAddress("px",&px,&bpx);
    
    TBranch *bpy = 0;
    gene->SetBranchAddress("py",&py,&bpy);
    
    TBranch *bpz = 0;
    gene->SetBranchAddress("pz",&pz,&bpz);
    
    TBranch *bvz = 0;
    gene->SetBranchAddress("vz",&vz,&bvz);
    
    
    // The root file contains a tree
    // this tree will contain the hit pads and the signal
    TFile *fout = TFile::Open("shaping.root","RECREATE");
    TTree *Rec = new TTree("RTPC","RTPC data");
    TTree *Gen = new TTree("generated","Generated data");
    
    int evn;
    std::vector<int>* Pad_daq = 0;
    std::vector<int>* pid_daq = 0;
    std::vector<int>* trackid_daq = 0;
    std::vector<double>* ADC_daq = 0;
    std::vector<double>* Time_daq = 0;
    std::vector<double>* TimeShift_daq = 0;
    std::vector<double>* vz_daq = 0;
    
    Rec->Branch("evn",&evn,"evn/I");
    Rec->Branch("pid",&pid_daq);
    Rec->Branch("TrackID",&trackid_daq);
    Rec->Branch("CellID",&Pad_daq);
    Rec->Branch("ADC",&ADC_daq);
    Rec->Branch("Time",&Time_daq);
    Rec->Branch("TimeShift",&TimeShift_daq);
    Rec->Branch("vz",&vz_daq);
    
    
    int evn_v;
    std::vector<int>* pid_v = 0;
    std::vector<double>* z_v = 0;
    std::vector<double>* p_v = 0;
    std::vector<double>* pt_v = 0;
    std::vector<double>* th_v = 0;
    std::vector<double>* phi_v = 0;
    
    Gen->Branch("evn_v",&evn_v,"evn/I");
    Gen->Branch("pid",&pid_v);
    Gen->Branch("z_v",&z_v);
    Gen->Branch("p_v",&p_v);
    Gen->Branch("pt_v",&pt_v);
    Gen->Branch("th_v",&th_v);
    Gen->Branch("phi_v",&phi_v);
    
    //______________________________________________________________________________________________
    //  __________________________________________ Readings __________________________________________
    //______________________________________________________________________________________________
    
    
    //--Creating the signal on pads with 1 ns steps.
    for(int eve=0; eve<tree->GetEntries(); eve++){
        // new z well outside RTPC
        
        cout << "Event: " << eve+1 << " of " << gene->GetEntries() << " events." << endl;
        
        // Initializations
        R_adc.clear(); // Raw depositions for CellID, adc
        PadN.clear();
        PadNum.clear();
        TimeNum.clear();
        ShiftNum.clear();
        VzNum.clear();
        
        flag_event=false;
        
        Pad_daq->clear();
        ADC_daq->clear();  // not reliable for now as the fit fails often
        Time_daq->clear();
        trackid_daq->clear();
        TimeShift_daq->clear();
        vz_daq->clear();
        pid_daq->clear();
        
        
        pid_v->clear();
        z_v->clear();
        p_v->clear();
        pt_v->clear();
        th_v->clear();
        phi_v->clear();
        
        Long64_t tentry2 = gene->LoadTree(eve);
        bpid_g->GetEntry(tentry2);
        bpx->GetEntry(tentry2);
        bpy->GetEntry(tentry2);
        bpz->GetEntry(tentry2);
        bvz->GetEntry(tentry2);
        
        Long64_t tentry = tree->LoadTree(eve);
        bcellID->GetEntry(tentry);
        btime->GetEntry(tentry);
        btotEdep->GetEntry(tentry);
        bpid->GetEntry(tentry);
        bz_track->GetEntry(tentry);
        btid->GetEntry(tentry);
        btime_shift->GetEntry(tentry);
        
        for(int k=0; k<pid_g->size(); k++){
            if(k==0) cout << "Storing vertex data..." << endl;
            evn_v = eve;
            pid_v->push_back(pid_g->at(k));
            z_v->push_back(vz->at(k));
            p_v->push_back(TMath::Sqrt(px->at(k)*px->at(k)+py->at(k)*py->at(k)+pz->at(k)*pz->at(k)));
            pt_v->push_back(TMath::Sqrt(px->at(k)*px->at(k)+py->at(k)*py->at(k)));
            th_v->push_back(TMath::ATan2(pt_v->at(k),pz->at(k)));
            phi_v->push_back(TMath::ATan2(py->at(k),px->at(k)));
            
        }
        
        if(CellID->size()!=0){
            
            for(int c = 0; c<CellID->size(); c++){
                if(c==0) {
                    cout << "Fitting signal on pads...CellID size= " << CellID->size() << endl;
                    TIDprev = tid->at(c);
                }
                
                //if(pid->at(c)!=2212) continue;
                
                // ensure fitting only occurs for same track
                if(tid->at(c) == TIDprev){
                    std::vector<int>::iterator it;
                    it = find(PadN.begin(), PadN.end(), CellID->at(c)); // searches in PadN if CellID already exists
                    
                    if(it!=PadN.end()){ // this pad has already seen signal
                        for(int t=0;t<TrigWindSize;t+=StepSize){
                            R_adc[CellID->at(c)][t] += EtoS(Time->at(c),t,totEdep->at(c));
                        }
                    }
                    else{ // first signal on this pad
                        for(int t=0;t<TrigWindSize;t+=StepSize){
                            R_adc[CellID->at(c)][t] = EtoS(Time->at(c),t,totEdep->at(c));
                        }
                        PadNum.push_back(CellID->at(c));
                        TidNum.push_back(tid->at(c));
                        PidNum.push_back(pid->at(c));
                        TimeNum.push_back(Time->at(c));
                        ShiftNum.push_back(time_shift->at(c));
                        VzNum.push_back(z_track->at(c));
                    }

                }
                else TIDprev = tid->at(c);
                
                PadN.push_back(CellID->at(c));
                
            } // c
            
            
            //--Signal created on pads with 1 ns steps
            
            
            //--For each pad
            // Read the signal on it
            // Integrates it into BinSize long bins
            // Keeps only 1 bin over 3
            
            inte=0;
            for(unsigned int p=0;p<PadNum.size();p++){
                time_inte=0.;
                inte_tot=0;
                for(int t=0;t<TrigWindSize;t+=StepSize){
                    if(p==0 && t==0) cout << "Integrating signal..." << endl;
                    
                    time_inte = t;//TimeNum[p]+t;
                    
                    // Add this step to the integral of the signal
                    if(t==0) inte = R_adc[PadNum[p]][t]*StepSize;
                    else inte+=0.5*(R_adc[PadNum[p]][t-1]+R_adc[PadNum[p]][t])*StepSize;
                    
                    inte_tot+=inte;
                    
                    if(t%BinSize==0 && 0<t){ // integration over BinSize
                        if(t%BinSize*NBinKept==0){ // one BinSize over NBinKept is read out, hence added to the histogram
                            g_pad_inte->SetPoint(time_inte/(BinSize*NBinKept),time_inte,inte);
                            if(max_inte<inte){
                                max_inte=inte;
                                max_t=time_inte;
                                //cout << "Max Time= " << max_t << endl;
                            }
                        }
                        inte=0; // the integral of 40 ns bin must be reset to 0
                    }
                }
                //max_t=TimeNum[p];
                mgaus->SetParameters(5*max_inte,max_t,165.);
                mgaus->SetRange(max_t-StepSize*4,max_t+StepSize*4);
                g_pad_inte->Fit("mgaus","QER");
                
                // Filling output tree
                evn = eve;
                if(0<mgaus->GetParameter(1) && mgaus->GetParameter(1)<10000){ // the fit is not robust for now
                    Pad_daq->push_back(PadNum[p]);
                    ADC_daq->push_back(inte_tot);  // use the signal integral for now
                    Time_daq->push_back(mgaus->GetParameter(1));
                    vz_daq->push_back(VzNum[p]);
                    pid_daq->push_back(PidNum[p]);
                    trackid_daq->push_back(TidNum[p]);
                    TimeShift_daq->push_back(ShiftNum[p]);
                    flag_event=true;
                }
                
                // Cleaning
                max_inte=0;
                max_t=0;
                if(eve!=tree->GetEntries()-1){
                    for(int ii=g_pad_inte->GetN();ii>-1;ii--){
                        g_pad_inte->RemovePoint(ii);
                    }
                }
                
            }
            //--each pad has its daq signal associated to it
            
            //______________________________________________________________________________________________
            //  _________________________________________ Cleaning ___________________________________________
            //______________________________________________________________________________________________
            
            if(eve!=tree->GetEntries()-1){
                for(int ii=g_pad_inte->GetN();ii>-1;ii--){
                    g_pad_inte->RemovePoint(ii);
                }
            }
            
        }
        else{ // CellID->size()==0
            evn = eve;
            Pad_daq->push_back(-999);
            ADC_daq->push_back(0);
            Time->push_back(-999);
            TimeShift_daq->push_back(-999);
            vz_daq->push_back(-999);
            pid_daq->push_back(-999);
            trackid_daq->push_back(-999);
        }
        cout << "CellID size= " << Pad_daq->size() << "Time size= " << Time_daq->size() << endl;
        Gen->Fill();
        Rec->Fill();
        
        
    } // eve
    
    fout->Write();
    
    delete fout;
    
    return;
    
    
}


double EtoS(double tini, double t, double e_tot){
    
    double sig;
    
    t = noise_elec(t);    // change t to simulate the electronics noise, also modifies the amplitude
    double p0 = 0.0;
    double p2 = 178.158;
    double p3 = 165.637;
    double p4 = 165.165;
    
    if(t<tini) sig = p0+e_tot*p2*exp(-(t-tini)*(t-tini)/(2*p3*p3))/(0.5*(p3+p4)*sqrt(2*TMath::Pi()));
    else       sig = p0+e_tot*p2*exp(-(t-tini)*(t-tini)/(2*p4*p4))/(0.5*(p3+p4)*sqrt(2*TMath::Pi()));
    
    return sig;
    
}

double noise_elec(double tim){
    
    return noise->Gaus(tim,sigTelec);
    
}

double drift_V(double tim){
    
    return noise->Gaus(tim,sigTVdrift);
    
}

