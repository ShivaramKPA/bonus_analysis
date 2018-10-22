///////////////////////////////////////////////////////////////////
//
// Nate Dzbenski
// Gas mixture analysis
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
#include "TLegend.h"
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <locale>
#include <sstream>
#include <math.h>

void bfieldmap_ana()
{
    gROOT->SetStyle("Plain");
    
    //________________________________________________________________________________________________
    // ____________________________________________Variables__________________________________________
    //________________________________________________________________________________________________

    float zneg19_phi[9] = {0.980011,
        0.837412,
        0.678139,
        0.536768,
        0.411291,
        0.295843,
        0.18943,
        0.0910976,
        0.0175907};
    
    float zneg19_t[9] = {4794.48,
        4486.93,
        4021.54,
        3510.23,
        2932.6,
        2287.31,
        1587.45,
        828.128,
        171.318};
    
    float zneg15_phi[9] = {0.990828,
        0.842372,
        0.685084,
        0.542275,
        0.41465,
        0.297506,
        0.19121,
        0.0917766,
        0.0179327};
    
    float zneg15_t[9] = {4848.16,
        4517.17,
        4062.49,
        3545.54,
        2959.08,
        2318.55,
        1608.88,
        833.545,
        171.688};
    
    float zneg10_phi[9] = {0.995662,
        0.848969,
        0.688194,
        0.545092,
        0.416746,
        0.300426,
        0.191803,
        0.0928831,
        0.0179778};
    
    float zneg10_t[9] = {4879.49,
        4554.65,
        4088.18,
        3571.7,
        2990.3,
        2329.17,
        1612.93,
        842.903,
        173.312};
    
    float zneg5_phi[9] = {0.997039,
        0.850915,
        0.689582,
        0.545335,
        0.417162,
        0.298644,
        0.191962,
        0.0929194,
        0.0181597};
    
    float zneg5_t[9] = {4888.32,
        4562.65,
        4096.14,
        3571.72,
        2984.51,
        2331.96,
        1615.02,
        844.293,
        172.865};
    
    float z0_phi[9] = {0.998838,
        0.85015,
        0.687601,
        0.545972,
        0.417705,
        0.298937,
        0.191559,
        0.0922974,
        0.0178475};
    
    float z0_t[9] = {4886,
        4561.83,
        4098.4,
        3570.77,
        2985.64,
        2330.14,
        1616.56,
        838.226,
        173.386};
    
    float z5_phi[9] = {0.995865,
        0.852157,
        0.68922,
        0.544917,
        0.416418,
        0.299752,
        0.192174,
        0.0927076,
        0.0179557};
    
    float z5_t[9] = {4878.8,
        4560.28,
        4095.34,
        3571.08,
        2977.79,
        2328.88,
        1613.33,
        838.968,
        173.207};
    
    float z10_phi[9] = {0.993505,
        0.847097,
        0.687927,
        0.544612,
        0.416464,
        0.29927,
        0.191271,
        0.0929344,
        0.0178747};
    
    float z10_t[9] = {4866.17,
        4534.6,
        4084.8,
        3556.08,
        2972.8,
        2327.71,
        1611.17,
        841.333,
        172.16};
    
    float z15_phi[9] = {0.993413,
        0.842892,
        0.681887,
        0.5423,
        0.413491,
        0.296665,
        0.190428,
        0.0914498,
        0.0179528};
    
    float z15_t[9] = {4868.64,
        4513.78,
        4047.89,
        3537.28,
        2952.37,
        2306.01,
        1606.77,
        833.152,
        173.03};
    
    float z19_phi[9] = {0.976629,
        0.833081,
        0.67597,
        0.536253,
        0.40962,
        0.294246,
        0.188567,
        0.0910907,
        0.0177384};
    
    float z19_t[9] = {4772.85,
        4447.88,
        4005.49,
        3493.07,
        2917.44,
        2280.85,
        1579.45,
        822.055,
        170.68};
    
    //errors
    float zneg19_phi_err[9] = {0.00900979,
        0.00929756,
        0.00694851,
        0.00583257,
        0.0051429,
        0.00463225,
        0.00335575,
        0.00333426,
        0.000974565};
    
    float zneg19_t_err[9] = {40.5679,
        42.8955,
        29.2347,
        34.5895,
        33.1942,
        40.4061,
        22.198,
        20.1889,
        9.51966};
    
    float zneg15_phi_err[9] = {0.00960376,
        0.0108475,
        0.00772402,
        0.00712173,
        0.00546197,
        0.00477844,
        0.00352452,
        0.0025788,
        0.00101236};
    
    float zneg15_t_err[9] = {35.1561,
        42.7189,
        31.8874,
        34.307,
        32.2444,
        38.4698,
        24.3454,
        17.3399,
        7.97082};
    
    float zneg10_phi_err[9] = {0.0103908,
        0.00940146,
        0.0072093,
        0.00910251,
        0.00498041,
        0.00412408,
        0.00444955,
        0.00187424,
        0.0016163};
    
    float zneg10_t_err[9] = {49.426,
        31.3864,
        38.9165,
        32.4095,
        43.9832,
        24.4538,
        24.4279,
        15.4337,
        9.87914};
    
    float zneg5_phi_err[9] = {0.0104792,
        0.00763382,
        0.00696928,
        0.00679426,
        0.00589499,
        0.00508229,
        0.00303503,
        0.00262334,
        0.00112457};
    
    float zneg5_t_err[9] = {36.8646,
        34.474,
        25.8173,
        37.8122,
        32.4846,
        32.031,
        27.1239,
        15.1818,
        7.94243};
    
    float z0_phi_err[9] = {0.0153287,
        0.0100367,
        0.0086736,
        0.00538724,
        0.00547198,
        0.00541904,
        0.00451889,
        0.00267296,
        0.00161708};
    
    float z0_t_err[9] = {42.0425,
        37.4574,
        24.8777,
        32.5755,
        36.8723,
        31.0018,
        21.2976,
        13.8519,
        8.54205};

    float z5_phi_err[9] = {0.0102535,
        0.00818188,
        0.0067975,
        0.00582456,
        0.00525655,
        0.00494794,
        0.00386793,
        0.00241614,
        0.00127802};
    
    float z5_t_err[9] = {35.8543,
        32.3121,
        39.6636,
        29.1931,
        38.4252,
        25.4011,
        27.0958,
        20.9858,
        9.74676};
    
    float z10_phi_err[9] = {0.00838353,
        0.00836322,
        0.00777414,
        0.00679242,
        0.00557131,
        0.00438086,
        0.00462493,
        0.00198503,
        0.00106549};
    
    float z10_t_err[9] = {38.2701,
        30.1537,
        33.3057,
        29.5365,
        28.0308,
        28.2062,
        35.7995,
        22.0911,
        9.8714};
    
    float z15_phi_err[9] = {0.00787928,
        0.0123061,
        0.00589799,
        0.0055731,
        0.00595455,
        0.00447543,
        0.00362025,
        0.00243182,
        0.000758088};
    
    float z15_t_err[9] = {45.2231,
        48.8116,
        26.0933,
        28.2662,
        34.8137,
        30.7452,
        29.6279,
        13.6175,
        5.15297};
    
    float z19_phi_err[9] = {0.00843516,
        0.00738115,
        0.00759956,
        0.0085201,
        0.00538106,
        0.00354116,
        0.00476369,
        0.00242365,
        0.00129273};
    
    float z19_t_err[9] = {25.1225,
        30.7316,
        32.8523,
        37.3773,
        32.413,
        23.1283,
        24.5472,
        17.7795,
        7.26676};
    
    
    // error on sigma
    float zneg19_sphi_err[9] = {0.00169665,
        0.00217279,
        0.00165936,
        0.00174571,
        0.000956303,
        0.000756981,
        0.000566057,
        8.09E-04,
        1.39E-04};
    
    float zneg15_sphi_err[9] = {0.00190519,
        0.00259419,
        0.00139945,
        0.00118923,
        0.000374621,
        0.000804529,
        0.000774409,
        6.73E-04,
        1.25E-04};
    
    float zneg10_sphi_err[9] = {0.00208574,
        0.00174064,
        0.000533065,
        0.00244649,
        0.00124295,
        0.00111156,
        0.00108335,
        1.65E-04,
        3.56E-04};
    
    float zneg5_sphi_err[9] = {0.00194006,
        0.00144419,
        0.00139117,
        0.00141667,
        0.000966202,
        0.000989933,
        0.00048698,
        2.95E-04,
        1.34E-04};
    
    float z0_sphi_err[9] = {0.0062012,
        0.0022705,
        0.00176774,
        0.00096853,
        0.00115863,
        0.00105666,
        0.000832015,
        6.25E-04,
        4.40E-04};
    
    float z5_sphi_err[9] = {0.00242292,
        0.00159916,
        0.00125963,
        0.000906184,
        0.00106825,
        0.00116369,
        0.00075764,
        3.68E-04,
        2.71E-04};
    
    float z10_sphi_err[9] = {0.00152226,
        0.00181872,
        0.00202506,
        0.00155767,
        0.000963057,
        0.00037309,
        0.000904094,
        3.67E-04,
        1.69E-04};
    
    float z15_sphi_err[9] = {0.00109183,
        0.00304634,
        0.000973293,
        0.00124964,
        0.0015755,
        0.00102381,
        0.000667755,
        4.71E-04,
        1.19E-04};
    
    float z19_sphi_err[9] = {0.00176278,
        0.0014858,
        0.00143059,
        0.00198016,
        0.000956183,
        0.000538472,
        0.00097437,
        4.88E-04,
        2.75E-04};
    
    float zneg19_st_err[9] = {10.1294,
        12.1273,
        4.95144,
        7.19185,
        7.6024,
        12.714,
        5.62251,
        4.20491,
        1.91055};
    
    float zneg15_st_err[9] = {5.56684,
        8.192,
        6.60044,
        8.75462,
        2.9584,
        10.3464,
        5.08767,
        3.00551,
        0.727335};
    
    float zneg10_st_err[9] = {12.1489,
        4.6737,
        2.27355,
        7.10417,
        11.4963,
        4.42385,
        5.84465,
        1.39053,
        2.40323};
    
    float zneg5_st_err[9] = {6.75859,
        6.39716,
        3.51514,
        7.78763,
        6.3207,
        9.08278,
        5.75799,
        1.26386,
        0.776315};
    
    float z0_st_err[9] = {8.9227,
        6.88496,
        3.58453,
        6.60128,
        8.13226,
        5.74539,
        3.99874,
        1.95276,
        1.75562};
    
    float z5_st_err[9] = {6.85307,
        5.23438,
        8.12424,
        6.02257,
        9.024,
        4.55003,
        5.03408,
        4.9939,
        2.1525};
    
    float z10_st_err[9] = {6.20501,
        5.69988,
        6.17566,
        5.32207,
        4.79438,
        2.24644,
        11.7752,
        7.20071,
        1.85747};
    
    float z15_st_err[9] = {9.73181,
        14.5631,
        4.74398,
        4.90346,
        7.66725,
        6.74596,
        5.37632,
        2.50897,
        3.18366};
    
    float z19_st_err[9] = {6.0066,
        6.54693,
        6.45825,
        9.80349,
        8.62312,
        4.04627,
        4.17741,
        3.1704,
        2.12024};
    
    
    float r[9] = {3.1, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 6.9};
    float r_err[9] = {0,0,0,0,0,0,0,0,0};
    
    float z[9] = {-19,-15,-10,-5,0,5,10,15,19};
    float z_err[9] = {0,0,0,0,0,0,0,0,0};
    
    float max_t[9]= {zneg19_t[0],zneg15_t[0],zneg10_t[0],zneg5_t[0],z0_t[0],z5_t[0],z10_t[0],z15_t[0],z19_t[0]};
    float max_t_err[9]= {zneg19_t_err[0],zneg15_t_err[0],zneg10_t_err[0],zneg5_t_err[0],z0_t_err[0],z5_t_err[0],z10_t_err[0],z15_t_err[0],z19_t_err[0]};
    
    float max_phi[9]= {zneg19_phi[0],zneg15_phi[0],zneg10_phi[0],zneg5_phi[0],z0_phi[0],z5_phi[0],z10_phi[0],z15_phi[0],z19_phi[0]};
    float max_phi_err[9]= {zneg19_phi_err[0],zneg15_phi_err[0],zneg10_phi_err[0],zneg5_phi_err[0],z0_phi_err[0],z5_phi_err[0],z10_phi_err[0],z15_phi_err[0],z19_phi_err[0]};
    
    //________________________________________________________________________________________________
    //_____________________________________ Canvas and histograms ____________________________________
    //________________________________________________________________________________________________
    

    TCanvas *c1 = new TCanvas("c1","t_d vs phi_d",800,600);
    TCanvas *c2 = new TCanvas("c2","phi_d vs r",800,600);
    TCanvas *c3 = new TCanvas("c3","t_d vs r",800,600);
    TCanvas *c4 = new TCanvas("c4","sphi_d vs r",800,600);
    TCanvas *c5= new TCanvas("c5","st_d vs r",800,600);
    TCanvas *c6= new TCanvas("c6","t_d vs z",800,600);
    TCanvas *c7= new TCanvas("c7","phi_d vs z",800,600);

    // for the t vs phi graph
    TGraphErrors *gr_zneg19 = new TGraphErrors(9,zneg19_phi,zneg19_t,zneg19_phi_err,zneg19_t_err);
    TGraphErrors *gr_zneg15 = new TGraphErrors(9,zneg15_phi,zneg15_t,zneg15_phi_err,zneg15_t_err);
    TGraphErrors *gr_zneg10 = new TGraphErrors(9,zneg10_phi,zneg10_t,zneg10_phi_err,zneg10_t_err);
    TGraphErrors *gr_zneg5 = new TGraphErrors(9,zneg5_phi,zneg5_t,zneg5_phi_err,zneg5_t_err);
    TGraphErrors *gr_z0 = new TGraphErrors(9,z0_phi,z0_t,z0_phi_err,z0_t_err);
    TGraphErrors *gr_z5 = new TGraphErrors(9,z5_phi,z5_t,z5_phi_err,z5_t_err);
    TGraphErrors *gr_z10 = new TGraphErrors(9,z10_phi,z10_t,z10_phi_err,z10_t_err);
    TGraphErrors *gr_z15 = new TGraphErrors(9,z15_phi,z15_t,z15_phi_err,z15_t_err);
    TGraphErrors *gr_z19 = new TGraphErrors(9,z19_phi,z19_t,z19_phi_err,z19_t_err);
    
    // phi vs r
    TGraphErrors *gr_phir19 = new TGraphErrors(9,r,z19_phi,r_err,z19_phi_err);
    TGraphErrors *gr_phir15 = new TGraphErrors(9,r,z15_phi,r_err,z15_phi_err);
    TGraphErrors *gr_phir10 = new TGraphErrors(9,r,z10_phi,r_err,z10_phi_err);
    TGraphErrors *gr_phir5 = new TGraphErrors(9,r,z5_phi,r_err,z5_phi_err);
    TGraphErrors *gr_phir0 = new TGraphErrors(9,r,z0_phi,r_err,z0_phi_err);
    TGraphErrors *gr_phirneg5 = new TGraphErrors(9,r,zneg5_phi,r_err,zneg5_phi_err);
    TGraphErrors *gr_phirneg10 = new TGraphErrors(9,r,zneg10_phi,r_err,zneg10_phi_err);
    TGraphErrors *gr_phirneg15 = new TGraphErrors(9,r,zneg15_phi,r_err,zneg15_phi_err);
    TGraphErrors *gr_phirneg19 = new TGraphErrors(9,r,zneg19_phi,r_err,zneg19_phi_err);
    
    // time vs r
    TGraphErrors *gr_tr19 = new TGraphErrors(9,r,z19_t,r_err,z19_t_err);
    TGraphErrors *gr_tr15 = new TGraphErrors(9,r,z15_t,r_err,z15_t_err);
    TGraphErrors *gr_tr10 = new TGraphErrors(9,r,z10_t,r_err,z10_t_err);
    TGraphErrors *gr_tr5 = new TGraphErrors(9,r,z5_t,r_err,z5_t_err);
    TGraphErrors *gr_tr0 = new TGraphErrors(9,r,z0_t,r_err,z0_t_err);
    TGraphErrors *gr_trneg5 = new TGraphErrors(9,r,zneg5_t,r_err,zneg5_t_err);
    TGraphErrors *gr_trneg10 = new TGraphErrors(9,r,zneg10_t,r_err,zneg10_t_err);
    TGraphErrors *gr_trneg15 = new TGraphErrors(9,r,zneg15_t,r_err,zneg15_t_err);
    TGraphErrors *gr_trneg19 = new TGraphErrors(9,r,zneg19_t,r_err,zneg19_t_err);
    
    // phi vs r
    TGraphErrors *gr_sphir19 = new TGraphErrors(9,r,z19_phi_err,r_err,z19_sphi_err);
    TGraphErrors *gr_sphir15 = new TGraphErrors(9,r,z15_phi_err,r_err,z15_sphi_err);
    TGraphErrors *gr_sphir10 = new TGraphErrors(9,r,z10_phi_err,r_err,z10_sphi_err);
    TGraphErrors *gr_sphir5 = new TGraphErrors(9,r,z5_phi_err,r_err,z5_sphi_err);
    TGraphErrors *gr_sphir0 = new TGraphErrors(9,r,z0_phi_err,r_err,z0_sphi_err);
    TGraphErrors *gr_sphirneg5 = new TGraphErrors(9,r,zneg5_phi_err,r_err,zneg5_sphi_err);
    TGraphErrors *gr_sphirneg10 = new TGraphErrors(9,r,zneg10_phi_err,r_err,zneg10_sphi_err);
    TGraphErrors *gr_sphirneg15 = new TGraphErrors(9,r,zneg15_phi_err,r_err,zneg15_sphi_err);
    TGraphErrors *gr_sphirneg19 = new TGraphErrors(9,r,zneg19_phi_err,r_err,zneg19_sphi_err);
    
    // time vs r
    TGraphErrors *gr_str19 = new TGraphErrors(9,r,z19_t_err,r_err,z19_st_err);
    TGraphErrors *gr_str15 = new TGraphErrors(9,r,z15_t_err,r_err,z15_st_err);
    TGraphErrors *gr_str10 = new TGraphErrors(9,r,z10_t_err,r_err,z10_st_err);
    TGraphErrors *gr_str5 = new TGraphErrors(9,r,z5_t_err,r_err,z5_st_err);
    TGraphErrors *gr_str0 = new TGraphErrors(9,r,z0_t_err,r_err,z0_st_err);
    TGraphErrors *gr_strneg5 = new TGraphErrors(9,r,zneg5_t_err,r_err,zneg5_st_err);
    TGraphErrors *gr_strneg10 = new TGraphErrors(9,r,zneg10_t_err,r_err,zneg10_st_err);
    TGraphErrors *gr_strneg15 = new TGraphErrors(9,r,zneg15_t_err,r_err,zneg15_st_err);
    TGraphErrors *gr_strneg19 = new TGraphErrors(9,r,zneg19_t_err,r_err,zneg19_st_err);
    
    TGraphErrors *gr_td_z = new TGraphErrors(9,z,max_t,z_err,max_t_err);
    TGraphErrors *gr_phid_z = new TGraphErrors(9,z,max_phi,z_err,max_phi_err);
    
    TMultiGraph *gr = new TMultiGraph();
    TMultiGraph *gr_phi_r = new TMultiGraph();
    TMultiGraph *gr_t_r = new TMultiGraph();
    
    TMultiGraph *gr_sphi_r = new TMultiGraph();
    TMultiGraph *gr_st_r = new TMultiGraph();
    
    TF1  *f1 = new TF1("f1","[0]*(7-x)+[1]*(7-x)*(7-x)",0,7);
    
    //________________________________________________________________________________________________
    // ___________________________________________ Openings __________________________________________
    //________________________________________________________________________________________________

    
    //________________________________________________________________________________________________
    // ___________________________________________ Readings __________________________________________
    //________________________________________________________________________________________________
    
    
    
    //________________________________________________________________________________________________
    // _____________________________________________ Displays ________________________________________
    //________________________________________________________________________________________________
    
    // time vs phi
    c1->cd();
    
    gr_zneg19->SetMarkerColor(1);
    gr_zneg19->SetMarkerStyle(21);
    gr_zneg19->SetMarkerSize(0.5);
    
    gr_zneg15->SetMarkerColor(1);
    gr_zneg15->SetMarkerStyle(21);
    gr_zneg15->SetMarkerSize(0.5);
    
    gr_zneg10->SetMarkerColor(2);
    gr_zneg10->SetMarkerStyle(21);
    gr_zneg10->SetMarkerSize(0.5);
    
    gr_zneg5->SetMarkerColor(3);
    gr_zneg5->SetMarkerStyle(21);
    gr_zneg5->SetMarkerSize(0.5);
    
    gr_z0->SetMarkerColor(4);
    gr_z0->SetMarkerStyle(21);
    gr_z0->SetMarkerSize(0.5);
    
    gr_z5->SetMarkerColor(5);
    gr_z5->SetMarkerStyle(21);
    gr_z5->SetMarkerSize(0.5);
    
    gr_z10->SetMarkerColor(6);
    gr_z10->SetMarkerStyle(21);
    gr_z10->SetMarkerSize(0.5);
    
    gr_z15->SetMarkerColor(7);
    gr_z15->SetMarkerStyle(21);
    gr_z15->SetMarkerSize(0.5);
    
    gr_z19->SetMarkerColor(8);
    gr_z19->SetMarkerStyle(21);
    gr_z19->SetMarkerSize(0.5);
    
    gr->SetTitle("t_{d} vs #phi_{d}");
    
    gr->Add(gr_zneg19);
    gr->Add(gr_zneg15);
    gr->Add(gr_zneg10);
    gr->Add(gr_zneg5);
    gr->Add(gr_z0);
    gr->Add(gr_z5);
    gr->Add(gr_z10);
    gr->Add(gr_z15);
    gr->Add(gr_z19);

    gr->Draw("AP");
    gr->GetXaxis()->SetTitle("#Delta #phi_{d} [rad]");
    gr->GetYaxis()->SetTitle("t_{d} [ns]");
    gr->GetYaxis()->SetTitleOffset(1.3);
    gr->SaveAs("figs/b_field_ana/t_phi.png");
    
    //gr->Fit("pol2", "F");
    //TF1 *fpol = gr->GetFunction("pol2");
    //fpol->SetLineWidth(1);
    
    TLegend *leg = new TLegend(0.8,0.2,0.9,0.5);
    leg->AddEntry(gr_zneg19,"z = -19 cm","lep");
    leg->AddEntry(gr_zneg15,"z = -15 cm","lep");
    leg->AddEntry(gr_zneg10,"z = -10 cm","lep");
    leg->AddEntry(gr_zneg5,"z = -5 cm","lep");
    leg->AddEntry(gr_z0,"z = 0 cm","lep");
    leg->AddEntry(gr_z5,"z = 5 cm","lep");
    leg->AddEntry(gr_z10,"z = 10 cm","lep");
    leg->AddEntry(gr_z15,"z = 15 cm","lep");
    leg->AddEntry(gr_z19,"z = 19 cm","lep");
    leg->Draw();
    
    
    // phi vs r
    c2->cd();
    
    gr_phir19->SetMarkerColor(1);
    gr_phir19->SetMarkerStyle(21);
    gr_phir19->SetMarkerSize(0.5);
    gr_phir15->SetMarkerColor(2);
    gr_phir15->SetMarkerStyle(21);
    gr_phir15->SetMarkerSize(0.5);
    gr_phir10->SetMarkerColor(3);
    gr_phir10->SetMarkerStyle(21);
    gr_phir10->SetMarkerSize(0.5);
    gr_phir5->SetMarkerColor(4);
    gr_phir5->SetMarkerStyle(21);
    gr_phir5->SetMarkerSize(0.5);
    gr_phir0->SetMarkerColor(5);
    gr_phir0->SetMarkerStyle(21);
    gr_phir0->SetMarkerSize(0.5);
    gr_phirneg5->SetMarkerColor(6);
    gr_phirneg5->SetMarkerStyle(21);
    gr_phirneg5->SetMarkerSize(0.5);
    gr_phirneg10->SetMarkerColor(7);
    gr_phirneg10->SetMarkerStyle(21);
    gr_phirneg10->SetMarkerSize(0.5);
    gr_phirneg15->SetMarkerColor(8);
    gr_phirneg15->SetMarkerStyle(21);
    gr_phirneg15->SetMarkerSize(0.5);
    gr_phirneg19->SetMarkerColor(9);
    gr_phirneg19->SetMarkerStyle(21);
    gr_phirneg19->SetMarkerSize(0.5);
    
    gr_phi_r->Add(gr_phir19);
    gr_phi_r->Add(gr_phir15);
    gr_phi_r->Add(gr_phir10);
    gr_phi_r->Add(gr_phir5);
    gr_phi_r->Add(gr_phir0);
    gr_phi_r->Add(gr_phirneg5);
    gr_phi_r->Add(gr_phirneg10);
    gr_phi_r->Add(gr_phirneg15);
    gr_phi_r->Add(gr_phirneg19);
    
    gr_phi_r->Fit("f1");
    
    gr_phi_r->SetTitle("#phi_{d} vs r");
    gr_phi_r->Draw("AP");
    gr_phi_r->GetXaxis()->SetTitle("r [cm]");
    gr_phi_r->GetYaxis()->SetTitle("#phi_{d} [rad]");
    gr_phi_r->SaveAs("figs/b_field_ana/phi_r.png");
    
    TLegend *leg1 = new TLegend(0.8,0.6,0.9,0.9);
    leg1->AddEntry(gr_phir19,"z = 19 cm","lep");
    leg1->AddEntry(gr_phir15,"z = 15 cm","lep");
    leg1->AddEntry(gr_phir10,"z = 10 cm","lep");
    leg1->AddEntry(gr_phir5,"z = 5 cm","lep");
    leg1->AddEntry(gr_phir0,"z = 0 cm","lep");
    leg1->AddEntry(gr_phirneg5,"z = -5 cm","lep");
    leg1->AddEntry(gr_phirneg10,"z = -10 cm","lep");
    leg1->AddEntry(gr_phirneg15,"z = -15 cm","lep");
    leg1->AddEntry(gr_phirneg19,"z = -19 cm","lep");
    leg1->Draw();
    
    // time vs r
    c3->cd();
    
    gr_tr19->SetMarkerColor(1);
    gr_tr19->SetMarkerStyle(21);
    gr_tr19->SetMarkerSize(0.5);
    gr_tr15->SetMarkerColor(2);
    gr_tr15->SetMarkerStyle(21);
    gr_tr15->SetMarkerSize(0.5);
    gr_tr10->SetMarkerColor(3);
    gr_tr10->SetMarkerStyle(21);
    gr_tr10->SetMarkerSize(0.5);
    gr_tr5->SetMarkerColor(4);
    gr_tr5->SetMarkerStyle(21);
    gr_tr5->SetMarkerSize(0.5);
    gr_tr0->SetMarkerColor(5);
    gr_tr0->SetMarkerStyle(21);
    gr_tr0->SetMarkerSize(0.5);
    gr_trneg5->SetMarkerColor(6);
    gr_trneg5->SetMarkerStyle(21);
    gr_trneg5->SetMarkerSize(0.5);
    gr_trneg10->SetMarkerColor(7);
    gr_trneg10->SetMarkerStyle(21);
    gr_trneg10->SetMarkerSize(0.5);
    gr_trneg15->SetMarkerColor(8);
    gr_trneg15->SetMarkerStyle(21);
    gr_trneg15->SetMarkerSize(0.5);
    gr_trneg19->SetMarkerColor(9);
    gr_trneg19->SetMarkerStyle(21);
    gr_trneg19->SetMarkerSize(0.5);
    
    gr_t_r->Add(gr_tr19);
    gr_t_r->Add(gr_tr15);
    gr_t_r->Add(gr_tr10);
    gr_t_r->Add(gr_tr5);
    gr_t_r->Add(gr_tr0);
    gr_t_r->Add(gr_trneg5);
    gr_t_r->Add(gr_trneg10);
    gr_t_r->Add(gr_trneg15);
    gr_t_r->Add(gr_trneg19);
    
    gr_t_r->Fit("f1");
    
    gr_t_r->SetTitle("t_{d} vs r");
    gr_t_r->Draw("AP");
    gr_t_r->GetXaxis()->SetTitle("r [cm]");
    gr_t_r->GetYaxis()->SetTitle("t_{d} [ns]");
    gr_t_r->GetYaxis()->SetTitleOffset(1.3);
    gr_t_r->SaveAs("figs/b_field_ana/t_r.png");
    
    TLegend *leg2 = new TLegend(0.8,0.6,0.9,0.9);
    leg2->AddEntry(gr_tr19,"z = 19 cm","lep");
    leg2->AddEntry(gr_tr15,"z = 15 cm","lep");
    leg2->AddEntry(gr_tr10,"z = 10 cm","lep");
    leg2->AddEntry(gr_tr5,"z = 5 cm","lep");
    leg2->AddEntry(gr_tr0,"z = 0 cm","lep");
    leg2->AddEntry(gr_trneg5,"z = -5 cm","lep");
    leg2->AddEntry(gr_trneg10,"z = -10 cm","lep");
    leg2->AddEntry(gr_trneg15,"z = -15 cm","lep");
    leg2->AddEntry(gr_trneg19,"z = -19 cm","lep");
    leg2->Draw();
    
    
    // sigma phi vs r
    c4->cd();
    
    gr_sphir19->SetMarkerColor(1);
    gr_sphir19->SetMarkerStyle(21);
    gr_sphir19->SetMarkerSize(0.5);
    gr_sphir15->SetMarkerColor(2);
    gr_sphir15->SetMarkerStyle(21);
    gr_sphir15->SetMarkerSize(0.5);
    gr_sphir10->SetMarkerColor(3);
    gr_sphir10->SetMarkerStyle(21);
    gr_sphir10->SetMarkerSize(0.5);
    gr_sphir5->SetMarkerColor(4);
    gr_sphir5->SetMarkerStyle(21);
    gr_sphir5->SetMarkerSize(0.5);
    gr_sphir0->SetMarkerColor(5);
    gr_sphir0->SetMarkerStyle(21);
    gr_sphir0->SetMarkerSize(0.5);
    gr_sphirneg5->SetMarkerColor(6);
    gr_sphirneg5->SetMarkerStyle(21);
    gr_sphirneg5->SetMarkerSize(0.5);
    gr_sphirneg10->SetMarkerColor(7);
    gr_sphirneg10->SetMarkerStyle(21);
    gr_sphirneg10->SetMarkerSize(0.5);
    gr_sphirneg15->SetMarkerColor(8);
    gr_sphirneg15->SetMarkerStyle(21);
    gr_sphirneg15->SetMarkerSize(0.5);
    gr_sphirneg19->SetMarkerColor(9);
    gr_sphirneg19->SetMarkerStyle(21);
    gr_sphirneg19->SetMarkerSize(0.5);
    
    gr_sphi_r->Add(gr_sphir19);
    gr_sphi_r->Add(gr_sphir15);
    gr_sphi_r->Add(gr_sphir10);
    gr_sphi_r->Add(gr_sphir5);
    gr_sphi_r->Add(gr_sphir0);
    gr_sphi_r->Add(gr_sphirneg5);
    gr_sphi_r->Add(gr_sphirneg10);
    gr_sphi_r->Add(gr_sphirneg15);
    gr_sphi_r->Add(gr_sphirneg19);
    
    gr_sphi_r->Fit("f1");
    
    gr_sphi_r->SetTitle("#sigma_{#phi}_{d} vs r");
    gr_sphi_r->Draw("AP");
    gr_sphi_r->GetXaxis()->SetTitle("r [cm]");
    gr_sphi_r->GetYaxis()->SetTitle("#sigma_{#phi}_{d} [rad]");
    gr_sphi_r->SaveAs("figs/b_field_ana/sphi_r.png");
    
    TLegend *leg3 = new TLegend(0.8,0.6,0.9,0.9);
    leg3->AddEntry(gr_phir19,"z = 19 cm","lep");
    leg3->AddEntry(gr_phir15,"z = 15 cm","lep");
    leg3->AddEntry(gr_phir10,"z = 10 cm","lep");
    leg3->AddEntry(gr_phir5,"z = 5 cm","lep");
    leg3->AddEntry(gr_phir0,"z = 0 cm","lep");
    leg3->AddEntry(gr_phirneg5,"z = -5 cm","lep");
    leg3->AddEntry(gr_phirneg10,"z = -10 cm","lep");
    leg3->AddEntry(gr_phirneg15,"z = -15 cm","lep");
    leg3->AddEntry(gr_phirneg19,"z = -19 cm","lep");
    leg3->Draw();
    
    // time vs r
    c5->cd();
    
    gr_str19->SetMarkerColor(1);
    gr_str19->SetMarkerStyle(21);
    gr_str19->SetMarkerSize(0.5);
    gr_str15->SetMarkerColor(2);
    gr_str15->SetMarkerStyle(21);
    gr_str15->SetMarkerSize(0.5);
    gr_str10->SetMarkerColor(3);
    gr_str10->SetMarkerStyle(21);
    gr_str10->SetMarkerSize(0.5);
    gr_str5->SetMarkerColor(4);
    gr_str5->SetMarkerStyle(21);
    gr_str5->SetMarkerSize(0.5);
    gr_str0->SetMarkerColor(5);
    gr_str0->SetMarkerStyle(21);
    gr_str0->SetMarkerSize(0.5);
    gr_strneg5->SetMarkerColor(6);
    gr_strneg5->SetMarkerStyle(21);
    gr_strneg5->SetMarkerSize(0.5);
    gr_strneg10->SetMarkerColor(7);
    gr_strneg10->SetMarkerStyle(21);
    gr_strneg10->SetMarkerSize(0.5);
    gr_strneg15->SetMarkerColor(8);
    gr_strneg15->SetMarkerStyle(21);
    gr_strneg15->SetMarkerSize(0.5);
    gr_strneg19->SetMarkerColor(9);
    gr_strneg19->SetMarkerStyle(21);
    gr_strneg19->SetMarkerSize(0.5);
   
    gr_st_r->Add(gr_str19);
    gr_st_r->Add(gr_str15);
    gr_st_r->Add(gr_str10);
    gr_st_r->Add(gr_str5);
    gr_st_r->Add(gr_str0);
    gr_st_r->Add(gr_strneg5);
    gr_st_r->Add(gr_strneg10);
    gr_st_r->Add(gr_strneg15);
    gr_st_r->Add(gr_strneg19);
    
    gr_st_r->Fit("f1");
    
    gr_st_r->SetTitle("#sigma_{t}_{d} vs r");
    gr_st_r->Draw("AP");
    gr_st_r->GetXaxis()->SetTitle("r [cm]");
    gr_st_r->GetYaxis()->SetTitle("#sigma_{t}_{d} [ns]");
    gr_st_r->GetYaxis()->SetTitleOffset(1.3);
    gr_st_r->SaveAs("figs/b_field_ana/st_r.png");
    
    TLegend *leg5 = new TLegend(0.8,0.6,0.9,0.9);
    leg5->AddEntry(gr_tr19,"z = 19 cm","lep");
    leg5->AddEntry(gr_tr15,"z = 15 cm","lep");
    leg5->AddEntry(gr_tr10,"z = 10 cm","lep");
    leg5->AddEntry(gr_tr5,"z = 5 cm","lep");
    leg5->AddEntry(gr_tr0,"z = 0 cm","lep");
    leg5->AddEntry(gr_trneg5,"z = -5 cm","lep");
    leg5->AddEntry(gr_trneg10,"z = -10 cm","lep");
    leg5->AddEntry(gr_trneg15,"z = -15 cm","lep");
    leg5->AddEntry(gr_trneg19,"z = -19 cm","lep");
    leg5->Draw();

    // drift time vs z
    c6->cd();
    
    gr_td_z->SetMarkerColor(1);
    gr_td_z->SetMarkerStyle(21);
    gr_td_z->SetMarkerSize(0.5);
    
    gr_td_z->SetTitle("t_{d} vs z");
    gr_td_z->Draw("AP");
    gr_td_z->GetXaxis()->SetTitle("z [cm]");
    gr_td_z->GetYaxis()->SetTitle("t_{d} [ns]");
    gr_td_z->GetYaxis()->SetTitleOffset(1.3);
    gr_td_z->SaveAs("figs/b_field_ana/td_z.png");
    
    // drift time vs z
    c7->cd();
    
    gr_phid_z->SetMarkerColor(1);
    gr_phid_z->SetMarkerStyle(21);
    gr_phid_z->SetMarkerSize(0.5);
    
    gr_phid_z->SetTitle("#phi_{d} vs z");
    gr_phid_z->Draw("AP");
    gr_phid_z->GetXaxis()->SetTitle("z [cm]");
    gr_phid_z->GetYaxis()->SetTitle("#phi_{d} [rad]");
    gr_phid_z->SaveAs("figs/b_field_ana/phid_z.png");

    return;
    
}
