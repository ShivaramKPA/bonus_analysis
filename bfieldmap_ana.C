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

    float zneg19_phi[9] = {0.979244,
        0.834884,
        0.678226,
        0.535752,
        0.410171,
        0.295019,
        0.188652,
        0.0910025,
        0.0177088};
    
    float zneg19_t[9] = {4786.47,
            4466.1,
            4018.76,
            3494.36,
            2925.06,
            2287.43,
            1586.39,
            822.569,
            169.395};
    
    float zneg15_phi[9] = {0.988219,
        0.843082,
        0.683793,
        0.541156,
        0.41388,
        0.297647,
        0.190756,
        0.0918548,
        0.0178318};
    
    float zneg15_t[9] = {4841.58,
            4518.43,
            4060.33,
            3531.92,
            2954.23,
            2309.71,
            1603.51,
            831.784,
            171.887};
    
    float zneg10_phi[9] = {0.994556,
        0.849106,
        0.68795,
        0.545733,
        0.416377,
        0.29937,
        0.19125,
        0.0923999,
        0.0180024};
    
    float zneg10_t[9] = {4875.19,
            4551.68,
            4087.48,
            3562.46,
            2976.59,
            2328.67,
            1612.73,
            836.586,
            173.9};
    
    float zneg5_phi[9] = {0.996141,
        0.850608,
        0.689877,
        0.546082,
        0.41678,
        0.299997,
        0.192103,
        0.0927311,
        0.0179201};
    
    float zneg5_t[9] = {4886.21,
            4562.93,
            4097.36,
            3573.92,
            2981.41,
            2332.51,
            1618.15,
            840.089,
            172.559};
    
    float z0_phi[9] = {0.996211,
        0.852161,
        0.688313,
        0.546201,
        0.417339,
        0.299948,
        0.191907,
        0.0927202,
        0.0179773};
    
    float z0_t[9] = {4885.6,
            4566.44,
            4091.69,
            3571.12,
            2985.53,
            2333.04,
            1616.46,
            841.034,
            173.272};
    
    float z5_phi[9] = {0.995114,
        0.850148,
        0.68858,
        0.545378,
        0.416608,
        0.299867,
        0.192076,
        0.0924266,
        0.017865};
    
    float z5_t[9] = {4881.83,
            4554.19,
            4091.26,
            3563.12,
            2979.04,
            2327.22,
            1614.84,
            836.69,
            171.505};
    
    float z10_phi[9] = {0.992428,
        0.846899,
        0.687275,
        0.544144,
        0.415107,
        0.298674,
        0.19119,
        0.09237,
        0.0179606};
    
    float z10_t[9] = {4865.08,
            4538.23,
            4078.53,
            3554.98,
            2967.41,
            2320.97,
            1608.09,
            837.274,
            172.176};
    
    float z15_phi[9] = {0.985556,
        0.846806,
        0.681562,
        0.540151,
        0.41257,
        0.296584,
        0.19015,
        0.0918789,
        0.0178904};
    
    float z15_t[9] = {4824.36,
            4538.86,
            4044.33,
            3524.2,
            2943.63,
            2301.8,
            1598.26,
            830.662,
            171.514};
    
    float z19_phi[9] = {0.975401,
        0.833244,
        0.675335,
        0.534832,
        0.408592,
        0.29378,
        0.188092,
        0.0908892,
        0.0176751};
    
    float z19_t[9] = {4768.56,
            4449.01,
            3998.29,
            3484.68,
            2909.41,
            2276.69,
            1578.53,
            821.215,
            169.258};
    
    //errors
    float zneg19_phi_err[9] = {0.00790143,
        0.00654047,
        0.00577262,
        0.00491193,
        0.0051883,
        0.00348415,
        0.00261958,
        0.00197608,
        0.000892642};
    
    float zneg19_t_err[9] = {28.486,
        28.9564,
        24.5626,
        24.3323,
        31.3933,
        20.0094,
        16.6143,
        13.3087,
        6.46456};
    
    float zneg15_phi_err[9] = {0.00709689,
        0.00605749,
        0.00784265,
        0.00660392,
        0.00413568,
        0.00364627,
        0.00314986,
        0.00192226,
        0.000843152};
    
    float zneg15_t_err[9] = {28.1514,
        29.4616,
        37.5309,
        34.9881,
        24.1082,
        22.609,
        22.3217,
        13.4994,
        5.62323};
    
    float zneg10_phi_err[9] = {0.0075145,
        0.00661811,
        0.00589645,
        0.00507366,
        0.0044142,
        0.00352366,
        0.00290565,
        0.00207782,
        0.000954755};
    
    float zneg10_t_err[9] = {30.4804,
        26.9286,
        28.7042,
        25.3424,
        24.7189,
        22.2829,
        19.2258,
        14.3192,
        6.80384};
    
    float zneg5_phi_err[9] = {0.00726524,
        0.00653824,
        0.00587723,
        0.00513663,
        0.00483408,
        0.00424706,
        0.00316663,
        0.00215564,
        0.00087729};
    
    float zneg5_t_err[9] = {28.551,
        28.4971,
        25.8622,
        26.7617,
        24.0515,
        27.4591,
        21.7022,
        14.7628,
        6.44206};
    
    float z0_phi_err[9] = {0.0074328,
        0.00627894,
        0.00569887,
        0.00513311,
        0.0044647,
        0.0037362,
        0.00327845,
        0.0021453,
        0.000915396};
    
    float z0_t_err[9] = {28.7466,
        28.5049,
        22.7742,
        25.0495,
        24.0274,
        24.0352,
        22.4133,
        15.9336,
        6.704};

    float z5_phi_err[9] = {0.00758277,
        0.00723578,
        0.00567169,
        0.00532652,
        0.00412484,
        0.00399051,
        0.00286755,
        0.00203164,
        0.000904946};
    
    float z5_t_err[9] = {25.8842,
        30.1976,
        26.8211,
        24.7646,
        24.6906,
        22.3199,
        19.1085,
        12.8893,
        6.69598};
    
    float z10_phi_err[9] = {0.0079517,
        0.00678811,
        0.00564569,
        0.00514464,
        0.00520969,
        0.00367976,
        0.00339461,
        0.00209676,
        0.000925309};
    
    float z10_t_err[9] = {31.596,
        27.8471,
        25.7786,
        25.0546,
        29.7535,
        21.3839,
        23.0089,
        14.7259,
        6.90231};
    
    float z15_phi_err[9] = {0.00789608,
        0.00674795,
        0.00599379,
        0.00525774,
        0.00466732,
        0.00332019,
        0.00298609,
        0.00190694,
        0.0009323};
    
    float z15_t_err[9] = {30.3422,
        29.4611,
        26.3213,
        26.0193,
        25.8689,
        21.8045,
        20.0273,
        12.7769,
        6.9929};
    
    float z19_phi_err[9] = {0.00743052,
        0.00651113,
        0.00628804,
        0.005549,
        0.0047017,
        0.0040629,
        0.0029009,
        0.00207914,
        0.000840564};
    
    float z19_t_err[9] = {29.2302,
        29.5884,
        28.6999,
        26.778,
        22.4365,
        25.8536,
        18.7653,
        14.5105,
        6.05247};
    
    float zneg19_sphi_err[9] = {0.000351997,
        0.000239236,
        0.000290456,
        0.000187707,
        0.000154306,
        0.000136742,
        0.000106177,
        8.74E-05,
        3.33E-05};
    
    float zneg15_sphi_err[9] = {0.000306157,
        0.000265493,
        0.000230052,
        0.000191892,
        0.000178926,
        0.000128785,
        0.000106578,
        6.91E-05,
        3.23E-05};
    
    float zneg10_sphi_err[9] = {0.000303399,
        0.000225117,
        0.000221807,
        0.000209913,
        0.000178391,
        0.000145892,
        0.000123428,
        7.53E-05,
        3.37E-05};
    
    float zneg5_sphi_err[9] = {0.000298755,
        0.00025374,
        0.000236981,
        0.00020235,
        0.000186904,
        0.00013438,
        0.000105513,
        8.12E-05,
        3.76E-05};
    
    float z0_sphi_err[9] = {0.000299118,
        0.000254479,
        0.000204695,
        0.000203234,
        0.000179085,
        0.000157334,
        0.000113575,
        6.72E-05,
        3.51E-05};
    
    float z5_sphi_err[9] = {0.000274701,
        0.00027257,
        0.00024733,
        0.000235656,
        0.000153885,
        0.000168294,
        0.000118222,
        8.62E-05,
        3.60E-05};
    
    float z10_sphi_err[9] = {0.000317614,
        0.000268006,
        0.000252839,
        0.000204963,
        0.000171743,
        0.000142802,
        0.000112857,
        6.75E-05,
        2.80E-05};
    
    float z15_sphi_err[9] = {0.000384659,
        0.000274441,
        0.000244776,
        0.000207162,
        0.000189574,
        0.000142777,
        0.000118588,
        7.67E-05,
        3.29E-05};
    
    float z19_sphi_err[9] = {0.00030282,
        0.000240911,
        0.000259822,
        0.000248542,
        0.000197171,
        0.00014452,
        0.000117841,
        5.80E-05,
        3.63E-05};
    
    float zneg19_st_err[9] = {1.19279,
        1.1964,
        0.958631,
        0.952715,
        0.99854,
        0.797354,
        0.590523,
        0.543354,
        0.252611};
    
    float zneg15_st_err[9] = {1.13925,
        1.40455,
        2.35018,
        1.16836,
        0.945168,
        0.866186,
        0.719946,
        0.519423,
        0.260578};
    
    float zneg10_st_err[9] = {1.38117,
        1.02439,
        1.08798,
        1.02371,
        0.95184,
        0.824445,
        0.684061,
        0.522268,
        0.229139};
    
    float zneg5_st_err[9] = {1.20666,
        1.28276,
        0.954335,
        1.12086,
        0.970093,
        0.986289,
        0.741688,
        0.438246,
        0.233441};
    
    float z0_st_err[9] = {1.21126,
        1.15066,
        0.79847,
        0.960503,
        0.876151,
        0.906762,
        0.700953,
        0.438141,
        0.257788};
    
    float z5_st_err[9] = {1.05898,
        1.31267,
        1.12558,
        1.06837,
        0.994015,
        0.847037,
        0.904993,
        0.560166,
        0.254771};
    
    float z10_st_err[9] = {1.25408,
        1.04339,
        1.10079,
        1.00323,
        1.56836,
        0.854004,
        0.640883,
        0.465881,
        0.220096};
    
    float z15_st_err[9] = {1.24873,
        1.23779,
        0.978658,
        1.06969,
        0.998819,
        0.914399,
        0.84109,
        0.527136,
        0.293209};
    
    float z19_st_err[9] = {1.24207,
        1.21537,
        1.34443,
        1.04935,
        0.871928,
        0.729927,
        0.734899,
        0.439147,
        0.251086};
    
    
    float r[9] = {3.1, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 6.9};
    float r_err[9] = {0,0,0,0,0,0,0,0,0};
    
    
    //________________________________________________________________________________________________
    //_____________________________________ Canvas and histograms ____________________________________
    //________________________________________________________________________________________________
    

    TCanvas *c1 = new TCanvas("c1","t_d vs phi_d",800,600);
    TCanvas *c2 = new TCanvas("c2","phi_d vs r",800,600);
    TCanvas *c3 = new TCanvas("c3","t_d vs r",800,600);
    TCanvas *c4 = new TCanvas("c4","sphi_d vs r",800,600);
    TCanvas *c5= new TCanvas("c5","st_d vs r",800,600);

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

    return;
    
}
