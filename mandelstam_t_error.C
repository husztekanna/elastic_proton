
// Plots the analytic uncertainty on Mandelstam t as a function of scattering 
// angle theta_scattering, propagated from a fixed angular resolution 
// (error_theta_scattering = 0.0004 rad), independent of any data file.
//
// Formula matches the err_mandelstam_t calculation used in 
// elastic_proton_scattering.C.
//
// Output: mandelstam_t_error.png, mandelstam_t_error.root (TF1 + TCanvas)
// Usage:  root -l mandelstam_t_error.C
//
// Note: constants here (mp, p_inc, e_inc) are hardcoded in GeV, whereas the 
// same quantities in elastic_proton_scattering.C are computed per-event in 
// MeV — check these match the specific beam/fragment case you're 
// illustrating before reusing this for a different energy.

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TStyle.h>

void mandelstam_t_error() {
    gStyle->SetOptStat(0);
    
    TCanvas *c1 = new TCanvas("c1", "Mandelstam t Error Function", 800, 600);
    c1->SetGrid();
    // Give the axis titles room so they aren't clipped
    c1->SetLeftMargin(0.15);
    c1->SetBottomMargin(0.15);
    c1->SetRightMargin(0.05);
    c1->SetTopMargin(0.08);

    // 1. Define Constants
    double error_theta_scattering = 0.0004; // [rad]
    double mp = 0.93827;  // GeV
    double p_inc = 31.0;  // GeV/c
    double e_inc = 37.0;  // GeV
    
    double d2r = 0.01745329251; // deg to rad factor (pi / 180)

    // 2. Build Formula
    TString formula = Form("( (2 * pow(%f, 2) * pow(%f, 2) * sin(x * %f)) / pow(%f + (2 * %f * pow(sin((x * %f) / 2.0), 2)), 2) ) * %f",
                           mp, p_inc, d2r, mp, e_inc, d2r, error_theta_scattering);

    double xmin_deg = 0.0;
    double xmax_deg = 3.0;
    
    TF1 *f_t_err = new TF1("f_t_err", formula.Data(), xmin_deg, xmax_deg);
    f_t_err->SetNpx(1000);
    f_t_err->SetTitle("");
    f_t_err->SetLineColor(kOrange + 1);
    f_t_err->SetLineWidth(2);
    f_t_err->GetXaxis()->SetTitleSize(0.06);
    f_t_err->GetYaxis()->SetTitleSize(0.06);
    f_t_err->GetXaxis()->SetTitleFont(62);
    f_t_err->GetYaxis()->SetTitleFont(62);
    f_t_err->GetXaxis()->SetTitle("#theta_{scattering} [deg]");
    f_t_err->GetYaxis()->SetTitle("Uncertainty of t [GeV^{2}/c^{2}]");
    f_t_err->GetXaxis()->SetLabelSize(0.05);
    f_t_err->GetYaxis()->SetLabelSize(0.05);
    f_t_err->GetXaxis()->SetTitleOffset(1.1);
    f_t_err->GetYaxis()->SetTitleOffset(1.2);

    c1->cd();
    f_t_err->Draw();
    c1->Modified();
    c1->Update();

    // Save PNG
    c1->SaveAs("mandelstam_t_error.png");

    // Save ROOT file
    TFile *outFile = new TFile("mandelstam_t_error.root", "RECREATE");
    f_t_err->Write();
    c1->Write();
    outFile->Close();
}
