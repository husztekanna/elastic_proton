#include <TCanvas.h>
#include <TString.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TPaveStats.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>


// Helper to apply canvas margins
void ApplyCanvasStyle(TCanvas *c) {
    c->SetMargin(0.15, 0.15, 0.15, 0.05); // [left, right, bottom, top]
}

// Helper to format 1D and 2D histogram axes uniformly
void FormatHistogram(TH1 *h) {
    h->SetTitle(""); // Remove top title
    
    // X-Axis
    //h->GetXaxis()->SetNdivisions(505, kTRUE); 
    h->GetXaxis()->SetTitleSize(0.04);
    h->GetXaxis()->SetLabelSize(0.04);
    h->GetXaxis()->SetTitleOffset(1.1);
    h->GetXaxis()->SetTitleFont(62); // Bold
    
    // Y-Axis
    //h->GetYaxis()->SetNdivisions(505, kTRUE);
    h->GetYaxis()->SetTitleSize(0.04);
    h->GetYaxis()->SetLabelSize(0.04);
    h->GetYaxis()->SetTitleOffset(1.1);
    h->GetYaxis()->SetTitleFont(62); // Bold
}

void applying_cuts() {

  // Data points for T2 (x) and CALIFA energy (y) in MeV
  const int nPoints = 17;
  double x_TKE[nPoints] = {31, 35, 40, 45, 50, 55, 60, 65, 70, 80, 90, 100, 120, 140, 160, 180, 200};
  double y_remaining[nPoints] = {4.214,  14.918, 23.393,  30.556,  37.103,  43.276, 49.209, 54.975, 60.614,
                                 71.633, 82.411, 93.026, 113.943, 134.596, 155.087, 175.469, 195.776};

  TGraph *g_calc = new TGraph(nPoints, x_TKE, y_remaining);
  g_calc->SetName("g_calc_T2_vs_E");
  g_calc->SetMarkerStyle(20);
  g_calc->SetMarkerSize(1);
  g_calc->SetMarkerColor(kBlack);
  g_calc->SetLineColor(kBlack);
  g_calc->SetLineWidth(2);

  gStyle->SetOptStat(0);


  // 1. Setup
  const int nFiles = 3;
  TString fragments[nFiles] = {"22O", "24F", "25F"};
  int colors[nFiles] = {kBlue + 1, kRed + 1, kGreen + 2};

  bool applyCuts = true; // Set to true to apply cuts, false to skip cuts

  TString thetaCut =
    "(theta_califa >= 70 && theta_califa <= 90)";

TString vertexCut =
    "(vertex_z >= -100.0 && vertex_z <= 200.0)";

TString phiCut =
    "(fabs(phi_califa - phi_proton) <= 14)";

TString trackCut =
    "(fabs(TX_diff) <= 0.0015 && fabs(TY_diff) <= 0.005)";

TString exclusionCut = ""
    "(!(phi_califa <= -55 && phi_califa >= -125))";

TString cuts = "";
if (applyCuts) {
  cuts = thetaCut + " && " +
         vertexCut + " && " +
         phiCut + " && " +
         trackCut  + " && " +
         exclusionCut;
}

  TCanvas *c1 = new TCanvas("c1", "Centre of mass angle Comparison", 800, 600);
  ApplyCanvasStyle(c1);
  c1->SetGrid();
  c1->SetLogy();

  TCanvas *c2 = new TCanvas("c2", "Mandelstam t comparison", 800, 600);
  ApplyCanvasStyle(c2);
  c2->SetGrid();
  c2->SetLogy();

  TCanvas *c3[nFiles];


  TLegend *legend = new TLegend(0.68, 0.72, 0.88, 0.88);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.06); // <-- Add this line (adjust between 0.04 and 0.06 as needed)
  legend->SetTextFont(62);    // <-- Optional: makes the text bold

  TH1F *theta_cm_histograms[nFiles];
  TH2F *px_py_histograms[nFiles];
  TH1F *vertex_z[nFiles];
  TH1F *frag_ptr[nFiles];
  TH2F *phi_califa_vs_phi_proton[nFiles];
  TH2F *TX_vs_TY_diff[nFiles];
  TH1F *mandelstam_t[nFiles];
  TH2F *t_error_vs_theta_scattering[nFiles];
  TH2F *califa_energy_vs_T2[nFiles];
  TH2F *vertex_x_y[nFiles];
  TH1F *phi_califa_minus_phi_proton[nFiles];
  TCanvas *c_phi_cal_vs_phi_p[nFiles];

  // 2. Process Files
  for (int i = 0; i < nFiles; i++) {

    TFile *file = TFile::Open(Form("analysis_output_%s.root", fragments[i].Data()), "READ");
    if (!file || file->IsZombie()) {
      Form("Error opening file for %s", fragments[i].Data());
      continue;
    }
    TTree *tree = (TTree *)file->Get("tree");

    TString theta_cm_histName = Form("h_theta_cm_%s", fragments[i].Data());
    TString px_py_histName = Form("h2_px_py_%s", fragments[i].Data());
    TString vertex_z_histName = Form("vertex_z_%s", fragments[i].Data());
    TString frag_ptr_histName = Form("frag_ptr_%s", fragments[i].Data());
    TString phi_califa_vs_phi_proton_histName = Form("phi_califa_vs_phi_proton_%s", fragments[i].Data());
    TString TX_vs_TY_diff_histName = Form("TX_vs_TY_diff_%s", fragments[i].Data());
    TString mandelstam_t_histName = Form("mandelstam_t_%s", fragments[i].Data());
    TString t_error_vs_theta_scattering_histName = Form("t_error_vs_theta_scattering_%s", fragments[i].Data());
    TString califa_energy_vs_T2_histName = Form("califa_energy_vs_T2_%s", fragments[i].Data());
    TString vertex_x_y_histName = Form("vertex_x_y_%s", fragments[i].Data());
    TString phi_califa_minus_phi_proton_histName = Form("phi_califa_minus_phi_proton_%s", fragments[i].Data());

    gROOT->cd(); // Ensure histograms belong to global ROOT directory

    theta_cm_histograms[i] = new TH1F(theta_cm_histName, "Center of Mass Angle Distribution; #theta_{cm} [deg]; Normalized Counts", 30, 15, 30);
    FormatHistogram(theta_cm_histograms[i]);
    theta_cm_histograms[i]->Sumw2();
    px_py_histograms[i] = new TH2F(px_py_histName, "; p_{x} [MeV/c]; p_{y} [MeV/c]", 200, -600, 600, 200, -600, 600);
    FormatHistogram(px_py_histograms[i]);
    vertex_z[i] = new TH1F(vertex_z_histName, Form("%s: Vertex-z; Vertex-z [mm]; Counts", fragments[i].Data()), 300, -300, 300);
    FormatHistogram(vertex_z[i]);
    frag_ptr[i] = new TH1F(frag_ptr_histName, Form("%s: frag_ptr [MeV/c]", fragments[i].Data()), 400, 0, 1500);
    FormatHistogram(frag_ptr[i]);
    phi_califa_vs_phi_proton[i] = new TH2F(phi_califa_vs_phi_proton_histName, "; #Phi of the nucleus [deg]; #Phi of proton in CALIFA [deg]", 360, -180, 180, 360, -180, 180);
    FormatHistogram(phi_califa_vs_phi_proton[i]);
    TX_vs_TY_diff[i] = new TH2F(TX_vs_TY_diff_histName, Form("%s: TX vs TY diff", fragments[i].Data()), 1000, -0.08, 0.08, 1000, -0.08, 0.08);
    FormatHistogram(TX_vs_TY_diff[i]);
    t_error_vs_theta_scattering[i] = new TH2F(t_error_vs_theta_scattering_histName, Form("%s: Error of mandelstam t vs theta_{scattering}", fragments[i].Data()), 180, 0, 10, 100, 0, 0.08);
    FormatHistogram(t_error_vs_theta_scattering[i]);
    mandelstam_t[i] = new TH1F(mandelstam_t_histName, "; t [GeV^{2}]; Normalized Counts", 30, 0.1, 0.4);
    FormatHistogram(mandelstam_t[i]);
    mandelstam_t[i]->Sumw2();
    califa_energy_vs_T2[i] = new TH2F(califa_energy_vs_T2_histName, "; Calculated energy of proton T_{R} [MeV]; Energy of high energy cluster in CALIFA [MeV]", 200, 0, 200, 200, 0, 200);
    FormatHistogram(califa_energy_vs_T2[i]);
    vertex_x_y[i] = new TH2F(vertex_x_y_histName, "; Vertex-x [mm]; Vertex-y [mm]", 200, -30, 30, 200, -30, 30);
    FormatHistogram(vertex_x_y[i]);
    phi_califa_minus_phi_proton[i] = new TH1F(phi_califa_minus_phi_proton_histName, Form("; #Phi_{nucleus} - #Phi_{proton} [deg]; Counts"), 300, -150, 150);
    FormatHistogram(phi_califa_minus_phi_proton[i]);

    // Fill histograms from tree branches
    TString thetaWeight = "(1.0 / sin(theta_cm * TMath::DegToRad()))";
    TString validityCut = "(theta_cm == theta_cm && theta_cm >= 0 && theta_cm <= 180 && TMath::Abs(sin(theta_cm*TMath::DegToRad())) > 1e-6)";
    TString thetaSelection = validityCut;
    if (applyCuts) {
      thetaSelection = Form("(%s) && (%s)", cuts.Data(), validityCut.Data());
    }
    TString weightedSelection = Form("(%s) ? (%s) : 0",
                                     thetaSelection.Data(), thetaWeight.Data());
    TString thetaExpression = Form("(%s) ? theta_cm : 0", validityCut.Data());

    tree->Draw(Form("%s >> %s", thetaExpression.Data(), theta_cm_histName.Data()),
           weightedSelection.Data(),
           "goff");
    tree->Draw(Form("frag_py:frag_px >> %s", px_py_histName.Data()), cuts.Data(), "goff");
    tree->Draw(Form("vertex_z >> %s", vertex_z_histName.Data()), cuts.Data(), "goff");
    tree->Draw(Form("frag_ptr >> %s", frag_ptr_histName.Data()), cuts.Data(), "goff");
    tree->Draw(Form("phi_califa:phi_proton >> %s", phi_califa_vs_phi_proton_histName.Data()), cuts.Data(), "goff");
    tree->Draw(Form("TX_diff:TY_diff >> %s", TX_vs_TY_diff_histName.Data()), cuts.Data(), "goff");
    tree->Draw(Form("- mandelstam_t >> %s", mandelstam_t_histName.Data()), cuts.Data(), "goff");
    tree->Draw(Form("err_mandelstam_t : theta_scattering * TMath::RadToDeg() >> %s", t_error_vs_theta_scattering_histName.Data()), cuts.Data(), "goff");
    tree->Draw(Form("califa_energy:T2 >> %s", califa_energy_vs_T2_histName.Data()), cuts.Data(), "goff");
    tree->Draw(Form("vertex_y : vertex_x >> %s", vertex_x_y_histName.Data()), cuts.Data(), "goff");
    tree->Draw(Form("phi_califa - phi_proton >> %s", phi_califa_minus_phi_proton_histName.Data()), cuts.Data(), "goff");

    file->Close();
    delete file;
  }

  // 3. Normalize & Calculate Peak Headroom
  double theta_cm_max = 0.0;
  double mandelstam_t_max = 0.0;

  for (int i = 0; i < nFiles; i++) {
    int bin_min = 1;
    int bin_max = theta_cm_histograms[i]->GetNbinsX();
    double integral_range = theta_cm_histograms[i]->Integral(bin_min, bin_max);
    if (integral_range > 0) {
      theta_cm_histograms[i]->Scale(1.0 / integral_range);
    }

    if (mandelstam_t[i]->Integral() > 0) {
      mandelstam_t[i]->Scale(1.0 / mandelstam_t[i]->Integral());
    }

    if (theta_cm_histograms[i]->GetMaximum() > theta_cm_max) {
      theta_cm_max = theta_cm_histograms[i]->GetMaximum();
    }

    if (mandelstam_t[i]->GetMaximum() > mandelstam_t_max) {
      mandelstam_t_max = mandelstam_t[i]->GetMaximum();
    }

    legend->AddEntry(theta_cm_histograms[i], fragments[i], "l");
  }

  // 4. Draw Plots on Canvases
  c1->cd();

  for (int i = 0; i < nFiles; i++) {
    theta_cm_histograms[i]->SetMinimum(1e-4);
    theta_cm_histograms[i]->SetMaximum(theta_cm_max * 2.0);
    TString option = (i == 0) ? "E1 P" : "E1 P SAMES";
    theta_cm_histograms[i]->SetLineColor(colors[i]);
    theta_cm_histograms[i]->SetLineWidth(2);
    theta_cm_histograms[i]->Draw(option);
  }
  legend->Draw();
  c1->Modified();
  c1->Update();
  c1->SaveAs("theta_cm_comparison.C");

  c2->cd();

  for (int i = 0; i < nFiles; i++) {
    mandelstam_t[i]->SetMaximum(mandelstam_t_max * 1.25);
    TString option = (i == 0) ? "E1 P" : "E1 P SAMES";
    mandelstam_t[i]->SetLineColor(colors[i]);
    mandelstam_t[i]->SetLineWidth(2);
    mandelstam_t[i]->Draw(option);
  }
  legend->Draw();
  c2->Modified();
  c2->Update();
  c2->SaveAs("mandelstam_t_comparison.C");

  // Draw per-isotope individual plots
  for (int i = 0; i < nFiles; i++) {
    TCanvas *c_pxpy = new TCanvas(Form("c_pxpy_%s", fragments[i].Data()), "", 800, 600);
    ApplyCanvasStyle(c_pxpy);
    c_pxpy->SetGrid();
    px_py_histograms[i]->Draw("COLZ");
    c_pxpy->SaveAs(Form("px_py_%s.png", fragments[i].Data()));

    TCanvas *c_xy = new TCanvas(Form("c_xy_%s", fragments[i].Data()), "", 800, 600);
    ApplyCanvasStyle(c_xy);
    c_xy->SetGrid();
    vertex_x_y[i]->Draw("COLZ");
    c_xy->SaveAs(Form("vertex_x_y_%s.png", fragments[i].Data()));

    TCanvas *c_z = new TCanvas(Form("c_z_%s", fragments[i].Data()), "", 800, 600);
    ApplyCanvasStyle(c_z);
    c_z->SetGrid();
    vertex_z[i]->Draw(); // Fixed: 1D Histogram drawing, removed "COLZ"
    c_z->SaveAs(Form("vertex_z_%s.png", fragments[i].Data()));

    // 2D Histogram: Phi CALIFA vs Phi Proton (Log Z)
    c_phi_cal_vs_phi_p[i] = new TCanvas(Form("phi_cal_vs_phi_p_%s", fragments[i].Data()), "", 800, 600);
    ApplyCanvasStyle(c_phi_cal_vs_phi_p[i]);
    c_phi_cal_vs_phi_p[i]->cd();
    c_phi_cal_vs_phi_p[i]->SetGrid();
    
    gPad->SetLogz(); // Set log z on the active canvas pad
    phi_califa_vs_phi_proton[i]->SetMinimum(1); // Set min positive value so 0-count bins don't break log scale
    phi_califa_vs_phi_proton[i]->Draw("COLZ");
    
    c_phi_cal_vs_phi_p[i]->Update();
    c_phi_cal_vs_phi_p[i]->SaveAs(Form("phi_califa_vs_phi_proton_%s.png", fragments[i].Data()));

    // 1D Histogram: Phi CALIFA minus Phi Proton
    TCanvas *c_phi_cal_minus_phi_p = new TCanvas(Form("phi_cal_minus_phi_p_%s", fragments[i].Data()), "", 800, 600);
    ApplyCanvasStyle(c_phi_cal_minus_phi_p);
    c_phi_cal_minus_phi_p->cd();
    c_phi_cal_minus_phi_p->SetGrid();
    
    phi_califa_minus_phi_proton[i]->Draw();
    c_phi_cal_minus_phi_p->SaveAs(Form("phi_califa_minus_phi_proton_%s.png", fragments[i].Data()));    

    TCanvas *c_mandelstam_t = new TCanvas(Form("mandelstam_t_%s", fragments[i].Data()), "", 800, 600);
    ApplyCanvasStyle(c_mandelstam_t);
    c_mandelstam_t->SetGrid();
    mandelstam_t[i]->Draw(); // Fixed: 1D Histogram drawing, removed "COLZ"
    c_mandelstam_t->SaveAs(Form("mandelstam_t_%s.png", fragments[i].Data()));

    TCanvas *c_theta_cm = new TCanvas(Form("theta_cm_%s", fragments[i].Data()), "", 800, 600);
    ApplyCanvasStyle(c_theta_cm);
    c_theta_cm->SetGrid();
    theta_cm_histograms[i]->Draw(); // Fixed: 1D Histogram drawing, removed "COLZ"
    c_theta_cm->SaveAs(Form("theta_cm_%s.png", fragments[i].Data()));

    TCanvas *c_vertex_z = new TCanvas(Form("vertex_z_%s", fragments[i].Data()), "", 800, 600);
    ApplyCanvasStyle(c_vertex_z);
    c_vertex_z->SetGrid();
    vertex_z[i]->Draw(); // Fixed: 1D Histogram drawing, removed "COLZ"
    c_vertex_z->SaveAs(Form("vertex_z_%s.png", fragments[i].Data()));

    TCanvas *c_vertex_x_y = new TCanvas(Form("vertex_x_y_%s", fragments[i].Data()), "", 800, 600);
    ApplyCanvasStyle(c_vertex_x_y);
    c_vertex_x_y->SetGrid();
    vertex_x_y[i]->Draw(); // Fixed: 1D Histogram drawing, removed "COLZ"
    c_vertex_x_y->SaveAs(Form("vertex_x_y_%s.png", fragments[i].Data()));

    c3[i] = new TCanvas(Form("c3_%s", fragments[i].Data()), Form("CALIFA Energy vs T2 (%s)", fragments[i].Data()), 800, 600);
    ApplyCanvasStyle(c3[i]); // Apply margin settings (0.14, 0.04, 0.14, 0.04)
    c3[i]->cd();
    c3[i]->SetGrid();

    // --- Legend styling ---
    TLegend *leg3 = new TLegend(0.18, 0.75, 0.48, 0.88);
    leg3->SetBorderSize(1);
    leg3->SetFillColor(kWhite);
    leg3->AddEntry(g_calc, "Calculated Values", "lp");

    califa_energy_vs_T2[i]->SetTitle(""); // Remove top title
    califa_energy_vs_T2[i]->GetXaxis()->SetTitle("Calculated energy of proton T_{R} [MeV]");
    califa_energy_vs_T2[i]->GetYaxis()->SetTitle("Energy of high energy cluster in CALIFA [MeV]");

    califa_energy_vs_T2[i]->Draw("COLZ");
    g_calc->SetMarkerColor(kRed);
    g_calc->SetLineColor(kRed);
    g_calc->Draw("LP SAME");
    leg3->Draw();
    c3[i]->Update();

    c3[i]->SaveAs(Form("califa_vs_t2_%s.png", fragments[i].Data()));
    
  }

  // 5. Save Output File
  TFile *outFile;

  if (applyCuts) {
    outFile = new TFile("plot_with_cuts.root", "RECREATE");
  } else {
    outFile = new TFile("plot_without_cuts.root", "RECREATE");
  }
  for (int i = 0; i < nFiles; i++) {
    theta_cm_histograms[i]->Write();
    px_py_histograms[i]->Write();
    vertex_z[i]->Write();
    frag_ptr[i]->Write();
    phi_califa_vs_phi_proton[i]->Write();
    c_phi_cal_vs_phi_p[i]->Write();
    TX_vs_TY_diff[i]->Write();
    mandelstam_t[i]->Write();
    t_error_vs_theta_scattering[i]->Write();
    califa_energy_vs_T2[i]->Write();
    vertex_x_y[i]->Write();
    phi_califa_minus_phi_proton[i]->Write();
  }

  c1->Write();
  c2->Write();
  g_calc->Write();
  for (int i = 0; i < nFiles; i++) {
    c3[i]->Write();
  }
  outFile->Close();
}
