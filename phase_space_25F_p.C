#include <TLorentzVector.h>
#include <TFile.h>
#include <TParameter.h>
#include <TGenPhaseSpace.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TROOT.h>
#include <iostream>
#include <TGraph.h>

void phase_space_25F_p() {
  gROOT->SetBatch(kTRUE);

  const Double_t m_p = 0.938272; // GeV
  const Double_t m_delta = 1.232;
  const Double_t exitation_energy1 = 100. / 1000.0;  // GeV
  const Double_t exitation_energy2 = 2.5 / 1000.0; // GeV
  const Double_t exitation_energy3 = 4. / 1000.0;  // GeV
  const Bool_t useDelta = kFALSE;

  Double_t m_out = useDelta ? m_delta : m_p;
  const char* outName = useDelta ? "delta" : "proton";
  const Double_t m_25F = 23.294089; // GeV

  const Double_t p_beam = 31.0;  // GeV/c
  TLorentzVector target(0.0, 0.0, 0.0, m_p);
  TLorentzVector beam(0.0, 0.0, p_beam, TMath::Sqrt(p_beam * p_beam + m_25F * m_25F));

  TLorentzVector W = beam + target;

  Double_t masses[2]     = {m_25F, m_out};
  Double_t inelastic1[2] = {m_25F + exitation_energy1, m_out};
  Double_t inelastic2[2] = {m_25F + exitation_energy2, m_out};
  Double_t inelastic3[2] = {m_25F + exitation_energy3, m_out};

  TGraph* g_ground = new TGraph(); g_ground->SetName("g_ground");
  TGraph* g_ex1    = new TGraph(); g_ex1->SetName("g_ex1");
  TGraph* g_ex2    = new TGraph(); g_ex2->SetName("g_ex2");
  TGraph* g_ex3    = new TGraph(); g_ex3->SetName("g_ex3");

  TGenPhaseSpace event;      event.SetDecay(W, 2, masses);
  TGenPhaseSpace event_ex1;  event_ex1.SetDecay(W, 2, inelastic1);
  TGenPhaseSpace event_ex2;  event_ex2.SetDecay(W, 2, inelastic2);
  TGenPhaseSpace event_ex3;  event_ex3.SetDecay(W, 2, inelastic3);

  TH2F* px_py_25F = new TH2F("px_py_25F", "px vs py of 25F; px (MeV/c); py (MeV/c)", 200, -1000.0, 1000.0, 200, -1000.0, 1000.0);

  // Histograms
  TH2F* phi_califa_vs_phi_proton_ground = new TH2F("phi_califa_vs_phi_proton_ground", "25F Ground State: Phi CALIFA vs Phi proton ", 360, -180, 180, 360, -180, 180);
  TH2F* phi_califa_vs_phi_proton_ex1    = new TH2F("phi_califa_vs_phi_proton_ex1",    "25F Ex1: Phi CALIFA vs Phi proton", 360, -180, 180, 360, -180, 180);
  TH2F* phi_califa_vs_phi_proton_ex2    = new TH2F("phi_califa_vs_phi_proton_ex2",    "25F Ex2: Phi CALIFA vs Phi proton", 360, -180, 180, 360, -180, 180);
  TH2F* phi_califa_vs_phi_proton_ex3    = new TH2F("phi_califa_vs_phi_proton_ex3",    "25F Ex3: Phi CALIFA vs Phi proton", 360, -180, 180, 360, -180, 180);
  

  for (Int_t i = 0; i < 1000000; ++i) {
    // --- Ground State ---
    Double_t weight = event.Generate();
    if (weight > 0.0) {

      TLorentzVector* p25F = event.GetDecay(0);
      TLorentzVector* pOut = event.GetDecay(1);

      Double_t phi_proton = pOut->Vect().Phi() * 180 / TMath::Pi();
      Double_t phi_califa = p25F->Vect().Phi() * 180 / TMath::Pi();
      Double_t px_25F = p25F->Px() * 1000.0;
      Double_t py_25F = p25F->Py() * 1000.0;
      Double_t theta_lab = pOut->Vect().Theta() * 180.0 / TMath::Pi();

        px_py_25F->Fill(px_25F, py_25F, weight);
        if (g_ground->GetN() < 5000) g_ground->SetPoint(g_ground->GetN(), phi_proton, phi_califa);
        phi_califa_vs_phi_proton_ground->Fill(phi_proton, phi_califa, weight);
        
    }

    // --- Excited State 1 ---
    Double_t weight_ex1 = event_ex1.Generate();
    if (weight_ex1 > 0.0) {
      TLorentzVector* p25F = event_ex1.GetDecay(0);
      TLorentzVector* pOut = event_ex1.GetDecay(1);

      Double_t phi_proton = pOut->Vect().Phi() * 180.0 / TMath::Pi();
      Double_t phi_califa = p25F->Vect().Phi() * 180.0 / TMath::Pi();

      if (g_ex1->GetN() < 5000) g_ex1->SetPoint(g_ex1->GetN(), phi_proton, phi_califa);
      phi_califa_vs_phi_proton_ex1->Fill(phi_proton, phi_califa, weight_ex1);
    }

    // --- Excited State 2 ---
    Double_t weight_ex2 = event_ex2.Generate();
    if (weight_ex2 > 0.0) {
      TLorentzVector* p25F = event_ex2.GetDecay(0);
      TLorentzVector* pOut = event_ex2.GetDecay(1);

      Double_t phi_proton = pOut->Vect().Phi() * 180.0 / TMath::Pi();
      Double_t phi_califa = p25F->Vect().Phi() * 180.0 / TMath::Pi();

      if (g_ex2->GetN() < 5000) g_ex2->SetPoint(g_ex2->GetN(), phi_proton, phi_califa);
      phi_califa_vs_phi_proton_ex2->Fill(phi_proton, phi_califa, weight_ex2);
    }

    // --- Excited State 3 ---
    Double_t weight_ex3 = event_ex3.Generate();
    if (weight_ex3 > 0.0) {
      TLorentzVector* p25F = event_ex3.GetDecay(0);
      TLorentzVector* pOut = event_ex3.GetDecay(1);

      Double_t phi_proton = pOut->Vect().Phi() * 180.0 / TMath::Pi();
      Double_t phi_califa = p25F->Vect().Phi() * 180.0 / TMath::Pi();

      if (g_ex3->GetN() < 5000) g_ex3->SetPoint(g_ex3->GetN(), phi_proton, phi_califa);
      phi_califa_vs_phi_proton_ex3->Fill(phi_proton, phi_califa, weight_ex3);
    }

    
  }

  TFile* out = new TFile(Form("phase_space_25F_p_%s.root", outName), "RECREATE");
  px_py_25F->Write();
  phi_califa_vs_phi_proton_ground->Write();
  phi_califa_vs_phi_proton_ex1->Write();
  phi_califa_vs_phi_proton_ex2->Write();
  phi_califa_vs_phi_proton_ex3->Write();

 // Set styling for TGraph scatter points
  g_ground->SetMarkerColor(kBlack);   g_ground->SetMarkerStyle(20); g_ground->SetMarkerSize(0.5);
  g_ex1->SetMarkerColor(kRed);       g_ex1->SetMarkerStyle(20); g_ex1->SetMarkerSize(0.5);
  g_ex2->SetMarkerColor(kBlue);      g_ex2->SetMarkerStyle(20); g_ex2->SetMarkerSize(0.5);
  g_ex3->SetMarkerColor(kGreen+2);   g_ex3->SetMarkerStyle(20); g_ex3->SetMarkerSize(0.5);

  TCanvas* c_overlay = new TCanvas("c_overlay", "Phi Correlation Comparison", 800, 700);

  // Draw empty TH2F frame to define axes
  TH2F* frame = new TH2F("frame", "Phi CALIFA vs Phi Proton Comparison; #phi_{proton} [deg]; #phi_{CALIFA} [deg]", 360, -180, 180, 360, -180, 180);
  frame->SetStats(kFALSE);
  frame->Draw();

  // Overlay graphs as scatter markers
  g_ground->Draw("P SAME");
  g_ex1->Draw("P SAME");
  g_ex2->Draw("P SAME");
  g_ex3->Draw("P SAME");

  // Legend
  TLegend* leg = new TLegend(0.15, 0.70, 0.40, 0.88);
  leg->AddEntry(g_ground, "Ground State", "p");
  leg->AddEntry(g_ex1,    "Ex1 (1.0 MeV)", "p");
  leg->AddEntry(g_ex2,    "Ex2 (2.5 MeV)", "p");
  leg->AddEntry(g_ex3,    "Ex3 (4.0 MeV)", "p");
  leg->Draw();

  c_overlay->Write();
  
  out->Close();
}