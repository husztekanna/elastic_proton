#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <iostream>

void ptotal_with_califa_hits()
{
    // Set nice drawing defaults
    gStyle->SetOptStat(0);
    gStyle->SetTextFont(42);

    // 1. Open the input file
    TFile *inFile = TFile::Open("analysis_output_califa_gt15MeV.root", "READ");
    if (!inFile || inFile->IsZombie())
    {
        std::cerr << "[ERROR] Could not open input file analysis_output_califa_gt15MeV.root\n";
        return;
    }

    // 2. Retrieve histograms created in analyse_all
    TH1F *h_all   = (TH1F*) inFile->Get("h1_frag_ptotal_lab_0");
    TH1F *h_1hit  = (TH1F*) inFile->Get("h1_frag_ptotal_lab_1hit_0");
    TH1F *h_2hit  = (TH1F*) inFile->Get("h1_frag_ptotal_lab_2hit_0");
    TH1F *h_3hit  = (TH1F*) inFile->Get("h1_frag_ptotal_lab_3hit_0");

    if (!h_all || !h_1hit || !h_2hit || !h_3hit)
    {
        std::cerr << "[ERROR] Could not retrieve one or more histograms from file.\n";
        std::cerr << "Check histogram names inside analysis_output_califa_gt15MeV.root!\n";
        inFile->Close();
        return;
    }

    // 3. Customize Visual Appearance (Colors & Styles)
    h_all->SetLineColor(kBlack);
    h_all->SetLineWidth(3);
    h_all->SetTitle("Fragment Total Momentum in Lab Frame;p_{total} (MeV/c);Counts");

    h_1hit->SetLineColor(kBlue + 1);
    h_1hit->SetLineWidth(2);

    h_2hit->SetLineColor(kRed + 1);
    h_2hit->SetLineWidth(2);

    h_3hit->SetLineColor(kGreen + 2);
    h_3hit->SetLineWidth(2);

    // 4. Create Canvas & Draw
    TCanvas *c1 = new TCanvas("c_ptotal", "Fragment Total Momentum Comparison", 800, 600);
    c1->SetGrid();

    // Determine max Y scale across all histograms so nothing gets clipped
    double max_y = h_all->GetMaximum();
    h_all->SetMaximum(max_y * 1.15); // Add 15% head room for legend

    h_1hit->Scale(h_all->Integral()/h_1hit->Integral());
    h_2hit->Scale(h_all->Integral()/h_2hit->Integral());
    h_3hit->Scale(h_all->Integral()/h_3hit->Integral());
    h_all->Draw("HIST");
    h_1hit->Draw("HIST SAME");
    h_2hit->Draw("HIST SAME");
    h_3hit->Draw("HIST SAME");

    // 5. Add Legend
    TLegend *legend = new TLegend(0.62, 0.68, 0.88, 0.88);
    legend->SetBorderSize(1);
    legend->SetFillColor(kWhite);
    legend->SetTextSize(0.035);
    legend->AddEntry(h_all,  "All Events",           "l");
    legend->AddEntry(h_1hit, "CALIFA 1 Hit",         "l");
    legend->AddEntry(h_2hit, "CALIFA 2 Hits",        "l");
    legend->AddEntry(h_3hit, "CALIFA 3 Hits",        "l");
    legend->Draw();

    // Update canvas
    c1->Update();

    // 6. Save Canvas and Histograms to output ROOT file
    TFile *outFile = new TFile("ptotal_with_califa_hits.root", "RECREATE");
    c1->Write();
    h_all->Write("h1_frag_ptotal_lab_all");
    h_1hit->Write("h1_frag_ptotal_lab_1hit");
    h_2hit->Write("h1_frag_ptotal_lab_2hit");
    h_3hit->Write("h1_frag_ptotal_lab_3hit");
    
    outFile->Close();
    inFile->Close();

    std::cout << "\n[INFO] Macro finished successfully. Canvas & histograms saved to ptotal_with_califa_hits.root\n";
}