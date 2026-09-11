#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <iostream>

void p_in_minus_p_out()
{
    // Set nice drawing defaults
    gStyle->SetOptStat(0);
    gStyle->SetTextFont(42);

    // 1. Open the input file
    TFile *inFile1 = TFile::Open("analysis_output_califa_gt15MeV.root", "READ");
    TFile *inFile2 = TFile::Open("analysis_output_califa_lt1MeV.root", "READ");
    if (inFile1->IsZombie())
    {
        std::cerr << "[ERROR] Could not open input file analysis_output_califa_gt15MeV.root\n";
        return;
    }
    if (inFile2->IsZombie())
    {
        std::cerr << "[ERROR] Could not open input file analysis_output_califa_lt1MeV.root\n";
        return;
    }

    // 2. Retrieve histograms created in analyse_all
    TH1F *h_0hit   = (TH1F*) inFile2->Get("h1_p_in_minus_p_out_0");
    TH1F *h_1hit  = (TH1F*) inFile1->Get("h1_p_in_minus_p_out_1hit_0");
    TH1F *h_2hit  = (TH1F*) inFile1->Get("h1_p_in_minus_p_out_2hit_0");
    TH1F *h_3hit  = (TH1F*) inFile1->Get("h1_p_in_minus_p_out_3hit_0"); 

    if (!h_0hit || !h_1hit || !h_2hit || !h_3hit)
    {
        std::cerr << "[ERROR] Could not retrieve one or more histograms from file.\n";
        inFile1->Close();
        inFile2->Close();
        return;
    }

    // 3. Customize Visual Appearance (Colors & Styles)
    h_0hit->SetLineColor(kBlack);
    h_0hit->SetLineWidth(3);
    h_0hit->Rebin(6); // Optional: Rebin for better visibility

    h_1hit->SetLineColor(kBlue + 1);
    h_1hit->SetLineWidth(2);
    h_1hit->Rebin(6); // Optional: Rebin for better visibility

    h_2hit->SetLineColor(kRed + 1);
    h_2hit->SetLineWidth(2);
    h_2hit->Rebin(6); // Optional: Rebin for better visibility

    h_3hit->SetLineColor(kGreen + 2);
    h_3hit->SetLineWidth(2);
    h_3hit->Rebin(6); // Optional: Rebin for better visibility

    // 4. Create Canvas & Draw
    TCanvas *c1 = new TCanvas("c_ptotal", "Incoming minus outgoing momentum comparison", 800, 600);
    c1->SetGrid();

    // Determine max Y scale across all histograms so nothing gets clipped
    double max_y = h_0hit->GetMaximum();
    h_0hit->SetMaximum(max_y * 1.15); // Add 15% head room for legend

    h_1hit->Scale(h_0hit->Integral()/h_1hit->Integral());
    h_2hit->Scale(h_0hit->Integral()/h_2hit->Integral());
    h_3hit->Scale(h_0hit->Integral()/h_3hit->Integral());
    h_0hit->Draw("HIST");
    h_1hit->Draw("HIST SAME");
    h_2hit->Draw("HIST SAME");
    h_3hit->Draw("HIST SAME");

    // 5. Add Legend
    TLegend *legend = new TLegend(0.62, 0.68, 0.88, 0.88);
    legend->SetBorderSize(1);
    legend->SetFillColor(kWhite);
    legend->SetTextSize(0.035);
    legend->AddEntry(h_0hit,  "No CALIFA Hits",           "l");
    legend->AddEntry(h_1hit, "CALIFA 1 Hit",         "l");
    legend->AddEntry(h_2hit, "CALIFA 2 Hits",        "l");
    legend->AddEntry(h_3hit, "CALIFA 3 Hits",        "l");
    legend->Draw();

    // Update canvas
    c1->Update();

    // 6. Save Canvas and Histograms to output ROOT file
    TFile *outFile = new TFile("p_in_minus_p_out.root", "RECREATE");
    c1->Write();
    h_0hit->Write("h1_p_in_minus_p_out_0");
    h_1hit->Write("h1_p_in_minus_p_out_1hit_0");
    h_2hit->Write("h1_p_in_minus_p_out_2hit_0");
    h_3hit->Write("h1_p_in_minus_p_out_3hit_0");

    outFile->Close();
    inFile1->Close();
    inFile2->Close();

    std::cout << "\n[INFO] Macro finished successfully. Canvas & histograms saved to p_in_minus_p_out.root\n";
}