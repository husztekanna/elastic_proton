#include <TCanvas.h>
#include <TPad.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TStyle.h>

void overlay_cuts() {
    gStyle->SetOptStat(0); // Turn off statistics box

    // 1. Open both ROOT files
    TFile *file_uncut = TFile::Open("plot_without_cuts.root", "READ");
    TFile *file_cut   = TFile::Open("plot_with_cuts.root", "READ");

    if (!file_uncut || file_uncut->IsZombie()) {
        printf("Error: Could not open plot_without_cuts.root\n");
        return;
    }
    if (!file_cut || file_cut->IsZombie()) {
        printf("Error: Could not open plot_with_cuts.root\n");
        return;
    }

    // 2. Retrieve histograms for 25F
    TH1F *h_uncut = (TH1F*)file_uncut->Get("phi_califa_minus_phi_proton_25F");
    TH1F *h_cut   = (TH1F*)file_cut->Get("phi_califa_minus_phi_proton_25F");

    // NEW: the 2D angle-correlation histogram for the inset (uncut version)
    TH2F *h2_corr = (TH2F*)file_uncut->Get("phi_califa_vs_phi_proton_25F");

    if (!h_uncut || !h_cut) {
        printf("Error: 1D histograms missing from ROOT files!\n");
        return;
    }
    if (!h2_corr) {
        printf("Error: 2D correlation histogram missing from plot_without_cuts.root!\n");
        return;
    }

    // Prevent deletion when closing files
    h_uncut->SetDirectory(0);
    h_cut->SetDirectory(0);
    h2_corr->SetDirectory(0);

    file_uncut->Close();
    file_cut->Close();

    // 3. Apply Styling
    // Uncut histogram: Blue outline
    h_uncut->SetLineColor(kBlue + 1);
    h_uncut->SetLineWidth(2);
    h_uncut->GetYaxis()->SetTitleSize(0.07);
    h_uncut->GetXaxis()->SetTitleSize(0.07);
    h_uncut->GetYaxis()->SetTitleFont(62);
    h_uncut->GetXaxis()->SetTitleFont(62);
    h_uncut->GetXaxis()->SetTitle("#Phi_{miss} - #Phi_{proton} [deg]");
    h_uncut->GetYaxis()->SetTitle("Counts");
    h_uncut->GetXaxis()->SetLabelSize(0.06);
    h_uncut->GetYaxis()->SetLabelSize(0.06);
    h_uncut->GetXaxis()->SetTitleOffset(1);
    h_uncut->GetYaxis()->SetTitleOffset(1);

    // Cut histogram: Red line
    h_cut->SetLineColor(kRed + 1);
    h_cut->SetLineWidth(2);
    //h_cut->SetFillColorAlpha(kRed + 1, 0.35); // Semi-transparent red fill

    // Adjust Y-axis scale to fit the uncut peak headroom
    h_uncut->SetMaximum(h_uncut->GetMaximum() * 1.15);

    // NEW: styling for the inset 2D histogram
    h2_corr->SetTitle("");
    h2_corr->GetXaxis()->SetTitle("#Phi of the missing momentum [deg]");
    h2_corr->GetYaxis()->SetTitle("#Phi of proton in CALIFA [deg]");
    h2_corr->GetXaxis()->SetTitleSize(0.06);
    h2_corr->GetYaxis()->SetTitleSize(0.06);
    h2_corr->GetXaxis()->SetLabelSize(0.06);
    h2_corr->GetYaxis()->SetLabelSize(0.06);
    h2_corr->GetXaxis()->SetTitleOffset(1.0);
    h2_corr->GetYaxis()->SetTitleOffset(1.1);

    // 4. Create Canvas and Draw main histograms
    TCanvas *c = new TCanvas("c_phi_comparison", "Phi Coplanarity Cut Comparison", 800, 600);
    c->SetMargin(0.14, 0.04, 0.15, 0.04); // Adjust margins [left, right, bottom, top]
    c->SetGrid();

    h_uncut->Draw("HIST");
    h_cut->Draw("HIST SAME");

    // 5. Add Legend
    TLegend *leg = new TLegend(0.65, 0.75, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.045);
    leg->AddEntry(h_uncut, "Without Cuts", "l");
    leg->AddEntry(h_cut, "With Cuts", "l");
    leg->Draw();

    TPaveText *cutStats = new TPaveText(0.65, 0.58, 0.90, 0.72, "NDC");
    cutStats->SetBorderSize(0);
    cutStats->SetFillStyle(0);
    cutStats->SetTextAlign(12);
    cutStats->SetTextSize(0.035);
    cutStats->SetTextColor(kRed + 1);
    cutStats->AddText("Mean = -0.314 #pm 0.052 deg");
    cutStats->AddText("Sigma = 4.604 #pm 0.041 deg");
    cutStats->Draw();

    // 6. NEW: draw the 2D correlation as an inset pad, top-left corner
    c->cd(); // make sure we're back on the main canvas before adding the pad
    TPad *insetPad = new TPad("insetPad", "insetPad", 0.17, 0.58, 0.48, 0.90);
    insetPad->SetMargin(0.16, 0.16, 0.15, 0.05); // room for z-axis palette on the right
    insetPad->Draw();
    insetPad->SetFrameLineWidth(2);
    insetPad->cd();
    insetPad->SetLogz(); // Set log scale for Z-axis

    gStyle->SetPalette(kInvertedDarkBodyRadiator); // Optional: set a color palette for the 2D histogram

    h2_corr->SetMinimum(1); // Set min positive value so 0-count bins don't break log scale
    h2_corr->Draw("COLZ");

    c->cd(); // return to main pad before saving
    c->Update();

    // 7. Save Image
    c->SaveAs("phi_califa_minus_phi_proton_25F_with_inset.png");
}