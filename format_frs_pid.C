void format_frs_pid() {
    // 1. Turn off the stats box globally
    gStyle->SetOptStat(0);

    // 2. Open the existing root file
    TFile *file = TFile::Open("frs_pid_cocktailbeam.root");
    if (!file || file->IsZombie()) return;

    // 3. Get the 2D histogram
    TH2F *h2 = (TH2F*)file->Get("h2_frs_pid");

    // 4. Create a clean canvas
    TCanvas *c1 = new TCanvas("c1", "FRS PID", 800, 600);
    c1->SetMargin(0.14, 0.04, 0.14, 0.04); // Adjust margins [left, right, bottom, top]

    // 5. Adjust title and axis styling
    h2->SetTitle(""); // Removes the main top title
    
    // X-Axis (A/Z)
    h2->GetXaxis()->SetTitle("A/Q");
    h2->GetXaxis()->SetTitleSize(0.07);
    h2->GetXaxis()->SetLabelSize(0.06);
    h2->GetXaxis()->SetTitleOffset(1);
    h2->GetXaxis()->SetRangeUser(2.65, 2.8);
    h2->GetXaxis()->SetTitleFont(62); // Bold X-Axis Title


    // Y-Axis (Z)
    h2->GetYaxis()->SetTitle("Z");
    h2->GetYaxis()->SetTitleSize(0.07);
    h2->GetYaxis()->SetLabelSize(0.06);
    h2->GetYaxis()->SetTitleOffset(1);
    h2->GetYaxis()->SetRangeUser(7.0, 10.5);
    h2->GetYaxis()->SetTitleFont(62); // Bold Y-Axis Title

    // 6. Draw and Save
    h2->Draw("COL"); // Draw with color palette
    c1->SaveAs("frs_pid_formatted.png");
}