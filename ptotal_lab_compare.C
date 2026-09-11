void ptotal_lab_compare()
{
    TFile *f1 = TFile::Open("analysis_output_califa_gt15MeV.root");
    TFile *f2 = TFile::Open("analysis_output_califa_lt1MeV.root");

    if (!f1 || !f2) {
        std::cout << "One of the files is missing.\n";
        return;
    }

    TH1 *h_gt = (TH1*)f1->Get("h1_frag_ptotal_lab_0");
    TH1 *h_lt = (TH1*)f2->Get("h1_frag_ptotal_lab_0");

    if (!h_gt) h_gt = (TH1*)f1->Get("h1_frag_ptotal_lab");
    if (!h_lt) h_lt = (TH1*)f2->Get("h1_frag_ptotal_lab");

    if (!h_gt || !h_lt) {
        std::cout << "Histogram not found. Check the exact name.\n";
        return;
    }

    h_gt = (TH1*)h_gt->Clone("h_gt_15MeV");
    h_lt = (TH1*)h_lt->Clone("h_lt_1MeV");

    h_gt->SetLineColor(kRed);
    h_gt->SetLineWidth(2);
    h_gt->SetStats(0);

    h_lt->SetLineColor(kBlue);
    h_lt->SetLineWidth(2);
    h_lt->SetStats(0);

    Double_t ymax = 1.15 * TMath::Max(h_gt->GetMaximum(), h_lt->GetMaximum());
    h_gt->SetMaximum(ymax);

    TCanvas *c = new TCanvas("c", "CALIFA comparison", 900, 600);
    c->SetGrid();

    h_gt->Scale(h_lt->Integral()/h_gt->Integral());
    h_gt->Draw("HIST");
    h_lt->Draw("HIST SAME");

    TLegend *leg = new TLegend(0.62, 0.72, 0.88, 0.88);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(h_gt, "gt 15 MeV", "l");
    leg->AddEntry(h_lt, "lt 1 MeV", "l");
    leg->Draw();

    TLatex *txt = new TLatex();
    txt->SetTextSize(0.035);
    txt->SetTextColor(kBlack);
    txt->DrawLatexNDC(0.15, 0.90, "Red = gt 15 MeV");
    txt->DrawLatexNDC(0.15, 0.85, "Blue = lt 1 MeV");

    c->Modified();
    c->Update();

    TFile *out = new TFile("ptotal_lab_comparison.root", "RECREATE");
    c->Write("califa_ptotal_comparison");
    out->Close();

    //c->SaveAs("ptotal_lab_comparison.pdf");
    //c->SaveAs("ptotal_lab_comparison.png");
}