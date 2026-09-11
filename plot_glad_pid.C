#include <TChain.h>

#include <TH2F.h>

#include <TH1F.h>

#include <TCanvas.h>

#include <TStyle.h>

#include <TClonesArray.h>

#include <TMath.h>

#include <TFile.h>

#include <TString.h>

#include <iostream>

#include "definitions.h"

//#include "/u/lndgst01/ahusztek/R3BRoot/tracking/R3BTrackingParticle.h"

//#include "/u/lndgst01/ahusztek/R3BRoot/r3bdata/R3BCalifaClusterData.h"

#include "R3BTrackingParticle.h"

#include "R3BCalifaClusterData.h"



using namespace std;



const double CALIFA_THRESHOLD = 5000.0; // keV



void plot_glad_pid(TString outgoing_fragment = "25F") {

    gStyle->SetOptStat(0);



    TChain *tree = new TChain("evt");

    TString input_file = Form("/lustre/r3b/vpanin/G249/rootfiles_neuland_clustering/filtered/filtered_%s_to_all.root", outgoing_fragment.Data());

    tree->Add(input_file);



    cout << "Processing " << tree->GetEntries() << " entries for GLAD PID..." << endl;



    // We use generic TClonesArray objects to bypass R3B class header inclusions

    TClonesArray *fragment_tracks = nullptr;

    TClonesArray *califa_data = nullptr;



    tree->SetBranchAddress("FragmentMDFTrack", &fragment_tracks);

    tree->SetBranchAddress("CalifaClusterData", &califa_data);



    TH2F *h2_glad_pid = new TH2F("h2_glad_pid", "Fragment PID (CALIFA Hits >= 1);A/Q;Z", 1000, 1.5, 3.8, 1000, 1, 12);



    Long64_t nEntries = tree->GetEntries();

   // nEntries = 1e5;

    for (Long64_t i = 0; i < nEntries; ++i) {

        if (fragment_tracks) fragment_tracks->Clear();

        if (califa_data) califa_data->Clear();



        tree->GetEntry(i);

        if (i % 10000 == 0) cout << "\rProcessing event " << i << " / " << nEntries << flush;



        // --- Check CALIFA Multiplicity (>= 1 hit above threshold) ---

        int califa_mult = 0;

        if (califa_data) {

            for (int j = 0; j < califa_data->GetEntriesFast(); ++j) {

                auto obj = dynamic_cast<R3BCalifaClusterData *>(califa_data->At(j));

                if (!obj) continue;

               

                // Read Energy directly via ROOT's dynamic method access

                double energy = obj->GetTitle() ? obj->GetEnergy() : 0.0;

                if (energy > CALIFA_THRESHOLD) {

                    califa_mult++;

                }

            }

        }



        if (califa_mult < 1) continue; // Require at least 1 CALIFA hit



        // --- Extract Fragment Tracks ---

        if (!fragment_tracks || fragment_tracks->GetEntriesFast() == 0) continue;

       

        auto frag_obj = dynamic_cast<R3BTrackingParticle *>(fragment_tracks->At(0));

        if (!frag_obj) continue;



        // Get Mass (A/Q) and Charge (Z) directly from the track object

        double frag_Z = frag_obj->GetCharge();

        double frag_AoQ = frag_obj->GetMass();
        frag_AoQ = frag_AoQ + 0.04;


        // Fill histogram without isotope cuts

        h2_glad_pid->Fill(frag_AoQ, frag_Z);

    }



    cout << "\nDraw and saving results..." << endl;



    TCanvas *c1 = new TCanvas("c1_pid", "GLAD PID", 800, 600);

    c1->SetMargin(0.14, 0.04, 0.14, 0.04);

// 5. Adjust title and axis styling

 h2_glad_pid->SetTitle(""); // Removes the main top title

// X-Axis (A/Z)

h2_glad_pid->GetXaxis()->SetTitle("A/Q");

h2_glad_pid->GetXaxis()->SetTitleSize(0.07);

h2_glad_pid->GetXaxis()->SetLabelSize(0.06);

h2_glad_pid->GetXaxis()->SetTitleOffset(1);

h2_glad_pid->GetXaxis()->SetRangeUser(2.4, 3.1);

h2_glad_pid->GetXaxis()->SetTitleFont(62); // Bold X-Axis Title

// Y-Axis (Z)

h2_glad_pid->GetYaxis()->SetTitle("Z");

h2_glad_pid->GetYaxis()->SetTitleSize(0.07);

h2_glad_pid->GetYaxis()->SetLabelSize(0.06);

h2_glad_pid->GetYaxis()->SetTitleOffset(1);

h2_glad_pid->GetYaxis()->SetRangeUser(7.5, 10);

h2_glad_pid->GetYaxis()->SetTitleFont(62); // Bold Y-Axis Title





    h2_glad_pid->Draw("COL");



    TString out_png = Form("glad_pid_califa_%s.png", outgoing_fragment.Data());

    c1->SaveAs(out_png);



    TFile *fOut = new TFile(Form("glad_pid_%s.root", outgoing_fragment.Data()), "RECREATE");

    h2_glad_pid->Write();

    fOut->Close();



    cout << "Finished! Saved to " << out_png << endl;

}



int main(int argc, char **argv)

{

    plot_glad_pid("25F"); // Default fragment

    return 0;

} 
