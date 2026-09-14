
// Fills a GLAD-spectrometer PID histogram (A/Q vs Z) from raw fragment 
// tracking data, unfiltered by any A/Q or Z gate — used to inspect the full 
// PID spectrum before deciding on isotope-selection cuts (compare against 
// FRS PID from format_frs_pid.C).
//
// Applies only a basic trigger filter (header->GetTrigger() == 1) and 
// requires at least one fragment track with positive Z and A/Q.
//
// Depends on: definitions.h
//
// Input:  /lustre/r3b/vpanin/G249/rootfiles_neuland_clustering/main0158_0008.root
//         (TChain "evt")
// Output: glad_pid_cocktailbeam.root (TH2F "h2_glad_pid")
// Usage:  compiled as a standalone executable (has main()): ./plot_glad_pid
//         (can also be run interactively as a ROOT macro via 
//         plot_glad_pid_unfiltered())

#include "R3BEventHeader.h"
#include "R3BTrackingParticle.h"
#include "definitions.h"

#include <TChain.h>
#include <TClonesArray.h>
#include <TFile.h>
#include <TH2F.h>
#include <TString.h>
#include <TStyle.h>
#include <iostream>

using namespace std;

void glad_pid()
{
    gStyle->SetOptStat(0);

    // Hardcoded single input file path (SAME AS FRS)
    TString input_file = "/lustre/r3b/vpanin/G249/rootfiles_neuland_clustering/main0158_0008.root";
    TString output_file = "glad_pid_cocktailbeam.root";

    TChain *tree = new TChain("evt");
    tree->Add(input_file);

    cout << "\n-- Analyzing GLAD PID across " << tree->GetEntries() << " entries --\n";
    cout << "Input:  " << input_file << "\n";
    cout << "Output: " << output_file << "\n\n";

    // GLAD PID Histogram: A/Q vs Z (matching FRS range for easy side-by-side comparison)
    TH2F *h2_glad_pid = new TH2F("h2_glad_pid", "GLAD PID; A/Q; Z", 1000, 2.55, 2.95, 1000, 6, 12);

    h2_glad_pid->GetXaxis()->SetTitleOffset(1.2);
    h2_glad_pid->GetYaxis()->SetTitleOffset(1.2);
    h2_glad_pid->GetXaxis()->SetLabelSize(0.025);
    h2_glad_pid->GetYaxis()->SetLabelSize(0.025);

    // Set branch address for GLAD fragment tracking
    R3BEventHeader *header = nullptr;
    TClonesArray *fragment_tracks = nullptr;

    tree->SetBranchAddress("EventHeader.", &header);
    tree->SetBranchAddress("FragmentMDFTrack", &fragment_tracks);

    Long64_t nEntries = tree->GetEntries();

    for (Long64_t i = 0; i < nEntries; ++i)
    {
        if (fragment_tracks) fragment_tracks->Clear();
        tree->GetEntry(i);

        if (i % 100000 == 0)
        {
            cout << "\rProcessed events: " << i << " / " << nEntries << flush;
        }

        // Optional basic trigger filter to drop noise
        if (header && header->GetTrigger() != 1) continue;

        if (!fragment_tracks || fragment_tracks->GetEntriesFast() == 0) continue;

        auto frag_obj = dynamic_cast<R3BTrackingParticle *>(fragment_tracks->At(0));
        if (!frag_obj) continue;

        double frag_Z = frag_obj->GetCharge();
        double frag_AoQ = frag_obj->GetMass(); // A/Q in R3BTrackingParticle

        if (frag_Z > 0 && frag_AoQ > 0)
        {
            h2_glad_pid->Fill(frag_AoQ, frag_Z);
        }
    }

    cout << "\nDone processing.\n";

    // Write to output ROOT file
    TFile *outputFile = new TFile(output_file, "RECREATE");
    h2_glad_pid->Write();
    outputFile->Close();

    cout << "GLAD PID successfully saved to " << output_file << "\n";
}

int main(int argc, char **argv)
{
    plot_glad_pid_unfiltered();
    return 0;
}
