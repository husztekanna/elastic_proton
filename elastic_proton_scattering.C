#include "definitions.h"
const double CALIFA_THRESHOLD = 5000; // keV: for no califa:1000, for califa: 15000
using namespace std;
const double m_proton = 938.272; //MeV
namespace
{

    bool configureFragmentParameters(const TString &outgoing_fragment, double &frag_Z, double &frag_mass, double &beam_mass, double &incoming_E_loss, double &outgoing_E_loss)
    {

    if (outgoing_fragment == "24O")
    {
        frag_Z = 8;
        frag_mass = 22.370272* 1000; // GeV
        beam_mass = frag_mass;
        incoming_E_loss = 0;
        outgoing_E_loss = 0;
    }
    else if (outgoing_fragment == "24F")
    {
        frag_Z = 9;
        frag_mass = 22.358806* 1000; // GeV
        beam_mass = frag_mass;
        incoming_E_loss = 10.02;
        outgoing_E_loss = 4.8;
    }
    else if (outgoing_fragment == "23F")
    {
        frag_Z = 9;
        frag_mass = 21.423052* 1000; // GeV
        beam_mass = frag_mass;
        incoming_E_loss = 0;
        outgoing_E_loss = 0;
    }
    else if (outgoing_fragment == "23O")
    {
        frag_Z = 8;
        frag_mass = 21.434899* 1000; // GeV
        beam_mass = frag_mass;
        incoming_E_loss = 0;
        outgoing_E_loss = 0;
    }
    else if (outgoing_fragment == "22O")
    {
        frag_Z = 8;
        frag_mass = 20.498066* 1000; // GeV
        beam_mass = frag_mass;
        incoming_E_loss = 8.89;
        outgoing_E_loss = 4.31;
    }
    else if (outgoing_fragment == "25F")
    {
        frag_Z = 9;
        frag_mass = 23.294089* 1000; // GeV
        beam_mass = frag_mass;
        incoming_E_loss = 10.02;
        outgoing_E_loss = 4.8;
    }
    else
    {
        std::cerr << "[ERROR] Unknown outgoing fragment: " << outgoing_fragment << "\n";
        return false;
    }
    
        return true;
    }

    void clearEventCollections(TClonesArray *tofd_array, TClonesArray *foot_array, TClonesArray *frs_array,
                               TClonesArray *fragment_tracks, TClonesArray *incoming_tracks, TClonesArray *califa_data)
    {
        if (tofd_array)
            tofd_array->Clear();
        if (foot_array)
            foot_array->Clear();
        if (frs_array)
            frs_array->Clear();
        if (fragment_tracks)
            fragment_tracks->Clear();
        if (incoming_tracks)
            incoming_tracks->Clear();
        if (califa_data)
            califa_data->Clear();
    }

    void setBranchAddresses(TChain *tree, R3BEventHeader **header, TClonesArray **frs_array, TClonesArray **tofd_array,
                            TClonesArray **foot_array, TClonesArray **fragment_tracks, TClonesArray **incoming_tracks,
                            TClonesArray **alpide_tracks, TClonesArray **califa_data)
    {
        tree->SetBranchAddress("EventHeader.", header);
        tree->SetBranchAddress("FrsData", frs_array);
        tree->SetBranchAddress("TofdHit", tofd_array);
        tree->SetBranchAddress("FootHitData", foot_array);
        tree->SetBranchAddress("FragmentMDFTrack", fragment_tracks);
        tree->SetBranchAddress("IncomingTrackFoot", incoming_tracks);
        tree->SetBranchAddress("AlpideTrack", alpide_tracks);
        tree->SetBranchAddress("CalifaClusterData", califa_data);
    }

    R3BFrsData *findFrsData(TClonesArray *frs_array)
    {
        if (!frs_array)
            return nullptr;

        for (int j = 0; j < frs_array->GetEntriesFast(); ++j)
        {
            auto item = dynamic_cast<R3BFrsData *>(frs_array->At(j));
            if (!item)
                continue;
            if (item->GetStaId() == 2)
                return item;
        }

        for (int j = 0; j < frs_array->GetEntriesFast(); ++j)
        {
            auto item = dynamic_cast<R3BFrsData *>(frs_array->At(j));
            if (!item)
                continue;
            if (item->GetStaId() == 1)
                return item;
        }

        return nullptr;
    }

    std::vector<R3BCalifaClusterData *> collectCalifaHits(TClonesArray *califa_data, double threshold)
    {
        std::vector<R3BCalifaClusterData *> hits;
        if (!califa_data)
            return hits;

        for (int j = 0; j < califa_data->GetEntriesFast(); ++j)
        {
            auto califa_item = dynamic_cast<R3BCalifaClusterData *>(califa_data->At(j));
            if (califa_item && califa_item->GetEnergy() > threshold)
            {
                hits.push_back(califa_item);
            }
        }
        return hits;
    }

} // namespace

void analyse_all(TString outgoing_fragment = "25F")
{
    TRandom3 rng(0);
    double frag_Z = 0;
    double frag_mass = 0;
    double beam_mass = 0;
    double incoming_E_loss = 0;
    double outgoing_E_loss = 0;
   
    if (!configureFragmentParameters(outgoing_fragment, frag_Z, frag_mass, beam_mass, incoming_E_loss, outgoing_E_loss))
    {
        return;
    }

    double beam_mass_amu = beam_mass / 931.494;
    double frag_mass_amu = frag_mass / 931.494;

    TChain *tree = new TChain("evt");
    // tree->Add("/lustre/r3b/vpanin/G249/rootfiles_neuland_clustering/filtered/filtered_24F_to_all.root");
    //tree->Add("/lustre/r3b/vpanin/G249/rootfiles_neuland_clustering/filtered/filtered_25F_to_all.root");
    //for 25F and 24F: TString input_file = Form("/lustre/r3b/vpanin/G249/rootfiles_25F_all_data/filtered/filtered_%s_to_all.root", outgoing_fragment.Data());
    //for 22O :
    TString input_file = Form("/lustre/r3b/vpanin/G249/rootfiles_neuland_clustering/filtered/filtered_%s_to_all.root", outgoing_fragment.Data());
    tree->Add(input_file);

    std::cout << "\n\n-- Analyse TTree with " << tree->GetEntries() << " entries \n";

    // FRS
    auto h2_frs_pid = CreateHistograms<TH2F>("frs_pid", 1, 1000, 2.4, 3.2, 1000, 3, 12);
    // Fragment
    auto h2_AoZ_vs_Z = CreateHistograms<TH2F>("h2_AoZ_vs_Z", 1, 1000, 1.5, 3.8, 1000, 1, 12);
    // Califa
    auto h1_califa_mult = CreateHistograms<TH1F>("h1_califa_mult", 1, 20, 0, 20);
    // Vertex
    auto h2_vertex_XY = CreateHistograms<TH2F>("h2_vertex_XY", 1, 400, -100, 100, 400, -100, 100);
    auto h2_vertex_XZ = CreateHistograms<TH2F>("h2_vertex_XZ", 1, 1000, -500, 500, 1000, -500, 500);
    auto h2_vertex_YZ = CreateHistograms<TH2F>("h2_vertex_YZ", 1, 1000, -500, 500, 1000, -500, 500);
    auto h1_DCA = CreateHistograms<TH1F>("h1_DCA", 1, 1000, 0, 10);
    auto h1_vertex_Z = CreateHistograms<TH1F>("h1_vertex_Z", 1, 500, -2000, 2000);
    auto h2_TX_TY_diff_foot_alpide = CreateHistograms<TH2F>("h2_TX_TY_diff_foot_alpide", 1, 1000, -0.08, 0.08, 1000, -0.08, 0.08);
     auto h2_charge_foot_vs_pos = CreateHistograms<TH2F>("h2_charge_foot_vs_pos", 4, 2000, -50, 50, 1000, 0, 15);
    // Beam momentum
    auto h1_beam_momentum = CreateHistograms<TH1F>("h1_beam_momentum", 1, 500, 20000, 40000);
    // Fragment momentum
    auto h1_frag_px = CreateHistograms<TH1F>("h1_frag_px", 1, 200, -1000, 1000);
    auto h1_frag_py = CreateHistograms<TH1F>("h1_frag_py", 1, 200, -1000, 1000);
    auto h1_frag_ptr = CreateHistograms<TH1F>("h1_frag_ptr", 1, 200, 0, 1000);
    auto h1_frag_pz_lab = CreateHistograms<TH1F>("h1_frag_pz_lab", 1, 500, 20000, 40000);
    auto h1_frag_ptotal_lab = CreateHistograms<TH1F>("h1_frag_ptotal_lab", 1, 500, 20000, 40000);
    auto h1_p_in_minus_p_out = CreateHistograms<TH1F>("h1_p_in_minus_p_out", 1, 1000, -1000, 2000);
    auto h2_frag_px_vs_py = CreateHistograms<TH2F>("h2_frag_px_vs_py", 1, 200, -1000, 1000, 200, -1000, 1000);
    auto h1_califa_energy_1hit = CreateHistograms<TH1F>("h1_califa_energy_1hit", 1, 2000, 0, 1000);
    auto h1_califa_theta_1hit = CreateHistograms<TH1F>("h1_califa_theta_1hit", 1, 180, 0, 90);
    auto h1_califa_phi_1hit = CreateHistograms<TH1F>("h1_califa_phi_1hit", 1, 360, -180, 180);
    auto h2_califa_energy_vs_theta_1hit = CreateHistograms<TH2F>("h2_califa_energy_vs_theta_1hit", 1, 300, 0, 90, 2000, 0, 1000);
    auto h2_califa_energy_vs_T = CreateHistograms<TH2F>("h2_califa_energy_vs_T", 1, 300, 0, 200, 300, 0, 200);
    auto h2_califa_energy_vs_T2 = CreateHistograms<TH2F>("h2_califa_energy_vs_T2", 1, 300, 0, 200, 300, 0, 200);
    auto h1_missing_mass = CreateHistograms<TH1F>("h1_missing_mass", 1, 200, 0, 2000);
    auto h2_T_vs_theta_proton = CreateHistograms<TH2F>("h2_T_vs_theta_proton", 1, 300, 0, 180, 2000, 0, 1000);
    auto h2_thetaCalifa_vs_theta_proton = CreateHistograms<TH2F>("h2_thetaCalifa_vs_theta_proton", 1, 300, 0, 180, 300, 0, 180);
    auto h2_phiCalifa_vs_phi_proton = CreateHistograms<TH2F>("h2_phiCalifa_vs_phi_proton", 1, 360, -180, 180, 360, -180, 180);
    auto h2_MissM_vs_MissP = CreateHistograms<TH2F>("h2_MissM_vs_MissP", 1, 300, 0, 2000, 300, 0, 2000);
    auto h1_phi_califa_minus_phi_proton = CreateHistograms<TH1F>("h1_phi_califa_minus_phi_proton", 1, 360, -180, 180);
    auto h1_mandelstam_t = CreateHistograms<TH1F>("h1_mandelstam_t", 1, 1000, -10, 0);
    auto h1_theta_cm = CreateHistograms<TH1F>("h1_theta_cm", 1, 180, 0, 90);


    R3BEventHeader *header = nullptr;
    TClonesArray *frs_array = nullptr;
    TClonesArray *tofd_array = nullptr;
    TClonesArray *foot_array = nullptr;
    TClonesArray *fragment_tracks = nullptr;
    TClonesArray *incoming_tracks = nullptr;
    TClonesArray *alpide_tracks = nullptr;
    TClonesArray *califa_data = nullptr;

    setBranchAddresses(tree, &header, &frs_array, &tofd_array, &foot_array, &fragment_tracks, &incoming_tracks,
                       &alpide_tracks, &califa_data);

    // Create tree to store simple variables
    TTree *outTree = new TTree("tree", "Output tree with simple variables");
    Double_t frag_px = 0, frag_py = 0, frag_pz_lab = 0, frag_ptr = 0, frag_beta = 0, frag_ptotal_lab = 0;
    Double_t frag_pz_rf = 0, beam_beta = 0, beam_gamma = 0, frag_ptot_lab = 0, frag_ptot_rf = 0;
    Double_t p1_theta = 0, p1_phi = 0, p2_theta = 0, p2_phi = 0;
    Double_t p_opang = 0, p_dif_phi = 0;
    Double_t n_beta = 0, n_gamma = 0, n_theta = 0, n_phi = 0, n_ptot_lab = 0, n_ptot_rf = 0, n_tof = 0, n_angle_farg = 0, n_flight_path = 0, n_first_energy = 0;
    Double_t Excitation_energy = 0, total_px = 0, total_py = 0, total_pz_lab = 0, total_pz_rf = 0;
    Double_t p_in_minus_p_out = 0;
    Double_t incoming_p = 0, outgoing_p = 0;
    Double_t frs_Z = 0, frs_AoZ = 0, frs_Beta = 0, frs_Brho = 0;
    Double_t frag_AoQ = 0;
    Double_t vertex_x = 0, vertex_y = 0, vertex_z = 0, dca = 0;
    Double_t califa_mult = 0, theta_califa = 0, phi_califa = 0, califa_energy = 0;
    Double_t missing_mass = 0, missing_momentum = 0, theta_proton = 0, phi_proton = 0;
    Double_t mandelstam_t = 0, T = 0, T2 = 0, theta_scattering = 0, phi_scattering = 0,theta_cm = 0;
    Double_t Nf = 0, Ns = 0, charge_tofd = 0, tofd_bar_id = 0;
    Double_t NbOfCrystalHits = 0, TX_diff = 0, TY_diff = 0;
    Double_t err_mandelstam_t = 0;

    outTree->Branch("p_in_minus_p_out", &p_in_minus_p_out, "p_in_minus_p_out/D");
    outTree->Branch("frag_px", &frag_px, "frag_px/D");
    outTree->Branch("frag_py", &frag_py, "frag_py/D");
    outTree->Branch("frag_ptotal_lab", &frag_ptotal_lab, "frag_ptotal_lab/D");
    outTree->Branch("frag_pz_lab", &frag_pz_lab, "frag_pz_lab/D");
    outTree->Branch("frag_ptr", &frag_ptr, "frag_ptr/D");
    outTree->Branch("p1_theta", &p1_theta, "p1_theta/D");
    outTree->Branch("p1_phi", &p1_phi, "p1_phi/D");
    outTree->Branch("p2_theta", &p2_theta, "p2_theta/D");
    outTree->Branch("p2_phi", &p2_phi, "p2_phi/D");
    outTree->Branch("p_opang", &p_opang, "p_opang/D");
    outTree->Branch("p_dif_phi", &p_dif_phi, "p_dif_phi/D");
    outTree->Branch("beam_beta", &beam_beta, "beam_beta/D");
    outTree->Branch("beam_gamma", &beam_gamma, "beam_gamma/D");
    outTree->Branch("frag_beta", &frag_beta, "frag_beta/D");
    outTree->Branch("frag_pz_rf", &frag_pz_rf, "frag_pz_rf/D");
    outTree->Branch("frag_ptot_rf", &frag_ptot_rf, "frag_ptot_rf/D");
    outTree->Branch("frag_ptot_lab", &frag_ptot_lab, "frag_ptot_lab/D");
    outTree->Branch("incoming_p", &incoming_p, "incoming_p/D");
    outTree->Branch("outgoing_p", &outgoing_p, "outgoing_p/D");
    outTree->Branch("frs_Z", &frs_Z, "frs_Z/D");
    outTree->Branch("frs_AoZ", &frs_AoZ, "frs_AoZ/D");
    outTree->Branch("frag_AoQ", &frag_AoQ, "frag_AoQ/D");
    outTree->Branch("vertex_x", &vertex_x, "vertex_x/D");
    outTree->Branch("vertex_y", &vertex_y, "vertex_y/D");
    outTree->Branch("vertex_z", &vertex_z, "vertex_z/D");
    outTree->Branch("dca", &dca, "dca/D");
    outTree->Branch("califa_mult", &califa_mult, "califa_mult/D");
    outTree->Branch("theta_califa", &theta_califa, "theta_califa/D");
    outTree->Branch("phi_califa", &phi_califa, "phi_califa/D");
    outTree->Branch("califa_energy", &califa_energy, "califa_energy/D");
    outTree->Branch("missing_mass", &missing_mass, "missing_mass/D");
    outTree->Branch("missing_momentum", &missing_momentum, "missing_momentum/D");
    outTree->Branch("theta_proton", &theta_proton, "theta_proton/D");
    outTree->Branch("phi_proton", &phi_proton, "phi_proton/D");
    outTree->Branch("mandelstam_t", &mandelstam_t, "mandelstam_t/D");
    outTree->Branch("T", &T, "T/D");
    outTree->Branch("T2", &T2, "T2/D");
    outTree->Branch("theta_cm", &theta_cm, "theta_cm/D");
    outTree->Branch("Nf", &Nf, "Nf/D");
    outTree->Branch("Ns", &Ns, "Ns/D");
    outTree->Branch("NbOfCrystalHits", &NbOfCrystalHits, "NbOfCrystalHits/D");
    outTree->Branch("charge_tofd", &charge_tofd, "charge_tofd/D");
    outTree->Branch("tofd_bar_id", &tofd_bar_id, "tofd_bar_id/D");
    outTree->Branch("phi_scattering", &phi_scattering, "phi_scattering/D");
    outTree->Branch("theta_scattering", &theta_scattering, "theta_scattering/D");
    outTree->Branch("TX_diff", &TX_diff, "TX_diff/D");
    outTree->Branch("TY_diff", &TY_diff, "TY_diff/D");
    outTree->Branch("err_mandelstam_t", &err_mandelstam_t, "err_mandelstam_t/D");

    h1_mandelstam_t[0]->Sumw2();

    Long64_t nEntries = tree->GetEntries(); // to analyse all data in the tree
    // Long64_t nEntries = 100000;
    for (Long64_t i = 0; i < nEntries; ++i)
    {
        clearEventCollections(tofd_array, foot_array, frs_array, fragment_tracks, incoming_tracks, califa_data);

        tree->GetEntry(i);
        cout << "\r Number of events processed: " << i << flush;

        frag_ptr = 0;

        if (header->GetTrigger() != 1)
            continue;
        if (header->GetTpat() > 32)
            continue;
        if (incoming_tracks->GetEntriesFast() == 0)
            continue;
        if (fragment_tracks->GetEntriesFast() == 0)
            continue;

        //-------- Get CALIFA data
        auto califa_hits = collectCalifaHits(califa_data, CALIFA_THRESHOLD);
        h1_califa_mult[0]->Fill(califa_hits.size());
        if (califa_hits.size() != 1)
        {
            califa_hits.clear();
            continue;
        }
        theta_califa = califa_hits[0]->GetTheta() * TMath::RadToDeg();
        phi_califa = califa_hits[0]->GetPhi() * TMath::RadToDeg();
        califa_energy = califa_hits[0]->GetEnergy() / 1000; // MeV
        Nf = califa_hits[0]->GetNf();
        Ns = califa_hits[0]->GetNs();
        califa_mult = califa_hits.size();
        NbOfCrystalHits  = califa_hits[0]->GetNbOfCrystalHits();
       // if (theta_califa < 50 || theta_califa > 100)
       // {
       //     califa_hits.clear();
       //     continue;
       // }
        califa_hits.clear();

        //-------- Get FRS data
        frs_Z = 0;
        frs_AoZ = 0;
        frs_Beta = 0;
        frs_Brho = 0;
        R3BFrsData *frs_item = findFrsData(frs_array);

        if (frs_item)
        {
            frs_Z = frs_item->GetZ();
            frs_AoZ = frs_item->GetAq();
            frs_Beta = frs_item->GetBeta();
            frs_Brho = frs_item->GetBrho();
        }
        if (frs_Z == 0 || frs_AoZ == 0)
            continue;

        // Convert MeV mass back to amu (~931.494 MeV/u) if required by momentum_after_energy_loss
        double frag_mass_amu = frag_mass / 931.494;

        double initial_incoming_p = frs_Brho * frag_Z / 3.33564 * 1000; // MeV/c
        incoming_p = momentum_after_energy_loss(initial_incoming_p, beam_mass_amu, incoming_E_loss, "in");
        //cout << "\nInitial incoming momentum: " << initial_incoming_p << " MeV/c, after energy loss: " << incoming_p << " MeV/c\n";

        h1_beam_momentum[0]->Fill(incoming_p);


        //-------- Get incoming outgoing tracks
        auto fragment = dynamic_cast<R3BTrackingParticle *>(fragment_tracks->At(0));
        auto intrack = dynamic_cast<R3BTrackingParticle *>(incoming_tracks->At(0));
        auto alpide_track = dynamic_cast<R3BTrackingParticle *>(alpide_tracks->At(0));
        if (!fragment || !intrack || !alpide_track)
            continue;

        auto out_track_mom = fragment->GetStartMomentum();
        auto out_track_pos = fragment->GetStartPosition();
        auto out_track_mom_unit = out_track_mom.Unit();
        frag_Z = fragment->GetCharge();
        frag_AoQ = fragment->GetMass();
        auto PoQ = out_track_mom.Mag();

        if (outgoing_fragment == "25F")
        {
            if (frag_AoQ < 2.68 || frag_AoQ > 2.82 || frag_Z > 9.6 || frag_Z < 8.55)
                continue;
        }
        if (outgoing_fragment == "24O")
        {
            if (frag_AoQ < 2.885 || frag_AoQ > 3.02 || frag_Z > 8.8 || frag_Z < 7.5)
                continue;
        }
        if (outgoing_fragment == "24F")
        {
            if (frag_AoQ < 2.57 || frag_AoQ > 2.68 || frag_Z > 9.59 || frag_Z < 8.48)
                continue;
        }
        if (outgoing_fragment == "23F")
        {
            if (frag_AoQ < 2.45 || frag_AoQ > 2.57 || frag_Z > 9.59 || frag_Z < 8.48)
                continue;
        }
        if (outgoing_fragment == "23O")
        {
            if (frag_AoQ < 2.8 || frag_AoQ > 2.88 || frag_Z > 8.8 || frag_Z < 7.5)
                continue;
        }
        if (outgoing_fragment == "22O")
        {
            if (frag_AoQ < 2.66 || frag_AoQ > 2.78 || frag_Z > 8.8 || frag_Z < 7.5)
                continue;
        }

        // Align outgoing momentum with incoming momentum in the center of the 50 mm target
        double initial_outgoing_p = PoQ * frag_Z * 1000; // MeV/c
        outgoing_p = momentum_after_energy_loss(initial_outgoing_p, frag_mass_amu, outgoing_E_loss, "out") ; //+ 48.28
        out_track_mom.SetMag(outgoing_p);
        // cout << "Initial outgoing momentum: " << initial_outgoing_p << " MeV/c, after energy loss: " << outgoing_p << " MeV/c\n";

        auto in_track_pos = intrack->GetStartPosition();
        auto in_track_mom = intrack->GetStartMomentum();
        auto in_track_mom_unit = in_track_mom.Unit();
        in_track_mom.SetMag(incoming_p);

        auto in_track_alpide_pos = alpide_track->GetStartPosition();
        auto in_track_alpide_mom = alpide_track->GetStartMomentum();
        in_track_alpide_mom.SetMag(incoming_p);

        TX_diff = in_track_alpide_mom.X() / in_track_alpide_mom.Z() - in_track_mom.X() / in_track_mom.Z();
        TY_diff = in_track_alpide_mom.Y() / in_track_alpide_mom.Z() - in_track_mom.Y() / in_track_mom.Z();
       // if (fabs(TX_diff) > 0.0015 || fabs(TY_diff) > 0.0015)
       //     continue;
        h2_TX_TY_diff_foot_alpide[0]->Fill(TX_diff, TY_diff);

        // Get fragment momentum in the frame of incident beam
        TVector3 v1p(0, 0, 1e-10);                   // helping vector for the final alignment
        TVector3 out_track_mom_aligned(0, 0, 1e-10); // aligned with incoming beam
        // TVector3 v_random(gRandom->Uniform(0., 0.1), gRandom->Uniform(0., 0.1), in_track_mom.Mag()); // to define axes
        TVector3 v_random(1, 0., 0.); // lab X axis to define other axes

        // Redefine coordinates relative to the incoming beam
        TVector3 Zaxis = in_track_mom;
        TVector3 Yaxis = (Zaxis.Unit()).Cross(v_random);
        TVector3 Xaxis = (Yaxis.Unit()).Cross(Zaxis.Unit());
        // And define unit vectors along each new axis in momentum units
        TVector3 Zu = Zaxis.Unit();
        TVector3 Xu = Xaxis.Unit();
        TVector3 Yu = Yaxis.Unit();
        out_track_mom_aligned.SetXYZ(out_track_mom.Dot(Xu), out_track_mom.Dot(Yu), out_track_mom.Dot(Zu));
        out_track_mom_aligned.SetMag(outgoing_p);

        //------ Target vertex
        auto target_vertex = get_vertex(in_track_pos, out_track_pos, in_track_mom_unit, out_track_mom_unit);
        vertex_z = target_vertex.vertexPos.Z();
        vertex_y = target_vertex.vertexPos.Y();
        vertex_x = target_vertex.vertexPos.X();
        dca = target_vertex.dca;
       // if (vertex_z < -80. || vertex_z > 180.)
       //     continue; // vertex z cut

        // fooot data
        double foot_charge = -1;
        for (int f = 0; f < 4; ++f)
        {
            auto foot_name = "foot" + to_string(f + 5);
            if (fragment->GetHitIndexByName(foot_name) < 0)
                continue;
            auto foot_item = dynamic_cast<R3BFootHitData *>(foot_array->At(fragment->GetHitIndexByName(foot_name)));
            if (!foot_item)
                continue;

            foot_charge = foot_item->GetZCharge();
            auto foot_position = foot_item->GetPos();
            auto foot_cluster_size = foot_item->GetMulStrip();

            h2_charge_foot_vs_pos[f]->Fill(foot_position, foot_charge);
            h2_charge_foot_vs_pos[f]->SetName(Form("foot %d", f + 5));
        }
        //-------- Get TOFD data
        if (fragment->GetHitIndexByName("tofd") < 0)
            continue;
        auto tofd_item = dynamic_cast<R3BTofdHitData *>(tofd_array->At(fragment->GetHitIndexByName("tofd")));
        if (!tofd_item)
            continue;
        charge_tofd = tofd_item->GetEloss();
        tofd_bar_id = tofd_item->GetBarId();

        //-------- Fragment analysis
        double ptotal_lab = out_track_mom_aligned.Mag();
        TLorentzVector incoming_p4(0, 0, in_track_mom.Mag(), sqrt((in_track_mom.Mag() * in_track_mom.Mag()) + (beam_mass * beam_mass)));
        TLorentzVector outgoing_p4(out_track_mom_aligned.X(), out_track_mom_aligned.Y(), out_track_mom_aligned.Z(), sqrt(out_track_mom_aligned.Mag2() + pow(frag_mass, 2)));
        TLorentzVector target_p4(0, 0, 0, 938.272);
        TVector3 missing_momentum_proton = incoming_p4.Vect() + target_p4.Vect() - outgoing_p4.Vect();
        auto missing_p4_proton = incoming_p4 + target_p4 - outgoing_p4;

        // beam_beta = in_track_mom.Mag() / sqrt(pow(in_track_mom.Mag(), 2) + pow(beam_mass, 2));
        // frag_beta = out_track_mom_aligned.Mag() / sqrt(pow(out_track_mom_aligned.Mag(), 2) + pow(frag_mass, 2));
        
        // Fill tree with momentum components
        frag_px = out_track_mom_aligned.X();
        frag_py = out_track_mom_aligned.Y();
        frag_pz_lab = out_track_mom_aligned.Z();
        frag_ptr = sqrt(frag_px * frag_px + frag_py * frag_py);
        //frag_ptotal_lab = sqrt(frag_px * frag_px + frag_py * frag_py + frag_pz_lab * frag_pz_lab);
        frag_ptotal_lab = ptotal_lab;
        p_in_minus_p_out = incoming_p - outgoing_p;
        beam_beta = incoming_p4.Beta();
        frag_beta = outgoing_p4.Beta();
        beam_gamma = 1 / sqrt(1 - beam_beta * beam_beta);
        frag_pz_rf = beam_gamma * (frag_pz_lab - beam_beta * sqrt(pow(frag_mass, 2) + pow(frag_ptot_lab, 2)));
        frag_ptot_rf = sqrt(frag_px * frag_px + frag_py * frag_py + frag_pz_rf * frag_pz_rf);

        //cout << "Beta beam" << beam_beta << ", Beta frag" << frag_beta << ", Gamma beam" << beam_gamma << ", Pz_rf frag" << frag_pz_rf << ", Ptot_rf frag" << frag_ptot_rf << "\n";


        
        missing_mass = missing_p4_proton.M();
        missing_momentum = missing_momentum_proton.Mag();
        theta_proton = missing_momentum_proton.Theta() * TMath::RadToDeg(); // degrees
        phi_proton = missing_momentum_proton.Phi() * TMath::RadToDeg(); // degrees
        double mandelstam_t1 = (incoming_p4 - outgoing_p4).M2();
        T = -mandelstam_t1 / 2 / m_proton ;// MeV
        TVector3 outgoing_p3 = outgoing_p4.Vect();
        TVector3 incoming_p3 = incoming_p4.Vect();
        theta_scattering = outgoing_p3.Angle(incoming_p3); // radians
        phi_scattering = out_track_mom_aligned.Phi(); //radians
        double incoming_energy = incoming_p4.E();
        T2 = ( 2 * pow(incoming_p, 2) * pow(sin(theta_scattering / 2), 2) ) / (m_proton + 2 * incoming_energy * pow(sin(theta_scattering / 2), 2));
        mandelstam_t = - (T2 * 2 * m_proton) / 1000000; //GeV 
        //cout << "mandelstam_t = " << mandelstam_t << "\n";

        double mandelstam_s = incoming_p4.M2() + pow(m_proton, 2) + 2 * incoming_p4.E() * m_proton;
        double sqrt_mandelstam_s = sqrt(mandelstam_s);
        double p_cm = 1/(2 * sqrt_mandelstam_s) * sqrt((mandelstam_s - (pow(m_proton - beam_mass,2))) * (mandelstam_s - (pow(m_proton + beam_mass,2))));
        double cos_theta_cm = 1 + (mandelstam_t / (2 * pow(p_cm / 1000, 2)));
        theta_cm = acos(cos_theta_cm) * TMath::RadToDeg(); //deg 
        double error_theta_scattering = 0.0004; //Rad
       // cout<< "incoming_p = " << incoming_p <<" , incoming energy= " << incoming_energy << "\n";
        err_mandelstam_t = ( (2 * pow(m_proton/1000, 2) * pow(incoming_p/1000, 2) * sin(theta_scattering))/pow(m_proton/1000 + (2 * incoming_energy/1000 * pow(sin(theta_scattering/2), 2)), 2) ) * error_theta_scattering;
        //cout << "Error in mandelstam_t = " << err_mandelstam_t << "\n";


       // if (T < 30)
        //    continue;

        
        // Align califa hits with the incoming beam direction
        TVector3 califa_vector;
        v1p.SetMagThetaPhi(califa_energy, theta_califa * TMath::DegToRad(), phi_califa * TMath::DegToRad());
        califa_vector.SetXYZ(v1p.Dot(Xu), v1p.Dot(Yu), v1p.Dot(Zu));

      //  if(fabs(califa_vector.Phi() * 180. / 3.14 - phi_proton) > 17 ) continue;
      //  if(outgoing_p3.Phi() * 180. / 3.14 > 30 && outgoing_p3.Phi() * 180. / 3.14 < 150) continue;

        //h2_T_vs_theta_proton[0]->Fill(theta_proton, T);
        h2_T_vs_theta_proton[0]->Fill(theta_califa, T);
        h2_thetaCalifa_vs_theta_proton[0]->Fill(theta_proton, califa_vector.Theta() * 180. / 3.14);
        h2_phiCalifa_vs_phi_proton[0]->Fill(phi_proton, califa_vector.Phi() * 180. / 3.14);

        h2_vertex_XY[0]->Fill(target_vertex.vertexPos.X(), target_vertex.vertexPos.Y());
        h2_vertex_XZ[0]->Fill(target_vertex.vertexPos.Z(), target_vertex.vertexPos.X());
        h2_vertex_YZ[0]->Fill(target_vertex.vertexPos.Z(), target_vertex.vertexPos.Y());
        h1_vertex_Z[0]->Fill(target_vertex.vertexPos.Z());
        h1_DCA[0]->Fill(target_vertex.dca);
        h2_MissM_vs_MissP[0]->Fill(missing_mass, missing_momentum);
        h1_califa_energy_1hit[0]->Fill(califa_energy);
        h1_califa_theta_1hit[0]->Fill(califa_vector.Theta() * 180. / 3.14);
        h1_califa_phi_1hit[0]->Fill(califa_vector.Phi() * 180. / 3.14);
        h2_califa_energy_vs_theta_1hit[0]->Fill(califa_vector.Theta() * 180. / 3.14, califa_energy);
        h2_califa_energy_vs_T[0]->Fill(califa_energy, T);
        h2_califa_energy_vs_T2[0]->Fill(califa_energy, T2);
        h1_missing_mass[0]->Fill(missing_mass);
        h1_phi_califa_minus_phi_proton[0]->Fill(califa_vector.Phi() * 180. / 3.14 - phi_proton);
        h1_mandelstam_t[0]->Fill(mandelstam_t);
        h1_theta_cm[0]->Fill(theta_cm);

        // Accumulate uncertainty in quadrature: sigma_total = sqrt(sigma_current^2 + error_new^2)
        int bin = h1_mandelstam_t[0]->FindBin(mandelstam_t);
        double current_err = h1_mandelstam_t[0]->GetBinError(bin);
        double combined_err = sqrt(current_err * current_err + err_mandelstam_t * err_mandelstam_t);
        h1_mandelstam_t[0]->SetBinError(bin, combined_err);

        h1_mandelstam_t[0]->Draw("E1");

        //------ fragmment histograms
        h1_frag_px[0]->Fill(out_track_mom_aligned.X());
        h1_frag_py[0]->Fill(out_track_mom_aligned.Y());
        h1_frag_ptr[0]->Fill(out_track_mom_aligned.Perp());
        h1_frag_pz_lab[0]->Fill(out_track_mom_aligned.Z());
        h1_frag_ptotal_lab[0]->Fill(out_track_mom_aligned.Mag());
        h2_frag_px_vs_py[0]->Fill(out_track_mom_aligned.X(), out_track_mom_aligned.Y());
       
        h1_p_in_minus_p_out[0]->Fill(p_in_minus_p_out);

        h2_frs_pid[0]->Fill(frs_AoZ, frs_Z);
        h2_AoZ_vs_Z[0]->Fill(frag_AoQ, frag_Z);
        outTree->Fill();
    } // eventloop

    
    //TString output_name = "analysis_output.root";
    TString output_name = Form("analysis_output_%s.root", outgoing_fragment.Data());
    TFile *outputFile = new TFile(output_name, "RECREATE");

    std::vector<TH1 *> all_hists;
    all_hists.clear();
    appendVector(all_hists, h2_frs_pid);
    appendVector(all_hists, h2_AoZ_vs_Z);
    appendVector(all_hists, h2_vertex_XY);
    appendVector(all_hists, h2_vertex_YZ);
    appendVector(all_hists, h2_vertex_XZ);
    appendVector(all_hists, h1_DCA);
    appendVector(all_hists, h1_vertex_Z);
    appendVector(all_hists, h2_TX_TY_diff_foot_alpide);
    appendVector(all_hists, h1_califa_mult);
    appendVector(all_hists, h1_frag_px);
    appendVector(all_hists, h1_frag_py);
    appendVector(all_hists, h1_frag_ptr);
    appendVector(all_hists, h1_frag_pz_lab);
    appendVector(all_hists, h1_frag_ptotal_lab);
    appendVector(all_hists, h2_frag_px_vs_py);
    appendVector(all_hists, h1_p_in_minus_p_out);
    appendVector(all_hists, h1_califa_energy_1hit);
    appendVector(all_hists, h1_califa_theta_1hit);
    appendVector(all_hists, h1_califa_phi_1hit);
    appendVector(all_hists, h2_califa_energy_vs_theta_1hit);
    appendVector(all_hists, h2_califa_energy_vs_T);
    appendVector(all_hists, h2_califa_energy_vs_T2);
    appendVector(all_hists, h1_missing_mass);
    appendVector(all_hists, h2_T_vs_theta_proton);
    appendVector(all_hists, h2_thetaCalifa_vs_theta_proton);
    appendVector(all_hists, h2_MissM_vs_MissP);
    appendVector(all_hists, h2_charge_foot_vs_pos);
    appendVector(all_hists, h2_phiCalifa_vs_phi_proton);
    appendVector(all_hists, h1_beam_momentum);
    appendVector(all_hists, h1_phi_califa_minus_phi_proton);
    appendVector(all_hists, h1_mandelstam_t);
    appendVector(all_hists, h1_theta_cm);

    for (auto hist : all_hists)
        hist->Write();
    outTree->Write();
    outputFile->Close();
    std::cout << "\n\nHistograms saved to " << output_name << "\n";
}

int main(int argc, char **argv)
{
    TString outgoing_fragment = "25F"; // default value
    bool is_Erel = false;
    if (argc > 1)
    {
        outgoing_fragment = argv[1];
        if (outgoing_fragment != "24O" && outgoing_fragment != "24F" && outgoing_fragment != "23O" && outgoing_fragment != "22O" && outgoing_fragment != "23F" && outgoing_fragment != "25F")
        {
            std::cout << "Usage: " << argv[0] << " [24O|24F|23F|23O|22O|25F]\n";
            std::cout << "  24O - analyze 24O fragment (default)\n";
            std::cout << "  23O - analyze 23O fragment\n";
            std::cout << "  22O - analyze 23O fragment\n";
            std::cout << "  24F - analyze 24F fragment\n";
            std::cout << "  23F - analyze 24F fragment\n";
            std::cout << "  25F - analyze 25F fragment\n";
            std::cout << "  --Erel - to include Erel analysis \n";
            return 1;
        }
    }
    std::cout << "Analyzing fragment: " << outgoing_fragment << "\n";
    analyse_all(outgoing_fragment);
    return 0;
}