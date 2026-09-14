
// Shared constants, helper functions, and templates used across the 
// analysis macros. Included via #include "definitions.h" in files that 
// need FRS/CALIFA/FOOT data structures and R3B tracking classes.
//
// Contents:
//   - Physical constants: SPEED_OF_LIGHT, LOS_TO_FOOT5_DISTANCE, NEUTRON_MASS
//   - momentum_after_energy_loss(): corrects momentum for energy loss/gain 
//     through a target (incoming vs outgoing particle)
//   - CreateHistograms<HistT>(): templated helper to batch-create a vector 
//     of histograms with sequential names
//   - PlotOnCanvas<THist>(): templated helper to draw a vector of histograms 
//     onto a divided canvas (auto-selects "colz" for 2D, log-y for 1D)
//   - extractFromTCA<T>(): extracts a member value from every object in a 
//     TClonesArray via a member-function pointer
//   - appendVector(): concatenates a vector of histograms into a TH1* vector
//   - vertex struct + get_vertex(): computes the distance-of-closest-approach 
//     vertex between two tracked particles from their positions/slopes
//   - IsInsideEllipse(): checks whether a (z, A/Z) point falls inside an 
//     elliptical PID gate
//   - Init_A_and_Z(): identifies fragment A/Z from FRS PID ellipse cuts 
//     (currently only Z=9, A=25 active; other isotope gates commented out)
//
// Depends on: R3BRoot classes (R3BFootHitData, R3BHit, R3BFrsData, 
// R3BCalifaClusterData, R3BTofdHitData, R3BNeulandHit/Cluster, 
// R3BTrackingParticle, etc.), ROOT Minuit2/GSL minimizer libraries

#include <TFile.h>
#include <TTree.h>
#include <TVector3.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <iostream>
#include "R3BFootHitData.h"
#include "R3BLogger.h"
#include "R3BHit.h"
#include "R3BWRData.h"
#include "R3BFrsData.h"
#include "R3BCalifaClusterData.h"
#include "R3BTofdHitData.h"
#include "R3BFootHitData.h"
#include "R3BNeulandHit.h"
#include "R3BNeulandCluster.h"
#include "R3BEventHeader.h"
#include <TClonesArray.h>
#include <TApplication.h>
#include <TString.h>
#include <TStyle.h>
#include <TLorentzVector.h>
#include <TRotation.h>
#include "R3BTrackingParticle.h"
#include "R3BTPropagator.h"
#include "FairIon.h"
#include <array>
#include "R3BTrackingDetector.h"
#include <TStopwatch.h>
#include <TRandom3.h>
#include <TRandomGen.h>
#include "TMath.h"
#include "TChain.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/GSLMinimizer.h"
#include "Math/Minimizer.h"
#include "Minuit2/Minuit2Minimizer.h"

const  double SPEED_OF_LIGHT = 29.9792458;//cm/ns
const  double LOS_TO_FOOT5_DISTANCE = 254.97;//cm
const  double NEUTRON_MASS = 0.939565420;//MeV/c^2

double momentum_after_energy_loss(double initial_momentum_MeV_c, double mass_AMU, double dE_MeV_u, TString in_out)
{
    double mass_MeV = mass_AMU * 0.92907667; // Convert mass from AMU to MeV/c^2
    double initial_TKE = sqrt(initial_momentum_MeV_c * initial_momentum_MeV_c + mass_MeV * mass_MeV);
    double initial_E_MeV_u = initial_TKE / mass_AMU; // Initial kinetic energy in MeV/u

    double final_E_MeV_u = 0.0;

    if (in_out == "in")
    {
        final_E_MeV_u = initial_E_MeV_u - dE_MeV_u; // Energy loss for incoming particle
    }
    else if (in_out == "out")
    {
        final_E_MeV_u = initial_E_MeV_u + dE_MeV_u; // Energy gain for outgoing particle
    }
    else
    {
        std::cerr << "[ERROR] Unknown in_out parameter: " << in_out << "\n";
        return -1;
    }

    double final_TKE = final_E_MeV_u * mass_AMU;
    double final_momentum_MeV_c = sqrt(pow((final_TKE + mass_MeV), 2) - mass_MeV * mass_MeV);
    return final_momentum_MeV_c;
}

//Template helper fucntion to create a vector of n histograms of any type
template <typename HistT, typename... Args>
std::vector<HistT*> CreateHistograms(const char* base, int n, Args&&... args)
{
    std::vector<HistT*> hists;
    hists.reserve(n);
    for (int i = 0; i < n; ++i) {
        TString name  = Form("%s_%d", base, i);
        auto* h = new HistT(name, name, std::forward<Args>(args)...);
        h->SetDirectory(nullptr);
        //h->Sumw2();
        hists.push_back(h);
    }
    return hists;
}

template <typename THist>
TCanvas* PlotOnCanvas(const char* cname,
                        const std::vector<THist*>& hists,
                        int nx = 4, int ny = 3)
{
    auto* c = new TCanvas(cname, cname, 1200, 1000);
    c->Divide(nx, ny);

    int n = static_cast<int>(hists.size());
    for (int i = 0; i < n && i < nx * ny; ++i)
    {
        c->cd(i + 1);
        if (!hists[i]) continue;
        // automatic draw option depending on histogram dimension
        if (hists[i]->InheritsFrom(TH2::Class()))
            hists[i]->Draw("colz");
        else
            gPad->SetLogy();
            hists[i]->Draw();
    }
    c->Update();
    return c;
}

template <typename T, typename F>
auto extractFromTCA(TClonesArray* arr, F getter)
{
    using RetT = decltype((std::declval<T*>()->*getter)());
    std::vector<RetT> result;
    result.reserve(arr->GetEntriesFast());

    for (int i = 0; i < arr->GetEntriesFast(); ++i)
    {
        if (auto obj = dynamic_cast<T*>(arr->At(i)); obj)
            result.push_back((obj->*getter)());
        else
            result.push_back(RetT{}); // or skip if you prefer
    }
    return result;
}

//Quick function to join two vectors of histograms
template <typename V>
void appendVector(std::vector<TH1*>& dest, const std::vector<V>& src)
{
    dest.insert(dest.end(), src.begin(), src.end());
}

struct vertex
{
    Double_t dca;
    TVector3 vertexPos;
    TVector3 normalVec;
};

std::vector<vertex> vertex_array;

vertex get_vertex(TVector3 pos1, TVector3 pos2, TVector3 slope1, TVector3 slope2)
{
    Double_t t = 100000;
    Double_t v = 100000;
    TVector3 normalVec = slope2.Cross(slope1).Unit();

    TVector3 posDifference = pos2 - pos1;
    t = (-posDifference.Dot(slope2) + posDifference.Dot(slope1) * slope2.Dot(slope1)) / (1 - (slope2.Dot(slope1)) * (slope2.Dot(slope1)));
    v = (posDifference.Dot(slope1) - (posDifference.Dot(slope2) * (slope2.Dot(slope1)))) / (1 - (slope2.Dot(slope1)) * (slope2.Dot(slope1)));

    TVector3 vertex1Global = pos2 + t * slope2;
    TVector3 vertex2Global = pos1 + v * slope1;

    TVector3 Vertex_calc;
    TVector3 Vertex_dir_helper;

    Vertex_dir_helper = vertex2Global - vertex1Global;
    Vertex_dir_helper.SetX(Vertex_dir_helper.X() / 2);
    Vertex_dir_helper.SetY(Vertex_dir_helper.Y() / 2);
    Vertex_dir_helper.SetZ(Vertex_dir_helper.Z() / 2);

    Vertex_calc = vertex1Global + Vertex_dir_helper;
    vertex recVertex;
    recVertex.vertexPos = Vertex_calc;
    recVertex.dca = Vertex_dir_helper.Mag();
    recVertex.normalVec = normalVec;

    return recVertex;
}

bool IsInsideEllipse(double z, double aoz, double MeanZ, double RadiusZ, double MeanAoZ, double RadiusAoZ)
{
    Double_t ell = pow(((z - MeanZ) / (RadiusZ)), 2) + pow(((aoz - MeanAoZ) / (RadiusAoZ)), 2);
    return (ell <= 1.);
}

bool Init_A_and_Z(double &A, double &Z, double frag_AoZ, double frag_Z, double footZ) 
{
    double Zradius = 0.3;
    double AoZradius = 0.03;
    
//    //===== Z=10 isotopes
//    if(IsInsideFrsEllipse(frag_Z, frag_AoZ, 9.87, Zradius, 3.0, AoZradius)){
//        A = ; Z = 10;  return true;
//    }

    //===== Z=9 isotopes
    //if(IsInsideEllipse(frag_Z, frag_AoZ, 8.93, Zradius, 2.95778, AoZradius)){
    //    if(footZ > 9.3 || footZ < 8.7) return false;
    //    A = 27; Z = 9;  return true;
    //}
    
    //if(IsInsideEllipse(frag_Z, frag_AoZ, 8.93, Zradius, 2.84376, AoZradius)){
    //    if(footZ > 9.3 || footZ < 8.7) return false;
    //    A = 26; Z = 9;  return true;
    //}

    if(IsInsideEllipse(frag_Z, frag_AoZ, 8.93, Zradius, 2.7359, AoZradius)){
        if(footZ > 9.3 || footZ < 8.7) return false;
        A = 25; Z = 9;  return true;
    }

    //if(IsInsideEllipse(frag_Z, frag_AoZ, 8.93, Zradius, 2.63112, AoZradius)){
    //    if(footZ > 9.3 || footZ < 8.7) return false;
    //    A = 24; Z = 9;  return true;
    //}

    //if(IsInsideEllipse(frag_Z, frag_AoZ, 8.93, Zradius, 2.5171, AoZradius)){
    //    if(footZ > 9.3 || footZ < 8.7) return false;
    //    A = 23; Z = 9;  return true;
    //}

    //if(IsInsideEllipse(frag_Z, frag_AoZ, 8.93, Zradius, 2.39692, AoZradius)){
    //    if(footZ > 9.3 || footZ < 8.7) return false;
    //    A = 22; Z = 9;  return true;
    //}

    //if(IsInsideEllipse(frag_Z, frag_AoZ, 8.93, Zradius, 2.27982, AoZradius)){
    //    if(footZ > 9.3 || footZ < 8.7) return false;
    //    A = 21; Z = 9;  return true;
    //}

    ////===== Z=8 isotopes

    //if(IsInsideEllipse(frag_Z, frag_AoZ, 8.02, Zradius, 2.96086, AoZradius)){
    //    if(footZ > 8.3 || footZ < 7.7) return false;
    //    A = 24; Z = 8;  return true;
    //}

    //if(IsInsideEllipse(frag_Z, frag_AoZ, 8.02, Zradius, 2.83762, AoZradius)){
    //    if(footZ > 8.3 || footZ < 7.7) return false;
    //    A = 23; Z = 8;  return true;
    //}
    //
    //if(IsInsideEllipse(frag_Z, frag_AoZ, 8.02, Zradius, 2.71741, AoZradius)){
    //    if(footZ > 8.3 || footZ < 7.7) return false;
    //    A = 22; Z = 8;  return true;
    //}

    //if(IsInsideEllipse(frag_Z, frag_AoZ, 8.02, Zradius, 2.59106, AoZradius)){
    //    if(footZ > 8.3 || footZ < 7.7) return false;
    //    A = 21; Z = 8;  return true;
    //}

    //if(IsInsideEllipse(frag_Z, frag_AoZ, 8.02, Zradius, 2.45855, AoZradius)){
    //    if(footZ > 8.3 || footZ < 7.7) return false;
    //    A = 20; Z = 8;  return true;
    //}

    //if(IsInsideEllipse(frag_Z, frag_AoZ, 8.02, Zradius, 2.3322, AoZradius)){
    //    if(footZ > 8.3 || footZ < 7.7) return false;
    //    A = 19; Z = 8;  return true;
    //}


    //if(IsInsideFrsEllipse(frag_Z, frag_AoZ, 7.04, Zradius, 3.0, AoZradius)){
    //    A = ; Z = 7;  return true;
    //}
    //
    //if(IsInsideFrsEllipse(frag_Z, frag_AoZ, 6.07, Zradius, 3.0, AoZradius)){
    //    A = ; Z = 6;  return true;
    //}

    return false;
}

