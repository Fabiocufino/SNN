#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TMath.h"

using namespace std;

void draw_histograms(const char* filename = "clusters_ntuple.root") {
    TFile* file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        cerr << "Error: Cannot open file " << filename << endl;
        return;
    }

    // Access the "clusterValidIT" directory
    TDirectoryFile* dirIT = dynamic_cast<TDirectoryFile*>(file->Get("clusterValidIT"));
    TDirectoryFile* dirOT = dynamic_cast<TDirectoryFile*>(file->Get("clusterValidOT"));

    if (!dirIT) {
        cerr << "Error: Cannot access directory clusterValidIT" << endl;
        file->Close();
        return;
    }

    if (!dirOT) {
        cerr << "Error: Cannot access directory clusterValidOT" << endl;
        file->Close();
        return;
    }

    // Access the "tree" in the "clusterValidIT" directory
    TTree* IT = dynamic_cast<TTree*>(dirIT->Get("tree"));
    TTree* OT = dynamic_cast<TTree*>(dirOT->Get("tree;1"));

    if (!IT) {
        cerr << "Error: Cannot access tree in clusterValidIT" << endl;
        file->Close();
        return;
    }

    if (!OT) {
        cerr << "Error: Cannot access tree in clusterValidOT" << endl;
        file->Close();
        return;
    }

    // Create a canvas for the histograms
    TCanvas* canvas = new TCanvas("histograms", "Histograms", 1500, 1200);
    canvas->Divide(2, 2); // Divide canvas into 2 rows and 2 columns

    // Create TH1F histograms for r, phi, and eta
    TH1F* rHistType1 = new TH1F("rHistType1", "r Distribution", 100, 0, 1200);
    rHistType1->SetLineColor(kRed);
    TH1F* rHistType2 = new TH1F("rHistType2", "r Distribution", 100, 0, 1200);
    rHistType2->SetLineColor(kBlue);

    TH1F* phiHistType1 = new TH1F("phiHistType1", "phi Distribution", 100, -TMath::Pi(), TMath::Pi());
    phiHistType1->SetLineColor(kRed);
    TH1F* phiHistType2 = new TH1F("phiHistType2", "phi Distribution ", 100, -TMath::Pi(), TMath::Pi());
    phiHistType2->SetLineColor(kBlue);

    TH1F* etaHistType1 = new TH1F("etaHistType1", "eta Distribution ", 100, -1.1, 1.5);
    etaHistType1->SetLineColor(kRed);
    TH1F* etaHistType2 = new TH1F("etaHistType2", "eta Distribution", 100, -1.1, 1.5);
    etaHistType2->SetLineColor(kBlue);

    

    // Define variables to hold leaf values
    Float_t cluster_R;
    Float_t cluster_phi;
    Float_t cluster_type;
    Float_t eventID;
    Float_t cluster_eta;

    // Set branch addresses
    IT->SetBranchAddress("cluster_R", &cluster_R);
    IT->SetBranchAddress("cluster_phi", &cluster_phi);
    IT->SetBranchAddress("cluster_type", &cluster_type);
    IT->SetBranchAddress("eventID", &eventID);
    IT->SetBranchAddress("cluster_eta", &cluster_eta);

    Float_t cluster_R_OT;
    Float_t cluster_phi_OT;
    Float_t cluster_type_OT;
    Float_t eventID_OT;
    Float_t cluster_eta_OT;

    // Set branch addresses
    OT->SetBranchAddress("cluster_R", &cluster_R_OT);
    OT->SetBranchAddress("cluster_phi", &cluster_phi_OT);
    OT->SetBranchAddress("cluster_type", &cluster_type_OT);
    OT->SetBranchAddress("eventID", &eventID_OT);
    OT->SetBranchAddress("cluster_eta", &cluster_eta_OT);

    // Loop through tree entries and fill the histograms
    Long64_t numEntries = IT->GetEntries();
    for (Long64_t entry = 0; entry < numEntries; ++entry) {
        IT->GetEntry(entry);
        if (eventID == 4) {
            if (cluster_type == 1) {
                rHistType1->Fill(cluster_R);
                phiHistType1->Fill(cluster_phi);
                etaHistType1->Fill(cluster_eta);
            } else if (cluster_type == 2) {
                rHistType2->Fill(cluster_R);
                phiHistType2->Fill(cluster_phi);
                etaHistType2->Fill(cluster_eta);
            }
        }
    }

    numEntries = OT->GetEntries();
    for (Long64_t entry = 0; entry < numEntries; ++entry) {
        OT->GetEntry(entry);
        if (eventID_OT == 4) {
            if (cluster_type_OT == 1) {
                rHistType1->Fill(cluster_R_OT);
                phiHistType1->Fill(cluster_phi_OT);
                etaHistType1->Fill(cluster_eta_OT);
            } else if (cluster_type_OT == 2) {
                rHistType2->Fill(cluster_R_OT);
                phiHistType2->Fill(cluster_phi_OT);
                etaHistType2->Fill(cluster_eta_OT);
            }
        }
    }

    // Draw the histograms in the sub-pads
    canvas->cd(1);
    rHistType1->Draw();
    rHistType2->Draw("SAME");
    canvas->cd(1)->SetTitle("r");

    canvas->cd(2);
    phiHistType1->Draw();
    phiHistType2->Draw("SAME");
    canvas->cd(2)->SetTitle("phi");

    canvas->cd(3);
    etaHistType1->Draw();
    etaHistType2->Draw("SAME");
    canvas->cd(3)->SetTitle("eta");

    // Customize plot axis labels and titles
    rHistType1->GetXaxis()->SetTitle("r");
    phiHistType1->GetXaxis()->SetTitle("phi");
    etaHistType1->GetXaxis()->SetTitle("eta");
    rHistType1->GetYaxis()->SetTitle("Frequency");

    gStyle->SetOptStat(0);


    // Draw the legend
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(rHistType1, "Type 1", "l");
    legend->AddEntry(rHistType2, "Type 2", "l");
    legend->Draw();

    canvas->Update();
    canvas->Modified();
    canvas->Update();
    
}
