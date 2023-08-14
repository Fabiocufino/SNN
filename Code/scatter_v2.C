#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TMath.h"

using namespace std;

void draw_scatter_plot(const char* filename = "clusters_ntuple.root") {
	int count_hit = 0;
	int count_bg = 0;
	
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
    // Create a canvas for the scatter plot
    TCanvas* canvas = new TCanvas("scatter_plot", "Scatter Plot", 800, 600);

    // Create TGraph for red points (cluster_type == 1)
    TGraph* redGraph = new TGraph();
    redGraph->SetMarkerStyle(20);
    redGraph->SetMarkerColor(kRed);

    // Create TGraph for blue points (cluster_type == 2)
    TGraph* blueGraph = new TGraph();
    blueGraph->SetMarkerStyle(20);
    blueGraph->SetMarkerColor(kBlue);

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

    // Loop through tree entries and fill the appropriate TGraph
    Long64_t numEntries = IT->GetEntries();
    for (Long64_t entry = 0; entry < numEntries; ++entry) {
        IT->GetEntry(entry);
        double theta = TMath::ATan(TMath::Exp(-cluster_eta_OT)) * 2.;
        Double_t y = cluster_R;
        Double_t x = cluster_R/TMath::Tan(theta);

        
        if(eventID == 4){
            
        cout << cluster_R << " " << cluster_phi << " " << cluster_eta << " " << cluster_type <<endl;
        if (cluster_type == 1) {
            redGraph->SetPoint(redGraph->GetN(), x, y);
            count_bg++;
        } else if (cluster_type == 2) {
            blueGraph->SetPoint(blueGraph->GetN(), x, y);
            count_hit++;
        }
        }
    }

    numEntries = OT->GetEntries();
    for (Long64_t entry = 0; entry < numEntries; ++entry) {
        OT->GetEntry(entry);
        double theta = TMath::ATan(TMath::Exp(-cluster_eta_OT)) * 2.;
        Double_t y = cluster_R_OT;
        Double_t x = cluster_R_OT/TMath::Tan(theta);

        if(eventID_OT == 4){
        cout << cluster_R_OT << " " << cluster_phi_OT << " " << cluster_eta_OT << " " << cluster_type_OT <<endl;
        if (cluster_type_OT == 1) {
            redGraph->SetPoint(redGraph->GetN(), x, y);
            count_bg++;
        } else if (cluster_type_OT == 2) {
            blueGraph->SetPoint(blueGraph->GetN(), x, y);
            count_hit++;
        }
        }
    }

    // Draw the TGraphs on the canvas
    canvas->cd();

	
    redGraph->Draw("AP");
    blueGraph->Draw("P SAME");

        redGraph -> GetYaxis()->SetRangeUser(0, 1300);
	redGraph -> GetXaxis()->SetRangeUser(-600, 600);
	
	blueGraph -> GetYaxis()->SetRangeUser(0, 1300);
	blueGraph -> GetXaxis()->SetRangeUser(-600, 600);

    // Customize plot axis labels and titles
    redGraph->GetYaxis()->SetTitle("r");
    redGraph->GetXaxis()->SetTitle("z");
    redGraph->SetTitle("Scatter Plot");

    // Draw the legend
    TLegend* legend = new TLegend(0.7, 0.8, 0.9, 0.9);
    legend->AddEntry(redGraph, "type 1", "p");
    legend->AddEntry(blueGraph, "type 2", "p");
    legend->Draw();

    canvas->Update();
    canvas->Modified();
    canvas->Update();
    
}
