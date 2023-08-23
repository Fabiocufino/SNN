#include "Snnt_constants.h"
#include <vector>
#include <TGraph.h>
#include <TCanvas.h>
#include <TArc.h>
#include <TMarker.h>

using namespace std;

struct Event
{
    float x, y, z;
    float r, phi, eta;
    float id_event;
    float type;
    float pclass;

    Event(float x_, float y_, float z_, float r_, float phi_, float eta_, float id_event_, float type_, float pclass_)
        : x(x_), y(y_), z(z_), r(r_), phi(phi_), eta(eta_), id_event(id_event_), type(type_), pclass(pclass_)
    {
    }
};

void createPNGFrames(const char *file_name)
{

    string rootInput;
    rootInput = file_name;

    TFile *file = TFile::Open(rootInput.c_str(), "READ");

    if (!file || file->IsZombie())
    {
        cerr << "Error: Cannot open file " << rootInput << endl;
        return;
    }

    cout << "Opening: " << rootInput << endl;

    // Access directories
    TDirectoryFile *dirIT = dynamic_cast<TDirectoryFile *>(file->Get("clusterValidIT"));
    TDirectoryFile *dirOT = dynamic_cast<TDirectoryFile *>(file->Get("clusterValidOT"));

    if (!dirIT)
    {
        cerr << "Error: Cannot access directory clusterValidIT" << endl;
        file->Close();
        return;
    }

    if (!dirOT)
    {
        cerr << "Error: Cannot access directory clusterValidOT" << endl;
        file->Close();
        return;
    }

    // Access the "tree" in the "clusterValidIT" directory
    TTree *IT = dynamic_cast<TTree *>(dirIT->Get("tree"));
    TTree *OT = dynamic_cast<TTree *>(dirOT->Get("tree"));

    if (!IT)
    {
        cerr << "Error: Cannot access tree in clusterValidIT" << endl;
        file->Close();
        return;
    }

    if (!OT)
    {
        cerr << "Error: Cannot access tree in clusterValidOT" << endl;
        file->Close();
        return;
    }

    //------------
    vector<Event> event = {};
    vector<Event> event_BKG = {};

    float x, y, z;
    float r, phi, eta;
    float id_event;

    float x_OT, y_OT, z_OT;
    float r_OT, phi_OT, eta_OT;
    float id_event_OT;

    // default value for background
    float pclass = -10;
    float type, type_OT;

    IT->SetBranchAddress("cluster_x", &x);
    IT->SetBranchAddress("cluster_y", &y);
    IT->SetBranchAddress("cluster_z", &z);
    IT->SetBranchAddress("cluster_R", &r);
    IT->SetBranchAddress("cluster_phi", &phi);
    IT->SetBranchAddress("cluster_eta", &eta);
    IT->SetBranchAddress("eventID", &id_event);
    IT->SetBranchAddress("cluster_type", &type);

    OT->SetBranchAddress("cluster_x", &x_OT);
    OT->SetBranchAddress("cluster_y", &y_OT);
    OT->SetBranchAddress("cluster_z", &z_OT);
    OT->SetBranchAddress("cluster_R", &r_OT);
    OT->SetBranchAddress("cluster_phi", &phi_OT);
    OT->SetBranchAddress("cluster_eta", &eta_OT);
    OT->SetBranchAddress("eventID", &id_event_OT);
    OT->SetBranchAddress("cluster_type", &type_OT);

    for (int i = 0; i <= IT->GetEntries()-1; i++)
    {
        IT->GetEntry(i);
        event.emplace_back(x, y, z, r, phi, eta, id_event, type, pclass);
    }

    for (int i = 0; i <= OT->GetEntries()-1; i++)
    {
        OT->GetEntry(i);
        event.emplace_back(x_OT, y_OT, z_OT, r_OT, phi_OT, eta_OT, id_event_OT, type_OT, pclass);
    }

    //------------
    TVectorF vec_x(event.size());
    TVectorF vec_y(event.size());

    for (int i = 0; i <= event.size(); i++)
    {
        if (event[i].type == 1)
        {
            vec_x(i) = event[i].r * sin(event[i].phi);
            vec_y(i) = event[i].r * cos(event[i].phi);
        }
        if (event[i].type == 2)
        {
           //vec_x(i) = event_BKG[q].x;
        }
    }

    //------------

    TCanvas *canvas = new TCanvas("canvas", "Animated GIF", 800, 800);
    //canvas->SetGrid();
    canvas->Range(-1200, -1200, 1200, 1200);

    TGraph *graph = new TGraph(vec_x, vec_y);
   // graph->SetTitle(Form("Frame %d", frameCount + 1));
    graph->GetXaxis()->SetTitle("X");
    graph->GetYaxis()->SetTitle("Y");
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.1);

    graph->Draw("P");

    // TGraph *bkg = new TGraph();
    // bkg->SetTitle(Form("Frame %d", frameCount + 1));
    // bkg->GetXaxis()->SetTitle("X");
    // bkg->GetYaxis()->SetTitle("Y");
    // bkg->SetMarkerStyle(20);
    // bkg->SetMarkerSize(2.2);

    // while (inputFile >> x >> y >> trackId)
    // {
    //     if (lineCount % 10 == 0)
    //     {
    //         canvas->Clear();

    //         // Draw all points in blue
    //         graph->SetMarkerColor(kBlue);
    //         graph->Draw("P");

    //         // Draw background points in red
    //         bkg->SetMarkerColor(kGreen);
    //         bkg->Draw("P");

    //         // Draw arcs for background points (y2=0)
    //         for (int i = 0; i < graph->GetN(); i++)
    //         {
    //             double radius = strip_r[i];
    //             TArc *arc = new TArc(0, 0, radius);
    //             arc->SetLineColor(kRed);
    //             arc->SetFillStyle(0);
    //             arc->Draw();
    //         }

    //         canvas->Update();
    //         canvas->Print(Form("frame%d.png", frameCount)); // Save the frame as an image PNG
    //         graph->Clear();
    //         frameCount++;
    //     }
    //     if (trackId == 0)
    //     {
    //         bkg->SetPoint(lineCount % 10, x, y);
    //     }
    //     else
    //         graph->SetPoint(lineCount % 10, x, y);
    //     lineCount++;
    // }

    // inputFile.close();
    // delete canvas;
}

void createAnimatedGIF()
{
    // Utilize ImageMagick to convert the PNG images into an animated GIF
    system("convert -delay 70 frame*.png output.gif");

    // Remove the temporary PNG images
    system("rm frame*.png");
}

void plotAnimatedGIF()
{
    createPNGFrames("plot.txt");
    createAnimatedGIF();
}
