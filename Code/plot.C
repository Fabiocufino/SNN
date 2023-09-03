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

void createPNGFrames(const char *file_name, int N_EV)
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

    for (int i = 0; i <= IT->GetEntries() - 1; i++)
    {
        IT->GetEntry(i);
        event.emplace_back(x, y, z, r, phi, eta, id_event, type, pclass);
    }

    for (int i = 0; i <= OT->GetEntries() - 1; i++)
    {
        OT->GetEntry(i);
        event.emplace_back(x_OT, y_OT, z_OT, r_OT, phi_OT, eta_OT, id_event_OT, type_OT, pclass);
    }

    //------------
    int frameCount = 0;

    for (int ev = 1; ev < N_EV; ev++)
    {
        TVectorF vec_x(400);
        TVectorF vec_y(400);

        TVectorF vec_x_BKG(400);
        TVectorF vec_y_BKG(400);

        int l = 0;
        for (int j = 0; j < event.size(); j++)
        {

            if (event[j].id_event == ev)
            {
                vec_x(l) = 0;
                vec_y(l) = 0;
                vec_x_BKG(l) = 0;
                vec_y_BKG(l) = 0;

                if (event[j].type == 1)
                {
                    vec_x(l) = event[j].r * sin(event[j].phi);
                    vec_y(l) = event[j].r * cos(event[j].phi);
                    l = l + 1;
                }
                if (event[j].type == 2)
                {

                    vec_x_BKG(l) = event[j].r * sin(event[j].phi);
                    vec_y_BKG(l) = event[j].r * cos(event[j].phi);
                    l = l + 1;
                }
            }
        }

        TCanvas *canvas = new TCanvas("canvas", "Animated GIF", 800, 800);
        // canvas->SetGrid();
        canvas->Range(-1200, -1200, 1200, 1200);
        TImage *img = TImage::Open("sfondo.jpg"); // put your own picture here
        img->Draw("x");

        TGraph *graph = new TGraph(vec_x, vec_y);

        graph->SetTitle(Form("Frame %d", frameCount + 1));
        graph->SetMarkerStyle(20);
        graph->SetMarkerColor(kBlue);
        graph->SetMarkerSize(0.7);
        graph->Draw("P");

        TGraph *bkg = new TGraph(vec_x_BKG, vec_y_BKG);

        bkg->SetTitle(Form("Frame %d", frameCount + 1));
        bkg->SetMarkerStyle(20);
        bkg->SetMarkerSize(0.7);
        bkg->SetMarkerColor(kRed);
        bkg->Draw("PSAME");

        TText *text = new TText(750, 1100, Form("Evento: %d", ev));
        text->SetTextSize(0.03);
        //text->Draw("SAME");
        // Create a transparent pad filling the full canvas
        TPad *p = new TPad("p", "p", 0, 0, 1, 1);
        p->SetFillStyle(4000);
        p->SetFrameFillStyle(4000);
        p->Draw("SAME");
        p->cd();


        canvas->Update();
        canvas->Print(Form("./img/frame%d.png", frameCount)); // Save the frame as an image PNG
        graph->Clear();
        frameCount++;
    }
    // if (trackId == 0)
    // {
    //     bkg->SetPoint(lineCount % 10, x, y);
    // }
    // else
    //     graph->SetPoint(lineCount % 10, x, y);
    // lineCount++;
    // delete canvas;
    // inputFile.close();
}

void createAnimatedGIF()
{
    // Utilize ImageMagick to convert the PNG images into an animated GIF
    system("convert -delay 70 frame*.png output.gif");

    // Remove the temporary PNG images
    system("rm frame*.png");
}

void plotAnimatedGIF(const char *file_name, int N_EV)
{
    createPNGFrames(file_name, N_EV);
    createAnimatedGIF();
}

void canvasimage(const char *file_name, int Nevent)
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

    for (int i = 0; i <= IT->GetEntries() - 1; i++)
    {
        if(id_event>Nevent) break;
        IT->GetEntry(i);
        event.emplace_back(x, y, z, r, phi, eta, id_event, type, pclass);
    }

    for (int i = 0; i <= OT->GetEntries() - 1; i++)
    {
        if(id_event_OT>Nevent) break;
        OT->GetEntry(i);
        event.emplace_back(x_OT, y_OT, z_OT, r_OT, phi_OT, eta_OT, id_event_OT, type_OT, pclass);
    }

    //------------
    int l = 0;

    TVectorF vec_x(event.size());
    TVectorF vec_y(event.size());

    TVectorF vec_x_BKG(event.size());
    TVectorF vec_y_BKG(event.size());

    for (int j = 0; j < event.size(); j++)
    {

        if (event[j].type == 1)
        {
            vec_x(l) = event[j].r * sin(event[j].phi);
            vec_y(l) = event[j].r * cos(event[j].phi);
            l = l + 1;
        }
        if (event[j].type == 2)
        {

            vec_x_BKG(l) = event[j].r * sin(event[j].phi);
            vec_y_BKG(l) = event[j].r * cos(event[j].phi);
            l = l + 1;
        }
    }
    TCanvas *canvas = new TCanvas("canvas", "Animated GIF", 800, 800);
    // canvas->SetGrid();
    canvas->Range(-1200, -1200, 1200, 1200);

    TGraph *graph = new TGraph(vec_x, vec_y);
    graph->SetTitle("Frame");
    graph->GetXaxis()->SetTitle("X");
    graph->GetYaxis()->SetTitle("Y");
    graph->SetMarkerStyle(24);
    // graph->SetMarkerColor(kBlue);
    graph->SetMarkerColorAlpha(kBlack, 0.95);

    graph->SetMarkerSize(0.1);

    graph->Draw("P");

    TGraph *bkg = new TGraph(vec_x_BKG, vec_y_BKG);

    // bkg->SetTitle(Form("Frame %d", frameCount + 1));
    bkg->SetTitle("Frame");

    bkg->GetXaxis()->SetTitle("X");
    bkg->GetYaxis()->SetTitle("Y");
    bkg->SetMarkerStyle(24);
    bkg->SetMarkerSize(0.1);
    // bkg->SetMarkerColor(kBlack);
    bkg->SetMarkerColorAlpha(kBlack, 0.95);
    bkg->Draw("PSAME");

    // canvas->Clear();

    // create a canvas with a picture (.png,.jpg,.gif, etc) in the background

    // TCanvas *c1 = new TCanvas("c1");
    // TImage *img = TImage::Open("sfondo.png"); // put your own picture here
    // img->Draw("x");

    // // Create a transparent pad filling the full canvas
    // TPad *p = new TPad("p", "p", 0, 0, 1, 1);
    // p->SetFillStyle(4000);
    // p->SetFrameFillStyle(4000);
    // p->Draw();
    // p->cd();
    // TH1F *h = new TH1F("h", "test", 100, -3, 3);
    // h->SetFillColor(kCyan);
    // h->FillRandom("gaus", 5000);
    // h->Draw();
}