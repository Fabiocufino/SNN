#include "Snnt_constants.h"
#include <vector>
#include <TGraph.h>
#include <TCanvas.h>
#include <TArc.h>
#include <TMarker.h>

using namespace std;

void createPNGFrames(const char* filename) {
    TCanvas *canvas = new TCanvas("canvas", "Animated GIF", 800, 800);
    canvas->SetGrid();
    canvas->Range(-120, -120, 120, 120);

    ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        cout << "Error: Could not open file " << filename << endl;
        return;
    }

    double x, y, trackId;
    int lineCount = 0;
    int frameCount = 0;

    TGraph *graph = new TGraph();
    graph->SetTitle(Form("Frame %d", frameCount + 1));
    graph->GetXaxis()->SetTitle("X");
    graph->GetYaxis()->SetTitle("Y");
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1.8);


    TGraph *bkg = new TGraph();
    bkg->SetTitle(Form("Frame %d", frameCount + 1));
    bkg->GetXaxis()->SetTitle("X");
    bkg->GetYaxis()->SetTitle("Y");
    bkg->SetMarkerStyle(20);
    bkg->SetMarkerSize(2.2);

    while (inputFile >> x >> y >> trackId) {
        if (lineCount % 10 == 0) {
            canvas->Clear();

            // Draw all points in blue
            graph->SetMarkerColor(kBlue);
            graph->Draw("P");

            // Draw background points in red
            bkg->SetMarkerColor(kGreen);
            bkg->Draw("P");


            // Draw arcs for background points (y2=0)
            for (int i = 0; i < graph->GetN(); i++) {
                double radius = strip_r[i];
                TArc *arc = new TArc(0, 0, radius);
                arc->SetLineColor(kRed);
                arc->SetFillStyle(0);
                arc->Draw();
            }

            canvas->Update();
            canvas->Print(Form("frame%d.png", frameCount)); // Save the frame as an image PNG
            graph->Clear();
            frameCount++;
        }
        if (trackId == 0){
            bkg->SetPoint(lineCount % 10, x, y);
        }
        else
            graph->SetPoint(lineCount % 10, x, y);
        lineCount++;
    }

    inputFile.close();
    delete canvas;
}

void createAnimatedGIF() {
    // Utilize ImageMagick to convert the PNG images into an animated GIF
    system("convert -delay 70 frame*.png output.gif");

    // Remove the temporary PNG images
    system("rm frame*.png");
}

void plotAnimatedGIF() {
    createPNGFrames("plot.txt");
    createAnimatedGIF();
}
