#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TTree.h"

void createRootFileFromText(const char* inputFileName, const char* outputFileName) {
    
    // Apri il file di testo in lettura
    std::ifstream inputFile(inputFileName);
    if (!inputFile.is_open()) {
        std::cerr << "Error: Could not open input file " << inputFileName << std::endl;
        return;
    }

    // Crea un nuovo file ROOT in modalità scrittura
    TFile* outputFile = TFile::Open(outputFileName, "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: Could not create output ROOT file " << outputFileName << std::endl;
        return;
    }

    // Crea l'albero
    TTree* tree = new TTree("tree", "HitsTree");

    double r, phi;
    int id;

    // Collega le variabili alle foglie (leaf) dell'albero
    tree->Branch("r", &r);
    tree->Branch("phi", &phi);
    tree->Branch("id", &id);

    while (inputFile >> r >> phi >> id) {
        tree->Fill();
    }

    // Scrivi l'albero nel file ROOT
    tree->Write();

    // Chiudi il file ROOT
    outputFile->Close();
}
