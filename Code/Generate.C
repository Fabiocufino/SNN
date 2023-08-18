#include "TH2.h"
#include "TH1.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "Riostream.h"

#include <math.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>

static TRandom3 * myRNG = new TRandom3(23);
static TFile *files[NFile];
static TTree *IT_list[NFile];
static TTree *OT_list[NFile];
static string P_name[3]= {"1", "3", "10"};
static float P_cum[];
long long int NIT;
static long int NOT;
static long int N_RandEvents;

//myRNG->Uniform()

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

long int Get_first_row_event(long int id_event_value)
{
    return 1;
}

vector<Event> GetEventFromMia(TTree *IT, TTree *OT, long int id_event_value, float pclass, long int new_id_event)
{
    vector<Event> event = {};

    float x, y, z;
    float r, phi, eta;
    float id_event;
    float type;

    IT->SetBranchAddress("cluster_x", &x);
    IT->SetBranchAddress("cluster_y", &y);
    IT->SetBranchAddress("cluster_z", &z);
    IT->SetBranchAddress("cluster_R", &r);
    IT->SetBranchAddress("cluster_phi", &phi);
    IT->SetBranchAddress("cluster_eta", &eta);
    IT->SetBranchAddress("eventID", &id_event);
    IT->SetBranchAddress("cluster_type", &type);

    long int first_row_event = Get_first_row_event(id_event_value);

    // Loop over entries and find rows with the specified id_event value
    for (long int i = first_row_event; i < IT->GetEntries(); ++i)
    {
        IT->GetEntry(i);
        if (static_cast<long int>(id_event) != id_event_value)
            break;

        event.emplace_back(x, y, z, r, phi, eta, new_id_event, type, pclass);
    }

    // OUT Tracker
    OT->SetBranchAddress("cluster_x", &x);
    OT->SetBranchAddress("cluster_y", &y);
    OT->SetBranchAddress("cluster_z", &z);
    OT->SetBranchAddress("cluster_R", &r);
    OT->SetBranchAddress("cluster_phi", &phi);
    OT->SetBranchAddress("cluster_eta", &eta);
    OT->SetBranchAddress("eventID", &id_event);
    OT->SetBranchAddress("cluster_type", &type);

    long int first_row_event_OT = Get_first_row_event(id_event_value);
    for (long int i = first_row_event_OT; i < OT->GetEntries(); ++i)
    {
        OT->GetEntry(i);
        if (static_cast<long int>(id_event) != id_event_value)
            break;

        event.emplace_back(x, y, z, r, phi, eta, new_id_event, type, pclass);
    }

    return event;
}

void ComputeCumulative(TTree *IT, TTree *OT){
    NIT = IT->GetEntries();
    NOT = OT->GetEntries();
    N_RandEvents = NIT + NOT;

    float x,y,z;


    //Inner Traker
    IT->SetBranchAddress("cluster_x", &x);
    IT->SetBranchAddress("cluster_y", &y);
    IT->SetBranchAddress("cluster_z", &z);

    P_cum = new float[Nevents];

    IT->GetEntry(0);
    P_cum[0] = x**2 + y**2 + z**2;
    for (long int i = 1; i < NIT; i++)
    {
        IT->GetEntry(i);
        P_cum[i] = P_cum[i-1] + (x**2 + y**2 + z**2);
    }
    
    //Outer Traker
    OT->SetBranchAddress("cluster_x", &x);
    OT->SetBranchAddress("cluster_y", &y);
    OT->SetBranchAddress("cluster_z", &z);

    for (long int i = 0; i < NOT; i++)
    {
        OT->GetEntry(i);
        P_cum[i] = P_cum[i-1] + (x**2 + y**2 + z**2);
    }
}

vector<Event> GetBackgroundFromMia(TTree *IT, TTree *OT, long int new_id_event, float bkg_rate = 50)
{
    vector<Event> event = {};

    float x, y, z;
    float r, phi, eta;
    float id_event;
    float type;

    float x_, y, z;
    float r, phi, eta;
    float id_event;
    float type;

    IT->SetBranchAddress("cluster_x", &x);
    IT->SetBranchAddress("cluster_y", &y);
    IT->SetBranchAddress("cluster_z", &z);
    IT->SetBranchAddress("cluster_R", &r);
    IT->SetBranchAddress("cluster_phi", &phi);
    IT->SetBranchAddress("cluster_eta", &eta);
    IT->SetBranchAddress("eventID", &id_event);
    IT->SetBranchAddress("cluster_type", &type);

    OT->SetBranchAddress("cluster_x", &x);
    OT->SetBranchAddress("cluster_y", &y);
    OT->SetBranchAddress("cluster_z", &z);
    OT->SetBranchAddress("cluster_R", &r);
    OT->SetBranchAddress("cluster_phi", &phi);
    OT->SetBranchAddress("cluster_eta", &eta);
    OT->SetBranchAddress("eventID", &id_event);
    OT->SetBranchAddress("cluster_type", &type);


    
    int N_gen = TRandom::Poisson(bkg_rate);	

    for (int i = 0; i < N_gen; i++)
    {
        int p_i = TRandom::Uniform(P_cum[N_RandEvents]);
    }
    
    
    return event;
}


void GenerateRootFromMia(const int NFile = 18, string folder = "/Users/Fabio/Desktop/DATA/MuGun/", string file_name = "clusters_ntuple_")
{
    for (int j=0; j < 3; j++)
    {
        for (int i = 0; i < 6; i++)
        {
            string rootInput = folder + P_name[j] + "GeV/SingleParticleEta0p4/" + file_name + to_string(i)+".root ";
            TFile *file = TFile::Open(rootInput.c_str(), "READ");
            if (!file || file->IsZombie())
            {
                cerr << "Error: Cannot open file " << rootInput << endl;
                return;
            }

            // Access the "clusterValidIT" directory
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
            TTree *OT = dynamic_cast<TTree *>(dirOT->Get("tree;1"));

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


            IT->SetMaxVirtualSize(250000000);
            IT->LoadBaskets();

            OT->SetMaxVirtualSize(250000000);
            OT->LoadBaskets();

            IT_list[i] = IT;
            OT_list[i] = OT;

            float id_event;
            IT->SetBranchAddress("eventID", &id_event);
            IT->GetEntry(IT->GetEntries()-1);


            
            
        }
    }
    //Computo la funzione cumulativa
    ComputeCumulative(IT_list[0], OT_list[0]);
    
    GetEventFromMia();
}


void 