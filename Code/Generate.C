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

static const int NFile = 18;
static TRandom3 * myRNG = new TRandom3(23);
static TFile *files[NFile];
static TTree *IT_list[NFile];
static TTree *OT_list[NFile];
static string P_name[3]= {"1", "3", "10"};
static double *P_cum;
long long int NIT;
static long int NOT;
static long int N_RandEvents;
double scale = 1.e3;

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

//Recursive search algorithm
long int recursive_binary_search(const double arr[],double target, long int low, long int high) {
    if (low > high) {
        return low;
    }

    long int mid = (low + high) / 2;

    if (arr[mid] == target) {
        return mid;
    } else if (arr[mid] < target) {
        return recursive_binary_search(arr, target, mid + 1, high);
    } else {
        return recursive_binary_search(arr, target, low, mid - 1);
    }
}

long int recursive_binary_search(TTree* tree, float *field, float target, long int low, long int high) {
    if (low > high) {
        return low;
    }

    long int mid = (low + high) / 2;
    tree->GetEntry(mid);
    if (*field == target) {
        return mid;
    } else if (*field < target) {
        return recursive_binary_search(tree, field, target, mid + 1, high);
    } else {
        return recursive_binary_search(tree, field, target, low, mid - 1);
    }
}

long int Get_first_row_event(TTree *tree, float* field, long int id_event_value)
{
    long int id = recursive_binary_search(tree, field, id_event_value, 0, tree->GetEntries());
    //loop back to find the first entry of the event
    if(id == 0) return id;
    tree -> GetEntry(id-1);
    while((long int) field == id && id != 0){
        id--;
        tree -> GetEntry(id-1);
    }
    return id;
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

    long int first_row_event = Get_first_row_event(IT, &id_event, id_event_value);

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

    long int first_row_event_OT = Get_first_row_event(OT, &id_event, id_event_value);
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

    P_cum = new double[N_RandEvents];

    IT->GetEntry(0);
    P_cum[0] = sqrt(x*x + y*y + z*z)/scale;
    for (long int i = 1; i < NIT; i++)
    {
        IT->GetEntry(i);
        P_cum[i] = P_cum[i-1] + sqrt(x*x + y*y + z*z)/scale;
    }
    
    //Outer Traker
    OT->SetBranchAddress("cluster_x", &x);
    OT->SetBranchAddress("cluster_y", &y);
    OT->SetBranchAddress("cluster_z", &z);

    for (long int i = 0; i < NOT; i++)
    {
        OT->GetEntry(i);
        P_cum[NIT+i] = P_cum[NIT+i-1] + sqrt(x*x + y*y + z*z)/scale;
    }
}

pair<std::vector<Event>, std::vector<Event>> GetBackgroundFromMia(TTree *IT, TTree *OT, long int new_id_event, float bkg_rate = 50)
{
    vector<Event> event_IT = {};
    vector<Event> event_OT = {};


    float x, y, z;
    float r, phi, eta;
    float id_event;

    float x_OT, y_OT, z_OT;
    float r_OT, phi_OT, eta_OT;
    float id_event_OT;

    //default value for background
    float pclass = -1;
    float type = 2;

    IT->SetBranchAddress("cluster_x", &x);
    IT->SetBranchAddress("cluster_y", &y);
    IT->SetBranchAddress("cluster_z", &z);
    IT->SetBranchAddress("cluster_R", &r);
    IT->SetBranchAddress("cluster_phi", &phi);
    IT->SetBranchAddress("cluster_eta", &eta);

    OT->SetBranchAddress("cluster_x", &x_OT);
    OT->SetBranchAddress("cluster_y", &y_OT);
    OT->SetBranchAddress("cluster_z", &z_OT);
    OT->SetBranchAddress("cluster_R", &r_OT);
    OT->SetBranchAddress("cluster_phi", &phi_OT);
    OT->SetBranchAddress("cluster_eta", &eta_OT);

 
    int N_gen = myRNG->Poisson(bkg_rate);	
    cout << "# clusters: " << N_gen << endl;

    for (int i = 1; i <= N_gen; i++)
    {   
        double p_i = myRNG->Uniform(P_cum[N_RandEvents-1]);
        //find the index corrisponding to p_i inside the Cumulative probability array
        int id = recursive_binary_search(P_cum, p_i, 0, N_RandEvents-1);
        //extracting that hit and saving it inside the event vector
        if(id < NIT){
            IT->GetEntry(id);
            event_IT.emplace_back(x, y, z, r, phi, eta, new_id_event, type, pclass);
        }
        else{
            OT->GetEntry(id - NIT);
            event_OT.emplace_back(x_OT, y_OT, z_OT, r_OT, phi_OT, eta_OT, new_id_event, type, pclass);
            if(y_OT < 2000) cout << "PROBLEMA: " << id << " " << y_OT << " " << new_id_event << endl;
            if(r_OT>2000) cout << "PROBLEMA: " << id << " " << r_OT << " " << new_id_event << endl; 
        }
    }
    return make_pair(event_IT, event_OT);
}

void test(int nEvents = 100, string rootInput="/Users/Fabio/Desktop/DATA/MuGun/1GeV/SingleParticleEta0p4/clusters_ntuple_0.root"){
    
    TFile* file = TFile::Open(rootInput.c_str());
    TDirectoryFile *dirIT = dynamic_cast<TDirectoryFile *>(file->Get("clusterValidIT"));
    TDirectoryFile *dirOT = dynamic_cast<TDirectoryFile *>(file->Get("clusterValidOT"));

    // Access the "tree" in the "clusterValidIT" directory
    TTree *IT = dynamic_cast<TTree *>(dirIT->Get("tree"));
    TTree *OT = dynamic_cast<TTree *>(dirOT->Get("tree"));

    IT->SetMaxVirtualSize(250000000);
    IT->LoadBaskets();

    OT->SetMaxVirtualSize(250000000);
    OT->LoadBaskets();

    ComputeCumulative(IT, OT);

    for (int i = 0; i < N_RandEvents; i++)
    {
        cout << i <<  " " << P_cum[i] << endl;
    }
    
    //initialize the output file
    TFile* out = new TFile("bkg.root", "RECREATE");

    TDirectory *dirIT_out = out->mkdir("clusterValidIT");
    TDirectory *dirOT_out = out->mkdir("clusterValidOT");    
    
    dirIT_out->cd();
    TTree *IT_out = new TTree("tree", "RECREATE");

    float x_IT, y_IT, z_IT;
    float r_IT, phi_IT, eta_IT;
    float id_event_IT;
    float type_IT;
    float pclass_IT;

    IT_out->Branch("cluster_x", &x_IT);
    IT_out->Branch("cluster_y", &y_IT);
    IT_out->Branch("cluster_z", &z_IT);
    IT_out->Branch("cluster_R", &r_IT);
    IT_out->Branch("cluster_phi", &phi_IT);
    IT_out->Branch("cluster_eta", &eta_IT);
    IT_out->Branch("eventID", &id_event_IT);
    IT_out->Branch("cluster_type", &type_IT);
    IT_out->Branch("pclass", &pclass_IT);

    dirOT_out->cd();
    TTree *OT_out = new TTree("tree", "RECREATE");

    float x_OT, y_OT, z_OT;
    float r_OT, phi_OT, eta_OT;
    float id_event_OT;
    float type_OT;
    float pclass_OT;

    OT_out->Branch("cluster_x", &x_OT);
    OT_out->Branch("cluster_y", &y_OT);
    OT_out->Branch("cluster_z", &z_OT);
    OT_out->Branch("cluster_R", &r_OT);
    OT_out->Branch("cluster_phi", &phi_OT);
    OT_out->Branch("cluster_eta", &eta_OT);
    OT_out->Branch("eventID", &id_event_OT);
    OT_out->Branch("cluster_type", &type_OT);
    OT_out->Branch("pclass", &pclass_OT);
    
    for (int i = 0; i < nEvents; i++)
    {
        cout << "Getting the background" << endl;
        pair <vector<Event>, vector<Event>> event = GetBackgroundFromMia(IT, OT, i, 50);
        vector<Event> event_IT = event.first;
        vector<Event> event_OT = event.second;

        int IT_size = event_IT.size();
        int OT_size = event_OT.size();

        cout << "IT_size " << IT_size << endl;
        cout << "OT_size " << OT_size << endl;


        cout << "Write IT" << endl;
        dirIT_out->cd();
        for(int j = 0; j < IT_size; j++){
                
            x_IT = event_IT[j].x;
            y_IT = event_IT[j].y;
            z_IT = event_IT[j].z;
            r_IT = event_IT[j].r;
            phi_IT = event_IT[j].phi;
            eta_IT = event_IT[j].eta;
            id_event_IT = event_IT[j].id_event;
            type_IT = event_IT[j].type;
            pclass_IT = event_IT[j].pclass;

            IT_out->Fill();
        } 

        cout << "Write OT" << endl;
        dirOT_out->cd();

        for(int j = 0; j < OT_size; j++){
            x_OT = event_OT[j].x;
            y_OT = event_OT[j].y;
            z_OT = event_OT[j].z;
            r_OT = event_OT[j].r;
            phi_OT = event_OT[j].phi;
            eta_OT = event_OT[j].eta;
            id_event_OT = event_OT[j].id_event;
            type_OT = event_OT[j].type;
            pclass_OT = event_OT[j].pclass;

            OT_out->Fill();
        }
    }
    dirIT_out->cd();
    IT_out->Write();
    
    dirOT_out->cd();
    OT_out->Write();

    file->Close();
    out->Close();

}

void GenerateRootFromMia(string folder = "/Users/Fabio/Desktop/DATA/MuGun/", string file_name = "clusters_ntuple_")
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
    
    //GetEventFromMia();
}
