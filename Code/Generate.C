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

//in the future 18
static const int NFile = 6;
static TRandom3 * myRNG = new TRandom3(23);
static TFile *files[NFile];
static TDirectory *dirIT_list[NFile];
static TDirectory *dirOT_list[NFile];
static TTree *IT_list[NFile];
static TTree *OT_list[NFile];
static int N_Events_list[NFile];
static vector<int> N_Events_list_id[NFile];
static string P_name[3]= {"1", "3", "10"};
static double *P_cum;
static int NIT;
static int NOT;
static int N_RandEvents;
double scale = 1.e3;
static float epsilon = 1.e-10;

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
int recursive_binary_search(const double arr[],double target, int low, int high) {
    if (low > high) {
        return low;
    }

    int mid = (low + high) / 2;

    if (arr[mid] == target) {
        return mid;
    } else if (arr[mid] < target) {
        return recursive_binary_search(arr, target, mid + 1, high);
    } else {
        return recursive_binary_search(arr, target, low, mid - 1);
    }
}

int recursive_binary_search(TTree* tree, float *field, float target, int low, int high) {
    if (low > high) {
        return low;
    }

    int mid = (low + high) / 2;
    tree->GetEntry(mid);
    if (*field == target) {
        return mid;
    } else if (*field < target) {
        return recursive_binary_search(tree, field, target, mid + 1, high);
    } else {
        return recursive_binary_search(tree, field, target, low, mid - 1);
    }
}

int Get_first_row_event(TTree *tree, float* field, int id_event_value)
{
    int id = recursive_binary_search(tree, field, id_event_value, 0, tree->GetEntries()-1);
    //loop back to find the first entry of the event
    if(id == 0) return id;
    tree -> GetEntry(id-1);
    while((int) (*field) == id_event_value && (id != 0)){
        id--;
        tree -> GetEntry(id-1);
    }
    return id;
}

pair<std::vector<Event>, std::vector<Event>> GetEventFromMia(TTree *IT, TTree *OT, int id_event_value, float pclass, int new_id_event)
{
    vector<Event> event_IT = {};
    vector<Event> event_OT = {};
    int count = 0;

    float x, y, z;
    float r, phi, eta;
    float id_event;
    float type;

    //cout << IT << endl;

    IT->SetBranchAddress("cluster_x", &x);
    IT->SetBranchAddress("cluster_y", &y);
    IT->SetBranchAddress("cluster_z", &z);
    IT->SetBranchAddress("cluster_R", &r);
    IT->SetBranchAddress("cluster_phi", &phi);
    IT->SetBranchAddress("cluster_eta", &eta);
    IT->SetBranchAddress("eventID", &id_event);
    IT->SetBranchAddress("cluster_type", &type);

    int first_row_event = Get_first_row_event(IT, &id_event, id_event_value);
    //cout << "Getting event: " << id_event_value << endl;
    //cout << "First Row: " << first_row_event << endl;

    // Loop over entries and find rows with the specified id_event value
    for (int i = first_row_event; i < IT->GetEntries(); i++)
    {
        IT->GetEntry(i);
        if (static_cast<int>(id_event) != id_event_value){
        //cout << "Last Row: " << i << endl;
            break;
        }

        if((int) (type)==1){
            event_IT.emplace_back(x, y, z, r, phi, eta, new_id_event, type, pclass);
        }
        else
            event_IT.emplace_back(x, y, z, r, phi, eta, new_id_event, type, -1.);

    }

    //cout << OT << endl;
    // OUT Tracker
    OT->SetBranchAddress("cluster_x", &x);
    OT->SetBranchAddress("cluster_y", &y);
    OT->SetBranchAddress("cluster_z", &z);
    OT->SetBranchAddress("cluster_R", &r);
    OT->SetBranchAddress("cluster_phi", &phi);
    OT->SetBranchAddress("cluster_eta", &eta);
    OT->SetBranchAddress("eventID", &id_event);
    OT->SetBranchAddress("cluster_type", &type);

    int first_row_event_OT = Get_first_row_event(OT, &id_event, id_event_value);
    //cout << "First Row: " << first_row_event_OT << endl;
    for (int i = first_row_event_OT; i < OT->GetEntries(); i++)
    {
        OT->GetEntry(i);
        if (static_cast<int>(id_event) != id_event_value)
            break;

        if((int) (type)==1){
            count++;
            event_OT.emplace_back(x, y, z, r, phi, eta, new_id_event, type, pclass);
        }
        else
            event_OT.emplace_back(x, y, z, r, phi, eta, new_id_event, type, -1.);
    }

    
    if(count<3){
        cout << "Problema: " << new_id_event << "  " << id_event << " " << first_row_event_OT << " " << id_event_value << endl;
        for (int i = first_row_event_OT; i < OT->GetEntries(); i++)
    {
        OT->GetEntry(i);
        cout << i << endl;
        cout << id_event << " " << id_event_value<<endl;
        if (static_cast<int>(id_event) != id_event_value)
            break;    
    }

    }
    return make_pair(event_IT, event_OT);
}

void ComputeCumulative(TTree *IT, TTree *OT){
    NIT = IT->GetEntries();
    NOT = OT->GetEntries();
    //cout << "NIT: " << NIT << endl;
    //cout << "NOT: " << NOT << endl;
    N_RandEvents = NIT + NOT;

    float x,y,z;


    //Inner Traker
    IT->SetBranchAddress("cluster_x", &x);
    IT->SetBranchAddress("cluster_y", &y);
    IT->SetBranchAddress("cluster_z", &z);

    P_cum = new double[N_RandEvents];

    IT->GetEntry(0);
    P_cum[0] = sqrt(x*x + y*y  )/scale;
    for (int i = 1; i < NIT; i++)
    {
        IT->GetEntry(i);
        P_cum[i] = P_cum[i-1] + sqrt(x*x + y*y  )/scale;
    }
    
    //Outer Traker
    OT->SetBranchAddress("cluster_x", &x);
    OT->SetBranchAddress("cluster_y", &y);
    OT->SetBranchAddress("cluster_z", &z);

    for (int i = 0; i < NOT; i++)
    {
        OT->GetEntry(i);
        P_cum[NIT+i] = P_cum[NIT+i-1] + sqrt(x*x + y*y  )/scale;
    }
}

pair<std::vector<Event>, std::vector<Event>> GetBackgroundFromMia(TTree *IT, TTree *OT, int new_id_event, float bkg_rate = 50)
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
    //cout << "# clusters: " << N_gen << endl;

    for (int i = 1; i <= N_gen; i++)
    {   
        double p_i = myRNG->Uniform(P_cum[N_RandEvents-1]);
        //find the index corrisponding to p_i inside the Cumulative probability array
        int id = recursive_binary_search(P_cum, p_i, 0, N_RandEvents-1);
        //cout << "id: " << id;
        //extracting that hit and saving it inside the event vector
        if(id < NIT){
            //cout << " INNER" << endl;
            IT->GetEntry(id);
            event_IT.emplace_back(x, y, z, r, phi, eta, new_id_event, type, pclass);
        }
        else{
            //cout << " OUTER" << endl;
            OT->GetEntry(id - NIT);
            event_OT.emplace_back(x_OT, y_OT, z_OT, r_OT, phi_OT, eta_OT, new_id_event, type, pclass);
        }
    }
    //cout << "Background generated" << endl;
    return make_pair(event_IT, event_OT);
}

void GenerateRootFromMia(int N_events = 100000, string outRoot="100k.root", float bkg_rate = 50,  string folder = "/home/ema/Documents/thesis/DATA/MuGun/", string file_name = "clusters_ntuple.root")
{   

    //momentaneamente j = 0 per gestire solo i file a 1GeV
    int combind = 0;
    for (int j=0; j < 3; j++)
    {   
        //open all root files and TTrees inside
        //momentaneamente solo 1
        for (int i = 0; i < 2; i++)
        {
            string rootInput;
            if(i%2==0) rootInput = folder + P_name[j] + "GeV/SingleParticleEta0p4/" + file_name;
            else rootInput = folder + P_name[j] + "GeV/SingleParticleEta0p4/plus/" + file_name;
            TFile *file = TFile::Open(rootInput.c_str(), "READ");
            if (!file || file->IsZombie())
            {
                cerr << "Error: Cannot open file " << rootInput << endl;
                return;
            }

            cout <<"Opening: " << rootInput << endl;

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

            IT_list[combind] = IT;
            OT_list[combind] = OT;
            
            float id_event;
            IT->SetBranchAddress("eventID", &id_event);
            int size = IT->GetEntries();
            IT->GetEntry(size-1);
            N_Events_list[combind] = (int) id_event;
            for(int i = 0; i < size; i++){
                IT->GetEntry(i);
                N_Events_list_id[combind].push_back(id_event);
            }
            combind++;
        }   
    }
    
    //opening output file
    TFile* out = new TFile(outRoot.c_str(), "RECREATE");

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
    
    //Computing the cumulative probability
    ComputeCumulative(IT_list[0], OT_list[0]);
    
    //loop on the number of events
    for (int i = 0; i < N_events; i++)
    {
        if(i%(N_events/10)==0)
            cout << i/(N_events/10)*10 << "%" <<endl;
        //cout << i << endl;
        //generate background
        //cout << "Getting the background" << endl;
        pair <vector<Event>, vector<Event>> event = GetBackgroundFromMia(IT_list[0], OT_list[0], i+1, bkg_rate);
        vector<Event> event_IT = event.first;
        vector<Event> event_OT = event.second;
        
        //cout << "Gen BKG" << endl;
        //generate only background or signal with 50% of probability
        bool signal = (myRNG->Uniform())<0.5;
        if (signal){
            //add signal to the event
            //select random an event from a random file
            int ID_file = (int) (myRNG->Uniform(NFile-epsilon));
            //int ID_file;
            /*
             do{
                ID_file =(int) (myRNG->Uniform(NFile-epsilon));

            }while(N_Events_list_id[ID_file].size()==0);
            */
           
            
            int ID_event =(int) (myRNG->Uniform(N_Events_list[ID_file]-epsilon))+1;
            //int ID_event_id = (int) (myRNG)->Uniform(N_Events_list_id[ID_file].size()-epsilon);
            //int ID_event = N_Events_list_id[ID_file][ID_event_id];
            //N_Events_list_id[ID_file].erase(N_Events_list_id[ID_file].begin()+ID_event_id);
            float pclass = ID_file;
            //cout << ID_file << " " << ID_event << endl;
            //cout << IT_list[ID_file] << " " << OT_list[ID_file] << endl;

            pair <vector<Event>, vector<Event>> event_sig = GetEventFromMia(IT_list[ID_file], OT_list[ID_file], ID_event, pclass, i+1);
            vector<Event> event_IT_sig = event_sig.first;
            vector<Event> event_OT_sig = event_sig.second;


            //merge vectors
            event_IT.insert(event_IT.end(), event_IT_sig.begin(), event_IT_sig.end());
            event_OT.insert(event_OT.end(), event_OT_sig.begin(), event_OT_sig.end());
        }
        
        //write the root file

        int IT_size = event_IT.size();
        int OT_size = event_OT.size();

        //cout << "IT_size " << IT_size << endl;
        //cout << "OT_size " << OT_size << endl;

        //cout << "Write IT" << endl;
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

        //cout << "Write OT" << endl;
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

    out->Close();
    
}
