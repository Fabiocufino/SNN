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
#include "Snnt_constants_eta.h"

#include <math.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

using namespace std;

// Constants and data used throughout the code
// -------------------------------------------
static double Occupancy;                                 // Probability of random hit firing
static int N_part;                                       // Number of generated particles in an event
static int KindOfSignal;                                 // 0 for normal single particle configs; 1 for V (+- pair); 2 for kinks (to be implemented)
static int StartLayer;                                   // Layer when first hit occurs (=0 for single particles, >=0 for V particles)
static double strip_phi[MaxStrips];               
static double First_angle;
static double Eff[MaxNeurons*MaxClasses];                // Efficiency of each neuron to signals of different classes
static double Efftot[MaxClasses];                        // Global efficiency to a different class
static double Weight[MaxNeurons][MaxStreams];            // Weight of synapse-neuron strength
static bool check_LTD[MaxNeurons][MaxStreams];           // checks to generate LTD after neuron discharge
static bool Void_weight[MaxNeurons][MaxStreams];         // These may be used to model disconnections
static bool bestVoid_weight[MaxNeurons][MaxStreams];
static double Weight_initial[MaxNeurons][MaxStreams];    // store to be able to return to initial conditions when optimizing
static double Delay[MaxNeurons][MaxStreams];             // Delay in incoming signals
static double bestDelay[MaxNeurons][MaxStreams];         // Opt delay in incoming signals
static vector<double> History_time[MaxNeurons];          // Time of signal events per each neuron
static vector<int> History_type[MaxNeurons];             // Type of signal
static vector<int> History_ID[MaxNeurons];               // ID of generating signal stream or neuron
static vector<double> Fire_time[MaxNeurons];             // Times of firing of each neuron
static int Neuron_layer[MaxNeurons];
static int N_neuronsL[2];                                // Number of neurons in layers 0 and 1
static int N_streams;
static int N_neurons;
static int N_classes;
static int N_events;
static int N_epochs;
static int NevPerEpoch;
static double ConnectedFraction_Input_L0;
static double ConnectedFraction_Input_L1;
static double ConnectedFraction_L0_L1;
static vector<double> PreSpike_Time;
static vector<int> PreSpike_Stream;
static vector<int> PreSpike_Signal;                      // 0 for background hit, 1 for signal hit, 2 for L1 neuron spike
static double tmax;                                      // t of max value for EPSP
static double Pmax_EPSP;                                 // maximum EPSP spike height
static double K;                                         // constant computed such that it sets the max of excitation spike at 1V
static bool update9;                                     // controls whether to optimize 7 network parameters
static bool updateDelays;                                // controls whether to optimize neuron delays
static bool updateConnections;                           // controls whether to optimize connections between streams and neurons
static bool anyHits            = true;                   // Whether to accept tracks with any number of hits <8 or not
static double Q_best;
static double SelL1_best;
static double Eff_best;
static double Acc_best;
static double T0_best;
static double T1_best;
static double A_best;
static double L1if_best;
static double K_best;
static double K1_best;
static double K2_best;
static double IEPC_best;
static double IPSPdf_best;
static int indfile;
static char progress[53]      = "[--10%--20%--30%--40%--50%--60%--70%--80%--90%-100%]"; // Progress bar
static long int ievent;
static long int NROOT = 100000;

// Sorting vector by column
bool sortcol(const vector<double>& v1, const vector<double>& v2)
{
    int id_col=1;
	return v1[id_col] < v2[id_col];
}

// New random number generator
// ---------------------------
static TRandom3 * myRNG = new TRandom3(23);


float getP(){
    if (pclass == 0 || pclass == 1)
        return 1.;
    if (pclass == 2 || pclass == 3)
        return 3.;
    return 10;
}

float getBeta(){
    return  myRNG->Uniform(2.*M_PI);
    //return  M_PI/2 - 0.1 - myRNG->Uniform(0.6);
}

short int getCharge(){
    return pow(-1.,pclass%2);
}


// Bisection method
// ---------------------------
float f(float r, float xc, float yc, float R, float phi){
    return r*r - 2*r * (xc*cos(phi) + yc*sin(phi)) + xc*xc + yc*yc - R*R;
}

float bisection(float a, float b, float s, float r, float xc, float yc, float R){
    float c;
    do{
        c = (b+a)*0.5;
        if((f(r,xc,yc,R,a) * f(r,xc,yc,R,c)) < 0)
            b = c;
        else
            a = c;

    }while((b-a) > s);
    return c;
}

// Define tracker geometry
// ---------------------------
void Define_tracker () {
    for (int is=0; is<N_strips; is++) {
        strip_phi[is] = pitch_rad*is;
    }
}

// Function that reads parameters from file
// ----------------------------------------
int Read_Parameters () {
    string Path  = "./MODE/SNNT/";
    ifstream tmpfile;
    indfile = -1;
    // Determine last available file number to read from, by attempting to open all files with same name and previous numbering
    do {
        if (indfile>-1) tmpfile.close();
        indfile++;
        stringstream tmpstring;
        tmpstring << "Params13_NL0=" << N_neuronsL[0] << "_NL1=" << N_neuronsL[1] << "_NCl=" << N_classes << "_" << indfile;
        string tmpfilename = Path + tmpstring.str() + ".txt";
        tmpfile.open(tmpfilename);
    } while (tmpfile.is_open());

    if (indfile==-1) {
        cout << "     Warning, no file to read parameters from. " << endl;
        return -1;
    }
    ifstream parfile;
    stringstream sstr;
    char num[40];
    sprintf (num, "NL0=%d_NL1=%d_NCl=%d_%d", N_neuronsL[0], N_neuronsL[1], N_classes, indfile-1); // we'll pick the last one in the list
    sstr << "Params13_";
    string nameparfile = Path  + sstr.str() + num + ".txt";
    parfile.open(nameparfile);
    double e;
    int ie;
    parfile >> ie;
    if (ie!=N_neurons) cout << "Warning, file " << nameparfile << " ie= " << ie << " N_neurons = " << N_neurons << " - input file not matching N_neurons" << endl;
    parfile >> ie;
    if (ie!=N_streams) cout << "Warning, file " << nameparfile << " ie= " << ie << " N_streams = " << N_streams << " - input file not matching N_streams" << endl;
    parfile >> e;
    Threshold[0] = e;
    parfile >> e; 
    Threshold[1] = e;
    parfile >> e;
    alpha = e;
    parfile >> e;
    L1inhibitfactor = e;
    parfile >> e;
    K = e;
    parfile >> e;
    K1 = e;
    parfile >> e;
    K2 = e;
    parfile >> e;
    IE_Pot_const = e;
    parfile >> e;
    IPSP_dt_dilation = e;
    for (int in=0; in<N_neurons; in++) {
        for (int is=0; is<N_streams; is++) {
            parfile >> e;
            Delay[in][is] = e;
        }
    }
    bool b;
    for (int in=0; in<N_neurons; in++) {
        for (int is=0; is<N_streams; is++) {
            parfile >> b;
            Void_weight[in][is] = b;
            if (in<N_neuronsL[0] && is>=N_InputStreams) Void_weight[in][is] = true; // just making sure
        }
    }
    parfile.close();
    return 0;
}

// Function that saves the layout data to file
// -------------------------------------------
void Write_Parameters () {

    string Path  = "./MODE/SNNT/";
    ifstream tmpfile;
    indfile = -1;
    // Determine first available file number to write, by attempting to open all files with same name and previous numbering
    do {
        if (indfile>-1) tmpfile.close();
        indfile++;
        stringstream tmpstring;
        tmpstring << "Params13_NL0=" << N_neuronsL[0] << "_NL1=" << N_neuronsL[1] << "_NCl=" << N_classes << "_" << indfile;
        string tmpfilename = Path + tmpstring.str() + ".txt";
        tmpfile.open(tmpfilename);
    } while (tmpfile.is_open());

    ofstream parfile;
    stringstream sstr;
    char num[40];
    sprintf (num,"NL0=%d_NL1=%d_NCl=%d_%d", N_neuronsL[0], N_neuronsL[1], N_classes, indfile);
    sstr << "Params13_";
    string nameparfile = Path  + sstr.str() + num + ".txt";
    parfile.open(nameparfile);
    parfile << N_neurons << " " << N_streams << endl;
    parfile << T0_best << " " << T1_best << " " << A_best << " " 
            << L1if_best << " " << K_best << " " << K1_best << " " << K2_best << " " << IEPC_best << " " << IPSP_dt_dilation << endl; 
    for (int in=0; in<N_neurons; in++) {
        for (int is=0; is<N_streams; is++) {
            parfile << bestDelay[in][is] << endl;
        }
    }
    for (int in=0; in<N_neurons; in++) {
        for (int is=0; is<N_streams; is++) {
            parfile << bestVoid_weight[in][is] << endl;
        }
    }

    // Also write optimization output
    parfile << Q_best << " " << " " << Eff_best << " " << Acc_best << " " << SelL1_best << endl;

    // Finally, write complete set of hyperparameters and settings
    parfile << "                       L0 neurons: " << N_neuronsL[0] << endl;
    parfile << "                       L1 neurons: " << N_neuronsL[1] << endl;
    parfile << "            Connected L0-L1 frac.: " << ConnectedFraction_L0_L1 << endl;
    parfile << "            Connected IN-L0 frac.: " << ConnectedFraction_Input_L0 << endl;
    parfile << "            Connected IN-L1 frac.: " << ConnectedFraction_Input_L1 << endl;
    parfile << "                  Noise per layer: " << Occupancy*N_strips << endl;
    parfile << "                    Track classes: " << N_classes << endl;
    parfile << "                     Total events: " << N_events << endl;
    parfile << "               Optimization loops: " << N_epochs << endl;
    parfile << "             Optimize SNN params.: ";
    if (update9) {
        parfile << "True" << endl; 
    } else {
        parfile << "False" << endl;
    }
    parfile << "                  Optimize delays: ";
    if (updateDelays) {
        parfile << "True" << endl; 
    } else {
        parfile << "False" << endl;
    }
    parfile << "             Optimize connections: ";
    if (updateConnections) {
        parfile << "True" << endl; 
    } else {
        parfile << "False" << endl;
    }
    parfile << "                  Max mod. factor: " << MaxFactor << endl;
    parfile << "                  Only " << N_TrackingLayers << "-hit tracks: ";
    if (!anyHits) {
        parfile << "True" << endl;
    } else {
        parfile << "False" << endl;
    }
    parfile.close();
    return;
}


// Initialize neuron potentials
// ----------------------------
void Init_neurons () {
    for (int in=0; in<N_neurons; in++) {
        // Set first event in history of this neuron
        History_time[in].push_back(0.);
        History_type[in].push_back(0); 
        History_ID[in].push_back(0);
        if (in<N_neuronsL[0]) {
            Neuron_layer[in] = 0;
        } else {
            Neuron_layer[in] = 1;
        }
    }
    return;
}

// Initialize synapse weights
// --------------------------
void Init_weights () {
    for (int in=0; in<N_neurons; in++) {
        for (int is=0; is<N_streams; is++) {
            check_LTD[in][is] = true; // flags used to see if we need to create a LTD signal after a neuron discharge
            Weight[in][is] = myRNG->Uniform();
            Weight_initial[in][is] = Weight[in][is];
        }
    }    
    return;
}

// Reset synapse weights
// ---------------------
void Reset_weights () {
    for (int in=0; in<N_neurons; in++) {
        for (int is=0; is<N_streams; is++) {
            Weight[in][is] = Weight_initial[in][is];
        }
    }
    return;
}

// Initialize synaptic delays
// --------------------------
void Init_delays () {
    // Define delays for IE signals
    for (int in=0; in<N_neurons; in++) {
        for (int is=0; is<N_streams; is++) {
            Delay[in][is] = 0.;
            if (learnDelays || updateDelays) Delay[in][is] = MaxDelay/2.;
//            if (is<N_TrackingLayers) { // no IE delay for neuron-originated spikes into L1
//                Delay[in][is] = myRNG->Uniform(MaxDelay); 
//            }
        }
    }
    return;
}

// Initialize connection map
// -------------------------
void Init_connection_map() {
    //Setting L0 input connections
    for (int in = 0; in < N_neuronsL[0]; in++)
    {
        //input connections Tracking layers -> L0
        for (int is = 0; is < N_InputStreams; is++)
        {
            Void_weight[in][is] = false;
            if(myRNG->Uniform()>ConnectedFraction_Input_L0) Void_weight[in][is] = true;
        }
        //input connections L0 -> L0
        for (int is = N_InputStreams; is < N_streams; is++)
        {
            Void_weight[in][is] = false;
        }
    }

    //Setting L1 input connections
    for (int in = N_neuronsL[0]; in < N_neurons; in++)
    {
        //input connections Tracking layers -> L1
        for (int is = 0; is < N_InputStreams; is++)
        {
            Void_weight[in][is] = false;
            if(myRNG->Uniform()>ConnectedFraction_Input_L1) Void_weight[in][is] = true;
        }
        //input connections L0 -> L1
        for (int is = N_InputStreams; is < N_streams; is++)
        {
            Void_weight[in][is] = false;
            if(myRNG->Uniform()>ConnectedFraction_L0_L1) Void_weight[in][is] = true;
        }
    }
    /*
        for (int in=0; in<N_neurons; in++) {
        for (int is=0; is<N_streams; is++) {
            if (in<N_neuronsL[0] && is>=N_InputStreams) {
                Void_weight[in][is] = true; 
            } else {
                Void_weight[in][is] = false; 
                if (myRNG->Uniform()>ConnectedFraction) Void_weight[in][is] = true;
            }
        }
    }
    return;
    */

}
//DOUBT: Here, it seems that the connections between the input layer and L1 are not removed by default.
//From the construction of the weight initialization function, it appears that such connection exists.

// Reset hits
// ----------
void Reset_hits () {
    hit_pos.clear();
    return;
}

// Simulate background hits
// ------------------------
void Simulate_bgrhits (double bgr_rate) {
    for (int itl=0; itl<N_TrackingLayers; itl++) {
        for (int is=0; is<N_strips; is++) {
            if (myRNG->Uniform()<bgr_rate){ 
                float r = strip_r[itl];
                float phi = strip_phi[is];
                hit_pos.emplace_back(r, 0, phi, BGR);
            }
        }   
    }
    return;
}

// Simulate hits from particle interactions
// This is a very crude model which does not account for any of the following (possible future implementation extensions):
// - charge deposition model in silicon sensors
// - multi-strip clustering
// - angle of incidence of track in silicon layer
// - electronic noise
// -----------------------------------------------------------------------------------------------------------------------
int Simulate_sighits () {
    // We first compute xc, yc, p, beta, charge
    // p is momentum, in GeV/c
    // beta is azimuthal angle from x axis

    int Nhits    = 0;
    StartLayer   = 0;
    
    //TODO implement V particles
    if (KindOfSignal>0) StartLayer = KindOfSignal*myRNG->Uniform(3.-epsilon); // Only simulate V-particles generated in first four layers
  
    float p     = getP(); 
    float beta  = getBeta(); 
    short int ch    = getCharge();
    float R     = p*abs(ch)/(0.3*Bfield)*100.;  //(in cm)
    
    float min_angle = max_angle; 

    //We compute the particle's ideal orbit: a circumference passing through the origin.
    //To determine the intersection between layers and the circumference, we employed the numerical bisection method.
    //For the first layer, we conducted the bisection within an angular interval centered around phi_orthogonal, defined as follows:

    float phi_ort;
    if(ch<0) 
        phi_ort = beta + (M_PI/2.);
    else
        phi_ort = beta + 3*(M_PI/2.);
    if (phi_ort > 2.*M_PI) phi_ort-= 2*M_PI;
    
    
    for (int itl=StartLayer; itl<N_TrackingLayers; itl++) { 
        float r_hit  = strip_r[itl];
        for (int ich=0; ich<N_part; ich++) { 
            //TODO  modify to include V-particles
            /*
            if (ich>0 && (KindOfSignal==1 || myRNG->Uniform()>0.5)) { // V-particles are +- pairs (of same momentum and beta, for now)
                ch = -ch;
                R  = -R;                
            }
            */
            
            //computing center position
            float yc = R*sin(beta);
            float xc = R*cos(beta);

            float phi_hit = bisection(phi_ort - bisection_window,phi_ort + bisection_window, bisection_precision, r_hit, xc, yc, R);

            
            if(phi_hit<0) phi_hit += 2*M_PI;
            else if (phi_hit > 2 * M_PI) phi_hit -= 2*M_PI;

            if(phi_hit < min_angle)
                min_angle = phi_hit;
            
            //updating the searching window for the next layer
            phi_ort = phi_hit;

            hit_pos.emplace_back(r_hit,0, phi_hit, SIG);
            Nhits++;
        }
    }
    First_angle = min_angle;
    return Nhits;
}

// Encode hits in spike stream / we do it in a single stream with code for layer
// -----------------------------------------------------------------------------
int GetITL(float r_hit){
    if (r_hit < 0) r_hit = 0;
    if (r_hit > max_R) r_hit = max_R-epsilon;

    return (int)(r_hit / r_strip); 
    
    /*
        for (int i = 0; i < N_TrackingLayers; i++){
        
        if((r_hit > strip_r[i] - confidence_r_left[i]) && (r_hit < strip_r[i] + confidence_r_right[i])){ 
            cout << i;
            return i;
        }
        
    }
    return N_TrackingLayers+1;

    */
}

int GetIEta(float eta){
    double tmp = (eta + eta_range/2.);
    if (tmp<0) tmp = 0;
    else if(tmp > eta_range) tmp = eta_range -epsilon;
    
    return (int) (tmp/ eta_strip);
}

int GetStreamID(int r, int eta){
    return r+eta*N_TrackingLayers; 
    //return r*N_etaLayers+z;
}
void GetR_Z(int combind){
    int r = combind%N_TrackingLayers;
    int z = (combind-r)/N_EtaLayers;
}

void Encode (double t_in) { 
    //sort by phi angle
    sort(hit_pos.begin(), hit_pos.end(), [](const Hit& h1, const Hit& h2) {
        return h1.phi < h2.phi;
    });

    for (auto &&row : hit_pos)
    {
        double time = t_in + row.phi/omega;
        //i'll need a row.eta
        //uncomment when implemented:
        int itl = GetStreamID(GetITL(row.r), GetIEta(row.eta));

        PreSpike_Time.push_back(time);
        PreSpike_Stream.push_back(itl);
        PreSpike_Signal.push_back(row.id-1); // 0,1,2 -> -1,0,1 respectively NoHit, Backgroung, Signal
    }

    //rescan from [0, delta]
    for (auto &&row : hit_pos)
    {
        if(row.phi > delta) break;
        double time = t_in + (row.phi + M_PI *2.)/omega;
        //uncomment when implemented:
        int itl = GetStreamID(GetITL(row.r), GetIEta(row.eta));
        
        PreSpike_Time.push_back(time);
        PreSpike_Stream.push_back(itl);
        PreSpike_Signal.push_back(row.id-1); // 0,1,2 -> -1,0,1 respectively NoHit, Backgroung, Signal
    }

}
// Model Excitatory Post-Synaptic Potential
// We take this as parametrized in T. Masquelier et al., "Competitive STDP-Based Spike Pattern Learning", DOI: 10.1162/neco.2008.06-08-804
// ---------------------------------------------------------------------------------------------------------------------------------------
double EPS_potential (double delta_t) {
    double psp = 0.;
    if (delta_t>=0. && delta_t<MaxDeltaT) psp = K * (exp(-delta_t/tau_m)-exp(-delta_t/tau_s));
    return psp;
}

// Model membrane potential after spike
// Also modeled as in paper cited above, like IPSP and other signals below
// -----------------------------------------------------------------------
double Spike_potential (double delta_t, int ilayer) {
    double sp = 0.;
    if (delta_t>=0. && delta_t<MaxDeltaT) sp = Threshold[ilayer] * (K1*exp(-delta_t/tau_m)-K2*(exp(-delta_t/tau_m)-exp(-delta_t/tau_s)));
    return sp;
}

// Model Inhibitory-Excitatory signal (IE) as combination of two EPSP-like shapes, a negative one followed by a positive one
// This is a crude model, loosely inspired by shapes in Fig.2 of F. Sandin and M. Nilsson, "Synaptic Delays for Insect-Inspired
// Feature Detection in Dynamic Neuromorphic Processors", doi.org/10.3389/fnins.2020.00150
// ----------------------------------------------------------------------------------------------------------------------------
double IE_potential (double delta_t, int in, int is) {
    double sp = 0.;
    if (delta_t>=0. && delta_t<Delay[in][is]) {
        sp = -IE_Pot_const * EPS_potential(delta_t);
    } else if (delta_t>=Delay[in][is] && delta_t<MaxDeltaT+Delay[in][is]) { 
        delta_t = delta_t - Delay[in][is];
        sp = IE_Pot_const * EPS_potential (delta_t); // So for zero delay, this is an EPSP
    }
    return sp;
}

// Model Inhibitory Post-Synaptic Potential (IPSP)
// -----------------------------------------------
double Inhibitory_potential (double delta_t, int ilayer) {
    // In order to dilate the inhibition to larger times (in the attempt at obtaining higher selectivity), 
    // we kludge it by multiplying delta_t and maxdeltat by a factor
    delta_t = delta_t * IPSP_dt_dilation;
    double ip = 0.;
    double thisalpha = alpha;
    if (ilayer>0) thisalpha = L1inhibitfactor*alpha; // Different inhibition in L1
    if (delta_t>=0. && delta_t<MaxDeltaT) ip = -thisalpha * Threshold[ilayer] * EPS_potential(delta_t); 
    return ip;
}

// Model Spike-Timing-Dependent Plasticity (LTP - Long-Term Potentiation) - modify weights 
// based on how close a previous synapse fired before the spike)
// ---------------------------------------------------------------------------------------
void LTP (int in, int is, int this_spike, double fire_time) {
    if (Void_weight[in][is]) return;
    // Use nearest-spike approximation: search for closest pre-spike
    bool no_prespikes = true;
    int isp = this_spike-1;
    do {
        if (PreSpike_Stream[isp]==is) {
            double delta_t = PreSpike_Time[isp]-fire_time;
            Weight[in][is] += a_plus*exp(delta_t/tau_plus);
            if (Weight[in][is]>1.) Weight[in][is] = 1.;
            no_prespikes = false;
        }
        isp--;
    } while (isp>=0 && PreSpike_Time[isp]>fire_time-7.*tau_plus && no_prespikes);

    // Also modify delays (experimental)
    if (learnDelays) {
        double dtbig = MaxDelay / NevPerEpoch;
        if (Delay[in][is]>=dtbig) Delay[in][is] -= dtbig;
        double dt = dtbig / N_neuronsL[Neuron_layer[in]];
        int inmin = 0;
        int inmax = N_neuronsL[0];
        if (Neuron_layer[in]==1) {
            inmin = N_neuronsL[0];
            inmax = N_neurons;
        }
        for (int in2=inmin; in2<inmax; in2++) { 
            if (in2!=in && Delay[in2][is]<=MaxDelay-dt) Delay[in2][is] += dt;
        }
    }
    return;
}

// Model Spike-Timing-Dependent Plasticity (LTD - Long-Term Depression)
// --------------------------------------------------------------------
void LTD (int in, int is, double spike_time) {
    // Use nearest-spike approximation: search for closest neuron spike to this input spike 
    if (Fire_time[in].size()==0) return;
    if (Void_weight[in][is]) return;
    double delta_t = spike_time-Fire_time[in].back();
    check_LTD[in][is] = false;
    if (delta_t>=0 && delta_t<7.*tau_minus) {
        Weight[in][is] -= a_minus*exp(-delta_t/tau_minus);
        if (Weight[in][is]<0.) Weight[in][is] = 0.;
    }
    // Also modify delay (experimental)
    if (learnDelays) {
        double dtbig = MaxDelay / NevPerEpoch;
        if (Delay[in][is]<MaxDelay-dtbig) Delay[in][is] += dtbig;
        double dt = dtbig/(N_neurons-1);
        for (int in2=0; in2<N_neurons; in2++) {
            if (in2!=in  && Neuron_layer[in]==Neuron_layer[in2] && Delay[in2][is]>dt) Delay[in2][is] -= dt;
        }
    }
    return;
}

// Compute collective effect of excitatory, post-spike, and inhibitory potentials on a neuron
// ------------------------------------------------------------------------------------------
double Neuron_firetime (int in, double t) {
    int ilayer = Neuron_layer[in];
    double P0 = 0.;
    double t0 = History_time[in][0];
    double delta_t = t - t0;
    if (t0>0. && delta_t>=0. && delta_t<MaxDeltaT) {
        int ilayer = Neuron_layer[in];
        P0 = Spike_potential(delta_t, ilayer); // the first event in the history sequence is a spike
    }
    double P = P0;

    // Now we extrapolate the effect of all past spikes and inhibitions to time t, to compute the potential when EPSP arrives
    int len = History_time[in].size();
    if (len>1) {
        for (int ih=1; ih<len-1; ih++) {
            delta_t = t - History_time[in][ih];
            if (History_type[in][ih]==1) { // EPSP
                if (delta_t>MaxDeltaT) { // Get rid of irrelevant events
                    History_time[in].erase(History_time[in].begin()+ih,History_time[in].begin()+ih+1);
                    History_type[in].erase(History_type[in].begin()+ih,History_type[in].begin()+ih+1);
                    History_ID[in].erase(History_ID[in].begin()+ih,History_ID[in].begin()+ih+1);
                    len = len-1;
                } else {
                    if (!Void_weight[in][History_ID[in][ih]]) // for type 1 or 3 signals, ID is the stream
                        P += Weight[in][History_ID[in][ih]]*EPS_potential(delta_t);
                }
            } else if (History_type[in][ih]==2) { // IPSP
                if (delta_t>MaxDeltaT) { // get rid of irrelevant events
                    History_time[in].erase(History_time[in].begin()+ih,History_time[in].begin()+ih+1);
                    History_type[in].erase(History_type[in].begin()+ih,History_type[in].begin()+ih+1);
                    History_ID[in].erase(History_ID[in].begin()+ih,History_ID[in].begin()+ih+1);
                    len = len-1;
                } else {
                    int ilayer = Neuron_layer[in];
                    P += Inhibitory_potential(delta_t, ilayer);
                }
            } else if (History_type[in][ih]==3) { // IE
                if (delta_t>MaxDeltaT) { // get rid of irrelevant events
                    History_time[in].erase(History_time[in].begin()+ih,History_time[in].begin()+ih+1);
                    History_type[in].erase(History_type[in].begin()+ih,History_type[in].begin()+ih+1);
                    History_ID[in].erase(History_ID[in].begin()+ih,History_ID[in].begin()+ih+1);
                    len = len-1;
                } else {
                    if (!Void_weight[in][History_ID[in][ih]]) // for type 1 or 3 signals, ID is the stream
                        P += IE_potential(delta_t, in, History_ID[in][ih]);
                }
            }
        }
    }
    if (P>Threshold[ilayer]) { // Neuron will fire as spike contribution will bring it above threshold
        // compute fire time by looping more finely from t to t+tmax (tmax is peak time of EPSP)
        double this_t = t;
        do {
            P = P0;
            for (int ih=1; ih<len; ih++) {
                double delta_t = this_t-History_time[in][ih];
                if (History_type[in][ih]==1) { // EPSP
                    if (!Void_weight[in][History_ID[in][ih]])
                        P += Weight[in][History_ID[in][ih]]*EPS_potential(delta_t);
                } else if (History_type[in][ih]==2) { // IPSP
                    int ilayer = Neuron_layer[in];
                    P += Inhibitory_potential(delta_t, ilayer);
                } else if (History_type[in][ih]==3) { // IE
                    if (!Void_weight[in][History_ID[in][ih]])
                        P += IE_potential(delta_t, in, History_ID[in][ih]);
                }
            }
            this_t += 1/(10000.*omega);
        } while (P<Threshold[ilayer] && this_t<=t+tmax);
        if (P>=Threshold[ilayer]) return this_t;
    }   
    return largenumber;
}

// For debugging purposes, compute neuron potential at given time, using same model above
// (warning, may be outdated)
// --------------------------------------------------------------------------------------
double Neuron_P (int in, double t) {
    double P0 = 0.;
    double t0 = History_time[in][0];
    double delta_t = t-t0;
    if (t0>0. && delta_t>=0. && delta_t<MaxDeltaT) {
        int ilayer = Neuron_layer[in];
        P0 = Spike_potential(delta_t, ilayer); // the first event in the history sequence is a spike
    }
    double P = P0;

    // Now we extrapolate the effect of all past spikes and inhibitions to time t, to compute the potential when EPSP arrives
    int len = History_time[in].size();
    if (len>1) {
        for (int ih=1; ih<len; ih++) {
            delta_t = t-History_time[in][ih];
            if (History_type[in][ih]==1) { // EPSP
                if (!Void_weight[in][History_ID[in][ih]])
                    P += Weight[in][History_ID[in][ih]]*EPS_potential(delta_t);
            } else if (History_type[in][ih]==2) { // IPSP
                int ilayer = Neuron_layer[in];
                P += Inhibitory_potential(delta_t, ilayer);
            } else if (History_type[in][ih]==3) { // IE
                if (!Void_weight[in][History_ID[in][ih]])
                    P += IE_potential (delta_t, in, History_ID[in][ih]);
            }
        }
    }
    return P;
}

// Learning rate scheduler - this returns an oscillating, dampened function as a function of the epoch
// ---------------------------------------------------------------------------------------------------
double LR_Scheduler (double LR0, int epoch, int Nepochs) {
    double par[3] = {-0.01,0.2,0.2};
    double x = 100.*epoch/Nepochs;
    return LR0*exp(par[0]*x)*(par[1]+(1.-par[1])*pow(cos(par[2]*x),2));
}

// Calculate selectivity of set of neurons
// --------------------------------------- 
double Compute_Selectivity(int level, int mode) {
    double S = 0.;
    int inmin, inmax;
    if (level==0) {
        inmin = 0;
        inmax = N_neuronsL[0];
    } else if (level==1) {
        inmin = N_neuronsL[0];
        inmax = N_neurons;
    }
    if (mode==0) { // Use additive rule
        for (int in=inmin; in<inmax; in++) {
        // select max efficiency class
            double maxeff = 0.;
            double sumeff = 0.;
            for (int ic=0; ic<N_classes; ic++) {
                double e = Eff[ic+N_classes*in];
                if (e>maxeff) maxeff = e; 
                sumeff = sumeff + e;
            }
            if (sumeff>0.) S += maxeff*N_classes/sumeff;
        }
        if (inmax-inmin>0) S = S/(inmax-inmin);
    } else if (mode==1) { // aim for collective effect (each neuron participates)
        S = 1.;
        for (int in=inmin; in<inmax; in++) {
        // select max efficiency class
            double maxeff = 0.;
            double sumeff = 0.;
            for (int ic=0; ic<N_classes; ic++) {
                double e = Eff[ic+N_classes*in];
                if (e>maxeff) maxeff = e; 
                sumeff = sumeff + e;
            }
            if (sumeff>0.) S *= maxeff/sumeff;
        }
    } else if (mode==2) { // compute mutual information
        // I(N_neurons,N_classes) = Sum_i^N_n Sum_j^N_c Eff(i,j) log_2 [Eff(i,j)/Eff_i Eff_j)]
        // where  is the average efficiency of neuron i over classes, and Eff_j is the
        // average efficiency on class j over neurons
        S = 0.;
        double Effn[MaxNeurons];
        double Effc[MaxClasses];
        double sumeff = 0.;
        for (int in=inmin; in<inmax; in++) {
            sumeff = 0.;
            for (int ic=0; ic<N_classes; ic++) {
                sumeff += Eff[ic+N_classes*in];
            }
            if (N_classes>0) Effn[in] = sumeff/N_classes;
        }
        for (int ic=0; ic<N_classes; ic++) {
            sumeff = 0.;
            for (int in=inmin; in<inmax; in++) {
                sumeff += Eff[ic+N_classes*in];
            }
            if (inmax>inmin) Effc[ic] = sumeff/(inmax-inmin);
        }
        for (int in=inmin; in<inmax; in++) {
            for (int ic=0; ic<N_classes; ic++) {
                if (Effc[ic]*Effn[in]>0.) S += Eff[ic+N_classes*in] * 
                                               (log2(Eff[ic+N_classes*in]+epsilon) - log2(Effc[ic]*Effn[in]));
            }
        }
    }
    if (inmax>inmin) S /= N_classes*(inmax-inmin);
    return S;
}

// Compute Q-value
// ---------------
double Compute_Q (double eff, double acc, double sel) {
    // efficiency saturates at eff_target
    if (eff>eff_target) eff = eff_target + pow(eff-eff_target,1.5);
    // acceptance saturates at acc_target
    if (acc<acc_target) acc = acc_target - pow(acc_target-acc,1.5);
    double Q0 = eff/sqrt(acc);
    double w  = 5.*(exp(Q0/4.)-1.)/(exp(1.)-1.);
    return Q0 + w * sel;
}

//TODO: optimize
void createRootFile(const char* nomeFile, int nev, int classes=6, double occ = 0.000390625) {
    N_classes = classes;
    Occupancy = occ;
    Define_tracker ();
    
    cout << "Create file" << endl;
    TFile* file = TFile::Open(nomeFile, "RECREATE");
    if (!file || file->IsZombie()) {
        cerr << "Errore nell'apertura o creazione del file ROOT." << endl;
        return;
    }

    cout << "Create variables" << endl;

    //Creating variables
    float r, phi;
    short int id;
    long int id_event;

    cout << "Create tree" << endl;
    TTree *tree = new TTree("tree", "tree");

    //Arrange the variables as pointers to branches of the tree 
    tree->Branch("r", &r);
    tree->Branch("phi", &phi);
    tree->Branch("id", &id);
    tree->Branch("id_event", &id_event);
    tree->Branch("pclass", &pclass);


    cout << "Generate events" << endl;

    for(id_event = 0; id_event<nev; id_event++){
        // Hit generation
        // --------------
        N_part     = 0;
        pclass = -1;
        int Nhits  = 0;            
        // Reset detector hits
        Reset_hits ();
        bool signal = false;
        if (myRNG->Uniform()>0.5) signal = true;
        if (signal) { // This is an event with signal in it (50%)
            N_part = 1; // + myRNG->Uniform(2.-epsilon); 
            KindOfSignal = 0; // turn off v-particles generation // myRNG->Uniform(2.-epsilon); 
            if (KindOfSignal>0) N_part = 2; 
            do {                
                // Reset detector hits if not first pass
                if (Nhits>0) Reset_hits ();
                First_angle = max_angle; // static that indicates leftmost strip hit by particle (or pair). Assigned in Simulate_sighits
                pclass = (int)myRNG->Uniform((double)N_classes-epsilon);                
                // Simulate hits from particle track
                Nhits = Simulate_sighits ();     
            } while ((!anyHits && Nhits<N_TrackingLayers) || Nhits<N_TrackingLayers); 
        }

        // Simulate background hits
        double bgr_rate = Occupancy;
        if (!signal) {
            Simulate_bgrhits (bgr_rate);
        } else {
            Simulate_bgrhits (bgr_rate+1./N_strips); // Make the occupancy identical in events with and without tracks
        }

        //Filling the TTree
        for (long int i = 0; i < hit_pos.size(); i++) {
            r = hit_pos[i].r;
            phi = hit_pos[i].phi;
            id = hit_pos[i].id;
            tree->Fill();
        }
    }
    tree->Write(); 

}

//To read simulated events

//Legge i dati GENERATI da noi
void ReadFromRoot(TTree* tree, long int id_event_value){
    // Open the ROOT file
    
    Reset_hits();
    int countPart = 0;

    // Attach the branches to variables
    float r;
    float phi;
    short int id;
    long int id_event;

    tree->SetBranchAddress("r", &r);
    tree->SetBranchAddress("phi", &phi);
    tree->SetBranchAddress("id", &id);
    tree->SetBranchAddress("id_event", &id_event);
    tree->SetBranchAddress("pclass", &pclass);

 
    // Loop over entries and find rows with the specified id_event value
    for (long int i = last_row_event; i < tree->GetEntries(); ++i) {
        tree->GetEntry(i);
        if (id_event!=id_event_value) {
            last_row_event = i;
            break;
        }
        if(static_cast<int>(id)==SIG){
            countPart++;
        }
        hit_pos.emplace_back(r,0, phi, id);
    }
    if(pclass==-1)
        pclass++;
    N_part =countPart/N_TrackingLayers;
}

//To read pure Monte Carlo events
//TODO work in progress
void ReadFromMia(TTree* IT, TTree* OT, long int id_event_value){
    Reset_hits();
  
    float r;
    float phi;
    float id;
    float id_event;
    float eta;

    IT->SetBranchAddress("cluster_R", &r);
    IT->SetBranchAddress("cluster_phi", &phi);
    IT->SetBranchAddress("cluster_type", &id);
    IT->SetBranchAddress("eventID", &id_event);
    IT->SetBranchAddress("cluster_eta", &eta);

    pclass = 0;

    // Loop over entries and find rows with the specified id_event value
    for (long int i = last_row_event; i < IT->GetEntries(); ++i) {
        IT->GetEntry(i);
        if (static_cast<long int>(id_event)!=id_event_value) {
            last_row_event = i;
            break;
        }
        if(static_cast<int>(id)==1){
            id = SIG;
        }
        else id=BGR;
        
        hit_pos.emplace_back(r,eta, phi, static_cast<int>(id));
    }
    
    //OUT Tracker
    
    OT->SetBranchAddress("cluster_R", &r);
    OT->SetBranchAddress("cluster_phi", &phi);
    OT->SetBranchAddress("cluster_type", &id);
    OT->SetBranchAddress("eventID", &id_event);
    OT->SetBranchAddress("cluster_eta", &eta);


    for (long int i = last_row_event_OT; i < OT->GetEntries(); ++i) {
        OT->GetEntry(i);
        if (static_cast<long int>(id_event)!=id_event_value) {
            last_row_event_OT = i;
            break;
        }
        if(static_cast<int>(id)==1){
            id = SIG;
        }
        else id=BGR;
        
        hit_pos.emplace_back(r,eta, phi, static_cast<int>(id));
    }

    N_part = 2;
}

//To read our preprocessed file
void ReadFromProcessed(TTree* IT, TTree* OT, long int id_event_value){
    Reset_hits();
    pclass = 0;
    N_part = 0;

    float z;
    float r, phi;
    float id_event;
    float type;
    float cluster_pclass;

    IT->SetBranchAddress("cluster_eta", &z);
    IT->SetBranchAddress("cluster_R", &r);
    IT->SetBranchAddress("cluster_phi", &phi);
    IT->SetBranchAddress("eventID", &id_event);
    IT->SetBranchAddress("cluster_type", &type);
    IT->SetBranchAddress("pclass", &cluster_pclass);

    // Loop over entries and find rows with the specified id_event value
    for (long int i = last_row_event; i < IT->GetEntries(); ++i) {
        IT->GetEntry(i);
        
        if (static_cast<long int>(id_event)!=id_event_value) {
            last_row_event = i;
            break;
        }
        phi+=M_PI;
        if(static_cast<int>(type)==1){
            type = SIG;
            pclass = (int) cluster_pclass;
            N_part=1;
            phi+=2.*M_PI*((int)(ievent/(NROOT))) * 1./ ((int)(N_events/NROOT) + 1 );
        }
        else type = BGR;

        if(phi>=2.*M_PI) phi-= 2.*M_PI;
        
        hit_pos.emplace_back(r,z, phi, static_cast<int>(type));
    }
    
    //OUT Tracker
    
    OT->SetBranchAddress("cluster_eta", &z);
    OT->SetBranchAddress("cluster_R", &r);
    OT->SetBranchAddress("cluster_phi", &phi);
    OT->SetBranchAddress("eventID", &id_event);
    OT->SetBranchAddress("cluster_type", &type);

    for (long int i = last_row_event_OT; i < OT->GetEntries(); ++i) {
        OT->GetEntry(i);   
        if (static_cast<long int>(id_event)!=id_event_value) {
            last_row_event_OT = i;
            break;
        }
        phi += M_PI;
        if(static_cast<int>(type)==1){
            phi+=2.*M_PI*((int)(ievent/(NROOT))) * 1./ ((int)(N_events/NROOT) + 1 );
            type = SIG;
        }
        else type=BGR;
        
        //cout << 1.*((int)(ievent/(NROOT))) * 1./ ((int)(N_events/NROOT)) << endl;
        if(phi>=2.*M_PI) phi-= 2.*M_PI;
        hit_pos.emplace_back(r,z, phi, static_cast<int>(type));
    }
}

// -----------------------------------------------------------------------------------------------------------------------------------------
// Main routine
// ------------
void SNN_Tracking (int N_ev, int N_ep, int NL0, int NL1, char* rootInput = nullptr, bool batch=false, double CFI0=1, double CFI1 = 1, double CF01 = 1,
                   double Thresh0=15, double Thresh1=10, double _MaxFactor = 0.2, double a=0.25, double l1if=1., double k=1., double k1=2., double k2=4., 
                   double IEPC=2.5, double ipspdf=1.0, double _MaxDelay = 0.1e-9, double _tau_m = 1e-09, double _tau_s = 0.25e-09, 
                   double _tau_plus = 1.68e-09, double _tau_minus = 3.37e-09, double _a_plus = 0.03125, double _a_minus = 0.03125, 
                   int N_cl=6, 
                   int TrainingCode=5,bool ReadPars=false, double Occ=0.000390625, long int _NROOT = 100000) {

    // Pass parameters:
    // ----------------__
    // N_ev:      total number of simulated events
    // N_ep:      number of weight-learning cycles divido N_ev in N_ep gruppi e faccio girare il learning 
    // NL0, NL1:  number of neurons performing track pattern recognition, organized in 2 layers
    // N_cl:      number of different signal classes (different particle momenta)
    // rootInput: name of the root file if you want to load the data. If not provided it will simulate the events
    // CF:        fraction of connected neurons between L0 and L1
    // Occ:       probability of random hit firing 
    // ipspdf:    IPSP time dilation factor (to increase effect of inhibition)
    // Trainingcode: binary code to turn on updates of 9 pars, delays, voids
    // ReadPars:  whether to read parameters in from file
    // Thresh0:   threshold for firing at Layer 0
    // Thresh1:   threshold for firing at Layer 1
    // alpha:     parameter of signal model
    // l1if:      L1 inhibit factor
    // k, k1, k2: parameters of signal model
    // IEPC:      inhibition-excitation potential constant
    // IPSPdf:    IPSP dt dilation factor

    // The routine works as follows:
    // -----------------------------
    // Define tracker geometry ok
    // Initialize neuron potentials and synapse weights ok
    // Simulate background hits 
    // Simulate hits from particle interactions 
    // Encode hits in spike streams
    // Loop on optimization cycles (N_epochs)
    //   Loop on time steps in event-based fashion and modify neuron and synapse potentials
    //     Determine when an output spike will occur and jump there
    //     If spike: 
    //       - record spike 
    //       - update synapses 
    //       - update neuron 
    //       - inhibit other neurons
    //       - check latency with respect to signal patterns, record it     
    // Dump statistics for run
    // ----------------------------------------------------------------------------------
    //
    // Every tracking layer provides an encoded stream of spikes into a dense layer of L0 neurons;
    // Every L0 neuron is densely connected to a L1 neuron, which also receives all tracking streams.
    // The topology is sketched below for 4 tracking layers, 3 L0 neurons, 2 L1 neurons
    //
    //  ------------|---- ....... -|---
    //                         O  ---|-
    //  ----|---|-------- ....... -----      O-|->
    //                     X   O  -----   X
    //  ------|----|----- ....... --|--      O--->
    //                         O  -|---
    //  --------------|-- ....... -----      
    //
    //  Notes:
    //  - The optimization strategy is naive and needs to be improved
    //  - In addition, rather than re-learning weights every time parameters are modified, some better
    //    way of handling this should be implemented
    //  - IE spike model is to be improved
    // -----------------------------------------------------------------------------------------------

    // Pass parameters can't update static values, so we need to reassign the latter
    if(batch) gROOT->SetBatch(kTRUE);
    
    N_events          = N_ev;
    N_epochs          = N_ep;
    N_neuronsL[0]     = NL0;
    N_neuronsL[1]     = NL1;
    ConnectedFraction_Input_L0  = CFI0;
    ConnectedFraction_Input_L1  = CFI1;
    ConnectedFraction_L0_L1     = CF01;
    Threshold[0]      = Thresh0;
    Threshold[1]      = Thresh1;
    alpha             = a;
    L1inhibitfactor   = l1if;
    K                 = k;
    K1                = k1; 
    K2                = k2;
    IE_Pot_const      = IEPC;
    IPSP_dt_dilation  = ipspdf;
    MaxDelay          = _MaxDelay;
    tau_m             = _tau_m;
    tau_s             = _tau_s;
    tau_plus          = _tau_plus;
    tau_minus         = _tau_minus;
    a_plus            = _a_plus;
    a_minus           = _a_minus;
    MaxFactor        = _MaxFactor;
    NROOT = _NROOT;

    N_neurons = N_neuronsL[0] + N_neuronsL[1];
    N_streams = N_InputStreams + N_neuronsL[0];
    NevPerEpoch = N_events/N_epochs;
    N_classes         = N_cl;
    Occupancy         = Occ;
    
    // Assign meta-learning booleans
    update9           = false;
    updateDelays      = false;
    updateConnections = false;

    if (TrainingCode/4>0) {
        update9 = true;
        TrainingCode-=4;
    }
    if (TrainingCode/2>0) {
        updateDelays = true;
        TrainingCode-=4;
    }
    if (TrainingCode>0) {
        updateConnections = true;
    }

    // Initial checks
    // --------------
    if (N_ev>MaxEvents) {
        cout << "  Sorry, max # of events is 10,000,000. Terminating." << endl;
        return;
    }
    if (N_ev/N_ep<10000) {
        cout << "  Too few events per epoch. Set to " << 10000*N_ep << endl;
        //N_ev = 10000*N_ep;
    }
    if (N_ep<1) {
        cout << "  Invalid N_epochs = " << N_ep << ". Set to 1." << endl;
        N_ep = 1;
    }
    if (NL0+NL1>MaxNeurons) {
        cout << "  Sorry, too many neurons. Terminating." << endl;
        return;
    }
    if (N_classes>MaxClasses) {
        cout << "  Sorry, too many classes (max is " << MaxClasses << "). Terminating." << endl;
        return;
    }

    // Welcome screen
    // --------------
    cout << endl << endl;
    cout << "                                 ------------------------------------" << endl;
    cout << endl;
    cout << "                                    S   N   N      T r a c k i n g"    << endl;
    cout << endl;
    cout << "                                 ------------------------------------" << endl;
    cout << endl << endl << endl << endl;
    cout << "         ------------------------------------------------------------------------------------    " << endl;
    cout << "         Unsupervised search for tracks in 8-layer strip detector with spiking neural network    " << endl;
    cout << "                                                                             T.Dorigo, 3/2023    " << endl;
    cout << "         ------------------------------------------------------------------------------------    " << endl;
    cout << endl;
    cout << "         Run parameters: " << endl;
    cout << "         -----------------------------------" << endl;
    cout << "                       L0 neurons: " << NL0 << endl;
    cout << "                       L1 neurons: " << NL1 << endl;
    cout << "            Connected L0-L1 frac.: " << CF01 << endl;
    cout << "            Connected IN-L0 frac.: " << CFI0 << endl;
    cout << "            Connected IN-L1 frac.: " << CFI1 << endl;
    cout << "                  Noise per layer: " << Occ*N_strips << endl;
    cout << "                    Track classes: " << N_cl << endl;
    cout << "                     Total events: " << N_ev << endl;
    cout << "               Optimization loops: " << N_ep << endl;
    cout << "             Optimize SNN params.: ";
    if (update9) {
        cout << "True" << endl; 
    } else {
        cout << "False" << endl;
    }
    cout << "                  Optimize delays: ";
    if (updateDelays) {
        cout << "True" << endl; 
    } else {
        cout << "False" << endl;
    }
    cout << "             Optimize connections: ";
    if (updateConnections) {
        cout << "True" << endl; 
    } else {
        cout << "False" << endl;
    }
    cout << "                  Max mod. factor: " << MaxFactor << endl;
    cout << "                Only " << N_TrackingLayers << "-hit tracks: ";
    if (!anyHits) {
        cout << "True" << endl;
    } else {
        cout << "False" << endl;
    }

    // Suppress root warnings
    gROOT->ProcessLine( "gErrorIgnoreLevel = 6001;");
    gROOT->ProcessLine( "gPrintViaErrorHandler = kTRUE;");

    // Histograms definition
    // ---------------------
    TH1F * SelectivityL0 = new TH1F ("SelectivityL0", "", N_epochs, 0.5, 0.5+N_epochs);
    TH1F * SelectivityL1 = new TH1F ("SelectivityL1", "", N_epochs, 0.5, 0.5+N_epochs);
    TH1F * Qvalue        = new TH1F ("Qvalue",        "", N_epochs, 0.5, 0.5+N_epochs);
    TH1F * Qmax          = new TH1F ("Qmax",          "", N_epochs, 0.5, 0.5+N_epochs);
    TH1F * HEff          = new TH1F ("HEff",          "", N_epochs, 0.5, 0.5+N_epochs);
    TH1F * HAcc          = new TH1F ("HAcc",          "", N_epochs, 0.5, 0.5+N_epochs);
    TH1F * HT0           = new TH1F ("HT0",           "", N_epochs, 0.5, 0.5+N_epochs);
    TH1F * HT1           = new TH1F ("HT1",           "", N_epochs, 0.5, 0.5+N_epochs);
    TH1F * HA            = new TH1F ("HA",            "", N_epochs, 0.5, 0.5+N_epochs);
    TH1F * HL1IF         = new TH1F ("HL1IF",         "", N_epochs, 0.5, 0.5+N_epochs);
    TH1F * HK            = new TH1F ("HK",            "", N_epochs, 0.5, 0.5+N_epochs);
    TH1F * HK1           = new TH1F ("HK1",           "", N_epochs, 0.5, 0.5+N_epochs);
    TH1F * HK2           = new TH1F ("HK2",           "", N_epochs, 0.5, 0.5+N_epochs);
    TH1F * HIEPC         = new TH1F ("HIEPC",         "", N_epochs, 0.5, 0.5+N_epochs);
    TH1F * HIPSPdf       = new TH1F ("HIPSPdf",       "", N_epochs, 0.5, 0.5+N_epochs);
    TH2F * EffMap        = new TH2F ("EffMap",        "", N_neurons, -0.5, N_neurons-0.5, N_classes, -0.5, N_classes-0.5);
    SelectivityL1->SetLineColor(kBlack);
    Qmax->SetLineColor(2);
    HEff->SetMaximum(1.1);
    HEff->SetMinimum(0.);
    HAcc->SetLineColor(kRed);
    HAcc->SetMaximum(1.1);
    HAcc->SetMinimum(0.);
    HA->SetLineColor(kBlack);
    HT1->SetLineColor(kRed);
    HK1->SetLineColor(kRed);
    HK2->SetLineColor(kBlue);
    HT0->SetMinimum(0.);
    HK2->SetMinimum(0.);
    Qvalue->SetMinimum(0.);
    Qmax->SetMinimum(0.);
    HA->SetMinimum(0.);
    HIEPC->SetMinimum(0.);
    HIPSPdf->SetMinimum(0.);
    HL1IF->SetMinimum(0.);
    TH1F * HDelays = new TH1F ("HDelays", "", 50, 0., MaxDelay);
    TH2F * HVoidWs = new TH2F ("HVoidWs", "", N_neurons, -0.5, -0.5+N_neurons, N_streams, -0.5, -0.5+N_streams);
    TH2F * Q_12 = new TH2F ("Q_12","", 20, -MaxFactor, MaxFactor, 20, -MaxFactor, MaxFactor);
    TH2F * Q_34 = new TH2F ("Q_34","", 20, -MaxFactor, MaxFactor, 20, -MaxFactor, MaxFactor);
    TH2F * Q_56 = new TH2F ("Q_56","", 20, -MaxFactor, MaxFactor, 20, -MaxFactor, MaxFactor);
    TH2F * Q_78 = new TH2F ("Q_78","", 20, -MaxFactor, MaxFactor, 20, -MaxFactor, MaxFactor);
    TH2F * Q_93 = new TH2F ("Q_93","", 20, -MaxFactor, MaxFactor, 20, -MaxFactor, MaxFactor);
    TH2F * Q_MV = new TH2F ("Q_MV","", 20, -MaxFactor, MaxFactor, 20, -MaxFactor, MaxFactor);
    TH2F * N_12 = new TH2F ("N_12","", 20, -MaxFactor, MaxFactor, 20, -MaxFactor, MaxFactor);
    TH2F * N_34 = new TH2F ("N_34","", 20, -MaxFactor, MaxFactor, 20, -MaxFactor, MaxFactor);
    TH2F * N_56 = new TH2F ("N_56","", 20, -MaxFactor, MaxFactor, 20, -MaxFactor, MaxFactor);
    TH2F * N_78 = new TH2F ("N_78","", 20, -MaxFactor, MaxFactor, 20, -MaxFactor, MaxFactor);
    TH2F * N_93 = new TH2F ("N_93","", 20, -MaxFactor, MaxFactor, 20, -MaxFactor, MaxFactor);
    TH2F * N_MV = new TH2F ("N_MV","", 20, -MaxFactor, MaxFactor, 20, -MaxFactor, MaxFactor);

    int N_bins = 100;
    TH2F * Latency[MaxNeurons*MaxClasses];
    char name[50];
    for (int i=0; i<N_neurons*N_classes; i++) {
        sprintf (name,"Latency%d", i);
        Latency[i] = new TH2F (name, name, N_bins, 0., (double)NevPerEpoch, max_angle+Empty_buffer, 0., (max_angle+Empty_buffer)/omega);
    }
    TH1F * HWeight[MaxNeurons*MaxStreams];
    TH1F * Efficiency[MaxNeurons*MaxClasses];
    TH1F * FakeRate[MaxNeurons];
    TH1F * Eff_totL0[MaxClasses];
    TH1F * Eff_totL1[MaxClasses];
    TH2F * StreamsS[10];
    TH2F * StreamsB[10];
    TH2F * StreamsN[10];
    TH1F * BestEff[MaxNeurons];
    TH1F * BestFR[MaxNeurons];
    TH1F * BestEtot[MaxNeurons];

    for (int i=0; i<N_neurons*N_streams; i++) {
        sprintf (name,"HWeight%d", i);
        HWeight[i] = new TH1F (name, name, N_bins, 0., (double)NevPerEpoch);
    }
    for (int i=0; i<N_neurons*N_classes; i++) {
        sprintf (name,"Efficiency%d", i);
        Efficiency[i] = new TH1F (name, name, N_epochs, 0.5, 0.5+N_epochs);
    }
    for (int in=0; in<N_neurons; in++) {
        sprintf (name,"FakeRate%d", in);
        FakeRate[in]  = new TH1F (name, name, N_epochs, 0.5, 0.5+N_epochs);
    }
    for (int in=0; in<N_neurons; in++) {
        sprintf (name, "BestEff%d", in);
        BestEff[in] = new TH1F (name, name, N_classes, -0.5, -0.5+N_classes);
        sprintf (name, "BestFR%d", in);
        BestFR[in]  = new TH1F (name, name, N_classes, -0.5, -0.5+N_classes);
        sprintf (name, "BestEtot%d", in);
        BestEtot[in]  = new TH1F (name, name, N_classes, -0.5, -0.5+N_classes);
    }
    for (int ic=0; ic<N_classes; ic++) {
        sprintf (name,"Eff_totL0%d", ic);
        Eff_totL0[ic] = new TH1F (name, name, N_epochs, 0.5, 0.5+N_epochs);
        sprintf (name,"Eff_totL1%d", ic);
        Eff_totL1[ic] = new TH1F (name, name, N_epochs, 0.5, 0.5+N_epochs);
    }

    for (int i=0; i<10; i++) {
        sprintf (name, "StreamsS%d", i);
        StreamsS[i] = new TH2F (name, name, (max_angle+Empty_buffer)*500, 0., (max_angle+Empty_buffer)*50./omega, N_InputStreams+N_neurons, 0.5, N_InputStreams+N_neurons+0.5);
        sprintf (name, "StreamsB%d", i);
        StreamsB[i] = new TH2F (name, name, (max_angle+Empty_buffer)*500, 0., (max_angle+Empty_buffer)*50./omega, N_InputStreams+N_neurons, 0.5, N_InputStreams+N_neurons+0.5);
        sprintf (name, "StreamsN%d", i);
        StreamsN[i] = new TH2F (name, name, (max_angle+Empty_buffer)*500, 0., (max_angle+Empty_buffer)*50./omega, N_InputStreams+N_neurons, 0.5, N_InputStreams+N_neurons+0.5);
    }

    // Calculation of constant in excitation spike, to make it max at 1
    // double deltat_max = (tau_m*tau_s)/(tau_m-tau_s)*log(tau_m/tau_s);
    // K = 1./(exp(-deltat_max/tau_m)-exp(-deltat_max/tau_s)); // Now an optimization parameter

    // Calculation of time at maximum of EPSP
    tmax = tau_s*tau_m/(tau_m-tau_s)*(log(tau_m)-log(tau_s));

    // Calculation of max spike height
    Pmax_EPSP = EPS_potential(tmax);

    // Define tracker geometry
    Define_tracker (); //note, this is currently not used - strip x and y are computed on the fly when needed

    // If requested, read in parameters
    int okfile;
    if (ReadPars) {
        okfile = Read_Parameters();
        if (okfile==-1) return;
    }

    // Final part of initial printout
    cout << "         -----------------------------------" << endl;
    cout << endl;
    cout << "         Starting values of parameters:" << endl;
    cout << "         -----------------------------------" << endl;
    cout << "                     L0 threshold: " << Threshold[0] << endl;
    cout << "                     L1 threshold: " << Threshold[1] << endl;
    cout << "                            alpha: " << alpha << endl;
    cout << "                        L1inhibit: " << L1inhibitfactor << endl;
    cout << "                                K: " << K << endl;
    cout << "                               K1: " << K1 << endl;
    cout << "                               K2: " << K2 << endl;
    cout << "                 IE pot. constant: " << IE_Pot_const << endl; 
    cout << "                 IPSP dt dilation: " << IPSP_dt_dilation << endl;
    cout << "         -----------------------------------" << endl;
    cout << endl;

    // Initialize neuron potentials
    Init_neurons ();

    // Initialize synapse activation
    Init_weights ();

    // Initialize delays
    if (!ReadPars) Init_delays ();

    // Initialize connection map
    if (!ReadPars) Init_connection_map();

    // Prime the event loop - we continuously sample detector readout and feed inputs to synapses
    // ------------------------------------------------------------------------------------------
    int N_fires[MaxNeurons];
    double LastP[MaxNeurons];
    for (int in=0; in<N_neurons; in++) { 
        N_fires[in] = 0.;
        LastP[in]   = 0.;
    }
    int fired_sum[MaxClasses][MaxNeurons];
    int random_fire[MaxNeurons];
    for (int in=0; in<N_neurons; in++) {
        random_fire[in] = 0;
        for (int ic=0; ic<N_classes; ic++) {
            fired_sum[ic][in] = 0;
        }
    }
    bool not_fired_bgr  = true;
    int atleastonefired = 0;
    int gen_sum[MaxClasses];
    int fired_anyL0[MaxClasses];
    int fired_anyL1[MaxClasses];
    for (int ic=0; ic<N_classes; ic++) {
        gen_sum[ic]              = 0;
        fired_anyL0[ic]          = 0;
        fired_anyL1[ic]          = 0;
    }
    bool doneL0[MaxClasses];
    bool doneL1[MaxClasses];
    bool Seen[MaxClasses][MaxNeurons];
    double selectivityL0 = 0.;
    double selectivityL1 = 0.;
    double averefftotL1  = 0.;
    double averacctotL1  = 0.;

    // Storing parameters subjected to random search
    double oldThresholdL0     = Threshold[0];
    double oldThresholdL1     = Threshold[1];
    double oldalpha           = alpha;
    double oldL1inhibitfactor = L1inhibitfactor;
    double oldK               = K;
    double oldK1              = K1;
    double oldK2              = K2;
    double oldIE_Pot_const    = IE_Pot_const;
    double oldIPSPdf          = IPSP_dt_dilation;
    double oldDelay[MaxNeurons][MaxStreams];
    for (int in=0; in<N_neurons; in++) {
        for (int is=0; is<N_streams; is++) {
            oldDelay[in][is]  = Delay[in][is];
            bestDelay[in][is] = Delay[in][is];
        }
    }
    bool oldVoid_weight[MaxNeurons][MaxStreams];
    for (int in=0; in<N_neurons; in++) {
        for (int is=0; is<N_streams; is++) {
            oldVoid_weight[in][is]  = Void_weight[in][is];
            bestVoid_weight[in][is] = Void_weight[in][is];
        }
    }
    double Q     = 0.;
    double Q_old = 0.;
    Q_best       = 0.;
    SelL1_best   = 0.;
    Eff_best     = 0.;
    Acc_best     = 0.;
    T0_best      = 0.;
    T1_best      = 0.;
    A_best       = 0.;
    L1if_best    = 0.;
    K_best       = 0.;
    K1_best      = 0.;
    K2_best      = 0.;
    IEPC_best    = 0.;
    IPSPdf_best  = 0.;

    double Optvar[9];
    double max_dx[9];
    double aver_dQ[9];
    Optvar[0] = 0.; // Threshold[0];
    Optvar[1] = 0.; // Threshold[1];
    Optvar[2] = 0.; // alpha;
    Optvar[3] = 0.; // L1inhibitfactor;
    Optvar[4] = 0.; // K;
    Optvar[5] = 0.; // K1;
    Optvar[6] = 0.; // K2;
    Optvar[7] = 0.; // IE_Pot_const;
    Optvar[8] = 0.; // IPSP_dt_dilation;
    for (int i=0; i<9; i++) {
        aver_dQ[i] = 0.;
        max_dx[i]  = MaxFactor; // max factor of change in parameter values during optimization
    }
    double OptvarD[MaxNeurons*MaxStreams];
    double max_dxD[MaxNeurons*MaxStreams];
    double aver_dQD[MaxNeurons*MaxStreams];
    for (int in=0; in<N_neurons; in++) {
        for (int is=0; is<N_streams; is++) {
            int id = in*N_streams+is;
            aver_dQD[id] = 0.;
            max_dxD[id]  = 0.01; // max delay change factor during optimization
            OptvarD[id]  = Delay[in][is];
        }
    }
    int MaxdQHist = 50;   // we do not want the full history, because we are moving in the par space; only last 10 points.
    vector <double> dQHist;
    vector <double> OptvarHist[9];
    vector <double> OptvarDHist[MaxNeurons*MaxStreams];
    int ibad   = 0;
    double LR  = MaxFactor;

    // Big loop on events
    // ------------------
    bool doprogress = true;
    int block = N_events/N_epochs/50;
    if (block<1) doprogress = false;
    if (doprogress) cout << "         " << progress[0];
    int currchar      = 1;
    ievent   = 0;
    int iev_thisepoch = 0;
    int iepoch        = 0;
    int ind_qbest     = 0;
    
    //Open the root data file if provided
    TFile* file = TFile::Open(rootInput, "READ");
    if (!file || file->IsZombie()) {
        cerr << "Error: Cannot open file " << rootInput << endl;
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
    TTree* OT = dynamic_cast<TTree*>(dirOT->Get("tree"));

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

    IT->SetMaxVirtualSize(250000000);
    IT->LoadBaskets();

    OT->SetMaxVirtualSize(250000000);
    OT->LoadBaskets();

    /*
    TFile *file = TFile::Open(rootInput, "READ");
    TTree *tree;
    if (!file || file->IsZombie()) {
        cerr << "Error opening file!" << endl;
        cout << "Falling back to events simulation" << endl;
        rootInput = nullptr;
    }
    else{
         // Get the TTree
        tree = (TTree*)(file->Get("tree;1"));
        if (!tree) {
            cerr << "Error getting TTree!" << endl;
            cout << "Falling back to events simulation" << endl;
            rootInput = nullptr;
            file->Close();
        }
        else{
            //set maximum amount of memory used for caching TBranches to 250MB
            tree->SetMaxVirtualSize(250000000);
            //load all the data blocks into memory
            tree->LoadBaskets();
        }
    }
    */
    
   
    do {
        iev_thisepoch++;
        if (doprogress) {
            if (ievent%block==0) {
	            cout << progress[currchar];
            	currchar++;
            }
        }

        //load data from the root file if provided
        if(rootInput){
            if(ievent%NROOT==0){
                last_row_event=0;
                last_row_event_OT = 0;
            }
            ReadFromProcessed(IT, OT, ievent%NROOT);
        }
        else{
        // Hit generation
        // --------------
        N_part     = 0;
        pclass = 0;
        int Nhits  = 0;            
        // Reset detector hits
        Reset_hits ();
        bool signal = false;
        if (myRNG->Uniform()>0.5) signal = true;
        if (signal) { // This is an event with signal in it (50%)
            N_part = 1; // + myRNG->Uniform(2.-epsilon); 
            KindOfSignal = 0; // turn off v-particles generation // myRNG->Uniform(2.-epsilon); 
            if (KindOfSignal>0) N_part = 2; 
            do {                
                // Reset detector hits if not first pass
                if (Nhits>0) Reset_hits ();
                First_angle = max_angle; // static that indicates leftmost strip hit by particle (or pair). Assigned in Simulate_sighits
                pclass = (int)myRNG->Uniform((double)N_classes-epsilon);                
                // Simulate hits from particle track
                Nhits = Simulate_sighits ();     
            } while ((!anyHits && Nhits<N_TrackingLayers) || Nhits<N_TrackingLayers); 
        }

        // Simulate background hits
        
        double bgr_rate = Occupancy;
        if (!signal) {
            Simulate_bgrhits (bgr_rate);
        } else {
            Simulate_bgrhits (bgr_rate+1./N_strips); // Make the occupancy identical in events with and without tracks
        }
        
    }

        // See if we find with track with positive latency by at least one neuron
        for (int in=0; in<N_neurons; in++) {
            Seen[pclass][in] = false;
        }
        doneL0[pclass] = false;
        doneL1[pclass] = false;
        not_fired_bgr  = true;

        // Encode hits in spike streams
        // Here we encode the position of hits through the timing of a spike,
        // In the future we might think at how the analog charge readout in the hits could also be added
        PreSpike_Time.clear();
        PreSpike_Stream.clear();
        PreSpike_Signal.clear();
        double t_in = ievent*(max_angle+Empty_buffer)/omega; // every event adds 0.256+maxdelay seconds (for 256 strips), allowing for L1 neurons to get L0 signal
        Encode (t_in); 

        // Keep track of latency calc for each neuron in this event
        bool not_filled[MaxNeurons];
        for (int in=0; in<N_neurons; in++) {
            not_filled[in] = true;
        }

        // Loop on spikes and modify neuron and synapse potentials
        // -------------------------------------------------------
        for (int ispike=0; ispike<PreSpike_Time.size(); ispike++) { 
            // By looping to size(), we can insert along the way and still make it to the end
            double t = PreSpike_Time[ispike];

            // Save information on hit-based streams for last 500 events to histograms
            if (ievent>=N_events-500.) {
                int is = (ievent-N_events+500)/50;
                //time = tin + thit - tin(primoevento)
                double time = PreSpike_Time[ispike]-(max_angle+Empty_buffer)/omega*(ievent/50)*50;
                if (PreSpike_Signal[ispike]==1) {
                    StreamsS[is]->Fill(time, PreSpike_Stream[ispike]+1);
                } else if (PreSpike_Signal[ispike]==0) {
                    StreamsB[is]->Fill(time, PreSpike_Stream[ispike]+1);
                } else if (PreSpike_Signal[ispike]==2) {
                    StreamsN[is]->Fill(time, PreSpike_Stream[ispike]+1);
                }
            }

            // Modify neuron potentials based on synapse weights
            // -------------------------------------------------
            double min_fire_time = largenumber-1.; // if no fire, neuron_firetime returns largenumber
            int in_first = -1;
            for (int in=0; in<N_neurons; in++) {
                // We implement a scheme where input streams produce an IE signal into L0, an EPS into L1, and L0 neurons EPS into L1
                // Add to neuron history, masking out L1 spikes for L0 neurons
                int is = PreSpike_Stream[ispike];
                if (is<N_InputStreams || Neuron_layer[in]>0) { // otherwise stream "is" does not lead to neuron "in"
                    History_time[in].push_back(t);
                    if (PreSpike_Signal[ispike]==2) { // L0 neuron-induced spike 
                        History_type[in].push_back(1); 
                    } else { // IE spike in input to neurons 
                        History_type[in].push_back(3);
                    }
                    History_ID[in].push_back(is);

                    // Model STDP: LTD. See if this spike depresses a neuron that fired earlier
                    if (check_LTD[in][is]) LTD(in,is,t);

                    // Compute future fire times of neurons and their order
                    double fire_time = Neuron_firetime(in,t); 
                    if (fire_time<min_fire_time) {
                        in_first = in;
                        min_fire_time = fire_time;
                    }
                }
            }
            if (in_first==-1) continue; // nothing happens, move on

            // Ok, neuron in_first is going to fire next.
            // Peek at next event in list, to see if it comes before in_first fires
            // --------------------------------------------------------------------
            if (ispike<PreSpike_Time.size()-1) {
                if (PreSpike_Time[ispike+1]>=min_fire_time) { // otherwise we go to next spike in list
                    // handle firing of neuron in_first
                    double latency = 0.;
                    N_fires[in_first]++;
                    Fire_time[in_first].push_back(min_fire_time);

                    // Reset history of this neuron
                    History_time[in_first].clear();
                    History_type[in_first].clear();
                    History_ID[in_first].clear();
                    History_time[in_first].push_back(min_fire_time);
                    History_type[in_first].push_back(0);
                    History_ID[in_first].push_back(0); // ID is not used for type 0 history events

                    // IPSP for all others at relevant layer
                    for (int in2=0; in2<N_neurons; in2++) {
                        if (in2!=in_first) {
                            if (Neuron_layer[in2]==Neuron_layer[in_first]) { // inhibitions within layer or across
                                History_time[in2].push_back(min_fire_time);
                                History_type[in2].push_back(2);
                                History_ID[in2].push_back(in_first); 
                            }
                        }
                    }

                    // Learn weights with spike-time-dependent plasticity: long-term synaptic potentiation
                    for (int is=0; is<N_streams; is++) { 
                        LTP(in_first, is, ispike, min_fire_time);
                        // Reset LTD check flags
                        check_LTD[in_first][is] = true;
                    }

                    // Create EPS signal in L0 neuron-originated streams 
                    if (Neuron_layer[in_first]==0) { // this is a Layer-0 neuron
                        PreSpike_Time.insert(PreSpike_Time.begin()+ispike+1, min_fire_time);
                        PreSpike_Stream.insert(PreSpike_Stream.begin()+ispike+1, N_InputStreams+in_first);
                        PreSpike_Signal.insert(PreSpike_Signal.begin()+ispike+1, 2);
                    }

                    // Fill spikes train histogram
                    if (ievent>=N_events-500.) {
                        int is = (ievent-N_events+500)/50;
                        double time = min_fire_time-(max_angle+Empty_buffer)/omega*(ievent/50)*50;
                        if (Neuron_layer[in_first]==1) StreamsN[is]->Fill(time, N_InputStreams+in_first+1);
                    }

                    // Fill latency histogram
                    if (N_part>0) {
                        //quanto tempo ci ha messo il primo neurone a sparare rispestto all'arrivo temporale della prima hit?
                        latency = min_fire_time-t_in-First_angle/omega;
                        if (latency>=0. && not_filled[in_first]) {
                            if (iepoch==N_epochs-1) Latency[in_first*N_classes+pclass]->Fill(0.5+iev_thisepoch,latency); 
                            Seen[pclass][in_first] = true;
                            not_filled[in_first] = false;
                        } 
                    } else {
                        if (not_filled[in_first]) {
                            random_fire[in_first]++;
                            not_filled[in_first] = false;
                        }
                        if (in_first>=N_neuronsL[0] && iev_thisepoch>NevPerEpoch*0.9) { // for Q-value calculations
                            if (not_fired_bgr) {
                                atleastonefired++;
                                not_fired_bgr = false;
                            }
                        }
                    }
                }
            } // end if in_first fires
        } // end ispike loop, ready to start over

        // Fill info for efficiency calculations
        if (N_part>0) {
            if (iev_thisepoch>NevPerEpoch*0.9) {
                gen_sum[pclass]++;
                for (int in=0; in<N_neurons; in++) {
                    if (Seen[pclass][in]) {
                        fired_sum[pclass][in]++;
                        if (in<N_neuronsL[0]) {
                            if (!doneL0[pclass]) {
                                doneL0[pclass] = true;
                                fired_anyL0[pclass]++;
                            } 
                        } else {
                            if (!doneL1[pclass]) {
                                doneL1[pclass] = true;
                                fired_anyL1[pclass]++;
                            }
                        }
                    }
                }
            }
        }

        // Write histograms of weights
        if (iepoch==N_epochs-1) {
            for (int in=0; in<N_neurons; in++) {
                for (int is=0; is<N_streams; is++) {
                    int bin = (int)(1000.*(double)iev_thisepoch/NevPerEpoch);
                    HWeight[in*N_streams+is]->SetBinContent(bin, Weight[in][is]);
                }
            }
        }

        //prespike_time.push_back(time) -> tempo della Hit o dello spike di L0
        //prespike_Stream -> id del tracking layer della hit o id del neurone di L0 che ha sparato
        //prespike_Signal -> spike type: 0 if BKG 1 if Track, 2 if NeuronFire
        // Fill efficiency histograms every NevPerEpoch events, compute Q value and Selectivity, modify parameters
        // ---------------------------------------------------------------------------------------------------
        if (iev_thisepoch==NevPerEpoch) { // we did NevPerEpoch events

            // Reset counter that inhibits efficiency and Q calculations until we reach steady state with weights
            iev_thisepoch = 0;
            iepoch++;
            // End of progress bar
            if (doprogress) cout << progress[51] << endl;

            for (int in=0; in<N_neurons; in++) {
                for (int ic=0; ic<N_classes; ic++) {
                    int combind = ic+N_classes*in;
                    Eff[combind] = fired_sum[ic][in];                
                    if (gen_sum[ic]>0) Eff[combind] /= gen_sum[ic];
                    Efficiency[combind]->SetBinContent(iepoch, Eff[combind]); 
                }
                double fakerate = random_fire[in]*2./NevPerEpoch; // there are NevPerEpoch/2 events with no tracks, where we compute random_fire per neuron
                FakeRate[in]->SetBinContent(iepoch, fakerate);
            }
            double Efftot[MaxClasses];
            for (int ic=0; ic<N_classes; ic++) {
                double etl0 = fired_anyL0[ic];
                if (gen_sum[ic]>0) etl0 /= gen_sum[ic];
                Eff_totL0[ic]->SetBinContent(iepoch, etl0);
                double etl1 = fired_anyL1[ic]; // L1 efficiency is what counts.
                if (gen_sum[ic]>0) etl1 /= gen_sum[ic];
                Eff_totL1[ic]->SetBinContent(iepoch, etl1);
                Efftot[ic] = etl1;
            }

            selectivityL0 = Compute_Selectivity(0,2);
            SelectivityL0->Fill(iepoch,selectivityL0);
            selectivityL1 = Compute_Selectivity(1,2);
            SelectivityL1->Fill(iepoch,selectivityL1);

            // Q value is average efficiency divided by sqrt (aver eff plus aver acceptance)
            // -----------------------------------------------------------------------------
            averacctotL1 = atleastonefired*(2./NevPerEpoch*10.); // total acceptance, computed with 0.1*NevPerEpoch/2 events with no tracks
            averefftotL1 = 0.;
            for (int ic=0; ic<N_classes; ic++) {
                averefftotL1 += Efftot[ic];
            }
            averefftotL1 /= N_classes;
            double den = 0.05 + averacctotL1; // deem 10% fake rate ok-ish
            Q = Compute_Q(averefftotL1, averacctotL1, selectivityL1); 
            
            // Fix maximum excursion of parameters with a schedule
            LR = LR_Scheduler(MaxFactor, iepoch, N_epochs);
            for (int i=0; i<9; i++) {
                max_dx[i] = LR;
            }
            for (int id=0; id<N_neurons*N_streams; id++) {
                max_dxD[id] = 0.1*LR;
            }

            // Re-initialize neurons
            Init_neurons();
            // Reset hits
            Reset_hits();
            // Reset weights to initial conditions before new investigation
            Reset_weights();
            // Init delays
            if (!updateDelays && !ReadPars && !learnDelays) Init_delays(); // This unlike void connections, because we can opt to learn these at each cycle too

            
            cout << "         Ev. # " << ievent+1 << " - LR = " << LR << "; Selectivity L0 = " << selectivityL0 << " L1 = " << selectivityL1 
                 << "; Eff = " << averefftotL1 << " Acc = " << averacctotL1 << "; Firings: ";
            
            for (int in=0; in<N_neurons; in++) {
                cout << N_fires[in] << " ";
            } 
            cout << endl;
            
            

            // Keep a history of recent delta Q values, so that we know how to sample next
            dQHist.push_back(Q);
            for (int i=0; i<9; i++) {
                OptvarHist[i].push_back(Optvar[i]);
            }
            for (int id=0; id<N_neurons*N_streams; id++) {
                OptvarDHist[id].push_back(OptvarD[id]);
            }
        
            // Keep only the last MaxdQHist values of dQ history
            if (dQHist.size()>MaxdQHist) {
                dQHist.erase(dQHist.begin());
                for (int i=0; i<9; i++) {
                    OptvarHist[i].erase(OptvarHist[i].begin());
                }
                for (int id=0; id<N_neurons*N_streams; id++) {
                    OptvarDHist[id].erase(OptvarDHist[id].begin());
                }
            }
            // Find running weighted average of dQ-values in par space
            double dQ_max = 0.;
            for (int i=0; i<9; i++) {
                double sum_dQdx = 0.;
                double sum_dx   = 0.;
                for (int j=1; j<dQHist.size(); j++) {
                    double dx = OptvarHist[i][j]-OptvarHist[i][j-1];
                    sum_dQdx += (dQHist[j]-dQHist[j-1])*dx;
                    sum_dx   += dx;
                }
                aver_dQ[i] = 0.;
                if (sum_dx!=0.) aver_dQ[i] = sum_dQdx/sum_dx;
                if (fabs(aver_dQ[i])>dQ_max) dQ_max = fabs(aver_dQ[i]);
            }
            // Same, for delays
            double dQ_maxD = 0.;
            for (int id=0; id<N_neurons*N_streams; id++) {
                double sum_dQdxD = 0.;
                double sum_dxD   = 0.;
                for (int j=1; j<dQHist.size(); j++) {
                    double dxD = OptvarDHist[id][j]-OptvarDHist[id][j-1];
                    sum_dQdxD += (dQHist[j]-dQHist[j-1])*dxD;
                    sum_dxD   += dxD;
                }
                aver_dQD[id] = 0.;
                if (sum_dxD!=0.) aver_dQD[id] = sum_dQdxD/sum_dxD;
                if (fabs(aver_dQD[id])>dQ_maxD) dQ_maxD = fabs(aver_dQD[id]);
            }

            // Fill debugging graphs of Q as a function of parameters
            int ibin, jbin;
            double n, cont, newcont;
            if (iepoch>1 || N_epochs==1) { // only do it from end of second epoch onwards, as we are filling delta values
                ibin = 1+(int)(20.*(Threshold[0]/oldThresholdL0-1.+MaxFactor)/(2.*MaxFactor));
                jbin = 1+(int)(20.*(Threshold[1]/oldThresholdL1-1.+MaxFactor)/(2.*MaxFactor));
                if (ibin>0 && ibin<21 && jbin>0 && jbin<21) {
                    cont    = Q_12->GetBinContent(ibin,jbin);
                    n       = N_12->GetBinContent(ibin,jbin);
                    newcont = (Q-Q_old)/(n+1) + cont*n/(n+1);
                    Q_12->SetBinContent(ibin,jbin,newcont);
                    N_12->SetBinContent(ibin,jbin,n+1);
                }
                ibin = 1+(int)(20.*(alpha/oldalpha-1.+MaxFactor)/(2*MaxFactor));
                jbin = 1+(int)(20.*(L1inhibitfactor/oldL1inhibitfactor-1.+MaxFactor)/(2.*MaxFactor));
                if (ibin>0 && ibin<21 && jbin>0 && jbin<21) {
                    cont    = Q_34->GetBinContent(ibin,jbin);
                    n       = N_34->GetBinContent(ibin,jbin);
                    newcont = (Q-Q_old)/(n+1) + cont*n/(n+1);
                    Q_34->SetBinContent(ibin,jbin,newcont);
                    N_34->SetBinContent(ibin,jbin,n+1);
                }
                ibin = 1+(int)(20.*(K/oldK-1.+MaxFactor)/(2.*MaxFactor));
                jbin = 1+(int)(20.*(K1/oldK1-1.+MaxFactor)/(2.*MaxFactor));
                if (ibin>0 && ibin<21 && jbin>0 && jbin<21) {
                    cont    = Q_56->GetBinContent(ibin,jbin);
                    n       = N_56->GetBinContent(ibin,jbin);
                    newcont = (Q-Q_old)/(n+1) + cont*n/(n+1);
                    Q_56->SetBinContent(ibin,jbin,newcont);
                    N_56->SetBinContent(ibin,jbin,n+1);
                }
                ibin = 1+(int)(20.*(K2/oldK2-1.+MaxFactor)/(2.*MaxFactor));
                jbin = 1+(int)(20.*(IE_Pot_const/oldIE_Pot_const-1.+MaxFactor)/(2.*MaxFactor));
                if (ibin>0 && ibin<21 && jbin>0 && jbin<21) {
                    cont    = Q_78->GetBinContent(ibin,jbin);
                    n       = N_78->GetBinContent(ibin,jbin);
                    newcont = (Q-Q_old)/(n+1) + cont*n/(n+1);
                    Q_78->SetBinContent(ibin,jbin,newcont);
                    N_78->SetBinContent(ibin,jbin,n+1);
                }
                jbin = 1+(int)(20.*(alpha/oldalpha-1.+MaxFactor)/(2.*MaxFactor));
                ibin = 1+(int)(20.*(IPSP_dt_dilation/oldIPSPdf-1.+MaxFactor)/(2.*MaxFactor));
                if (ibin>0 && ibin<21 && jbin>0 && jbin<21) {
                    cont    = Q_93->GetBinContent(ibin,jbin);
                    n       = N_93->GetBinContent(ibin,jbin);
                    newcont = (Q-Q_old)/(n+1) + cont*n/(n+1);
                    Q_93->SetBinContent(ibin,jbin,newcont);
                    N_93->SetBinContent(ibin,jbin,n+1);
                }
                // Now graph of mean vs sqm of Delay distribution
                double meanDelay    = 0.;
                double sqmDelay     = 0.;
                double meanOldDelay = 0.;
                double sqmOldDelay  = 0.;
                for (int in=0; in<N_neurons; in++) {
                    for (int is=0; is<N_streams; is++) {
                        meanDelay    += Delay[in][is];
                        sqmDelay     += pow(Delay[in][is],2);
                        meanOldDelay += oldDelay[in][is];
                        sqmOldDelay  += pow(oldDelay[in][is],2);
                    }
                }
                meanDelay /= N_neurons*N_streams;
                sqmDelay = sqmDelay/(N_neurons*N_streams) - meanDelay*meanDelay;
                sqmDelay = sqrt(sqmDelay);
                meanOldDelay /= N_neurons*N_streams;
                sqmOldDelay = sqmOldDelay/(N_neurons*N_streams) - meanOldDelay*meanOldDelay;
                sqmOldDelay = sqrt(sqmOldDelay);
                ibin = 1+(int)(20.*(meanDelay/(meanOldDelay+epsilon)-1.+MaxFactor)/(2.*MaxFactor));
                jbin = 1+(int)(20.*(sqmDelay/(sqmOldDelay+epsilon)-1.+MaxFactor)/(2.*MaxFactor));
                if (ibin>0 && ibin<21 && jbin>0 && jbin<21) {
                    cont    = Q_MV->GetBinContent(ibin,jbin);
                    n       = N_MV->GetBinContent(ibin,jbin);
                    newcont = (Q-Q_old)/(n+1) + cont*n/(n+1);
                    Q_MV->SetBinContent(ibin,jbin,newcont);
                    N_MV->SetBinContent(ibin,jbin,n+1);
                }
                TCanvas * QQ = new TCanvas ("QQ","", 800,600);
                QQ->Divide(3,2);
                QQ->cd(1);
                Q_12->Draw("COL4");
                QQ->cd(2);
                Q_34->Draw("COL4");
                QQ->cd(3);
                Q_56->Draw("COL4");
                QQ->cd(4);
                Q_78->Draw("COL4");
                QQ->cd(5);
                Q_93->Draw("COL4");
                QQ->cd(6);
                Q_MV->Draw("COL4");
                QQ->Update();
            }

            // Fill histograms with delays
            HDelays->Reset();
            HVoidWs->Reset();
            for (int in=0; in<N_neurons; in++) {
                for (int is=0; is<N_streams; is++) {
                    HDelays->Fill(Delay[in][is]);
                    if (Void_weight[in][is]) {
                        HVoidWs->SetBinContent(in+1,is+1,0.);
                    } else {
                        HVoidWs->SetBinContent(in+1,is+1,1.);
                    }
                }
            }

            // Is this Q factor not larger than before?
            cout << "         Q = " << Q << " Old = " << Q_old << " Best = " << Q_best << " ib = " << ibad;

            // Update histograms with current parameter values and optimization metrics
            Qvalue->SetBinContent(iepoch, Q);
            if (Q>Q_best) {
                ind_qbest   = iepoch;
                Q_best      = Q;
                SelL1_best  = selectivityL1;
                Eff_best    = averefftotL1;
                Acc_best    = averacctotL1;
                T0_best     = Threshold[0];
                T1_best     = Threshold[1];
                A_best      = alpha;
                L1if_best   = L1inhibitfactor;
                K_best      = K;
                K1_best     = K1;
                K2_best     = K2;
                IEPC_best   = IE_Pot_const;
                IPSPdf_best = IPSP_dt_dilation;
                for (int in=0; in<N_neurons; in++) {
                    for (int is=0; is<N_streams; is++) {
                        bestDelay[in][is] = Delay[in][is];
                        bestVoid_weight[in][is] = Void_weight[in][is];
                    }
                }
            }
            Qmax->SetBinContent(iepoch, Q_best);
            HEff->SetBinContent(iepoch, averefftotL1);
            HAcc->SetBinContent(iepoch, averacctotL1);
            HT0->SetBinContent(iepoch, Threshold[0]);
            HT1->SetBinContent(iepoch, Threshold[1]);
            HA->SetBinContent(iepoch, alpha);
            HL1IF->SetBinContent(iepoch, L1inhibitfactor);
            HK->SetBinContent(iepoch, K);
            HK1->SetBinContent(iepoch, K1);
            HK2->SetBinContent(iepoch, K2);
            HIEPC->SetBinContent(iepoch, IE_Pot_const);
            HIPSPdf->SetBinContent(iepoch,IPSP_dt_dilation);
            TCanvas * CU = new TCanvas ("CU","",1600,700);
            CU->Divide(5,2);
            CU->cd(1);
            Qvalue->Draw();
            Qmax->Draw("SAME");
            CU->cd(2);
            HEff->Draw();
            HAcc->Draw("SAME");
            double h;
            double hmax = -largenumber;
            for (int ibin=1; ibin<=N_epochs; ibin++) {
                h = HT0->GetBinContent(ibin);
                if (hmax<h) hmax = h;
                h = HT1->GetBinContent(ibin);
                if (hmax<h) hmax = h;
            }
            HT0->SetMaximum(hmax+0.1*fabs(hmax));
            HT1->SetMaximum(hmax+0.1*fabs(hmax));
            CU->cd(3);
            HT0->Draw();
            HT1->Draw("SAME");
            HA->Draw("SAME");
            CU->cd(4);
            HL1IF->Draw();
            hmax = -largenumber;
            for (int ibin=1; ibin<=N_epochs; ibin++) {
                h = HK->GetBinContent(ibin);
                if (hmax<h) hmax = h;
                h = HK1->GetBinContent(ibin);
                if (hmax<h) hmax = h;
                h = HK2->GetBinContent(ibin);
                if (hmax<h) hmax = h;
            }
            HK->SetMaximum(hmax+0.1*fabs(hmax));
            HK1->SetMaximum(hmax+0.1*fabs(hmax));
            HK2->SetMaximum(hmax+0.1*fabs(hmax));            
            CU->cd(5);
            HK2->Draw();
            HK->Draw("SAME");
            HK1->Draw("SAME");
            CU->cd(6);
            HIEPC->Draw();
            CU->cd(7);
            HIPSPdf->Draw();
            CU->cd(8);
            HDelays->Draw();
            CU->cd(9); 
            HVoidWs->Draw("COL4");
            CU->cd(10);
            hmax = -largenumber;
            for (int ibin=1; ibin<=N_epochs; ibin++) {
                h = SelectivityL0->GetBinContent(ibin);
                if (hmax<h) hmax = h;
                h = SelectivityL1->GetBinContent(ibin);
                if (hmax<h) hmax = h;
            }
            SelectivityL0->SetMaximum(hmax+0.1*fabs(hmax));
            SelectivityL1->SetMaximum(hmax+0.1*fabs(hmax));
            SelectivityL1->Draw();
            SelectivityL0->Draw("SAME");
            CU->Update();

            if (iepoch==1) { // The first time we modify at random the parameters

                // Store previous values
                if (update9) {
                    oldThresholdL0     = Threshold[0];
                    oldThresholdL1     = Threshold[1];
                    oldalpha           = alpha;
                    oldL1inhibitfactor = L1inhibitfactor;
                    oldK               = K;
                    oldK1              = K1;
                    oldK2              = K2;
                    oldIE_Pot_const    = IE_Pot_const;
                    oldIPSPdf          = IPSP_dt_dilation;
                    for (int i=0; i<9; i++) {
                        Optvar[i] = myRNG->Uniform(-max_dx[i],max_dx[i]);
                    }
                    Threshold[0]    *= 1.+Optvar[0]; 
                    Threshold[1]    *= 1.+Optvar[1];
                    alpha           *= 1.+Optvar[2];
                    L1inhibitfactor *= 1.+Optvar[3];
                    K               *= 1.+Optvar[4];
                    K1              *= 1.+Optvar[5];
                    K2              *= 1.+Optvar[6];
                    IE_Pot_const    *= 1.+Optvar[7];
                    IPSP_dt_dilation*= 1.+Optvar[8];
                }
                if (updateDelays) {
                    for (int in=0; in<N_neurons; in++) {
                        for (int is=0; is<N_streams; is++) {
                            int id = in*N_streams+is;
                            oldDelay[in][is] = Delay[in][is];
                            OptvarD[id]   = myRNG->Uniform(-max_dxD[id],max_dxD[id]);
                            Delay[in][is] += OptvarD[id];
                            if (Delay[in][is]>MaxDelay) Delay[in][is] = MaxDelay;
                            if (Delay[in][is]<0.) Delay[in][is] = 0.;
                        }
                    }
                }
                if (updateConnections) {
                    for (int in=0; in<N_neurons; in++) {
                        for (int is=0; is<N_streams; is++) {
                            if (in>=N_neuronsL[0] || is<N_InputStreams) {
                                oldVoid_weight[in][is] = Void_weight[in][is];
                                if (!Void_weight[in][is]) {
                                    if (myRNG->Uniform()<ProbWSwitchDown) { 
                                        Void_weight[in][is] = !Void_weight[in][is];
                                    }
                                } else {
                                    if (myRNG->Uniform()<ProbWSwitchUp) { 
                                        Void_weight[in][is] = !Void_weight[in][is];
                                    }
                                }
                            }
                        }
                    }
                }
                Q_old = Q;
                ibad  = 0;

            } else { // We are in second or larger epoch, can look at history of improvement for directions

                if (Q<=Q_old) { // We did a step in the wrong direction, need to take it back and modify at random, sampling wisely...
                               // ibad counts how many trials we make from last improved Q. After 3 unsuccessful trials we allow for a step away
                    if (ibad<2) {
                        // Store previous values
                        if (update9) {
                            Threshold[0]     = oldThresholdL0;
                            Threshold[1]     = oldThresholdL1;
                            alpha            = oldalpha;
                            L1inhibitfactor  = oldL1inhibitfactor;
                            K                = oldK;
                            K1               = oldK1;
                            K2               = oldK2;
                            IE_Pot_const     = oldIE_Pot_const;
                            IPSP_dt_dilation = oldIPSPdf;
                        }
                        if (updateDelays) {
                            for (int in=0; in<N_neurons; in++) {
                                for (int is=0; is<N_streams; is++) {
                                    Delay[in][is] = oldDelay[in][is];
                                }
                            }
                        }
                        if (updateConnections) {
                            for (int in=0; in<N_neurons; in++) {
                                for (int is=0; is<N_streams; is++) {
                                    Void_weight[in][is] = oldVoid_weight[in][is];
                                }
                            }                            
                        }
                        // The calculations below do the following:
                        // - find a number between exp(-2.) and exp(+2.) to fix the slope of the dx cumulant distribution
                        // - find the multiplier of each parameter as a function of how much dx increases Q
                        //   The random number is converted by the function pow(r,lambda) to a number between -max_dx and max_dx
                        //   which distributes uniformly for no slope of dq vs dx, and peaky at the extrema for larger correlation 
                        for (int i=0; i<9; i++) {
                            double lambda = 1.;
                            if (dQ_max!=0.) lambda = exp(2.*aver_dQ[i]/dQ_max);
                            double r = myRNG->Uniform();
                            Optvar[i] = -max_dx[i] + 2.*max_dx[i]*pow(r,lambda);
                        }
                        if (update9) { 
                            Threshold[0]     *= 1.+Optvar[0]; 
                            Threshold[1]     *= 1.+Optvar[1];
                            alpha            *= 1.+Optvar[2];
                            L1inhibitfactor  *= 1.+Optvar[3];
                            K                *= 1.+Optvar[4];
                            K1               *= 1.+Optvar[5];
                            K2               *= 1.+Optvar[6];
                            IE_Pot_const     *= 1.+Optvar[7];
                            IPSP_dt_dilation *= 1.+Optvar[8];
                        }
                        if (updateDelays) {
                            double lambda;
                            for (int in=0; in<N_neurons; in++) {
                                for (int is=0; is<N_streams; is++) {
                                    int id = in*N_streams+is;
                                    lambda = 1.;
                                    if (dQ_maxD>0.) lambda = exp(2.*aver_dQD[id]/dQ_maxD);
                                    double r = myRNG->Uniform();
                                    OptvarD[id] = -max_dxD[id] + 2.*max_dxD[id]*pow(r,lambda);
                                    Delay[in][is] += OptvarD[id];
                                    if (Delay[in][is]>MaxDelay) Delay[in][is] = MaxDelay;
                                    if (Delay[in][is]<0.) Delay[in][is] = 0.;
                                }
                            }
                        }
                        if (updateConnections) {
                            for (int in=0; in<N_neurons; in++) {
                                for (int is=0; is<N_streams; is++) {
                                    if (in>=N_neuronsL[0] || is<N_InputStreams) {
                                        if (!Void_weight[in][is]) {
                                            if (myRNG->Uniform()<ProbWSwitchDown) { 
                                                Void_weight[in][is] = !Void_weight[in][is];
                                            }
                                        } else {
                                            if (myRNG->Uniform()<ProbWSwitchUp) { 
                                                Void_weight[in][is] = !Void_weight[in][is];
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        ibad++;
                    } else { // ibad=2, need to reset

                        if (update9) {
                            Threshold[0]     = T0_best;
                            Threshold[1]     = T1_best;
                            alpha            = A_best;
                            L1inhibitfactor  = L1if_best;
                            K                = K_best;
                            K1               = K1_best;
                            K2               = K2_best;
                            IE_Pot_const     = IEPC_best;
                            IPSP_dt_dilation = IPSPdf_best;
                            oldThresholdL0     = Threshold[0];
                            oldThresholdL1     = Threshold[1];
                            oldalpha           = alpha;
                            oldL1inhibitfactor = L1inhibitfactor;
                            oldK               = K;
                            oldK1              = K1;
                            oldK2              = K2;
                            oldIE_Pot_const    = IE_Pot_const;
                            oldIPSPdf          = IPSP_dt_dilation;
                            for (int i=0; i<9; i++) {
                                Optvar[i] = myRNG->Uniform(-max_dx[i],max_dx[i]);
                            }
                            Threshold[0]    *= 1.+Optvar[0]; 
                            Threshold[1]    *= 1.+Optvar[1];
                            alpha           *= 1.+Optvar[2];
                            L1inhibitfactor *= 1.+Optvar[3];
                            K               *= 1.+Optvar[4];
                            K1              *= 1.+Optvar[5];
                            K2              *= 1.+Optvar[6];
                            IE_Pot_const    *= 1.+Optvar[7];
                            IPSP_dt_dilation*= 1.+Optvar[8];
                        }
                        if (updateDelays) {
                            for (int in=0; in<N_neurons; in++) {
                                for (int is=0; is<N_streams; is++) {
                                    int id = in*N_streams+is;
                                    Delay[in][is]    = bestDelay[in][is];
                                    oldDelay[in][is] = Delay[in][is];
                                    OptvarD[id] = myRNG->Uniform(-max_dxD[id],max_dxD[id]);
                                    Delay[in][is] += OptvarD[id];
                                    if (Delay[in][is]>MaxDelay) Delay[in][is] = MaxDelay;
                                    if (Delay[in][is]<0.) Delay[in][is] = 0.;
                                }
                            }
                        }
                        if (updateConnections) {
                            for (int in=0; in<N_neurons; in++) {
                                for (int is=0; is<N_streams; is++) {
                                    if (in>=N_neuronsL[0] || is<N_InputStreams) {
                                        Void_weight[in][is]    = bestVoid_weight[in][is];
                                        oldVoid_weight[in][is] = Void_weight[in][is];
                                        if (!Void_weight[in][is]) {
                                            if (myRNG->Uniform()<ProbWSwitchDown) { 
                                                Void_weight[in][is] = !Void_weight[in][is];
                                            }
                                        } else {
                                            if (myRNG->Uniform()<ProbWSwitchUp) { 
                                                Void_weight[in][is] = !Void_weight[in][is];
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        Q_old = Q_best;
                        ibad  = 0;
                    }

                } else if (Q>Q_old) { // The direction was ok 

                    // Find parameter multipliers, to continue in the same direction that improved Q (with some added momentum and stochasticity)
                    if (update9) {
                        double RT0 = myRNG->Gaus(1.1,0.1)*Threshold[0]/oldThresholdL0;
                        double RT1 = myRNG->Gaus(1.1,0.1)*Threshold[1]/oldThresholdL1;
                        double Ra  = myRNG->Gaus(1.1,0.1)*alpha/oldalpha;
                        double RI  = myRNG->Gaus(1.1,0.1)*L1inhibitfactor/oldL1inhibitfactor;
                        double RK  = myRNG->Gaus(1.1,0.1)*K/oldK;
                        double RK1 = myRNG->Gaus(1.1,0.1)*K1/oldK1;
                        double RK2 = myRNG->Gaus(1.1,0.1)*K2/oldK2;
                        double RIE = myRNG->Gaus(1.1,0.1)*IE_Pot_const/oldIE_Pot_const; 
                        double RID = myRNG->Gaus(1.1,0.1)*IPSP_dt_dilation/oldIPSPdf;
                        // Store previous values
                        oldThresholdL0     = Threshold[0];
                        oldThresholdL1     = Threshold[1];
                        oldalpha           = alpha;
                        oldL1inhibitfactor = L1inhibitfactor;
                        oldK               = K;
                        oldK1              = K1;
                        oldK2              = K2;
                        oldIE_Pot_const    = IE_Pot_const;
                        oldIPSPdf          = IPSP_dt_dilation;
                        // And update parameters
                        Optvar[0] = RT0-1.;
                        Optvar[1] = RT1-1.;
                        Optvar[2] = Ra-1.;
                        Optvar[3] = RI-1.;
                        Optvar[4] = RK-1.;
                        Optvar[5] = RK1-1.;
                        Optvar[6] = RK2-1.;
                        Optvar[7] = RIE-1.;
                        Optvar[8] = RID-1.;
                        Threshold[0]    *= RT0;
                        Threshold[1]    *= RT1;
                        alpha           *= Ra;
                        L1inhibitfactor *= RI;
                        K               *= RK;
                        K1              *= RK1;
                        K2              *= RK2;
                        IE_Pot_const    *= RIE;
                        IPSP_dt_dilation*= RID;
                    }
                    // Same story, for delays
                    if (updateDelays) {
                        for (int in=0; in<N_neurons; in++) {
                            for (int is=0; is<N_streams; is++) {
                                double R = myRNG->Gaus(1.1,0.1)*(Delay[in][is]-oldDelay[in][is]);
                                oldDelay[in][is] = Delay[in][is];
                                if (R>0.) {
                                    Delay[in][is] += R;
                                    OptvarD[in*N_streams+is] = R;
                                }
                                if (Delay[in][is]>MaxDelay) Delay[in][is] = MaxDelay;
                                if (Delay[in][is]<0.) Delay[in][is] = 0.;
                            }
                        }
                    }
                    if (updateConnections) {
                        // Need to do nothing, no "direction" to go further in
                    }
                    Q_old = Q;
                }                 

            } // if not iepoch==1

            if (ievent<N_events-1) { // Otherwise we graciously exit loop
                if (update9) {
                    cout << " - Try TL0 = " << Threshold[0] 
                        << " TL1 = " << Threshold[1] << " a = " << alpha << " L1inh = " << L1inhibitfactor 
                        << " K = " << K << " K1 = " << K1 << " K2 = " << K2 << " IEPC = " << IE_Pot_const << " IPSPdf = " << IPSP_dt_dilation << endl << endl;
                } else {
                    cout << endl << endl;
                }

                // Reset a few counters
                for (int in=0; in<N_neurons; in++) { 
                    N_fires[in] = 0.;
                }
                for (int in=0; in<N_neurons; in++) {
                    random_fire[in] = 0;
                    for (int ic=0; ic<N_classes; ic++) {
                        fired_sum[ic][in] = 0;
                    }
                }
                for (int ic=0; ic<N_classes; ic++) {
                    gen_sum[ic]              = 0;
                    fired_anyL0[ic]          = 0;
                    fired_anyL1[ic]          = 0;
                }
                atleastonefired = 0;

                // Reset progress bar
                if (doprogress) {
                    cout << "         " << progress[0];
                    currchar = 1;
                }
            }

        } // if ievent+1%NevPerEpoch = 0

        ievent++; // only go to next event if we did a backward pass too

    } while (ievent<N_events);

    //closing the input file
    delete IT;
    delete OT;
    delete dirIT;
    delete dirOT;

    file->Close();
    delete file;

    // Draw histograms
    // ---------------
    cout << "Drawing histos" << endl;
    TCanvas * S = new TCanvas ("S","",1000,1000);
    S->Divide(2,5);
    for (int i=0; i<10; i++) {
        S->cd(i+1);
        StreamsB[i]->SetLineColor(kRed);
        StreamsB[i]->Draw("BOX");

        
        StreamsS[i]->SetLineColor(kBlue);
        StreamsS[i]->Draw("BOXSAME");


        StreamsN[i]->SetLineColor(kGreen);
        StreamsN[i]->Draw("BOXSAME");
    }

    TCanvas * C = new TCanvas ("C", "", 1000,1000);
    C->Divide(N_classes, N_neurons);
    for (int i=0; i<N_neurons*N_classes; i++) {
        C->cd(i+1);
        Latency[i]->Draw("COL4");
    }

    TCanvas * E0 = new TCanvas ("E0","", 800,800);
    E0->Divide(N_classes, N_neuronsL[0]);
    for (int i=0; i<N_neuronsL[0]*N_classes; i++) {
        E0->cd(i+1);
        Efficiency[i]->SetMaximum(1.1);
        Efficiency[i]->SetMinimum(0.);
        Efficiency[i]->Draw("");
        int in = i/N_classes;
        FakeRate[in]->SetMarkerColor(2);
        FakeRate[in]->SetLineColor(2);
        FakeRate[in]->Draw("SAME");
        int ic = i%N_classes;
        Eff_totL0[ic]->SetMarkerColor(3);
        Eff_totL0[ic]->SetLineColor(3);
        Eff_totL0[ic]->Draw("SAME");
        Efficiency[i]->Draw("SAME");
    }
 
    TCanvas * E1 = new TCanvas ("E1","", 800,800);
    E1->Divide(N_classes, N_neuronsL[1]);
    for (int i=N_neuronsL[0]*N_classes; i<N_neurons*N_classes; i++) {
        E1->cd(i+1-N_neuronsL[0]*N_classes);
        Efficiency[i]->SetMaximum(1.1);
        Efficiency[i]->SetMinimum(0.);
        Efficiency[i]->Draw("");
        int in = i/N_classes;
        FakeRate[in]->SetMarkerColor(2);
        FakeRate[in]->SetLineColor(2);
        FakeRate[in]->Draw("SAME");
        int ic = i%N_classes;
        Eff_totL1[ic]->SetMarkerColor(3);
        Eff_totL1[ic]->SetLineColor(3);
        Eff_totL1[ic]->Draw("SAME");
        Efficiency[i]->Draw("SAME");
    }

    // Plot the efficiencies and acceptances for the best q-value run
    // --------------------------------------------------------------
    for (int i=0; i<N_neurons*N_classes; i++) {
        int in = i/N_classes;
        int ic = i%N_classes;
        BestEff[in]->SetBinContent(ic+1,Efficiency[i]->GetBinContent(ind_qbest));
        BestFR[in]->SetBinContent(ic+1,FakeRate[in]->GetBinContent(ind_qbest));
        if (in<N_neuronsL[0]) {
            BestEtot[in]->SetBinContent(ic+1,Eff_totL0[ic]->GetBinContent(ind_qbest));
        } else {
            BestEtot[in]->SetBinContent(ic+1,Eff_totL1[ic]->GetBinContent(ind_qbest));
        }
    }
    TCanvas * BE = new TCanvas ("BE","",400,800);
    BE->Divide(1,N_neurons);
    for (int in=0; in<N_neurons; in++) {
        BE->cd(in+1);
        BestEff[in]->SetMaximum(1.1);
        BestEff[in]->SetMinimum(0.);
        BestEff[in]->Draw("");
        BestFR[in]->SetMarkerColor(2);
        BestFR[in]->SetLineColor(2);
        BestFR[in]->Draw("SAME");
        BestEtot[in]->SetMarkerColor(3);
        BestEtot[in]->SetLineColor(3);
        BestEtot[in]->Draw("SAME");
        BestEff[in]->Draw("SAME");
    }

    TCanvas * SE = new TCanvas ("SE","", 800, 400);
    SE->Divide(3,1);
    SE->cd(1);
    SelectivityL0->Draw();
    SE->cd(2);
    SelectivityL1->Draw();
    SE->cd(3);
    Qvalue->Draw();
    Qmax->Draw("SAME");

    TCanvas * W = new TCanvas ("W", "", 500, 800);
    int nrow = N_neurons/2;
    if (N_neurons%2!=0) nrow += 1;
    W->Divide(2, nrow);
    for (int in=0; in<N_neurons; in++) {
        W->cd(in+1);
        for (int is=0; is<N_streams; is++) {
            HWeight[in*N_streams+is]->SetMaximum(1.1);
            HWeight[in*N_streams+is]->SetMinimum(0.);
            int color = is+1;
            if (color==8) color=9;
            HWeight[in*N_streams+is]->SetLineColor(color);
            if (is==0) {
                HWeight[in*N_streams+is]->Draw();
            } else {
                HWeight[in*N_streams+is]->Draw("SAME");
            }
        }
    }

    // Draw final Efficiency and acceptance maps
    for (int in=0; in<N_neurons; in++) {
        for (int ic=0; ic<N_classes; ic++) {
            EffMap->SetBinContent(in+1,ic+1,Efficiency[ic+in*N_classes]->GetBinContent(ind_qbest));
        }
    }
    TCanvas * Y = new TCanvas ("Y","",600,900);
    Y->cd();
    EffMap->Draw("COL4");

    // Final Statistics
    // ----------------
    cout << endl << endl;
    cout << "         Run parameters" << endl;
    cout << "         -----------------------------------" << endl;
    cout << "                       L0 neurons: " << NL0 << endl;
    cout << "                       L1 neurons: " << NL1 << endl;
    cout << "            Connected L0-L1 frac.: " << CF01 << endl;
    cout << "            Connected IN-L0 frac.: " << CFI0 << endl;
    cout << "            Connected IN-L1 frac.: " << CFI1 << endl;
   
    cout << "                  Noise per layer: " << Occupancy*N_strips << endl;
    cout << "                    Track classes: " << N_cl << endl;
    cout << "                Only 8-hit tracks: ";
    if (!anyHits) {
        cout << "True" << endl;
    } else {
        cout << "False" << endl;
    }
    cout << endl;
    cout << "         Optimization results" << endl;
    cout << "         -----------------------------------" << endl;
    cout << "               Average efficiency: " << Eff_best<< endl;
    cout << "                Average fake rate: " << Acc_best << endl;
    cout << "                  Maximum Q value: " << Q_best << endl;
    cout << "                   L1 selectivity: " << SelL1_best << endl;
    cout << endl;
    cout << "         Optimized parameter values" << endl;
    cout << "         -----------------------------------" << endl;
    cout << "                     L0 threshold: " << T0_best<< endl;
    cout << "                     L1 threshold: " << T1_best << endl;
    cout << "                            alpha: " << A_best << endl;
    cout << "                        L1inhibit: " << L1if_best << endl;
    cout << "                                K: " << K_best << endl;
    cout << "                               K1: " << K1_best << endl;
    cout << "                               K2: " << K2_best << endl;
    cout << "                 IE pot. constant: " << IEPC_best << endl; 
    cout << "                 IPSP dt dilation: " << IPSP_dt_dilation << endl;
    cout << "         -----------------------------------" << endl;
    cout << endl;

    // Dump to file optimized parameters and results
    // ---------------------------------------------
    if (N_epochs>1) {
        Write_Parameters(); // This also defines indfile, used below
    }

    // Dump histograms to root file
    string Path = "./MODE/SNNT/";
    std::stringstream sstr;
    char num[40];
    sprintf (num,"NL0=%d_NL1=%d_NCl=%d_%d", N_neuronsL[0], N_neuronsL[1], N_classes, indfile);
    sstr << "Histos13_";
    string namerootfile = Path  + sstr.str() + num + ".root";
    TFile * rootfile = new TFile (namerootfile.c_str(),"RECREATE");
    rootfile->cd();

    // Write canvases first
    cout << "Writing canvas" << endl;
    S->Write();
    C->Write();
    E0->Write();
    E1->Write();
    BE->Write();
    SE->Write();
    W->Write();
    Y->Write();

    cout << "Saving pdf" << endl;
    S->SaveAs("./pdf/S.pdf" , "pdf");
    C->SaveAs("./pdf/C.pdf" , "pdf");
    E0->SaveAs("./pdf/E0.pdf", "pdf");
    E1->SaveAs("./pdf/E1.pdf", "pdf");
    BE->SaveAs("./pdf/BE.pdf", "pdf");
    SE->SaveAs("./pdf/SE.pdf", "pdf");
    W->SaveAs("./pdf/W.pdf" , "pdf");
    Y->SaveAs("./pdf/Y.pdf" , "pdf");

    // Then histograms
    SelectivityL0->Write();
    SelectivityL1->Write();
    Qvalue->Write();
    Qmax->Write();
    HEff->Write();
    HAcc->Write();
    HT0->Write();
    HT1->Write();
    HA->Write();
    HL1IF->Write();
    HK->Write();
    HK1->Write();
    HK2->Write();
    HIEPC->Write();
    HVoidWs->Write();
    HDelays->Write();
    Q_12->Write();
    Q_34->Write();
    Q_56->Write();
    Q_78->Write();
    Q_MV->Write();
    N_12->Write();
    N_34->Write();
    N_56->Write();
    N_78->Write();
    N_MV->Write();
    for (int in=0; in<N_neurons; in++) {
        for (int ic=0; ic<N_classes; ic++) {
            int id = in*N_classes+ic;
            Latency[id]->Write();
            Efficiency[id]->Write();
        }
        for (int is=0; is<N_streams; is++) {
            int id = in*N_streams+is;
            HWeight[id]->Write();
        }
        FakeRate[in]->Write();
        BestEff[in]->Write();
        BestFR[in]->Write();
        BestEtot[in]->Write();
    }
    for (int ic=0; ic<N_classes; ic++) {
        Eff_totL0[ic]->Write();
        Eff_totL1[ic]->Write();
    }
    for (int i=0; i<10; i++) {
        StreamsS[i]->Write();
        StreamsB[i]->Write();
        StreamsN[i]->Write();
    }
    EffMap->Write();

    rootfile->Write();

    // End of program
    rootfile->Close();
    gROOT->Time();

    return;
}
