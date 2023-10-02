#ifndef SNNT_CONSTANTS_H
#define SNNT_CONSTANTS_H
#include <math.h>

// hard coded constants and data used throughout the code 
// ------------------------------------------- 

// ------------------------------------------- 
//Tracker geometry and physical constants
// ------------------------------------------- 
static float Bfield                  = 3.8;                              // Tesla
static const short int N_TrackingLayers     = 50;
static float strip_r[N_TrackingLayers] = {32, 64, 107, 149, 250, 370, 525, 670, 710, 845, 876, 1068, 1099};//radial positions of the layers from the origin
//static float strip_r[N_TrackingLayers] = {3.0, 6.15, 10.45, 14.65, 24.9, 37.1678, 52.27, 68.7, 86.0, 108.3};//radial positions of the layers from 
//static float confidence_r_left[N_TrackingLayers] = {0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1};
//static float confidence_r_right[N_TrackingLayers] = {0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1};
static float confidence_r_left[N_TrackingLayers] = {7, 6, 7, 8, 40, 40, 30, 10, 10, 5, 6, 5, 4};
static float confidence_r_right[N_TrackingLayers] = {7, 6, 7, 8, 40, 40, 30, 10, 10, 5, 6, 5, 4};
static const int MaxStrips            = 10000;                            //Discretization of the phi angle (retained only for noise calculation)
static int N_strips                   = 10000;
static float pitch_rad               = 2.0 * M_PI / N_strips;            //[rad]
static float confidence_r            = 0.1;                              // Confidence used to convert the radial position to the id of the strip in cm 
static const int N_EtaLayers              = 50;
static const int N_InputStreams = (N_TrackingLayers) * N_EtaLayers ;
// ------------------------------------------- 
//Time encoding constants
// ------------------------------------------- 
static int Empty_buffer               = 0;                                // (old) We process events every 300 TimeSteps, leaving time for L0 neurons to pass delayed signal to L1 ones
static double delta     = 0.7;                                            //max delta for 1Gev is 0.66rad
static double max_angle = 2.0 * M_PI + delta;                             //max angle to scan
static double frequency = 40e6;                                           //CMS tracker reading frequency [Hz]
static double omega = max_angle * frequency;                              //reading angular velocity

// ------------------------------------------- 
//neural network constants
// ------------------------------------------- 
static const int MaxEvents            = 10000000;
static const double largenumber       = 999999999.;
static const double epsilon           = 1./largenumber;
static const int MaxNeurons           = 20;
static double ProbWSwitchUp           = 0.5;
static double ProbWSwitchDown         = 0.05;
static double MaxDelay                = 0.1e-9;                          // Determines shape of IE signal[s]
static double tau_m             = 1e-09;                           // membrane time constant[s]
static  double tau_s             = 0.25e-09;                        // synapse time constant[s]
static double K1                      = 2.;                              // constants to tune post-synaptic shape
static double K2                      = 4.;                              // see above
static double IE_Pot_const            = 2.5;                             // Constant for IE modeling
static double Threshold[2]            = {15.,10.};                       // Neuron threshold in arbitrary units; in paper it is 550V but they have 1000 channels, 100x ours
static double alpha                   = 0.25;                            // 0.25; // factor tuning inhibition strength
static double L1inhibitfactor         = 1.;                              // multiplier for L1 inhibits
static double MaxDeltaT         = 7.*tau_m;                        // time window wherein pre-synaptic, inhibition, and post-synaptic kernels affect neuron potential
static double tau_plus          = 1.68e-09;                        // [s]
static double tau_minus         = 3.37e-09;                        // [s]
static double IPSP_dt_dilation        = 1.;                              // shape factor of exponential IPSP signal
static double a_plus            = 0.03125;                         // for model of EPSP
static double a_minus           = 0.0265625;                       // 0.85*a_plus;
static double MaxFactor               = 0.2;                             // Initial factor of excursion of parameters for optimization it was 0.2
static double eff_target              = 0.9;
static double acc_target              = 0.05;
static bool learnDelays = false;
static const int MaxStreams           = MaxNeurons + N_TrackingLayers*N_EtaLayers;

// ------------------------------------------- 
//Other constants
// -------------------------------------------     
struct Hit
{
    float eta;
    float r, phi;
    short int id;

    Hit(float r_, float eta_, float phi_, short int id_)
        : r(r_), eta(eta_), phi(phi_), id(id_)
    {
    }
};

static vector<Hit> hit_pos;

static float bisection_window        = (M_PI/5);                        //interval of confidence used in the bisection method to find the interceptions between layers and tracks
static float bisection_precision     = 1.e-04;                          //[rad]
static long int last_row_event        = 1;                               //last row associated to the previous event read
static long int last_row_event_OT     = 1;
static const int MaxClasses           = 20;
static int pclass;

static const short int BGR = 1;
static const short int SIG = 2; 

ofstream fout_n("n.txt");
ofstream fout_h_bg("h_bg.txt");
ofstream fout_h_sig("h_sig.txt");


static double eta_range = 4;
static double eta_strip = eta_range / N_EtaLayers;
static double max_R = 1200.;
static double r_bin_length = max_R/N_TrackingLayers;

#endif // SNNT_CONSTANTS_H
