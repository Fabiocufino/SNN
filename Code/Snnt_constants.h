#ifndef SNNT_CONSTANTS_H
#define SNNT_CONSTANTS_H
#include <math.h>

// hard coded constants and data used throughout the code

// -------------------------------------------
// Time encoding constants
// -------------------------------------------

static int Empty_buffer = 0;                  // (old) We process events every 300 TimeSteps, leaving time for L0 neurons to pass delayed signal to L1 ones
static double delta = 0.7;                    // max delta for 1Gev is 0.66rad
static double max_angle = 2.0 * M_PI + delta; // max angle to scan
static double frequency = 40e6;               // CMS tracker reading frequency [Hz]
static double omega = max_angle * frequency;  // reading angular velocity

static double z_range = 1200.;                 // bin z in [-z_range/2, z_range/2]          
static double max_R = 1200.;                  // bin r in [0, max_R]
static const short int N_bin_r = 50;
static const int N_bin_z = 50;
static const int N_InputStreams = (N_bin_r)*N_bin_z;
static double z_bin_length = z_range / N_bin_z;    
static double r_bin_length = max_R / N_bin_r;       

// -------------------------------------------
// Hits and signals related constants
// -------------------------------------------
struct Hit //holds info about one hit (position of a cluster)
{
    float z;
    float r, phi;
    short int id;

    Hit(float r_, float z_, float phi_, short int id_)
        : r(r_), z(z_), phi(phi_), id(id_)
    {
    }
};

static vector<Hit> hit_pos;

static long int last_row_event_IT = 1;         // last row associated to the previous event read in IT
static long int last_row_event_OT = 1;         // last row associated to the previous event read in OT
static const int MaxClasses = 20;
static int pclass;

static const short int BGR = 1;
static const short int SIG = 2;

static long int NROOT = 100000;                //number of events inside the root file

// -------------------------------------------
// neural network constants
// -------------------------------------------
static const int MaxEvents = 10000000;
static const double largenumber = 999999999.;
static const double epsilon = 1. / largenumber;
static const int MaxNeurons = 20;
static double ProbWSwitchUp = 0.5;
static double ProbWSwitchDown = 0.05;
static double MaxDelay = 0.1e-9;         // Determines shape of IE signal[s]
static double tau_m = 1e-09;             // membrane time constant[s]
static double tau_s = 0.25e-09;          // synapse time constant[s]
static double K1 = 2.;                   // constants to tune post-synaptic shape
static double K2 = 4.;                   // see above
static double IE_Pot_const = 2.5;        // Constant for IE modeling
static double Threshold[2] = {15., 10.}; // Neuron threshold in arbitrary units; in paper it is 550V but they have 1000 channels, 100x ours
static double alpha = 0.25;              // 0.25; // factor tuning inhibition strength
static double L1inhibitfactor = 1.;      // multiplier for L1 inhibits
static double MaxDeltaT = 7. * tau_m;    // time window wherein pre-synaptic, inhibition, and post-synaptic kernels affect neuron potential
static double tau_plus = 1.68e-09;       // [s]
static double tau_minus = 3.37e-09;      // [s]
static double IPSP_dt_dilation = 1.;     // shape factor of exponential IPSP signal
static double a_plus = 0.03125;          // for model of EPSP
static double a_minus = 0.0265625;       // 0.85*a_plus;
static double MaxFactor = 0.2;           // Initial factor of excursion of parameters for optimization
static double eff_target = 0.9;
static double acc_target = 0.05;
static bool learnDelays = false;
static const int MaxStreams = MaxNeurons + N_bin_r * N_bin_z;

#endif // SNNT_CONSTANTS_H
