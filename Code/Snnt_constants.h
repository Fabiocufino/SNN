#ifndef SNNT_CONSTANTS_H
#define SNNT_CONSTANTS_H

 
// hard coded constants and data used throughout the code 
// ------------------------------------------- 
static const double largenumber       = 999999999.;
static const double epsilon           = 1./largenumber;
static const int MaxNeurons           = 20;
static const int N_TrackingLayers     = 10;
static const int MaxStreams           = MaxNeurons + N_TrackingLayers;
static const int MaxStrips            = 1024;
static const int MaxEvents            = 10000000;
static const int MaxClasses           = 20;
static double TimeStep                = 0.001; // s
static int Empty_buffer               = 0;    // (old) We process events every 300 TimeSteps, leaving time for L0 neurons to pass delayed signal to L1 ones
static int N_strips                   = 1024;
static double pitch                   = 0.02;  // cm
static double pitch_rad               = 2*M_PI / N_strips;  // rad
static double width                   = 2.0;   // cm
static double Bfield                  = 3.8; // Tesla
static double ProbWSwitchUp           = 0.5;
static double ProbWSwitchDown         = 0.05;
static double MaxDelay                = 0.100;     // Determines shape of IE signal
static const double tau_m             = 0.01;      // membrane time constant
static const double tau_s             = 0.0025;    // synapse time constant
static double K1                      = 2.;        // constants to tune post-synaptic shape
static double K2                      = 4.;        // see above
static double IE_Pot_const            = 2.5;       // Constant for IE modeling
static double Threshold[2]            = {5.5,3.5}; // Neuron threshold in arbitrary units; in paper it is 550V but they have 1000 channels, 100x ours
static double alpha                   = 0.25;      // 0.25; // factor tuning inhibition strength
static double L1inhibitfactor         = 1.;        // multiplier for L1 inhibits
static const double MaxDeltaT         = 7.*tau_m;  // time window wherein pre-synaptic, inhibition, and post-synaptic kernels affect neuron potential
static const double tau_plus          = 0.0168;    // s
static const double tau_minus         = 0.0337;    // s
static double IPSP_dt_dilation        = 1.;        // shape factor of exponential IPSP signal
static const double a_plus            = 0.03125;   // for model of EPSP
static const double a_minus           = 0.0265625; // 0.85*a_plus;
static double MaxFactor               = 0.2;       // Initial factor of excursion of parameters for optimization
static double eff_target              = 0.9;
static double acc_target              = 0.05;
static double strip_r[N_TrackingLayers] = {3.0, 6.15, 10.45, 14.65, 24.9, 37.1678, 52.27, 68.7, 86.0, 108.3};//in cm (Sono posizioni assolute rispetto allo 0)
static double confidence_r              = 0.1;     // Confidence used to convert the radial position to the id of the strip in cm 
#endif // SNNT_CONSTANTS_H
