#include <math.h>
#include <stdlib.h>

#include <stdint.h>

// Export the DLL functions as "C" and not C++
// #define DLLEXPORT extern "C" __declspec(dllexport)
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define DIM(x, y) (((x) > (y)) ? (x - y) : (0))

#define PI                                      3.1415926535897932384
#define SQRT2                                   sqrt(2)
#define a_0__meter                              6370e3
#define a_9000__meter                           9000e3
#define THIRD                                   1.0 / 3.0

#define MODE__P2P                               0
#define MODE__AREA                              1


/////////////////////////////
// Data Structures

typedef struct
{
    double theta_hzn[2];        // Terminal horizon angles
    double d_hzn__meter[2];     // Terminal horizon distances, in meters
    double h_e__meter[2];       // Terminal effective heights, in meters
    double N_s;                 // Surface refractivity, in N-Units
    double delta_h__meter;      // Terrain irregularity parameter, in meters
    double A_ref__db;           // Reference attenuation, in dB
    double A_fs__db;            // Free space basic transmission loss, in dB
    double d__km;               // Path distance, in km
    int mode;                   // Mode of propagation value
} IntermediateValues;

/////////////////////////////
// Custom Functions
double ITM_P2P_TLS_PL(double h_tx__meter, double h_rx__meter, double pfl[], int8_t climate, double N_0, double f__mhz,
    int8_t pol, double epsilon, double sigma, int8_t mdvar, double time, double location, double situation);

/////////////////////////////
// Main ITM Functions

int ITM_P2P_TLS(double h_tx__meter, double h_rx__meter, double pfl[], int climate, double N_0, double f__mhz,
    int pol, double epsilon, double sigma, int mdvar, double time, double location, double situation,
    double *A__db, long *warnings);
int ITM_P2P_TLS_Ex(double h_tx__meter, double h_rx__meter, double pfl[], int climate, double N_0, double f__mhz,
    int pol, double epsilon, double sigma, int mdvar, double time, double location, double situation,
    double *A__db, long *warnings, IntermediateValues *interValues);
int ITM_P2P_CR(double h_tx__meter, double h_rx__meter, double pfl[], int climate, double N_0, double f__mhz,
    int pol, double epsilon, double sigma, int mdvar, double confidence, double reliability,
    double *A__db, long *warnings);
int ITM_P2P_CR_Ex(double h_tx__meter, double h_rx__meter, double pfl[], int climate, double N_0, double f__mhz,
    int pol, double epsilon, double sigma, int mdvar, double confidence, double reliability,
    double *A__db, long *warnings, IntermediateValues *interValues);
int ITM_AREA_TLS(double h_tx__meter, double h_rx__meter, int tx_site_criteria, int rx_site_criteria, double d__km,
    double delta_h__meter, int climate, double N_0, double f__mhz, int pol, double epsilon, double sigma,
    int mdvar, double time, double location, double situation, double *A__db, long *warnings);
int ITM_AREA_TLS_Ex(double h_tx__meter, double h_rx__meter, int tx_site_criteria, int rx_site_criteria, double d__km,
    double delta_h__meter, int climate, double N_0, double f__mhz, int pol, double epsilon, double sigma,
    int mdvar, double time, double location, double situation, double *A__db, long *warnings, IntermediateValues *interValues);
int ITM_AREA_CR(double h_tx__meter, double h_rx__meter, int tx_site_criteria, int rx_site_criteria, double d__km,
    double delta_h__meter, int climate, double N_0, double f__mhz, int pol, double epsilon, double sigma,
    int mdvar, double confidence, double reliability, double *A__db, long *warnings);
int ITM_AREA_CR_Ex(double h_tx__meter, double h_rx__meter, int tx_site_criteria, int rx_site_criteria, double d__km,
    double delta_h__meter, int climate, double N_0, double f__mhz, int pol, double epsilon, double sigma,
    int mdvar, double confidence, double reliability, double *A__db, long *warnings, IntermediateValues *interValues);