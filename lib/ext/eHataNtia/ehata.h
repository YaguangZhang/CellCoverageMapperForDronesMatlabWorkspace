// #pragma once
#ifdef _WIN32
// Export the DLL functions as "C" and not C++
// #define DLLEXPORT extern "C" __declspec(dllexport)
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#endif
#ifdef __linux__
// #define DLLEXPORT extern "C"
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#endif

#include <stdint.h>

typedef struct
{
    double d_bp__km;
    double att_1km;
    double att_100km;

    double h_b_eff__meter;
    double h_m_eff__meter;

    // Terrain Stats
    double pfl10__meter;
    double pfl50__meter;
    double pfl90__meter;
    double deltah__meter;

    // Path Geometry
    double d__km;
    double d_hzn__meter[2];
    double h_avg__meter[2];
    double theta_m__mrad;
    double beta;
    int iend_ov_sea;
    double hedge_tilda;
    bool single_horizon;

    // Misc
    double slope_max;
    double slope_min;
} InterValues;

#define PI 3.14159265358979323846

// public
double ExtendedHata(double pfl[], double f__mhz, double h_b__meter, double h_m__meter, int8_t enviro_code, double reliability);
void ExtendedHata_DBG(double pfl[], double f__mhz, double h_b__meter, double h_m__meter, int environment, double reliability, double *plb, InterValues *interValues);
