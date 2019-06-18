#pragma once
#ifdef _WIN32
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#endif
#ifdef __linux__
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#endif

#include <stdint.h>

struct InterValues
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
};

#define PI 3.14159265358979323846

// public
double ExtendedHata(double pfl[], double f__mhz, double h_b__meter, double h_m__meter, int8_t environment, double reliability);
void ExtendedHata_DBG(double pfl[], double f__mhz, double h_b__meter, double h_m__meter, int environment, double reliability, double *plb, struct InterValues *interValues);

// private
void FindAverageGroundHeight(double *pfl, double h_m__meter, double h_b__meter, struct InterValues *interValues);
void MobileTerrainSlope(double *pfl, struct InterValues *interValues);
void LeastSquares(double *pfl_segment, double x1, double x2, double *z0, double *zn);
void AnalyzeSeaPath(double* pfl, struct InterValues *interValues);
void FindHorizons(double *pfl, double gme, double d__meter, double h_1__meter, double h_2__meter, double *d_hzn__meter);
void SingleHorizonTest(double *pfl, double h_m__meter, double h_b__meter, struct InterValues *interValues);
void ComputeTerrainStatistics(double *pfl, struct InterValues *interValues);
double FindQuantile(const int nn, double *apfl, const int ir);
void PreprocessTerrainPath(double *pfl, double h_b__meter, double h_m__meter, struct InterValues *interValues);
double AverageTerrainHeight(double *pfl);
double GeneralSlopeCorrectionFactor(double theta_m__mrad, double d__km);
double FineRollingHillyTerrainCorrectionFactor(struct InterValues *interValues, double h_m_gnd__meter);
double MixedPathCorrectionFactor(double d__km, struct InterValues *interValues);
double MedianRollingHillyTerrainCorrectionFactor(double deltah);
void MedianBasicPropLoss(double f__mhz, double h_b__meter, double h_m__meter, double d__km, int environment, double* plb_med__db, struct InterValues *interValues);
double IsolatedRidgeCorrectionFactor(double d1_hzn__km, double d2_hzn__km, double h_edge__meter);
double Variability(double plb_med__db, double f__mhz, int enviro_code, double reliability);
double Sigma_u(double f__mhz);
double Sigma_r(double f__mhz);
double ierf(double q);

