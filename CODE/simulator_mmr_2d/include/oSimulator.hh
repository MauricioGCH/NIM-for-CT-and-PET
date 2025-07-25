#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <ctime>
#include <vector>
#include <algorithm>
#include "oScanner.hh"
#include "oTableBiograph.hh"
#include "oTableBiograph2D.hh"
#include "oTableHRRT.hh"
#include "oTableHRplus.hh"
#include "oTableInveon.hh"
#include "oTableMMR2D.hh"
#include "oTableVision600.hh"
#include <mmintrin.h>
using namespace std;

#ifndef OSIMULATOR_HH
#define OSIMULATOR_HH 1

#define EVENT_TRUE 0
#define EVENT_RAND 1
#define EVENT_SCAT 2

class oSimulator
{
  // Constructor & Destructor
  public:
    oSimulator(int f_NbThreads, int f_Seed, int f_Verbose);
    ~oSimulator();

  // Member's functions
  public:
    int InitScannerStuff(const string& f_ScannerName, const string& f_FileCrystalMap, int f_Mash, int f_Span, int f_MaxRingDiff);
    int InitInputImages(const string& f_FileImage, const string& f_FileAttenuation, float f_OffsetX, float f_OffsetY, float f_OffsetZ, int f_Projector);
    int ApplyCounts(float f_ScatterFraction, float f_RandomFraction, float f_RandomFractionLSO, long int f_NbCounts, float f_ECF, bool f_ListMode, bool f_FillEqualLORs, int f_NbReplicates);
    int InitPSF(float f_PsfTransFWHM, float f_PsfAxialFWHM);
    int Project(bool f_SaveImage, int f_Projector);
    int SaveSinograms();

  // Member's private functions
  private:
    void SiddonForwardProjection( float x1, float y1, float z1, float x2, float y2, float z2,
                                  float* img, float* mumap,
                                  int dimX, int dimY, int dimZ, float voxX, float voxY, float voxZ,
                                  float* emission, float* transmission );
    void SiddonDidierForwardProjection( float x1, float y1, float z1, float x2, float y2, float z2,
                                        float* img, float* mumap,
                                        int dimX, int dimY, int dimZ, float voxX, float voxY, float voxZ,
                                        float* emission, float* transmission );
    void ComputeRandomFanSum();
    void ConvolveScatterComponent();
    void Make3DGaussianKernel();
    void Convolve3D(float* fp_Input, float* fp_Result);
    long int PoissonSampleLittle(float f_Mean, int f_Random);
    long int PoissonSampleLittle(float f_Mean);
    long int PoissonSampleBig(float f_Mean);

  // Member's data
  private:
    // ----------------------------------------------------------------------------------------------------
    // For constructor
    int m_NbThreads;
    int m_Verbose;
    unsigned int m_Seed;
    // For crystal map and scanner
    oScanner* mp_Scanner;
    oTableBiograph* mp_TableBiograph;
    oTableBiograph2D* mp_TableBiograph2D;
    oTableHRRT* mp_TableHRRT;
    oTableHRplus* mp_TableHRplus;
    oTableInveon* mp_TableInveon;
    oTableMMR2D* mp_TableMMR2D;
    oTableVision600* mp_TableVision600;
    // ----------------------------------------------------------------------------------------------------
    // For the image
    float* mp_InputImage; // The input image
    float* mp_MuMapImage; // The mumap in mm-1
    float* mp_ProjectImage; // Image to be projected
    float* mp_FanSum; // Singles counts in fan sum
    int m_DimX;
    int m_DimY;
    int m_DimZ;
    int m_DimXY;
    int m_DimTot;
    float m_VoxSizeX;
    float m_VoxSizeY;
    float m_VoxSizeZ;
    float m_FOVSizeX;
    float m_FOVSizeY;
    float m_FOVSizeZ;
    float m_OffsetX;
    float m_OffsetY;
    float m_OffsetZ;
    // ----------------------------------------------------------------------------------------------------
    // Counts
    double m_TotalTrue;
    double m_TotalScat;
    double m_TotalRand;
    long int m_CountTrue;
    long int m_CountScat;
    long int m_CountRand;
    long int m_CountPrompt;
    int m_NbReplicates;
    // ----------------------------------------------------------------------------------------------------
    // Corrections
    bool m_ScatCorr;
    bool m_RandCorr;
    bool m_AttnCorr;
    double m_ECF;
    // ----------------------------------------------------------------------------------------------------
    // For the PSF
    float m_PsfTransFWHM; // The PSF transaxial FWHM (gaussian 3D symetric)
    float m_PsfAxialFWHM; // The PSF axial FWHM (gaussian 3D symetric)
    int m_PsfKernSizeX; // The PSF kernel size along axis X
    int m_PsfKernSizeY; // The PSF kernel size along axis Y
    int m_PsfKernSizeZ; // The PSF kernel size along axis Z
    float *** mp_PsfKernel; // The kernel allocated in 3D
    // ----------------------------------------------------------------------------------------------------
    // For the sinograms
    bool m_FloatBool;     // Say if we stay in float (do not add Poisson noise)
    float*** mp_SinoForw; // The projected sinogram
    float*** mp_SinoScat; // The scatter sinogram
    float*** mp_SinoRand; // The random sinogram
    float*** mp_SinoNorm; // The normalization sinogram
    float*** mp_SinoAttn; // The attenuation sinogram
    int**** mp_SinoShoot; // Random numbers (for thread-safe and repeatable shooting)
    short int*** mp_SinoPrompt; // The final sinogram
    int m_NbSino;
    int m_NbView;
    int m_NbElem;
    int m_Mash;
    int m_Span;
    int m_MaxRingDiff;
    // ----------------------------------------------------------------------------------------------------
    // For the correct algorithm initialization
    bool m_HaveInitScannerStuff;
    bool m_HaveInitInputImage;
    bool m_HaveInitPSF;
    bool m_HaveProject;
};

#endif

