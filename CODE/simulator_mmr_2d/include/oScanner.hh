#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <ctime>
using namespace std;

#ifndef OSCANNER_HH
#define OSCANNER_HH 1

// Some predefined scanners
#define SCANNER_HRRT        0
#define SCANNER_BIOGRAPH    1
#define SCANNER_HRPLUS      2
#define SCANNER_BIOGRAPH2D  3
#define SCANNER_INVEON      4
#define SCANNER_FOCUS       5
#define SCANNER_CASTOR      6
#define SCANNER_SIGNA       7
#define SCANNER_MMR2D       8
#define SCANNER_VISION600   9
class oScanner
{
  // Constructor & Destructor
  public:
    oScanner( const string& f_ScannerName, int f_Verbose );
    oScanner( int f_ScannerModel, int f_Verbose );
    void GenericConstructor( const string& f_ScannerName, int f_Verbose );
    ~oScanner();
    void Initialize( int f_NbHeads, int f_NbAxialBlocks, int f_NbTransBlocks,
                     int f_NbAxialCrystals, int f_NbTransCrystals, int f_NbLayers,
                     int f_NbGapsBetweenHeads, int f_NbAxialGapsBetweenBlocks, int f_NbTransGapsBetweenBlocks,
                     float f_Radius, float f_MeanDOI,
                     float f_AxialCrystalSize, float f_TransCrystalSize, float f_DepthCrystalSize,
                     float f_AxialCrystalGap, float f_TransCrystalGap,
                     float f_AxialBlockGap, float f_TransBlockGap,
                     int f_NbElem, int f_NbView, float f_RadialBinSize,
                     float f_FOVSizeX, float f_FOVSizeY, float f_FOVSizeZ,
                     float f_BedOffset,
                     int f_NbVoxelX, int f_NbVoxelY, int f_NbVoxelZ,
                     int f_ScannerModel, int f_BitMode, int f_Verbose
                   );
    int ComputeNbSinoFromSpanAndMaxRingDiff(int f_Span, int f_MaxRingDiff);

  // Get & Set functions
  public:
    inline float GetCornerX1(int f_CrystalID) {return mp_CornerX1[f_CrystalID];}
    inline float GetCornerY1(int f_CrystalID) {return mp_CornerY1[f_CrystalID];}
    inline float GetCornerZ1(int f_CrystalID) {return mp_CornerZ1[f_CrystalID];}
    inline float GetCornerX2(int f_CrystalID) {return mp_CornerX2[f_CrystalID];}
    inline float GetCornerY2(int f_CrystalID) {return mp_CornerY2[f_CrystalID];}
    inline float GetCornerZ2(int f_CrystalID) {return mp_CornerZ2[f_CrystalID];}
    inline float GetEfficiency(int f_CrystalID) {return mp_Efficiency[f_CrystalID];}
    inline int GetScannerModel() {return m_ScannerModel;}
    inline int GetNbView() {return m_NbView;}
    inline int GetNbElem() {return m_NbElem;}
    inline int GetNbSinoBin() {return m_NbSinoBin;}
    inline int GetNbPlanes() {return m_NbPlanes;}
    inline int GetNbAxialCrystals() {return m_NbAxialCrystals;}
    inline int GetNbTotalAxialCrystals() {return m_NbTotalAxialCrystals;}
    inline int GetNbTotalAxialCrystalsWithGaps() {return m_NbTotalAxialCrystalsWithGaps;}
    inline int GetNbTotalTransCrystals() {return m_NbTotalTransCrystals;}
    inline int GetNbTotalTransCrystalsWithGaps() {return m_NbTotalTransCrystalsWithGaps;}
    inline int GetNbTransCrystalsInHead() {return m_NbTransCrystalsInHead;}
    inline int GetNbTransCrystals() {return m_NbTransCrystals;}
    inline int GetNbTransBlocks() {return m_NbTransBlocks;}
    inline int GetNbAxialGapsBetweenBlocks() {return m_NbAxialGapsBetweenBlocks;}
    inline int GetNbTransGapsBetweenBlocks() {return m_NbTransGapsBetweenBlocks;}
    inline int GetNbGapsBetweenHeads() {return m_NbGapsBetweenHeads;}
    inline float GetAxialScannerSize() {return m_AxialScannerSize;}
    inline int GetBitMode() {return m_BitMode;}
    inline int GetNbTotalCrystals() {return m_NbTotalCrystals;}
    inline float GetAxialCrystalSize() {return m_AxialCrystalSize;}
    inline int GetNbLayers() {return m_NbLayers;}
    inline float GetFOVSizeX() {return m_FOVSizeX;}
    inline float GetFOVSizeY() {return m_FOVSizeY;}
    inline float GetFOVSizeZ() {return m_FOVSizeZ;}
    inline int GetNbVoxelX() {return m_NbVoxelX;}
    inline int GetNbVoxelY() {return m_NbVoxelY;}
    inline int GetNbVoxelZ() {return m_NbVoxelZ;}
    inline float GetRadialBinSize() {return m_RadialBinSize;}
    inline float GetBedOffset() {return m_BedOffset;}

  // Member's functions
  private:
    int ComputeDiagCoordinates( int f_LayerIndex, int f_AxialBlockIndex, int f_AxialCrystalIndex,
                                int f_HeadIndex, int f_TransBlockIndex, int f_TransCrystalIndex,
                                float* fp_Corner1X, float* fp_Corner1Y, float* fp_Corner1Z,
                                float* fp_Corner2X, float* fp_Corner2Y, float* fp_Corner2Z, bool f_Castor
                              );
    int ComputeOrientationVector( int f_LayerIndex, int f_AxialBlockIndex, int f_AxialCrystalIndex,
                                  int f_HeadIndex, int f_TransBlockIndex, int f_TransCrystalIndex,
                                  float* fp_OrientationX, float* fp_OrientationY, float* fp_OrientationZ,
                                  bool f_Castor
                                );
  public:
    int ProcessCrystalEfficiency( bool f_Uniform, bool f_Random, float f_RandomEfficiencyPercentDeviation );
    int ProcessDeadBlocks(int f_NbDeadBlocks, int* f_DeadBlockHeadIndex, int* f_DeadBlockAxialIndex, int* f_DeadBlockTransIndex);
    int WriteCrystalEfficiency( const string& f_FileBaseOut );
    int WriteCrystalMap( const string& f_FileBaseOut, bool f_Castor );
    int ReadCrystalMap( const string& f_MapFileName );
    int ComputeUniformCrystalMap(bool f_Castor);

  // Member's data
  private:
    // Scanner model
    int m_ScannerModel;
    int m_BitMode; // 32 or 64 bits
    // Dimensions
    int m_NbHeads;         // Number of heads in the ring
    int m_NbAxialBlocks;   // Number of blocks axially inside a head
    int m_NbTransBlocks;   // Number of blocks transaxially inside a head
    int m_NbAxialCrystals; // Number of crystals axially inside a block
    int m_NbTransCrystals; // Number of crystals transaxially inside a block
    int m_NbLayers;        // Number of layers
    int m_NbTotalTransCrystals; // Total number of crystals in the transaxial direction (along the ring)
    int m_NbTotalAxialCrystals; // Total number of crystals axially (along the tomograph axis)
    int m_NbTotalCrystals;      // Total number of crystals
    int m_NbGapsBetweenHeads;       // Number of transaxial gaps between heads
    int m_NbAxialGapsBetweenBlocks; // Number of axial gaps between blocks
    int m_NbTransGapsBetweenBlocks; // Number of transaxial gaps between blocks in a head
    int m_NbTotalTransCrystalsWithGaps;    // Number of crystals per ring including gaps
    int m_NbTotalAxialCrystalsWithGaps;    // Number of axial crystals including gaps
    int m_NbTotalCrystalsWithGaps; // Total number of crystals including gaps
    int m_NbTransCrystalsInHead;        // Number of crystals transaxially in a head
    float m_Radius;           // Radius from one head to another, face to face
    float m_MeanDOI;          // Mean depth of interaction in the crystals
    float m_AxialCrystalSize; // Size of a crystal along the tomograph axis
    float m_TransCrystalSize; // Size of a crystal transaxially
    float m_DepthCrystalSize; // Depth of a crystal (including all layers)
    float m_AxialCrystalGap;  // Axial gap between two crystals
    float m_TransCrystalGap;  // Transaxial gap between two crystals
    float m_AxialBlockSize;   // Size of a block along the tomograph axis
    float m_TransBlockSize;   // Size of a block transaxially
    float m_AxialBlockGap;    // Axial gap between two blocks
    float m_TransBlockGap;    // Transaxial gap between two blocks
    float m_TransHeadSize;    // Total transaxial size of the head
    float m_AxialScannerSize; // Total axial scanner size
    float m_RadialBinSize;    // Size in mm of a radial bin
    // Bed offset for whole-body imaging
    float m_BedOffset;
    // Image default size
    float m_FOVSizeX;
    float m_FOVSizeY;
    float m_FOVSizeZ;
    int m_NbVoxelX;
    int m_NbVoxelY;
    int m_NbVoxelZ;
    // Sinogram dimensions
    int m_NbElem;      // Number of radial elements
    int m_NbView;      // Number of views (angles) in a 2D sinogram
    int m_NbSinoBin;   // Number of sinogram bins (elem*view)
    int m_NbPlanes;    // Number of direct planes
    // Crystal efficiency
    float*** mp_CrystalEfficiency; // Table of crystal efficiencies for calculation
    // Variables for reading
    string m_MapFileName;
    float* mp_CornerX1;
    float* mp_CornerY1;
    float* mp_CornerZ1;
    float* mp_CornerX2;
    float* mp_CornerY2;
    float* mp_CornerZ2;
    float* mp_Efficiency;
    // Verbosity
    int m_Verbose;
};

#endif

