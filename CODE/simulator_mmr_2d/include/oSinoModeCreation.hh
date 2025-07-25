#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <ctime>
#include "oTableBiograph.hh"
#include "oTableBiograph2D.hh"
#include "oTableHRplus.hh"
#include "oTableHRRT.hh"
#include "oTableInveon.hh"
#include "oTableFocus.hh"
#include "oTableCastor.hh"
#include "oTableMMR2D.hh"
#include "oTableVision600.hh"
#include "oScanner.hh"
#include "oMiscellaneous.hh"
using namespace std;

#ifndef OSINOMODECREATION_HH
#define OSINOMODECREATION_HH 1

#define ALL_WEIGHTING 1
#define POS_WEIGHTING 0

#define SINO_BIN_TAG -32000

#define ATTENUATION_FROM_UMAP 0
#define ATTENUATION_FROM_SINO 1

// Different word type
#define WORD_EOF -1
#define WORD_EVENT 0
#define WORD_OTHER 1

class oSinoModeCreation
{
  // Constructor & Destructor
  public:
    oSinoModeCreation( const string& f_ScannerName, int f_Weighting, bool f_SiddonDidier, int f_Seed, int f_Verbose );
    ~oSinoModeCreation();

  // Member's functions (unique)
  public:
    // For corrections
    int InitInputFiles( const string& f_FileFloat, const string& f_FileTrue, const string& f_FilePrompt, const string& f_FileListMode,
                        const string& f_FileRandom, const string& f_FileScatter,
                        const string& f_FileNormalization, const string& f_FileAttenuation, int f_AttenuationMode,
                        PRECISION f_TimeDelay, bool f_Mode2D );
    int InitFloat(const string& f_FileFloat);
    int InitNetTrue(const string& f_FileTrue);
    int InitPrompt(const string& f_FilePrompt);
    int InitListMode(const string& f_FileListMode);
    int InitScatter(string f_FileScat); // String name can be modified because two files separated by a comma can be provided
    int ApplyScatterFraction();
    int InitRandom(const string& f_FileRand);
    int InitNormalization(const string& f_FileNorm);
    int InitAttenuation(const string& f_FileAttn);
    int ComputeTotalScatter();
    int ComputeTotalRandom();
    float ComputeScatterRate(int f_Elem, int f_View, int f_Sino);
    float ComputeRandomRate(int f_Elem, int f_View, int f_Sino);
    float ComputeNormalizationFactor(int f_Elem, int f_View, int f_Sino);
    float ComputeACF(int f_CrystalID1, int f_CrystalID2);

    // General
    int InitConversionTable();
    int InitCrystalMap(const string& f_FileMap, bool f_Castor);
    int ProcessListMode(const string& f_FileBaseOut, bool f_FillEqualLORs, int f_Threads, bool f_Castor, PRECISION f_LimitedAnglePercent);

  // Private member's functions
  private:
    double SiddonAttenuationProjection( double x1, double y1, double z1, double x2, double y2, double z2);
    int ReadWord32bits(FILE* fin, int* f_Word);

  // Member's data
  private:
    // ----------------------------------------------------------------------------------------------------
    // Trues, scatters and normalisation sinograms
    bool m_Mode2D;
    bool m_InitBool;
    bool m_FloatBool;
    bool m_TrueBool;
    bool m_PromptBool;
    bool m_ScatCorr;
    bool m_RandCorr;
    bool m_NormCorr;
    float *** mp_FloatSino;
    short int *** mp_TrueSino;
    short int *** mp_PromptSino;
    float *** mp_ScatSino;
    float *** mp_RandSino;
    float *** mp_NormSino;
    int m_Weighting;
    int m_ReadSpan;
    int m_ReadNbSino;
    int m_ReadNbView;
    int m_ReadNbElem;
    int m_ReadSinoSize;
    int m_ReadMaxRingDiff;
    double m_ECF;
    float m_DeadTimeCorrectionFactor;
    int m_NbReplicates;
    string m_Isotope;
    int m_Frame;
    // ----------------------------------------------------------------------------------------------------
    // Attenuation (u-map or sino)
    bool m_AttnCorr;
    bool m_AttnHalfVoxelShift;
    int m_AttnDimX;
    int m_AttnDimY;
    int m_AttnDimZ;
    int m_AttnDimXY;
    int m_AttnDimTot;
    double m_AttnVoxSizeX;
    double m_AttnVoxSizeY;
    double m_AttnVoxSizeZ;
    float* mp_AttnUMap;
    float *** mp_AttnSino;
    int m_AttenuationMode;
    double m_AxisOrientationX;
    double m_AxisOrientationY;
    string m_Orientation;
    bool m_SiddonDidier;
    // ----------------------------------------------------------------------------------------------------
    // Scanner
    oScanner* mp_Scanner;
    // ----------------------------------------------------------------------------------------------------
    // Span and mash decompression
    oTableBiograph* mp_BiographTable;
    oTableBiograph2D* mp_BiographTable2D;
    oTableHRplus* mp_HRplusTable;
    oTableHRRT* mp_HrrtTable;
    oTableInveon* mp_InveonTable;
    oTableFocus* mp_FocusTable;
    oTableCastor* mp_CastorTable;
    oTableMMR2D* mp_MMRTable2D;
    oTableVision600* mp_Vision600Table;
    // ----------------------------------------------------------------------------------------------------
    // Counters
    double m_NbFloats;
    long int m_NbPrompts;
    long int m_NbDelays;
    long int m_NbNetTrues;
    double m_NbScatters;
    double m_ScatterFraction; // For the focus special case
    double m_NbRandoms;
    // ----------------------------------------------------------------------------------------------------
    // Time management
    int m_TimeDelay; // in milliseconds
    unsigned int m_StartTime; // in seconds
    unsigned int m_Duration;  // in seconds
    unsigned int m_RelativeStartTime; // in milliseconds
    unsigned int m_RelativeStopTime; // in milliseconds
    // ----------------------------------------------------------------------------------------------------
    // Verbosity
    int m_Verbose;
    unsigned int m_Seed;
};

#endif

