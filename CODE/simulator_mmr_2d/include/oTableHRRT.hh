#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <ctime>
#include "oScanner.hh"
using namespace std;

#ifndef OTABLEHRRT_HH
#define OTABLEHRRT_HH 1

typedef struct
{
  int nsino;
  float d;
  float z;
} SOL;

#define NDOIS   2
#define NXCRYS  72
#define NYCRYS  104
#define NHEADS  8
#define NCRYS 119808
#define NTRANS 576
#define NCRYSINLAYER 59904

class oTableHRRT
{
  // Constructor & Destructor
  public:
    oTableHRRT( oScanner* fp_Scanner, int f_Verbose );
    ~oTableHRRT();

  // Member's functions
  public:
    int ReadLUTTable(string f_ScannerTable="");
    int ComputeSpanTables();
    int GetSinogramBinFromCrystalIndices( int f_Mp, int f_Layer1, int f_Layer2, int f_Trans1, int f_Trans2, int f_Axial1, int f_Axial2, int* f_BinSino, int* f_BinView, int* f_BinElem );
    int GetHeadIndicesFromHeadPairIndex( int f_Mp, int* f_Head1, int* f_Head2 );
    int GetHeadPairIndexFromHeadIndices( int f_Head1, int f_Head2, int* f_Mp );

    int ComputeFinalTables();
    int ComputeElemViewRingTables();
    int ComputeSino2PlaneSegmentTables();
    int GetLORsFromSinoBin( int f_BinSino, int f_BinView, int f_BinElem, unsigned int* fp_ID1, unsigned int* fp_ID2 );
    int ComputeFullTable();


  // Get & Set functions
  public:
    inline int GetMash()        {return m_Mash;}
    inline int GetSpan()        {return m_Span;}
    inline int GetMaxRingDiff() {return m_MaxRingDiff;}
    inline int GetNbElem()      {return m_NbElem;}
    inline int GetNbView()      {return m_NbView;}
    inline int GetNbSino()      {return m_NbSino;}

  // Member's data
  private:
    int** mp_HeadPairs;
    int** mp_HeadPairsReversed;
    int* mp_Span3to9ConversionTable;
    SOL*** mp_LUTsino;
    float* mp_LUTzpos;
    float* mp_LUTzpos2;
    int* m_segz0;
    int* m_segzmax;
    int* m_segzoffset;
    int m_nsegs;
    short **m_segplane;
    int m_MaxRingDiff;
    int m_NbElem;
    int m_NbView;
    int m_NbSino;
    int m_Span;
    int m_Mash;
    string m_TableFileName;
    oScanner* mp_Scanner;
    int m_Verbose;
    // Final ring conversion tables
    int*** mp_ElemView2Mp;
    int*** mp_ElemView2XX1;
    int*** mp_ElemView2XX2;
    int**  mp_ElemView2NbLors;
    // Sino 2 plane segment tables
    int* mp_Sino2Segment;
    int* mp_Sino2Plane;
    // Full table
    int**** mp_Bin2CrystalID1;
    int**** mp_Bin2CrystalID2;
    int***  mp_Bin2NbLors;
};

#endif

