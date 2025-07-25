#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <ctime>
#include "oScanner.hh"
using namespace std;

#ifndef OTABLEBIOGRAPH2D_HH
#define OTABLEBIOGRAPH2D_HH 1

class oTableBiograph2D
{
  // Constructor & Destructor
  public:
    oTableBiograph2D( oScanner* fp_Scanner, int f_Verbose );
    ~oTableBiograph2D();

  // Member's functions
  public:
    int ComputeRingTables();
    void ElemView2CrystalIDs(int f_Elem, int f_View, int* f_Crystal1, int* f_Crystal2);

  // Get & Set functions
  public:
    inline int GetSpan()        {return m_Span;}
    inline int GetNbSino()      {return m_NbSino;}
    inline int GetMash()        {return m_Mash;}
    inline int GetMaxRingDiff() {return m_MaxRingDiff;}
    inline int GetNbElem()      {return m_NbElem;}
    inline int GetNbView()      {return m_NbView;}

    inline unsigned short int GetNbRingPairsBySinoIndex(int f_Sino) {return 1;};
    inline void GetRingPairBySinoIndex(int f_Sino, int f_RingPair, unsigned int* f_Ring1, unsigned int* f_Ring2) {f_Ring1[0]=0; f_Ring2[0]=0;};
    inline void GetRingPairsBySinoIndex(int f_Sino, unsigned int* f_Ring1, unsigned int* f_Ring2) {f_Ring1[0]=0; f_Ring2[0]=0;};
    inline int GetSinoIndexFromRingPair(unsigned int f_Ring1, unsigned int f_Ring2) {return 0;};
    void GetCrystalIDsFromElemView(int f_Elem, int f_View, int* f_Crystal1, int* f_Crystal2);
    void GetElemViewFromCrystalIDs(int* f_Elem, int* f_View, int f_Crystal1, int f_Crystal2);

    int GetLORsFromSinoBin( int f_BinSino, int f_BinView, int f_BinElem, unsigned int* fp_ID1, unsigned int* fp_ID2 );

  // Member's data
  private:
    // Lists for crystals indices to elem,view and vice-versa
    int*** mp_ElemView2CrystalID1;
    int*** mp_ElemView2CrystalID2;
    int** mp_CrystalIDs2Elem;
    int** mp_CrystalIDs2View;
    // Dimensions
    int m_NbElem;
    int m_NbView;
    int m_NbCrystalsPerRingWithGaps;
    int m_NbSino;
    int m_Span;
    int m_Mash;
    int m_MaxRingDiff;
    int m_Verbose;
    oScanner* mp_Scanner;
};

#endif

