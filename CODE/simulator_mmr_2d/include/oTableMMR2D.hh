#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <ctime>
#include "oScanner.hh"
using namespace std;

#ifndef OTABLEMMR2D_HH
#define OTABLEMMR2D_HH 1

class oTableMMR2D
{
  // Constructor & Destructor
  public:
    oTableMMR2D( int f_NbElem, int f_NbView, int f_NbSino, int f_MaxRingDiff, int f_Span, oScanner* fp_Scanner, int f_Verbose );
    oTableMMR2D( oScanner* f_Scanner, int f_Verbose );
    ~oTableMMR2D();

  // Member's functions
  public:
    inline int ComputeSpanTables() {return 0;};
    int ComputeRingTables();
    void ElemView2CrystalIDs(int f_Elem, int f_View, int* f_Crystal1, int* f_Crystal2);

  // Get & Set functions
  public:
    inline int GetMash()   {return m_Mash;}
    inline int GetSpan()   {return 1;}
    inline int GetNbSino() {return 1;}
    inline int GetNbView() {return m_NbView;}
    inline int GetNbElem() {return m_NbElem;}

    void GetCrystalIDsFromElemView(int f_Elem, int f_View, int* f_Crystal1, int* f_Crystal2);
    void GetCrystalIDFromElemView(int f_Index, int f_Elem, int f_View, int* f_Crystal1, int* f_Crystal2);
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
    int m_Mash;
    int m_Span;
    int m_MaxRingDiff;
    int m_NbSino;
    int m_NbCalculatedSino;
    int m_NbSegments;
    int m_NbAbsoluteSegments;
    oScanner* mp_Scanner;
    int m_Verbose;
};

#endif

