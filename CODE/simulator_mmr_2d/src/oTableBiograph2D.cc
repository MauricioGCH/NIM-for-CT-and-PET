#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <ctime>
#include <fstream>
#include "oTableBiograph2D.hh"
#include "oOutputManager.hh"
#include "oScanner.hh"
using namespace std;

// ==========================================================================================================================================
// Constructor
// ==========================================================================================================================================
oTableBiograph2D::oTableBiograph2D( oScanner* fp_Scanner, int f_Verbose )
{
  // Verbose
  m_Verbose = f_Verbose;
  if (m_Verbose>=1) LogCout ("oTableBiograph2D::Constructor() -> Initialize biograph 2D conversion table" << endl);

  // Sinogram dimensions and characteristics
  m_NbSino = 1; // Imposed by an eletronic card inside the scanner !
  m_Span = 1;    // Imposed by an eletronic card inside the scanner !
  m_NbElem = 336; // Imposed by an eletronic card inside the scanner !
  m_NbView = 336; // Imposed by an eletronic card inside the scanner !
  m_NbCrystalsPerRingWithGaps = 14*48;
  m_Mash = 1;     // Imposed by an eletronic card inside the scanner !
  m_MaxRingDiff = 1; // Imposed by an eletronic card inside the scanner !

  // Scanner model
  mp_Scanner = fp_Scanner;
}
// ==========================================================================================================================================
// Destructor
// ==========================================================================================================================================
oTableBiograph2D::~oTableBiograph2D()
{
  if (m_Verbose>=1) LogCout ("oTableBiograph2D::Destructor() -> Destroy conversion table" << endl);

}
// ==========================================================================================================================================
// int ComputeRingTables()
//   --> This function computes the ring conversion tables (crystals to elem,view) and returns 0 upon success and other value otherwise.
// ==========================================================================================================================================
int oTableBiograph2D::ComputeRingTables()
{
  if (m_Verbose>=1) LogCout ("oTableBiograph2D::ComputeRingTables() -> Compute the ring conversion tables" << endl);

  // Allocations
  mp_ElemView2CrystalID1 = (int***)malloc(m_NbElem*sizeof(int**));
  mp_ElemView2CrystalID2 = (int***)malloc(m_NbElem*sizeof(int**));
  for (int e=0; e<m_NbElem; e++)
  {
    mp_ElemView2CrystalID1[e] = (int**)malloc(m_NbView*sizeof(int*));
    mp_ElemView2CrystalID2[e] = (int**)malloc(m_NbView*sizeof(int*));
    for (int v=0; v<m_NbView; v++)
    {
      mp_ElemView2CrystalID1[e][v] = (int*)malloc(m_Mash*sizeof(int));
      mp_ElemView2CrystalID2[e][v] = (int*)malloc(m_Mash*sizeof(int));
    }
  }
  int nb_total_trans_crystals = mp_Scanner->GetNbTotalTransCrystalsWithGaps();
  mp_CrystalIDs2Elem = (int**)malloc(nb_total_trans_crystals*sizeof(int*));
  mp_CrystalIDs2View = (int**)malloc(nb_total_trans_crystals*sizeof(int*));
  for (int c1=0; c1<nb_total_trans_crystals; c1++)
  {
    mp_CrystalIDs2Elem[c1] = (int*)malloc(nb_total_trans_crystals*sizeof(int));
    mp_CrystalIDs2View[c1] = (int*)malloc(nb_total_trans_crystals*sizeof(int));
    for (int c2=0; c2<nb_total_trans_crystals; c2++)
    {
      mp_CrystalIDs2Elem[c1][c2] = -1;
      mp_CrystalIDs2View[c1][c2] = -1;
    }
  }

  // Loop on all elem,view indices
  int* crystals1 = (int*)malloc(m_Mash*sizeof(int));
  int* crystals2 = (int*)malloc(m_Mash*sizeof(int));
  for (int e=0; e<m_NbElem; e++) for (int v=0; v<m_NbView; v++)
  {
    ElemView2CrystalIDs(e, v, crystals1, crystals2);
    for (int m=0; m<m_Mash; m++)
    {
      mp_ElemView2CrystalID1[e][v][m] = crystals1[m];
      mp_ElemView2CrystalID2[e][v][m] = crystals2[m];
      mp_CrystalIDs2Elem[crystals1[m]][crystals2[m]] = e;
      mp_CrystalIDs2View[crystals1[m]][crystals2[m]] = v;
      mp_CrystalIDs2Elem[crystals2[m]][crystals1[m]] = e;
      mp_CrystalIDs2View[crystals2[m]][crystals1[m]] = v;
    }
  }
  free(crystals1);
  free(crystals2);

  // Ending
  return 0;
}
// ==========================================================================================================================================
// Set & Get functions to access the table informations
// ==========================================================================================================================================
void oTableBiograph2D::GetCrystalIDsFromElemView(int f_Elem, int f_View, int* f_Crystal1, int* f_Crystal2)
{
  *f_Crystal1 = mp_ElemView2CrystalID1[f_Elem][f_View][0];
  *f_Crystal2 = mp_ElemView2CrystalID2[f_Elem][f_View][0];
}
void oTableBiograph2D::GetElemViewFromCrystalIDs(int* f_Elem, int* f_View, int f_Crystal1, int f_Crystal2)
{
  *f_Elem = mp_CrystalIDs2Elem[f_Crystal1][f_Crystal2];
  *f_View = mp_CrystalIDs2View[f_Crystal1][f_Crystal2];
}
// ==========================================================================================================================================
// void Sino2Crystal()
//   --> This function computes the crystal coordinates (inside a ring) for given element and view.
// ==========================================================================================================================================
void oTableBiograph2D::ElemView2CrystalIDs(int f_Elem, int f_View, int* f_Crystal1, int* f_Crystal2)
{
  // Compute coordinates
  int opti = m_NbElem/4;
  int det1_c = f_View;
  int det2_c = f_View + (m_NbCrystalsPerRingWithGaps/2);
  *f_Crystal1 = det1_c + f_Elem/2     - opti;
  *f_Crystal2 = det2_c - (f_Elem+1)/2 + opti;

  // Apply modulos
  if (*f_Crystal1 < 0) *f_Crystal1 += m_NbCrystalsPerRingWithGaps;
  else if (*f_Crystal1 >= m_NbCrystalsPerRingWithGaps) *f_Crystal1 -= m_NbCrystalsPerRingWithGaps;
  if (*f_Crystal2 < 0) *f_Crystal2 += m_NbCrystalsPerRingWithGaps;
  else if (*f_Crystal2 >= m_NbCrystalsPerRingWithGaps) *f_Crystal2 -= m_NbCrystalsPerRingWithGaps;
}
// ==========================================================================================================================================
// int GetLORsFromSinoBin()
//   --> This function computes the number of LORs contributing to a sinogram bin, and associated crystal indices.
// ==========================================================================================================================================
int oTableBiograph2D::GetLORsFromSinoBin( int f_BinSino, int f_BinView, int f_BinElem, unsigned int* fp_ID1, unsigned int* fp_ID2 )
{
  // Get xpos
  int xpos1 = mp_ElemView2CrystalID1[f_BinElem][f_BinView][0];
  int xpos2 = mp_ElemView2CrystalID2[f_BinElem][f_BinView][0];

  // Check if we are in a gap, we skip it but decrement the number of events stored
  if ((xpos1+1)%(mp_Scanner->GetNbTransCrystalsInHead()+1)==0) return 0;
  if ((xpos2+1)%(mp_Scanner->GetNbTransCrystalsInHead()+1)==0) return 0;

  // Compute the crystal index on ring without gaps
  xpos1 -= xpos1/(mp_Scanner->GetNbTransCrystalsInHead()+1);
  xpos2 -= xpos2/(mp_Scanner->GetNbTransCrystalsInHead()+1);

  // Compute IDs and affect
  fp_ID1[0] = xpos1;
  fp_ID2[0] = xpos2;

  // Return number of LORs
  return 1;
}

