#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <ctime>
#include <fstream>
#include "oTableMMR2D.hh"
#include "oScanner.hh"
#include "oOutputManager.hh"
using namespace std;

// ==========================================================================================================================================
// Constructors
// ==========================================================================================================================================
oTableMMR2D::oTableMMR2D( int f_NbElem, int f_NbView, int f_NbSino, int f_MaxRingDiff, int f_Span, oScanner* fp_Scanner, int f_Verbose )
{
  // Verbose
  m_Verbose = f_Verbose;
  if (m_Verbose>=1) LogCout ("oTableMMR2D::Constructor() -> Initialize MMR 2D conversion table" << endl);

  // Sinogram dimensions
  m_NbSino = 1;
  m_NbView = f_NbView;
  m_NbElem = f_NbElem;
  m_Span = 1;
  m_MaxRingDiff = 0;
  mp_Scanner = fp_Scanner;

  // Check dimensions validity and calculate mashing factor
  if (m_NbElem!=mp_Scanner->GetNbElem())
  {
    LogCerr ("***** oTableMMR2D::Constructor() -> Number of elements read in sinogram is different from the normal number of the scanner !" << endl);
    exit(1);
  }
  if (mp_Scanner->GetNbView()%m_NbView!=0)
  {
    LogCerr ("***** oTableMMR2D::Constructor() -> Number of views read in sinogram is uncompatible with the scanner !" << endl);
    exit(1);
  }
  m_Mash = mp_Scanner->GetNbView()/m_NbView;

  m_NbCalculatedSino = 1;
}
oTableMMR2D::oTableMMR2D( oScanner* f_Scanner, int f_Verbose )
{
  // Verbose
  m_Verbose = f_Verbose;
  if (m_Verbose>=1) LogCout ("oTableMMR2D::Constructor() -> Initialize MMR 2D conversion table" << endl);

  // Sinogram characteristics
  m_Mash = 1;
  m_Span = 1;
  m_MaxRingDiff = 0;
  mp_Scanner = f_Scanner;
  
  m_NbCalculatedSino = 1;

  // Affect sinogram dimensions
  m_NbView = mp_Scanner->GetNbView()/m_Mash;
  m_NbElem = mp_Scanner->GetNbElem();
}
// ==========================================================================================================================================
// Destructor
// ==========================================================================================================================================
oTableMMR2D::~oTableMMR2D()
{
  if (m_Verbose>=1) LogCout ("oTableMMR2D::Destructor() -> Destroy conversion table" << endl);
}
// ==========================================================================================================================================
// int ComputeRingTables()
//   --> This function computes the ring conversion tables (crystals to elem,view) and returns 0 upon success and other value otherwise.
// ==========================================================================================================================================
int oTableMMR2D::ComputeRingTables()
{
  if (m_Verbose>=1) LogCout ("oTableMMR2D::ComputeRingTables() -> Compute the ring conversion tables" << endl);

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
void oTableMMR2D::GetCrystalIDsFromElemView(int f_Elem, int f_View, int* f_Crystal1, int* f_Crystal2)
{
  for (int m=0; m<m_Mash; m++)
  {
    f_Crystal1[m] = mp_ElemView2CrystalID1[f_Elem][f_View][m];
    f_Crystal2[m] = mp_ElemView2CrystalID2[f_Elem][f_View][m];
  }
}
void oTableMMR2D::GetCrystalIDFromElemView(int f_Index, int f_Elem, int f_View, int* f_Crystal1, int* f_Crystal2)
{
  *f_Crystal1 = mp_ElemView2CrystalID1[f_Elem][f_View][f_Index];
  *f_Crystal2 = mp_ElemView2CrystalID2[f_Elem][f_View][f_Index];
}
void oTableMMR2D::GetElemViewFromCrystalIDs(int* f_Elem, int* f_View, int f_Crystal1, int f_Crystal2)
{
  *f_Elem = mp_CrystalIDs2Elem[f_Crystal1][f_Crystal2];
  *f_View = mp_CrystalIDs2View[f_Crystal1][f_Crystal2];
}
// ==========================================================================================================================================
// void ElemView2CrystalIDs()
//   --> This function computes the crystal coordinates (inside a ring) for given element and view.
// ==========================================================================================================================================
void oTableMMR2D::ElemView2CrystalIDs(int f_Elem, int f_View, int* f_Crystal1, int* f_Crystal2)
{
  // Loop on mash
  for (int m=0; m<m_Mash; m++)
  {
    // Compute coordinates
    int opti = m_NbElem/4;
    int det1_c = f_View*m_Mash + m;
    int det2_c = f_View*m_Mash + m + (mp_Scanner->GetNbTotalTransCrystalsWithGaps()/2);
    f_Crystal1[m] = det1_c + f_Elem/2     - opti;
    f_Crystal2[m] = det2_c - (f_Elem+1)/2 + opti;
    // Apply modulos
    if (f_Crystal1[m] < 0) f_Crystal1[m] += mp_Scanner->GetNbTotalTransCrystalsWithGaps();
    else if (f_Crystal1[m] >= mp_Scanner->GetNbTotalTransCrystalsWithGaps()) f_Crystal1[m] -= mp_Scanner->GetNbTotalTransCrystalsWithGaps();
    if (f_Crystal2[m] < 0) f_Crystal2[m] += mp_Scanner->GetNbTotalTransCrystalsWithGaps();
    else if (f_Crystal2[m] >= mp_Scanner->GetNbTotalTransCrystalsWithGaps()) f_Crystal2[m] -= mp_Scanner->GetNbTotalTransCrystalsWithGaps();
  }
}
// ==========================================================================================================================================
// int GetLORsFromSinoBin()
//   --> This function computes the number of LORs contributing to a sinogram bin, and associated crystal indices.
// ==========================================================================================================================================
int oTableMMR2D::GetLORsFromSinoBin( int f_BinSino, int f_BinView, int f_BinElem, unsigned int* fp_ID1, unsigned int* fp_ID2 )
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

