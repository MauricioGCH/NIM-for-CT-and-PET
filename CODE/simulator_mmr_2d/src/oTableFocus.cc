#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <ctime>
#include <fstream>
#include "oTableFocus.hh"
#include "oScanner.hh"
#include "oOutputManager.hh"
using namespace std;

// ==========================================================================================================================================
// Constructors
// ==========================================================================================================================================
oTableFocus::oTableFocus( int f_NbElem, int f_NbView, int f_NbSino, int f_MaxRingDiff, int f_Span, oScanner* fp_Scanner, int f_Verbose )
{
  // Verbose
  m_Verbose = f_Verbose;
  if (m_Verbose>=1) LogCout ("oTableFocus::Constructor() -> Initialize Focus conversion table" << endl);

  // Sinogram dimensions
  m_NbSino = f_NbSino;
  m_NbView = f_NbView;
  m_NbElem = f_NbElem;
  m_Span = f_Span;
  m_MaxRingDiff = f_MaxRingDiff;
  mp_Scanner = fp_Scanner;

  // Check dimensions validity and calculate mashing factor
  if (m_NbElem!=mp_Scanner->GetNbElem())
  {
    LogCerr ("***** oTableFocus::Constructor() -> Number of elements read in sinogram is different from the normal number of the scanner !" << endl);
    exit(1);
  }
  if (mp_Scanner->GetNbView()%m_NbView!=0)
  {
    LogCerr ("***** oTableFocus::Constructor() -> Number of views read in sinogram is uncompatible with the scanner !" << endl);
    exit(1);
  }
  m_Mash = mp_Scanner->GetNbView()/m_NbView;

  // Calculate the number of sinograms given the span and maximum ring difference
  if (CalculateNbSino())
  {
    LogCerr ("***** oTableFocus::Constructor() -> Problem while calculating the number of sinograms !" << endl);
    exit(1);
  }
  // Check consistency or affect
  if (m_NbSino<1)
  {
    m_NbSino = m_NbCalculatedSino;
  }
  else if (m_NbSino!=m_NbCalculatedSino)
  {
    LogCerr ("***** oTableFocus::Constructor() -> Theoretical number of sinograms (" << m_NbCalculatedSino << ") differs from the given one (" << m_NbSino << ") !" << endl);
    exit(1);
  }
}
oTableFocus::oTableFocus( int f_Mash, int f_MaxRingDiff, int f_Span, oScanner* f_Scanner, int f_Verbose )
{
  // Verbose
  m_Verbose = f_Verbose;
  if (m_Verbose>=1) LogCout ("oTableFocus::Constructor() -> Initialize Focus conversion table" << endl);

  // Sinogram characteristics
  m_Mash = f_Mash;
  m_Span = f_Span;
  m_MaxRingDiff = f_MaxRingDiff;
  mp_Scanner = f_Scanner;
  
  // Calculate the number of sinograms given the span and maximum ring difference
  if (CalculateNbSino())
  {
    LogCerr ("***** oTableFocus::Constructor() -> Problem while calculating the number of sinograms !" << endl);
    exit(1);
  }
  m_NbSino = m_NbCalculatedSino;

  // Affect sinogram dimensions
  m_NbView = mp_Scanner->GetNbView()/m_Mash;
  m_NbElem = mp_Scanner->GetNbElem();
}
// ==========================================================================================================================================
// Destructor
// ==========================================================================================================================================
oTableFocus::~oTableFocus()
{
  if (m_Verbose>=1) LogCout ("oTableFocus::Destructor() -> Destroy conversion table" << endl);

  if (mp_ListRingPairsBySinoIndex && mp_NbRingPairsBySinoIndex)
  {
    for (int s=0; s<m_NbSino; s++)
    {
      if (mp_ListRingPairsBySinoIndex[s])
      {
        for (int r=0; r<mp_NbRingPairsBySinoIndex[s]; r++) if (mp_ListRingPairsBySinoIndex[s][r]) free(mp_ListRingPairsBySinoIndex[s][r]);
        free(mp_ListRingPairsBySinoIndex[s]);
      }
    }
    free(mp_NbRingPairsBySinoIndex);
  }
  if (mp_RingPair2SinogramIndex)
  {
    for (int r=0; r<mp_Scanner->GetNbTotalAxialCrystals(); r++) if (mp_RingPair2SinogramIndex[r]) free(mp_RingPair2SinogramIndex[r]);
    free(mp_RingPair2SinogramIndex);
  }
}
// ==========================================================================================================================================
// int ComputeSpanTables()
//   --> This function compute the span conversion tables and returns 0 upon success and other value otherwise.
// ==========================================================================================================================================
int oTableFocus::ComputeSpanTables()
{
  if (m_Verbose>=1) LogCout ("oTableFocus::ComputeSpanTables() -> Compute the span conversion tables" << endl);

  // Get the number of rings
  int nb_rings = mp_Scanner->GetNbTotalAxialCrystals();

  // Allocations
  mp_NbRingPairsBySinoIndex = (unsigned short int*)calloc(m_NbSino,sizeof(unsigned short int));
  mp_ListRingPairsBySinoIndex = (unsigned int***)malloc(m_NbSino*sizeof(unsigned int**));
  for (int s=0; s<m_NbSino; s++) mp_ListRingPairsBySinoIndex[s] = (unsigned int**)malloc(1*sizeof(unsigned int*));
  mp_RingPair2SinogramIndex = (int**)malloc(nb_rings*sizeof(int*));
  for (int r1=0; r1<nb_rings; r1++)
  {
    mp_RingPair2SinogramIndex[r1] = (int*)malloc(nb_rings*sizeof(int));
    for (int r2=0; r2<nb_rings; r2++) mp_RingPair2SinogramIndex[r1][r2] = -1;
  }

  // Check if the max ring difference is consistent with the span (and calculate number of segments at the same time)
  int check_span = m_MaxRingDiff;
  bool correct = false;
  m_NbAbsoluteSegments = 0;
  while (check_span >= 0)
  {
    if (check_span==m_Span/2)
    {
      correct = true;
      break;
    }
    check_span -= m_Span;
    m_NbAbsoluteSegments++;
  }
/* Apparemment meme si ca coince ici ça peut passer
  if (!correct)
  {
    LogCerr ("***** oTableFocus::ComputeSpanTables() -> Inconsistence between the span (" << m_Span << ") and the maximum ring difference (" << m_MaxRingDiff << ") !" << endl);
    return 1;
  }
*/
  m_NbSegments = 2*m_NbAbsoluteSegments+1;
  if (m_Verbose>=1) LogCout ("  --> Number of segments: " << m_NbSegments << endl);

  // Main loop on segments and increasing sino
  int segment_number = 0;
  int sino_bin = 0;
  for (int s=0; s<m_NbSegments; s++)
  {
    // Calculate the min and max allowed ring differences for that segment
    int min_allowed_ring_difference = -m_Span/2 + segment_number*m_Span;
    int max_allowed_ring_difference = m_Span/2 + segment_number*m_Span;

    // Do a loop on all possible ring indices sum
    for (int sum=0; sum<mp_Scanner->GetNbTotalAxialCrystals()*2; sum++)
    {
      // Loop on all ring pairs
      for (int r1=0; r1<nb_rings; r1++) for (int r2=0; r2<nb_rings; r2++)
      {
        // First test if this ring pairs is in the good segment
        int ring_diff = r1-r2;
        if (ring_diff<min_allowed_ring_difference || ring_diff>max_allowed_ring_difference) continue;

        // Then test if the ring sum match the current seeked sum
        int ring_sum = r1+r2;
        if (ring_sum!=sum) continue;

        // Increment the number of ring pairs for this sinogram
        mp_NbRingPairsBySinoIndex[sino_bin]++;
        // Realloc the list of ring pairs for this sinogram
        mp_ListRingPairsBySinoIndex[sino_bin] = (unsigned int**)realloc(mp_ListRingPairsBySinoIndex[sino_bin], mp_NbRingPairsBySinoIndex[sino_bin]*sizeof(unsigned int*));
        // Alloc and affect the new ring pair
        mp_ListRingPairsBySinoIndex[sino_bin][mp_NbRingPairsBySinoIndex[sino_bin]-1] = (unsigned int*)malloc(2*sizeof(unsigned int));
        mp_ListRingPairsBySinoIndex[sino_bin][mp_NbRingPairsBySinoIndex[sino_bin]-1][0] = r1;
        mp_ListRingPairsBySinoIndex[sino_bin][mp_NbRingPairsBySinoIndex[sino_bin]-1][1] = r2;
        // Set the sinogram indices for this ring pair in the table
        mp_RingPair2SinogramIndex[r1][r2] = sino_bin;
      }

      // Increment the sinogram index if some ring pairs were found
      if (mp_NbRingPairsBySinoIndex[sino_bin]>0)
      {
        if (m_Verbose>=3) LogCout ("        --> Sino index: " << sino_bin << " | Nb ring pairs: " << mp_NbRingPairsBySinoIndex[sino_bin] << endl);
        sino_bin++;
      }
    }

    // Update the segment number
    int plus_or_minus = (1 - (s%2)*2);
    segment_number = segment_number + plus_or_minus*(s+1);
  }

  // Ending
  return 0;
}
// ==========================================================================================================================================
// Function: CalculateNbSino
//   --> This function calculates the number of sinograms given the span and maximum ring difference
// ==========================================================================================================================================
int oTableFocus::CalculateNbSino()
{
  // Get the number of rings
  int nb_rings = mp_Scanner->GetNbTotalAxialCrystals();

  // Check if the max ring difference is consistent with the span (and calculate number of segments at the same time)
  int check_span = m_MaxRingDiff;
  bool correct = false;
  m_NbAbsoluteSegments = 0;
  while (check_span >= 0)
  {
    if (check_span==m_Span/2)
    {
      correct = true;
      break;
    }
    check_span -= m_Span;
    m_NbAbsoluteSegments++;
  }
/* Apparemment meme si ca coince ici ça peut passer
  if (!correct)
  {
    LogCerr ("***** oTableFocus::CalculateNbSino() -> Inconsistence between the span (" << m_Span << ") and the maximum ring difference (" << m_MaxRingDiff << ") !" << endl);
    return 1;
  }
*/
  m_NbSegments = 2*m_NbAbsoluteSegments+1;

  // Main loop on segments and increasing sino
  int segment_number = 0;
  int sino_bin = 0;
  int nb_ring_pairs_by_sino = 0;
  for (int s=0; s<m_NbSegments; s++)
  {
    // Calculate the min and max allowed ring differences for that segment
    int min_allowed_ring_difference = -m_Span/2 + segment_number*m_Span;
    int max_allowed_ring_difference = m_Span/2 + segment_number*m_Span;

    // Do a loop on all possible ring indices sum
    for (int sum=0; sum<mp_Scanner->GetNbTotalAxialCrystals()*2; sum++)
    {
      // Loop on all ring pairs
      for (int r1=0; r1<nb_rings; r1++) for (int r2=0; r2<nb_rings; r2++)
      {
        // First test if this ring pairs is in the good segment
        int ring_diff = r1-r2;
        if (ring_diff<min_allowed_ring_difference || ring_diff>max_allowed_ring_difference) continue;

        // Then test if the ring sum match the current seeked sum
        int ring_sum = r1+r2;
        if (ring_sum!=sum) continue;

        // Increment the number of ring pairs for this sinogram
        nb_ring_pairs_by_sino++;
      }

      // Increment the sinogram index if some ring pairs were found
      if (nb_ring_pairs_by_sino>0)
      {
        sino_bin++;
        nb_ring_pairs_by_sino = 0;
      }
    }

    // Update the segment number
    int plus_or_minus = (1 - (s%2)*2);
    segment_number = segment_number + plus_or_minus*(s+1);
  }

  // Affect
  m_NbCalculatedSino = sino_bin;

  // End
  return 0;
}
// ==========================================================================================================================================
// int ComputeRingTables()
//   --> This function computes the ring conversion tables (crystals to elem,view) and returns 0 upon success and other value otherwise.
// ==========================================================================================================================================
int oTableFocus::ComputeRingTables()
{
  if (m_Verbose>=1) LogCout ("oTableFocus::ComputeRingTables() -> Compute the ring conversion tables" << endl);

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
  int nb_total_trans_crystals = mp_Scanner->GetNbTotalTransCrystals();
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
unsigned short int oTableFocus::GetNbRingPairsBySinoIndex(int f_Sino)
{
  return mp_NbRingPairsBySinoIndex[f_Sino];
}
void oTableFocus::GetRingPairBySinoIndex(int f_Sino, int f_RingPair, unsigned int* f_Ring1, unsigned int* f_Ring2)
{
  *f_Ring1 = mp_ListRingPairsBySinoIndex[f_Sino][f_RingPair][0];
  *f_Ring2 = mp_ListRingPairsBySinoIndex[f_Sino][f_RingPair][1];
}
void oTableFocus::GetRingPairsBySinoIndex(int f_Sino, unsigned int* f_Ring1, unsigned int* f_Ring2)
{
  for (int r=0; r<mp_NbRingPairsBySinoIndex[f_Sino]; r++)
  {
    f_Ring1[r] = mp_ListRingPairsBySinoIndex[f_Sino][r][0];
    f_Ring2[r] = mp_ListRingPairsBySinoIndex[f_Sino][r][1];
  }
}
int oTableFocus::GetSinoIndexFromRingPair(unsigned int f_Ring1, unsigned int f_Ring2)
{
  return mp_RingPair2SinogramIndex[f_Ring1][f_Ring2];
}
void oTableFocus::GetCrystalIDsFromElemView(int f_Elem, int f_View, int* f_Crystal1, int* f_Crystal2)
{
  for (int m=0; m<m_Mash; m++)
  {
    f_Crystal1[m] = mp_ElemView2CrystalID1[f_Elem][f_View][m];
    f_Crystal2[m] = mp_ElemView2CrystalID2[f_Elem][f_View][m];
  }
}
void oTableFocus::GetCrystalIDFromElemView(int f_Index, int f_Elem, int f_View, int* f_Crystal1, int* f_Crystal2)
{
  *f_Crystal1 = mp_ElemView2CrystalID1[f_Elem][f_View][f_Index];
  *f_Crystal2 = mp_ElemView2CrystalID2[f_Elem][f_View][f_Index];
}
void oTableFocus::GetElemViewFromCrystalIDs(int* f_Elem, int* f_View, int f_Crystal1, int f_Crystal2)
{
  *f_Elem = mp_CrystalIDs2Elem[f_Crystal1][f_Crystal2];
  *f_View = mp_CrystalIDs2View[f_Crystal1][f_Crystal2];
}
// ==========================================================================================================================================
// void ElemView2CrystalIDs()
//   --> This function computes the crystal coordinates (inside a ring) for given element and view.
// ==========================================================================================================================================
void oTableFocus::ElemView2CrystalIDs(int f_Elem, int f_View, int* f_Crystal1, int* f_Crystal2)
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
int oTableFocus::GetLORsFromSinoBin( int f_BinSino, int f_BinView, int f_BinElem, unsigned int* fp_ID1, unsigned int* fp_ID2 )
{
  // The number of LORs
  int nb_lors = 0;

  // Loop on span
  for (int r=0; r<mp_NbRingPairsBySinoIndex[f_BinSino]; r++)
  {
    // Get rings
    int ring1 = mp_ListRingPairsBySinoIndex[f_BinSino][r][0];
    int ring2 = mp_ListRingPairsBySinoIndex[f_BinSino][r][1];

    // Loop on mash
    for (int m=0; m<m_Mash; m++)
    {
      // Get xpos
      int xpos1 = mp_ElemView2CrystalID1[f_BinElem][f_BinView][m];
      int xpos2 = mp_ElemView2CrystalID2[f_BinElem][f_BinView][m];

      // Compute IDs and affect
      fp_ID1[nb_lors] = ring1 * mp_Scanner->GetNbTotalTransCrystals() + xpos1;
      fp_ID2[nb_lors] = ring2 * mp_Scanner->GetNbTotalTransCrystals() + xpos2;
      // Increment number of LORs
      nb_lors++;
    }
  }

  // Return number of LORs
  return nb_lors;
}

