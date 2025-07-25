#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <ctime>
#include <fstream>
#include "oTableCastor.hh"
#include "oOutputManager.hh"
#include "oScanner.hh"
using namespace std;

// ==========================================================================================================================================
// Constructor
// ==========================================================================================================================================
oTableCastor::oTableCastor( oScanner* fp_Scanner, int f_Verbose )
{
  // Verbose
  m_Verbose = f_Verbose;
  if (m_Verbose>=1) LogCout ("oTableCastor::Constructor() -> Initialize biograph conversion table" << endl);

  // Sinogram dimensions and characteristics
  m_NbSino = 109; // Imposed by an eletronic card inside the scanner !
  m_Span = 11;    // Imposed by an eletronic card inside the scanner !
  m_NbElem = 160; // Imposed by an eletronic card inside the scanner !
  m_NbView = 112; // Imposed by an eletronic card inside the scanner !
  m_NbCrystalsPerRingWithGaps = 14*48;
  m_Mash = 3;     // Imposed by an eletronic card inside the scanner !
  m_MaxRingDiff = 5; // Imposed by an eletronic card inside the scanner !

  // Scanner model
  mp_Scanner = fp_Scanner;
}
// ==========================================================================================================================================
// Destructor
// ==========================================================================================================================================
oTableCastor::~oTableCastor()
{
  if (m_Verbose>=1) LogCout ("oTableCastor::Destructor() -> Destroy conversion table" << endl);

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
}
// ==========================================================================================================================================
// int ReadTable()
//   --> This function read the conversion table from the given file and returns 0 upon success and other value otherwise.
// ==========================================================================================================================================
int oTableCastor::ReadSpanTable(string f_ScannerTable)
{
  // Affect table file name
  m_TableFileName = f_ScannerTable;

  // Get table file name if no table name is given
  if (m_TableFileName=="")
  {
    char* tmp_string = getenv("RECON_TABLES");
    if (tmp_string==NULL)
    {
      LogCerr ("***** oTableCastor::Constructor() -> Environment variable 'RECON_TABLES' is not set !" << endl);
      LogCerr ("                                       Please set it to point to the appropriate tables." << endl);
      exit(1);
    }
    string table_path = (string)tmp_string;
    m_TableFileName = table_path+"/biograph_4rings.txt";
  }

  if (m_Verbose>=1) LogCout ("oTableCastor::ReadSpanTable() -> Read table from file '" << m_TableFileName << "'" << endl);

  // Allocation
  mp_NbRingPairsBySinoIndex = (unsigned short int*)calloc(m_NbSino,sizeof(unsigned short int));
  mp_ListRingPairsBySinoIndex = (unsigned int***)malloc(m_NbSino*sizeof(unsigned int**));
  for (int s=0; s<m_NbSino; s++) mp_ListRingPairsBySinoIndex[s] = (unsigned int**)malloc(1*sizeof(unsigned int*));
  int nb_rings = mp_Scanner->GetNbTotalAxialCrystalsWithGaps();
  mp_RingPair2SinogramIndex = (int**)malloc(nb_rings*sizeof(int*));
  for (int r=0; r<nb_rings; r++) mp_RingPair2SinogramIndex[r] = (int*)malloc(nb_rings*sizeof(int));

  // Open the file
  ifstream ftab(m_TableFileName.c_str());
  if (!ftab)
  {
    LogCerr ("***** oTableCastor::ReadSpanTable() -> Input table file '" << m_TableFileName << "' is missing or corrupted !" << endl);
    return 1;
  }

  // Read it now
  string buffer;
  ftab >> buffer;
  for (int s=0; s<m_NbSino; s++)
  {
    // First read a new buffer to get a ring pair in it
    ftab >> buffer;

    // Then do a while loop as the ':' character means we get a valid ring pair
    while (buffer.find(":")==string::npos && buffer.find("<")!=string::npos)
    {
      // Increment the number of ring pairs for this sinogram
      mp_NbRingPairsBySinoIndex[s]++;
      // Realloc the list of ring pairs for this sinogram
      mp_ListRingPairsBySinoIndex[s] = (unsigned int**)realloc(mp_ListRingPairsBySinoIndex[s], mp_NbRingPairsBySinoIndex[s]*sizeof(unsigned int*));
      // Decode the ring indices
      unsigned int ring1, ring2;
      if (sscanf(buffer.c_str(),"<%u,%u>",&ring1,&ring2) != 2)
      {
        LogCerr ("***** oTableCastor::ReadSpanTable() -> Error decoding a ring pair \"" << buffer << "\" !" << endl);
        ftab.close();
        return 1;
      }
      // Alloc and affect the new ring pair
      mp_ListRingPairsBySinoIndex[s][mp_NbRingPairsBySinoIndex[s]-1] = (unsigned int*)malloc(2*sizeof(unsigned int));
      mp_ListRingPairsBySinoIndex[s][mp_NbRingPairsBySinoIndex[s]-1][0] = ring1;
      mp_ListRingPairsBySinoIndex[s][mp_NbRingPairsBySinoIndex[s]-1][1] = ring2;
      // Set the sinogram indices for this ring pair in the table
      mp_RingPair2SinogramIndex[ring1][ring2] = s;
      // Read a new buffer
      ftab >> buffer;
    }
  }

  // Close the file
  ftab.close();

  // Ending
  return 0;
}
// ==========================================================================================================================================
// int ComputeRingTables()
//   --> This function computes the ring conversion tables (crystals to elem,view) and returns 0 upon success and other value otherwise.
// ==========================================================================================================================================
int oTableCastor::ComputeRingTables()
{
  if (m_Verbose>=1) LogCout ("oTableCastor::ComputeRingTables() -> Compute the ring conversion tables" << endl);

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
//  mp_NbCrystalPairsByElemView = (int**)malloc(m_NbElem*sizeof(int*));
//  for (int e=0; e<m_nbElem; e++) mp_NbCrystalPairsByElemView[e] = (int*)malloc(m_NbView*sizeof(int));

  // Loop on all elem,view indices
  int* crystals1 = (int*)malloc(m_Mash*sizeof(int));
  int* crystals2 = (int*)malloc(m_Mash*sizeof(int));
  for (int e=0; e<m_NbElem; e++) for (int v=0; v<m_NbView; v++)
  {
//    mp_NbCrystalPairsByElemView[e][v] = ElemView2CrystalIDs(e, v, crystals1, crystals2);
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
unsigned short int oTableCastor::GetNbRingPairsBySinoIndex(int f_Sino)
{
  return mp_NbRingPairsBySinoIndex[f_Sino];
}
void oTableCastor::GetRingPairBySinoIndex(int f_Sino, int f_RingPair, unsigned int* f_Ring1, unsigned int* f_Ring2)
{
  *f_Ring1 = mp_ListRingPairsBySinoIndex[f_Sino][f_RingPair][0];
  *f_Ring2 = mp_ListRingPairsBySinoIndex[f_Sino][f_RingPair][1];
}
void oTableCastor::GetRingPairsBySinoIndex(int f_Sino, unsigned int* f_Ring1, unsigned int* f_Ring2)
{
  for (int r=0; r<mp_NbRingPairsBySinoIndex[f_Sino]; r++)
  {
    f_Ring1[r] = mp_ListRingPairsBySinoIndex[f_Sino][r][0];
    f_Ring2[r] = mp_ListRingPairsBySinoIndex[f_Sino][r][1];
  }
}
int oTableCastor::GetSinoIndexFromRingPair(unsigned int f_Ring1, unsigned int f_Ring2)
{
  return mp_RingPair2SinogramIndex[f_Ring1][f_Ring2];
}
void oTableCastor::GetCrystalIDsFromElemView(int f_Elem, int f_View, int* f_Crystal1, int* f_Crystal2)
{
  for (int m=0; m<m_Mash; m++)
  {
    f_Crystal1[m] = mp_ElemView2CrystalID1[f_Elem][f_View][m];
    f_Crystal2[m] = mp_ElemView2CrystalID2[f_Elem][f_View][m];
  }
/*
  *f_Crystal1 = mp_ElemView2CrystalID1[f_Elem][f_View][0];
  *f_Crystal2 = mp_ElemView2CrystalID2[f_Elem][f_View][0];
*/
}
void oTableCastor::GetElemViewFromCrystalIDs(int* f_Elem, int* f_View, int f_Crystal1, int f_Crystal2)
{
  *f_Elem = mp_CrystalIDs2Elem[f_Crystal1][f_Crystal2];
  *f_View = mp_CrystalIDs2View[f_Crystal1][f_Crystal2];
}
// ==========================================================================================================================================
// void Sino2Crystal()
//   --> This function computes the crystal coordinates (inside a ring) for given element and view.
// ==========================================================================================================================================
void oTableCastor::ElemView2CrystalIDs(int f_Elem, int f_View, int* f_Crystal1, int* f_Crystal2)
{
  // The number of crystal pairs
//  int nb_crystal_pairs = 0;
  // Loop on mash
  for (int m=0; m<m_Mash; m++)
  {
    // Compute coordinates
    int opti = m_NbElem/4;
    int det1_c = f_View*m_Mash + m;
    int det2_c = f_View*m_Mash + m + (mp_Scanner->GetNbTotalTransCrystalsWithGaps()/2);
    int crystal1 = det1_c + f_Elem/2     - opti;
    int crystal2 = det2_c - (f_Elem+1)/2 + opti;
    // Apply modulos
    if (crystal1 < 0) crystal1 += mp_Scanner->GetNbTotalTransCrystalsWithGaps();
    else if (crystal1 >= mp_Scanner->GetNbTotalTransCrystalsWithGaps()) crystal1 -= mp_Scanner->GetNbTotalTransCrystalsWithGaps();
    if (crystal2 < 0) crystal2 += mp_Scanner->GetNbTotalTransCrystalsWithGaps();
    else if (crystal2 >= mp_Scanner->GetNbTotalTransCrystalsWithGaps()) crystal2 -= mp_Scanner->GetNbTotalTransCrystalsWithGaps();
/*
    // Check if we are in a gap, we skip it
    if ((crystals1[t][0]+1)%(nb_trans_crystals_in_head+1)==0) continue;
    if ((crystals2[t][0]+1)%(nb_trans_crystals_in_head+1)==0) continue;
    // Compute the crystal indices on ring without gaps
        crystals1[t][0] -= crystals1[t][0]/(nb_trans_crystals_in_head+1);
        crystals2[t][0] -= crystals2[t][0]/(nb_trans_crystals_in_head+1);
*/
    f_Crystal1[m] = crystal1;
    f_Crystal2[m] = crystal2;
  }


/*
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
*/
}
// ==========================================================================================================================================
// int GetLORsFromSinoBin()
//   --> This function computes the number of LORs contributing to a sinogram bin, and associated crystal indices.
// ==========================================================================================================================================
int oTableCastor::GetLORsFromSinoBin( int f_BinSino, int f_BinView, int f_BinElem, unsigned int* fp_ID1, unsigned int* fp_ID2 )
{
  // The number of LORs
  int nb_lors = 0;

  // Loop on mashing
  for (int m=0; m<m_Mash; m++)
  {
    // Get xpos of the crystal pair
    int xpos1 = mp_ElemView2CrystalID1[f_BinElem][f_BinView][m];
    int xpos2 = mp_ElemView2CrystalID2[f_BinElem][f_BinView][m];

    // Check if we are in a gap, we skip it
    if ((xpos1+1)%(mp_Scanner->GetNbTransCrystalsInHead()+1)==0) continue;
    if ((xpos2+1)%(mp_Scanner->GetNbTransCrystalsInHead()+1)==0) continue;

    // Compute the crystal index on ring without gaps
    xpos1 -= xpos1/(mp_Scanner->GetNbTransCrystalsInHead()+1);
    xpos2 -= xpos2/(mp_Scanner->GetNbTransCrystalsInHead()+1);

    // Loop on span
    for (int r=0; r<mp_NbRingPairsBySinoIndex[f_BinSino]; r++)
    {
      // Get rings
      int ring1 = mp_ListRingPairsBySinoIndex[f_BinSino][r][0];
      int ring2 = mp_ListRingPairsBySinoIndex[f_BinSino][r][1];

      // Compute the ring indices without gaps (no need to check wether we are in a gap because gaps are already removed from the span table)
      ring1 -= ring1/(mp_Scanner->GetNbAxialCrystals()+1);
      ring2 -= ring2/(mp_Scanner->GetNbAxialCrystals()+1);

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

