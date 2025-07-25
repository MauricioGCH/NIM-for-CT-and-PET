#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <ctime>
#include <fstream>
#include "oTableHRRT.hh"
#include "oOutputManager.hh"
#include "oScanner.hh"
using namespace std;

// ==========================================================================================================================================
// Constructor
// ==========================================================================================================================================
oTableHRRT::oTableHRRT( oScanner* fp_Scanner, int f_Verbose )
{
  // Verbose
  m_Verbose = f_Verbose;
  if (m_Verbose>=1) LogCout ("oTableHRRT::Constructor() -> Initialize hrrt conversion table" << endl);

  // Scanner model
  mp_Scanner = fp_Scanner;

  // Fixed table for possible head pairs
  mp_HeadPairs = (int**)malloc(8*sizeof(int*));
  for (int h1=0; h1<8; h1++)
  {
    mp_HeadPairs[h1] = (int*)malloc(8*sizeof(int));
    for (int h2=0; h2<8; h2++) mp_HeadPairs[h1][h2] = -1;
  }
  mp_HeadPairs[0][2] = 1;
  mp_HeadPairs[0][3] = 2;
  mp_HeadPairs[0][4] = 3;
  mp_HeadPairs[0][5] = 4;
  mp_HeadPairs[0][6] = 5;
  mp_HeadPairs[1][3] = 6;
  mp_HeadPairs[1][4] = 7;
  mp_HeadPairs[1][5] = 8;
  mp_HeadPairs[1][6] = 9;
  mp_HeadPairs[1][7] = 10;
  mp_HeadPairs[2][4] = 11;
  mp_HeadPairs[2][5] = 12;
  mp_HeadPairs[2][6] = 13;
  mp_HeadPairs[2][7] = 14;
  mp_HeadPairs[3][5] = 15;
  mp_HeadPairs[3][6] = 16;
  mp_HeadPairs[3][7] = 17;
  mp_HeadPairs[4][6] = 18;
  mp_HeadPairs[4][7] = 19;
  mp_HeadPairs[5][7] = 20;

  // Same fixed table but reversed
  mp_HeadPairsReversed = (int**)malloc(21*sizeof(int*));
  for (int mp=0; mp<21; mp++) mp_HeadPairsReversed[mp] = (int*)calloc(2,sizeof(int));
  mp_HeadPairsReversed[0][0]  = -1;  mp_HeadPairsReversed[0][1]  = -1;
  mp_HeadPairsReversed[1][0]  = 0;   mp_HeadPairsReversed[1][1]  = 2;
  mp_HeadPairsReversed[2][0]  = 0;   mp_HeadPairsReversed[2][1]  = 3;
  mp_HeadPairsReversed[3][0]  = 0;   mp_HeadPairsReversed[3][1]  = 4;
  mp_HeadPairsReversed[4][0]  = 0;   mp_HeadPairsReversed[4][1]  = 5;
  mp_HeadPairsReversed[5][0]  = 0;   mp_HeadPairsReversed[5][1]  = 6;
  mp_HeadPairsReversed[6][0]  = 1;   mp_HeadPairsReversed[6][1]  = 3;
  mp_HeadPairsReversed[7][0]  = 1;   mp_HeadPairsReversed[7][1]  = 4;
  mp_HeadPairsReversed[8][0]  = 1;   mp_HeadPairsReversed[8][1]  = 5;
  mp_HeadPairsReversed[9][0]  = 1;   mp_HeadPairsReversed[9][1]  = 6;
  mp_HeadPairsReversed[10][0] = 1;   mp_HeadPairsReversed[10][1] = 7;
  mp_HeadPairsReversed[11][0] = 2;   mp_HeadPairsReversed[11][1] = 4;
  mp_HeadPairsReversed[12][0] = 2;   mp_HeadPairsReversed[12][1] = 5;
  mp_HeadPairsReversed[13][0] = 2;   mp_HeadPairsReversed[13][1] = 6;
  mp_HeadPairsReversed[14][0] = 2;   mp_HeadPairsReversed[14][1] = 7;
  mp_HeadPairsReversed[15][0] = 3;   mp_HeadPairsReversed[15][1] = 5;
  mp_HeadPairsReversed[16][0] = 3;   mp_HeadPairsReversed[16][1] = 6;
  mp_HeadPairsReversed[17][0] = 3;   mp_HeadPairsReversed[17][1] = 7;
  mp_HeadPairsReversed[18][0] = 4;   mp_HeadPairsReversed[18][1] = 6;
  mp_HeadPairsReversed[19][0] = 4;   mp_HeadPairsReversed[19][1] = 7;
  mp_HeadPairsReversed[20][0] = 5;   mp_HeadPairsReversed[20][1] = 7;

  // Fixed table for span3 to span9 conversion
  mp_Span3to9ConversionTable = (int*)calloc(45,sizeof(int));
  mp_Span3to9ConversionTable[0]  = 0;    mp_Span3to9ConversionTable[1]  = 0;    mp_Span3to9ConversionTable[2]  = 0;
  mp_Span3to9ConversionTable[3]  = 1;    mp_Span3to9ConversionTable[4]  = 2;    mp_Span3to9ConversionTable[5]  = 1;
  mp_Span3to9ConversionTable[6]  = 2;    mp_Span3to9ConversionTable[7]  = 1;    mp_Span3to9ConversionTable[8]  = 2;
  mp_Span3to9ConversionTable[9]  = 3;    mp_Span3to9ConversionTable[10] = 4;    mp_Span3to9ConversionTable[11] = 3;
  mp_Span3to9ConversionTable[12] = 4;    mp_Span3to9ConversionTable[13] = 3;    mp_Span3to9ConversionTable[14] = 4;
  mp_Span3to9ConversionTable[15] = 5;    mp_Span3to9ConversionTable[16] = 6;    mp_Span3to9ConversionTable[17] = 5;
  mp_Span3to9ConversionTable[18] = 6;    mp_Span3to9ConversionTable[19] = 5;    mp_Span3to9ConversionTable[20] = 6;
  mp_Span3to9ConversionTable[21] = 7;    mp_Span3to9ConversionTable[22] = 8;    mp_Span3to9ConversionTable[23] = 7;
  mp_Span3to9ConversionTable[24] = 8;    mp_Span3to9ConversionTable[25] = 7;    mp_Span3to9ConversionTable[26] = 8;
  mp_Span3to9ConversionTable[27] = 9;    mp_Span3to9ConversionTable[28] = 10;   mp_Span3to9ConversionTable[29] = 9;
  mp_Span3to9ConversionTable[30] = 10;   mp_Span3to9ConversionTable[31] = 9;    mp_Span3to9ConversionTable[32] = 10;
  mp_Span3to9ConversionTable[33] = 11;   mp_Span3to9ConversionTable[34] = 12;   mp_Span3to9ConversionTable[35] = 11;
  mp_Span3to9ConversionTable[36] = 12;   mp_Span3to9ConversionTable[37] = 11;   mp_Span3to9ConversionTable[38] = 12;
  mp_Span3to9ConversionTable[39] = 13;   mp_Span3to9ConversionTable[40] = 14;   mp_Span3to9ConversionTable[41] = 13;
  mp_Span3to9ConversionTable[42] = 14;   mp_Span3to9ConversionTable[43] = 13;   mp_Span3to9ConversionTable[44] = 14;

  // LUT allocation
  mp_LUTsino = (SOL***)calloc(21,sizeof(SOL**));
  for (int mp=1; mp<=20; mp++)
  {
    mp_LUTsino[mp] = (SOL**)calloc(NXCRYS*NDOIS,sizeof(SOL *));
    for (int ax=0; ax<NXCRYS*NDOIS; ax++) mp_LUTsino[mp][ax] = (SOL*)calloc(NXCRYS*NDOIS,sizeof(SOL));
  }
  mp_LUTzpos  = (float*)calloc(NYCRYS,sizeof(float));
  mp_LUTzpos2 = (float*)calloc(NYCRYS,sizeof(float));

  // Fixed values
  m_MaxRingDiff = 67;
  m_NbElem = 256;
  m_NbView = 288;
  m_NbSino = 2207;
  m_Span = 9;
  m_Mash = 1;

  // Default values
  m_segz0 = NULL;
  m_segzmax = NULL;
  m_segzoffset = NULL;
  m_nsegs = -1;
  m_segplane = NULL;
}
// ==========================================================================================================================================
// Destructor
// ==========================================================================================================================================
oTableHRRT::~oTableHRRT()
{
  if (m_Verbose>=1) LogCout ("oTableHRRT::Destructor() -> Destroy conversion table" << endl);

}
// ==========================================================================================================================================
// int ReadLUTTable()
//   --> This function reads the LUT conversion table from the given file and returns 0 upon success and other value otherwise.
// ==========================================================================================================================================
int oTableHRRT::ReadLUTTable(string f_ScannerTable)
{
  // Affect table file name
  m_TableFileName = f_ScannerTable;

  // Get table file name if no table name is given
  if (m_TableFileName=="")
  {
    char* tmp_string = getenv("RECON_TABLES");
    if (tmp_string==NULL)
    {
      LogCerr ("***** oTableHRRT::Constructor() -> Environment variable 'RECON_TABLES' is not set !" << endl);
      LogCerr ("                                   Please set it to point to the appropriate tables." << endl);
      exit(1);
    }
    string table_path = (string)tmp_string;
    m_TableFileName = table_path+"/hrrt_rebinner.lut";
  }

  if (m_Verbose>=1) LogCout ("oTableHRRT::ReadRingTable() -> Read table from file '" << m_TableFileName << "'" << endl);

  // Open the file
  FILE* flut = fopen(m_TableFileName.c_str(),"rb");
  if (flut==NULL)
  {
    LogCerr ("***** oTableHRRT::ReadRingTable() -> Input table file '" << m_TableFileName << "' is missing or corrupted !" << endl);
    return 1;
  }

  // Read it now
  if (m_Verbose>=2) LogCout ("  --> Read sino2D component ..." << endl);
  for (int mp=1; mp<=20; mp++)
  {
    for (int ax=0; ax<NXCRYS*NDOIS; ax++)
    {
      if (fread(mp_LUTsino[mp][ax],sizeof(SOL),NXCRYS*NDOIS,flut) != NXCRYS*NDOIS)
      {
        LogCerr ("***** oTableHRRT::ReadRingTable() -> Failed to read some data into LUT file !" << endl);
        fclose(flut);
        return 1;
      }
    }
  }
  if (m_Verbose>=2) LogCout ("  --> Read zpos component ..." << endl);
  if (fread(mp_LUTzpos,sizeof(float),NYCRYS,flut) != NYCRYS)
  {
    LogCerr ("***** oTableHRRT::ReadRingTable() -> Failed to read some data into LUT file !" << endl);
    fclose(flut);
    return 1;
  }
  if (m_Verbose>=2) LogCout ("  --> Read zpos2 component ..." << endl);
  if (fread(mp_LUTzpos2,sizeof(float),NYCRYS,flut) != NYCRYS)
  {
    LogCerr ("***** oTableHRRT::ReadRingTable() -> Failed to read some data into LUT file !" << endl);
    fclose(flut);
    return 1;
  }

  // Close the file
  fclose(flut);

  /*
  // Show the content of LUT tables
  for (int mp=1; mp<=20; mp++) for (int xx1=0; xx1<NXCRYS*NDOIS; xx1++) for (int xx2=0; xx2<NXCRYS*NDOIS; xx2++)
  {
    cout << "mp: " << mp << " | xx1: " << xx1 << " | xx2: " << xx2 << " | .z: " << mp_LUTsino[mp][xx1][xx2].z << " | .d: " << mp_LUTsino[mp][xx1][xx2].d << endl;
    getchar();
  }
  */

  // Ending
  return 0;
}
// ==========================================================================================================================================
// int ComputeSpanTables()
//   --> This function computes the span tables and returns 0 upon success and other value otherwise.
//       (code copied from segment_info class from gen_delays_lib)
// ==========================================================================================================================================
int oTableHRRT::ComputeSpanTables()
{
  if (m_Verbose>=1) LogCout ("oTableHRRT::ComputeSpanTables() -> Compute the span conversion tables" << endl);

  // Everything is copied from the original Merence's code

  int nbRings = 104;
  int span = 9;
  double plane_sep = 1.21875;
  double crystal_radius = 234.5;
  // First function
  m_nsegs=2*(m_MaxRingDiff/span)+1;
  int *m_segzoffset_span9=(int *)calloc(m_nsegs,sizeof(int));
  int np = 2*nbRings-1;
  int sp2 = (span+1)/2;
  m_segz0   = (int*) malloc(m_nsegs*sizeof(int));
  m_segzmax = (int*) malloc(m_nsegs*sizeof(int));
  int* segnz   = (int*) malloc(m_nsegs*sizeof(int));
  int* segzoff = (int*) malloc(m_nsegs*sizeof(int));
  int nplanes = 0;
  for (int i=0; i<m_nsegs; i++)
  {
    int segnum = (1-2*(i%2))*(i+1)/2;
    if (i==0) m_segz0[0]=0;
    else m_segz0[i]=sp2+span*((i-1)/2);
    segnz[i]=np-2*m_segz0[i];
    if (i==0) segzoff[0]=0;
    else segzoff[i] = segzoff[i-1] + segnz[i-1];
    nplanes += segnz[i];
    m_segzmax[i]=m_segz0[i]+segnz[i]-1;
    m_segzoffset_span9[i]=-m_segz0[i]+segzoff[i];
  }
  free(m_segz0);
  free(segnz);
  free(segzoff);
  free(m_segzmax);
  // Second function
  m_nsegs=2*(m_MaxRingDiff/3)+1;
  m_segzoffset=(int *)calloc(m_nsegs,sizeof(int));
  sp2=(3+1)/2;
  m_segz0=(int*) malloc(m_nsegs*sizeof(int));
  segnz=(int*) malloc(m_nsegs*sizeof(int));
  segzoff=(int*) malloc(m_nsegs*sizeof(int));
  m_segzmax=(int*) malloc(m_nsegs*sizeof(int));
  for (int i=0; i<m_nsegs; i++)
  {
    int segnum=(1-2*(i%2))*(i+1)/2;
    if (i==0) m_segz0[0]=0;
    else m_segz0[i]=sp2+3*((i-1)/2);
    segnz[i]=np-2*m_segz0[i];
    if (i==0) segzoff[0]=0;
    else segzoff[i] = segzoff[i-1] + segnz[i-1];
    m_segzmax[i]=m_segz0[i]+segnz[i]-1;
    m_segzoffset[i]=-m_segz0[i]+segzoff[i];
  }
  free(segnz);
  free(segzoff);
  for (int i=0;i<m_nsegs;i++) m_segzoffset[i]=m_segzoffset_span9[mp_Span3to9ConversionTable[i]];
  free(m_segzoffset_span9);
  // LUT span part
  m_segplane = (short **) calloc(63,sizeof(short *));
  for (int i=0; i<63; i++)
  {
    m_segplane[i]=(short *) calloc(NYCRYS*2-1,sizeof(short));
  }
  for (int i=0; i<63; i++)
  {
    for(int plane=0; plane<NYCRYS*2-1; plane++)
    {
      if (i>m_nsegs-1)
      {
        m_segplane[i][plane]=-1;
        continue;
      }
      if (plane < m_segz0[i]) m_segplane[i][plane]=-1;
      if (plane > m_segzmax[i]) m_segplane[i][plane]=-1;
      if (m_segplane[i][plane]!=-1) m_segplane[i][plane] = plane+m_segzoffset[i];
    }
  }
  // Ending
  return 0;
}
// ==========================================================================================================================================
// int ComputeFinalTables()
//   --> This function computes the final tables, namely elem, view, sino to friendly stuff (reverse order with respect to rebinnin approach)
// ==========================================================================================================================================
int oTableHRRT::ComputeFinalTables()
{
  if (ComputeElemViewRingTables())
  {
    LogCerr ("***** oTableHRRT::ComputeFinalTables() -> An error occured while calling function ComputeElemViewRingTables() !" << endl);
    return 1;
  }
  if (ComputeSino2PlaneSegmentTables())
  {
    LogCerr ("***** oTableHRRT::ComputeFinalTables() -> An error occured while calling function ComputeSino2PlaneSegmentTables() !" << endl);
    return 1;
  }

  // End
  return 0;
}
// ==========================================================================================================================================
// int ComputeElemViewRingTables()
//   --> This function computes tables for reverse binning (elem, view 2 indices)
// ==========================================================================================================================================
int oTableHRRT::ComputeElemViewRingTables()
{
  if (m_Verbose>=1) LogCout ("oTableHRRT::ComputeElemViewRingTables() -> Compute the conversion tables" << endl);

  // Allocate tables
  mp_ElemView2Mp  = (int***)malloc(m_NbElem*sizeof(int**));
  mp_ElemView2XX1 = (int***)malloc(m_NbElem*sizeof(int**));
  mp_ElemView2XX2 = (int***)malloc(m_NbElem*sizeof(int**));
  mp_ElemView2NbLors = (int**)malloc(m_NbElem*sizeof(int*));
  for (int e=0; e<m_NbElem; e++)
  {
    mp_ElemView2Mp [e] = (int**)malloc(m_NbView*sizeof(int*));
    mp_ElemView2XX1[e] = (int**)malloc(m_NbView*sizeof(int*));
    mp_ElemView2XX2[e] = (int**)malloc(m_NbView*sizeof(int*));
    mp_ElemView2NbLors[e] = (int*)calloc(m_NbView,sizeof(int));
    for (int v=0; v<m_NbView; v++)
    {
      mp_ElemView2Mp [e][v] = (int*)calloc(1,sizeof(int));
      mp_ElemView2XX1[e][v] = (int*)calloc(1,sizeof(int));
      mp_ElemView2XX2[e][v] = (int*)calloc(1,sizeof(int));
    }
  }

  // Loop on heads
  for (int mp=1; mp<=20; mp++)
  {
    // Loop on Layer and Trans 1
    for (int l1=0; l1<NDOIS; l1++) for (int t1=0; t1<NXCRYS; t1++)
    {
      // Loop on Layer and Trans 2
      for (int l2=0; l2<NDOIS; l2++) for (int t2=0; t2<NXCRYS; t2++)
      {
        // Get the 2D sinogram bin
        int xx1 = t1 + NXCRYS*l1;
        int xx2 = t2 + NXCRYS*l2;
        int addr_2D = mp_LUTsino[mp][xx1][xx2].nsino;
        if (addr_2D==-1) continue;
        // Decompose the 2D bin into elem, view
        int elem = addr_2D % m_NbElem;
        int view = addr_2D / m_NbElem;
        // Realloc corresponding table entries
        mp_ElemView2NbLors[elem][view]++;
        mp_ElemView2Mp [elem][view] = (int*)realloc(mp_ElemView2Mp[elem][view],mp_ElemView2NbLors[elem][view]*sizeof(int));
        mp_ElemView2XX1[elem][view] = (int*)realloc(mp_ElemView2XX1[elem][view],mp_ElemView2NbLors[elem][view]*sizeof(int));
        mp_ElemView2XX2[elem][view] = (int*)realloc(mp_ElemView2XX2[elem][view],mp_ElemView2NbLors[elem][view]*sizeof(int));
        // Add the corresponding values
        mp_ElemView2Mp [elem][view][mp_ElemView2NbLors[elem][view]-1] = mp;
        mp_ElemView2XX1[elem][view][mp_ElemView2NbLors[elem][view]-1] = xx1;
        mp_ElemView2XX2[elem][view][mp_ElemView2NbLors[elem][view]-1] = xx2;
      }
    }
  }

  // End
  return 0;
}
// ==========================================================================================================================================
// int ComputeSino2PlaneSegmentTables()
//   --> This function computes tables for reverse binning (sino 2 plane and segment numbers)
// ==========================================================================================================================================
int oTableHRRT::ComputeSino2PlaneSegmentTables()
{
  if (m_Verbose>=1) LogCout ("oTableHRRT::ComputeSino2PlaneSegmentTables() -> Compute the conversion tables" << endl);

  // Allocate tables
  mp_Sino2Segment = (int*)malloc(m_NbSino*sizeof(int));
  mp_Sino2Plane   = (int*)malloc(m_NbSino*sizeof(int));

  // Number of segments
  int nb_segments = 2 * (m_MaxRingDiff/m_Span) + 1;
  int nb_planes   = NYCRYS;

  // Compute them
  for (int segment=0; segment<nb_segments; segment++)
  {
    // Relative segment number
    
    for (int plane=0; plane<nb_planes; plane++)
    {
      // If this combination exists (and is subsequentely unique)
      if (m_segplane[segment][plane] != -1)
      {
        int offset = m_segzoffset[segment];
        int sino   = plane + offset;
        mp_Sino2Plane[sino]   = plane;
        mp_Sino2Segment[sino] = segment;
      }
    }
  }

  // End
  return 0;
}
// ==========================================================================================================================================
// int ComputeSino2RingsTables()
//   --> This function computes the ring indices corresponding to the sinogram bins, at the FOV center.
// ==========================================================================================================================================
/*
int oTableHRRT::ComputeSino2RingsTables()
{
  if (m_Verbose>=1) LogCout ("oTableHRRT::ComputeSino2RingsTables() -> Compute the conversion tables" << endl);

  // Allocate tables
  mp_Sino2Ring1 = (int*)malloc(m_NbSino*sizeof(int));
  mp_Sino2Ring2 = (int*)malloc(m_NbSino*sizeof(int));



}
*/

// ==========================================================================================================================================
// int GetLORsFromSinoBin()
//   --> This function computes using brut-force, the number of LORs contributing to a sinogram bin, and associated crystal indices.
// ==========================================================================================================================================
int oTableHRRT::GetLORsFromSinoBin( int f_BinSino, int f_BinView, int f_BinElem, unsigned int* fp_ID1, unsigned int* fp_ID2 )
{
  // Number of LORs
  int nb_lors = 0;

  // Calculate max ring diff for this sino bin
  int segment = mp_Sino2Segment[f_BinSino];

//2  int max_ring_diff_in_segment0 = m_Span / 2 + 1;
//2  int max_ring_diff = segment/2 * m_Span + max_ring_diff_in_segment0 + 1;

  // Loop on transaxial possibilities
  for (int i=0; i<mp_ElemView2NbLors[f_BinElem][f_BinView]; i++)
  {
    // Local variables on stack
    int mp  = mp_ElemView2Mp[f_BinElem][f_BinView][i];
    int xx1 = mp_ElemView2XX1[f_BinElem][f_BinView][i];
    int xx2 = mp_ElemView2XX2[f_BinElem][f_BinView][i];

    // Get z and d from original tables
    float z = mp_LUTsino[mp][xx1][xx2].z;
    float d = mp_LUTsino[mp][xx1][xx2].d;

    // Brut-force loop on axial possibilities
    for (int ax1=0; ax1<NYCRYS; ax1++)
    {
      float cay = mp_LUTzpos2[ax1];
      float dz2_orig = -mp_LUTzpos[ax1];

/*
      // Opt 1
      int loop_start = max( 0 , ax1 - 1 - m_MaxRingDiff );
      int loop_stop  = min( NYCRYS , ax1 + 2 + m_MaxRingDiff );
*/
/*
      // Opt 2
      int loop_start = max( 0 , ax1 - 1 - max_ring_diff );
      int loop_stop  = min( NYCRYS , ax1 + 2 + max_ring_diff );
*/
      // Opt 0
      int loop_start = 0;
      int loop_stop  = NYCRYS;

      for (int ax2=loop_start; ax2<loop_stop; ax2++)
      {
        float dz2 = dz2_orig + mp_LUTzpos[ax2];
        int plane = (int)(cay + z*dz2);
        if (plane!=mp_Sino2Plane[f_BinSino]) continue;
        float seg = 0.5 + dz2*d;
        int segnum = (int)seg;
        if (seg<0) segnum = 1 - (segnum<<1);
        else segnum = segnum<<1;
        // Do we finally get one ?
        if (segnum==segment)
        {
          if (nb_lors==100)
          {
            LogCerr ("!!!!!!!!!!!!!!!!!! LORS 100 !!!!!!!!!!!!!!!!!!!" << endl << endl << endl << flush);
            exit(1);
          }
          // Increment number of LORs and affect
          fp_ID1[nb_lors] = (xx1/NXCRYS)*NCRYSINLAYER + ax1*NTRANS + mp_HeadPairsReversed[mp][0]*NXCRYS + xx1%NXCRYS;
          fp_ID2[nb_lors] = (xx2/NXCRYS)*NCRYSINLAYER + ax2*NTRANS + mp_HeadPairsReversed[mp][1]*NXCRYS + xx2%NXCRYS;
          nb_lors++;
        }
      }
    }
  }

cout << "Nb ring LOR: " << mp_ElemView2NbLors[f_BinElem][f_BinView] << " | Final number of LORs: " << nb_lors << endl;
  // Return the number of LORs
  return nb_lors;
}
// ==========================================================================================================================================
// int ComputeFullTable()
//   --> This function computes the full conversion table from sinogram bin to contributing LORs (need a huge amount of RAM !).
// ==========================================================================================================================================
int oTableHRRT::ComputeFullTable()
{
  if (m_Verbose>=1) LogCout ("oTableHRRT::ComputeFullTables() -> Compute the conversion tables" << endl);

  // Allocate tables
  mp_Bin2CrystalID1 = (int****)malloc(m_NbSino*sizeof(int***));
  mp_Bin2CrystalID2 = (int****)malloc(m_NbSino*sizeof(int***));
  mp_Bin2NbLors = (int***)malloc(m_NbSino*sizeof(int**));
  for (int s=0; s<m_NbSino; s++)
  {
    mp_Bin2CrystalID1[s] = (int***)malloc(m_NbView*sizeof(int**));
    mp_Bin2CrystalID2[s] = (int***)malloc(m_NbView*sizeof(int**));
    mp_Bin2NbLors[s] = (int**)malloc(m_NbView*sizeof(int*));
    for (int v=0; v<m_NbView; v++)
    {
      mp_Bin2CrystalID1[s][v] = (int**)malloc(m_NbElem*sizeof(int*));
      mp_Bin2CrystalID2[s][v] = (int**)malloc(m_NbElem*sizeof(int*));
      mp_Bin2NbLors[s][v] = (int*)calloc(m_NbElem,sizeof(int*));
      for (int e=0; e<m_NbElem; e++)
      {
        mp_Bin2CrystalID1[s][v][e] = (int*)malloc(1*sizeof(int));
        mp_Bin2CrystalID2[s][v][e] = (int*)malloc(1*sizeof(int));
      }
    }
  }

  // Loop on first crystal index
  for (int id1=0; id1<NCRYS; id1++)
  {
    // Decompose the crystal index into usual indices (head 0 is on top-flat)
    int l1 = id1 / NCRYSINLAYER;
    int rest1 = id1 % NCRYSINLAYER;
    int a1 = rest1 / NTRANS;
    rest1 = rest1 % NTRANS;
    int h1 = rest1 / NXCRYS;
    int t1 = rest1 % NXCRYS;

    // Loop on second crystal index
    for (int id2=id1+1; id2<NCRYS; id2++)
    {
      // Decompose the crystal index into usual indices (head 0 is on top-flat)
      int l2 = id2 / NCRYSINLAYER;
      int rest2 = id2 % NCRYSINLAYER;
      int a2 = rest2 / NTRANS;
      rest2 = rest2 % NTRANS;
      int h2 = rest2 / NXCRYS;
      int t2 = rest2 % NXCRYS;

      // Get the head pair index
      int mp = abs(h1-h2);
      if (mp<2) continue;
      if (h1<h2) mp = mp_HeadPairs[h1][h2];
      else mp = mp_HeadPairs[h2][h1];

      // Get the sinogram bin
      int elem, view, sino;
      if (GetSinogramBinFromCrystalIndices( mp, l1, l2, t1, t2, a1, a2, &sino, &view, &elem )) continue;

      // Add it
      mp_Bin2NbLors[sino][view][elem]++;
      mp_Bin2CrystalID1[sino][view][elem] = (int*)realloc(mp_Bin2CrystalID1[sino][view][elem],mp_Bin2NbLors[sino][view][elem]*sizeof(int));
      mp_Bin2CrystalID2[sino][view][elem] = (int*)realloc(mp_Bin2CrystalID1[sino][view][elem],mp_Bin2NbLors[sino][view][elem]*sizeof(int));
      mp_Bin2CrystalID1[sino][view][elem][mp_Bin2NbLors[sino][view][elem]] = id1;
      mp_Bin2CrystalID2[sino][view][elem][mp_Bin2NbLors[sino][view][elem]] = id2;
    }
  }

  // End
  return 0;
}
// ==========================================================================================================================================
// int GetSinogramBinFromCrystalIndices()
//   --> This function computes the sinogram bin indices from two crystal indices.
//       It returns 0 upon success, 1 if not in maxrd range, 2 if the bin 2D sinogram index is wrong, 3 if the segment is out of range,
//       and 4 if the plane is out of range.
//       (the code is copied from the rebin_event() from LM_Rebinner.cpp
// ==========================================================================================================================================
int oTableHRRT::GetSinogramBinFromCrystalIndices( int f_Mp, int f_Layer1, int f_Layer2,
                                                  int f_Trans1, int f_Trans2, int f_Axial1, int f_Axial2,
                                                  int* f_BinSino, int* f_BinView, int* f_BinElem )
{
  // If the head pair index is negative, just call the fonction inversing the 1 and 2
  if (f_Mp<0) return GetSinogramBinFromCrystalIndices( -f_Mp, f_Layer2, f_Layer1, f_Trans2, f_Trans1, f_Axial2, f_Axial1, f_BinSino, f_BinView, f_BinElem );

  // Calculate the global indices in the 2D sinogram
  int xx1 = f_Trans1 + NXCRYS*f_Layer1;
  int xx2 = f_Trans2 + NXCRYS*f_Layer2;

  // Get the sino 2D addr
  int addr_2D = mp_LUTsino[f_Mp][xx1][xx2].nsino;
  if (addr_2D==-1) return 2;
  *f_BinElem = addr_2D % m_NbElem;
  *f_BinView = addr_2D / m_NbElem;

  // Things copied from the rebin_event() from LM_Rebinner.cpp
  float cay = mp_LUTzpos2[f_Axial1];
  float dz2 = mp_LUTzpos[f_Axial2] - mp_LUTzpos[f_Axial1];
  float z = mp_LUTsino[f_Mp][xx1][xx2].z;
  float d = mp_LUTsino[f_Mp][xx1][xx2].d;
  int plane = (int)(cay + z*dz2);
  float seg = 0.5 + dz2*d;
  int segnum = (int)seg;
  if (seg<0) segnum = 1 - (segnum<<1);
  else segnum = segnum<<1;
  if (segnum >= m_nsegs) return 3;
  if (m_segplane[segnum][plane] != -1)
  {
    int offset = m_segzoffset[segnum];
    *f_BinSino = plane + offset;
    return 0;
  }
  else return 4;

  // Ending
  return 0;
}
// ==========================================================================================================================================
// int GetHeadIndicesFromHeadPairIndex()
//   --> This function gets the head indices from a head pair index.
//       It returns 0 upon success, 1 if the pair index does not exist.
// ==========================================================================================================================================
int oTableHRRT::GetHeadIndicesFromHeadPairIndex( int f_Mp, int* f_Head1, int* f_Head2 )
{
  // Check if the pair index is out of range
  if (f_Mp<1 || f_Mp>20) return 1;
  // Else just get the head indices
  *f_Head1 = mp_HeadPairsReversed[f_Mp][0];
  *f_Head2 = mp_HeadPairsReversed[f_Mp][1];
  return 0;
}
// ==========================================================================================================================================
// int GetHeadPairIndexFromHeadIndices()
//   --> This function gets the head pair index from head indices.
//       If the index 1 is higher than the second, the Mp returned is set negative.
//       It returns 0 upon success, 1 if the pair index does not exist.
// ==========================================================================================================================================
int oTableHRRT::GetHeadPairIndexFromHeadIndices( int f_Head1, int f_Head2, int* f_Mp )
{
  // Check the head order
  if (f_Head1<=f_Head2)
  {
    *f_Mp = mp_HeadPairs[f_Head1][f_Head2];
    // If the pair does not exist
    if (*f_Mp==-1)
    {
      *f_Mp = 0;
      return 1;
    }
    else return 0;
  }
  else
  {
    *f_Mp = mp_HeadPairs[f_Head2][f_Head1];
    if (*f_Mp==-1)
    {
      *f_Mp = 0;
      return 1;
    }
    else
    {
      *f_Mp = -(*f_Mp);
      return 0;
    }
  }
}

