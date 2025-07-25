#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <ctime>
#include <fstream>
#include "oTableVision600.hh"
#include "oOutputManager.hh"
#include "oScanner.hh"
using namespace std;

// ==========================================================================================================================================
// Constructor
// ==========================================================================================================================================
oTableVision600::oTableVision600( oScanner* fp_Scanner, int f_Verbose )
{
  // Verbose
  m_Verbose = f_Verbose;
  if (m_Verbose>=1) LogCout ("oTableVision600::Constructor() -> Initialize vision600 conversion table" << endl);

  // Sinogram dimensions and characteristics
  m_NbSino = 815; 
  m_Span = 19;    
  m_NbElem = 520; 
  m_NbView = 50; 
  m_NbCrystalsPerRingWithGaps = 38*21; 
  m_Mash = 8;   
  m_MaxRingDiff = 79; 

  // Scanner model
  mp_Scanner = fp_Scanner;
}
// ==========================================================================================================================================
// Destructor
// ==========================================================================================================================================
oTableVision600::~oTableVision600()
{

  //Destruction of the lists for sino index to ring pairs and vice-versa
    
  if (m_Verbose>=1) LogCout ("oTableVision::Destructor() -> Destroy conversion table" << endl);

  if (mp_ListRingPairsBySinoIndex && mp_NbRingPairsBySinoIndex)
  {
    for (int s(0); s<m_NbSino; s++)
    {
      if (mp_ListRingPairsBySinoIndex[s])
      {
        for (int r(0); r<mp_NbRingPairsBySinoIndex[s]; r++)
        {
          if (mp_ListRingPairsBySinoIndex[s][r])
          {
            free(mp_ListRingPairsBySinoIndex[s][r]);
          }
        }
        free(mp_ListRingPairsBySinoIndex[s]);
      }
    }
    free(mp_ListRingPairsBySinoIndex);
    free(mp_NbRingPairsBySinoIndex);
  }
  
  int nb_rings = mp_Scanner->GetNbTotalAxialCrystalsWithGaps();
  
  if (mp_RingPair2SinogramIndex)
  {
    for(int s(0);s<nb_rings;s++)
    {
      if (mp_RingPair2SinogramIndex[s])
      {
        free(mp_RingPair2SinogramIndex[s]);
      }
    }
    free(mp_RingPair2SinogramIndex);
  }
  
  //Destruction of the lists for crystals indices to elem,view and vice-versa
  if (mp_ElemView2CrystalID1)
  {
    for (int e(0); e<m_NbElem; e++)
    {
      if (mp_ElemView2CrystalID1[e])
      {
        for (int v(0); v<m_NbView; v++)
        {
          if (mp_ElemView2CrystalID1[e][v])
          {
            free(mp_ElemView2CrystalID1[e][v]);
          }
        }
        free(mp_ElemView2CrystalID1[e]);
      }
    }
    free(mp_ElemView2CrystalID1);
  }


  if (mp_ElemView2CrystalID2)
  {
    for (int e(0); e<m_NbElem; e++)
    {
      if (mp_ElemView2CrystalID2[e])
      {
        for (int v(0); v<m_NbView; v++)
        {
          if (mp_ElemView2CrystalID2[e][v])
          {
            free(mp_ElemView2CrystalID2[e][v]);
          }
        }
        free(mp_ElemView2CrystalID2[e]);
      }
    }
    free(mp_ElemView2CrystalID2);
  }


  int nbTotalTransCrystalsWithGaps(mp_Scanner->GetNbTotalTransCrystalsWithGaps());

  if (mp_CrystalIDs2Elem)
  {
    for(int s(0);s<nbTotalTransCrystalsWithGaps;s++)
    {
      if (mp_CrystalIDs2Elem[s])
      {
        free(mp_CrystalIDs2Elem[s]);
      }
    }
    free(mp_CrystalIDs2Elem);
  }

  if (mp_CrystalIDs2View)
  {
    for(int s(0);s<nbTotalTransCrystalsWithGaps;s++)
    {
      if (mp_CrystalIDs2View[s])
      {
        free(mp_CrystalIDs2View[s]);
      }
    }
    free(mp_CrystalIDs2View);
  }


  //We do not delete the scanner because it's already deleted in oSinoModeCreation
  
  //if (mp_Scanner) delete mp_Scanner;
  
}
// ==========================================================================================================================================
// int ReadTable()
//   --> This function read the conversion table from the given file and returns 0 upon success and other value otherwise.
// ==========================================================================================================================================
int oTableVision600::ReadSpanTable(string f_ScannerTable)
{
  // Affect table file name
  m_TableFileName = f_ScannerTable;

/* Obsolete as the tab is built in below in the code
  // Get table file name if no table name is given
  if (m_TableFileName=="")
  {
    char* tmp_string = getenv("RECON_TABLES");
    if (tmp_string==NULL)
    {
      LogCerr ("***** oTableVision600::Constructor() -> Environment variable 'RECON_TABLES' is not set !" << endl);
      LogCerr ("                                       Please set it to point to the appropriate tables." << endl);
      exit(1);
    }
    string table_path = (string)tmp_string;
    m_TableFileName = table_path+"/Vision600_4rings.txt";
  }
*/

  if (m_Verbose>=1) LogCout ("oTableVision600::ReadSpanTable() -> Read table from file '" << m_TableFileName << "'" << endl);

  // Allocation
  mp_NbRingPairsBySinoIndex = (unsigned short int*)calloc(m_NbSino,sizeof(unsigned short int));
  mp_ListRingPairsBySinoIndex = (unsigned int***)calloc(m_NbSino,sizeof(unsigned int**));

  //for (int s=0; s<m_NbSino; s++) mp_ListRingPairsBySinoIndex[s] = (unsigned int**)malloc(1*sizeof(unsigned int*));
  int nb_rings = mp_Scanner->GetNbTotalAxialCrystalsWithGaps();
  mp_RingPair2SinogramIndex = (int**)calloc(nb_rings,sizeof(int*));

  for (int r=0; r<nb_rings; r++) mp_RingPair2SinogramIndex[r] = (int*)calloc(nb_rings,sizeof(int));

  // -------------------------------------------------------------------------------------------------------------------
  // Read the file if not empty
  // -------------------------------------------------------------------------------------------------------------------
  if (m_TableFileName!="")
  {
    // Open the file
    ifstream ftab(m_TableFileName.c_str());
    if (!ftab)
    {
      LogCerr ("***** oTableVision600::ReadSpanTable() -> Input table file '" << m_TableFileName << "' is missing or corrupted !" << endl);
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
        // alloc the list of ring pairs for this sinogram
        mp_ListRingPairsBySinoIndex[s] = (unsigned int**)calloc(mp_NbRingPairsBySinoIndex[s],sizeof(unsigned int*));
        // Decode the ring indices
        unsigned int ring1, ring2;
        if (sscanf(buffer.c_str(),"<%u,%u>",&ring1,&ring2) != 2)
        {
          LogCerr ("***** oTableVision600::ReadSpanTable() -> Error decoding a ring pair \"" << buffer << "\" !" << endl);
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
  }
  // -------------------------------------------------------------------------------------------------------------------
  // Otherwise, build it! This part is a generic part allowing to build the lists for sino index to ring pairs and vice-versa from the caracteristics of the scanner (span, max ring diff...).
  // -------------------------------------------------------------------------------------------------------------------
  
  else
  {
  
    //The general idea is to build the lists by a loop over the segments, inside this loop, there is a loop over the sinograms, and inside the sinograms, one over the ring pairs.
    //We first allow any ring index, even if this index doesn't really exists (we can have negative indices, ...). This makes the code much easier to write. Before wrighting a component in a list, we check that this index really exists, like this we forbid at the same time negative and too large ring indeces, as well as taking care of gaps.
    
    
    int sino(0); //number of sinograms already done
    int segment(0); //number of segments already done
    // The four following indices represent the ring pair limits of a specific sinogram ie In the loop over the ring pairs inside a specific sinogram, the first ring pair is (index1,index2), and the last one is (index3,index4).
    int index1;
    int index2;
    int index3;
    int index4;
    int ring1;
    int ring2;
    bool isNotEmpty(false);
    int const nbTotalAxialCrystalsWithGaps(mp_Scanner -> GetNbTotalAxialCrystalsWithGaps());
    int const nbAxialCrystals(mp_Scanner -> GetNbAxialCrystals());
    int const nbAxialGapsBetweenBlocks(mp_Scanner -> GetNbAxialGapsBetweenBlocks());
    int const nbAxialCrystalsPerBlockWithGaps(nbAxialCrystals + nbAxialGapsBetweenBlocks);
    
    int oddSpanComponent; // eg 5 if the span is 11, or 3 if the span is 5
    int evenSpanComponent; // eg 6 if the span is 11, or 2 if the span is 5
    
    
    // computation of the "span components"
    if ((m_Span % 4) == 1)
    {
        oddSpanComponent = (m_Span+1)/2;
        evenSpanComponent = (m_Span-1)/2;
    }
    else if ((m_Span % 4) == 3)
    {
        oddSpanComponent = (m_Span-1)/2;
        evenSpanComponent = (m_Span+1)/2;
    }
    
    
    
    //Loop over the segments. We first cover the segment 0, then 1, then -1, then 2, etc.
    while(abs(segment)*m_Span-(m_Span-1)/2<=m_MaxRingDiff)
    {
    
    
        //Initialization of the indices
        if (segment>0)
        {
            index1 = segment*m_Span-(m_Span-1)/2;
            index2 = 0;
            if ((index1+index2+segment)%2==1)
            {
                index3=index1+evenSpanComponent-1;
                index4=index2-evenSpanComponent+1;
            }
            else
            {
                index3=index1+oddSpanComponent-1;
                index4=index2-oddSpanComponent+1;
            }
            
        }
        else if (segment<0)
        {
            index3 = 0;
            index4 = -segment*m_Span-(m_Span-1)/2;
            if ((index3+index4+segment)%2==1)
            {
                index1=index3-evenSpanComponent+1;
                index2=index4+evenSpanComponent-1;
            }
            else
            {
                index1=index3-oddSpanComponent+1;
                index2=index4+oddSpanComponent-1;
            }
        }
        else
        {
            index1 = -(oddSpanComponent-1)/2;
            index2 = (oddSpanComponent-1)/2;
            if ((index1+index2+segment)%2==1)
            {
                index3=index1+evenSpanComponent-1;
                index4=index2-evenSpanComponent+1;
            }
            else
            {
                index3=index1+oddSpanComponent-1;
                index4=index2-oddSpanComponent+1;
            }
        }

        

        // Loop over the sinograms inside a segment (we first cover the sinograms whose ring pair indices are the lowest).
        while ((segment>0 && index1<nbTotalAxialCrystalsWithGaps) || (segment<0 && index4 < nbTotalAxialCrystalsWithGaps) || (segment==0 && index1+index2<=2*(nbTotalAxialCrystalsWithGaps-1)))
        {
            ring1 = index1;
            ring2 = index2;
            isNotEmpty=false; //IsNotEmpty checks if the sinogram really exists. ie if at least one ring pair belong to it.
            
            mp_NbRingPairsBySinoIndex[sino]=0;
            
            // Loop over the ring pairs inside a sinogram. We start by the ring pairs whose first ring index is the lowest
            while (ring1<=index3)
            {
                if (ring1>=0 && ring1<nbTotalAxialCrystalsWithGaps && ring1 % nbAxialCrystalsPerBlockWithGaps < nbAxialCrystals && ring2>=0 && ring2<nbTotalAxialCrystalsWithGaps && ring2 % nbAxialCrystalsPerBlockWithGaps < nbAxialCrystals && abs(ring2-ring1)<=m_MaxRingDiff)
                {
                    mp_NbRingPairsBySinoIndex[sino]++;
                    mp_RingPair2SinogramIndex[ring1][ring2] = sino;
                    isNotEmpty=true;
                }
                
                
                
                ring1++;
                ring2--;
                
            }
            if (isNotEmpty)
            {
                mp_ListRingPairsBySinoIndex[sino] = (unsigned int**)calloc((mp_NbRingPairsBySinoIndex[sino]),sizeof(unsigned int*));
               
                ring1=index1;
                ring2=index2;
                int nbPair(0); //integer representing the index of the pair inside the sinogram (form 0 to max 5)
            
                while (ring1<=index3) // We need this second loop because that's only after the first loop that we can know how much we have to allocate in mp_ListRingPairsBySinoIndex[sino]
                {
                    if (ring1>=0 && ring1<nbTotalAxialCrystalsWithGaps && ring1 % nbAxialCrystalsPerBlockWithGaps < nbAxialCrystals && ring2>=0 && ring2<nbTotalAxialCrystalsWithGaps && ring2 % nbAxialCrystalsPerBlockWithGaps < nbAxialCrystals && abs(ring2-ring1)<=m_MaxRingDiff)
                    {

                        mp_ListRingPairsBySinoIndex[sino][nbPair] = (unsigned int*)calloc(2,sizeof(unsigned int));
                        mp_ListRingPairsBySinoIndex[sino][nbPair][0] = ring1;
                        mp_ListRingPairsBySinoIndex[sino][nbPair][1] = ring2;
                        nbPair++;
                    }
                
                
                
                    ring1++;
                    ring2--;
                
                }
            }
            
            //incrementations of the indices, and of the number of sinograms already done
            
            if (isNotEmpty)
            {
                sino++;
            }
            
            if (index3-index1==m_Span/2-1)
            {
                index2++;
                index3++;
            }
            else
            {
                index1++;
                index4++;
            }
        }
        
        
        
        
        //incrementation of the segments
        if (segment>0)
        {
            segment=-segment;
        }
        else if (segment<0)
        {
            segment=(-segment)+1;
        }
        else
        {
            segment++;
        }
    }    
  }
  
  
  
 

  // Ending
  return 0;
}
// ==========================================================================================================================================
// int ComputeRingTables()
//   --> This function computes the ring conversion tables (crystals to elem,view) and returns 0 upon success and other value otherwise.
// ==========================================================================================================================================
int oTableVision600::ComputeRingTables()
{
  if (m_Verbose>=1) LogCout ("oTableVision600::ComputeRingTables() -> Compute the ring conversion tables" << endl);

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
unsigned short int oTableVision600::GetNbRingPairsBySinoIndex(int f_Sino)
{
  return mp_NbRingPairsBySinoIndex[f_Sino];
}
void oTableVision600::GetRingPairBySinoIndex(int f_Sino, int f_RingPair, unsigned int* f_Ring1, unsigned int* f_Ring2)
{
  *f_Ring1 = mp_ListRingPairsBySinoIndex[f_Sino][f_RingPair][0];
  *f_Ring2 = mp_ListRingPairsBySinoIndex[f_Sino][f_RingPair][1];
}
void oTableVision600::GetRingPairsBySinoIndex(int f_Sino, unsigned int* f_Ring1, unsigned int* f_Ring2)
{
  for (int r=0; r<mp_NbRingPairsBySinoIndex[f_Sino]; r++)
  {
    f_Ring1[r] = mp_ListRingPairsBySinoIndex[f_Sino][r][0];
    f_Ring2[r] = mp_ListRingPairsBySinoIndex[f_Sino][r][1];
  }
}
int oTableVision600::GetSinoIndexFromRingPair(unsigned int f_Ring1, unsigned int f_Ring2)
{
  return mp_RingPair2SinogramIndex[f_Ring1][f_Ring2];
}
void oTableVision600::GetCrystalIDsFromElemView(int f_Elem, int f_View, int* f_Crystal1, int* f_Crystal2)
{
  for (int m=0; m<m_Mash; m++)
  {
    *f_Crystal1 = mp_ElemView2CrystalID1[f_Elem][f_View][m];
    *f_Crystal2 = mp_ElemView2CrystalID2[f_Elem][f_View][m];
  }
}
void oTableVision600::GetElemViewFromCrystalIDs(int* f_Elem, int* f_View, int f_Crystal1, int f_Crystal2)
{
  *f_Elem = mp_CrystalIDs2Elem[f_Crystal1][f_Crystal2];
  *f_View = mp_CrystalIDs2View[f_Crystal1][f_Crystal2];
}
// ==========================================================================================================================================
// void Sino2Crystal()
//   --> This function computes the crystal coordinates (inside a ring) for given element and view.
// ==========================================================================================================================================
void oTableVision600::ElemView2CrystalIDs(int f_Elem, int f_View, int* f_Crystal1, int* f_Crystal2)
{

  // Loop on mash
  for (int m=0; m<m_Mash; m++)
  {
  
    // Compute coordinates
    int opti = m_NbElem/4;
    int det1_c = f_View*m_Mash + m;
    int det2_c = f_View*m_Mash + m + (m_NbCrystalsPerRingWithGaps/2);
    f_Crystal1[m] = det1_c + f_Elem/2     - opti;
    f_Crystal2[m] = det2_c - (f_Elem+1)/2 + opti;

    // Apply modulos
    
    f_Crystal1[m] = f_Crystal1[m] % m_NbCrystalsPerRingWithGaps;
    f_Crystal2[m] = f_Crystal2[m] % m_NbCrystalsPerRingWithGaps;
    
    if (f_Crystal1[m] < 0) f_Crystal1[m] += m_NbCrystalsPerRingWithGaps;
    else if (f_Crystal1[m] >= m_NbCrystalsPerRingWithGaps) f_Crystal1[m] -= m_NbCrystalsPerRingWithGaps;
    if (f_Crystal2[m] < 0) f_Crystal2[m] += m_NbCrystalsPerRingWithGaps;
    else if (f_Crystal2[m] >= m_NbCrystalsPerRingWithGaps) f_Crystal2[m] -= m_NbCrystalsPerRingWithGaps;
  }
}
// ==========================================================================================================================================
// int GetLORsFromSinoBin()
//   --> This function computes the number of LORs contributing to a sinogram bin, and associated crystal indices.
// ==========================================================================================================================================
int oTableVision600::GetLORsFromSinoBin( int f_BinSino, int f_BinView, int f_BinElem, unsigned int* fp_ID1, unsigned int* fp_ID2 )
{
  // The number of LORs
  int nb_lors = 0;
  
  // Loop on mashing
  for (int m=0; m<m_Mash; m++)
  {  
      
      // Get xpos of the crystal pair
      int xpos1 = mp_ElemView2CrystalID1[f_BinElem][f_BinView][m];
      int xpos2 = mp_ElemView2CrystalID2[f_BinElem][f_BinView][m];
    
      // Check if we are in a gap, we skip it but decrement the number of events stored
      int nbTransCrystalsInHead(mp_Scanner->GetNbTransCrystalsInHead());
      int nbGapsBetweenHeads(mp_Scanner->GetNbGapsBetweenHeads());
      int nbTransCrystalsInHeadWithGaps(nbTransCrystalsInHead+nbGapsBetweenHeads);
      int nbTransCrystals(mp_Scanner->GetNbTransCrystals());
      int nbTransGapsBetweenBlocks(mp_Scanner->GetNbTransGapsBetweenBlocks());
      int nbTransCrystalsInBlockWithGaps(nbTransCrystals+nbTransGapsBetweenBlocks);
      
      if ( (xpos1 % nbTransCrystalsInHeadWithGaps >= nbTransCrystalsInHead) || (xpos1 % nbTransCrystalsInHeadWithGaps) %  nbTransCrystalsInBlockWithGaps >= nbTransCrystals) continue;
      if ( (xpos2 % nbTransCrystalsInHeadWithGaps >= nbTransCrystalsInHead) || (xpos2 % nbTransCrystalsInHeadWithGaps) %  nbTransCrystalsInBlockWithGaps >= nbTransCrystals) continue;

      // Compute the crystal index on ring without gaps
      
      int NbTransBlocks(mp_Scanner->GetNbTransBlocks());
      
      xpos1 -= (xpos1 / nbTransCrystalsInHeadWithGaps)*(nbGapsBetweenHeads+nbTransGapsBetweenBlocks*(NbTransBlocks-1));
      xpos1 -= ( (xpos1 % nbTransCrystalsInHeadWithGaps)/nbTransCrystalsInBlockWithGaps )*nbTransGapsBetweenBlocks;
      xpos2 -= (xpos2 / nbTransCrystalsInHeadWithGaps)*(nbGapsBetweenHeads+nbTransGapsBetweenBlocks*(NbTransBlocks-1));
      xpos2-= ( (xpos2 % nbTransCrystalsInHeadWithGaps)/nbTransCrystalsInBlockWithGaps )*nbTransGapsBetweenBlocks;


      

      
      // Loop on span
      //for (int r=0; r<mp_NbRingPairsBySinoIndex[f_BinSino]; r++)
      //f_BinSino=0;
      
      for (int r=0; r<mp_NbRingPairsBySinoIndex[f_BinSino]; r++)
      {
        


        
        int ring1 = mp_ListRingPairsBySinoIndex[f_BinSino][r][0];

        int ring2 = mp_ListRingPairsBySinoIndex[f_BinSino][r][1];

        // Compute the ring indices without gaps (no need to check wether we are in a gap because gaps are already removed from the span table)
        int nbAxialCrystals(mp_Scanner->GetNbAxialCrystals());
        int nbAxialGapsBetweenBlocks(mp_Scanner->GetNbAxialGapsBetweenBlocks());
        int nbAxialCrystalsInBlockWithGaps(nbAxialCrystals+nbAxialGapsBetweenBlocks);
    
        ring1 -= ( ring1/(nbAxialCrystalsInBlockWithGaps) )*nbAxialGapsBetweenBlocks;
        ring2 -= ( ring2/(nbAxialCrystalsInBlockWithGaps) )*nbAxialGapsBetweenBlocks;

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

