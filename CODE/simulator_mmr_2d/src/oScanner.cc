#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <ctime>
#include <float.h>
#include "oScanner.hh"
#include "oMiscellaneous.hh"
#include "oOutputManager.hh"
using namespace std;

// ==========================================================================================================================================
// Constructor
// ==========================================================================================================================================
oScanner::oScanner(int f_ScannerModel, int f_Verbose)
{
  if (f_ScannerModel==SCANNER_HRPLUS)          GenericConstructor("hrplus"     ,f_Verbose);
  else if (f_ScannerModel==SCANNER_HRRT)       GenericConstructor("hrrt"       ,f_Verbose);
  else if (f_ScannerModel==SCANNER_BIOGRAPH)   GenericConstructor("biograph"   ,f_Verbose);
  else if (f_ScannerModel==SCANNER_BIOGRAPH2D) GenericConstructor("biograph2D" ,f_Verbose);
  else if (f_ScannerModel==SCANNER_VISION600) GenericConstructor("vision600" ,f_Verbose);
}
oScanner::oScanner(const string& f_ScannerName, int f_Verbose)
{
  GenericConstructor(f_ScannerName,f_Verbose);
}
void oScanner::GenericConstructor(const string& f_ScannerName, int f_Verbose)
{
  if (f_Verbose) LogCout ("oScanner::Constructor() -> Construct" << endl);

  if (f_ScannerName == "hrrt")
  {
    float mean_depth_of_interaction = 5.;
    Initialize( 8, 13, 9,
                8, 8, 2,
                1/*?*/, 0, 0,
                234.5, mean_depth_of_interaction,
                2.2, 2.2, 20.,
                0.2, 0.2, // 0.2375, 0.2375,
                0.5, 0.5,
                256, 288, 1.218750,
                312.0128, 312.0128, 252.2916,
                0., // No bed offset !
                256, 256, 207,
                SCANNER_HRRT, 64, f_Verbose
              );
  }
  else if (f_ScannerName == "vision600")
  {
    cerr << "***** scanner correct" << endl;
    float mean_depth_of_interaction = 8;
//    float mean_depth_of_interaction = 0.;
    Initialize( 38, 8, 1,
                10, 20, 1,
                1, 0, 0,
                410., mean_depth_of_interaction,
                3.2, 3.2, 20.,
                0., 0.,
                1, 0,
                520, 400, 1.6,
                726., 726., 263. ,
                0.,
                440, 440, 159,
                SCANNER_VISION600, 32, f_Verbose );
  }
  else if (f_ScannerName == "hrplus")
  {
    float mean_depth_of_interaction = 15.;
    Initialize( 72, 4, 1,
                8, 8, 1,
                0, 0, 0,
                412., mean_depth_of_interaction,
                4.75, 4.4, 30.,
                0.025, 0.0942,
                0.025, 0.,
                288, 288, 2.25,
                659.0208, 659.0208, 152.775,
                126.1, // 11 slices of overlap, so bed offset of 52 x 2.425mm slices
                256, 256, 63,
                SCANNER_HRPLUS, 32, f_Verbose );
  }
  
  else if (f_ScannerName == "biograph")
  {
//    float mean_depth_of_interaction = 10.;
    float mean_depth_of_interaction = 9.6;
/*
    Initialize( 48, 4, 1,
                13, 13, 1,
                1, 1, 0,
                427.6, mean_depth_of_interaction,
                4., 4., 20.,
//LAST USED                0.0, 0.01,
                0.017462963, 0.01,
//                0.054, 0.01, // From claire's macro
//                0.166666667, 0.166666667, // computed from official biograph publication TNS
//                4.108, 0.,
//LASTUSED                4., 0.,
                4.034925926, 0.,
                336, 336, 2.005,
                684.23712, 684.23712, 218.,
                148., // Bed offset of 74 x 2mm slices
//                149.998, // Bed offset of 74 x 2.027mm slices (this is what is written in the official biograph headers)
                256, 256, 109,
                SCANNER_BIOGRAPH, 32, f_Verbose );
*/
    Initialize( 48, 4, 1,
                13, 13, 1,
                1, 1, 0,
                427.6, mean_depth_of_interaction,
                4., 4., 20.,
                0.017462963, 0.01,
                4.034925926, 0.,
                336, 336, 2.005,
                684.23712, 684.23712, 220.943,
                149.998, // Bed offset of 74 x 2.027mm slices (this is what is written in the official biograph headers)
                256, 256, 109,
                SCANNER_BIOGRAPH, 32, f_Verbose );

  }
  else if (f_ScannerName == "castor")
  {
    float mean_depth_of_interaction = 10.;
    Initialize( 48, 4, 1,
                13, 13, 1,
                1, 1, 0,
                427.6, mean_depth_of_interaction,
                4., 4., 20.,
                0.017462963, 0.01,
                4.034925926, 0.,
                160, 112, 2.005,
                672., 672., 220.9437,
                149.998, // Bed offset of 74 x 2.027mm slices
                256, 256, 109,
                SCANNER_CASTOR, 32, f_Verbose );
  }
  else if (f_ScannerName == "biograph2D")
  {
    float mean_depth_of_interaction = 9.6;
    Initialize( 48, 1, 1,
                1, 13, 1,
                1, 0, 0,
                427.6, mean_depth_of_interaction,
                4., 4., 20.,
                0., 0.01, // From claire's macro
//                0.166666667, 0.166666667, // computed from official biograph publication TNS
                0., 0.,
                336, 336, 2.005,
                684.23712, 684.23712, 4.,
                4., // Bed offset of 1 x 4mm slices
                256, 256, 1,
                SCANNER_BIOGRAPH2D, 32, f_Verbose );
  }
  else if (f_ScannerName == "mmr2d")
  {
    float mean_depth_of_interaction = 10.;
    Initialize( 56, 1, 1,
                1, 8, 1,
                1, 0, 0,
                328, mean_depth_of_interaction,
                4., 4., 20.,
                0., 0.17252,
                0., 0.0,
                344, 252, 2.03125,
                717.67344, 717.67344, 4.,
                0., // Bed offset
                344, 344, 1,
                SCANNER_MMR2D, 32, f_Verbose );
  }
  else if (f_ScannerName == "inveon")
  {
    float mean_depth_of_interaction = 4.584;
    Initialize( 16, 4, 1,
                20, 20, 1,
                0, 0, 0,
                80.54, mean_depth_of_interaction,
                1.5, 1.5, 10.,
                0.092, 0.13,
                0., 0.,
                128, 160, 0.815,
                99.377152, 99.377152, 126.564,
                0.,
                256, 256, 159,
                SCANNER_INVEON, 32, f_Verbose );
  }
  else if (f_ScannerName == "focus")
  {
    float mean_depth_of_interaction = 4.584;
    Initialize( 42, 4, 1,
                12, 12, 1,
                0, 0, 0,
                130.673, mean_depth_of_interaction,
                1.5, 1.5, 10.,
                0.092, 0.13,
                0., 0.,
                288, 252, 0.815,
                108., 108., 75.62,
                0.,
                256, 256, 95,
                SCANNER_FOCUS, 32, f_Verbose );
  }
  else if (f_ScannerName == "signa")
  {
    float mean_depth_of_interaction = 3.3;
//    float mean_depth_of_interaction = 0.;
    Initialize( 28, 5, 4,
                9, 4, 1,
                0, 0, 0,
                317., mean_depth_of_interaction,
                5.316, 3.975, 25.,
                0., 0.,
                2.8, 0.3,
                357, 224, 2.,
                600., 600., 250.4,
                0.,
                192, 192, 109,
                SCANNER_SIGNA, 32, f_Verbose );
  }
  
  else
  {
    cerr << "***** a problem" << endl;
    LogCerr ("***** oScanner::oScanner() -> Unknown scanner name '" << f_ScannerName << "' or not yet supported oo!" << endl);
    exit(1);
  }
}
// ==========================================================================================================================================
// void Initialize()
//    --> This function is only called by the constructor to initialize all members.
// ==========================================================================================================================================
void oScanner::Initialize( int f_NbHeads, int f_NbAxialBlocks, int f_NbTransBlocks,
                           int f_NbAxialCrystals, int f_NbTransCrystals, int f_NbLayers,
                           int f_NbGapsBetweenHeads, int f_NbAxialGapsBetweenBlocks, int f_NbTransGapsBetweenBlocks,
                           float f_Radius, float f_MeanDOI,
                           float f_AxialCrystalSize, float f_TransCrystalSize, float f_DepthCrystalSize,
                           float f_AxialCrystalGap, float f_TransCrystalGap,
                           float f_AxialBlockGap, float f_TransBlockGap,
                           int f_NbElem, int f_NbView, float f_RadialBinSize,
                           float f_FOVSizeX, float f_FOVSizeY, float f_FOVSizeZ,
                           float f_BedOffset,
                           int f_NbVoxelX, int f_NbVoxelY, int f_NbVoxelZ,
                           int f_ScannerModel, int f_BitMode, int f_Verbose
                         )
{
  // Verbosity
  m_Verbose = f_Verbose;
  if (m_Verbose>=1) LogCout ("oScanner::Initialize() ..." << endl);

  // Scanner model
  m_ScannerModel = f_ScannerModel;
  m_BitMode = f_BitMode;

  // Dimensions affectations (sorry there is no check here ...)
  m_NbHeads = f_NbHeads;
  m_NbAxialBlocks = f_NbAxialBlocks;
  m_NbTransBlocks = f_NbTransBlocks;
  m_NbAxialCrystals = f_NbAxialCrystals;
  m_NbTransCrystals = f_NbTransCrystals;
  m_NbLayers = f_NbLayers;
  m_NbGapsBetweenHeads = f_NbGapsBetweenHeads;
  m_NbAxialGapsBetweenBlocks = f_NbAxialGapsBetweenBlocks;
  m_NbTransGapsBetweenBlocks = f_NbTransGapsBetweenBlocks;
  m_Radius = f_Radius;
  m_MeanDOI = f_MeanDOI;
  m_AxialCrystalSize = f_AxialCrystalSize;
  m_TransCrystalSize = f_TransCrystalSize;
  m_DepthCrystalSize = f_DepthCrystalSize;
  m_AxialCrystalGap = f_AxialCrystalGap;
  m_TransCrystalGap = f_TransCrystalGap;
  m_AxialBlockGap = f_AxialBlockGap;
  m_TransBlockGap = f_TransBlockGap;
  m_AxialBlockSize = ((float)m_NbAxialCrystals)*m_AxialCrystalSize + ((float)(m_NbAxialCrystals-1))*m_AxialCrystalGap;
  m_TransBlockSize = ((float)m_NbTransCrystals)*m_TransCrystalSize + ((float)(m_NbTransCrystals-1))*m_TransCrystalGap;
  m_TransHeadSize  = ((float)m_NbTransBlocks)*m_TransBlockSize + ((float)(m_NbTransBlocks-1))*m_TransBlockGap;
  m_AxialScannerSize = m_AxialBlockSize*((float)m_NbAxialBlocks) + m_AxialBlockGap*((float)(m_NbAxialBlocks-1));

  // Total number of crystals
  m_NbTransCrystalsInHead = m_NbTransBlocks * m_NbTransCrystals;
  m_NbTotalTransCrystals = m_NbHeads * m_NbTransBlocks * m_NbTransCrystals;
  m_NbTotalAxialCrystals = m_NbAxialBlocks * m_NbAxialCrystals;
  m_NbTotalCrystals = m_NbLayers * m_NbTotalTransCrystals * m_NbTotalAxialCrystals;
  m_NbTotalTransCrystalsWithGaps = m_NbHeads*(m_NbTransBlocks*m_NbTransCrystals+(m_NbTransBlocks-1)*m_NbTransGapsBetweenBlocks)+m_NbHeads*m_NbGapsBetweenHeads;
  m_NbTotalAxialCrystalsWithGaps = m_NbAxialBlocks*m_NbAxialCrystals + (m_NbAxialBlocks-1)*m_NbAxialGapsBetweenBlocks;
  m_NbTotalCrystalsWithGaps = m_NbLayers * m_NbTotalTransCrystalsWithGaps * m_NbTotalAxialCrystalsWithGaps;

  // Sinogram dimensions affectations
  m_NbElem = f_NbElem;
  m_NbView = f_NbView;
  m_RadialBinSize = f_RadialBinSize;
  m_NbSinoBin = m_NbElem * m_NbView;
  m_NbPlanes = (m_NbAxialCrystals*m_NbAxialBlocks + (m_NbAxialBlocks-1)*m_NbAxialGapsBetweenBlocks)*2 -1;

  // Image default dimensions
  m_FOVSizeX = f_FOVSizeX;
  m_FOVSizeY = f_FOVSizeY;
  m_FOVSizeZ = f_FOVSizeZ;
  m_NbVoxelX = f_NbVoxelX;
  m_NbVoxelY = f_NbVoxelY;
  m_NbVoxelZ = f_NbVoxelZ;
  m_BedOffset = f_BedOffset;

  // Crystal efficiency
  mp_CrystalEfficiency = NULL;

  // Default all pointers
  mp_CornerX1 = NULL;
  mp_CornerY1 = NULL;
  mp_CornerZ1 = NULL;
  mp_CornerX2 = NULL;
  mp_CornerY2 = NULL;
  mp_CornerZ2 = NULL;
  mp_Efficiency = NULL;

  // Verbosity
  if (m_Verbose>=1)
  {
    LogCout ("  --> Total number of crystals: " << m_NbTotalCrystals << "  [" << m_NbTotalTransCrystals << ";" << m_NbTotalAxialCrystals << "]" << endl);
  }

}

// ==========================================================================================================================================
// Destructor
// ==========================================================================================================================================
oScanner::~oScanner()
{
  if (m_Verbose>=1) LogCout ("oScanner::Destructor() -> Destroy" << endl);
  if (mp_CrystalEfficiency) free(mp_CrystalEfficiency);
  if (mp_CornerX1) free(mp_CornerX1);
  if (mp_CornerY1) free(mp_CornerY1);
  if (mp_CornerZ1) free(mp_CornerZ1);
  if (mp_CornerX2) free(mp_CornerX2);
  if (mp_CornerY2) free(mp_CornerY2);
  if (mp_CornerZ2) free(mp_CornerZ2);
  if (mp_Efficiency) free(mp_Efficiency);
}
// ==========================================================================================================================================
// int ProcessCrystalEfficiency()
//   --> This function allocates the crystal efficiency map and reads the crystal efficiency file given in argument. It returns 0 if success
//       and another value otherwise.
// ==========================================================================================================================================
int oScanner::ProcessCrystalEfficiency( bool f_Uniform, bool f_Random, float f_RandomEfficiencyPercentDeviation)
{
  if (m_Verbose>=1) LogCout ("oScanner::ProcessCrystalEfficiency() ..." << endl);

  if (!f_Uniform && !f_Random)
  {
    LogCerr ("***** oScanner::ProcessCrystalEfficiency() -> Must choose between uniform and random efficiencies !" << endl);
    return 1;
  }

  // Allocate memory for the crystal efficiency
  mp_CrystalEfficiency = (float***)malloc(m_NbLayers*sizeof(float**));
  for (int l=0; l<m_NbLayers; l++)
  {
    mp_CrystalEfficiency[l] = (float**)malloc(m_NbTotalAxialCrystals*sizeof(float*));
    for (int a=0; a<m_NbTotalAxialCrystals; a++)
    {
      mp_CrystalEfficiency[l][a] = (float*)calloc(m_NbTotalTransCrystals,sizeof(float));
      if (f_Uniform) for (int t=0; t<m_NbTotalTransCrystals; t++) mp_CrystalEfficiency[l][a][t] = 1.;
    }
  }

  if (f_Random)
  {
    srand(time(NULL)*getpid());
    if (m_Verbose>=2) LogCout ("  --> Shoot random defficiencies at " << f_RandomEfficiencyPercentDeviation << "% around 1" << endl);
    for (int l=0; l<m_NbLayers; l++)
    {
      for (int a=0; a<m_NbTotalAxialCrystals; a++)
      {
        for (int t=0; t<m_NbTotalTransCrystals; t++)
        {
          int random_number1 = rand();
          int random_number2 = rand();
          float random_float1 = ((float)random_number1) / ((float)RAND_MAX);
          float random_float2 = ((float)random_number2) / ((float)RAND_MAX);
          float plus_minus = -1.;
          if (random_float2>0.5) plus_minus = 1.;
          mp_CrystalEfficiency[l][a][t] = 1. + plus_minus * random_float1 * f_RandomEfficiencyPercentDeviation/100.;
        }
      }
    }
  }

  return 0;
}
// ==========================================================================================================================================
// int ProcessDeadBlocks()
//   --> This function set to 0 the crystal efficiencies for crystals contained within the given list of dead blocks.
//       It returns 0 if success and 1 otherwise.
// ==========================================================================================================================================
int oScanner::ProcessDeadBlocks(int f_NbDeadBlocks, int* f_DeadBlockHeadIndex, int* f_DeadBlockAxialIndex, int* f_DeadBlockTransIndex)
{
  // If no dead block, we exit
  if (f_NbDeadBlocks==0) return 0;

  // Verbosity
  if (m_Verbose>=1) LogCout ("oScanner::ProcessDeadBlocks() ..." << endl);

  // Process the block list
  for (int b=0; b<f_NbDeadBlocks; b++)
  {
    // Check if block exists
    if (f_DeadBlockHeadIndex[b]<0 || f_DeadBlockHeadIndex[b]>=m_NbHeads)
    {
      LogCerr ("***** oScanner::ProcessDeadBlocks() -> Head index (" << f_DeadBlockHeadIndex[b] << ") is out of range for dead block n°" << b << " !" << endl);
      return 1;
    }
    if (f_DeadBlockAxialIndex[b]<0 || f_DeadBlockAxialIndex[b]>=m_NbAxialBlocks)
    {
      LogCerr ("***** oScanner::ProcessDeadBlocks() -> Axial block index (" << f_DeadBlockAxialIndex[b] << ") is out of range for dead block n°" << b << " !" << endl);
      return 1;
    }
    if (f_DeadBlockTransIndex[b]<0 || f_DeadBlockTransIndex[b]>=m_NbTransBlocks)
    {
      LogCerr ("***** oScanner::ProcessDeadBlocks() -> Transaxial block index (" << f_DeadBlockTransIndex[b] << ") is out of range for dead block n°" << b << " !" << endl);
      return 1;
    }
    // Verbosity
    if (m_Verbose>=2) LogCout ("oScanner::ProcessDeadBlocks() -> Dead block n° " << b << " for head " << f_DeadBlockHeadIndex[b] << " at position [" << f_DeadBlockAxialIndex[b] << ";" << f_DeadBlockTransIndex[b] << "]" << endl);
    // Process
    for (int l=0; l<m_NbLayers; l++)
    {
      int axial_begin_index = f_DeadBlockAxialIndex[b]*m_NbAxialCrystals;
      for (int a=axial_begin_index; a<axial_begin_index+m_NbAxialCrystals; a++)
      {
        int trans_begin_index = f_DeadBlockHeadIndex[b]*m_NbTransBlocks*m_NbTransCrystals + f_DeadBlockTransIndex[b]*m_NbTransCrystals;
        for (int t=trans_begin_index; t<trans_begin_index+m_NbTransCrystals; t++)
        {
          mp_CrystalEfficiency[l][a][t] = 0.;
        }
      }
    }
  }

  return 0;
}
// ==========================================================================================================================================
// int WriteCrystalEfficiency()
//   --> This functions write the crystal efficiency in a file organized as layers, then axial crystals, and finally transaxial crystals.
//       It returns 0 if success and other value otherwise.
// ==========================================================================================================================================
int oScanner::WriteCrystalEfficiency( const string& f_FileBaseOut )
{
  if (m_Verbose>=1) LogCout ("oScanner::WriteCrystalEfficiency() ..." << endl);
  // Open file
  string f_FileName = f_FileBaseOut+"/"+f_FileBaseOut+"_eff.dat";
  FILE* fout = fopen(f_FileName.c_str(),"wb");
  if (fout==NULL)
  {
    LogCerr ("***** oScanner::WriteCrystalEfficiency() -> Failed to open output file '" << f_FileName << "' !" << endl);
    return 1;
  }
  // Write data
  int nb_data_written = 0;
  for (int l=0; l<m_NbLayers; l++) for (int a=0; a<m_NbTotalAxialCrystals; a++) for (int t=0; t<m_NbTotalTransCrystals; t++)
  {
    nb_data_written += fwrite(&mp_CrystalEfficiency[l][a][t],sizeof(float),1,fout);
  }
  // Close file
  fclose(fout);
  // Check if all data were written
  int all_data = m_NbLayers*m_NbTotalAxialCrystals*m_NbTotalTransCrystals;
  if (nb_data_written != all_data)
  {
    LogCerr ("***** oScanner::WriteCrystalEfficiency() -> Failed to write all data (" << all_data << ") in file '" << f_FileName << "' (" << nb_data_written << " written) !" << endl);
    return 1;
  }
  // Verbosity
  if (m_Verbose>=1)
  {
    LogCout ("  --> Written in file '" << f_FileName << "'." << endl);
    LogCout ("  --> Dimensions: [" << m_NbTotalTransCrystals << ";" << m_NbTotalAxialCrystals << ";" << m_NbLayers << "] in float 32 bits raw format." << endl);
  }
  return 0;
}
// ==========================================================================================================================================
// int WriteCrystalMap()
//
// ==========================================================================================================================================
int oScanner::WriteCrystalMap( const string& f_FileBaseOut, bool f_Castor )
{
  // Verbosity
  if (m_Verbose>=1) LogCout ("oScanner::WriteCrystalMap() ..." << endl);
  // Open file
  string f_FileName = "";
  if (!f_Castor) f_FileName = f_FileBaseOut+"/"+f_FileBaseOut+".ecm";
  else f_FileName = f_FileBaseOut+"/"+f_FileBaseOut+".lut";
  FILE* fout = fopen(f_FileName.c_str(),"wb");
  if (fout==NULL)
  {
    LogCerr ("***** oScanner::WriteCrystalMap() -> Failed to open output file '" << f_FileName << "' !" << endl);
    return 1;
  }
  // Loop on all crystals (generic for all scanners)
  int nb_data_written = 0;
  for (int l=0; l<m_NbLayers; l++)
  {
    for (int ab=0; ab<m_NbAxialBlocks; ab++)
    {
      for (int ac=0; ac<m_NbAxialCrystals; ac++)
      {
        for (int h=0; h<m_NbHeads; h++)
        {
          for (int tb=0; tb<m_NbTransBlocks; tb++)
          {
            for (int tc=0; tc<m_NbTransCrystals; tc++)
            {
              // Deal first with coordinates
              float corner_1X, corner_1Y, corner_1Z, corner_2X, corner_2Y, corner_2Z;
              ComputeDiagCoordinates(l,ab,ac,h,tb,tc,&corner_1X,&corner_1Y,&corner_1Z,&corner_2X,&corner_2Y,&corner_2Z,f_Castor);
              if (!f_Castor)
              {
                nb_data_written += fwrite(&corner_1X,sizeof(float),1,fout);
                nb_data_written += fwrite(&corner_1Y,sizeof(float),1,fout);
                nb_data_written += fwrite(&corner_1Z,sizeof(float),1,fout);
                nb_data_written += fwrite(&corner_2X,sizeof(float),1,fout);
                nb_data_written += fwrite(&corner_2Y,sizeof(float),1,fout);
                nb_data_written += fwrite(&corner_2Z,sizeof(float),1,fout);
              }
              else
              {
                float posX = (corner_1X+corner_2X) / 2.;
                float posY = (corner_1Y+corner_2Y) / 2.;
                float posZ = (corner_1Z+corner_2Z) / 2.;
                nb_data_written += fwrite(&posX,sizeof(float),1,fout);
                nb_data_written += fwrite(&posY,sizeof(float),1,fout);
                nb_data_written += fwrite(&posZ,sizeof(float),1,fout);
                float orientationX, orientationY, orientationZ;
                ComputeOrientationVector(l,ab,ac,h,tb,tc,&orientationX,&orientationY,&orientationZ,f_Castor);
                nb_data_written += fwrite(&orientationX,sizeof(float),1,fout);
                nb_data_written += fwrite(&orientationY,sizeof(float),1,fout);
                nb_data_written += fwrite(&orientationZ,sizeof(float),1,fout);
              }
              // Then with efficiencies
              int axial_index = ab*m_NbAxialCrystals+ac;
              int trans_index = h*m_NbTransBlocks*m_NbTransCrystals+tb*m_NbTransCrystals+tc;
              float efficiency = mp_CrystalEfficiency[l][axial_index][trans_index];
              if (!f_Castor) nb_data_written += fwrite(&efficiency,sizeof(float),1,fout);
              // Verbose
              if (m_Verbose>=3) LogCout ("  --> Layer: " << l << " | Axial block: " << ab << " | Axial crystal: " << ac << " | Head: " << h << " | Trans block: " << tb << " | Trans crystal: " << tc << " | Axial index: " << axial_index << " | Trans index: " << trans_index << endl << "      Corner1: [" << corner_1X << ";" << corner_1Y << ";" << corner_1Z << "] | Corner2: [" << corner_2X << ";" << corner_2Y << ";" << corner_2Z << "] | Efficiency: " << efficiency << endl);
            }
          }
        }
      }
    }
  }
  // Close file
  fclose(fout);
  // Check if all data were written
  int all_data = 0;
  if (!f_Castor) all_data = m_NbLayers*m_NbTotalAxialCrystals*m_NbTotalTransCrystals*7;
  else all_data = m_NbLayers*m_NbTotalAxialCrystals*m_NbTotalTransCrystals*6;
  if (nb_data_written != all_data)
  {
    LogCerr ("***** oScanner::WriteCrystalMap() -> Failed to write all data (" << all_data << " in file '" << f_FileName << "' (" << nb_data_written << " written) !" << endl);
    return 1;
  }
  // Verbosity
  if (m_Verbose>=1) LogCout ("  --> Written in file '" << f_FileName << "'." << endl);

  // ------------------------------------------------------------
  // For castor we need to write a header
  // ------------------------------------------------------------

  if (f_Castor)
  {
    // Open file
    string f_HeaderName = f_FileBaseOut+"/"+f_FileBaseOut+".hscan";
    ofstream fhdr(f_HeaderName.c_str());
    if (!fhdr)
    {
      LogCerr ("***** oScanner::WriteCrystalMap() -> Failed to open output header file '" << f_HeaderName << "' !" << endl);
      return 1;
    }
    // Write infos
    if (m_ScannerModel==SCANNER_BIOGRAPH) fhdr << "scanner name: PET_SIEMENS_BIOGRAPH6_TRUEPOINT_TRUEV" << endl;
    else if (m_ScannerModel==SCANNER_FOCUS) fhdr << "scanner name: PET_SIEMENS_FOCUS" << endl;
    fhdr << "modality: PET" << endl;
    fhdr << "scanner radius: " << m_Radius + m_MeanDOI << endl;
    fhdr << "number of rings in scanner: " << m_NbAxialBlocks * m_NbAxialCrystals << endl;
    fhdr << "number of elements: " << m_NbTotalTransCrystals * m_NbTotalAxialCrystals * m_NbLayers << endl;
    fhdr << "number of layers: " << m_NbLayers << endl;
    fhdr << "number of crystals in layer: " << m_NbTotalTransCrystals * m_NbTotalAxialCrystals * m_NbLayers << endl;
    fhdr << "layers size depth: " << m_DepthCrystalSize << endl;
    fhdr << "layers size trans: " << m_TransCrystalSize << endl;
    fhdr << "layers size axial: " << m_AxialCrystalSize << endl;
    fhdr << "voxels number trans: 256" << endl;
    fhdr << "voxels number axial: 128" << endl;
    fhdr << "field of view trans: 768." << endl;
    fhdr << "field of view axial: 256." << endl;
    fhdr << "mean depth of interaction: -1" << endl;
    fhdr << "min angle difference: 90" << endl;
    fhdr << "description: " << m_ScannerModel << endl;
    // Close file
    fhdr.close();
  }

  // End
  return 0;
}
// ==========================================================================================================================================
// int ComputeDiagCoordinates()
//   --> This functions computes the coordinates of two corners of the given crystal. The two corners must be diagonal opposites and the
//       first one has the smallest Z coordinate and the first transaxial position clock-wise.
//       The resulting 3D coordinates of the two corners are transferred in the arguments, while the function returns 0 if crystal exists
//       and othre value otherwise.
// ==========================================================================================================================================
int oScanner::ComputeDiagCoordinates( int f_LayerIndex, int f_AxialBlockIndex, int f_AxialCrystalIndex,
                                      int f_HeadIndex, int f_TransBlockIndex, int f_TransCrystalIndex,
                                      float* fp_Corner1X, float* fp_Corner1Y, float* fp_Corner1Z,
                                      float* fp_Corner2X, float* fp_Corner2Y, float* fp_Corner2Z, bool f_Castor )
{
  // Coordinates calculated for an arbitrary block without any rotation (cylindrical)
  float Xr = ((float)f_TransBlockIndex)   * (m_TransBlockSize+m_TransBlockGap)
            + ((float)f_TransCrystalIndex) * (m_TransCrystalSize+m_TransCrystalGap)
            - m_TransHeadSize/2. + m_TransCrystalSize/2.;
  // OLD VERSION float Yr = m_Radius + m_MeanDOI + ((float)f_LayerIndex)*m_DepthCrystalSize/((float)m_NbLayers);
  float Yr = 0.;
  if (f_Castor) Yr = m_Radius + m_DepthCrystalSize/2.;
  else Yr = m_Radius + m_MeanDOI;

  float Zr = ((float)f_AxialBlockIndex)   * (m_AxialBlockSize+m_AxialBlockGap)
            + ((float)f_AxialCrystalIndex) * (m_AxialCrystalSize+m_AxialCrystalGap);
  if (f_Castor) Zr += m_AxialCrystalSize/2. - m_AxialScannerSize/2.;

  // Global coordinates in the cartesian grid
  float head_angle = 0.;
  if ( m_ScannerModel == SCANNER_BIOGRAPH || m_ScannerModel == SCANNER_INVEON ||
       m_ScannerModel == SCANNER_FOCUS || m_ScannerModel == SCANNER_CASTOR ||
       m_ScannerModel == SCANNER_BIOGRAPH2D )
    head_angle = 2.*M_PI*((float)f_HeadIndex)/((float)m_NbHeads) + 1.*M_PI / ((float)m_NbHeads);
  else head_angle = 2.*M_PI*((float)f_HeadIndex)/((float)m_NbHeads);

  float pos_x =  Xr * cos(head_angle) + Yr * sin(head_angle);
  float pos_y = -Xr * sin(head_angle) + Yr * cos(head_angle);
  float pos_z = Zr;

  // Determine the 3D deltas to find the opposite corners' coordinates
  float alpha = M_PI - head_angle;
  float delta_x = m_TransCrystalSize * cos(alpha) / 2.;
  float delta_y = m_TransCrystalSize * sin(alpha) / 2.;
  float delta_z = m_AxialCrystalSize / 2.;

  // Results
  *fp_Corner1X = pos_x - delta_x;
  *fp_Corner1Y = pos_y - delta_y;
  *fp_Corner1Z = pos_z - delta_z;
  *fp_Corner2X = pos_x + delta_x;
  *fp_Corner2Y = pos_y + delta_y;
  *fp_Corner2Z = pos_z + delta_z;
  // Note: the zero along Z axis is at half the first crystal

  // End
  return 0;
}
// ==========================================================================================================================================
// int ComputeOrientationVector()
//   --> This functions computes orientation of a crystal, ie the normalized vector that follow the direction of the crystal's axis, from
//       its entrance to its back.
// ==========================================================================================================================================
int oScanner::ComputeOrientationVector( int f_LayerIndex, int f_AxialBlockIndex, int f_AxialCrystalIndex,
                                        int f_HeadIndex, int f_TransBlockIndex, int f_TransCrystalIndex,
                                        float* fp_OrientationX, float* fp_OrientationY, float* fp_OrientationZ,
                                        bool f_Castor )
{
  // Coordinates calculated for an arbitrary block without any rotation (cylindrical)
  float Xr = ((float)f_TransBlockIndex)   * (m_TransBlockSize+m_TransBlockGap)
            + ((float)f_TransCrystalIndex) * (m_TransCrystalSize+m_TransCrystalGap)
            - m_TransHeadSize/2. + m_TransCrystalSize/2.;
/*OLD VERSION
  float Yr1 = m_Radius + 1. + ((float)f_LayerIndex)*m_DepthCrystalSize/((float)m_NbLayers);
  float Yr2 = m_Radius - 1. + ((float)f_LayerIndex)*m_DepthCrystalSize/((float)m_NbLayers);
*/
  float Yr1 = m_Radius + 1.;
  float Yr2 = m_Radius - 1.;
  if (f_Castor)
  {
    Yr1 += m_DepthCrystalSize/2.;
    Yr2 += m_DepthCrystalSize/2.;
  }
  else
  {
    Yr1 += m_MeanDOI;
    Yr2 += m_MeanDOI;
  }
  float Zr = ((float)f_AxialBlockIndex)   * (m_AxialBlockSize+m_AxialBlockGap)
            + ((float)f_AxialCrystalIndex) * (m_AxialCrystalSize+m_AxialCrystalGap);
  if (f_Castor) Zr += m_AxialCrystalSize/2. - m_AxialScannerSize/2.;

  // Global coordinates in the cartesian grid
  float head_angle = 0.;
  if ( m_ScannerModel == SCANNER_BIOGRAPH || m_ScannerModel == SCANNER_INVEON ||
       m_ScannerModel == SCANNER_FOCUS || m_ScannerModel == SCANNER_CASTOR ||
       m_ScannerModel == SCANNER_BIOGRAPH2D )
    head_angle = 2.*M_PI*((float)f_HeadIndex)/((float)m_NbHeads) + 1.*M_PI / ((float)m_NbHeads);
  else head_angle = 2.*M_PI*((float)f_HeadIndex)/((float)m_NbHeads);

  float pos_x1 =  Xr * cos(head_angle) + Yr1 * sin(head_angle);
  float pos_y1 = -Xr * sin(head_angle) + Yr1 * cos(head_angle);
  float pos_x2 =  Xr * cos(head_angle) + Yr2 * sin(head_angle);
  float pos_y2 = -Xr * sin(head_angle) + Yr2 * cos(head_angle);
  float pos_z = Zr;

  // Compute the orientation vector from points 1 and 2
  *fp_OrientationZ = 0.;
  // We subtract point2 from point 1 to go from crystal entrance to crystal back
  float vectorX = pos_x1 - pos_x2;
  float vectorY = pos_y1 - pos_y2;
  float length = sqrt(vectorX*vectorX + vectorY*vectorY);
  *fp_OrientationX = vectorX/length;
  *fp_OrientationY = vectorY/length;

  // End
  return 0;
}
// ==========================================================================================================================================
// int ReadCrystalMap()
//   --> This function read the crystal map from the given file.
// ==========================================================================================================================================
int oScanner::ReadCrystalMap( const string& f_MapFileName )
{
  m_MapFileName = f_MapFileName;
  if (m_Verbose>=1) LogCout ("oScanner::ReadCrystalMap() -> Read esteban crystal map from file '" << m_MapFileName << "'" << endl);

  // Allocations
  mp_CornerX1 = (float*)malloc(m_NbTotalCrystals*sizeof(float));
  mp_CornerY1 = (float*)malloc(m_NbTotalCrystals*sizeof(float));
  mp_CornerZ1 = (float*)malloc(m_NbTotalCrystals*sizeof(float));
  mp_CornerX2 = (float*)malloc(m_NbTotalCrystals*sizeof(float));
  mp_CornerY2 = (float*)malloc(m_NbTotalCrystals*sizeof(float));
  mp_CornerZ2 = (float*)malloc(m_NbTotalCrystals*sizeof(float));
  mp_Efficiency = (float*)malloc(m_NbTotalCrystals*sizeof(float));

  // Open the file
  FILE* fmap = fopen(m_MapFileName.c_str(),"rb");
  if (fmap==NULL)
  {
    LogCerr ("***** oScanner::ReadCrystalMap() -> Input map file '" << m_MapFileName << "' is missing or corrupted !" << endl);
    return 1;
  }

  // Read file
  int nb_data_read = 0;
  for (int l=0, i=0; l<m_NbLayers; l++) for (int r=0; r<m_NbTotalAxialCrystals; r++) for (int c=0; c<m_NbTotalTransCrystals; c++, i++)
  {
    nb_data_read += fread(&mp_CornerX1[i],sizeof(float),1,fmap);
    nb_data_read += fread(&mp_CornerY1[i],sizeof(float),1,fmap);
    nb_data_read += fread(&mp_CornerZ1[i],sizeof(float),1,fmap);
    nb_data_read += fread(&mp_CornerX2[i],sizeof(float),1,fmap);
    nb_data_read += fread(&mp_CornerY2[i],sizeof(float),1,fmap);
    nb_data_read += fread(&mp_CornerZ2[i],sizeof(float),1,fmap);
    nb_data_read += fread(&mp_Efficiency[i],sizeof(float),1,fmap);
  }

  // Close file
  fclose(fmap);

  // Check
  int nb_data_to_be_read = 7*m_NbLayers*m_NbTotalAxialCrystals*m_NbTotalTransCrystals;
  if (nb_data_read!=nb_data_to_be_read)
  {
    LogCerr ("***** oScanner::ReadCrystalMap() -> Failed to read all data (" << nb_data_to_be_read << ") in input map '" << m_MapFileName << "' (" << nb_data_read << " read) !" << endl);
    return 1;
  }

  // Ending
  return 0;
}
// ==========================================================================================================================================
// int ComputeUniformCrystalMap()
//   --> This function allocate and compute a crystal map with uniform efficiencies
// ==========================================================================================================================================
int oScanner::ComputeUniformCrystalMap(bool f_Castor)
{
  if (m_Verbose>=1) LogCout ("oScanner::ComputeUniformCrystalMap() -> Start computation" << endl);

  // Allocations
  mp_CornerX1 = (float*)malloc(m_NbTotalCrystals*sizeof(float));
  mp_CornerY1 = (float*)malloc(m_NbTotalCrystals*sizeof(float));
  mp_CornerZ1 = (float*)malloc(m_NbTotalCrystals*sizeof(float));
  mp_CornerX2 = (float*)malloc(m_NbTotalCrystals*sizeof(float));
  mp_CornerY2 = (float*)malloc(m_NbTotalCrystals*sizeof(float));
  mp_CornerZ2 = (float*)malloc(m_NbTotalCrystals*sizeof(float));
  mp_Efficiency = (float*)malloc(m_NbTotalCrystals*sizeof(float));

  // Loop on all crystals (generic for all scanners)
  int nb_data_written = 0;
  for (int l=0, index=0; l<m_NbLayers; l++)
  {
    for (int ab=0; ab<m_NbAxialBlocks; ab++)
    {
      for (int ac=0; ac<m_NbAxialCrystals; ac++)
      {
        for (int h=0; h<m_NbHeads; h++)
        {
          for (int tb=0; tb<m_NbTransBlocks; tb++)
          {
            for (int tc=0; tc<m_NbTransCrystals; tc++, index++)
            {
              // Deal first with coordinates
              float corner_1X, corner_1Y, corner_1Z, corner_2X, corner_2Y, corner_2Z;
              ComputeDiagCoordinates(l,ab,ac,h,tb,tc,&corner_1X,&corner_1Y,&corner_1Z,&corner_2X,&corner_2Y,&corner_2Z,f_Castor);
              // Affect tables
              mp_CornerX1[index] = corner_1X;
              mp_CornerY1[index] = corner_1Y;
              mp_CornerZ1[index] = corner_1Z;
              mp_CornerX2[index] = corner_2X;
              mp_CornerY2[index] = corner_2Y;
              mp_CornerZ2[index] = corner_2Z;
              // Uniform efficiency
              mp_Efficiency[index] = 1.;
            }
          }
        }
      }
    }
  }

  // Ending
  return 0;
}

// ==========================================================================================================================================
// int ComputeNbSinoFromSpanAndMaxRingDiff()
//   --> This function can be used to compute the number of sinograms from a given span and maximum ring difference
//       It returns -1 if there is a inconsistency betwenn the given parameters.
// ==========================================================================================================================================
int oScanner::ComputeNbSinoFromSpanAndMaxRingDiff(int f_Span, int f_MaxRingDiff)
{
  // Get the number of rings
  int nb_rings = m_NbTotalAxialCrystals;

  // Check if the max ring difference is consistent with the span (and calculate number of segments at the same time)
  int check_span = f_MaxRingDiff;
  bool correct = false;
  int nb_absolute_segments = 0;
  while (check_span >= 0)
  {
    if (check_span==f_Span/2)
    {
      correct = true;
      break;
    }
    check_span -= f_Span;
    nb_absolute_segments++;
  }
/*
  if (!correct)
  {
    LogCerr ("***** oScanner::ComputeNbSinoFromSpanAndMaxRingDiff() -> Inconsistence between the span (" << f_Span << ") and the maximum ring difference (" << f_MaxRingDiff << ") !" << endl);
    return -1;
  }
*/
  int nb_segments = 2*nb_absolute_segments+1;
  // Main loop on segments and increasing sino
  int segment_number = 0;
  int sino_bin = 0;
  int nb_ring_pairs_by_sino = 0;
  for (int s=0; s<nb_segments; s++)
  {
    // Calculate the min and max allowed ring differences for that segment
    int min_allowed_ring_difference = -f_Span/2 + segment_number*f_Span;
    int max_allowed_ring_difference = f_Span/2 + segment_number*f_Span;

    // Do a loop on all possible ring indices sum
    for (int sum=0; sum<m_NbTotalAxialCrystals*2; sum++)
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
  // Return
  return sino_bin;
}

