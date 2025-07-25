#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <ctime>
#ifdef OMP_MODE
#include <omp.h>
#endif
#include "oSinoModeCreation.hh"
#include "oTableBiograph.hh"
#include "oTableBiograph2D.hh"
#include "oTableMMR2D.hh"
#include "oTableVision600.hh"
#include "oTableHRRT.hh"
#include "oTableHRplus.hh"
#include "oTableCastor.hh"
#include "oScanner.hh"
#include "oMiscellaneous.hh"
#include "oOutputManager.hh"
using namespace std;

// ==========================================================================================================================================
// Constructor
// ==========================================================================================================================================
oSinoModeCreation::oSinoModeCreation(const string& f_ScannerName, int f_Weighting, bool f_SiddonDidier, int f_Seed, int f_Verbose)
{
  if (f_Verbose) LogCout ("oSinoModeCreation::Constructor() -> Construct" << endl);

  // Weighting
  m_Weighting = f_Weighting;

  // Mode 2D
  m_Mode2D = false;

  // Zeroing counters
  m_NbPrompts  = 0;
  m_NbDelays   = 0;
  m_NbNetTrues = 0;
  m_NbScatters = 0.;
  m_NbRandoms  = 0.;
  m_ScatterFraction = -1.;

  // Initializing the random generator
  m_Seed = Misc_InitRandomGenerator(f_Seed);

  // Init time
  m_StartTime = 0.;
  m_Duration = 0.;

  // Corrections
  m_InitBool = false;
  m_TrueBool = false;
  mp_TrueSino = NULL;
  m_PromptBool = false;
  mp_PromptSino = NULL;
  m_RandCorr = false;
  mp_RandSino = NULL;
  m_ScatCorr = false;
  mp_ScatSino = NULL;
  m_NormCorr = false;
  mp_NormSino = NULL;
  m_FloatBool = false;
  mp_FloatSino = NULL;
  m_AttnCorr = false;
  mp_AttnSino = NULL;
  mp_AttnUMap = NULL;
  m_AttenuationMode = -1;
  m_AttnHalfVoxelShift = true;
  m_SiddonDidier = f_SiddonDidier;
  m_ReadSpan = -1;
  m_ReadNbSino = -1;
  m_ReadNbView = -1;
  m_ReadNbElem = -1;
  m_ReadMaxRingDiff = -1;
  m_ECF = 1.;
  m_AxisOrientationX = 0.;
  m_AxisOrientationY = 0.;

  // Dead time correction factor: we initialize it at 1 and it will be included into the norm factor (only needed for inveon, focus and hrrt)
  m_DeadTimeCorrectionFactor = 1.;

  // Patient orientation is default to unknown (do nothing)
  m_Orientation = "UNKNOWN";

  // Tables
  mp_BiographTable = NULL;
  mp_CastorTable = NULL;
  mp_BiographTable2D = NULL;
  mp_MMRTable2D = NULL;
  mp_HRplusTable = NULL;
  mp_HrrtTable = NULL;
  mp_InveonTable = NULL;
  mp_FocusTable = NULL;
  mp_Vision600Table = NULL;

  // Isotope
  m_Isotope = "INF";

  // Frame number
  m_Frame = 0;

  // Verbose
  m_Verbose = f_Verbose;

  // Initialize the scanner's related parameters
  mp_Scanner = new oScanner(f_ScannerName, f_Verbose);
}
// ==========================================================================================================================================
// Destructor
// ==========================================================================================================================================
oSinoModeCreation::~oSinoModeCreation()
{
  if (m_Verbose>=1) LogCout ("oSinoModeCreation::Destructor() -> Destroy" << endl);

  // CHECK OTHER FREE OR DELETE TO DO

  // Print out counters
  if (m_Verbose>=1)
  {
  }

  // Delete objects
  if (mp_BiographTable) delete mp_BiographTable;
  if (mp_CastorTable) delete mp_CastorTable;
  if (mp_BiographTable2D) delete mp_BiographTable2D;
  if (mp_MMRTable2D) delete mp_MMRTable2D;
  if (mp_Vision600Table) delete mp_Vision600Table;
  if (mp_HRplusTable) delete mp_HRplusTable;
  if (mp_HrrtTable) delete mp_HrrtTable;
  if (mp_InveonTable) delete mp_InveonTable;
  if (mp_FocusTable) delete mp_FocusTable;
  if (mp_Scanner) delete mp_Scanner;

}
// ==========================================================================================================================================
// int InitConversionTable()
//    --> This function initialize the conversion table. It returns 0 upon success and other value otherwise.
// ==========================================================================================================================================
int oSinoModeCreation::InitConversionTable()
{
  if (m_Verbose>=1) LogCout ("oSinoModeCreation::InitConversionTable() -> Initialize tables" << endl);

  // Security: if net trues have not been read, stop
  if (!m_InitBool)
  {
    LogCerr ("***** oSinoModeCreation::InitConversionTable() -> Cannot initialize conversion table before net trues !" << endl);
    return 1;
  }

  // Get the scanner model
  int scanner_model = mp_Scanner->GetScannerModel();

  // Switch on different scanner
  if (scanner_model==SCANNER_BIOGRAPH)
  {
    // Create the object
    mp_BiographTable = new oTableBiograph(mp_Scanner, m_Verbose);
    // Read the span table
    if (mp_BiographTable->ReadSpanTable())
    {
      LogCerr ("***** oSinoModeCreation::InitConversionTable() -> An error occured when reading the span table !" << endl);
      return 1;
    }
    // Compute ring tables
    if (mp_BiographTable->ComputeRingTables())
    {
      LogCerr ("***** oSinoModeCreation::InitConversionTable() -> An error occured when computing the ring conversion tables !" << endl);
      return 1;
    }
  }
  else if (scanner_model==SCANNER_CASTOR)
  {
    // Create the object
    mp_CastorTable = new oTableCastor(mp_Scanner, m_Verbose);
    // Read the span table
    if (mp_CastorTable->ReadSpanTable())
    {
      LogCerr ("***** oSinoModeCreation::InitConversionTable() -> An error occured when reading the span table !" << endl);
      return 1;
    }
    // Compute ring tables
    if (mp_CastorTable->ComputeRingTables())
    {
      LogCerr ("***** oSinoModeCreation::InitConversionTable() -> An error occured when computing the ring conversion tables !" << endl);
      return 1;
    }
  }
  else if (scanner_model==SCANNER_BIOGRAPH2D)
  {
    // Create the object
    mp_BiographTable2D = new oTableBiograph2D(mp_Scanner, m_Verbose);
    // Compute ring tables
    if (mp_BiographTable2D->ComputeRingTables())
    {
      LogCerr ("***** oSinoModeCreation::InitConversionTable() -> An error occured when computing the ring conversion tables !" << endl);
      return 1;
    }
  }
  else if (scanner_model==SCANNER_MMR2D)
  {
    // Create the object
    mp_MMRTable2D = new oTableMMR2D(mp_Scanner, m_Verbose);
    // Compute ring tables
    if (mp_MMRTable2D->ComputeRingTables())
    {
      LogCerr ("***** oSinoModeCreation::InitConversionTable() -> An error occured when computing the ring conversion tables !" << endl);
      return 1;
    }
  }
  else if (scanner_model==SCANNER_VISION600)
  {
    // Create the object
    mp_Vision600Table = new oTableVision600(mp_Scanner, m_Verbose);
    // Compute ring tables
    // Read the span table
    if (mp_Vision600Table->ReadSpanTable())
    {
      LogCerr ("***** oSinoModeCreation::InitConversionTable() -> An error occured when reading the span table !" << endl);
      return 1;
    }
    // Compute ring tables
    if (mp_Vision600Table->ComputeRingTables())
    {
      LogCerr ("***** oSinoModeCreation::InitConversionTable() -> An error occured when computing the ring conversion tables !" << endl);
      return 1;
    }
  }
  else if (scanner_model==SCANNER_HRPLUS)
  {
    // Create object
    mp_HRplusTable = new oTableHRplus( m_ReadNbElem, m_ReadNbView, m_ReadNbSino, m_ReadMaxRingDiff, m_ReadSpan, mp_Scanner, m_Verbose );
    // Compute span tables
    if (mp_HRplusTable->ComputeSpanTables())
    {
      LogCerr ("***** oSinoModeCreation::InitConversionTable() -> An error occured when computing the span conversion tables !" << endl);
      return 1;
    }
    // Compute ring tables
    if (mp_HRplusTable->ComputeRingTables())
    {
      LogCerr ("***** oSinoModeCreation::InitConversionTable() -> An error occured when computing the ring conversion tables !" << endl);
      return 1;
    }
  }
  else if (scanner_model==SCANNER_HRRT)
  {
    // Create the object
    mp_HrrtTable = new oTableHRRT(mp_Scanner, m_Verbose);
    // Read the ring table
    if (mp_HrrtTable->ReadLUTTable())
    {
      LogCerr ("***** oSinoModeCreation::InitConversionTable() -> An error occured when reading the ring table !" << endl);
      return 1;
    }
    // Compute span tables
    if (mp_HrrtTable->ComputeSpanTables())
    {
      LogCerr ("***** oSinoModeCreation::InitConversionTable() -> An error occured when computing the span conversion tables !" << endl);
      return 1;
    }
    // Compute final conversion tables (user friendly !)
    if (mp_HrrtTable->ComputeFinalTables())
    {
      LogCerr ("***** oSinoModeCreation::InitConversionTable() -> An error occured when computing the final conversion tables !" << endl);
      return 1;
    }
  }
  else if (scanner_model==SCANNER_INVEON)
  {
    // Create object
    mp_InveonTable = new oTableInveon( m_ReadNbElem, m_ReadNbView, m_ReadNbSino, m_ReadMaxRingDiff, m_ReadSpan, mp_Scanner, m_Verbose );
    // Compute span tables
    if (mp_InveonTable->ComputeSpanTables())
    {
      LogCerr ("***** oSinoModeCreation::InitConversionTable() -> An error occured when computing the span conversion tables !" << endl);
      return 1;
    }
    // Compute ring tables
    if (mp_InveonTable->ComputeRingTables())
    {
      LogCerr ("***** oSinoModeCreation::InitConversionTable() -> An error occured when computing the ring conversion tables !" << endl);
      return 1;
    }
  }
  else if (scanner_model==SCANNER_FOCUS)
  {
    // Create object
    mp_FocusTable = new oTableFocus( m_ReadNbElem, m_ReadNbView, m_ReadNbSino, m_ReadMaxRingDiff, m_ReadSpan, mp_Scanner, m_Verbose );
    // Compute span tables
    if (mp_FocusTable->ComputeSpanTables())
    {
      LogCerr ("***** oSinoModeCreation::InitConversionTable() -> An error occured when computing the span conversion tables !" << endl);
      return 1;
    }
    // Compute ring tables
    if (mp_FocusTable->ComputeRingTables())
    {
      LogCerr ("***** oSinoModeCreation::InitConversionTable() -> An error occured when computing the ring conversion tables !" << endl);
      return 1;
    }
  }
  else
  {
    LogCerr ("***** oSinoModeCreation::InitConversionTable() -> Unknown scanner model !" << endl);
    return 1;
  }

  // Ending
  return 0;
}
// ==========================================================================================================================================
// int InitCrystalMap()
//    --> This function initialize the crystal map. It returns 0 upon success and other value otherwise.
// ==========================================================================================================================================
int oSinoModeCreation::InitCrystalMap(const string& f_FileMap, bool f_Castor)
{
  // Read the map if file is provided
  if (f_FileMap!="")
  {
    if (mp_Scanner->ReadCrystalMap( f_FileMap ))
    {
      LogCerr ("***** oSinoModeCreation::InitCrystalMap() -> An error occured when reading the crystal map !" << endl);
      return 1;
    }
  }
  // Otherwise, compute it assuming an uniform crystal efficiency
  else
  {
    if (mp_Scanner->ComputeUniformCrystalMap(false))
    {
      LogCerr ("***** oSinoModeCreation::InitCrystalMap() -> An error occured while computing the crystal map !" << endl);
      return 1;
    }
  }

  // Ending
  return 0;
}
// ==========================================================================================================================================
// int ProcessListMode()
//    --> This function is the general one that process the input sinograms and create the output list-mode file.
// ==========================================================================================================================================
int oSinoModeCreation::ProcessListMode(const string& f_FileBaseOut, bool f_FillEqualLORs, int f_Threads, bool f_Castor, PRECISION f_LimitedAnglePercent)
{
  if (m_Verbose>=1) LogCout ("oSinoModeCreation::ProcessListMode() -> Begin processing ..." << endl << flush);

  // ======================================================================================================================
  // Limited angle tomography
  int keptNbView  = m_ReadNbView;
  int shiftView = 0;
  if (f_LimitedAnglePercent!=100.)
  {
    // Check validity of the provided value
    if (f_LimitedAnglePercent<=0. || f_LimitedAnglePercent>100.)
    {
      LogCerr ("***** oSinoModeCreation::ProcessListMode() -> The provided percentage of azymuthal angles must be in the range ]0:100] !" << endl);
      return 1;
    }
    // Verbose
    if (m_Verbose>=2) LogCout ("  --> Doing limited angle tomography" << endl);
    // Compute the end angle from the percentage, making it even
//    keptNbView = (int)( ((PRECISION)m_ReadNbView)*f_LimitedAnglePercent/100. );
    PRECISION keptNbViewFloat = ((PRECISION)m_ReadNbView)*f_LimitedAnglePercent/100.;
    PRECISION rest = keptNbViewFloat - ((PRECISION)( (((int)(keptNbViewFloat))/2) ))*2.;
    if (rest<1.) keptNbView = ((int)(keptNbViewFloat));
    else keptNbView = ((int)(keptNbViewFloat)) + 1;
    if (keptNbView%2!=0)
    {
      LogCerr("***** oSinoModeCreation::ProcessListMode() -> Error in the implementation of the calculation of the number of views kept for limited angle tomography !" << endl);
      LogCerr("                                              This number should be even, and it is not ! Please revise the code." << endl);
      return 1;
    }
    if (m_Verbose>=2) LogCout ("      Keep " << keptNbView << " covering " << 1.8*f_LimitedAnglePercent << " degrees" << endl);
    // For some scanners, shift the views from half a block, depending on the definition of the first view
    if (mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH) shiftView = mp_Scanner->GetNbTransCrystals()/2;
    else if (mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH2D) shiftView = mp_Scanner->GetNbTransCrystals()/2;
    else if (mp_Scanner->GetScannerModel()==SCANNER_MMR2D) shiftView = mp_Scanner->GetNbTransCrystals()/2;
    else if (mp_Scanner->GetScannerModel()==SCANNER_VISION600) shiftView = mp_Scanner->GetNbTransCrystals()/2;
  }

  // ======================================================================================================================
  // Verbose
  if (m_Verbose>=2) LogCout ("  --> Initialization (writing mode: " << m_Weighting << ") ..." << endl);

  // ======================================================================================================================
  // Check special case
  if (mp_Scanner->GetScannerModel()==SCANNER_HRRT && f_FillEqualLORs)
  {
    LogCerr ("***** oSinoModeCreation::ProcessListMode() -> Cannot fill empty LORs for equal numbering with HRRT scanner !" << endl);
    return 1;
  }

  // ======================================================================================================================
  // Open output listmode files

  // Get path and root names
  string baseName = oOutputManager::GetInstance()->GetRootName();
  string pathName = oOutputManager::GetInstance()->GetPathName();

  // Create file for each thread
  FILE** fout = (FILE**)malloc(f_Threads*sizeof(FILE*));
  for (int th=0; th<f_Threads; th++)
  {
    char tmp_th[100]; sprintf(tmp_th,"%d",th);
    string file_name = pathName + baseName + "/" + baseName + "_th" + ((string)tmp_th);
    fout[th] = fopen(file_name.c_str(),"wb");
    if (fout[th]==NULL)
    {
      LogCerr ("***** oSinoModeCreation::ProcessListMode() -> Failed to create output file '" << file_name << "' !" << endl);
      return 1;
    }
  }

  // ======================================================================================================================
  // Get some sinogram infos	
  int nb_axial_crystals = mp_Scanner->GetNbAxialCrystals();
  int nb_total_trans_crystals = mp_Scanner->GetNbTotalTransCrystals();
  int nb_trans_crystals_in_head = mp_Scanner->GetNbTransCrystalsInHead();
  int scanner_model = mp_Scanner->GetScannerModel();

  // Get the mashing power
  int mash = 1;
  if (scanner_model==SCANNER_BIOGRAPH) mash = mp_BiographTable->GetMash();
  else if (scanner_model==SCANNER_CASTOR) mash = mp_CastorTable->GetMash();
  else if (scanner_model==SCANNER_BIOGRAPH2D) mash = mp_BiographTable2D->GetMash();
  else if (scanner_model==SCANNER_MMR2D) mash = mp_MMRTable2D->GetMash();
  else if (scanner_model==SCANNER_VISION600) mash = mp_Vision600Table->GetMash();
  else if (scanner_model==SCANNER_HRPLUS) mash = mp_HRplusTable->GetMash();
  else if (scanner_model==SCANNER_HRRT) mash = mp_HrrtTable->GetMash();
  else if (scanner_model==SCANNER_INVEON) mash = mp_InveonTable->GetMash();
  else if (scanner_model==SCANNER_FOCUS) mash = mp_FocusTable->GetMash();

  // ======================================================================================================================
  // Do a first pass in the sinogram to remove all neglected bins
  int rest_bins = 0;
  for (int s=0; s<m_ReadNbSino; s++)
//    for (int v=0; v<m_ReadNbView; v++)
//    for (int v=0; v<keptNbView; v++)
    for (int view=m_ReadNbView+shiftView-keptNbView/2; view<m_ReadNbView+shiftView+keptNbView/2; view++)
    {
      int v = view%m_ReadNbView;
      for (int e=0; e<m_ReadNbElem; e++)
      {
        if (m_Weighting!=ALL_WEIGHTING)
        {
          if (m_TrueBool && mp_TrueSino[s][v][e]<=0) mp_TrueSino[s][v][e] = SINO_BIN_TAG;
          else if (m_PromptBool && mp_PromptSino[s][v][e]<=0) mp_PromptSino[s][v][e] = SINO_BIN_TAG;
        }
        else rest_bins++;
      }
    }

  // ======================================================================================================================
  // Build a vector containing all valid sino bin indices
  int* indices = (int*)malloc(rest_bins*sizeof(int));
//  int the_sino_index = 0;
  int the_index = 0;
  int* indices_to_be_placed = (int*)malloc(rest_bins*sizeof(int));
  for (int s=0; s<m_ReadNbSino; s++)
//    for (int v=0; v<m_ReadNbView; v++)
//    for (int v=0; v<keptNbView; v++)
    for (int view=m_ReadNbView+shiftView-keptNbView/2; view<m_ReadNbView+shiftView+keptNbView/2; view++)
    {
      int v = view%m_ReadNbView;
      for (int e=0; e<m_ReadNbElem; e++)
      {
        if ( m_FloatBool ||
             (m_TrueBool && mp_TrueSino[s][v][e]!=SINO_BIN_TAG) ||
             (m_PromptBool && mp_PromptSino[s][v][e]!=SINO_BIN_TAG) )
        {
          int the_sino_index = s*m_ReadNbView*m_ReadNbElem + v*m_ReadNbElem + e;
          indices_to_be_placed[the_index] = the_sino_index;
          the_index++;
        }
//        the_sino_index++;
      }
    }

  // ======================================================================================================================
  // Shuffle the vector (to avoid pattern artifacts...)
  if (m_Verbose>=2) LogCout ("  --> Shuffle sinogram indices ..." << endl);
  int rest_bins_to_be_placed = rest_bins-1;
  for (int i=0; i<rest_bins-1; i++)
  {
    if (i%1000000==0 && m_Verbose>=2) cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
                                           << "      " << ((double)i)*100./((double)rest_bins) << " %                     " << flush;
    int shoot = rand()%rest_bins_to_be_placed;
    indices[i] = indices_to_be_placed[shoot];
    indices_to_be_placed[shoot] = indices_to_be_placed[rest_bins_to_be_placed];
    rest_bins_to_be_placed--;
  }
  indices[rest_bins-1] = indices_to_be_placed[0];
  if (m_Verbose>=2) cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
                         << "      100 %                     " << endl << flush;
  free(indices_to_be_placed);

  // ======================================================================================================================
  // Loop to create the whole list-mode

  // Verbose
  if (m_Verbose>=2) LogCout ("  --> Process sinograms ..." << endl);

  #ifdef OMP_MODE
  // Threads with OpenMP
  omp_set_num_threads(f_Threads);
  #endif

  // For eventual failure
  bool* failed = (bool*)malloc(f_Threads*sizeof(bool));
  for (int th=0; th<f_Threads; th++) failed[th] = false;

  // Crystal indices buffers for mult-thread
  unsigned int** crystalID1 = (unsigned int**)malloc(f_Threads*sizeof(unsigned int*));
  unsigned int** crystalID2 = (unsigned int**)malloc(f_Threads*sizeof(unsigned int*));
  int nb_crystals_in_buffer = 100;
  for (int th=0; th<f_Threads; th++)
  {
    crystalID1[th] = (unsigned int*)malloc(nb_crystals_in_buffer*sizeof(unsigned int));
    crystalID2[th] = (unsigned int*)malloc(nb_crystals_in_buffer*sizeof(unsigned int));
  }

  // The main loop
  int* rejected_bins = (int*)calloc(f_Threads,sizeof(int));
  int i;

  #pragma omp parallel for private(i)
  for (i=0; i<rest_bins; i++)
  {
    // Get the thread number
    int th = 0;
    #ifdef OMP_MODE
    th = omp_get_thread_num();
    #endif

    // Verbose
    if (m_Verbose>=2 && th==0 && i%5000==0) cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
                           << "      bins to be treated: " << rest_bins/f_Threads-1-i << "             " << flush;

    // The number of written data to be checked later
    int nb_data_written = 0;

    // -------------------------------------------------------------------------------------------------------------------------
    // 1. Determine the sinogram bin to be store and contributing LORs
    // -------------------------------------------------------------------------------------------------------------------------

    // Compute sinogram bin
    int sino_bin = indices[i] / (m_ReadSinoSize);
    int view_bin = indices[i] % (m_ReadSinoSize) / m_ReadNbElem;
    int elem_bin = indices[i] % m_ReadNbElem;

    // Get norm factor and continue if null or negative
    float norm_factor  = ComputeNormalizationFactor(elem_bin, view_bin, sino_bin);
    if (norm_factor<=0.)
    {
      rejected_bins[th]++;
      continue;
    }

    // Get the number of contributing LORs to this sinogram bin
    int nb_lors = 0;
    if (scanner_model==SCANNER_BIOGRAPH)
    {
      nb_lors = mp_BiographTable->GetLORsFromSinoBin( sino_bin, view_bin, elem_bin, crystalID1[th], crystalID2[th] );
    }
    else if (scanner_model==SCANNER_CASTOR)
    {
      nb_lors = mp_CastorTable->GetLORsFromSinoBin( sino_bin, view_bin, elem_bin, crystalID1[th], crystalID2[th] );
    }
    else if (scanner_model==SCANNER_BIOGRAPH2D)
    {
      nb_lors = mp_BiographTable2D->GetLORsFromSinoBin( sino_bin, view_bin, elem_bin, crystalID1[th], crystalID2[th] );
    }
    else if (scanner_model==SCANNER_MMR2D)
    {
      nb_lors = mp_MMRTable2D->GetLORsFromSinoBin( sino_bin, view_bin, elem_bin, crystalID1[th], crystalID2[th] );
    }
    else if (scanner_model==SCANNER_VISION600)
    {

      
      nb_lors = mp_Vision600Table->GetLORsFromSinoBin( sino_bin, view_bin, elem_bin, crystalID1[th], crystalID2[th] );


    }
    else if (scanner_model==SCANNER_HRPLUS)
    {
      nb_lors = mp_HRplusTable->GetLORsFromSinoBin( sino_bin, view_bin, elem_bin, crystalID1[th], crystalID2[th] );
    }
    else if (scanner_model==SCANNER_HRRT)
    {
      nb_lors = mp_HrrtTable->GetLORsFromSinoBin( sino_bin, view_bin, elem_bin, crystalID1[th], crystalID2[th] );
    }
    else if (scanner_model==SCANNER_INVEON)
    {
      nb_lors = mp_InveonTable->GetLORsFromSinoBin( sino_bin, view_bin, elem_bin, crystalID1[th], crystalID2[th] );
    }
    else if (scanner_model==SCANNER_FOCUS)
    {
      nb_lors = mp_FocusTable->GetLORsFromSinoBin( sino_bin, view_bin, elem_bin, crystalID1[th], crystalID2[th] );
    }
    else
    {
      LogCerr ("***** oSinoModeCreation::ProcessListMode()  -> Unknown scanner model !" << endl);
      failed[th] = true;
    }

    // Eventually rejected bin
    if (nb_lors==0)
    {
      rejected_bins[th]++;
      continue;
    }

    // -------------------------------------------------------------------------------------------------------------------------
    // 2. Determine and write the counts and time
    // -------------------------------------------------------------------------------------------------------------------------

    // Number of counts
    short int nb_counts = 0;
    if (m_TrueBool) nb_counts = mp_TrueSino[sino_bin][view_bin][elem_bin];
    else if (m_PromptBool) nb_counts = mp_PromptSino[sino_bin][view_bin][elem_bin];
    else if (m_FloatBool) nb_counts = 1; // By default, even if it will be ignored

    // Write the time (same for all event, except for the first one)
    unsigned int the_time = 0;
    if (i==0) the_time = m_RelativeStartTime;
    else the_time = m_RelativeStopTime-0.1; // We subtract 0.1 millisecond from the stop time

    // -------------------------------------------------------------------------------------------------------------------------
    // 3. Determine the scatter and random rates, and the normalization factor
    // -------------------------------------------------------------------------------------------------------------------------

    // Scatters, randoms and normalization
    float scatter_rate = ComputeScatterRate(elem_bin, view_bin, sino_bin);
    float random_rate  = ComputeRandomRate(elem_bin, view_bin, sino_bin);

    // Scatter rate have to be un-normalized
    if (norm_factor>0.) scatter_rate /= norm_factor;
    else scatter_rate = 0.;

    // Special case if we are in float mode (we put the count value in the scatter rate variable)
    if (m_FloatBool) scatter_rate = mp_FloatSino[sino_bin][view_bin][elem_bin];

    // -------------------------------------------------------------------------------------------------------------------------
    // 4. Write and loop on the number of LORs contributing to the associated bin
    // -------------------------------------------------------------------------------------------------------------------------

    if (!f_Castor)
    {
      // Write as the prompt/delay flag, the number of counts
      nb_data_written += fwrite(&nb_counts,sizeof(short int),1,fout[th]);
      nb_data_written += fwrite(&the_time,sizeof(unsigned int),1,fout[th]);
      // Write rates
      nb_data_written += fwrite(&scatter_rate,sizeof(float),1,fout[th]);
      nb_data_written += fwrite(&random_rate,sizeof(float),1,fout[th]);
      nb_data_written += fwrite(&norm_factor,sizeof(float),1,fout[th]);
      // Write this number into the output list-mode
      nb_data_written += fwrite(&nb_lors,sizeof(unsigned short int),1,fout[th]);
    }
    else
    {
      // CASTOR: write time in uint32
      nb_data_written += fwrite(&the_time,sizeof(unsigned int),1,fout[th]);
      // CASTOR: write ACF in float32
      if (m_AttnCorr)
      {
        float acf = 0.;
        if (m_AttenuationMode==ATTENUATION_FROM_SINO) acf = mp_AttnSino[sino_bin][view_bin][elem_bin];
        else
        {
/* OLD VERSION WHERE EACH CRYSTAL PAIR IS PROJECTED
          for (unsigned short int l=0; l<nb_lors; l++)
            acf += ComputeACF(crystalID1[th][l], crystalID2[th][l]);
          acf /= ((float)nb_lors);
*/
          double x1 = 0.; double y1 = 0.; double z1 = 0.;
          double x2 = 0.; double y2 = 0.; double z2 = 0.;
          for (unsigned short int l=0; l<nb_lors; l++)
          {
            x1 += (mp_Scanner->GetCornerX1(crystalID1[th][l])+mp_Scanner->GetCornerX2(crystalID1[th][l]))/2.;
            y1 += (mp_Scanner->GetCornerY1(crystalID1[th][l])+mp_Scanner->GetCornerY2(crystalID1[th][l]))/2.;
            z1 += (mp_Scanner->GetCornerZ1(crystalID1[th][l])+mp_Scanner->GetCornerZ2(crystalID1[th][l]))/2.;
            x2 += (mp_Scanner->GetCornerX1(crystalID2[th][l])+mp_Scanner->GetCornerX2(crystalID2[th][l]))/2.;
            y2 += (mp_Scanner->GetCornerY1(crystalID2[th][l])+mp_Scanner->GetCornerY2(crystalID2[th][l]))/2.;
            z2 += (mp_Scanner->GetCornerZ1(crystalID2[th][l])+mp_Scanner->GetCornerZ2(crystalID2[th][l]))/2.;
          }
          x1 /= ((double)nb_lors);
          y1 /= ((double)nb_lors);
          z1 /= ((double)nb_lors);
          x2 /= ((double)nb_lors);
          y2 /= ((double)nb_lors);
          z2 /= ((double)nb_lors);
          // Project the ACF
          acf = SiddonAttenuationProjection(x1,y1,z1,x2,y2,z2);
        }
        nb_data_written += fwrite(&acf,sizeof(float),1,fout[th]);
      }
      else nb_data_written++;
      // CASTOR: write random rate in float32
      if (m_RandCorr) nb_data_written += fwrite(&random_rate,sizeof(float),1,fout[th]);
      else nb_data_written++;
      // CASTOR: write norm factor in float32
      if (m_NormCorr) nb_data_written += fwrite(&norm_factor,sizeof(float),1,fout[th]);
      else nb_data_written++;
      // CASTOR: write prompts in float 32
      float prompt = ((float)nb_counts);
      nb_data_written += fwrite(&prompt,sizeof(float),1,fout[th]);
      // CASTOR: write scatter rate in float 32
      if (m_ScatCorr) nb_data_written += fwrite(&scatter_rate,sizeof(float),1,fout[th]);
      else nb_data_written++;
      // CASTOR: write nb_lors as uint16, only if more than 1
      if ( (m_ReadSpan/2+1)*mash > 1 ) nb_data_written += fwrite(&nb_lors,sizeof(unsigned short int),1,fout[th]);
    }

    // Launch the loop on this number
    for (unsigned short int l=0; l<nb_lors; l++)
    {
      // -------------------------------------------------------------------------------------------------------------------------
      // Project the attenuation correction factor using SIDDON
      // -------------------------------------------------------------------------------------------------------------------------

      if (!f_Castor)
      {
        // Compute ACF
        float acf = 1.;
        // For inveon or focus, the sinogram can directly be provided
        if (m_AttenuationMode==ATTENUATION_FROM_SINO) acf = mp_AttnSino[sino_bin][view_bin][elem_bin];
        else acf = ComputeACF(crystalID1[th][l], crystalID2[th][l]);
        // Check
        if (acf<1.)
        {
          LogCerr ("***** oSinoModeCreation::ProcessListMode() -> Have computed an invalid ACF (" << acf << ") !" << endl);
         failed[th] = true;
        }
        // Write crystal IDs
        nb_data_written += fwrite(&crystalID1[th][l],sizeof(unsigned int),1,fout[th]);
        nb_data_written += fwrite(&crystalID2[th][l],sizeof(unsigned int),1,fout[th]);
        // Write ACF
        nb_data_written += fwrite(&acf,sizeof(float),1,fout[th]);
      }
      else
      {
        // CASTOR: write crystal IDs as uint32
        nb_data_written += fwrite(&crystalID1[th][l],sizeof(unsigned int),1,fout[th]);
        nb_data_written += fwrite(&crystalID2[th][l],sizeof(unsigned int),1,fout[th]);
      }
    }

    // We add a second loop here to fill the number of LORs included in the event so that there are
    // always the same number of LORs stored in the file (they are filled with -1 values)
    if (f_FillEqualLORs)
    {
      unsigned short int nb_rest_lors = (m_ReadSpan/2+1)*mash - nb_lors;
      for (unsigned short int r=0; r<nb_rest_lors; r++)
      {
        // We write -1 as crystal indices
        int crystal_factice_index = 0;
        // We write -1 as the ACF
        float factice_acf = -1.;
        if (!f_Castor)
        {
          nb_data_written += fwrite(&crystal_factice_index,sizeof(unsigned int),1,fout[th]);
          nb_data_written += fwrite(&crystal_factice_index,sizeof(unsigned int),1,fout[th]);
          nb_data_written += fwrite(&factice_acf,sizeof(float),1,fout[th]);
        }
        else
        {
          // CASTOR: write garbage ids
          nb_data_written += fwrite(&crystal_factice_index,sizeof(unsigned int),1,fout[th]);
          nb_data_written += fwrite(&crystal_factice_index,sizeof(unsigned int),1,fout[th]);
        }
      }
      nb_lors += nb_rest_lors;
    }

    // Check if all data were written
    int nb_data_to_be_written = 0;
    if (!f_Castor) nb_data_to_be_written = 6 + nb_lors*3;
    else {nb_data_to_be_written = nb_lors*2 + 6; if (nb_lors>1) nb_data_to_be_written++;}
    if (nb_data_written!=nb_data_to_be_written)
    {
      LogCerr ("***** oSinoModeCreation::ProcessListMode() -> Failed to write all data (" << nb_data_to_be_written << " in output esteban file (" << nb_data_written << " written) for sinogram bin [" << sino_bin << ";" << view_bin << ";" << elem_bin << "] !" << endl);
      failed[th] = true;
    }
  }

  // Verbose
  if (m_Verbose>=2) cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
                         << "      bins to be treated: 0             " << endl;

  // Close file
  for (int th=0; th<f_Threads; th++) fclose(fout[th]);

  // Check failure
  for (int th=0; th<f_Threads; th++)
  {
    if (failed[th])
    {
      LogCerr ("***** One or multiple threads have failed in writing list-mode file !" << endl);
      return 1;
    }
  }

  // Correct the real number of bins
  for (int th=0; th<f_Threads; th++) rest_bins -= rejected_bins[th];

  // ======================================================================================================================
  // Merge files created by each thread
  if (m_Verbose>=2) LogCout ("  --> Merge files ..." << endl);

  // Really merge if more than 1 thread
  if (f_Threads>1)
  {
    // Merge
    string instruction = "cat ";
    for (int th=0; th<f_Threads; th++)
    {
      char tmp_th[100]; sprintf(tmp_th,"%d",th);
      string file_name = pathName + baseName + "/" + baseName + "_th" + ((string)tmp_th);
      instruction += file_name + " ";
    }
    if (!f_Castor) instruction += "> " + pathName + baseName + "/" + baseName + ".elm";
    else instruction += "> " + pathName + baseName + "/" + baseName + ".cdf";
    int err = system(instruction.c_str());
    // Remove
    for (int th=0; th<f_Threads; th++)
    {
      char tmp_th[100]; sprintf(tmp_th,"%d",th);
      string file_name = pathName + baseName + "/" + baseName + "_th" + ((string)tmp_th);
      instruction = "rm -f " + file_name;
      int err = system(instruction.c_str());
    }
  }
  // Otherwise, just rename the file
  else
  {
    string instruction = "";
    if (!f_Castor) instruction = "mv " + pathName + baseName + "/" + baseName + "_th0 " + pathName + baseName + "/" + baseName + ".elm";
    else instruction = "mv " + pathName + baseName + "/" + baseName + "_th0 " + pathName + baseName + "/" + baseName + ".cdf";
    int err = system(instruction.c_str());
  }

  // ======================================================================================================================
  // Create the header file
  string header_file = "";
  if (!f_Castor) header_file = pathName + baseName + "/" + baseName + ".elm.hdr";
  else header_file = pathName + baseName + "/" + baseName + ".cdh";
  if (m_Verbose>=2) LogCout ("  --> Write header file '" << header_file << "' ..." << endl);
  // Open the file
  ofstream fhdr(header_file.c_str());
  if (!fhdr)
  {
    LogCerr ("***** oSinoModeCreation::ProcessListMode() -> Failed to create output header file '" << header_file << "' !" << endl);
    return 1;
  }
  if (!f_Castor)
  {
    // Write some infos in the header file
    fhdr << "List-mode := " << baseName << ".elm" << endl;
    fhdr << "Frame := " << m_Frame << endl;
    fhdr << "Nb events := " << rest_bins << endl;
    fhdr << "Start time (sec) := " << m_StartTime << endl;
    fhdr << "Stop time (sec) := " << m_StartTime+m_Duration << endl;
    if (mp_Scanner->GetScannerModel()==SCANNER_HRPLUS) fhdr << "PET system := HR+" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_INVEON) fhdr << "PET system := Inveon" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_FOCUS) fhdr << "PET system := Focus" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_HRRT) fhdr << "PET system := HRRT" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH) fhdr << "PET system := Biograph" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_CASTOR) fhdr << "PET system := Castor" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH2D) fhdr << "PET system := Biograph2D" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_MMR2D) fhdr << "PET system := MMR2D" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_VISION600) fhdr << "PET system := Vision600" << endl;
    fhdr << "Input := sinogram" << endl;
    fhdr << "Span := " << m_ReadSpan << endl;
    fhdr << "MaxRingDiff := " << m_ReadMaxRingDiff << endl;
    fhdr << "Mashing := " << mash << endl;
    fhdr << "ECF := " << m_ECF * ((double)m_NbReplicates) << endl;
    fhdr << "Isotope := " << m_Isotope << endl;
    fhdr << "Random generator seed := " << m_Seed << endl;
    if (f_FillEqualLORs) fhdr << "Equal number of LORs := YES" << endl;
    else fhdr << "Equal number of LORs := NO" << endl;
    if (m_FloatBool) fhdr << "Float precorrected data := YES" << endl;
    else fhdr << "Float precorrected data := NO" << endl;
  }
  else
  {
    // Write castor header file
    if (scanner_model==SCANNER_BIOGRAPH || scanner_model==SCANNER_BIOGRAPH2D || scanner_model==SCANNER_MMR2D  || scanner_model==SCANNER_VISION600  ||scanner_model==SCANNER_HRPLUS || scanner_model==SCANNER_CASTOR)
    {
      fhdr << "Data filename: " << baseName << ".cdf" << endl;
      fhdr << "Number of events: " << rest_bins << endl;
      fhdr << "Data mode: histogram" << endl;
      fhdr << "Data type: pet" << endl;
      fhdr << "Start time (s): " << m_StartTime << endl;
      fhdr << "Duration (s): " << m_Duration << endl;
      if (scanner_model==SCANNER_BIOGRAPH) fhdr << "Scanner name: PET_SIEMENS_BIOGRAPH6_TRUEPOINT_TRUEV" << endl;
      else if (scanner_model==SCANNER_CASTOR) fhdr << "Scanner name: PET_CASTOR_BENCHMARK" << endl;
      else if (scanner_model==SCANNER_BIOGRAPH2D) fhdr << "Scanner name: PET_SIEMENS_BIOGRAPH6_TRUEPOINT_2D" << endl;
      else if (scanner_model==SCANNER_MMR2D) fhdr << "Scanner name: PET_SIEMENS_MMR_2D" << endl;
      else if (scanner_model==SCANNER_VISION600) fhdr << "Scanner name: PET_SIEMENS_VISION600" << endl;
      else if (scanner_model==SCANNER_HRPLUS) fhdr << "Scanner name: PET_SIEMENS_HRPLUS" << endl;
      fhdr << "Maximum number of lines per event: " << (m_ReadSpan/2+1)*mash << endl;
      fhdr << "Axial compression: " << m_ReadSpan << endl;
      fhdr << "Azymutal compression: " << mash << endl;
      fhdr << "Max ring diff: " << m_ReadMaxRingDiff << endl;
      fhdr << "Calibration factor: " << m_ECF << endl;
      fhdr << "Isotope: " << m_Isotope << endl;
      if (m_AttnCorr) fhdr << "Attenuation correction flag: 1" << endl;
      else fhdr << "Attenuation correction flag: 0" << endl;
      if (m_NormCorr) fhdr << "Normalization correction flag: 1" << endl;
      else fhdr << "Normalization correction flag: 0" << endl;
      if (m_ScatCorr) fhdr << "Scatter correction flag: 1" << endl;
      else fhdr << "Scatter correction flag: 0" << endl;
      if (m_RandCorr) fhdr << "Random correction flag: 1" << endl;
      else fhdr << "Random correction flag: 0" << endl;
    }
    else
    {
      LogCerr("!!!!! Unsupported scanner for CASTOR... sorry." << endl);
    }
  }
  // Close the file
  fhdr.close();

  // Ending
  return 0;
}

