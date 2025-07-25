#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <ctime>
#include <fstream>
#ifdef OMP_MODE
#include <omp.h>
#endif
#include "oSimulator.hh"
#include "oTableHRplus.hh"
#include "oTableBiograph.hh"
#include "oTableHRRT.hh"
#include "oTableInveon.hh"
#include "oTableMMR2D.hh"
#include "oTableVision600.hh"
#include "oMiscellaneous.hh"
#include "oOutputManager.hh"
using namespace std;

// ==========================================================================================================================================
// Constructor
// ==========================================================================================================================================
oSimulator::oSimulator(int f_NbThreads, int f_Seed, int f_Verbose)
{
  // Boolean
  m_HaveInitScannerStuff = false;
  m_HaveInitInputImage = false;
  m_HaveInitPSF = false;

  // Default
  mp_Scanner = NULL;
  mp_TableBiograph = NULL;
  mp_TableBiograph2D = NULL;
  mp_TableMMR2D = NULL;
  mp_TableVision600 = NULL;
  mp_TableHRRT = NULL;
  mp_TableHRplus = NULL;
  mp_TableInveon = NULL;
  mp_InputImage = NULL;
  mp_MuMapImage = NULL;
  mp_ProjectImage = NULL;
  mp_MuMapImage = NULL;
  mp_SinoForw = NULL;
  mp_SinoScat = NULL;
  mp_SinoRand = NULL;
  mp_SinoNorm = NULL;
  mp_SinoAttn = NULL;
  m_RandCorr = false;
  m_ScatCorr = false;
  m_AttnCorr = false;
  m_FloatBool = false;
  m_NbSino = -1;
  m_NbView = -1;
  m_NbElem = -1;
  m_Mash = -1;
  m_Span = -1;
  m_MaxRingDiff = -1;

  // Affectations
  m_NbThreads = f_NbThreads;
  m_Verbose = f_Verbose;

  // Random generator seed
  m_Seed = Misc_InitRandomGenerator(f_Seed);

  // Checks
  if (m_NbThreads<1)
  {
    LogCerr ("***** oSimulator::Constructor() -> Cannot set less than 1 thread !" << endl);
    exit(1);
  }
  if (m_Verbose<0) m_Verbose = 0;
}
// ==========================================================================================================================================
// Destructor
// ==========================================================================================================================================
oSimulator::~oSimulator()
{

}
// ==========================================================================================================================================
// Function InitScannerStuff
//   This function initialize the scanner, the crystal map and other scanner related stuff.
// ==========================================================================================================================================
int oSimulator::InitScannerStuff(const string& f_ScannerName, const string& f_FileCrystalMap, int f_Mash, int f_Span, int f_MaxRingDiff)
{
  // Verbose
  if (m_Verbose>=1) LogCout ("oSimulator::InitScannerStuff() -> Initialize all scanner related stuff" << endl);

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // For correct initialization
  m_HaveInitScannerStuff = true;

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Affectations
  m_Mash = f_Mash;
  m_Span = f_Span;
  m_MaxRingDiff = f_MaxRingDiff;

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Scanner itself
  mp_Scanner = new oScanner( f_ScannerName, m_Verbose );

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Crystal map
  if (mp_Scanner->ReadCrystalMap(f_FileCrystalMap))
  {
    LogCerr ("***** oSimulator::InitScannerStuff() -> An error occured while reading the crystal map !" << endl);
    return 1;
  }

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Scanner table
  if (mp_Scanner->GetScannerModel()==SCANNER_HRPLUS)
  {
    // Check the mashing consistency
    if (mp_Scanner->GetNbView()%f_Mash!=0)
    {
      LogCerr ("***** oSimulator::InitScannerStuff() -> Mashing level (" << f_Mash << ") inconsistent with default number of views (" << mp_Scanner->GetNbView() << ") !" << endl);
      return 1;
    }
    // Build the tables
    mp_TableHRplus = new oTableHRplus(mp_Scanner->GetNbElem(), mp_Scanner->GetNbView()/f_Mash, -1, f_MaxRingDiff, f_Span, mp_Scanner, m_Verbose);
    if (mp_TableHRplus->ComputeSpanTables())
    {
      LogCerr ("***** oSimulator::InitScannerStuff() -> Problem while computing the span tables !" << endl);
      return 1;
    }
    if (mp_TableHRplus->ComputeRingTables())
    {
      LogCerr ("***** oSimulator::InitScannerStuff() -> Problem while computing the ring tables !" << endl);
      return 1;
    }
    // Get mashing power and sinogram dimensions
    m_Mash   = mp_TableHRplus->GetMash();
    m_NbSino = mp_TableHRplus->GetNbSino();
    m_NbView = mp_TableHRplus->GetNbView();
    m_NbElem = mp_TableHRplus->GetNbElem();
  }
  else if (mp_Scanner->GetScannerModel()==SCANNER_INVEON)
  {
    // Check the mashing consistency
    if (mp_Scanner->GetNbView()%f_Mash!=0)
    {
      LogCerr ("***** oSimulator::InitScannerStuff() -> Mashing level (" << f_Mash << ") inconsistent with default number of views (" << mp_Scanner->GetNbView() << ") !" << endl);
      return 1;
    }
    // Build the tables
    mp_TableInveon = new oTableInveon(mp_Scanner->GetNbElem(), mp_Scanner->GetNbView()/f_Mash, -1, f_MaxRingDiff, f_Span, mp_Scanner, m_Verbose);
    if (mp_TableInveon->ComputeSpanTables())
    {
      LogCerr ("***** oSimulator::InitScannerStuff() -> Problem while computing the span tables !" << endl);
      return 1;
    }
    if (mp_TableInveon->ComputeRingTables())
    {
      LogCerr ("***** oSimulator::InitScannerStuff() -> Problem while computing the ring tables !" << endl);
      return 1;
    }
    // Get mashing power and sinogram dimensions
    m_Mash   = mp_TableInveon->GetMash();
    m_NbSino = mp_TableInveon->GetNbSino();
    m_NbView = mp_TableInveon->GetNbView();
    m_NbElem = mp_TableInveon->GetNbElem();
  }
  else if (mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH)
  {
    // Build the tables first
    mp_TableBiograph = new oTableBiograph(mp_Scanner, m_Verbose);
    // Read span table
    if (mp_TableBiograph->ReadSpanTable())
    {
      LogCerr ("***** oSimulator::InitScannerStuff() -> Problem while reading span table !" << endl);
      return 1;
    }
    // Compute ring tables
    if (mp_TableBiograph->ComputeRingTables())
    {
      LogCerr ("***** oSimulator::InitScannerStuff() -> Problem while computing ring tables !" << endl);
      return 1;
    }
    // Force given mash, span and maxringdiff values
    if (m_Span!=mp_TableBiograph->GetSpan())
    {
      LogCerr ("!!!!! oSimulator::InitScannerStuff() -> Force span to the system's one (" << mp_TableBiograph->GetSpan() << ") !" << endl);
      m_Span = mp_TableBiograph->GetSpan();
    }
    if (m_MaxRingDiff!=mp_TableBiograph->GetMaxRingDiff())
    {
      LogCerr ("!!!!! oSimulator::InitScannerStuff() -> Force maximum ring difference to the system's one (" << mp_TableBiograph->GetMaxRingDiff() << ") !" << endl);
      m_MaxRingDiff = mp_TableBiograph->GetMaxRingDiff();
    }
    if (m_Mash!=mp_TableBiograph->GetMash())
    {
      LogCerr ("!!!!! oSimulator::InitScannerStuff() -> Force mashing to the system's one (" << mp_TableBiograph->GetMash() << ") !" << endl);
      m_Mash = mp_TableBiograph->GetMash();
    }
    // Get sinogram dimensions
    m_NbSino = mp_TableBiograph->GetNbSino();
    m_NbView = mp_TableBiograph->GetNbView();
    m_NbElem = mp_TableBiograph->GetNbElem();
  }
  else if (mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH2D)
  {
    // Build the tables first
    mp_TableBiograph2D = new oTableBiograph2D(mp_Scanner, m_Verbose);
    // Compute ring tables
    if (mp_TableBiograph2D->ComputeRingTables())
    {
      LogCerr ("***** oSimulator::InitScannerStuff() -> Problem while computing ring tables !" << endl);
      return 1;
    }
    // Force given mash, span and maxringdiff values
    if (m_Span!=mp_TableBiograph2D->GetSpan())
    {
      LogCerr ("!!!!! oSimulator::InitScannerStuff() -> Force span to the system's one (" << mp_TableBiograph2D->GetSpan() << ") !" << endl);
      m_Span = mp_TableBiograph2D->GetSpan();
    }
    if (m_MaxRingDiff!=mp_TableBiograph2D->GetMaxRingDiff())
    {
      LogCerr ("!!!!! oSimulator::InitScannerStuff() -> Force maximum ring difference to the system's one (" << mp_TableBiograph2D->GetMaxRingDiff() << ") !" << endl);
      m_MaxRingDiff = mp_TableBiograph2D->GetMaxRingDiff();
    }
    if (m_Mash!=mp_TableBiograph2D->GetMash())
    {
      LogCerr ("!!!!! oSimulator::InitScannerStuff() -> Force mashing to the system's one (" << mp_TableBiograph2D->GetMash() << ") !" << endl);
      m_Mash = mp_TableBiograph2D->GetMash();
    }
    // Get sinogram dimensions
    m_NbSino = mp_TableBiograph2D->GetNbSino();
    m_NbView = mp_TableBiograph2D->GetNbView();
    m_NbElem = mp_TableBiograph2D->GetNbElem();
  }
  else if (mp_Scanner->GetScannerModel()==SCANNER_MMR2D)
  {
    // Build the tables first
    mp_TableMMR2D = new oTableMMR2D(mp_Scanner, m_Verbose);
    // Compute ring tables
    if (mp_TableMMR2D->ComputeRingTables())
    {
      LogCerr ("***** oSimulator::InitScannerStuff() -> Problem while computing ring tables !" << endl);
      return 1;
    }
    // Force given mash, span and maxringdiff values
    if (m_Span!=mp_TableMMR2D->GetSpan())
    {
      LogCerr ("!!!!! oSimulator::InitScannerStuff() -> Force span to the system's one (" << mp_TableMMR2D->GetSpan() << ") !" << endl);
      m_Span = mp_TableMMR2D->GetSpan();
    }
    m_MaxRingDiff = 0;
    if (m_Mash!=mp_TableMMR2D->GetMash())
    {
      LogCerr ("!!!!! oSimulator::InitScannerStuff() -> Force mashing to the system's one (" << mp_TableMMR2D->GetMash() << ") !" << endl);
      m_Mash = mp_TableMMR2D->GetMash();
    }
    // Get sinogram dimensions
    m_NbSino = mp_TableMMR2D->GetNbSino();
    m_NbView = mp_TableMMR2D->GetNbView();
    m_NbElem = mp_TableMMR2D->GetNbElem();
  }
  else if (mp_Scanner->GetScannerModel()==SCANNER_VISION600)
  {
    // Build the tables first
    mp_TableVision600 = new oTableVision600(mp_Scanner, m_Verbose);
    // Compute ring tables
    if (mp_TableVision600->ComputeRingTables())
    {
      LogCerr ("***** oSimulator::InitScannerStuff() -> Problem while computing ring tables !" << endl);
      return 1;
    }
    // Force given mash, span and maxringdiff values
    if (m_Span!=mp_TableVision600->GetSpan())
    {
      LogCerr ("!!!!! oSimulator::InitScannerStuff() -> Force span to the system's one (" << mp_TableVision600->GetSpan() << ") !" << endl);
      m_Span = mp_TableVision600->GetSpan();
    }
    m_MaxRingDiff = 0;
    if (m_Mash!=mp_TableVision600->GetMash())
    {
      LogCerr ("!!!!! oSimulator::InitScannerStuff() -> Force mashing to the system's one (" << mp_TableVision600->GetMash() << ") !" << endl);
      m_Mash = mp_TableVision600->GetMash();
    }
    // Get sinogram dimensions
    m_NbSino = mp_TableVision600->GetNbSino();
    m_NbView = mp_TableVision600->GetNbView();
    m_NbElem = mp_TableVision600->GetNbElem();
  }
  else if (mp_Scanner->GetScannerModel()==SCANNER_HRRT)
  {
    LogCerr ("***** oSimulator::InitScannerStuff() -> Not yet implemented for this scanner !" << endl);
    return 1;
  }

  // End
  return 0;
}
// ==========================================================================================================================================
// Function InitInputImages
//   This function initialize and read the input images to be projected (emission and transmission).
// ==========================================================================================================================================
int oSimulator::InitInputImages(const string& f_FileImage, const string& f_FileAttenuation, float f_OffsetX, float f_OffsetY, float f_OffsetZ, int f_Projector)
{
  // Verbose
  if (m_Verbose>=1) LogCout ("oSimulator::InitInputImages() -> Initialize input images" << endl);

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // For correct initialization
  m_HaveInitInputImage = true;

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Emission image
  // ----------------------------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------------------------

  // Affect the offset
  m_OffsetX = f_OffsetX;
  m_OffsetY = f_OffsetY;
  m_OffsetZ = f_OffsetZ;

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // We have to read the interfile header, so we first open it
  ifstream fhead(f_FileImage.c_str());
  if (!fhead)
  {
    LogCerr ("***** oSimulator::InitInputImages() -> Input image header file '" << f_FileImage << "' is missing or corruped !" << endl);
    return 1;
  }

  // Verbose
  if (m_Verbose>=1) LogCout ("  --> Read input image header file '" << f_FileImage << "' ..." << endl);

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Now we have to find the data file name, the dimensions and voxel sizes
  string data_name = "";
  int data_found = 0;
  int data_to_be_found = 7;
  char line[1024];
  fhead.getline(line,1024);

  while (!fhead.eof() && data_found<data_to_be_found)
  {
    size_t found;
    // Transform it to string to benefit from useful c++ functions
    string test = (string)line;
    // Test 1: for the data file name
    found = test.find("name of data file");
    if (found!=string::npos)
    {
      // Get the file name
      found = test.find("=");
      data_name = test.substr(found+1);
      // Remove blancks and cariage return
      while ( (found=data_name.find(" ")) != string::npos) data_name.replace(found,1,"");
      while ( (found=data_name.find("\r")) != string::npos) data_name.replace(found,1,"");
      // Increment the number of data found
      data_found++;
    }
    // Test 2,3,4: for the dimXYZ
    found = test.find("matrix size");
    if (found!=string::npos)
    {
      // Get the number as a string
      found = test.find("=");
      string number = test.substr(found+1);
      // Remove blancks
      while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
      // Test the dimension and get the value
      if (test.find("[1]")!=string::npos)      m_DimX = atoi(number.c_str());
      else if (test.find("[2]")!=string::npos) m_DimY = atoi(number.c_str());
      else if (test.find("[3]")!=string::npos) m_DimZ = atoi(number.c_str());
      // Increment the number of data found
      data_found++;
    }
    // Test 5,6,7: for the voxel size XYZ
    if (test.find("scale factor (mm/pixel)")!=string::npos || test.find("scaling factor (mm/pixel)")!=string::npos)
    {
      // Get the number as a string
      found = test.find("=");
      string number = test.substr(found+1);
      // Remove blancks
      while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
      // Test the dimension and get the value
      if (test.find("[1]")!=string::npos)      m_VoxSizeX = atof(number.c_str());
      else if (test.find("[2]")!=string::npos) m_VoxSizeY = atof(number.c_str());
      else if (test.find("[3]")!=string::npos) m_VoxSizeZ = atof(number.c_str());
      // Increment the number of data found
      data_found++;
    }
    // Read a new line
    fhead.getline(line,1024);
  }

  // Close the header file
  fhead.close();

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Check if all data were found
  if (data_found!=data_to_be_found)
  {
    LogCerr ("***** oSimulator::InitInputImages() -> Failed to get all data in header file '" << f_FileImage << "' !" << endl);
    return 1;
  }

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Allocate image
  m_DimXY  = m_DimX * m_DimY;
  m_DimTot = m_DimXY * m_DimZ;
  mp_InputImage = (float*)malloc(m_DimTot*sizeof(float));

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Calculate FOV size
  m_FOVSizeX = ((float)m_DimX)*m_VoxSizeX;
  m_FOVSizeY = ((float)m_DimY)*m_VoxSizeY;
  m_FOVSizeZ = ((float)m_DimZ)*m_VoxSizeZ;

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Add the relative path to the data file name
  int pos; if ((pos=f_FileImage.find_last_of("/"))!=string::npos) data_name = f_FileImage.substr(0,pos) + "/" + data_name;

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Open data file
  FILE* fdata = fopen(data_name.c_str(),"rb");
  if (fdata==NULL)
  {
    LogCerr ("***** oSimulator::InitInputImages() -> Input data file '" << data_name << "' is missing or corrupted !" << endl);
    return 1;
  }

  // Verbose
  if (m_Verbose>=1)
  {
    LogCout ("  --> Dimensions [" << m_DimX << ";" << m_DimY << ";" << m_DimZ << "]" << endl);
    LogCout ("  --> Voxel size [" << m_VoxSizeX << ";" << m_VoxSizeY << ";" << m_VoxSizeZ << "] mm" << endl);
    LogCout ("  --> FOV size [" << m_FOVSizeX << ";" << m_FOVSizeY << ";" << m_FOVSizeZ << "] mm" << endl);
    LogCout ("  --> Reading input image from file '" << data_name << "' ..." << endl);
  }

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Read input image
  int nb_data_read = fread(&mp_InputImage[0],sizeof(float),m_DimTot,fdata);

  // Close data file
  fclose(fdata);

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Check
  if (nb_data_read!=m_DimTot)
  {
    LogCerr ("***** oSimulator::InitInputImages() -> Failed to read all data (" << m_DimTot << ") in image (" << nb_data_read << " read) !" << endl);
    return 1;
  }

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Apply projector dependent correction factor

  // Here we always use Siddon. As this factor will be applied in the reconstruction, we have to apply it here to be exact
  // This factor does not depend on the voxel size (Siddon gives unit per volume)
/* No need to do that
  if (f_Projector==0)
  {
    PRECISION siddon_projector_factor = 0.986004159;
    for (int v=0; v<m_DimTot; v++) mp_InputImage[v] /= siddon_projector_factor;
  }
  else if (f_Projector==1)
  {
//    PRECISION siddon_projector_factor = 0.986004159;
    PRECISION siddon_projector_factor = 1.;
    for (int v=0; v<m_DimTot; v++) mp_InputImage[v] /= siddon_projector_factor;
  }
*/

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Transmission image
  // ----------------------------------------------------------------------------------------------------------------------------------------
  // ----------------------------------------------------------------------------------------------------------------------------------------

  if (f_FileAttenuation!="")
  {
    // ----------------------------------------------------------------------------------------------------------------------------------------
    // We have to read the interfile header, so we first open it
    ifstream fhead(f_FileAttenuation.c_str());
    if (!fhead)
    {
      LogCerr ("***** oSimulator::InitInputImages() -> Input mumap header file '" << f_FileAttenuation << "' is missing or corruped !" << endl);
      return 1;
    }

    // Verbose
    if (m_Verbose>=1) LogCout ("  --> Read input mumap header file '" << f_FileAttenuation << "' ..." << endl);

    // ----------------------------------------------------------------------------------------------------------------------------------------
    // Now we have to find the data file name, the dimensions and voxel sizes
    string data_name = "";
    int data_found = 0;
    int data_to_be_found = 7;
    char line[1024];
    fhead.getline(line,1024);

    // Parameters to be find and compared to the emission image
    int mu_dimX = -1;
    int mu_dimY = -1;
    int mu_dimZ = -1;
    float mu_voxX = -1.;
    float mu_voxY = -1.;
    float mu_voxZ = -1.;

    while (!fhead.eof() && data_found<data_to_be_found)
    {
      size_t found;
      // Transform it to string to benefit from useful c++ functions
      string test = (string)line;
      // Test 1: for the data file name
      found = test.find("name of data file");
      if (found!=string::npos)
      {
        // Get the file name
        found = test.find("=");
        data_name = test.substr(found+1);
        // Remove blancks and cariage return
        while ( (found=data_name.find(" ")) != string::npos) data_name.replace(found,1,"");
        while ( (found=data_name.find("\r")) != string::npos) data_name.replace(found,1,"");
        // Increment the number of data found
        data_found++;
      }
      // Test 2,3,4: for the dimXYZ
      found = test.find("matrix size");
      if (found!=string::npos)
      {
        // Get the number as a string
        found = test.find("=");
        string number = test.substr(found+1);
        // Remove blancks
        while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
        // Test the dimension and get the value
        if (test.find("[1]")!=string::npos)      mu_dimX = atoi(number.c_str());
        else if (test.find("[2]")!=string::npos) mu_dimY = atoi(number.c_str());
        else if (test.find("[3]")!=string::npos) mu_dimZ = atoi(number.c_str());
        // Increment the number of data found
        data_found++;
      }
      // Test 5,6,7: for the voxel size XYZ
      if (test.find("scale factor (mm/pixel)")!=string::npos || test.find("scaling factor (mm/pixel)")!=string::npos)
      {
        // Get the number as a string
        found = test.find("=");
        string number = test.substr(found+1);
        // Remove blancks
        while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
        // Test the dimension and get the value
        if (test.find("[1]")!=string::npos)      mu_voxX = atof(number.c_str());
        else if (test.find("[2]")!=string::npos) mu_voxY = atof(number.c_str());
        else if (test.find("[3]")!=string::npos) mu_voxZ = atof(number.c_str());
        // Increment the number of data found
        data_found++;
      }
      // Read a new line
      fhead.getline(line,1024);
    }

    // Close the header file
    fhead.close();

    // ----------------------------------------------------------------------------------------------------------------------------------------
    // Check if all data were found
    if (data_found!=data_to_be_found)
    {
      LogCerr ("***** oSimulator::InitInputImages() -> Failed to get all data in header file '" << f_FileAttenuation << "' !" << endl);
      return 1;
    }

    // ----------------------------------------------------------------------------------------------------------------------------------------
    // Check all dimensions and allocate image
    if ( mu_dimX!=m_DimX || mu_dimY!=m_DimY || mu_dimZ!=m_DimZ ||
         mu_voxX!=m_VoxSizeX || mu_voxY!=m_VoxSizeY || mu_voxZ!=m_VoxSizeZ )
    {
      LogCerr ("***** oSimulator::InitInputImages() -> Dimensions of the mumap differ from those of the emission image !" << endl);
      return 1;
    }
    mp_MuMapImage = (float*)malloc(m_DimTot*sizeof(float));

    // ----------------------------------------------------------------------------------------------------------------------------------------
    // Add the relative path to the data file name
    int pos; if ((pos=f_FileAttenuation.find_last_of("/"))!=string::npos) data_name = f_FileAttenuation.substr(0,pos) + "/" + data_name;

    // ----------------------------------------------------------------------------------------------------------------------------------------
    // Open data file
    FILE* fdata = fopen(data_name.c_str(),"rb");
    if (fdata==NULL)
    {
      LogCerr ("***** oSimulator::InitInputImages() -> Input data file '" << data_name << "' is missing or corrupted !" << endl);
      return 1;
    }

    // Verbose
    if (m_Verbose>=1) LogCout ("  --> Reading input mumap from file '" << data_name << "' ..." << endl);

    // ----------------------------------------------------------------------------------------------------------------------------------------
    // Read input image
    int nb_data_read = fread(&mp_MuMapImage[0],sizeof(float),m_DimTot,fdata);

    // Close data file
    fclose(fdata);

    // ----------------------------------------------------------------------------------------------------------------------------------------
    // Check
    if (nb_data_read!=m_DimTot)
    {
      LogCerr ("***** oSimulator::InitInputImages() -> Failed to read all data (" << m_DimTot << ") in mumap (" << nb_data_read << " read) !" << endl);
      return 1;
    }

    // Valid
    m_AttnCorr = true;
  }

  // End
  return 0;
}

// ==========================================================================================================================================
// Function InitPSF
//   This function initialize and apply the PSF onto the input image.
// ==========================================================================================================================================
int oSimulator::InitPSF(float f_PsfTransFWHM, float f_PsfAxialFWHM)
{
  // Verbose
  if (m_Verbose>=1) LogCout ("oSimulator::InitPSF() -> Initialize PSF stuff" << endl);

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Check for correct initialization
  if (!m_HaveInitInputImage)
  {
    LogCerr ("***** oSimulator::InitPSF() -> Must call this function after having initialized the input image !" << endl);
    return 1;
  }
  m_HaveInitPSF = true;

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Affectations
  m_PsfTransFWHM = f_PsfTransFWHM;
  m_PsfAxialFWHM = f_PsfAxialFWHM;

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Kernel computation
  if (m_PsfTransFWHM>0. || m_PsfAxialFWHM>0.)
  {
    // Verbose
    if (m_Verbose>=1) LogCout ("  --> Dealing with PSF stuffs ..." << endl);
    // PSF kernel size
    m_PsfKernSizeX = (int) ( m_PsfTransFWHM*CONVOLUTION_SIGMAS_TRUNCATION/(TWO_SQRT_TWO_LN_2*m_VoxSizeX) );
    m_PsfKernSizeY = (int) ( m_PsfTransFWHM*CONVOLUTION_SIGMAS_TRUNCATION/(TWO_SQRT_TWO_LN_2*m_VoxSizeY) );
    m_PsfKernSizeZ = (int) ( m_PsfAxialFWHM*CONVOLUTION_SIGMAS_TRUNCATION/(TWO_SQRT_TWO_LN_2*m_VoxSizeZ) );
    m_PsfKernSizeX = m_PsfKernSizeX*2 + 1;
    m_PsfKernSizeY = m_PsfKernSizeY*2 + 1;
    m_PsfKernSizeZ = m_PsfKernSizeZ*2 + 1;
    // Verbose
    if (m_Verbose >= 1) LogCout ("  --> PSF kernel of [" << m_PsfKernSizeX << ";" << m_PsfKernSizeY << ";" << m_PsfKernSizeZ << "] voxels width for [" << m_PsfTransFWHM << ";" << m_PsfAxialFWHM << "] mm FWHMs (truncation at " << CONVOLUTION_SIGMAS_TRUNCATION << " sigmas) ..." << endl);
    // Allocating kernel
    mp_PsfKernel = (float***)malloc(m_PsfKernSizeX*sizeof(float**));
    for (int p1=0; p1<m_PsfKernSizeX; p1++)
    {
      mp_PsfKernel[p1] = (float**)malloc(m_PsfKernSizeY*sizeof(float*));
      for (int p2=0; p2<m_PsfKernSizeY; p2++) mp_PsfKernel[p1][p2] = (float*)calloc(m_PsfKernSizeZ,sizeof(float));
    }
    // Calculating kernel
    Make3DGaussianKernel();
  }

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Allocating buffer image
  mp_ProjectImage = (float*)malloc(m_DimTot*sizeof(float));

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Either copy the input into the image to be projected
  if (m_PsfTransFWHM<=0. && m_PsfAxialFWHM<=0.)
  {
    for (int v=0; v<m_DimTot; v++) mp_ProjectImage[v] = mp_InputImage[v];
  }
  // Or convolve the input into the one to be projected
  else
  {
    Convolve3D(mp_InputImage,mp_ProjectImage);
  }

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Same for the mumap if any
  if (m_AttnCorr)
  {
    string baseName = oOutputManager::GetInstance()->GetRootName();
    string pathName = oOutputManager::GetInstance()->GetPathName();
    // -----------------------------------------------------------
    // Save mumap before convolving
    // -----------------------------------------------------------
    // Mumap data file name
    string file_data = pathName + baseName + "/" + baseName + "_mumap.i";
    // Open file
    FILE* fmap = fopen(file_data.c_str(),"wb");
    if (fmap==NULL)
    {
      LogCerr("***** oSimulatorLabel::Project() -> Failed to create output file '" << file_data << "' !" << endl);
      return 1;
    }
    // Write file
    int nb_data_written = fwrite(mp_MuMapImage,sizeof(float),m_DimTot,fmap);
    // Close file
    fclose(fmap);
    // Check data
    if (nb_data_written!=m_DimTot)
    {
      LogCerr("***** oSimulatorLabel::Project() -> Failed to write all data in output file '" << file_data << "' !" << endl);
      return 1;
    }
    // Mumap header file name
    string file_head = pathName + baseName + "/" + baseName + "_mumap.i.hdr";
    // Open file
    ofstream fhead(file_head.c_str());
    if (!fhead)
    {
      LogCerr("***** oSimulatorLabel::Project() -> Failed to create output file '" << file_head << "' !" << endl);
      return 1;
    }
    // Write header
    fhead << "!INTERFILE" << endl;
    fhead << "!name of data file := " << baseName << "_mumap.i" << endl;
    if (mp_Scanner->GetScannerModel()==SCANNER_HRPLUS) fhead << "!originating system := HR+" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_HRRT) fhead << "!originating system := HRRT" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH) fhead << "!originating system := Biograph" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH2D) fhead << "!originating system := Biograph2D" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_MMR2D) fhead << "!originating system := MMR2D" << endl;
    fhead << "!PET data type := transmission" << endl;
    fhead << "data format := image" << endl;
    fhead << "number format := float" << endl;
    fhead << "number of bytes per pixel := 4" << endl;
    fhead << "Patient name := " << baseName << endl;
    fhead << "number of dimensions := 3" << endl;
    fhead << "matrix size [1] := " << m_DimX << endl;
    fhead << "matrix size [2] := " << m_DimY << endl;
    fhead << "matrix size [3] := " << m_DimZ << endl;
    fhead << "scaling factor (mm/pixel) [1] := " << m_VoxSizeX << endl;
    fhead << "scaling factor (mm/pixel) [2] := " << m_VoxSizeY << endl;
    fhead << "scaling factor (mm/pixel) [3] := " << m_VoxSizeZ << endl;
    fhead << "FOV size (mm) [1] := " << ((float)m_DimX)*m_VoxSizeX << endl;
    fhead << "FOV size (mm) [2] := " << ((float)m_DimY)*m_VoxSizeY << endl;
    fhead << "FOV size (mm) [3] := " << ((float)m_DimZ)*m_VoxSizeZ << endl;
    fhead << "patient orientation := NOPITCH" << endl;
    fhead << "PSF smoothing := NO" << endl;
    // Close file
    fhead.close();
    // -----------------------------------------------------------
    // Smooth it
    if (m_PsfTransFWHM>0. || m_PsfAxialFWHM>0.)
    {
      // Convolve the mumap into the input image (now empty)
      Convolve3D(mp_MuMapImage,mp_InputImage);
      // Copy back in the mumap matrix
      for (int v=0; v<m_DimTot; v++) mp_MuMapImage[v] = mp_InputImage[v];
      // -----------------------------------------------------------
      // Save mumap after convolving
      // -----------------------------------------------------------
      // Mumap data file name
      string file_data = pathName + baseName + "/" + baseName + "_mumap_psf.i";
      // Open file
      FILE* fmap = fopen(file_data.c_str(),"wb");
      if (fmap==NULL)
      {
        LogCerr("***** oSimulatorLabel::Project() -> Failed to create output file '" << file_data << "' !" << endl);
        return 1;
      }
      // Write file
      int nb_data_written = fwrite(mp_MuMapImage,sizeof(float),m_DimTot,fmap);
      // Close file
      fclose(fmap);
      // Check data
      if (nb_data_written!=m_DimTot)
      {
        LogCerr("***** oSimulatorLabel::Project() -> Failed to write all data in output file '" << file_data << "' !" << endl);
        return 1;
      }
      // Mumap header file name
      string file_head = pathName + baseName + "/" + baseName + "_mumap_psf.i.hdr";
      // Open file
      ofstream fhead(file_head.c_str());
      if (!fhead)
      {
        LogCerr("***** oSimulatorLabel::Project() -> Failed to create output file '" << file_head << "' !" << endl);
        return 1;
      }
      // Write header
      fhead << "!INTERFILE" << endl;
      fhead << "!name of data file := " << baseName << "_mumap_psf.i" << endl;
      if (mp_Scanner->GetScannerModel()==SCANNER_HRPLUS) fhead << "!originating system := HR+" << endl;
      else if (mp_Scanner->GetScannerModel()==SCANNER_HRRT) fhead << "!originating system := HRRT" << endl;
      else if (mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH) fhead << "!originating system := Biograph" << endl;
      else if (mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH2D) fhead << "!originating system := Biograph2D" << endl;
      else if (mp_Scanner->GetScannerModel()==SCANNER_MMR2D) fhead << "!originating system := MMR2D" << endl;
      fhead << "!PET data type := transmission" << endl;
      fhead << "data format := image" << endl;
      fhead << "number format := float" << endl;
      fhead << "number of bytes per pixel := 4" << endl;
      fhead << "Patient name := " << baseName << endl;
      fhead << "number of dimensions := 3" << endl;
      fhead << "matrix size [1] := " << m_DimX << endl;
      fhead << "matrix size [2] := " << m_DimY << endl;
      fhead << "matrix size [3] := " << m_DimZ << endl;
      fhead << "scaling factor (mm/pixel) [1] := " << m_VoxSizeX << endl;
      fhead << "scaling factor (mm/pixel) [2] := " << m_VoxSizeY << endl;
      fhead << "scaling factor (mm/pixel) [3] := " << m_VoxSizeZ << endl;
      fhead << "FOV size (mm) [1] := " << ((float)m_DimX)*m_VoxSizeX << endl;
      fhead << "FOV size (mm) [2] := " << ((float)m_DimY)*m_VoxSizeY << endl;
      fhead << "FOV size (mm) [3] := " << ((float)m_DimZ)*m_VoxSizeZ << endl;
      fhead << "patient orientation := NOPITCH" << endl;
      fhead << "PSF smoothing := YES" << endl;
      // Close file
      fhead.close();
    }
  }

  // Free the input image
  free(mp_InputImage);

  // End
  return 0;
}

// ==========================================================================================================================================
// Function Project
//   This function project the input image into the output sinogram.
// ==========================================================================================================================================
int oSimulator::Project(bool f_SaveImage, int f_Projector)
{
  // Verbose
  if (m_Verbose>=1) LogCout ("oSimulator::Project() -> Project input image into output sinogram" << endl);

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Check for correct initialization
  if (!m_HaveInitScannerStuff)
  {
    LogCerr ("***** oSimulator::Project() -> Must call this function after having initialized the scanner related stuff !" << endl);
    return 1;
  }
  if (!m_HaveInitInputImage)
  {
    LogCerr ("***** oSimulator::Project() -> Must call this function after having initialized the input image !" << endl);
    return 1;
  }
  if (!m_HaveInitPSF)
  {
    LogCerr ("***** oSimulator::Project() -> Must call this function after having initialized the PSF !" << endl);
    return 1;
  }
  m_HaveProject = true;

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Save forward image if asked for
  if (f_SaveImage)
  {
    if (oOutputManager::GetInstance()->SaveOrdinaryImage( mp_ProjectImage, m_DimX, m_DimY, m_DimZ, m_VoxSizeX, m_VoxSizeY, m_VoxSizeZ, mp_Scanner->GetScannerModel(), 1., 0., "INF" ))
    {
      LogCerr ("***** oSimulator::Project() -> An error occured while saving the image to be projected !" << endl);
      return 1;
    }
/*OBS
    if (m_AttnCorr)
    {
      if (oOutputManager::GetInstance()->SaveOrdinaryMuMap( mp_MuMapImage, m_DimX, m_DimY, m_DimZ, m_VoxSizeX, m_VoxSizeY, m_VoxSizeZ, mp_Scanner->GetScannerModel(), 1., 0., "INF" ))
      {
        LogCerr ("***** oSimulator::Project() -> An error occured while saving the mumap to be projected !" << endl);
        return 1;
      }
    }
*/
  }

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Clock start
  clock_t clock_start = clock();
  time_t time_start = time(NULL);

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Allocate output sinograms
  if (m_Verbose>=1)
  {
    LogCout ("  --> Sinogram dimensions: [" << m_NbElem << ";" << m_NbView << ";" << m_NbSino << "]" << endl);
    LogCout ("  --> Span: " << m_Span << " | Maximum ring difference: " << m_MaxRingDiff << " | Mash: " << m_Mash << endl);
  }
  mp_SinoForw = (float***)malloc(m_NbSino*sizeof(float**));
  for (int s=0; s<m_NbSino; s++)
  {
    mp_SinoForw[s] = (float**)malloc(m_NbView*sizeof(float*));
    for (int v=0; v<m_NbView; v++) mp_SinoForw[s][v] = (float*)calloc(m_NbElem,sizeof(float));
  }
  mp_SinoNorm = (float***)malloc(m_NbSino*sizeof(float**));
  for (int s=0; s<m_NbSino; s++)
  {
    mp_SinoNorm[s] = (float**)malloc(m_NbView*sizeof(float*));
    for (int v=0; v<m_NbView; v++) mp_SinoNorm[s][v] = (float*)calloc(m_NbElem,sizeof(float));
  }
  if (m_AttnCorr)
  {
    mp_SinoAttn = (float***)malloc(m_NbSino*sizeof(float**));
    for (int s=0; s<m_NbSino; s++)
    {
      mp_SinoAttn[s] = (float**)malloc(m_NbView*sizeof(float*));
      for (int v=0; v<m_NbView; v++) mp_SinoAttn[s][v] = (float*)calloc(m_NbElem,sizeof(float));
    }
  }

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Allocate fan sum
  mp_FanSum = (float*)calloc(mp_Scanner->GetNbTotalCrystals(),sizeof(float));

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Allocate buffers and counters
  unsigned long int* nb_useful_lors = (unsigned long int*)calloc(m_NbThreads,sizeof(unsigned long int));

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Get some usefull stuff
  int scanner_model = mp_Scanner->GetScannerModel();
  int nb_total_trans_crystals = mp_Scanner->GetNbTotalTransCrystals();
  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Set number of threads
  #ifdef OMP_MODE
  omp_set_num_threads(m_NbThreads);
  #endif
  // Allocate buffers for crystals and rings IDs (multi-thread)
  int **crystals1 = (int**)malloc(m_NbThreads*sizeof(int*));
  int **crystals2 = (int**)malloc(m_NbThreads*sizeof(int*));
  for (int t=0; t<m_NbThreads; t++)
  {
    crystals1[t] = (int*)malloc(m_Mash*sizeof(int));
    crystals2[t] = (int*)malloc(m_Mash*sizeof(int));
  }
  unsigned int** rings1 = (unsigned int**)malloc(m_NbThreads*sizeof(unsigned int*));
  unsigned int** rings2 = (unsigned int**)malloc(m_NbThreads*sizeof(unsigned int*));
  int max_nb_ring_pairs = m_Span/2 + 1;
  for (int t=0; t<m_NbThreads; t++)
  {
    rings1[t] = (unsigned int*)malloc(max_nb_ring_pairs*sizeof(unsigned int));
    rings2[t] = (unsigned int*)malloc(max_nb_ring_pairs*sizeof(unsigned int));
  }
  int nb_axial_crystals = mp_Scanner->GetNbAxialCrystals();
  int nb_trans_crystals_in_head = mp_Scanner->GetNbTransCrystalsInHead();

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Loop on sinogram bins using openMP
  if (m_Verbose>=1) LogCout ("  --> Project true component ... " << endl);
  int sino_index;
  #ifdef OMP_MODE
  #pragma omp parallel for private(sino_index) schedule(static, 1)
  #endif
  for (sino_index=0; sino_index<m_NbSino; sino_index++)
  {
    // Get the thread number
    int t = 0;
    #ifdef OMP_MODE
    t = omp_get_thread_num();
    #endif

    // Verbose
    if (m_Verbose>=1 && t==0)
    {
      float percent = ((float)sino_index)*100./((float)m_NbSino);
      cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
           << "      " << percent << " %              " << flush;
    }

    // Get the number of ring pairs for this sinogram index
    int nb_ring_pairs;
    if (scanner_model==SCANNER_HRPLUS)
    {
      nb_ring_pairs = mp_TableHRplus->GetNbRingPairsBySinoIndex(sino_index);
      mp_TableHRplus->GetRingPairsBySinoIndex(sino_index,rings1[t],rings2[t]);
    }
    else if (scanner_model==SCANNER_INVEON)
    {
      nb_ring_pairs = mp_TableInveon->GetNbRingPairsBySinoIndex(sino_index);
      mp_TableInveon->GetRingPairsBySinoIndex(sino_index,rings1[t],rings2[t]);
    }
    else if (scanner_model==SCANNER_HRRT)
    {
      // To be implemented
      cerr << endl << endl << "*!*!*!*!*  Part of code not yet implemented  *!*!*!*!*" << endl << endl;
      exit(-245);
    }
    else if (scanner_model==SCANNER_BIOGRAPH)
    {
      nb_ring_pairs = mp_TableBiograph->GetNbRingPairsBySinoIndex(sino_index);
      mp_TableBiograph->GetRingPairsBySinoIndex(sino_index,rings1[t],rings2[t]);
      // Compute indices without gaps (no need to check wether we are in a gap because gaps are already removed from the span table)
      for (int r=0; r<nb_ring_pairs; r++)
      {
        rings1[t][r] -= rings1[t][r]/(nb_axial_crystals+1);
        rings2[t][r] -= rings2[t][r]/(nb_axial_crystals+1);
      }
    }
    else if (scanner_model==SCANNER_BIOGRAPH2D || scanner_model==SCANNER_MMR2D)
    {
      nb_ring_pairs = 1;
      rings1[t][0] = 0;
      rings2[t][0] = 0;
    }

    // Compute the number of LORs
    int nb_lors = nb_ring_pairs * m_Mash;

    // Loops on view and elem
    for (int view_index=0; view_index<m_NbView; view_index++) for (int elem_index=0; elem_index<m_NbElem; elem_index++)
    {
      // The the transaxial crystal indices taking mashing into account
      if (scanner_model==SCANNER_HRPLUS)
      {
        // Crystal indices
        mp_TableHRplus->GetCrystalIDsFromElemView(elem_index, view_index, crystals1[t], crystals2[t]);
      }
      else if (scanner_model==SCANNER_INVEON)
      {
        // Crystal indices
        mp_TableInveon->GetCrystalIDsFromElemView(elem_index, view_index, crystals1[t], crystals2[t]);
      }
      else if (scanner_model==SCANNER_HRRT)
      {
        // To be implemented
        cerr << endl << endl << "*!*!*!*!*  Part of code not yet implemented  *!*!*!*!*" << endl << endl;
        exit(-245);
      }
      else if (scanner_model==SCANNER_BIOGRAPH)
      {
        // Crystal indices
        mp_TableBiograph->GetCrystalIDsFromElemView(elem_index, view_index, crystals1[t], crystals2[t]);
        // Check if we are in a gap, we skip it
        if ((crystals1[t][0]+1)%(nb_trans_crystals_in_head+1)==0) continue;
        if ((crystals2[t][0]+1)%(nb_trans_crystals_in_head+1)==0) continue;
        // Compute the crystal indices on ring without gaps
        crystals1[t][0] -= crystals1[t][0]/(nb_trans_crystals_in_head+1);
        crystals2[t][0] -= crystals2[t][0]/(nb_trans_crystals_in_head+1);
      }
      else if (scanner_model==SCANNER_BIOGRAPH2D)
      {
        // Crystal indices
        mp_TableBiograph2D->GetCrystalIDsFromElemView(elem_index, view_index, crystals1[t], crystals2[t]);
        // Check if we are in a gap, we skip it
        if ((crystals1[t][0]+1)%(nb_trans_crystals_in_head+1)==0) continue;
        if ((crystals2[t][0]+1)%(nb_trans_crystals_in_head+1)==0) continue;
        // Compute the crystal indices on ring without gaps
        crystals1[t][0] -= crystals1[t][0]/(nb_trans_crystals_in_head+1);
        crystals2[t][0] -= crystals2[t][0]/(nb_trans_crystals_in_head+1);
      }
      else if (scanner_model==SCANNER_MMR2D)
      {
        // Crystal indices
        mp_TableMMR2D->GetCrystalIDsFromElemView(elem_index, view_index, crystals1[t], crystals2[t]);
        // Check if we are in a gap, we skip it
        if ((crystals1[t][0]+1)%(nb_trans_crystals_in_head+1)==0) continue;
        if ((crystals2[t][0]+1)%(nb_trans_crystals_in_head+1)==0) continue;
        // Compute the crystal indices on ring without gaps
        crystals1[t][0] -= crystals1[t][0]/(nb_trans_crystals_in_head+1);
        crystals2[t][0] -= crystals2[t][0]/(nb_trans_crystals_in_head+1);
      }
      // Loop on all ring pairs
      for (int ring=0; ring<nb_ring_pairs; ring++)
      {
        // Loop on all mashing LORs
        for (int mash=0; mash<m_Mash; mash++)
        {
          // Compute the global crystal IDs
          int id1 = rings1[t][ring]*nb_total_trans_crystals + crystals1[t][mash];
          int id2 = rings2[t][ring]*nb_total_trans_crystals + crystals2[t][mash];
          // Get crystal efficiencies
          float eff1 = mp_Scanner->GetEfficiency(id1);
          float eff2 = mp_Scanner->GetEfficiency(id2);
          float efficiency = eff1 * eff2;
          if (efficiency<=0.) continue;
          // Get crystal corners' coordinates (for the 2 given)
          float cx1 = mp_Scanner->GetCornerX1(id1);
          float cx2 = mp_Scanner->GetCornerX2(id1);
          float cy1 = mp_Scanner->GetCornerY1(id1);
          float cy2 = mp_Scanner->GetCornerY2(id1);
          float cz1 = mp_Scanner->GetCornerZ1(id1);
          float cz2 = mp_Scanner->GetCornerZ2(id1);
          float cx3 = mp_Scanner->GetCornerX1(id2);
          float cx4 = mp_Scanner->GetCornerX2(id2);
          float cy3 = mp_Scanner->GetCornerY1(id2);
          float cy4 = mp_Scanner->GetCornerY2(id2);
          float cz3 = mp_Scanner->GetCornerZ1(id2);
          float cz4 = mp_Scanner->GetCornerZ2(id2);
          // Compute crystals positions
          float x1 = (cx1 + cx2) / 2.;
          float y1 = (cy1 + cy2) / 2.;
          float z1 = (cz1 + cz2) / 2.;
          float x2 = (cx3 + cx4) / 2.;
          float y2 = (cy3 + cy4) / 2.;
          float z2 = (cz3 + cz4) / 2.;
          // Apply offsets
          x1 += m_OffsetX;
          x2 += m_OffsetX;
          y1 += m_OffsetY;
          y2 += m_OffsetY;
          z1 += m_OffsetZ;
          z2 += m_OffsetZ;
          // Update normalization factor
          mp_SinoNorm[sino_index][view_index][elem_index] += efficiency;
          // Project using Siddon
          float emission, transmission;
          if (f_Projector==0) SiddonForwardProjection( x1, y1, z1, x2, y2, z2, mp_ProjectImage, mp_MuMapImage, m_DimX, m_DimY, m_DimZ, m_VoxSizeX, m_VoxSizeY, m_VoxSizeZ, &emission, &transmission);
          else if (f_Projector==1) SiddonDidierForwardProjection( x1, y1, z1, x2, y2, z2, mp_ProjectImage, mp_MuMapImage, m_DimX, m_DimY, m_DimZ, m_VoxSizeX, m_VoxSizeY, m_VoxSizeZ, &emission, &transmission);
          float value = efficiency * emission / transmission;
          // Register sinogram attenuation
          if (m_AttnCorr) mp_SinoAttn[sino_index][view_index][elem_index] += transmission;
          // Add contribution to forward sinogram
          mp_SinoForw[sino_index][view_index][elem_index] += value;
          // Add contribution to fan sum
          mp_FanSum[id1] += value;
          mp_FanSum[id2] += value;
          // Increment counter
          nb_useful_lors[t]++;
        } // End loop on all mashing LORs
      } // End loop on all ring pairs

      // Mean the attenuation sinogram bin
      if (m_AttnCorr) mp_SinoAttn[sino_index][view_index][elem_index] /= ((float)nb_lors);

      // Inverse normalization factor
      if (mp_SinoNorm[sino_index][view_index][elem_index]>0.) mp_SinoNorm[sino_index][view_index][elem_index] = 1./mp_SinoNorm[sino_index][view_index][elem_index];
      else mp_SinoNorm[sino_index][view_index][elem_index] = 0.;

    } // End loop on view and elem
  }
  // Verbose
  if (m_Verbose>=1) cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
                         << "      100 %              " << endl;

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Add lor counters
  for (int t=1; t<m_NbThreads; t++) nb_useful_lors[0] += nb_useful_lors[t];
  if (m_Verbose>=1) LogCout ("  --> Useful projected LORs: " << nb_useful_lors[0] << endl);

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Clock stop
  clock_t clock_stop = clock();
  time_t time_stop = time(NULL);
  if (m_Verbose>=1)
  {
    LogCout ("  --> Time spent | User: " << time_stop-time_start << " sec | CPU: " << (clock_stop-clock_start)/((float)CLOCKS_PER_SEC) << " sec" << endl);
  }

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Free some buffers
  free(nb_useful_lors);

  // End
  return 0;
}

// ==========================================================================================================================================
// Function ApplyCounts
//   This function applies scatter and random components and computes accurate total number of counts
// ==========================================================================================================================================
int oSimulator::ApplyCounts(float f_ScatterFraction, float f_RandomFraction, float f_RandomFractionLSO, long int f_NbCounts, float f_ECF, bool f_ListMode, bool f_FillEqualLORs, int f_NbReplicates)
{
  // Verbose
  if (m_Verbose>=1) LogCout ("oSimulator::ApplyCounts() -> Apply scatter, random and poisson effects" << endl);

  // Check if list-mode must give counts
  if (f_ListMode && f_NbCounts<=0 && f_ECF<=0.)
  {
    LogCerr ("***** oSimulator::ApplyCounts() -> Cannot perform list-mode simulation if no counts or ECF are given !" << endl);
    return 1;
  }

  // Set the number of replicates
  m_NbReplicates = f_NbReplicates;

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // True part

  // Verbose
  if (m_Verbose>=1) LogCout ("  --> Compute true component ..." << endl);

  // Divide true sinogram by sinogram bin size
  if (mp_Scanner->GetScannerModel()==SCANNER_HRPLUS)
  {
    if (m_Verbose>=2) LogCout ("      Divide by radial bin size" << endl);
    for (int s=0; s<m_NbSino; s++) for (int v=0; v<m_NbView; v++) for (int e=0; e<m_NbElem; e++)
      mp_SinoForw[s][v][e] /= mp_Scanner->GetRadialBinSize();
  }

  // Compute total true
  if (m_Verbose>=2) LogCout ("      Compute total true" << endl);
  m_TotalTrue = 0.;
  for (int s=0; s<m_NbSino; s++) for (int v=0; v<m_NbView; v++) for (int e=0; e<m_NbElem; e++)
  {
    m_TotalTrue += (double)mp_SinoForw[s][v][e];
  }

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Scatter part
  m_TotalScat = 0.;
  if (f_ScatterFraction>0.)
  {
    // Check that scatter fraction is less than 1
    if (f_ScatterFraction>=1.)
    {
      LogCerr ("***** oSimulator::ApplyCounts() -> Scatter fraction must be in [0.;1.[ !" << endl);
      return 1;
    }

    // Verbose
    if (m_Verbose>=1) LogCout ("  --> Compute scatter component ..." << endl);

    // Allocate sinogram and copy trues in it
    if (m_Verbose>=2) LogCout ("      Allocate sinogram" << endl);
    mp_SinoScat = (float***)malloc(m_NbSino*sizeof(float**));
    for (int s=0; s<m_NbSino; s++)
    {
      mp_SinoScat[s] = (float**)malloc(m_NbView*sizeof(float*));
      for (int v=0; v<m_NbView; v++)
      {
        mp_SinoScat[s][v] = (float*)malloc(m_NbElem*sizeof(float));
        for (int e=0; e<m_NbElem; e++) mp_SinoScat[s][v][e] = mp_SinoForw[s][v][e];
      }
    }
    // Correct for normalization
    if (m_Verbose>=2) LogCout ("      Correct for normalization" << endl);
    for (int s=0; s<m_NbSino; s++) for (int v=0; v<m_NbView; v++) for (int e=0; e<m_NbElem; e++)
      mp_SinoScat[s][v][e] *= mp_SinoNorm[s][v][e];
    // Convolve
    if (m_Verbose>=2) LogCout ("      Convolve" << endl);
    ConvolveScatterComponent();
    // Affect by normalization
    if (m_Verbose>=2) LogCout ("      Affect by normalization" << endl);
    for (int s=0; s<m_NbSino; s++) for (int v=0; v<m_NbView; v++) for (int e=0; e<m_NbElem; e++)
    {
      if (mp_SinoNorm[s][v][e]>0.) mp_SinoScat[s][v][e] /= mp_SinoNorm[s][v][e];
      else mp_SinoScat[s][v][e] = 0.;
    }
    // Compute total scatter and set to desired scale
    if (m_Verbose>=2) LogCout ("      Scale to desired scatter fraction (" << f_ScatterFraction << ")" << endl);
    double tmp_total_scat = 0.;
    for (int s=0; s<m_NbSino; s++) for (int v=0; v<m_NbView; v++) for (int e=0; e<m_NbElem; e++)
      tmp_total_scat += (double)mp_SinoScat[s][v][e];
    m_TotalScat = m_TotalTrue * ((double)f_ScatterFraction) / (1.-((double)f_ScatterFraction));
    double ratio = m_TotalScat / tmp_total_scat;
    for (int s=0; s<m_NbSino; s++) for (int v=0; v<m_NbView; v++) for (int e=0; e<m_NbElem; e++)
      mp_SinoScat[s][v][e] = (float)(((double)mp_SinoScat[s][v][e])*ratio);

    // Valid
    m_ScatCorr = true;
  }

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Random part
  m_TotalRand = 0.;
  if (f_RandomFraction>0.)
  {
    // Check that random fraction is less than 1
    if (f_RandomFraction>=1.)
    {
      LogCerr ("***** oSimulator::ApplyCounts() -> Random fraction must be in [0.;1.[ !" << endl);
      return 1;
    }
    // Check that LSO random fraction is between 0 and 1
    if (f_RandomFractionLSO<0. || f_RandomFractionLSO>1.)
    {
      LogCerr ("***** oSimulator::ApplyCounts() -> LSO random fraction must be in [0.;1.] !" << endl);
      return 1;
    }

    // Verbose
    if (m_Verbose>=1) LogCout ("  --> Compute random component ..." << endl);

    // Allocate sinogram
    if (m_Verbose>=2) LogCout ("      Allocate sinogram" << endl);
    mp_SinoRand = (float***)malloc(m_NbSino*sizeof(float**));
    for (int s=0; s<m_NbSino; s++)
    {
      mp_SinoRand[s] = (float**)malloc(m_NbView*sizeof(float*));
      for (int v=0; v<m_NbView; v++) mp_SinoRand[s][v] = (float*)calloc(m_NbElem,sizeof(float));
    }

    // Fan sum
    if (m_Verbose>=2) LogCout ("      Compute random from fan sum" << endl);
    ComputeRandomFanSum();

    // Add LSO flat contribution
    if (m_Verbose>=2) LogCout ("      Add LSO flat contribution (" << f_RandomFractionLSO*100. << "% of total randoms)" << endl);
    double tmp_total_rand = 0.;
    unsigned long int nb_contributing_bins = 0;
    for (int s=0; s<m_NbSino; s++) for (int v=0; v<m_NbView; v++) for (int e=0; e<m_NbElem; e++)
    {
      if (mp_SinoNorm[s][v][e]>0.)
      {
        tmp_total_rand += (double)mp_SinoRand[s][v][e];
        nb_contributing_bins++;
      }
    }
    // The following is the contribution of LSO for each sinogram bin
    float lso_random_per_bin = tmp_total_rand * f_RandomFractionLSO / ((1.-f_RandomFractionLSO) * ((float)nb_contributing_bins));
    for (int s=0; s<m_NbSino; s++) for (int v=0; v<m_NbView; v++) for (int e=0; e<m_NbElem; e++)
    {
      if (mp_SinoNorm[s][v][e]>0.) mp_SinoRand[s][v][e] += lso_random_per_bin;
    }
    
    // Compute total random and set to desired scale
    if (m_Verbose>=2) LogCout ("      Scale to desired random fraction (" << f_RandomFraction << ")" << endl);
    tmp_total_rand = 0.;
    for (int s=0; s<m_NbSino; s++) for (int v=0; v<m_NbView; v++) for (int e=0; e<m_NbElem; e++)
      tmp_total_rand += (double)mp_SinoRand[s][v][e];
    m_TotalRand = (m_TotalTrue+m_TotalScat) * ((double)f_RandomFraction) / (1.-((double)f_RandomFraction));
    double ratio = m_TotalRand / tmp_total_rand;
    for (int s=0; s<m_NbSino; s++) for (int v=0; v<m_NbView; v++) for (int e=0; e<m_NbElem; e++)
      mp_SinoRand[s][v][e] = (float)(((double)mp_SinoRand[s][v][e])*ratio);

    // Valid
    m_RandCorr = true;
  }

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Verbose total after projection

  // Verbose
  if (m_Verbose>=1)
  {
    LogCout ("  --> Total projections:" << endl);
    LogCout ("      Trues   : " << m_TotalTrue << endl);
    LogCout ("      Scatters: " << m_TotalScat << endl);
    LogCout ("      Randoms : " << m_TotalRand << endl);
  }

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Prompt sinogram

  // Verbose
  if (m_Verbose>=1)
  {
    if (f_ListMode) LogCout ("  --> Compute list-mode ..." << endl);
    else LogCout ("  --> Compute prompt component ..." << endl);
  }

  // Allocate final sinogram
  if (m_Verbose>=2) LogCout ("      Allocate prompt sinogram" << endl);
  mp_SinoPrompt = (short int***)malloc(m_NbSino*sizeof(short int**));
  for (int s=0; s<m_NbSino; s++)
  {
    mp_SinoPrompt[s] = (short int**)malloc(m_NbView*sizeof(short int*));
    for (int v=0; v<m_NbView; v++)
    {
      mp_SinoPrompt[s][v] = (short int*)calloc(m_NbElem,sizeof(short int));
    }
  }

  // Reset counts
  m_CountTrue = 0;
  m_CountScat = 0;
  m_CountRand = 0;
  m_CountPrompt = 0;

  // _____________________________________________________________________________________________________________________________________
  // List-Mode
  if (f_ListMode)
  {
    // Compute total prompt and search maximum
    if (m_Verbose>=2) LogCout ("      Compute total prompt" << flush);
    double tmp_total_prompt = 0.;
    for (int s=0; s<m_NbSino; s++) for (int v=0; v<m_NbView; v++) for (int e=0; e<m_NbElem; e++)
    {
      double prompt = (double)mp_SinoForw[s][v][e];
      if (m_ScatCorr) prompt += (double)mp_SinoScat[s][v][e];
      if (m_RandCorr) prompt += (double)mp_SinoRand[s][v][e];
      tmp_total_prompt += prompt;
    }
    LogCout (": " << tmp_total_prompt << endl);

    // Determine ECF and scale each sinogram and calculate sum of trues, scat and rand
    if (f_ECF<=0.)
    {
      LogCout ("      Number of counts provided, so compute associated ECF" << endl);
      m_ECF = tmp_total_prompt / ((double)f_NbCounts);
    }
    else if (f_NbCounts<=0)
    {
      LogCout ("      ECF provided, so compute associated number of counts" << endl);
      m_ECF = f_ECF;
      f_NbCounts = ((long int)( tmp_total_prompt / m_ECF ));
    }
    else
    {
      LogCerr ("***** oSimulator::ApplyCounts() -> Internal problem, reason unknown ... Sorry." << endl);
      return 1;
    }
    if (m_ECF==0.)
    {
      LogCerr ("***** oSimulator::ApplyCounts() -> Total prompt data is null ! Emission phantom is probably empty... Abort." << endl);
      return 1;
    }
    double tmp_total_true = 0.;
    double tmp_total_scat = 0.;
    double tmp_total_rand = 0.;
    if (m_Verbose>=2) LogCout ("      Scale sinograms (ECF: " << m_ECF << " for total counts of " << f_NbCounts << ")" << endl);
    for (int s=0; s<m_NbSino; s++) for (int v=0; v<m_NbView; v++) for (int e=0; e<m_NbElem; e++)
    {
      // Trues
      mp_SinoForw[s][v][e] = (float)(((double)mp_SinoForw[s][v][e])/m_ECF);
      tmp_total_true += ((double)mp_SinoForw[s][v][e]);
      // Scatters
      if (m_ScatCorr)
      {
        mp_SinoScat[s][v][e] = (float)(((double)mp_SinoScat[s][v][e])/m_ECF);
        tmp_total_scat += ((double)mp_SinoScat[s][v][e]);
      }
      // Randoms
      if (m_RandCorr)
      {
        mp_SinoRand[s][v][e] = (float)(((double)mp_SinoRand[s][v][e])/m_ECF);
        tmp_total_rand += ((double)mp_SinoRand[s][v][e]);
      }
    }

    // Compute trues PDF for each sinogram and view (for optimization)
    double*  sino_pdf_true = NULL;
    double** view_pdf_true = NULL;
    if (m_Verbose>=2) LogCout ("      Compute trues PDF for optimization" << endl);
    sino_pdf_true = (double*) calloc(((unsigned short int)m_NbSino),sizeof(double));
    view_pdf_true = (double**)malloc(m_NbSino*sizeof(double*));
    for (int s=0; s<m_NbSino; s++)
    {
      view_pdf_true[s] = (double*)calloc(m_NbView,sizeof(double));
      for (int v=0; v<m_NbView; v++)
      {
        for (int e=0; e<m_NbElem; e++) view_pdf_true[s][v] += ((double)mp_SinoForw[s][v][e]);
        sino_pdf_true[s] += view_pdf_true[s][v];
      }
    }

    // Compute scatters PDF for each sinogram and view (for optimization)
    double*  sino_pdf_scat = NULL;
    double** view_pdf_scat = NULL;
    if (m_ScatCorr)
    {
      if (m_Verbose>=2) LogCout ("      Compute scatters PDF for optimization" << endl);
      sino_pdf_scat = (double*) calloc(m_NbSino,sizeof(double));
      view_pdf_scat = (double**)malloc(m_NbSino*sizeof(double*));
      for (int s=0; s<m_NbSino; s++)
      {
        view_pdf_scat[s] = (double*)calloc(m_NbView,sizeof(double));
        for (int v=0; v<m_NbView; v++)
        {
          for (int e=0; e<m_NbElem; e++) view_pdf_scat[s][v] += ((double)mp_SinoScat[s][v][e]);
          sino_pdf_scat[s] += view_pdf_scat[s][v];
        }
      }
    }

    // Compute randoms PDF for each sinogram and view (for optimization)
    double*  sino_pdf_rand = NULL;
    double** view_pdf_rand = NULL;
    if (m_RandCorr)
    {
      if (m_Verbose>=2) LogCout ("      Compute scatters PDF for optimization" << endl);
      sino_pdf_rand = (double*) calloc(m_NbSino,sizeof(double));
      view_pdf_rand = (double**)malloc(m_NbSino*sizeof(double*));
      for (int s=0; s<m_NbSino; s++)
      {
        view_pdf_rand[s] = (double*)calloc(m_NbView,sizeof(double));
        for (int v=0; v<m_NbView; v++)
        {
          for (int e=0; e<m_NbElem; e++) view_pdf_rand[s][v] += ((double)mp_SinoRand[s][v][e]);
          sino_pdf_rand[s] += view_pdf_rand[s][v];
        }
      }
    }

    // Compute relative event fractions
    double relative_random_fraction = f_RandomFraction;
    double relative_scatter_fraction = f_RandomFraction + (1.-f_RandomFraction)*f_ScatterFraction;

    // Multi-thread sinogram bin list
    long int* thread_nb_events = (long int*)malloc(m_NbThreads*sizeof(long int));
    int** thread_elem_bin = (int**)malloc(m_NbThreads*sizeof(int*));
    int** thread_view_bin = (int**)malloc(m_NbThreads*sizeof(int*));
    int** thread_sino_bin = (int**)malloc(m_NbThreads*sizeof(int*));
    #ifdef OMP_MODE
    omp_set_num_threads(m_NbThreads);
    #endif

    // ----------------------------------------------------------------
    // Start the replicates loop here
    for (int replicate=1; replicate<=m_NbReplicates; replicate++)
    {
      // Replicate number as a string
      char tmp_rep[100]; sprintf(tmp_rep,"%d",replicate);
      string str_rep = (string)tmp_rep;

      // Compute poisson total number of counts
      long int real_number_of_counts = PoissonSampleBig( ((float)f_NbCounts) );
      long int real_number_of_counts_per_thread = real_number_of_counts / ((long int)m_NbThreads) + 1;
      m_CountPrompt = real_number_of_counts;
      if (m_Verbose>=2) LogCout ("      Compute list-mode #" << replicate << " with " << real_number_of_counts << " events ..." << endl);

/*
Verifier les randoms multi-thread.
Verifier que les comptes sont bons sur les trues+randoms+scatters à tous les niveaux.
Verifier que le multi-bed dans la recon donne bien exactement les memes offset que dans les weights et dans la simu.
*/

      // Reset global counters
      m_CountScat = 0;
      m_CountRand = 0;
      m_CountTrue = 0;

      // Thread safe counters
      int* thread_true_counter = (int*)calloc(m_NbThreads,sizeof(int));
      int* thread_scat_counter = (int*)calloc(m_NbThreads,sizeof(int));
      int* thread_rand_counter = (int*)calloc(m_NbThreads,sizeof(int));

      // Allocate thread buffers and reset thread counters
      for (int th=0; th<m_NbThreads; th++)
      {
        thread_elem_bin[th] = (int*)calloc(real_number_of_counts_per_thread,sizeof(int));
        thread_view_bin[th] = (int*)calloc(real_number_of_counts_per_thread,sizeof(int));
        thread_sino_bin[th] = (int*)malloc(real_number_of_counts_per_thread*sizeof(int));
        for (int e=0; e<real_number_of_counts_per_thread; e++) thread_sino_bin[th][e] = -1;
        thread_nb_events[th] = 0;
      }

/*
ofstream** fout = (ofstream**)malloc(m_NbThreads*sizeof(ofstream*));
for (int th=0; th<m_NbThreads; th++)
{
  char tmp_th[100]; sprintf(tmp_th,"%d",th);
  string file_name = "seeds_from_thread"+((string)tmp_th)+".txt";
  fout[th] = new ofstream(file_name.c_str());
}
*/

      // Will precompute random numbers to be thread safe
      int* rand_type_shoot = (int*)malloc(real_number_of_counts*sizeof(int));
      int* rand_sino_shoot = (int*)malloc(real_number_of_counts*sizeof(int));
      for (int e=0; e<real_number_of_counts; e++)
      {
        rand_type_shoot[e] = rand();
        rand_sino_shoot[e] = rand();
      }

      // Loop to create list-mode
      long int event;
      long int event_th0 = 0;
      #pragma omp parallel for private(event) schedule(static, 1)
      for (event=0; event<real_number_of_counts; event++)
      {
        // Get the thread number
        int th = 0;
        #ifdef OMP_MODE
        th = omp_get_thread_num();
        #endif
        // Verbose
        if (th==0 && m_Verbose>=1)
        {
          if (event_th0==1000)
          {
            float percent = ((float)(event))*100./((float)(real_number_of_counts));
            cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
            cout << "          " << percent << " %                 " << flush;
            event_th0 = 0;
          }
          else event_th0++;
        }
        // Shoot to choose between true, scatter and random event
        double type_shoot = ((double)rand_type_shoot[event])/((double)RAND_MAX);
        // Shoot to select in the sinograms PDF
        double sino_shoot = ((double)rand_sino_shoot[event])/((double)RAND_MAX);
//*(fout[th]) << event << "\t" << type_shoot << "\t" << sino_shoot << endl;
        // Set type
        int event_type = EVENT_TRUE;
        if (type_shoot <= relative_random_fraction) event_type = EVENT_RAND;
        else if (type_shoot <= relative_scatter_fraction) event_type = EVENT_SCAT;
        // Sinogram to shoot in
        float*** sinogram = NULL;
        double*  sino_pdf = NULL;
        double** view_pdf = NULL;
        // Switch on event type
        switch (event_type)
        {
          case EVENT_TRUE:
            sinogram = mp_SinoForw;
            thread_true_counter[th]++;
            sino_shoot *= tmp_total_true;
            sino_pdf = sino_pdf_true;
            view_pdf = view_pdf_true;
            break;
          case EVENT_RAND:
            sinogram = mp_SinoRand;
            thread_rand_counter[th]++;
            sino_shoot *= tmp_total_rand;
            sino_pdf = sino_pdf_rand;
            view_pdf = view_pdf_rand;
            break;
          case EVENT_SCAT:
            sinogram = mp_SinoScat;
            thread_scat_counter[th]++;
            sino_shoot *= tmp_total_scat;
            sino_pdf = sino_pdf_scat;
            view_pdf = view_pdf_scat;
            break;
        }
        // Search the bin corresponding to the shot value
        double tmp_pdf = 0.;
        bool break_all = false;
        for (int s=0; s<m_NbSino; s++)
        {
          // Test if the shoot is in this sinogram
          if (sino_shoot-tmp_pdf > sino_pdf[s])
          {
            tmp_pdf += sino_pdf[s];
            continue;
          }
          // Loop on views
          for (int v=0; v<m_NbView; v++)
          {
            // Test if the shoot is in this view
            if (sino_shoot-tmp_pdf > view_pdf[s][v])
            {
              tmp_pdf += view_pdf[s][v];
              continue;
            }
            // Loop on elems
            for (int e=0; e<m_NbElem; e++)
            {
              // Update PDF summation
              tmp_pdf += sinogram[s][v][e];
              // Compare
              if (sino_shoot<tmp_pdf)
              {
                thread_elem_bin[th][thread_nb_events[th]] = e;
                thread_view_bin[th][thread_nb_events[th]] = v;
                thread_sino_bin[th][thread_nb_events[th]] = s;
                break_all = true;
                break;
              }
              // At this step, if this is the last elem, this means that we must choose it but we passed because of rounding errors, so we force it
              if (e==m_NbElem-1)
              {
                thread_elem_bin[th][thread_nb_events[th]] = e;
                thread_view_bin[th][thread_nb_events[th]] = v;
                thread_sino_bin[th][thread_nb_events[th]] = s;
                break_all = true;
                break;
              }
            }
            if (break_all) break;
          }
          if (break_all) break;
        }
        // Increment number of events
        thread_nb_events[th]++;
      }
      if (m_Verbose>=1)
       {
       cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
        cout << "          100 %                       " << endl;
      }

      // Free thread-safe tables for random numbers
      free(rand_type_shoot);
      free(rand_sino_shoot);

      // Add thread-safe counters
      for (int th=0; th<m_NbThreads; th++)
      {
        m_CountTrue += thread_true_counter[th];
        m_CountScat += thread_scat_counter[th];
        m_CountRand += thread_rand_counter[th];
      }

      // Free thread-safe counters
      free(thread_true_counter);
      free(thread_scat_counter);
      free(thread_rand_counter);

      // Open list-mode file
      string baseName = oOutputManager::GetInstance()->GetRootName();
      string pathName = oOutputManager::GetInstance()->GetPathName();
      string file_listmode = "";
      if (m_NbReplicates==1) file_listmode = pathName + baseName+"/"+baseName+"_lm.elm";
      else file_listmode = pathName + baseName+"/"+baseName+"_lm_rep"+str_rep+".elm";
      FILE* flm = fopen(file_listmode.c_str(),"wb");
      if (flm==NULL)
      {
        LogCerr ("***** oSimulator::ApplyCounts() -> Failed to create output data list-mode file '" << file_listmode << "' !" << endl);
        return 1;
      }
      if (m_Verbose>=2) LogCout ("      Write list-mode file '" << file_listmode << "' ..." << endl);

      // Buffers
      int *crystals1 = (int*)malloc(m_Mash*sizeof(int));
      int *crystals2 = (int*)malloc(m_Mash*sizeof(int));
      unsigned int *rings1 = (unsigned int*)malloc((m_Span/2+1)*sizeof(unsigned int));
      unsigned int *rings2 = (unsigned int*)malloc((m_Span/2+1)*sizeof(unsigned int));
      unsigned short int max_nb_lors = m_Mash * (m_Span/2+1);
      // Get some usefull stuff
      int scanner_model = mp_Scanner->GetScannerModel();
      int nb_total_trans_crystals = mp_Scanner->GetNbTotalTransCrystals();
      int nb_axial_crystals = mp_Scanner->GetNbAxialCrystals();
      int nb_trans_crystals_in_head = mp_Scanner->GetNbTransCrystalsInHead();
      // Write list-mode
      long int total_event_print = 0;
      for (long int th=0, total_event=0; th<m_NbThreads; th++) for (int event=0; event<thread_nb_events[th]; event++, total_event++)
      {
        int nb_data_written = 0;
        // Verbose
        if (m_Verbose>=1)
        {
          if (total_event_print==1000)
          {
            float percent = ((float)(total_event))*100./((float)(real_number_of_counts));
            cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
            cout << "          " << percent << " %                 " << flush;
            total_event_print = 0;
          }
          else total_event_print++;
        }
        // Elem, view and sino
        int elem_bin = thread_elem_bin[th][event];
        int view_bin = thread_view_bin[th][event];
        int sino_bin = thread_sino_bin[th][event];
        // Increment the prompt sinogram corresponding bin
        mp_SinoPrompt[sino_bin][view_bin][elem_bin]++;
        // Get the number of ring pairs for this sinogram index and the crystal indices
        int nb_ring_pairs;
        if (scanner_model==SCANNER_HRPLUS)
        {
          nb_ring_pairs = mp_TableHRplus->GetNbRingPairsBySinoIndex(sino_bin);
          mp_TableHRplus->GetRingPairsBySinoIndex(sino_bin,rings1,rings2);
          mp_TableHRplus->GetCrystalIDsFromElemView(elem_bin, view_bin, crystals1, crystals2);
        }
        else if (scanner_model==SCANNER_INVEON)
        {
          nb_ring_pairs = mp_TableInveon->GetNbRingPairsBySinoIndex(sino_bin);
          mp_TableInveon->GetRingPairsBySinoIndex(sino_bin,rings1,rings2);
          mp_TableInveon->GetCrystalIDsFromElemView(elem_bin, view_bin, crystals1, crystals2);
        }
        else if (scanner_model==SCANNER_HRRT)
        {
          // To be implemented
          cerr << endl << endl << "*!*!*!*!*  Part of code not yet implemented  *!*!*!*!*" << endl << endl;
          exit(-245);
        }
        else if (scanner_model==SCANNER_BIOGRAPH)
        {
          nb_ring_pairs = mp_TableBiograph->GetNbRingPairsBySinoIndex(sino_bin);
          mp_TableBiograph->GetRingPairsBySinoIndex(sino_bin,rings1,rings2);
          // Compute indices without gaps (no need to check wether we are in a gap because gaps are already removed from the span table)
          for (int r=0; r<nb_ring_pairs; r++)
          {
            rings1[r] -= rings1[r]/(nb_axial_crystals+1);
            rings2[r] -= rings2[r]/(nb_axial_crystals+1);
          }
          // Crystal indices
          mp_TableBiograph->GetCrystalIDsFromElemView(elem_bin, view_bin, crystals1, crystals2);
          // Check if we are in a gap, we skip it
          if ((crystals1[0]+1)%(nb_trans_crystals_in_head+1)==0)
          {
            LogCerr ("***** oSimulator::ApplyCounts() -> Critical internal error while generating list-mode: events into a gap !! Stop." << endl);
            return 1;
          }
          if ((crystals2[0]+1)%(nb_trans_crystals_in_head+1)==0)
          {
            LogCerr ("***** oSimulator::ApplyCounts() -> Critical internal error while generating list-mode: events into a gap !! Stop." << endl);
            return 1;
          }
          // Compute the crystal indices on ring without gaps
          crystals1[0] -= crystals1[0]/(nb_trans_crystals_in_head+1);
          crystals2[0] -= crystals2[0]/(nb_trans_crystals_in_head+1);
        }
        else if (scanner_model==SCANNER_BIOGRAPH2D)
        {
          nb_ring_pairs = 1;
          rings1[0] = 0;
          rings2[0] = 0;
          // Crystal indices
          mp_TableBiograph2D->GetCrystalIDsFromElemView(elem_bin, view_bin, crystals1, crystals2);
          // Check if we are in a gap, we skip it
          if ((crystals1[0]+1)%(nb_trans_crystals_in_head+1)==0)
          {
            LogCerr ("***** oSimulator::ApplyCounts() -> Critical internal error while generating list-mode: events into a gap !! Stop." << endl);
            return 1;
          }
          if ((crystals2[0]+1)%(nb_trans_crystals_in_head+1)==0)
          {
            LogCerr ("***** oSimulator::ApplyCounts() -> Critical internal error while generating list-mode: events into a gap !! Stop." << endl);
            return 1;
          }
          // Compute the crystal indices on ring without gaps
          crystals1[0] -= crystals1[0]/(nb_trans_crystals_in_head+1);
          crystals2[0] -= crystals2[0]/(nb_trans_crystals_in_head+1);
        }
        else if (scanner_model==SCANNER_MMR2D)
        {
          nb_ring_pairs = 1;
          rings1[0] = 0;
          rings2[0] = 0;
          // Crystal indices
          mp_TableMMR2D->GetCrystalIDsFromElemView(elem_bin, view_bin, crystals1, crystals2);
          // Check if we are in a gap, we skip it
          if ((crystals1[0]+1)%(nb_trans_crystals_in_head+1)==0)
          {
            LogCerr ("***** oSimulator::ApplyCounts() -> Critical internal error while generating list-mode: events into a gap !! Stop." << endl);
            return 1;
          }
          if ((crystals2[0]+1)%(nb_trans_crystals_in_head+1)==0)
          {
            LogCerr ("***** oSimulator::ApplyCounts() -> Critical internal error while generating list-mode: events into a gap !! Stop." << endl);
            return 1;
          }
          // Compute the crystal indices on ring without gaps
          crystals1[0] -= crystals1[0]/(nb_trans_crystals_in_head+1);
          crystals2[0] -= crystals2[0]/(nb_trans_crystals_in_head+1);
        }
        // Event type
        unsigned short int event_type_prompt = 1;
        nb_data_written += fwrite(&event_type_prompt,sizeof(unsigned short int),1,flm);
        // Event time
        unsigned int event_time = 1; // 1 ms
        nb_data_written += fwrite(&event_time,sizeof(unsigned int),1,flm);
        // Scatter rate
        float scatter_rate = 0.;
        if (m_ScatCorr) scatter_rate = mp_SinoScat[sino_bin][view_bin][elem_bin];
        nb_data_written += fwrite(&scatter_rate,sizeof(float),1,flm);
        // Random rate
        float random_rate = 0.;
        if (m_RandCorr) random_rate = mp_SinoRand[sino_bin][view_bin][elem_bin];
        nb_data_written += fwrite(&random_rate,sizeof(float),1,flm);
        // Normalization factor
        float norm_factor = mp_SinoNorm[sino_bin][view_bin][elem_bin];
        nb_data_written += fwrite(&norm_factor,sizeof(float),1,flm);
        // Number of contributing LORs
        unsigned short int nb_lors = nb_ring_pairs * m_Mash;
        nb_data_written += fwrite(&nb_lors,sizeof(unsigned short int),1,flm);
        // Loop on all ring pairs
        for (int ring=0; ring<nb_ring_pairs; ring++)
        {
          // Loop on all mashing LORs
          for (int mash=0; mash<m_Mash; mash++)
          {
            // Crystal IDs
            unsigned int id1 = rings1[ring]*nb_total_trans_crystals + crystals1[mash];
            unsigned int id2 = rings2[ring]*nb_total_trans_crystals + crystals2[mash];
            nb_data_written += fwrite(&id1,sizeof(unsigned int),1,flm);
            nb_data_written += fwrite(&id2,sizeof(unsigned int),1,flm);
            // Attenuation correction factor
            float acf = 1.;
            if (m_AttnCorr) acf = mp_SinoAttn[sino_bin][view_bin][elem_bin];
            nb_data_written += fwrite(&acf,sizeof(float),1,flm);
          }
        }
        // Fill equal number of LORs
        if (f_FillEqualLORs)
        {
          unsigned int id_factice = 0;
          float acf_factice = -1.;
          for (int l=0; l<max_nb_lors-nb_lors; l++)
          {
            nb_data_written += fwrite(&id_factice,sizeof(unsigned int),1,flm);
            nb_data_written += fwrite(&id_factice,sizeof(unsigned int),1,flm);
            nb_data_written += fwrite(&acf_factice,sizeof(float),1,flm);
          }
        }
        // Check writting
        int nb_data_to_be_written = 0;
        if (f_FillEqualLORs) nb_data_to_be_written = 6 + max_nb_lors*3;
        else nb_data_to_be_written = 6 + nb_lors*3;
        if (nb_data_written!=nb_data_to_be_written)
        {
          LogCerr ("***** oSimulator::ApplyCounts() -> Failed to write all data (" << nb_data_to_be_written << " in output esteban file (" << nb_data_written << " written) !" << endl);
          fclose(flm);
          return 1;
        }
      }

      // Close list-mode
      fclose(flm);
      if (m_Verbose>=1)
      {
        cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
        cout << "          100 %                       " << endl;
      }

      // Free buffers
      free(crystals1);
      free(crystals2);
      free(rings1);
      free(rings2);

      // Free thread buffers
      for (int th=0; th<m_NbThreads; th++)
      {
        free(thread_elem_bin[th]);
        free(thread_view_bin[th]);
        free(thread_sino_bin[th]);
      }

      // Save list-mode header file
      string header_file = "";
      if (m_NbReplicates==1) header_file = pathName + baseName+"/"+baseName+"_lm.elm.hdr";
      else header_file = pathName + baseName + "/" + baseName + "_lm_rep" + str_rep + ".elm.hdr";
      if (m_Verbose>=2) LogCout ("      Write header file '" << header_file << "' ..." << endl);
      // Open the file
      ofstream fhdr(header_file.c_str());
      if (!fhdr)
      {
        LogCerr ("***** oSimulator::ApplyCounts() -> Failed to create output header file '" << header_file << "' !" << endl);
        return 1;
      }
      // Write some infos in the header file
      if (m_NbReplicates==1) fhdr << "List-mode := " << baseName << "_lm.elm" << endl;
      else fhdr << "List-mode := " << baseName << "_lm_rep" << str_rep << ".elm" << endl;
      fhdr << "Nb events := " << m_CountPrompt << endl;
      fhdr << "Start time (sec) := 0." << endl;
      fhdr << "Stop time (sec) := 1." << endl;
      if (mp_Scanner->GetScannerModel()==SCANNER_HRRT) fhdr << "PET system := HRRT" << endl;
      else if (mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH) fhdr << "PET system := Biograph" << endl;
      else if (mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH2D) fhdr << "PET system := Biograph2D" << endl;
      else if (mp_Scanner->GetScannerModel()==SCANNER_MMR2D) fhdr << "PET system := MMR2D" << endl;
      else if (mp_Scanner->GetScannerModel()==SCANNER_HRPLUS) fhdr << "PET system := HR+" << endl;
      fhdr << "Input := list-mode" << endl;
      fhdr << "Span := " << m_Span << endl;
      fhdr << "MaxRingDiff := " << m_MaxRingDiff << endl;
      fhdr << "Mashing := " << m_Mash << endl;
      fhdr << "ECF := " << m_ECF << endl;
      fhdr << "Isotope := INF" << endl;
      fhdr << "Dose fraction := 1." << endl;
      fhdr << "Random generator seed := " << m_Seed << endl;
      if (f_FillEqualLORs) fhdr << "Equal number of LORs := YES" << endl;
      else fhdr << "Equal number of LORs := NO" << endl;
      fhdr << "Prompts := " << m_CountPrompt << endl;
      fhdr << "Delays := 0" << endl;
      fhdr << "Trues := " << m_CountTrue << endl;
      fhdr << "Scatters := " << m_CountScat << endl;
      fhdr << "Randoms := " << m_CountRand << endl;
      fhdr << "Trues ratio := " << ((float)m_CountTrue)*100./((float)m_CountPrompt) << endl;
      // Close the file
      fhdr.close();
    }
    // End the replicates loop
    // ----------------------------------------------------------------
  }
  // _____________________________________________________________________________________________________________________________________
  // Prompt with poisson noise
  else if (f_NbCounts>=0 || f_ECF>0.)
  {
    // Compute total prompt and search maximum
    if (m_Verbose>=2) LogCout ("      Compute total prompt" << flush);
    double tmp_total_prompt = 0.;
    for (int s=0; s<m_NbSino; s++) for (int v=0; v<m_NbView; v++) for (int e=0; e<m_NbElem; e++)
    {
      double prompt = (double)mp_SinoForw[s][v][e];
      if (m_ScatCorr) prompt += (double)mp_SinoScat[s][v][e];
      if (m_RandCorr) prompt += (double)mp_SinoRand[s][v][e];
      tmp_total_prompt += prompt;
    }
    LogCout (": " << tmp_total_prompt << endl);

    // Determine ECF and scale each sinogram
    if (f_ECF<=0.)
    {
      LogCout ("      Number of counts provided, so compute associated ECF" << endl);
      m_ECF = tmp_total_prompt / ((double)f_NbCounts);
    }
    else if (f_NbCounts<=0)
    {
      LogCout ("      ECF provided, so compute associated number of counts" << endl);
      m_ECF = f_ECF;
      f_NbCounts = ((long int)( tmp_total_prompt / m_ECF ));
    }
    else
    {
      LogCerr ("***** oSimulator::ApplyCounts() -> Internal problem, reason unknown ... Sorry." << endl);
      return 1;
    }
    if (m_ECF==0.)
    {
      LogCerr ("***** oSimulator::ApplyCounts() -> Total prompt data is null ! Emission phantom is probably empty... Abort." << endl);
      return 1;
    }
    if (m_Verbose>=2) LogCout ("      Scale sinograms (ECF: " << m_ECF << " for total counts of " << f_NbCounts << ")" << endl);
    for (int s=0; s<m_NbSino; s++) for (int v=0; v<m_NbView; v++) for (int e=0; e<m_NbElem; e++)
    {
      mp_SinoForw[s][v][e] = (float)(((double)mp_SinoForw[s][v][e])/m_ECF);
      if (m_ScatCorr) mp_SinoScat[s][v][e] = (float)(((double)mp_SinoScat[s][v][e])/m_ECF);
      if (m_RandCorr) mp_SinoRand[s][v][e] = (float)(((double)mp_SinoRand[s][v][e])/m_ECF);
    }

    // Shoot all random numbers to be thread-safe and repeatable later
    if (m_Verbose>=2) LogCout ("      Shoot random numbers" << endl);
    int nb_shoot_per_bin = 1;
    if (m_ScatCorr) nb_shoot_per_bin++;
    if (m_RandCorr) nb_shoot_per_bin++;
    mp_SinoShoot = (int****)malloc(m_NbSino*sizeof(int***));
    for (int s=0; s<m_NbSino; s++)
    {
      mp_SinoShoot[s] = (int***)malloc(m_NbView*sizeof(int**));
      for (int v=0; v<m_NbView; v++)
      {
        mp_SinoShoot[s][v] = (int**)malloc(m_NbElem*sizeof(int*));
        for (int e=0; e<m_NbElem; e++)
        {
          mp_SinoShoot[s][v][e] = (int*)malloc(nb_shoot_per_bin*sizeof(int));
          for (int z=0; z<nb_shoot_per_bin; z++)
            mp_SinoShoot[s][v][e][z] = rand();
        }
      }
    }

    // Temporary thread-safe counters
    long int* count_true  = (long int*)calloc(m_NbThreads,sizeof(long int));
    long int* count_scat  = (long int*)calloc(m_NbThreads,sizeof(long int));
    long int* count_rand  = (long int*)calloc(m_NbThreads,sizeof(long int));
    long int* count_prompt = (long int*)calloc(m_NbThreads,sizeof(long int));

    // Apply Poisson noise to compute prompt sinogram
    if (m_Verbose>=2) LogCout ("      Apply Poisson noise to compute prompt sinogram" << endl);
    int s;
    #ifdef OMP_MODE
    omp_set_num_threads(m_NbThreads);
    #endif
    #pragma omp parallel for private(s) schedule(static, 1)
    for (s=0; s<m_NbSino; s++)
    {
      // Get the thread number
      int th = 0;
      #ifdef OMP_MODE
      th = omp_get_thread_num();
      #endif

      for (int v=0; v<m_NbView; v++) for (int e=0; e<m_NbElem; e++)
      {
        // The random shoot rank
        int random_rank = 0;
        // Get Poisson sample from trues
        long int prompt_long = PoissonSampleLittle(mp_SinoForw[s][v][e],mp_SinoShoot[s][v][e][random_rank]);
        count_true[th] += prompt_long;
        // Get Poisson sample from scatters
        if (m_ScatCorr)
        {
          random_rank++;
          long int tmp = PoissonSampleLittle(mp_SinoScat[s][v][e],mp_SinoShoot[s][v][e][random_rank]);
          prompt_long += tmp;
          count_scat[th] += tmp;
        }
        // Get Poisson sample from randoms
        if (m_RandCorr)
        {
          random_rank++;
          long int tmp = PoissonSampleLittle(mp_SinoRand[s][v][e],mp_SinoShoot[s][v][e][random_rank]);
          prompt_long += tmp;
          count_rand[th] += tmp;
        }
        // Check for short int overflow
        if (prompt_long>32767)
        {
          LogCerr ("!!!!! Sinogram bin [" << e << ";" << v << ";" << s << "] exceeds short int limit (32767) !" << endl);
          LogCerr ("      Set this bin to 32767." << endl);
          prompt_long = 32767;
        }
        // Affect prompt sinogram
        count_prompt[th] += prompt_long;
        mp_SinoPrompt[s][v][e] = (short int)prompt_long;
      }
    }

    // Add counters
    for (int th=0; th<m_NbThreads; th++)
    {
      m_CountTrue  += count_true[th];
      m_CountScat  += count_scat[th];
      m_CountRand  += count_rand[th];
      m_CountPrompt += count_prompt[th];
    }
  }
/* This was the old version of kind of float mode, in this case we scaled the prompt integer sinogram to have the minimal noise possible...
   Now we can stay in real float mode
  // _____________________________________________________________________________________________________________________________________
  // Prompt scaled to the max value (i.e. almost without noise)
  else
  {
    // Search the maximum value in the forward sinogram
    if (m_Verbose>=2) LogCout ("      Search maximum and calculate ECF scaling" << endl);
    double max_value = 0.;
    for (int s=0; s<m_NbSino; s++) for (int v=0; v<m_NbView; v++) for (int e=0; e<m_NbElem; e++)
    {
      double prompt = (double)mp_SinoForw[s][v][e];
      if (m_ScatCorr) prompt += (double)mp_SinoScat[s][v][e];
      if (m_RandCorr) prompt += (double)mp_SinoRand[s][v][e];
      if (max_value<prompt) max_value = prompt;
    }

    // Deduce the ECF to fit in signed short format
    double max_signed_short_value = 10000.;
    m_ECF = max_value / max_signed_short_value;
    if (m_Verbose>=2) LogCout ("      Max value: " << max_value << " | ECF: " << m_ECF << endl);

    // Scale each sinogram
    if (m_Verbose>=2) LogCout ("      Scale sinograms and compute prompt one" << endl);
    for (int s=0; s<m_NbSino; s++) for (int v=0; v<m_NbView; v++) for (int e=0; e<m_NbElem; e++)
    {
      // The final prompt value
      double prompt = 0.;
      // True component
      prompt += (double)mp_SinoForw[s][v][e];
      mp_SinoForw[s][v][e] = (float)(((double)mp_SinoForw[s][v][e])/m_ECF);
      m_CountTrue += (long int)mp_SinoForw[s][v][e];
      // Scatter component
      if (m_ScatCorr)
      {
        prompt += (double)mp_SinoScat[s][v][e];
        mp_SinoScat[s][v][e] = (float)(((double)mp_SinoScat[s][v][e])/m_ECF);
        m_CountScat += (long int)mp_SinoScat[s][v][e];
      }
      // Random component
      if (m_RandCorr)
      {
        prompt += (double)mp_SinoRand[s][v][e];
        mp_SinoRand[s][v][e] = (float)(((double)mp_SinoRand[s][v][e])/m_ECF);
        m_CountRand += (long int)mp_SinoRand[s][v][e];
      }
      // Affect prompt sinogram
      mp_SinoPrompt[s][v][e] = (short int)(prompt/m_ECF);
      m_CountPrompt += (long int)mp_SinoPrompt[s][v][e];
    }
  }
*/
  // _____________________________________________________________________________________________________________________________________
  // Float mode, we add all contribution to the forward sinogram
  else
  {
    // Set the float mode (this is used when saving sinograms to know that we have to save the forward sino
    m_FloatBool = true;

    // Set the ECF to 1
    m_ECF = 1.;

    // Compute total prompt in the forward sino and compute statistics
    if (m_Verbose>=2) LogCout ("      Compute total float sinogram" << endl);
    double tmp_total_prompt = 0.;
    double tmp_total_true   = 0.;
    double tmp_total_scat   = 0.;
    double tmp_total_rand   = 0.;
    for (int s=0; s<m_NbSino; s++) for (int v=0; v<m_NbView; v++) for (int e=0; e<m_NbElem; e++)
    {
      tmp_total_true += (double)mp_SinoForw[s][v][e];
      if (m_ScatCorr)
      {
        double scatter = (double)mp_SinoScat[s][v][e];
        tmp_total_scat += scatter;
        mp_SinoForw[s][v][e] += ((float)scatter);
      }
      if (m_RandCorr)
      {
        double random = (double)mp_SinoRand[s][v][e];
        tmp_total_rand += random;
        mp_SinoForw[s][v][e] += ((float)random);
      }
      tmp_total_prompt += (double)mp_SinoForw[s][v][e];
    }

    // Verbose
    if (m_Verbose>=1)
    {
      LogCout ("  --> Prompt counts  : " << tmp_total_prompt << endl);
      LogCout ("  --> True counts    : " << tmp_total_true << endl);
      LogCout ("  --> Scatter counts : " << tmp_total_scat << endl);
      LogCout ("  --> Random counts  : " << tmp_total_rand << endl);
    }

    // Return here because we don't have any integer counters here
    return 0;
  }

  // Verbose
  if (m_Verbose>=1)
  {
    LogCout ("  --> Prompt counts  : " << m_CountPrompt << endl);
    LogCout ("  --> True counts    : " << m_CountTrue << endl);
    LogCout ("  --> Scatter counts : " << m_CountScat << endl);
    LogCout ("  --> Random counts  : " << m_CountRand << endl);
  }

  // End
  return 0;
}

// ==========================================================================================================================================
// Function Save sinograms
//   This function saves the output sinograms.
// ==========================================================================================================================================
int oSimulator::SaveSinograms()
{
  // Verbose
  if (m_Verbose>=1) LogCout ("oSimulator::SaveSinograms() -> Save all sinograms" << endl);

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Check for correct initialization
  if (!m_HaveProject)
  {
    LogCerr ("***** oSimulator::SaveSinograms() -> Must call this function after having projected the image !" << endl);
    return 1;
  }

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Get the output root filename
  string baseName = oOutputManager::GetInstance()->GetRootName();
  string pathName = oOutputManager::GetInstance()->GetPathName();

  // ________________________________________________________________________________________________________________________________________
  // ________________________________________________________________________________________________________________________________________
  //                                                    F L O A T   S I N O
  // ________________________________________________________________________________________________________________________________________
  // ________________________________________________________________________________________________________________________________________
  if (m_FloatBool)
  {
    // ----------------------------------------------------------------------------------------------------------------------------------------
    // Write the forward sinogram data

    // Open data file
    string file_float = pathName + baseName+"/"+baseName+"_fw.s";
    FILE* ffloat = fopen(file_float.c_str(),"wb");
    if (ffloat==NULL)
    {
      LogCerr ("***** oSimulator::SaveSinograms() -> Failed to create output data float sinogram '" << file_float << "' !" << endl);
      return 1;
    }

    // Verbose
    if (m_Verbose>=1) LogCout ("  --> Write float sinogram data as '" << file_float << "'" << endl);

    // Write it
    int nb_data_written = 0;
    int nb_data_to_be_written = m_NbSino * m_NbView * m_NbElem;
    for (int s=0; s<m_NbSino; s++) for (int v=0; v<m_NbView; v++) for (int e=0; e<m_NbElem; e++)
    {
      nb_data_written += fwrite(&mp_SinoForw[s][v][e],sizeof(float),1,ffloat);
    }

    // Close and check integrity
    fclose(ffloat);
    if (nb_data_written!=nb_data_to_be_written)
    {
      LogCerr ("***** oSimulator::SaveSinograms() -> Failed to write all data (" << nb_data_to_be_written << ") in output file (" << nb_data_written << " data written) for float sinogram !" << endl);
      return 1;
    }

    // ----------------------------------------------------------------------------------------------------------------------------------------
    // Write the float sinogram header

    // Open header file
    string head_float = pathName + baseName+"/"+baseName+"_fw.s.hdr";
    ofstream hfloat(head_float.c_str());
    if (!hfloat)
    {
      LogCerr ("***** oSimulator::SaveSinograms() -> Failed to create output header float sinogram '" << head_float << "' !" << endl);
      return 1;
    }

    // Verbose
    if (m_Verbose>=1) LogCout ("  --> Write float sinogram header as '" << head_float << "'" << endl);

    // Write file
    hfloat << "!INTERFILE" << endl;
    hfloat << "!name of data file := " << baseName << "_fw.s" << endl;
    if (mp_Scanner->GetScannerModel()==SCANNER_HRPLUS) hfloat << "!originating system := HR+" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_INVEON) hfloat << "!originating system := Inveon" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH) hfloat << "!originating system := Biograph" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH2D) hfloat << "!originating system := Biograph2D" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_MMR2D) hfloat << "!originating system := MMR2D" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_HRRT) hfloat << "!originating system := HRRT" << endl;
    hfloat << "!PET data type := emission" << endl;
    hfloat << "data format := sinogram" << endl;
    hfloat << "number format := float" << endl;
    hfloat << "number of bytes per pixel := 4" << endl;
    hfloat << "axial compression := " << m_Span << endl;
    hfloat << "maximum ring difference := " << m_MaxRingDiff << endl;
    hfloat << "mashing power := " << m_Mash << endl;
    hfloat << "image duration := 1." << endl;
    hfloat << "image start time := 0." << endl;
    hfloat << "Dose type := INF" << endl;
    hfloat << "isotope halflife := 1.e+99" << endl;
    hfloat << "branching factor := 1." << endl;
    hfloat << "Patient ID := " << baseName << endl;
    hfloat << "Patient name := ******************************* " << endl;
    hfloat << "number of dimensions := 3" << endl;
    hfloat << "matrix size [1] := " << m_NbElem << endl;
    hfloat << "matrix size [2] := " << m_NbView << endl;
    hfloat << "matrix size [3] := " << m_NbSino << endl;
    hfloat << "scaling factor (mm/pixel) [1] := 1." << endl;
    hfloat << "scaling factor (mm/pixel) [2] := 1." << endl;
    hfloat << "scaling factor (mm/pixel) [3] := 1." << endl;
    hfloat << "frame := 0" << endl;
    hfloat << "bed := 0" << endl;
    hfloat << "Total Prompts := " << m_TotalTrue+m_TotalScat+m_TotalRand << endl;
    hfloat << "Total Randoms := " << m_TotalRand << endl;
    hfloat << "Total Net Trues := " << m_TotalTrue+m_TotalScat << endl;
    hfloat << "Total Scatters := " << m_TotalScat << endl;
    hfloat << "decay correction factor  := 1." << endl;
    hfloat << "decay correction factor2 := 1." << endl;
    hfloat << "Dead time correction factor := 1." << endl;
    hfloat << "start horizontal bed position (mm) := " << m_OffsetZ << endl;
    hfloat << "Simulated sinogram with random seed := " << m_Seed << endl;

    // Close file
    hfloat.close();
  }
  // ________________________________________________________________________________________________________________________________________
  // ________________________________________________________________________________________________________________________________________
  //                                                       P R O M P T
  // ________________________________________________________________________________________________________________________________________
  // ________________________________________________________________________________________________________________________________________
  else
  {
    // ----------------------------------------------------------------------------------------------------------------------------------------
    // Write the prompt sinogram data

    // Open data file
    string file_prompt = pathName + baseName+"/"+baseName+"_pt.s";
    FILE* fprompt = fopen(file_prompt.c_str(),"wb");
    if (fprompt==NULL)
    {
      LogCerr ("***** oSimulator::SaveSinograms() -> Failed to create output data prompt sinogram '" << file_prompt << "' !" << endl);
      return 1;
    }

    // Verbose
    if (m_Verbose>=1) LogCout ("  --> Write prompt sinogram data as '" << file_prompt << "'" << endl);

    // Write it
    int nb_data_written = 0;
    int nb_data_to_be_written = m_NbSino * m_NbView * m_NbElem;
    for (int s=0; s<m_NbSino; s++) for (int v=0; v<m_NbView; v++) for (int e=0; e<m_NbElem; e++)
    {
      nb_data_written += fwrite(&mp_SinoPrompt[s][v][e],sizeof(short int),1,fprompt);
    }

    // Close and check integrity
    fclose(fprompt);
    if (nb_data_written!=nb_data_to_be_written)
    {
      LogCerr ("***** oSimulator::SaveSinograms() -> Failed to write all data (" << nb_data_to_be_written << ") in output file (" << nb_data_written << " data written) for prompt sinogram !" << endl);
      return 1;
    }

    // ----------------------------------------------------------------------------------------------------------------------------------------
    // Write the prompt sinogram header

    // Open header file
    string head_prompt = pathName + baseName+"/"+baseName+"_pt.s.hdr";
    ofstream hprompt(head_prompt.c_str());
    if (!hprompt)
    {
      LogCerr ("***** oSimulator::SaveSinograms() -> Failed to create output header prompt sinogram '" << head_prompt << "' !" << endl);
      return 1;
    }

    // Verbose
    if (m_Verbose>=1) LogCout ("  --> Write prompt sinogram header as '" << head_prompt << "'" << endl);

    // Write file
    hprompt << "!INTERFILE" << endl;
    hprompt << "!name of data file := " << baseName << "_pt.s" << endl;
    if (mp_Scanner->GetScannerModel()==SCANNER_HRPLUS) hprompt << "!originating system := HR+" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_INVEON) hprompt << "!originating system := Inveon" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH) hprompt << "!originating system := Biograph" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH2D) hprompt << "!originating system := Biograph2D" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_MMR2D) hprompt << "!originating system := MMR2D" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_HRRT) hprompt << "!originating system := HRRT" << endl;
    hprompt << "!PET data type := emission" << endl;
    hprompt << "data format := sinogram" << endl;
    hprompt << "number format := signed integer" << endl;
    hprompt << "number of bytes per pixel := 2" << endl;
    hprompt << "axial compression := " << m_Span << endl;
    hprompt << "maximum ring difference := " << m_MaxRingDiff << endl;
    hprompt << "mashing power := " << m_Mash << endl;
    hprompt << "image duration := 1." << endl;
    hprompt << "image start time := 0." << endl;
    hprompt << "Dose type := INF" << endl;
    hprompt << "isotope halflife := 1.e+99" << endl;
    hprompt << "branching factor := 1." << endl;
    hprompt << "Patient ID := " << baseName << endl;
    hprompt << "Patient name := ******************************* " << endl;
    hprompt << "number of dimensions := 3" << endl;
    hprompt << "matrix size [1] := " << m_NbElem << endl;
    hprompt << "matrix size [2] := " << m_NbView << endl;
    hprompt << "matrix size [3] := " << m_NbSino << endl;
    hprompt << "scaling factor (mm/pixel) [1] := 1." << endl;
    hprompt << "scaling factor (mm/pixel) [2] := 1." << endl;
    hprompt << "scaling factor (mm/pixel) [3] := 1." << endl;
    hprompt << "frame := 0" << endl;
    hprompt << "bed := 0" << endl;
    hprompt << "Total Prompts := " << m_CountPrompt << endl;
    hprompt << "Total Randoms := " << m_CountRand << endl;
    hprompt << "Total Net Trues := " << m_CountPrompt-m_CountRand << endl;
    hprompt << "Total Scatters := " << m_CountScat << endl;
    hprompt << "decay correction factor  := 1." << endl;
    hprompt << "decay correction factor2 := 1." << endl;
    hprompt << "Dead time correction factor := 1." << endl;
    hprompt << "start horizontal bed position (mm) := " << m_OffsetZ << endl;
    hprompt << "Simulated sinogram with random seed := " << m_Seed << endl;

    // Close file
    hprompt.close();
  }

  // ________________________________________________________________________________________________________________________________________
  // ________________________________________________________________________________________________________________________________________
  //                                                       S C A T T E R
  // ________________________________________________________________________________________________________________________________________
  // ________________________________________________________________________________________________________________________________________

  if (m_ScatCorr)
  {
    // ----------------------------------------------------------------------------------------------------------------------------------------
    // Normalize scatter sinogram

    // Verbose
    if (m_Verbose>=1) LogCout ("  --> Normalize scatter sinogram" << endl);

    // Normalize it
    for (int s=0; s<m_NbSino; s++) for (int v=0; v<m_NbView; v++) for (int e=0; e<m_NbElem; e++)
      mp_SinoScat[s][v][e] *= mp_SinoNorm[s][v][e];

    // ----------------------------------------------------------------------------------------------------------------------------------------
    // Write the scatter sinogram data

    // Open data file
    string file_scatter = pathName + baseName+"/"+baseName+"_sc.s";
    FILE* fscatter = fopen(file_scatter.c_str(),"wb");
    if (fscatter==NULL)
    {
      LogCerr ("***** oSimulator::SaveSinograms() -> Failed to create output data scatter sinogram '" << file_scatter << "' !" << endl);
      return 1;
    }

    // Verbose
    if (m_Verbose>=1) LogCout ("  --> Write scatter sinogram data as '" << file_scatter << "'" << endl);

    // Write it
    int nb_data_written = 0;
    int nb_data_to_be_written = m_NbSino * m_NbView * m_NbElem;
    for (int s=0; s<m_NbSino; s++) for (int v=0; v<m_NbView; v++) for (int e=0; e<m_NbElem; e++)
    {
      nb_data_written += fwrite(&mp_SinoScat[s][v][e],sizeof(float),1,fscatter);
    }

    // Close and check integrity
    fclose(fscatter);
    if (nb_data_written!=nb_data_to_be_written)
    {
      LogCerr ("***** oSimulator::SaveSinograms() -> Failed to write all data (" << nb_data_to_be_written << ") in output file (" << nb_data_written << " data written) for scatter sinogram !" << endl);
      return 1;
    }

    // ----------------------------------------------------------------------------------------------------------------------------------------
    // Write the scatter sinogram header

    // Open header file
    string head_scatter = pathName + baseName+"/"+baseName+"_sc.s.hdr";
    ofstream hscatter(head_scatter.c_str());
    if (!hscatter)
    {
      LogCerr ("***** oSimulator::SaveSinograms() -> Failed to create output header scatter sinogram '" << head_scatter << "' !" << endl);
      return 1;
    }

    // Verbose
    if (m_Verbose>=1) LogCout ("  --> Write scatter sinogram header as '" << head_scatter << "'" << endl);

    // Write file
    hscatter << "!INTERFILE" << endl;
    hscatter << "!name of data file := " << baseName << "_sc.s" << endl;
    if (mp_Scanner->GetScannerModel()==SCANNER_HRPLUS) hscatter << "!originating system := HR+" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_INVEON) hscatter << "!originating system := Inveon" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH) hscatter << "!originating system := Biograph" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH2D) hscatter << "!originating system := Biograph2D" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_MMR2D) hscatter << "!originating system := MMR2D" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_HRRT) hscatter << "!originating system := HRRT" << endl;
    hscatter << "!PET data type := scatter" << endl;
    hscatter << "data format := sinogram" << endl;
    hscatter << "number format := float" << endl;
    hscatter << "number of bytes per pixel := 4" << endl;
    hscatter << "axial compression := " << m_Span << endl;
    hscatter << "maximum ring difference := " << m_MaxRingDiff << endl;
    hscatter << "mashing power := " << m_Mash << endl;
    hscatter << "image duration := 1." << endl;
    hscatter << "image start time := 0." << endl;
    hscatter << "branching factor := 1." << endl;
    hscatter << "Patient ID := " << baseName << endl;
    hscatter << "Patient name := ******************************* " << endl;
    hscatter << "number of dimensions := 3" << endl;
    hscatter << "matrix size [1] := " << m_NbElem << endl;
    hscatter << "matrix size [2] := " << m_NbView << endl;
    hscatter << "matrix size [3] := " << m_NbSino << endl;
    hscatter << "scaling factor (mm/pixel) [1] := 1." << endl;
    hscatter << "scaling factor (mm/pixel) [2] := 1." << endl;
    hscatter << "scaling factor (mm/pixel) [3] := 1." << endl;
    hscatter << "frame := 0" << endl;
    hscatter << "bed := 0" << endl;
    hscatter << "start horizontal bed position (mm) := " << m_OffsetZ << endl;

    // Close file
    hscatter.close();
  }

  // ________________________________________________________________________________________________________________________________________
  // ________________________________________________________________________________________________________________________________________
  //                                                       R A N D O M
  // ________________________________________________________________________________________________________________________________________
  // ________________________________________________________________________________________________________________________________________

  if (m_RandCorr)
  {
    // ----------------------------------------------------------------------------------------------------------------------------------------
    // Write the random sinogram data

    // Open data file
    string file_random = pathName + baseName+"/"+baseName+"_rd.s";
    FILE* frandom = fopen(file_random.c_str(),"wb");
    if (frandom==NULL)
    {
      LogCerr ("***** oSimulator::SaveSinograms() -> Failed to create output data random sinogram '" << file_random << "' !" << endl);
      return 1;
    }

    // Verbose
    if (m_Verbose>=1) LogCout ("  --> Write random sinogram data as '" << file_random << "'" << endl);

    // Write it
    int nb_data_written = 0;
    int nb_data_to_be_written = m_NbSino * m_NbView * m_NbElem;
    for (int s=0; s<m_NbSino; s++) for (int v=0; v<m_NbView; v++) for (int e=0; e<m_NbElem; e++)
    {
      nb_data_written += fwrite(&mp_SinoRand[s][v][e],sizeof(float),1,frandom);
    }

    // Close and check integrity
    fclose(frandom);
    if (nb_data_written!=nb_data_to_be_written)
    {
      LogCerr ("***** oSimulator::SaveSinograms() -> Failed to write all data (" << nb_data_to_be_written << ") in output file (" << nb_data_written << " data written) for random sinogram !" << endl);
      return 1;
    }

    // ----------------------------------------------------------------------------------------------------------------------------------------
    // Write the random sinogram header

    // Open header file
    string head_random = pathName + baseName+"/"+baseName+"_rd.s.hdr";
    ofstream hrandom(head_random.c_str());
    if (!hrandom)
    {
      LogCerr ("***** oSimulator::SaveSinograms() -> Failed to create output header random sinogram '" << head_random << "' !" << endl);
      return 1;
    }

    // Verbose
    if (m_Verbose>=1) LogCout ("  --> Write random sinogram header as '" << head_random << "'" << endl);

    // Write file
    hrandom << "!INTERFILE" << endl;
    hrandom << "!name of data file := " << baseName << "_rd.s" << endl;
    if (mp_Scanner->GetScannerModel()==SCANNER_HRPLUS) hrandom << "!originating system := HR+" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_INVEON) hrandom << "!originating system := Inveon" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH) hrandom << "!originating system := Biograph" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH2D) hrandom << "!originating system := Biograph2D" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_MMR2D) hrandom << "!originating system := MMR2D" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_HRRT) hrandom << "!originating system := HRRT" << endl;
    hrandom << "!PET data type := random" << endl;
    hrandom << "data format := sinogram" << endl;
    hrandom << "number format := float" << endl;
    hrandom << "number of bytes per pixel := 4" << endl;
    hrandom << "axial compression := " << m_Span << endl;
    hrandom << "maximum ring difference := " << m_MaxRingDiff << endl;
    hrandom << "mashing power := " << m_Mash << endl;
    hrandom << "image duration := 1." << endl;
    hrandom << "image start time := 0." << endl;
    hrandom << "branching factor := 1." << endl;
    hrandom << "Patient ID := " << baseName << endl;
    hrandom << "Patient name := ******************************* " << endl;
    hrandom << "number of dimensions := 3" << endl;
    hrandom << "matrix size [1] := " << m_NbElem << endl;
    hrandom << "matrix size [2] := " << m_NbView << endl;
    hrandom << "matrix size [3] := " << m_NbSino << endl;
    hrandom << "scaling factor (mm/pixel) [1] := 1." << endl;
    hrandom << "scaling factor (mm/pixel) [2] := 1." << endl;
    hrandom << "scaling factor (mm/pixel) [3] := 1." << endl;
    hrandom << "frame := 0" << endl;
    hrandom << "bed := 0" << endl;
    hrandom << "start horizontal bed position (mm) := " << m_OffsetZ << endl;

    // Close file
    hrandom.close();
  }

  // ________________________________________________________________________________________________________________________________________
  // ________________________________________________________________________________________________________________________________________
  //                                                       A T T E N U A T I O N
  // ________________________________________________________________________________________________________________________________________
  // ________________________________________________________________________________________________________________________________________

  if (m_AttnCorr)
  {
    // ----------------------------------------------------------------------------------------------------------------------------------------
    // Write the attenuation sinogram data

    // Open data file
    string file_attn = pathName + baseName+"/"+baseName+"_at.s";
    FILE* fattn = fopen(file_attn.c_str(),"wb");
    if (fattn==NULL)
    {
      LogCerr ("***** oSimulator::SaveSinograms() -> Failed to create output data attn sinogram '" << file_attn << "' !" << endl);
      return 1;
    }

    // Verbose
    if (m_Verbose>=1) LogCout ("  --> Write attn sinogram data as '" << file_attn << "'" << endl);

    // Write it
    int nb_data_written = 0;
    int nb_data_to_be_written = m_NbSino * m_NbView * m_NbElem;
    for (int s=0; s<m_NbSino; s++) for (int v=0; v<m_NbView; v++) for (int e=0; e<m_NbElem; e++)
    {
      nb_data_written += fwrite(&mp_SinoAttn[s][v][e],sizeof(float),1,fattn);
    }

    // Close and check integrity
    fclose(fattn);
    if (nb_data_written!=nb_data_to_be_written)
    {
      LogCerr ("***** oSimulator::SaveSinograms() -> Failed to write all data (" << nb_data_to_be_written << ") in output file (" << nb_data_written << " data written) for attn sinogram !" << endl);
      return 1;
    }

    // ----------------------------------------------------------------------------------------------------------------------------------------
    // Write the attn sinogram header

    // Open header file
    string head_attn = pathName + baseName+"/"+baseName+"_at.s.hdr";
    ofstream hattn(head_attn.c_str());
    if (!hattn)
    {
      LogCerr ("***** oSimulator::SaveSinograms() -> Failed to create output header attn sinogram '" << head_attn << "' !" << endl);
      return 1;
    }

    // Verbose
    if (m_Verbose>=1) LogCout ("  --> Write attn sinogram header as '" << head_attn << "'" << endl);

    // Write file
    hattn << "!INTERFILE" << endl;
    hattn << "!name of data file := " << baseName << "_at.s" << endl;
    if (mp_Scanner->GetScannerModel()==SCANNER_HRPLUS) hattn << "!originating system := HR+" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_INVEON) hattn << "!originating system := Inveon" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH) hattn << "!originating system := Biograph" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH2D) hattn << "!originating system := Biograph2D" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_MMR2D) hattn << "!originating system := MMR2D" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_HRRT) hattn << "!originating system := HRRT" << endl;
    hattn << "!PET data type := attn" << endl;
    hattn << "data format := sinogram" << endl;
    hattn << "number format := float" << endl;
    hattn << "number of bytes per pixel := 4" << endl;
    hattn << "axial compression := " << m_Span << endl;
    hattn << "maximum ring difference := " << m_MaxRingDiff << endl;
    hattn << "mashing power := " << m_Mash << endl;
    hattn << "image duration := 1." << endl;
    hattn << "image start time := 0." << endl;
    hattn << "branching factor := 1." << endl;
    hattn << "Patient ID := " << baseName << endl;
    hattn << "Patient name := ******************************* " << endl;
    hattn << "number of dimensions := 3" << endl;
    hattn << "matrix size [1] := " << m_NbElem << endl;
    hattn << "matrix size [2] := " << m_NbView << endl;
    hattn << "matrix size [3] := " << m_NbSino << endl;
    hattn << "scaling factor (mm/pixel) [1] := 1." << endl;
    hattn << "scaling factor (mm/pixel) [2] := 1." << endl;
    hattn << "scaling factor (mm/pixel) [3] := 1." << endl;
    hattn << "frame := 0" << endl;
    hattn << "bed := 0" << endl;
    hattn << "start horizontal bed position (mm) := " << m_OffsetZ << endl;

    // Close file
    hattn.close();
  }

  // ________________________________________________________________________________________________________________________________________
  // ________________________________________________________________________________________________________________________________________
  //                                                   N O R M A L I Z A T I O N
  // ________________________________________________________________________________________________________________________________________
  // ________________________________________________________________________________________________________________________________________
  {
    // ----------------------------------------------------------------------------------------------------------------------------------------
    // Write the norm sinogram data

    // Open data file
    string file_norm = pathName + baseName+"/"+baseName+"_nm.s";
    FILE* fnorm = fopen(file_norm.c_str(),"wb");
    if (fnorm==NULL)
    {
      LogCerr ("***** oSimulator::SaveSinograms() -> Failed to create output data normalization sinogram '" << file_norm << "' !" << endl);
      return 1;
    }

    // Verbose
    if (m_Verbose>=1) LogCout ("  --> Write normalization sinogram data as '" << file_norm << "'" << endl);

    // Write it
    int nb_data_written = 0;
    int nb_data_to_be_written = m_NbSino * m_NbView * m_NbElem;
    for (int s=0; s<m_NbSino; s++) for (int v=0; v<m_NbView; v++) for (int e=0; e<m_NbElem; e++)
    {
      float value = (float)mp_SinoNorm[s][v][e];
      nb_data_written += fwrite(&value,sizeof(float),1,fnorm);
    }

    // Close and check integrity
    fclose(fnorm);
    if (nb_data_written!=nb_data_to_be_written)
    {
      LogCerr ("***** oSimulator::SaveSinograms() -> Failed to write all data (" << nb_data_to_be_written << ") in output file (" << nb_data_written << " data written) !" << endl);
      return 1;
    }

    // ----------------------------------------------------------------------------------------------------------------------------------------
    // Write the norm sinogram header

    // Open header file
    string head_norm = pathName + baseName+"/"+baseName+"_nm.s.hdr";
    ofstream hnorm(head_norm.c_str());
    if (!hnorm)
    {
      LogCerr ("***** oSimulator::SaveSinograms() -> Failed to create output header normalization sinogram '" << head_norm << "' !" << endl);
      return 1;
    }

    // Verbose
    if (m_Verbose>=1) LogCout ("  --> Write normalization sinogram header as '" << head_norm << "'" << endl);

    // Write it
    hnorm << "!INTERFILE" << endl;
    hnorm << "!name of data file := " << baseName << "_nm.s" << endl;
    if (mp_Scanner->GetScannerModel()==SCANNER_HRPLUS) hnorm << "!originating system := HR+" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_INVEON) hnorm << "!originating system := Inveon" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH) hnorm << "!originating system := Biograph" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH2D) hnorm << "!originating system := Biograph2D" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_MMR2D) hnorm << "!originating system := MMR2D" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_HRRT) hnorm << "!originating system := HRRT" << endl;
    hnorm << "!PET data type := normalization" << endl;
    hnorm << "Patient ID := " << baseName << endl;
    hnorm << "data format := sinogram" << endl;
    hnorm << "number format := float" << endl;
    hnorm << "number of bytes per pixel := 4" << endl;
    hnorm << "axial compression := " << m_Span << endl;
    hnorm << "maximum ring difference := " << m_MaxRingDiff << endl;
    hnorm << "mashing power := " << m_Mash << endl;
    hnorm << "number of dimensions := 3" << endl;
    hnorm << "matrix size [1] := " << m_NbElem << endl;
    hnorm << "matrix size [2] := " << m_NbView << endl;
    hnorm << "matrix size [3] := " << m_NbSino << endl;

    // Close file
    hnorm.close();

    // ----------------------------------------------------------------------------------------------------------------------------------------
    // Write the norm sinogram calib file

    // Open calib file
    string calib_norm = pathName + baseName+"/"+baseName+"_nm.calib.hdr";
    ofstream calib(calib_norm.c_str());
    if (!calib)
    {
      LogCerr ("***** oSimulator::SaveSinograms() -> Failed to create output calibration normalization file '" << calib_norm << "' !" << endl);
      return 1;
    }

    // Verbose
    if (m_Verbose>=1) LogCout ("  --> Write normalization calibration header as '" << calib_norm << "'" << endl);

    // Write file
    calib << "!INTERFILE" << endl;
    if (mp_Scanner->GetScannerModel()==SCANNER_HRPLUS) calib << "!originating system := HR+" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_INVEON) calib << "!originating system := Inveon" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH) calib << "!originating system := Biograph" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH2D) calib << "!originating system := Biograph2D" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_MMR2D) calib << "!originating system := MMR2D" << endl;
    else if (mp_Scanner->GetScannerModel()==SCANNER_HRRT) calib << "!originating system := HRRT" << endl;
    calib << "calibration factor := " << m_ECF << endl;
    calib << "calibration unit := Bq/cc" << endl;
    calib << "branching factor corrected := true" << endl;
    for (int p=0; p<mp_Scanner->GetNbPlanes(); p++) calib << "efficient factor for plane " << p << " := 1" << endl;

    // Close file
    calib.close();
  }

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // End
  return 0;
}

// __________________________________________________________________________________________________________________________________________
// __________________________________________________________________________________________________________________________________________
// __________________________________________________________________________________________________________________________________________
//
//                                                  PRIVATE SECTION OF ROUTINES IMPLEMENTATIONS
// __________________________________________________________________________________________________________________________________________
// __________________________________________________________________________________________________________________________________________
// __________________________________________________________________________________________________________________________________________

// ==========================================================================================================================================
// Function Make3DGaussianKernel
//   This function calculate the kernel components for the convolution of a 3D gaussian.
// ==========================================================================================================================================
void oSimulator::Make3DGaussianKernel()
{
  // Verbose
  if (m_Verbose>=1) LogCout ("oSimulator::Make3DGaussianKernel() -> Calculating kernel components" << endl);

  // Variables
  float sigma_x = m_PsfTransFWHM/(2.*sqrt(2.*log(2.)));
  float sigma_y = m_PsfTransFWHM/(2.*sqrt(2.*log(2.)));
  float sigma_z = m_PsfAxialFWHM/(2.*sqrt(2.*log(2.)));
  float mu_x = (float)(m_PsfKernSizeX / 2);
  float mu_y = (float)(m_PsfKernSizeY / 2);
  float mu_z = (float)(m_PsfKernSizeZ / 2);

  // Compute kernel
  float sum_kernel = 0.;
  for (int l = 0; l < m_PsfKernSizeZ; l++)
  {
    float kern_orig1 = 0.5 * ( erf(  ( (((float)l)-mu_z+0.5)*m_VoxSizeZ) / (sqrt(2.)*sigma_z)  )
                              - erf(  ( (((float)l)-mu_z-0.5)*m_VoxSizeZ) / (sqrt(2.)*sigma_z)  ) );
    for (int k = 0; k < m_PsfKernSizeY; k++)
    {
      float kern_orig2 = 0.5 * ( erf(  ( (((float)k)-mu_y+0.5)*m_VoxSizeY) / (sqrt(2.)*sigma_y)  )
                                - erf(  ( (((float)k)-mu_y-0.5)*m_VoxSizeY) / (sqrt(2.)*sigma_y)  ) )
                              * kern_orig1;
      for (int j = 0; j < m_PsfKernSizeX; j++)
      {
        mp_PsfKernel[j][k][l] = 0.5 * ( erf(  ( (((float)j)-mu_x+0.5)*m_VoxSizeX) / (sqrt(2.)*sigma_x)  )
                                      - erf(  ( (((float)j)-mu_x-0.5)*m_VoxSizeX) / (sqrt(2.)*sigma_x)  ) )
                                    * kern_orig2;
        sum_kernel += mp_PsfKernel[j][k][l];
      }
    }
  }

  // Normalize kernel
  for (int j=0; j<m_PsfKernSizeX; j++) for (int k=0; k<m_PsfKernSizeY; k++) for (int l=0; l<m_PsfKernSizeZ; l++) mp_PsfKernel[j][k][l] /= sum_kernel;
}

// ==========================================================================================================================================
// Function Convolve3D
//   This function convolves the input image with the PSF and put it into the image to be projected.
// ==========================================================================================================================================
void oSimulator::Convolve3D(float* fp_Input, float* fp_Result)
{
  // Verbose
  if (m_Verbose>=1) LogCout ("oSimulator::Convolve3D() -> Convolving input image with the PSF ..." << endl);

  int xy_dim2 = m_DimX*m_DimY;
  int npointsX_on_2 = m_PsfKernSizeX/2;
  int npointsY_on_2 = m_PsfKernSizeY/2;
  int npointsZ_on_2 = m_PsfKernSizeZ/2;
  // Reset output image
  for (int pp=0; pp<m_DimTot; pp++) fp_Result[pp] = 0.;
  // Set number of threads
  #ifdef OMP_MODE
  omp_set_num_threads(m_NbThreads);
  #endif
  // Convolve
  int x;
  #ifdef OMP_MODE
  #pragma omp parallel for private(x)
  #endif
  for (x = 0; x < m_DimX; x++)
  {
    int z, y, zz, xx , yy, zzz, xxx, yyy;
    for (z = 0; z < m_DimZ; z++) for (y = 0; y < m_DimY; y++)
    {
      int index1d_out = z*xy_dim2 + y*m_DimX + x;
      float kernel_sum = 0.; // Will be used for integral density conservation
      for (zz=-npointsZ_on_2; zz<=npointsZ_on_2; zz++) for (yy=-npointsY_on_2; yy<=npointsY_on_2; yy++) for (xx=-npointsX_on_2; xx<=npointsX_on_2; xx++)
      {
        int z_in = z+zz;
        int y_in = y+yy;
        int x_in = x+xx;
        if ( x_in>=0 && x_in<m_DimX && y_in>=0 && y_in<m_DimY && z_in>=0 && z_in<m_DimZ )
        {
          int index1d_in = z_in*xy_dim2 + y_in*m_DimX + x_in;
          int x_kern = xx+npointsX_on_2;
          int y_kern = yy+npointsY_on_2;
          int z_kern = zz+npointsZ_on_2;
          fp_Result[index1d_out] += ((float)fp_Input[index1d_in]) * mp_PsfKernel[x_kern][y_kern][z_kern];
          kernel_sum += mp_PsfKernel[x_kern][y_kern][z_kern];
        }
      }
      fp_Result[index1d_out] /= kernel_sum;
    }
  }
}

// ==========================================================================================================================================
// Function SiddonForwardProjection
//   This function forward projects an image along a line given by its end points, using SIDDON algorithm.
// ==========================================================================================================================================
void oSimulator::SiddonForwardProjection( float x1, float y1, float z1, float x2, float y2, float z2,
                                          float* img, float* mumap,
                                          int dimX, int dimY, int dimZ, float voxX, float voxY, float voxZ,
                                          float* emission, float* transmission )
{
  int k1, ck;

  float xp, dl;
  float r1xyz[3], dxyz[3], pn_xyz[3];
  float alpha_1[3]={0.0, 0.0, 0.0}, alpha_n[3]={1.0, 1.0, 1.0}, alpha_min, alpha_max, oldvalue, almin[3], almax[3];
  float d_inc[3], d_xyz[3],i_min[3], i_max[3];
  int step_ijk[3];
  long i_c[3];
  int mxyz [] = {dimX, dimY, dimZ};
  int dimXY = dimX*dimY;

  // Reset values
  *emission = 0.;
  *transmission = 0.;

  // Shift in z half the detector size and add half crystal size
  z1 = z1 + mp_Scanner->GetAxialCrystalSize()/2. - mp_Scanner->GetAxialScannerSize()/2.;
  z2 = z2 + mp_Scanner->GetAxialCrystalSize()/2. - mp_Scanner->GetAxialScannerSize()/2.;

  // Voxel sizes
  float vox[3];
  vox[0] = voxX;
  vox[1] = voxY;
  vox[2] = voxZ;
  dxyz[0] = (x1-x2);
  dxyz[1] = (y1-y2);
  dxyz[2] = (z1-z2);
  r1xyz[0] = x2;
  r1xyz[1] = y2;
  r1xyz[2] = z2;
  dl = sqrt((dxyz[0]*dxyz[0])+(dxyz[1]*dxyz[1])+(dxyz[2]*dxyz[2]));
  pn_xyz[0] = ((float) dimX)/2.*voxX;
  pn_xyz[1] = ((float) dimY)/2.*voxY;
  pn_xyz[2] = ((float) dimZ)/2.*voxZ;


  for (k1=0;k1<3;k1++)
  {
    if (fabs(dxyz[k1]) > 1.0e-5)
    {
      alpha_1[k1] = ( -pn_xyz[k1] - r1xyz[k1]) / dxyz[k1];
      alpha_n[k1] = (  pn_xyz[k1] - r1xyz[k1]) / dxyz[k1];
    }
  }

  alpha_min=0.;
  alpha_max=1.;
  for (k1=0;k1<3;k1++)
  {
    if (alpha_1[k1] < alpha_n[k1])
    {
      almin[k1]=alpha_1[k1];
      almax[k1]=alpha_n[k1];
    }
    else
    {
      almin[k1]=alpha_n[k1];
      almax[k1]=alpha_1[k1];
    }
  }
  for(k1=0;k1<3;k1++)
  {
    if (almin[k1] > alpha_min) alpha_min = almin[k1];
    if (almax[k1] < alpha_max) alpha_max = almax[k1];
  }
  if (alpha_min >= alpha_max)
  {
    // No voxels crossed
    *emission = 0.;
    *transmission = 1.;
    // Corrected recently: must do a return here (otherwise, strange attenuation values where observed in the sinogram 1,105171
    return;
  }

  for(k1 = 0; k1<3; k1++)
  {
    if (fabs(dxyz[k1]) > 1.0e-5)
    {
      if (dxyz[k1] > 0.)
      {
        i_min[k1] = (mxyz[k1]-((int) ((pn_xyz[k1] - r1xyz[k1] - alpha_min*dxyz[k1])/vox[k1]+1.e-4)));
        d_xyz[k1] = ((-pn_xyz[k1] + vox[k1] * i_min[k1] - r1xyz[k1])/dxyz[k1]);
        d_inc[k1] = dl/(dxyz[k1]/vox[k1]);
        if (i_min[k1] == 0)
        {
          i_c[k1]=0;
          d_xyz[k1] = (dl * d_xyz[k1])+d_inc[k1];
        }
        else
        {
//          if (fabs(d_xyz[k1] - alpha_min > 1.0e-6))
          if (fabs(d_xyz[k1] - alpha_min) > 1.0e-6)
          {
            i_c[k1] = i_min[k1] - 1;
            d_xyz[k1] = dl * d_xyz[k1];
          }
          else
          {
            i_c[k1] = i_min[k1];
            d_xyz[k1] = (dl * d_xyz[k1] ) + d_inc[k1];
          }
        }
        step_ijk[k1] = 1;
      }
      else
      {
        i_max[k1] = (int)((pn_xyz[k1] + r1xyz[k1] + alpha_min*dxyz[k1])/vox[k1]+1.0e-4);
        d_xyz[k1] = ((-pn_xyz[k1] + vox[k1]*i_max[k1] - r1xyz[k1])/dxyz[k1]);
        d_inc[k1] = - dl/(dxyz[k1]/vox[k1]);      
        if (i_max[k1] == mxyz[k1])
        {
          i_c[k1] = mxyz[k1] - 1;
          d_xyz[k1] = (dl * d_xyz[k1] + d_inc[k1]);
        }
        else
        {
          if (fabs(d_xyz[k1] - alpha_min)>1.0e-6)
          {
            i_c[k1] = i_max[k1];
            d_xyz[k1] = dl * d_xyz[k1];
          }
          else
          {
            i_c[k1] = i_max[k1] - 1;
            d_xyz[k1] = (dl* d_xyz[k1] + d_inc[k1]);
          }
        }
        step_ijk[k1] = -1;
      }
    }
    else
    {
      step_ijk[k1] = 0;
      d_xyz[k1] = dl*2.;
      i_c[k1] = (int) ((pn_xyz[k1] + r1xyz[k1] + 1.0e-4)/vox[k1]);
    }
  }

  ck = 0;
  for (oldvalue = alpha_min * dl; (i_c[0] < dimX && i_c[1] < dimY && i_c[2] < dimZ && i_c[0] >-1 && i_c[1]>-1 && i_c[2]>-1);)
  {
    ck = (d_xyz[0] < d_xyz[1] ) ? 0 : 1;
    ck = (d_xyz[ck] < d_xyz[2] ) ? ck : 2;
    {
      long int index = i_c[2]*dimXY + i_c[1]*dimX + i_c[0];
      float alpha = d_xyz[ck]-oldvalue;
      // Update emission and transmission values
      *emission += alpha * img[index];
      if (m_AttnCorr) *transmission += alpha * mumap[index];
    }
    oldvalue = d_xyz[ck];
    if (fabs(d_xyz[0] - oldvalue) < 1.0e-6) {d_xyz[0] += d_inc[0];i_c[0] += step_ijk[0];}
    if (fabs(d_xyz[1] - oldvalue) < 1.0e-6) {d_xyz[1] += d_inc[1];i_c[1] += step_ijk[1];}
    if (fabs(d_xyz[2] - oldvalue) < 1.0e-6) {d_xyz[2] += d_inc[2];i_c[2] += step_ijk[2];}
  }

  // Finish the transmission calculation
  if (m_AttnCorr) *transmission = max(exp((*transmission)/10.),1.);
  else *transmission = 1.;
}

// ==========================================================================================================================================
// Function SiddonDidierForwardProjection
//   This function forward projects an image along a line given by its end points, using SIDDON algorithm (Didier's version).
// ==========================================================================================================================================
void oSimulator::SiddonDidierForwardProjection( float p1x, float p1y, float p1z, float p2x, float p2y, float p2z,
                                                float* img, float* mumap,
                                                int dimX, int dimY, int dimZ, float voxX, float voxY, float voxZ,
                                                float* emission, float* transmission )
{
  // Reset values
  *emission = 0.;
  *transmission = 0.;

  // Shift in z half the detector size and add half crystal size
  p1z = p1z + mp_Scanner->GetAxialCrystalSize()/2. - mp_Scanner->GetAxialScannerSize()/2.;
  p2z = p2z + mp_Scanner->GetAxialCrystalSize()/2. - mp_Scanner->GetAxialScannerSize()/2.;

  // Didier's stuff
  long int nPlane_[3];
  nPlane_[0] = dimX+1;
  nPlane_[1] = dimY+1;
  nPlane_[2] = dimZ+1;
  PRECISION fovSizeX = ((PRECISION)dimX)*voxX;
  PRECISION fovSizeY = ((PRECISION)dimY)*voxY;
  PRECISION fovSizeZ = ((PRECISION)dimZ)*voxZ;
  PRECISION xPlane_[2]; xPlane_[1] = fovSizeX / 2.; xPlane_[0] = -xPlane_[1];
  PRECISION yPlane_[2]; yPlane_[1] = fovSizeY / 2.; yPlane_[0] = -yPlane_[1];
  PRECISION zPlane_[2]; zPlane_[1] = fovSizeZ / 2.; zPlane_[0] = -zPlane_[1];

	// Computing the distance between axis
	PRECISION const dX21 = p2x - p1x;
	PRECISION const dY21 = p2y - p1y;
	PRECISION const dZ21 = p2z - p1z;

	// Computing distance between the two points
	PRECISION const d21 = sqrt( dX21 * dX21 + dY21 * dY21 + dZ21 * dZ21 );

	// Computing the parametric value alphaMin and alphaMax
	PRECISION alphaMin = 0.0, alphaMax = 1.0;

	// First step is to compute the extrem alpha for each axis
	// For the X-axis
	if( dX21 != 0.0 )
	{
		PRECISION const alphaX_1  = - ( fovSizeX/2. + p1x ) / dX21;
		PRECISION const alphaX_NX =   ( fovSizeX/2. - p1x ) / dX21;
		alphaMin = max( alphaMin, min( alphaX_1, alphaX_NX ) );
		alphaMax = min( alphaMax, max( alphaX_1, alphaX_NX ) );
	}

	// For the Y-axis
	if( dY21 != 0.0 )
	{
		PRECISION const alphaY_1  = - ( fovSizeY/2. + p1y ) / dY21;
		PRECISION const alphaY_NY =   ( fovSizeY/2. - p1y ) / dY21;
		alphaMin = max( alphaMin, min( alphaY_1, alphaY_NY ) );
		alphaMax = min( alphaMax, max( alphaY_1, alphaY_NY ) );
	}

	// For the Z-axis
	if( dZ21 != 0.0 )
	{
		PRECISION const alphaZ_1  = - ( fovSizeZ/2. + p1z ) / dZ21;
		PRECISION const alphaZ_NZ =   ( fovSizeZ/2. - p1z ) / dZ21;
		alphaMin = max( alphaMin, min( alphaZ_1, alphaZ_NZ ) );
		alphaMax = min( alphaMax, max( alphaZ_1, alphaZ_NZ ) );
	}

	// if alphaMax is less than or equal to alphaMin no intersection
	// and return an empty buffer
	if( alphaMax <= alphaMin )
  {
    // No voxels crossed
    *emission = 0.;
    *transmission = 1.;
    return;
  }

	// Now we have to find the indices of the particular plane
	// (iMin,iMax), (jMin,jMax), (kMin,kMax)
	long int iMin = 0, iMax = 0;
	long int jMin = 0, jMax = 0;
	long int kMin = 0, kMax = 0;

	// For the X-axis
	if( dX21 > 0.0 )
	{
		iMin = ceil( nPlane_[ 0 ] - ( xPlane_[ 1 ] - alphaMin * dX21 - p1x ) / voxX );
		iMax = floor( 1 + ( p1x + alphaMax * dX21 - xPlane_[ 0 ] ) / voxX );
	}
	else if( dX21 < 0.0 )
	{
		iMin = ceil( nPlane_[ 0 ] - ( xPlane_[ 1 ] - alphaMax * dX21 - p1x ) / voxX );
		iMax = floor( 1 + ( p1x + alphaMin * dX21 - xPlane_[ 0 ] ) / voxX );
	}
  else
	{
		iMin = 1, iMax = 0;
	}

	// For the Y-axis
	if( dY21 > 0.0 )
	{
		jMin = ceil( nPlane_[ 1 ] - ( yPlane_[ 1 ] - alphaMin * dY21 - p1y ) / voxY );
		jMax = floor( 1 + ( p1y + alphaMax * dY21 - yPlane_[ 0 ] ) / voxY );
	}
	else if( dY21 < 0.0 )
	{
		jMin = ceil( nPlane_[ 1 ] - ( yPlane_[ 1 ] - alphaMax * dY21 - p1y ) / voxY );
		jMax = floor( 1 + ( p1y + alphaMin * dY21 - yPlane_[ 0 ] ) / voxY );
	}
  else
	{
		jMin = 1, jMax = 0;
	}

	// For the Z-axis
	if( dZ21 > 0.0 )
	{
		kMin = ceil( nPlane_[ 2 ] - ( zPlane_[ 1 ] - alphaMin * dZ21 - p1z ) / voxZ );
		kMax = floor( 1 + ( p1z + alphaMax * dZ21 - zPlane_[ 0 ] ) / voxZ );
	}
	else if( dZ21 < 0.0 )
	{
		kMin = ceil( nPlane_[ 2 ] - ( zPlane_[ 1 ] - alphaMax * dZ21 - p1z ) / voxZ );
		kMax = floor( 1 + ( p1z + alphaMin * dZ21 - zPlane_[ 0 ] ) / voxZ );
	}
  else
	{
		kMin = 1, kMax = 0;
	}

	// Computing the last term n number of intersection
	long int n = ( iMax - iMin + 1 ) + ( jMax - jMin + 1 )
		+ ( kMax - kMin + 1 ) + 1;

	// We create a buffer storing the merging data
	// We merge alphaMin, alphaMax, alphaX, alphaY and alphaZ
	vector< PRECISION > alpha;
	vector< PRECISION > alphaMinMax( 2 );
	alphaMinMax[ 0 ] = alphaMin; alphaMinMax[ 1 ] = alphaMax;

	long int iElement = iMax - iMin + 1;
	vector< PRECISION > alphaX;
	if( iElement > 0 )
	{
		alphaX.resize( iElement );
		vector< PRECISION >::iterator idx = alphaX.begin();
		if( dX21 > 0 )
		{
			for( long int i = iMin; i <= iMax; ++i )
			*idx++ = ( ( xPlane_[ 0 ] + ( i - 1 ) * voxX ) - p1x ) / dX21;
		}
		else if( dX21 < 0 )
		{
			for( long int i = iMax; i >= iMin; --i )
				*idx++ = ( ( xPlane_[ 0 ] + ( i - 1 ) * voxX ) - p1x ) / dX21;
		}
	}

	// For alphaY
	long int jElement = jMax - jMin + 1;
	vector< PRECISION > alphaY;
	if( jElement > 0 )
	{
		alphaY.resize( jElement );
		vector< PRECISION >::iterator idx = alphaY.begin();
		if( dY21 > 0 )
		{
			for( long int j = jMin; j <= jMax; ++j )
				*idx++ = ( ( yPlane_[ 0 ] + ( j - 1 ) * voxY ) - p1y ) / dY21;
		}
		else if( dY21 < 0 )
		{
			for( long int j = jMax; j >= jMin; --j )
				*idx++ = ( ( yPlane_[ 0 ] + ( j - 1 ) * voxY ) - p1y ) / dY21;
		}
	}

	// For alphaZ
	long int kElement = kMax - kMin + 1;
	vector< PRECISION > alphaZ;
	if( kElement > 0 )
	{
		alphaZ.resize( kElement );
		vector< PRECISION >::iterator idx = alphaZ.begin();
		if( dZ21 > 0 )
		{
			for( long int k = kMin; k <= kMax; ++k )
				*idx++ = ( ( zPlane_[ 0 ] + ( k - 1 ) * voxZ ) - p1z ) / dZ21;
		}
		else if( dZ21 < 0 )
		{
			for( long int k = kMax; k >= kMin; --k )
				*idx++ = ( ( zPlane_[ 0 ] + ( k - 1 ) * voxZ ) - p1z ) / dZ21;
		}
	}

	vector< PRECISION > tmpAlpha; // Temporary vector

	// Merging buffer each buffer
	merge(
		alphaMinMax.begin(), alphaMinMax.end(),
		alphaX.begin(), alphaX.end(),
		std::back_inserter( alpha ) );

	tmpAlpha = alpha;
	alpha.clear();

	merge(
		alphaY.begin(), alphaY.end(),
		tmpAlpha.begin(), tmpAlpha.end(),
		std::back_inserter( alpha ) );

	tmpAlpha = alpha;
	alpha.clear();

	merge(
		alphaZ.begin(), alphaZ.end(),
		tmpAlpha.begin(), tmpAlpha.end(),
		std::back_inserter( alpha ) );

		// Computing the index of the voxels
	PRECISION alphaMid = 0.0;
	PRECISION length = 0.0; // voxel intersection length
	long int index = 0; // Index of the voxel
	long int i = 0, j = 0, k = 0; // indices of the voxel

	// Loop over the number of crossed planes
	for( long int nP = 1; nP <= n; ++nP )
	{
		alphaMid = ( alpha[ nP ] + alpha[ nP - 1 ] ) * 0.5;
		i = 1 + ( p1x + alphaMid * dX21 - xPlane_[ 0 ] ) / voxX;
		j = 1 + ( p1y + alphaMid * dY21 - yPlane_[ 0 ] ) / voxY;
		k = 1 + ( p1z + alphaMid * dZ21 - zPlane_[ 0 ] ) / voxZ;

		if( i < 1 || i > ( (int)nPlane_[ 0 ] - 1 ) )
		{
			continue;
		}

		if( j < 1 || j > ( (int)nPlane_[ 1 ] - 1 ) )
		{
			continue;
		}

		if( k < 1 || k > ( (int)nPlane_[ 2 ] - 1 ) )
		{
			continue;
		}

		// Computing the length
		length = d21 * ( alpha[ nP ] - alpha[ nP - 1 ] );
		// Computing the indices
    index = (i-1) + (j-1)*m_DimX + (k-1)*m_DimXY;
/*
		index = ( i - 1 )
			+ ( ( ( nPlane_[ 1 ] - 2 ) - ( j - 1 ) ) * ( nPlane_[ 0 ] - 1 ) )
			+ ( ( ( nPlane_[ 2 ] - 2 ) - ( k - 1 ) ) * ( nPlane_[ 0 ] - 1 )
			* ( nPlane_[ 1 ] - 1 ) );
*/

    // Update emission and transmission values
    *emission += length * img[index];
    if (m_AttnCorr) *transmission += length * mumap[index];
	}

  // Finish the transmission calculation
  if (m_AttnCorr) *transmission = max(exp((*transmission)/10.),1.);
  else *transmission = 1.;
}
// ==========================================================================================================================================
// Function ComputeRandomFanSum
//   This function computes the random amount for each sinogram bin, based on the fan sum.
// ==========================================================================================================================================
void oSimulator::ComputeRandomFanSum()
{
  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Get some usefull stuff
  int scanner_model = mp_Scanner->GetScannerModel();
  int nb_total_trans_crystals = mp_Scanner->GetNbTotalTransCrystals();
  int nb_axial_crystals = mp_Scanner->GetNbAxialCrystals();
  int nb_trans_crystals_in_head = mp_Scanner->GetNbTransCrystalsInHead();

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Set number of threads
  #ifdef OMP_MODE
  omp_set_num_threads(m_NbThreads);
  #endif

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Loop on sinogram bins using openMP
  int sino_index;
  #ifdef OMP_MODE
  #pragma omp parallel for private(sino_index) schedule(static, 1)
  #endif
  for (sino_index=0; sino_index<m_NbSino; sino_index++)
  {
    // Get the thread number
    int t = 0;
    #ifdef OMP_MODE
    t = omp_get_thread_num();
    #endif

    // Verbose
    if (m_Verbose>=2 && t==0)
    {
      float percent = ((float)sino_index)*100./((float)m_NbSino);
      cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
           << "      " << percent << " %              " << flush;
    }

    // Buffers
    int *crystals1 = (int*)malloc(m_Mash*sizeof(int));
    int *crystals2 = (int*)malloc(m_Mash*sizeof(int));

    // Get the number of ring pairs for this sinogram index
    int nb_ring_pairs;
    unsigned int *rings1, *rings2;
    if (scanner_model==SCANNER_HRPLUS)
    {
      nb_ring_pairs = mp_TableHRplus->GetNbRingPairsBySinoIndex(sino_index);
      rings1 = (unsigned int*)malloc(nb_ring_pairs*sizeof(unsigned int));
      rings2 = (unsigned int*)malloc(nb_ring_pairs*sizeof(unsigned int));
      mp_TableHRplus->GetRingPairsBySinoIndex(sino_index,rings1,rings2);
    }
    else if (scanner_model==SCANNER_INVEON)
    {
      nb_ring_pairs = mp_TableInveon->GetNbRingPairsBySinoIndex(sino_index);
      rings1 = (unsigned int*)malloc(nb_ring_pairs*sizeof(unsigned int));
      rings2 = (unsigned int*)malloc(nb_ring_pairs*sizeof(unsigned int));
      mp_TableInveon->GetRingPairsBySinoIndex(sino_index,rings1,rings2);
    }
    else if (scanner_model==SCANNER_HRRT)
    {
      // To be implemented
      cerr << endl << endl << "*!*!*!*!*  Part of code not yet implemented  *!*!*!*!*" << endl << endl;
      exit(-245);
    }
    else if (scanner_model==SCANNER_BIOGRAPH)
    {
      nb_ring_pairs = mp_TableBiograph->GetNbRingPairsBySinoIndex(sino_index);
      rings1 = (unsigned int*)malloc(nb_ring_pairs*sizeof(unsigned int));
      rings2 = (unsigned int*)malloc(nb_ring_pairs*sizeof(unsigned int));
      mp_TableBiograph->GetRingPairsBySinoIndex(sino_index,rings1,rings2);
      // Compute indices without gaps (no need to check wether we are in a gap because gaps are already removed from the span table)
      for (int r=0; r<nb_ring_pairs; r++)
      {
        rings1[r] -= rings1[r]/(nb_axial_crystals+1);
        rings2[r] -= rings2[r]/(nb_axial_crystals+1);
      }
    }
    else if (scanner_model==SCANNER_BIOGRAPH2D || scanner_model==SCANNER_MMR2D)
    {
      nb_ring_pairs = 1;
      rings1 = (unsigned int*)malloc(nb_ring_pairs*sizeof(unsigned int));
      rings2 = (unsigned int*)malloc(nb_ring_pairs*sizeof(unsigned int));
      rings1[0] = 0;
      rings2[0] = 0;
    }

    // Compute the number of LORs
    int nb_lors = nb_ring_pairs * m_Mash;

    // Loops on view and elem
    for (int view_index=0; view_index<m_NbView; view_index++) for (int elem_index=0; elem_index<m_NbElem; elem_index++)
    {
      // The the transaxial crystal indices taking mashing into account
      if (scanner_model==SCANNER_HRPLUS)
      {
        // Crystal indices
        mp_TableHRplus->GetCrystalIDsFromElemView(elem_index, view_index, crystals1, crystals2);
      }
      else if (scanner_model==SCANNER_INVEON)
      {
        // Crystal indices
        mp_TableInveon->GetCrystalIDsFromElemView(elem_index, view_index, crystals1, crystals2);
      }
      else if (scanner_model==SCANNER_HRRT)
      {
        // To be implemented
        cerr << endl << endl << "*!*!*!*!*  Part of code not yet implemented  *!*!*!*!*" << endl << endl;
        exit(-245);
      }
      else if (scanner_model==SCANNER_BIOGRAPH)
      {
        // Crystal indices
        mp_TableBiograph->GetCrystalIDsFromElemView(elem_index, view_index, crystals1, crystals2);
        // Check if we are in a gap, we skip it
        if ((crystals1[0]+1)%(nb_trans_crystals_in_head+1)==0) continue;
        if ((crystals2[0]+1)%(nb_trans_crystals_in_head+1)==0) continue;
        // Compute the crystal indices on ring without gaps
        crystals1[0] -= crystals1[0]/(nb_trans_crystals_in_head+1);
        crystals2[0] -= crystals2[0]/(nb_trans_crystals_in_head+1);
      }
      else if (scanner_model==SCANNER_BIOGRAPH2D)
      {
        // Crystal indices
        mp_TableBiograph2D->GetCrystalIDsFromElemView(elem_index, view_index, crystals1, crystals2);
        // Check if we are in a gap, we skip it
        if ((crystals1[0]+1)%(nb_trans_crystals_in_head+1)==0) continue;
        if ((crystals2[0]+1)%(nb_trans_crystals_in_head+1)==0) continue;
        // Compute the crystal indices on ring without gaps
        crystals1[0] -= crystals1[0]/(nb_trans_crystals_in_head+1);
        crystals2[0] -= crystals2[0]/(nb_trans_crystals_in_head+1);
      }
      else if (scanner_model==SCANNER_MMR2D)
      {
        // Crystal indices
        mp_TableMMR2D->GetCrystalIDsFromElemView(elem_index, view_index, crystals1, crystals2);
        // Check if we are in a gap, we skip it
        if ((crystals1[0]+1)%(nb_trans_crystals_in_head+1)==0) continue;
        if ((crystals2[0]+1)%(nb_trans_crystals_in_head+1)==0) continue;
        // Compute the crystal indices on ring without gaps
        crystals1[0] -= crystals1[0]/(nb_trans_crystals_in_head+1);
        crystals2[0] -= crystals2[0]/(nb_trans_crystals_in_head+1);
      }
      // Loop on all ring pairs
      for (int ring=0; ring<nb_ring_pairs; ring++)
      {
        // Loop on all mashing LORs
        for (int mash=0; mash<m_Mash; mash++)
        {
          // Compute the global crystal IDs
          int id1 = rings1[ring]*nb_total_trans_crystals + crystals1[mash];
          int id2 = rings2[ring]*nb_total_trans_crystals + crystals2[mash];
          // Add contribution to fan sum
          mp_SinoRand[sino_index][view_index][elem_index] += mp_FanSum[id1]*mp_FanSum[id2];
        } // End loop on all mashing LORs
      } // End loop on all ring pairs
    } // End loop on view and elem
    if (rings1) free(rings1);
    if (rings2) free(rings2);
    free(crystals1);
    free(crystals2);
  }
  // Verbose
  if (m_Verbose>=2) cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
                         << "      100 %              " << endl;
}

// ==========================================================================================================================================
// Function ConvolveScatterComponent
//   This function convolves the scatter sinogram.
//   The convolution kernel is arbitrary set, which works quite well for hr+ and biograph.
//   The convolution is 1D along the radial dimension.
// ==========================================================================================================================================
void oSimulator::ConvolveScatterComponent()
{
  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Build arbitrary scatter kernel
  int kernel_size = 0;
  float* kernel = NULL;
  if (mp_Scanner->GetScannerModel()==SCANNER_HRPLUS || mp_Scanner->GetScannerModel()==SCANNER_INVEON)
  {
    kernel_size = 73;
    kernel = (float*)malloc(kernel_size*sizeof(float));
    kernel[0]  = 0.000100;
    kernel[1]  = 0.00012;
    kernel[2]  = 0.00015;
    kernel[3]  = 0.00019;
    kernel[4]  = 0.00024;
    kernel[5]  = 0.00030;
    kernel[6]  = 0.00037;
    kernel[7]  = 0.00045;
    kernel[8]  = 0.00054;
    kernel[9]  = 0.00064;
    kernel[10] = 0.00075;
    kernel[11] = 0.00087;
    kernel[12] = 0.00100;
    kernel[13] = 0.012;
    kernel[14] = 0.015;
    kernel[15] = 0.019;
    kernel[16] = 0.024;
    kernel[17] = 0.030;
    kernel[18] = 0.037;
    kernel[19] = 0.045;
    kernel[20] = 0.054;
    kernel[21] = 0.064;
    kernel[22] = 0.075;
    kernel[23] = 0.087;
    kernel[24] = 0.100;
    kernel[25] = 0.12;
    kernel[26] = 0.15;
    kernel[27] = 0.19;
    kernel[28] = 0.24;
    kernel[29] = 0.30;
    kernel[30] = 0.37;
    kernel[31] = 0.45;
    kernel[32] = 0.54;
    kernel[33] = 0.64;
    kernel[34] = 0.75;
    kernel[35] = 0.87;
    kernel[36] = 1.;
    kernel[37] = 0.87;
    kernel[38] = 0.75;
    kernel[39] = 0.64;
    kernel[40] = 0.54;
    kernel[41] = 0.45;
    kernel[42] = 0.37;
    kernel[43] = 0.30;
    kernel[44] = 0.24;
    kernel[45] = 0.19;
    kernel[46] = 0.15;
    kernel[47] = 0.12;
    kernel[48] = 0.100;
    kernel[49] = 0.087;
    kernel[50] = 0.075;
    kernel[51] = 0.064;
    kernel[52] = 0.054;
    kernel[53] = 0.045;
    kernel[54] = 0.037;
    kernel[55] = 0.030;
    kernel[56] = 0.024;
    kernel[57] = 0.019;
    kernel[58] = 0.015;
    kernel[59] = 0.012;
    kernel[60] = 0.00100;
    kernel[61] = 0.00087;
    kernel[62] = 0.00075;
    kernel[63] = 0.00064;
    kernel[64] = 0.00054;
    kernel[65] = 0.00045;
    kernel[66] = 0.00037;
    kernel[67] = 0.00030;
    kernel[68] = 0.00024;
    kernel[69] = 0.00019;
    kernel[70] = 0.00015;
    kernel[71] = 0.00012;
    kernel[72] = 0.000100;
  }
  else if (mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH || mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH2D || mp_Scanner->GetScannerModel()==SCANNER_MMR2D)
  {
    kernel_size = 169;
    kernel = (float*)malloc(kernel_size*sizeof(float));
    for (int k=0; k<24; k++) kernel[k] = 0.0001;
    kernel[24] = 0.000100;
    kernel[25] = 0.000115;
    kernel[26] = 0.000135;
    kernel[27] = 0.000155;
    kernel[28] = 0.000180;
    kernel[29] = 0.000220;
    kernel[30] = 0.000270;
    kernel[31] = 0.000330;
    kernel[32] = 0.000395;
    kernel[33] = 0.000465;
    kernel[34] = 0.000530;
    kernel[35] = 0.000610;
    kernel[36] = 0.000700;
    kernel[37] = 0.000790;
    kernel[38] = 0.000890;
    kernel[39] = 0.00100;
    kernel[40] = 0.00115;
    kernel[41] = 0.00135;
    kernel[42] = 0.00155;
    kernel[43] = 0.00180;
    kernel[44] = 0.00220;
    kernel[45] = 0.00270;
    kernel[46] = 0.00330;
    kernel[47] = 0.00395;
    kernel[48] = 0.00465;
    kernel[49] = 0.00530;
    kernel[50] = 0.00610;
    kernel[51] = 0.00700;
    kernel[52] = 0.00790;
    kernel[53] = 0.00890;
    kernel[54] = 0.0100;
    kernel[55] = 0.0115;
    kernel[56] = 0.0135;
    kernel[57] = 0.0155;
    kernel[58] = 0.0180;
    kernel[59] = 0.0220;
    kernel[60] = 0.0270;
    kernel[61] = 0.0330;
    kernel[62] = 0.0395;
    kernel[63] = 0.0465;
    kernel[64] = 0.0530;
    kernel[65] = 0.0610;
    kernel[66] = 0.0700;
    kernel[67] = 0.0790;
    kernel[68] = 0.0890;
    kernel[69] = 0.100;
    kernel[70] = 0.115;
    kernel[71] = 0.135;
    kernel[72] = 0.155;
    kernel[73] = 0.180;
    kernel[74] = 0.220;
    kernel[75] = 0.270;
    kernel[76] = 0.330;
    kernel[77] = 0.395;
    kernel[78] = 0.465;
    kernel[79] = 0.530;
    kernel[80] = 0.610;
    kernel[81] = 0.700;
    kernel[82] = 0.790;
    kernel[83] = 0.890;
    kernel[84] = 1.;
    kernel[85] = 0.890;
    kernel[86] = 0.790;
    kernel[87] = 0.700;
    kernel[88] = 0.610;
    kernel[89] = 0.530;
    kernel[90] = 0.465;
    kernel[91] = 0.395;
    kernel[92] = 0.330;
    kernel[93] = 0.270;
    kernel[94] = 0.220;
    kernel[95] = 0.180;
    kernel[96] = 0.155;
    kernel[97] = 0.135;
    kernel[98] = 0.115;
    kernel[100] = 0.100;
    kernel[101] = 0.0890;
    kernel[102] = 0.0790;
    kernel[103] = 0.0700;
    kernel[104] = 0.0610;
    kernel[105] = 0.0530;
    kernel[106] = 0.0465;
    kernel[107] = 0.0395;
    kernel[108] = 0.0330;
    kernel[109] = 0.0270;
    kernel[110] = 0.0220;
    kernel[111] = 0.0180;
    kernel[112] = 0.0155;
    kernel[113] = 0.0135;
    kernel[114] = 0.0115;
    kernel[115] = 0.0100;
    kernel[116] = 0.00890;
    kernel[117] = 0.00790;
    kernel[118] = 0.00700;
    kernel[119] = 0.00610;
    kernel[120] = 0.00530;
    kernel[121] = 0.00465;
    kernel[122] = 0.00395;
    kernel[123] = 0.00330;
    kernel[124] = 0.00270;
    kernel[125] = 0.00220;
    kernel[126] = 0.00180;
    kernel[127] = 0.00155;
    kernel[128] = 0.00135;
    kernel[129] = 0.00115;
    kernel[130] = 0.00100;
    kernel[131] = 0.000890;
    kernel[132] = 0.000790;
    kernel[133] = 0.000700;
    kernel[134] = 0.000610;
    kernel[135] = 0.000530;
    kernel[136] = 0.000465;
    kernel[137] = 0.000395;
    kernel[138] = 0.000330;
    kernel[139] = 0.000270;
    kernel[140] = 0.000220;
    kernel[141] = 0.000180;
    kernel[142] = 0.000155;
    kernel[143] = 0.000135;
    kernel[144] = 0.000115;
    kernel[145] = 0.000100;
    for (int k=146; k<kernel_size; k++) kernel[k] = 0.0001;
  }
  else if (mp_Scanner->GetScannerModel()==SCANNER_HRRT)
  {
    cerr << "!!!***!!!  Scatter convolve not yet implemented for hrrt !" << endl;
    exit(-235);
  }

  float kernel_sum = 0.;
  for (int k=0; k<kernel_size; k++) kernel_sum += kernel[k];
  for (int k=0; k<kernel_size; k++) kernel[k] /= kernel_sum;

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Allocate convolution buffers
  float** buffers = (float**)malloc(m_NbThreads*sizeof(float*));
  for (int th=0; th<m_NbThreads; th++) buffers[th] = (float*)malloc(m_NbElem*sizeof(float));

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Set number of threads
  #ifdef OMP_MODE
  omp_set_num_threads(m_NbThreads);
  #endif

  // ----------------------------------------------------------------------------------------------------------------------------------------
  // Loop on sinogram bins using openMP
  int sino_index;
  #ifdef OMP_MODE
  #pragma omp parallel for private(sino_index) schedule(static, 1)
  #endif
  for (sino_index=0; sino_index<m_NbSino; sino_index++)
  {
    // Get the thread number
    int th = 0;
    #ifdef OMP_MODE
    th = omp_get_thread_num();
    #endif

    // Verbose
    if (m_Verbose>=2 && th==0)
    {
      float percent = ((float)sino_index)*100./((float)m_NbSino);
      cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
           << "      " << percent << " %              " << flush;
    }

    // Loop on view
    for (int view_index=0; view_index<m_NbView; view_index++)
    {
      // Reset buffer
      for (int elem_out=0; elem_out<m_NbElem; elem_out++) buffers[th][elem_out] = 0.;
      // Loop on elem in
      for (int elem_in=0; elem_in<m_NbElem; elem_in++)
      {
        // Loop on convolution kernel
        for (int k=0; k<kernel_size; k++)
        {
          // Compute output index
          int elem_out = elem_in + k - kernel_size/2;
          // Convolve if output index is valid
          if (elem_out>=0 && elem_out<m_NbElem)
          {
            buffers[th][elem_out] += mp_SinoScat[sino_index][view_index][elem_in] * kernel[k];
          }
        }
      }
      // Copy buffer
      for (int elem_out=0; elem_out<m_NbElem; elem_out++) mp_SinoScat[sino_index][view_index][elem_out] = buffers[th][elem_out];
    } // End loop on view

  } // End loop on sino
  // Verbose
  if (m_Verbose>=2) cout << "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
                         << "      100 %              " << endl;
}

// ==========================================================================================================================================
// Function PoissonSampleLittle
//   This function return a sample of the Poisson law of the given mean. The mean cannot be higher than 10000 (does not converge else)
//   It is implemented in long double precision to avoid rounding error stucking the loop.
// ==========================================================================================================================================
long int oSimulator::PoissonSampleLittle(float f_Mean, int f_Random)
{
  // Impose a limit in shooting near 1. (this impose a bais, but limited thought)
  long double shoot_limit = 0.999999;
  // Shooting between 0 and 1 [0;shoot_limit[
  long double shoot = ((long double)f_Random)/((long double)(RAND_MAX));
  if (shoot>shoot_limit) shoot = shoot_limit;
  // This variable contains the poisson sample (incremented along loops)
  long int poisson_value = 0;
  // This variable contains the probability of the current poisson lambda, knowing the mean
  long double proba_poi = exp(-((long double)f_Mean));
  // This variable contains the integrated probability which can be compared to the random uniform value
  long double proba_sum = proba_poi;
  // The loop (stopped when the integrated density is bigger or equal to the shoot
  while (proba_sum <= shoot)
  {
    // Increment Poisson sample
    poisson_value++;
    // Update the poisson probability for this sample
    proba_poi *= ((long double)f_Mean) / ((long double)poisson_value);
    // Update the integrated probability
    proba_sum += proba_poi;
  }
  // It is done
  return poisson_value;
}
long int oSimulator::PoissonSampleLittle(float f_Mean)
{
  // Impose a limit in shooting near 1. (this impose a bais, but limited thought)
  long double shoot_limit = 0.999999;
  // Shooting between 0 and 1 [0;shoot_limit[
  long double shoot = ((long double)rand())/((long double)(RAND_MAX));
  if (shoot>shoot_limit) shoot = shoot_limit;
  // This variable contains the poisson sample (incremented along loops)
  long int poisson_value = 0;
  // This variable contains the probability of the current poisson lambda, knowing the mean
  long double proba_poi = exp(-((long double)f_Mean));
  // This variable contains the integrated probability which can be compared to the random uniform value
  long double proba_sum = proba_poi;
  // The loop (stopped when the integrated density is bigger or equal to the shoot
  while (proba_sum <= shoot)
  {
    // Increment Poisson sample
    poisson_value++;
    // Update the poisson probability for this sample
    proba_poi *= ((long double)f_Mean) / ((long double)poisson_value);
    // Update the integrated probability
    proba_sum += proba_poi;
  }
  // It is done
  return poisson_value;
}
// ==========================================================================================================================================
// Function PoissonSampleBig
//   This function return a sample of the Poisson law of the given mean. The computation is cut in part of 10000.
//   It is implemented in long double precision to avoid rounding error stucking the loop.
// ==========================================================================================================================================
long int oSimulator::PoissonSampleBig(float f_Mean)
{
  // Impose a limit in shooting near 1. (this impose a bais, but limited thought)
  long double shoot_limit = 0.999999;

  // The mean limit to be randomize
  long double mean_limit = 10000.;

  // Final result
  long int final_poisson_value = 0;

  // Rest to be randomize
  long double rest = ((long double)f_Mean);

  // While loop on rest
  while (rest>0.)
  {
    // Mean value to be randomize
    long double mean_value = -1.;
    if (rest>mean_limit) mean_value = mean_limit;
    else mean_value = rest;
    rest -= mean_limit;
    // Shooting between 0 and 1 [0;shoot_limit[
    long double shoot = ((long double)rand())/((long double)(RAND_MAX));
    if (shoot>shoot_limit) shoot = shoot_limit;
    // This variable contains the poisson sample (incremented along loops)
    long int poisson_value = 0;
    // This variable contains the probability of the current poisson lambda, knowing the mean
    long double proba_poi = exp(-mean_value);
    // This variable contains the integrated probability which can be compared to the random uniform value
    long double proba_sum = proba_poi;
    // The loop (stopped when the integrated density is bigger or equal to the shoot
    while (proba_sum <= shoot)
    {
      // Increment Poisson sample
      poisson_value++;
      // Update the poisson probability for this sample
      proba_poi *= mean_value / ((long double)poisson_value);
      // Update the integrated probability
      proba_sum += proba_poi;
    }
    // Add result
    final_poisson_value += poisson_value;
  }
  // It is done
  return final_poisson_value;
}

