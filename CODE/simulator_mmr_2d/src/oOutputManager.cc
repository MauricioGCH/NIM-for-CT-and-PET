#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <ctime>
//#include <chrono>
#include "oOutputManager.hh"
#include "oMiscellaneous.hh"
#include "oScanner.hh"
using namespace std;

// ==========================================================================================================================================
// Static instance for singleton implementation
// ==========================================================================================================================================
oOutputManager* oOutputManager::mp_Instance = 0;

// ==========================================================================================================================================
// Constructor
// ==========================================================================================================================================
oOutputManager::oOutputManager( const string& f_RootName, int f_Verbose )
{
  // Default
  m_Directory = false;

  // Affectations
  m_RootName = f_RootName;
  m_Verbose = f_Verbose;

  // Checks
  if (m_RootName=="")
  {
    cerr << "***** oOutputManager::Constructor() -> Root name is empty !" << endl;
    exit(1);
  }

  // The logging is done only if the name differs from the no-log keyword
  if (m_RootName!=NO_LOG_KEYWORD)
  {
    // Dissociate given path from the real root name
    int pos;
    if ((pos=m_RootName.find_last_of("/"))!=string::npos)
    {
      m_PathName = m_RootName.substr(0,pos) + "/";
      m_RootName = m_RootName.substr(pos+1);
    }
    else
    {
      m_PathName = "";
    }

    // Init output directory
    if (InitOutputDirectory())
    {
      cerr << "***** oOutputManager::Constructor() -> Failed to initialize output directory !" << endl;
      exit(1);
    }

    // Open log file
    if (InitLogFile())
    {
      cerr << "***** oOutputManager::Constructor() -> Failed to initialize log file !" << endl;
      exit(1);
    }
  }

  // Instance
  mp_Instance = this;
}

// ==========================================================================================================================================
// Destructor
// ==========================================================================================================================================
oOutputManager::~oOutputManager()
{
  if (m_Verbose>=3) cout << "oOutputManager::Destructor() -> Destroy" << endl;

  // Close the log file
  if (m_Log) m_Log.close();
}

// ==========================================================================================================================================
// Function InitOutputDirectory
//    This function initializes the directory where output files are stored.
//    It returns the returned value of the system call.
// ==========================================================================================================================================
int oOutputManager::InitOutputDirectory()
{
  // Build the system call
//  #ifdef WINDOWS
  #ifdef _WIN32
  //string instruction = "md " + m_pathName;
  string instruction = "if not exist " + m_RootName + " mkdir " + m_RootName;
  #else
  string instruction = "mkdir -p " + m_PathName + m_RootName;
  #endif

  // Send instruction
  int error = system(instruction.c_str());
  // Check
  if (error)
  {
    cerr << "***** oOutputManager::InitOutputDirectory() -> Failed to correctly create the output directory '" << m_PathName << m_RootName << "' !" << endl;
    return error;
  }
  // Boolean to say that the directory is OK
  m_Directory = true;
  // End
  return 0;
}

// ==========================================================================================================================================
// Function InitLogFile
//    This function initializes the log file by simply opening it.
//    It returns 0 upon success and other value otherwise.
// ==========================================================================================================================================
int oOutputManager::InitLogFile()
{
  // Check
  if (!m_Directory)
  {
    cerr << "***** oOutputManager::InitLogFile() -> Output directory '" << m_PathName << m_RootName << "' is not valid !" << endl;
    return 1;
  }

  // Open file
  string log_file = m_PathName+m_RootName+"/"+m_RootName+".log";
  m_Log.open(log_file.c_str());
  if (!m_Log)
  {
    cerr << "***** oOutputManager::InitLogFile() -> Failed to create log file '" << log_file << "' !" << endl;
    return 1;
  }

  // End
  return 0;
}

// ==========================================================================================================================================
// Function LogCommandLine
//    This function logs the command line.
//    It returns 0 upon success and other value otherwise.
// ==========================================================================================================================================
int oOutputManager::LogCommandLine(int argc, char** argv)
{
  // Check
  if (m_Log)
  {
    m_Log << "==================================================================================================" << endl;
    m_Log << "                                      COMMAND LINE CONTEXT" << endl;
    m_Log << "==================================================================================================" << endl;
    // Print command line
    for (int i=0; i<argc; i++) m_Log << argv[i] << " ";
    m_Log << endl;
    // Get and print date of execution start
    time_t end_time = time(NULL);
// Fuck me !! the ctime function allocates a static char* that shift memory and produces segmentation fault in a very incomprehensible way ! Do not use it !!
//    m_Log << "Date of execution: " << ctime(&end_time) << endl;
    m_Log << "Date of execution (since Epoch): " << end_time << endl;
    m_Log << "Code precision in bytes: " << sizeof(PRECISION) << endl;
    m_Log << "==================================================================================================" << endl;
  }
  else return 1;

  // End
  return 0;
}

// ==========================================================================================================================================
// Function LogCommandLine
//    This function logs the CPU time.
//    It returns 0 upon success and other value otherwise.
// ==========================================================================================================================================
int oOutputManager::LogCPUTime()
{
  // Check
  if (m_Log)
  {
    m_Log << "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=" << endl;
    m_Log << "Current CPU time: " << time(NULL) << endl;
    m_Log << "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=" << endl;
  }
  else return 1;

  // End
  return 0;
}

// ==========================================================================================================================================
// Function SaveWeightMatrix
//    This function saves a weight matrix and associated header file.
//    It returns 0 upon success and other value otherwise.
// ==========================================================================================================================================
int oOutputManager::SaveWeightMatrix( PRECISION* fp_Matrix, int f_DimX, int f_DimY, int f_DimZ, PRECISION f_VoxX, PRECISION f_VoxY, PRECISION f_VoxZ, PRECISION f_OrigX, PRECISION f_OrigY, PRECISION f_OrigZ,
                                      PRECISION f_TimeStart, PRECISION f_Duration, int f_ScannerModel )
{
  // Check
  if (!m_Directory)
  {
    LogCerr ("***** oOutputManager::SaveWeightMatrix() -> Output directory '" << m_PathName << m_RootName << "' is not valid !" << endl);
    return 1;
  }

  // Flip weight matrix (should not be flipped in order to be properly used with the BNP code)
//  Misc_FlipXYZ(fp_Matrix,f_DimX,f_DimY,f_DimZ);

  // Pre-stuffs
  int dim_tot = f_DimX * f_DimY * f_DimZ;

  // -----------------------------------------------------------------------------------------------------------------------------
  //                            Save data
  // -----------------------------------------------------------------------------------------------------------------------------

  // Open file
  string file_weights = m_PathName+m_RootName+"/"+m_RootName+"_weights.i";
  FILE* fout = fopen(file_weights.c_str(),"wb");
  if (fout==NULL)
  {
    LogCerr ("***** oOutputManager::SaveWeightMatrix() -> Failed to create data file '" << file_weights << "' !" << endl);
    return 1;
  }
  // Verbose
  if (m_Verbose>=1) LogCout ("  --> Save weights in file '" << file_weights << "' ..." << endl);
  // Write data
  float* buffer = (float*)malloc(dim_tot*sizeof(float));
  for (int v=0; v<dim_tot; v++) buffer[v] = ((float)fp_Matrix[v]);
  int nb_data_written = fwrite(buffer,sizeof(float),dim_tot,fout);
  // Close file
  fclose(fout);
  free(buffer);
  // Check integrity
  if (nb_data_written!=dim_tot)
  {
    LogCerr ("***** oOutputManager::SaveWeightMatrix() -> Failed to write all data (" << dim_tot << ") in file '" << file_weights << "' (" << nb_data_written << " written) !" << endl);
    return 1;
  }

  // -----------------------------------------------------------------------------------------------------------------------------
  //                            Save header
  // -----------------------------------------------------------------------------------------------------------------------------

  // Open file
  string file_head = m_PathName+m_RootName+"/"+m_RootName+"_weights.i.hdr";
  ofstream fhead(file_head.c_str());
  if (!fhead)
  {
    LogCerr ("***** oOoutputManager::SaveWeightMatrix() -> Failed to create header file '" << file_head << "' !" << endl);
    return 1;
  }
  // Write info
  fhead << "!INTERFILE" << endl;
  fhead << "!name of data file := " << m_RootName << "_weights.i" << endl;
  if (f_ScannerModel==SCANNER_HRPLUS) fhead << "!originating system := HR+" << endl;
  else if (f_ScannerModel==SCANNER_HRRT) fhead << "!originating system := HRRT" << endl;
  else if (f_ScannerModel==SCANNER_BIOGRAPH) fhead << "!originating system := Biograph" << endl;
  else if (f_ScannerModel==SCANNER_BIOGRAPH2D) fhead << "!originating system := Biograph2D" << endl;
  else if (f_ScannerModel==SCANNER_INVEON) fhead << "!originating system := Inveon" << endl;
  else if (f_ScannerModel==SCANNER_FOCUS) fhead << "!originating system := Focus" << endl;
  fhead << "!PET data type := weights matrix used in reconstruction" << endl;
  fhead << "data format := image" << endl;
  fhead << "number format := float" << endl;
  fhead << "number of bytes per pixel := 4" << endl;
  fhead << "Patient name := " << m_RootName << endl;
  fhead << "number of dimensions := 3" << endl;
  fhead << "matrix size [1] := " << f_DimX << endl;
  fhead << "matrix size [2] := " << f_DimY << endl;
  fhead << "matrix size [3] := " << f_DimZ << endl;
  fhead << "scaling factor (mm/pixel) [1] := " << f_VoxX << endl;
  fhead << "scaling factor (mm/pixel) [2] := " << f_VoxY << endl;
  fhead << "scaling factor (mm/pixel) [3] := " << f_VoxZ << endl;
  fhead << "bounding origin (mm) [1] := " << f_OrigX-f_VoxX/2. << endl;
  fhead << "bounding origin (mm) [2] := " << f_OrigY-f_VoxY/2. << endl;
  fhead << "bounding origin (mm) [3] := " << f_OrigZ-f_VoxZ/2. << endl;
  fhead << "FOV size (mm) [1] := " << ((PRECISION)f_DimX)*f_VoxX << endl;
  fhead << "FOV size (mm) [2] := " << ((PRECISION)f_DimY)*f_VoxY << endl;
  fhead << "FOV size (mm) [3] := " << ((PRECISION)f_DimZ)*f_VoxZ << endl;
  fhead << "image duration := " << f_Duration << endl;
  fhead << "image start time := " << f_TimeStart << endl;
  // Close file
  fhead.close();

  // Flip-back weight matrix
//  Misc_FlipXYZ(fp_Matrix,f_DimX,f_DimY,f_DimZ);

  // End
  return 0;
}

// ==========================================================================================================================================
// Function SaveFrame
//    This function saves a frame image coming from the end of an iteration.
//    It returns 0 upon success and other value otherwise.
// ==========================================================================================================================================
int oOutputManager::SaveFrame( PRECISION* fp_Image, int f_DimX, int f_DimY, int f_DimZ, PRECISION f_VoxX, PRECISION f_VoxY, PRECISION f_VoxZ, int f_Iter, int f_Frame,
                               int f_ScannerModel, PRECISION f_Duration, PRECISION f_TimeStart, const string& f_Isotope, PRECISION f_DCF1, PRECISION f_DCF2, int f_NbUpdates,
                               string f_Message )
{
  // Check
  if (!m_Directory)
  {
    LogCerr ("***** oOutputManager::SaveFrame() -> Output directory '" << m_PathName << m_RootName << "' is not valid !" << endl);
    return 1;
  }

  // Pre-stuffs
  int dim_tot = f_DimX * f_DimY * f_DimZ;

  // -----------------------------------------------------------------------------------------------------------------------------
  //                            Save data
  // -----------------------------------------------------------------------------------------------------------------------------

  // Open file
  char tmp_it[100]; if (f_Iter>=0) sprintf(tmp_it,"it%d",f_Iter+1); else sprintf(tmp_it,"bkprj");
  char tmp_fr[100]; sprintf(tmp_fr,"%d",f_Frame+1);
  string file_data = m_PathName+m_RootName+"/"+m_RootName+"_"+((string)tmp_it)+"_fr"+((string)tmp_fr)+".i";
  FILE* fout = fopen(file_data.c_str(),"wb");
  if (fout==NULL)
  {
    LogCerr ("***** oOutputManager::SaveFrame() -> Failed to create data file '" << file_data << "' !" << endl);
    return 1;
  }
  // Verbose
  if (m_Verbose>=1) LogCout ("          Write in file '" << file_data << "' ..." << endl);
  // Write data
  int nb_data_written = 0;
  for (int i=0; i<dim_tot; i++)
  {
    float buffer = (float)fp_Image[i];
    nb_data_written += fwrite(&buffer,sizeof(float),1,fout);
  }
  // Close file
  fclose(fout);
  // Check integrity
  if (nb_data_written!=dim_tot)
  {
    LogCerr ("***** oOutputManager::SaveFrame() -> Failed to write all data (" << dim_tot << ") in file '" << file_data << "' (" << nb_data_written << " written) !" << endl);
    return 1;
  }

  // -----------------------------------------------------------------------------------------------------------------------------
  //                            Save header as INTERFILE format
  // -----------------------------------------------------------------------------------------------------------------------------

  // Open file
  string file_head = m_PathName+m_RootName+"/"+m_RootName+"_"+((string)tmp_it)+"_fr"+((string)tmp_fr)+".i.hdr";
  ofstream fhead(file_head.c_str());
  if (!fhead)
  {
    LogCerr ("***** oOutputManager::SaveFrame() -> Failed to create header file '" << file_head << "' !" << endl);
    return 1;
  }
  // Write info
  fhead << "!INTERFILE" << endl;
  fhead << "!name of data file := " << m_RootName << "_" << tmp_it << "_fr" << tmp_fr << ".i" << endl;
  if (f_ScannerModel==SCANNER_HRPLUS) fhead << "!originating system := HR+" << endl;
  else if (f_ScannerModel==SCANNER_HRRT) fhead << "!originating system := HRRT" << endl;
  else if (f_ScannerModel==SCANNER_BIOGRAPH) fhead << "!originating system := Biograph" << endl;
  else if (f_ScannerModel==SCANNER_BIOGRAPH2D) fhead << "!originating system := Biograph2D" << endl;
  else if (f_ScannerModel==SCANNER_INVEON) fhead << "!originating system := Inveon" << endl;
  else if (f_ScannerModel==SCANNER_FOCUS) fhead << "!originating system := Focus" << endl;
  fhead << "!PET data type := emission" << endl;
  fhead << "data format := image" << endl;
  fhead << "number format := float" << endl;
  fhead << "number of bytes per pixel := 4" << endl;
  fhead << "image duration := " << f_Duration << endl;
  fhead << "image start time := " << f_TimeStart << endl;
  fhead << "Data units := Bq/cc" << endl;
  fhead << "Dose type := " << f_Isotope << endl;
  PRECISION half_life, branching_ratio;
  Misc_GetIsotopeCharacteristics(f_Isotope, &half_life, &branching_ratio);
  fhead << "isotope halflife := " << half_life << endl;
  fhead << "branching ratio := " << branching_ratio << endl;
  fhead << "Patient ID := " << f_Frame << endl;
  fhead << "Patient name := " << m_RootName << endl;
  fhead << "number of dimensions := 3" << endl;
  fhead << "matrix size [1] := " << f_DimX << endl;
  fhead << "matrix size [2] := " << f_DimY << endl;
  fhead << "matrix size [3] := " << f_DimZ << endl;
  fhead << "scaling factor (mm/pixel) [1] := " << f_VoxX << endl;
  fhead << "scaling factor (mm/pixel) [2] := " << f_VoxY << endl;
  fhead << "scaling factor (mm/pixel) [3] := " << f_VoxZ << endl;
  fhead << "FOV size (mm) [1] := " << ((PRECISION)f_DimX)*f_VoxX << endl;
  fhead << "FOV size (mm) [2] := " << ((PRECISION)f_DimY)*f_VoxY << endl;
  fhead << "FOV size (mm) [3] := " << ((PRECISION)f_DimZ)*f_VoxZ << endl;
  fhead << "frame := " << f_Frame << endl;
  fhead << "bed := 0" << endl;
  fhead << "decay correction factor  := " << f_DCF1 << endl;
  fhead << "decay correction factor2 := " << f_DCF2 << endl;
  if (f_Iter>=0) fhead << "reconstruction method := 3D static OSEM " << f_NbUpdates << " updates" << endl;
  else fhead << "reconstruction method := backprojection" << endl;
  fhead << f_Message << endl;
  // Close file
  fhead.close();

  // End
  return 0;
}
// ==========================================================================================================================================
// Function SaveDynFrame
//    This function saves an ascii file that contains the whole list of frames.
//    It returns 0 upon success and other value otherwise.
// ==========================================================================================================================================
int oOutputManager::SaveDynFrame(int f_NbFrames, int f_Iter, bool f_Smoothed)
{
  // Check
  if (!m_Directory)
  {
    LogCerr ("***** oOutputManager::SaveDynFrame() -> Output directory '" << m_PathName << m_RootName << "' is not valid !" << endl);
    return 1;
  }

  // Create file name
  char tmp_it[100]; if (f_Iter>=0) sprintf(tmp_it,"it%d",f_Iter+1); else sprintf(tmp_it,"bkprj");
  string file_name = "";
  if (f_Smoothed) file_name = m_PathName + m_RootName + "/" + m_RootName + "_" + ((string)tmp_it) + "_smth.dyn";
  else file_name = m_PathName + m_RootName + "/" + m_RootName + "_" + ((string)tmp_it) + ".dyn";

  // Open file
  ofstream fout(file_name.c_str());
  if (!fout)
  {
    LogCerr ("***** oOutputManager::SaveDynFrame() -> Failed to create output file '" << file_name << "' !" << endl);
    return 1;
  }

  // First line is the number of frames
  fout << f_NbFrames << endl;

  // Following lines are the list of image file names
  for (int fr=0; fr<f_NbFrames; fr++)
  {
    if (f_Smoothed) fout << m_RootName << "_" << tmp_it << "_fr" << fr+1 << "_smth.i.hdr" << endl;
    else fout << m_RootName << "_" << tmp_it << "_fr" << fr+1 << ".i.hdr" << endl;
  }

  // Close file
  fout.close();

  // End
  return 0;
}
// ==========================================================================================================================================
// Function SaveCoefficients
//    This function saves a coefficients image of a given time function, coming from the end of an iteration.
//    It returns 0 upon success and other value otherwise.
// ==========================================================================================================================================
int oOutputManager::SaveCoefficients( PRECISION* fp_Image, int f_DimX, int f_DimY, int f_DimZ, PRECISION f_VoxX, PRECISION f_VoxY, PRECISION f_VoxZ, int f_Iter, int f_TimeFunction,
                                      int f_ScannerModel, int f_NbUpdates, bool f_Smoothed, string f_Message )
{
  // Check
  if (!m_Directory)
  {
    LogCerr ("***** oOutputManager::SaveCoefficients() -> Output directory '" << m_PathName << m_RootName << "' is not valid !" << endl);
    return 1;
  }

  // Pre-stuffs
  int dim_tot = f_DimX * f_DimY * f_DimZ;

  // -----------------------------------------------------------------------------------------------------------------------------
  //                            Save data
  // -----------------------------------------------------------------------------------------------------------------------------

  // Open file
  char tmp_it[100]; sprintf(tmp_it,"%d",f_Iter+1);
  char tmp_tf[100]; sprintf(tmp_tf,"%d",f_TimeFunction+1);
  string file_data = "";
  if (f_Smoothed) file_data = m_PathName+m_RootName+"/"+m_RootName+"_it"+((string)tmp_it)+"_tf"+((string)tmp_tf)+"_smth.i";
  else file_data = m_PathName+m_RootName+"/"+m_RootName+"_it"+((string)tmp_it)+"_tf"+((string)tmp_tf)+".i";
  FILE* fout = fopen(file_data.c_str(),"wb");
  if (fout==NULL)
  {
    LogCerr ("***** oOutputManager::SaveCoefficients() -> Failed to create data file '" << file_data << "' !" << endl);
    return 1;
  }
  // Verbose
  if (m_Verbose>=1) LogCout ("          Write in file '" << file_data << "' ..." << endl);
  // Write data
  int nb_data_written = 0;
  for (int i=0; i<dim_tot; i++)
  {
    float buffer = (float)fp_Image[i];
    nb_data_written += fwrite(&buffer,sizeof(float),1,fout);
  }
  // Close file
  fclose(fout);
  // Check integrity
  if (nb_data_written!=dim_tot)
  {
    LogCerr ("***** oOutputManager::SaveCoefficients() -> Failed to write all data (" << dim_tot << ") in file '" << file_data << "' (" << nb_data_written << " written) !" << endl);
    return 1;
  }

  // -----------------------------------------------------------------------------------------------------------------------------
  //                            Save header as INTERFILE format
  // -----------------------------------------------------------------------------------------------------------------------------

  // Open file
  string file_head = "";
  if (f_Smoothed) file_head = m_PathName+m_RootName+"/"+m_RootName+"_it"+((string)tmp_it)+"_tf"+((string)tmp_tf)+"_smth.i.hdr";
  else file_head = m_PathName+m_RootName+"/"+m_RootName+"_it"+((string)tmp_it)+"_tf"+((string)tmp_tf)+".i.hdr";
  ofstream fhead(file_head.c_str());
  if (!fhead)
  {
    LogCerr ("***** oOutputManager::SaveCoefficients() -> Failed to create header file '" << file_head << "' !" << endl);
    return 1;
  }
  // Write info
  fhead << "!INTERFILE" << endl;
  if (f_Smoothed) fhead << "!name of data file := " << m_RootName << "_it" << tmp_it << "_tf" << tmp_tf << "_smth.i" << endl;
  else fhead << "!name of data file := " << m_RootName << "_it" << tmp_it << "_tf" << tmp_tf << ".i" << endl;
  if (f_ScannerModel==SCANNER_HRPLUS) fhead << "!originating system := HR+" << endl;
  else if (f_ScannerModel==SCANNER_HRRT) fhead << "!originating system := HRRT" << endl;
  else if (f_ScannerModel==SCANNER_BIOGRAPH) fhead << "!originating system := Biograph" << endl;
  else if (f_ScannerModel==SCANNER_BIOGRAPH2D) fhead << "!originating system := Biograph2D" << endl;
  else if (f_ScannerModel==SCANNER_INVEON) fhead << "!originating system := Inveon" << endl;
  else if (f_ScannerModel==SCANNER_FOCUS) fhead << "!originating system := Focus" << endl;
  fhead << "!PET data type := emission" << endl;
  fhead << "data format := image" << endl;
  fhead << "number format := float" << endl;
  fhead << "number of bytes per pixel := 4" << endl;
  fhead << "Data units := Arbitrary" << endl;
  fhead << "Patient ID := " << f_TimeFunction << endl;
  fhead << "Patient name := " << m_RootName << endl;
  fhead << "number of dimensions := 3" << endl;
  fhead << "matrix size [1] := " << f_DimX << endl;
  fhead << "matrix size [2] := " << f_DimY << endl;
  fhead << "matrix size [3] := " << f_DimZ << endl;
  fhead << "scaling factor (mm/pixel) [1] := " << f_VoxX << endl;
  fhead << "scaling factor (mm/pixel) [2] := " << f_VoxY << endl;
  fhead << "scaling factor (mm/pixel) [3] := " << f_VoxZ << endl;
  fhead << "FOV size (mm) [1] := " << ((PRECISION)f_DimX)*f_VoxX << endl;
  fhead << "FOV size (mm) [2] := " << ((PRECISION)f_DimY)*f_VoxY << endl;
  fhead << "FOV size (mm) [3] := " << ((PRECISION)f_DimZ)*f_VoxZ << endl;
  fhead << "time function := " << f_TimeFunction << endl;
  fhead << "bed := 0" << endl;
  fhead << "reconstruction method := 4D dynamic OSEM " << f_NbUpdates << " updates" << endl;
  fhead << "input function file name := " << m_RootName << "_input_function.txt" << endl;
  fhead << "time-basis functions file name := " << m_RootName << "_basis_functions.txt" << endl;
  fhead << "convolved time-basis functions file name := " << m_RootName << "_convolved_basis_functions.txt" << endl;
  fhead << f_Message << endl;
  // Close file
  fhead.close();

  // End
  return 0;
}
// ==========================================================================================================================================
// Function SaveDynCoefficients
//    This function saves an ascii file that contains the whole list of time functions.
//    It returns 0 upon success and other value otherwise.
// ==========================================================================================================================================
int oOutputManager::SaveDynCoefficients(int f_NbTimeFunctions, int f_Iter, bool f_Smoothed)
{
  // Check
  if (!m_Directory)
  {
    LogCerr ("***** oOutputManager::SaveDynCoefficients() -> Output directory '" << m_PathName << m_RootName << "' is not valid !" << endl);
    return 1;
  }

  // Create file name
  char tmp_it[100]; sprintf(tmp_it,"%d",f_Iter+1);
  string file_name = "";
  if (f_Smoothed) file_name = m_PathName + m_RootName + "/" + m_RootName + "_it" + ((string)tmp_it) + "_smth.tbf";
  else file_name = m_PathName + m_RootName + "/" + m_RootName + "_it" + ((string)tmp_it) + ".tbf";

  // Open file
  ofstream fout(file_name.c_str());
  if (!fout)
  {
    LogCerr ("***** oOutputManager::SaveDynCoefficients() -> Failed to create output file '" << file_name << "' !" << endl);
    return 1;
  }

  // First line is the number of time functions
  fout << f_NbTimeFunctions << endl;

  // Following lines are the list of image file names
  for (int tf=0; tf<f_NbTimeFunctions; tf++)
  {
    if (f_Smoothed) fout << m_RootName << "_it" << tmp_it << "_tf" << tf+1 << "_smth.i.hdr" << endl;
    else fout << m_RootName << "_it" << tmp_it << "_tf" << tf+1 << ".i.hdr" << endl;
  }

  // Close file
  fout.close();

  // End
  return 0;
}
// ==========================================================================================================================================
// Function SaveOrdinaryImage
//    This function saves an image without too much informations.
//    It returns 0 upon success and other value otherwise.
// ==========================================================================================================================================
int oOutputManager::SaveOrdinaryImage( float* fp_Image, int f_DimX, int f_DimY, int f_DimZ, float f_VoxX, float f_VoxY, float f_VoxZ,
                                       int f_ScannerModel, float f_Duration, float f_TimeStart, const string& f_Isotope )
{
  // Check
  if (!m_Directory)
  {
    LogCerr ("***** oOutputManager::SaveOrdinaryImage() -> Output directory '" << m_PathName << m_RootName << "' is not valid !" << endl);
    return 1;
  }

  // Pre-stuffs
  int dim_tot = f_DimX * f_DimY * f_DimZ;

  // -----------------------------------------------------------------------------------------------------------------------------
  //                            Save data
  // -----------------------------------------------------------------------------------------------------------------------------

  // Open file
  string file_data = m_PathName+m_RootName+"/"+m_RootName+"_emission.i";
  FILE* fout = fopen(file_data.c_str(),"wb");
  if (fout==NULL)
  {
    LogCerr ("***** oOutputManager::SaveOrdinaryImage() -> Failed to create data file '" << file_data << "' !" << endl);
    return 1;
  }
  // Verbose
  if (m_Verbose>=1) LogCout ("          Write in file '" << file_data << "' ..." << endl);
  // Write data
  int nb_data_written = 0;
  for (int i=0; i<dim_tot; i++)
  {
    float buffer = (float)fp_Image[i];
    nb_data_written += fwrite(&buffer,sizeof(float),1,fout);
  }
  // Close file
  fclose(fout);
  // Check integrity
  if (nb_data_written!=dim_tot)
  {
    LogCerr ("***** oOutputManager::SaveOrdinaryImage() -> Failed to write all data (" << dim_tot << ") in file '" << file_data << "' (" << nb_data_written << " written) !" << endl);
    return 1;
  }

  // -----------------------------------------------------------------------------------------------------------------------------
  //                            Save header as INTERFILE format
  // -----------------------------------------------------------------------------------------------------------------------------

  // Open file
  string file_head = m_PathName+m_RootName+"/"+m_RootName+"_emission.i.hdr";
  ofstream fhead(file_head.c_str());
  if (!fhead)
  {
    LogCerr ("***** oOutputManager::SaveFrame() -> Failed to create header file '" << file_head << "' !" << endl);
    return 1;
  }
  // Write info
  fhead << "!INTERFILE" << endl;
  fhead << "!name of data file := " << m_RootName << "_emission.i" << endl;
  if (f_ScannerModel==SCANNER_HRPLUS) fhead << "!originating system := HR+" << endl;
  else if (f_ScannerModel==SCANNER_HRRT) fhead << "!originating system := HRRT" << endl;
  else if (f_ScannerModel==SCANNER_BIOGRAPH) fhead << "!originating system := Biograph" << endl;
  else if (f_ScannerModel==SCANNER_BIOGRAPH2D) fhead << "!originating system := Biograph2D" << endl;
  else if (f_ScannerModel==SCANNER_INVEON) fhead << "!originating system := Inveon" << endl;
  else if (f_ScannerModel==SCANNER_FOCUS) fhead << "!originating system := Focus" << endl;
  fhead << "!PET data type := emission" << endl;
  fhead << "data format := image" << endl;
  fhead << "number format := float" << endl;
  fhead << "number of bytes per pixel := 4" << endl;
  fhead << "image duration := " << f_Duration << endl;
  fhead << "image start time := " << f_TimeStart << endl;
  fhead << "Data units := Bq/cc" << endl;
  fhead << "Dose type := " << f_Isotope << endl;
  PRECISION half_life, branching_ratio;
  Misc_GetIsotopeCharacteristics(f_Isotope, &half_life, &branching_ratio);
  fhead << "isotope halflife := " << half_life << endl;
  fhead << "branching ratio := " << branching_ratio << endl;
  fhead << "Patient name := " << m_RootName << endl;
  fhead << "number of dimensions := 3" << endl;
  fhead << "matrix size [1] := " << f_DimX << endl;
  fhead << "matrix size [2] := " << f_DimY << endl;
  fhead << "matrix size [3] := " << f_DimZ << endl;
  fhead << "scaling factor (mm/pixel) [1] := " << f_VoxX << endl;
  fhead << "scaling factor (mm/pixel) [2] := " << f_VoxY << endl;
  fhead << "scaling factor (mm/pixel) [3] := " << f_VoxZ << endl;
  fhead << "FOV size (mm) [1] := " << ((float)f_DimX)*f_VoxX << endl;
  fhead << "FOV size (mm) [2] := " << ((float)f_DimY)*f_VoxY << endl;
  fhead << "FOV size (mm) [3] := " << ((float)f_DimZ)*f_VoxZ << endl;
  // Close file
  fhead.close();

  // End
  return 0;
}
// ==========================================================================================================================================
// Function SaveOrdinaryMuMap
//    This function saves a mumap without too much informations.
//    It returns 0 upon success and other value otherwise.
// ==========================================================================================================================================
int oOutputManager::SaveOrdinaryMuMap( float* fp_Image, int f_DimX, int f_DimY, int f_DimZ, float f_VoxX, float f_VoxY, float f_VoxZ,
                                       int f_ScannerModel, float f_Duration, float f_TimeStart, const string& f_Isotope )
{
  // Check
  if (!m_Directory)
  {
    LogCerr ("***** oOutputManager::SaveOrdinaryImage() -> Output directory '" << m_PathName << m_RootName << "' is not valid !" << endl);
    return 1;
  }

  // Pre-stuffs
  int dim_tot = f_DimX * f_DimY * f_DimZ;

  // -----------------------------------------------------------------------------------------------------------------------------
  //                            Save data
  // -----------------------------------------------------------------------------------------------------------------------------

  // Open file
  string file_data = m_PathName+m_RootName+"/"+m_RootName+"_transmission.i";
  FILE* fout = fopen(file_data.c_str(),"wb");
  if (fout==NULL)
  {
    LogCerr ("***** oOutputManager::SaveOrdinaryImage() -> Failed to create data file '" << file_data << "' !" << endl);
    return 1;
  }
  // Verbose
  if (m_Verbose>=1) LogCout ("          Write in file '" << file_data << "' ..." << endl);
  // Write data
  int nb_data_written = 0;
  for (int i=0; i<dim_tot; i++)
  {
    float buffer = (float)fp_Image[i];
    nb_data_written += fwrite(&buffer,sizeof(float),1,fout);
  }
  // Close file
  fclose(fout);
  // Check integrity
  if (nb_data_written!=dim_tot)
  {
    LogCerr ("***** oOutputManager::SaveOrdinaryImage() -> Failed to write all data (" << dim_tot << ") in file '" << file_data << "' (" << nb_data_written << " written) !" << endl);
    return 1;
  }

  // -----------------------------------------------------------------------------------------------------------------------------
  //                            Save header as INTERFILE format
  // -----------------------------------------------------------------------------------------------------------------------------

  // Open file
  string file_head = m_PathName+m_RootName+"/"+m_RootName+"_transmission.i.hdr";
  ofstream fhead(file_head.c_str());
  if (!fhead)
  {
    LogCerr ("***** oOutputManager::SaveFrame() -> Failed to create header file '" << file_head << "' !" << endl);
    return 1;
  }
  // Write info
  fhead << "!INTERFILE" << endl;
  fhead << "!name of data file := " << m_RootName << "_transmission.i" << endl;
  if (f_ScannerModel==SCANNER_HRPLUS) fhead << "!originating system := HR+" << endl;
  else if (f_ScannerModel==SCANNER_HRRT) fhead << "!originating system := HRRT" << endl;
  else if (f_ScannerModel==SCANNER_BIOGRAPH) fhead << "!originating system := Biograph" << endl;
  else if (f_ScannerModel==SCANNER_BIOGRAPH2D) fhead << "!originating system := Biograph2D" << endl;
  else if (f_ScannerModel==SCANNER_INVEON) fhead << "!originating system := Inveon" << endl;
  else if (f_ScannerModel==SCANNER_FOCUS) fhead << "!originating system := Focus" << endl;
  fhead << "!PET data type := transmission" << endl;
  fhead << "data format := image" << endl;
  fhead << "number format := float" << endl;
  fhead << "number of bytes per pixel := 4" << endl;
  fhead << "image duration := " << f_Duration << endl;
  fhead << "image start time := " << f_TimeStart << endl;
  fhead << "Data units := Bq/cc" << endl;
  fhead << "Dose type := " << f_Isotope << endl;
  PRECISION half_life, branching_ratio;
  Misc_GetIsotopeCharacteristics(f_Isotope, &half_life, &branching_ratio);
  fhead << "isotope halflife := " << half_life << endl;
  fhead << "branching ratio := " << branching_ratio << endl;
  fhead << "Patient name := " << m_RootName << endl;
  fhead << "number of dimensions := 3" << endl;
  fhead << "matrix size [1] := " << f_DimX << endl;
  fhead << "matrix size [2] := " << f_DimY << endl;
  fhead << "matrix size [3] := " << f_DimZ << endl;
  fhead << "scaling factor (mm/pixel) [1] := " << f_VoxX << endl;
  fhead << "scaling factor (mm/pixel) [2] := " << f_VoxY << endl;
  fhead << "scaling factor (mm/pixel) [3] := " << f_VoxZ << endl;
  fhead << "FOV size (mm) [1] := " << ((float)f_DimX)*f_VoxX << endl;
  fhead << "FOV size (mm) [2] := " << ((float)f_DimY)*f_VoxY << endl;
  fhead << "FOV size (mm) [3] := " << ((float)f_DimZ)*f_VoxZ << endl;
  fhead << "patient orientation := NOPITCH" << endl;
  // Close file
  fhead.close();

  // End
  return 0;
}
// ==========================================================================================================================================
// Function SaveInputFunction
//    This function saves an input function.
//    It returns 0 upon success and other value otherwise.
// ==========================================================================================================================================
int oOutputManager::SaveInputFunction(int f_NbFrames, PRECISION* fp_FrameTimeStart, PRECISION* fp_FrameTimeStop, PRECISION* fp_FrameDuration, PRECISION* fp_InputFunction)
{
  // Check
  if (!m_Directory)
  {
    LogCerr ("***** oOutputManager::SaveInputFunction() -> Output directory '" <<m_PathName << m_RootName << "' is not valid !" << endl);
    return 1;
  }

  // Build file name and open file
  string file_name = m_PathName + m_RootName + "/" + m_RootName + "_input_function.txt";
  ofstream fout(file_name.c_str());
  if (!fout)
  {
    LogCerr ("***** oOutputManager::SaveInputFunction() -> Failed to create output file '" << file_name << "' !" << endl);
    return 1;
  }
  // Verbose
  if (m_Verbose>=1) LogCout ("  --> Saving input function as '" << file_name << "' ..." << endl);

  // Register values
  for (int fr=0; fr<f_NbFrames; fr++)
  {
    fout << fr+1 << "\t" << fp_FrameTimeStart[fr] << "\t" << fp_FrameTimeStop[fr] << "\t" << fp_FrameDuration[fr] << "\t" << fp_InputFunction[fr] << endl;
  }

  // Register column labels
  fout << endl;
  fout << "Frame #\tStart (sec)\tStop (sec)\tDuration (sec)\tValue" << endl;
  // Close file
  fout.close();

  // End
  return 0;
}
// ==========================================================================================================================================
// Function SaveBasisFunctions
//    This function saves the basis functions
//    It returns 0 upon success and other value otherwise.
// ==========================================================================================================================================
int oOutputManager::SaveBasisFunctions(int f_NbFrames, int f_NbTimeFunctions, PRECISION* fp_FrameTimeStart, PRECISION* fp_FrameTimeStop, PRECISION* fp_FrameDuration, PRECISION** fp_BasisFunctions)
{
  // Check
  if (!m_Directory)
  {
    LogCerr ("***** oOutputManager::SaveBasisFunctions() -> Output directory '" << m_PathName << m_RootName << "' is not valid !" << endl);
    return 1;
  }

  // Build file name and open file
  string file_name = m_PathName + m_RootName + "/" + m_RootName + "_basis_functions.txt";
  ofstream fout(file_name.c_str());
  if (!fout)
  {
    LogCerr ("***** oOutputManager::SaveBasisFunctions() -> Failed to create output file '" << file_name << "' !" << endl);
    return 1;
  }
  // Verbose
  if (m_Verbose>=1) LogCout ("  --> Saving time-basis functions as '" << file_name << "' ..." << endl);

  // Register values
  for (int fr=0; fr<f_NbFrames; fr++)
  {
    fout << fr+1 << "\t" << fp_FrameTimeStart[fr] << "\t" << fp_FrameTimeStop[fr] << "\t" << fp_FrameDuration[fr];
    for (int tf=0; tf<f_NbTimeFunctions; tf++) fout << "\t" << fp_BasisFunctions[tf][fr];
    fout << endl;
  }

  // Register column labels
  fout << endl;
  fout << "Frame #\tStart (sec)\tStop (sec)\tDuration (sec)";
  for (int tf=0; tf<f_NbTimeFunctions; tf++) fout << "\tValue" << tf+1;
  fout << endl;

  // Close file
  fout.close();

  // End
  return 0;
}
// ==========================================================================================================================================
// Function SaveConvolvedBasisFunctions
//    This function saves the basis functions
//    It returns 0 upon success and other value otherwise.
// ==========================================================================================================================================
int oOutputManager::SaveConvolvedBasisFunctions(int f_NbFrames, int f_NbTimeFunctions, PRECISION* fp_FrameTimeStart, PRECISION* fp_FrameTimeStop, PRECISION* fp_FrameDuration, PRECISION** fp_ConvolvedBasisFunctions)
{
  // Check
  if (!m_Directory)
  {
    LogCerr ("***** oOutputManager::SaveConvolvedBasisFunctions() -> Output directory '" << m_PathName << m_RootName << "' is not valid !" << endl);
    return 1;
  }

  // Build file name and open file
  string file_name = m_PathName + m_RootName + "/" + m_RootName + "_convolved_basis_functions.txt";
  ofstream fout(file_name.c_str());
  if (!fout)
  {
    LogCerr ("***** oOutputManager::SaveConvolvedBasisFunctions() -> Failed to create output file '" << file_name << "' !" << endl);
    return 1;
  }
  // Verbose
  if (m_Verbose>=1) LogCout ("  --> Saving convolved time-basis functions as '" << file_name << "' ..." << endl);

  // Register values
  for (int fr=0; fr<f_NbFrames; fr++)
  {
    fout << fr+1 << "\t" << fp_FrameTimeStart[fr] << "\t" << fp_FrameTimeStop[fr] << "\t" << fp_FrameDuration[fr];
    for (int tf=0; tf<f_NbTimeFunctions; tf++) fout << "\t" << fp_ConvolvedBasisFunctions[tf][fr];
    fout << endl;
  }

  // Register column labels
  fout << endl;
  fout << "Frame #\tStart (sec)\tStop (sec)\tDuration (sec)";
  for (int tf=0; tf<f_NbTimeFunctions; tf++) fout << "\tValue" << tf+1;
  fout << endl;

  // Close file
  fout.close();

  // End
  return 0;
}

