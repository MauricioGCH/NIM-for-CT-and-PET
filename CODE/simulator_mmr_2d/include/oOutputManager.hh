#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <ctime>
#include <sstream>
#ifdef MPI_MODE
#include <mpi.h>
#endif
#include "oMiscellaneous.hh"
using namespace std;

#ifndef OOUTPUTMANAGER_HH
#define OOUTPUTMANAGER_HH 1

#ifdef MPI_MODE

// Macros for logging messages on standard output (with MPI)
#define LogCout(MESSAGE)                                        \
  do                                                            \
  {                                                             \
    int mpi_rank; MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);     \
    if (mpi_rank==0)                                            \
    {                                                           \
      std::cout << MESSAGE;                                     \
      oOutputManager* instance = oOutputManager::GetInstance(); \
      if (instance!=NULL)                                       \
      {                                                         \
        ofstream& logMac = instance->GetLog();                  \
        if (logMac) logMac << MESSAGE;                          \
      }                                                         \
    }                                                           \
  } while(0)
// Macros for logging messages on standard error (with MPI)
#define LogCerr(MESSAGE)                                        \
  do                                                            \
  {                                                             \
    int mpi_rank; MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);     \
    if (mpi_rank==0)                                            \
    {                                                           \
      std::cerr << MESSAGE;                                     \
      oOutputManager* instance = oOutputManager::GetInstance(); \
      if (instance!=NULL)                                       \
      {                                                         \
        ofstream& logMac = instance->GetLog();                  \
        if (logMac) logMac << MESSAGE;                          \
      }                                                         \
    }                                                           \
  } while(0)

#else

// Macros for logging messages on standard output (without MPI)
#define LogCout(MESSAGE)                                      \
  do                                                          \
  {                                                           \
    std::cout << MESSAGE;                                     \
    oOutputManager* instance = oOutputManager::GetInstance(); \
    if (instance!=NULL)                                       \
    {                                                         \
      ofstream& logMac = instance->GetLog();                  \
      if (logMac) logMac << MESSAGE;                          \
    }                                                         \
  } while(0)
// Macros for logging messages on standard error (without MPI)
#define LogCerr(MESSAGE)                                      \
  do                                                          \
  {                                                           \
    std::cerr << MESSAGE;                                     \
    oOutputManager* instance = oOutputManager::GetInstance(); \
    if (instance!=NULL)                                       \
    {                                                         \
      ofstream& logMac = instance->GetLog();                  \
      if (logMac) logMac << MESSAGE;                          \
    }                                                         \
  } while(0)

#endif

#define NO_LOG_KEYWORD "___NO_LOG___"

class oOutputManager
{
  // Constructor & Destructor
  public:
    oOutputManager( const string& f_RootName, int f_Verbose );
    ~oOutputManager();
    // Singleton implementation (but a first initialization is required)
    static oOutputManager* GetInstance()
    {
      if (!mp_Instance)
      {
        /* Maybe too fashist, just return NULL instead, and assume one tests its existency before using it...
        cerr << "***** oOutputManager() -> Try to get an instance while never initialized !" << endl;
        cerr << "                          The output manager must be used." << endl;
        exit(1);
        */
        return NULL;
      }
      return mp_Instance;
    }
    static oOutputManager* mp_Instance;

  // Get & Set functions
  public:
    ofstream& GetLog() {return m_Log;}
    inline const string& GetRootName() {return m_RootName;}
    inline const string& GetPathName() {return m_PathName;}

  // Member's functions
  public:
    int InitOutputDirectory();
    int InitLogFile();
    int LogCommandLine(int argc, char** argv);
    int LogCPUTime();
    int SaveWeightMatrix( PRECISION* fp_Matrix, int f_DimX, int f_DimY, int f_DimZ, PRECISION f_VoxX, PRECISION f_VoxY, PRECISION f_VoxZ, PRECISION f_OrigX, PRECISION f_OrigY, PRECISION f_OrigZ, PRECISION f_TimeStart, PRECISION f_Duration, int f_ScannerModel );
    int SaveOrdinaryImage( float* fp_Image, int f_DimX, int f_DimY, int f_DimZ, float f_VoxX, float f_VoxY, float f_VoxZ,
                           int f_ScannerModel, float f_Duration, float f_TimeStart, const string& f_Isotope );
    int SaveOrdinaryMuMap( float* fp_Image, int f_DimX, int f_DimY, int f_DimZ, float f_VoxX, float f_VoxY, float f_VoxZ,
                           int f_ScannerModel, float f_Duration, float f_TimeStart, const string& f_Isotope );
    int SaveFrame( PRECISION* fp_Image, int f_DimX, int f_DimY, int f_DimZ, PRECISION f_VoxX, PRECISION f_VoxY, PRECISION f_VoxZ, int f_Iter, int f_Frame,
                   int f_ScannerModel, PRECISION f_Duration, PRECISION f_TimeStart, const string& f_Isotope, PRECISION f_DCF1, PRECISION f_DCF2, int f_NbUpdates,
                   string f_Message = "" );
    int SaveCoefficients( PRECISION* fp_Image, int f_DimX, int f_DimY, int f_DimZ, PRECISION f_VoxX, PRECISION f_VoxY, PRECISION f_VoxZ, int f_Iter, int f_TimeFunction,
                            int f_ScannerModel, int f_NbUpdates, bool f_Smoothed, string f_Message = "" );
    int SaveInputFunction(int f_NbFrames, PRECISION* fp_FrameTimeStart, PRECISION* fp_FrameTimeStop, PRECISION* fp_FrameDuration, PRECISION* fp_InputFunction);
    int SaveBasisFunctions(int f_NbFrames, int f_NbTimeFunctions, PRECISION* fp_FrameTimeStart, PRECISION* fp_FrameTimeStop, PRECISION* fp_FrameDuration, PRECISION** fp_BasisFunctions);
    int SaveConvolvedBasisFunctions(int f_NbFrames, int f_NbTimeFunctions, PRECISION* fp_FrameTimeStart, PRECISION* fp_FrameTimeStop, PRECISION* fp_FrameDuration, PRECISION** fp_ConvolvedBasisFunctions);
    int SaveDynFrame(int f_NbFrames, int f_Iter, bool f_Smoothed);
    int SaveDynCoefficients(int f_NbTimeFunctions, int f_Iter, bool f_Smoothed);

  // Member's data
  private:
    ofstream m_Log;
    bool m_Directory;
    string m_RootName;
    string m_PathName;
    // Verbosity
    int m_Verbose;
};

#endif

