#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <ctime>
#include "oSinoModeCreation.hh"
#include "oOutputManager.hh"
using namespace std;

void showHelp(int code)
{
  cout << endl;
  cout << "Usage:  SMmaker  -m scannerName  -t sino_tr.s.hdr  -o output_root" << endl;
  cout << "                 [-l listMode.l.hdr] [-p sino_pt.s.hdr] [-r sino_de.s.hdr] [-f sino_float.s.hdr]" << endl;
  cout << "                 [-s sino_sc.s.hdr] [-n sino_nm.s.hdr] [-a mumap.i.hdr] [-c crystal_map.ecm]" << endl;
  cout << "                 [-w weightMode] [-S seed] [-v verbose] [-F] [-2D] [-D]" << endl;
  cout << "                 [-la percent_angle]" << endl;
  #ifdef OMP_MODE
  cout << "                 [-T nbThread]" << endl;
  #endif
  cout << endl;
  cout << "   [Options]:" << endl;
  cout << "   -m  value: give the scanner name ('hrplus', 'biograph', 'inveon', 'focus', 'mmr2d' or 'vision600' for the moment)" << endl;
  cout << "   -t  value: give the net trues sinogram header file name" << endl;
  cout << "   -p  value: give the prompts sinogram header file name (alternative to -t, -l or -f)" << endl;
  cout << "   -r  value: give the randoms sinogram header file name (useful when -p or -l is used)" << endl;
  cout << "   -f  value: give a float sinogram assumed to be precorrected, except for normalization (alternative to -t or -l or -p)" << endl;
  cout << "   -l  value: give the list-mode header file name (alternative to -t or -p)" << endl;
  cout << "   -s  value: give the scatter sinogram header file name (also requires -n)" << endl;
  cout << "              (for the focus scanner, can give also the log file of the e7_sino program to get the actual scatter fraction, give this file separated by a comma)" << endl;
  cout << "   -n  value: give the normalisation sinogram header file name" << endl;
  cout << "   -a  value: give the file name of the mu-map interfile header (have to provide a crystal map -c) (alternative to -A)" << endl;
  cout << "   -A  value: give the attenuation sinogram data file (supposed to have same dimensions as other sinograms) (alternative to -a)" << endl;
  cout << "   -c  value: give the file name of the esteban crystal map corresponding to the scanner model" << endl;
  cout << "   -o  value: give the root output file name" << endl;
  cout << "   -w  value: give the writing mode (0: only positive counts, 1: all counts) (default: 1)" << endl;
  cout << "   -F       : do not fill empty LORs to get the same number of LORs in each event" << endl;
  cout << "   -d  value: give the delay post-acquisition in seconds, if any (default: 0)" << endl;
  cout << "   -2D      : set in 2D mode = assume span 1, ring diff 1" << endl;
  cout << "   -castor  : write the output in castor format" << endl;
  cout << "   -D       : use Didier's implementation of siddon for forward projection of the mumap" << endl;
  cout << "   -la value: limit the azymuthal covering to the provided percent of angles (to simulate limited angle tomography)" << endl;
  #ifdef OMP_MODE
  cout << "   -T  value: give the number of threads for computation (default: 1)" << endl;
  #endif
  cout << "   -S  value: random generator seed (negative seeds will take time and pid)" << endl;
  cout << "   -v  value: give the verbose level" << endl;
  cout << "   -h       : print this help page and quit" << endl;
  cout << endl;
  cout << "  This program is used to create generic list-mode files from manufacturer sinograms." << endl;
  cout << "  The generic list-mode files can then be reconstructed with the MLrecon program." << endl;
  cout << endl;
  #ifdef OMP_MODE
  cout << "  Compiled with OpenMP" << endl;
  cout << endl;
  #endif
  #ifdef BUILD_DATE
  cout << "  Build date: " << BUILD_DATE << endl;
  cout << endl;
  #endif
  exit(code);
}

int main(int argc, char** argv)
{

  if (argc==1) showHelp(0);
  cout << endl;

  // ============================================================================================================
  //                                             INITIALIZATION
  // ============================================================================================================

  // Parameters
  string fileBaseOut = "";
  string scannerName = "";
  string fileFloat   = "";
  string fileNetTrue  = "";
  string filePrompt   = "";
  string fileRandom   = "";
  string fileScatter = "";
  string fileListMode = "";
  string fileAttenuation = "";
  int attenuationMode = -1;
  string fileCrystalMap = "";
  string fileNormalization = "";
  bool fillEqualLORs = true;
  int threads = 1;
  int weighting = ALL_WEIGHTING;
  int seed = -1;
  int verbose = 0;
  PRECISION timeDelay = 0.;
  bool mode2D = false;
  bool castor = false;
  bool siddonDidier = false;
  PRECISION limited_angle_percent = 100.;

  // Read command-line arguments
  for (int i=1; i<argc; i++)
  {
    string option = (string)argv[i];
    if (option=="-m")
    {
      if (i>=argc-1)
      {
        cerr << "***** Not enaugh arguments after option " << option << " !" << endl;
        showHelp(1);
      }
      scannerName = (string)argv[i+1];
      i++;
    }
    else if (option=="-o")
    {
      if (i>=argc-1)
      {
        cerr << "***** Not enaugh arguments after option " << option << " !" << endl;
        showHelp(1);
      }
      fileBaseOut = (string)argv[i+1];
      i++;
    }
    else if (option=="-s")
    {
      if (i>=argc-1)
      {
        cerr << "***** Not enaugh arguments after option " << option << " !" << endl;
        showHelp(1);
      }
      fileScatter = (string)argv[i+1];
      i++;
    }
    else if (option=="-t")
    {
      if (i>=argc-1)
      {
        cerr << "***** Not enaugh arguments after option " << option << " !" << endl;
        showHelp(1);
      }
      fileNetTrue = (string)argv[i+1];
      i++;
    }
    else if (option=="-f")
    {
      if (i>=argc-1)
      {
        cerr << "***** Not enaugh arguments after option " << option << " !" << endl;
        showHelp(1);
      }
      fileFloat = (string)argv[i+1];
      i++;
    }
    else if (option=="-p")
    {
      if (i>=argc-1)
      {
        cerr << "***** Not enaugh arguments after option " << option << " !" << endl;
        showHelp(1);
      }
      filePrompt = (string)argv[i+1];
      i++;
    }
    else if (option=="-r")
    {
      if (i>=argc-1)
      {
        cerr << "***** Not enaugh arguments after option " << option << " !" << endl;
        showHelp(1);
      }
      fileRandom = (string)argv[i+1];
      i++;
    }
    else if (option=="-l")
    {
      if (i>=argc-1)
      {
        cerr << "***** Not enaugh arguments after option " << option << " !" << endl;
        showHelp(1);
      }
      fileListMode = (string)argv[i+1];
      i++;
    }
    else if (option=="-n")
    {
      if (i>=argc-1)
      {
        cerr << "***** Not enaugh arguments after option " << option << " !" << endl;
        showHelp(1);
      }
      fileNormalization = (string)argv[i+1];
      i++;
    }
    else if (option=="-a")
    {
      if (i>=argc-1)
      {
        cerr << "***** Not enaugh arguments after option " << option << " !" << endl;
        showHelp(1);
      }
      fileAttenuation = (string)argv[i+1];
      attenuationMode = ATTENUATION_FROM_UMAP;
      i++;
    }
    else if (option=="-A")
    {
      if (i>=argc-1)
      {
        cerr << "***** Not enaugh arguments after option " << option << " !" << endl;
        showHelp(1);
      }
      fileAttenuation = (string)argv[i+1];
      attenuationMode = ATTENUATION_FROM_SINO;
      i++;
    }
    else if (option=="-c")
    {
      if (i>=argc-1)
      {
        cerr << "***** Not enaugh arguments after option " << option << " !" << endl;
        showHelp(1);
      }
      fileCrystalMap = (string)argv[i+1];
      i++;
    }
    else if (option=="-v")
    {
      if (i>=argc-1)
      {
        cerr << "***** Not enaugh arguments after option " << option << " !" << endl;
        showHelp(1);
      }
      verbose = atoi(argv[i+1]);
      i++;
    }
    else if (option=="-w")
    {
      if (i>=argc-1)
      {
        cerr << "***** Not enaugh arguments after option " << option << " !" << endl;
        showHelp(1);
      }
      weighting = atoi(argv[i+1]);
      i++;
    }
    else if (option=="-d")
    {
      if (i>=argc-1)
      {
        cerr << "***** Not enaugh arguments after option " << option << " !" << endl;
        showHelp(1);
      }
      timeDelay = (PRECISION)atof(argv[i+1]);
      i++;
    }
    else if (option=="-la")
    {
      if (i>=argc-1)
      {
        cerr << "***** Not enaugh arguments after option " << option << " !" << endl;
        showHelp(1);
      }
      limited_angle_percent = (PRECISION)atof(argv[i+1]);
      i++;
    }
    #ifdef OMP_MODE
    else if (option=="-T")
    {
      if (i>=argc-1)
      {
        cerr << "***** Not enaugh arguments after option " << option << " !" << endl;
        showHelp(1);
      }
      threads = atoi(argv[i+1]);
      i++;
    }
    #endif
    else if (option=="-S")
    {
      if (i>=argc-1)
      {
        cerr << "***** Not enaugh arguments after option " << option << " !" << endl;
        showHelp(1);
      }
      seed = atoi(argv[i+1]);
      i++;
    }
    else if (option=="-F") fillEqualLORs = false;
    else if (option=="-2D") mode2D = true;
    else if (option=="-castor") castor = true;
    else if (option=="-D") siddonDidier = true;
    else if (option=="-h" || option=="--help") showHelp(0);
    else
    {
      cerr << "***** Unknown option " << option << " !" << endl;
      showHelp(1);
    }
  }

  // Checks
  if (scannerName=="")
  {
    cerr << "***** Please provide a scanner name !" << endl;
    showHelp(1);
  }
  if (fileBaseOut=="")
  {
    cerr << "***** Please provide a root output file name !" << endl;
    showHelp(1);
  }
  if (fileNetTrue=="" && fileListMode=="" && filePrompt=="" && fileFloat=="")
  {
    cerr << "***** Please provide an input sinogram or list-mode file !" << endl;
    showHelp(1);
  }
  if (fileAttenuation!="" && fileCrystalMap=="")
  {
    cerr << "***** Cannot correct for attenuation (have to provide an esteban crystal map) !" << endl;
    showHelp(1);
  }
  if (fileNormalization=="" && fileFloat=="")
  {
    cerr << "***** Cannot correct for scatter or get the ECF (have to provide a normalisation sinogram) !" << endl;
    showHelp(1);
  }
  if (weighting!=POS_WEIGHTING && weighting!=ALL_WEIGHTING)
  {
    cerr << "***** Must choose an appropriate writing mode !" << endl;
    showHelp(1);
  }
  // For inveon or focus, this is not possible now to use prompts and delays, only net trues
  if ((scannerName=="inveon" || scannerName=="focus") && (fileNetTrue=="" &&fileFloat==""))
  {
    cerr << "***** For inveon and focus scanners, only the net true sinogram is supported for the moment !" << endl;
    return 1;
  }
  // For inveonTR or focus, 2D mode is possible but not for other scanners
  if ((scannerName!="inveon" && scannerName!="focus") && mode2D)
  {
    cerr << "***** 2D mode is only available for inveon and focus scanners !" << endl;
    return 1;
  }

  // ============================================================================================================
  //                                           CREATION OF THE GENERIC LIST MODE
  // ============================================================================================================

  // Create the output manager and log command
  oOutputManager* output = new oOutputManager(fileBaseOut, verbose);
  output->LogCommandLine(argc,argv);

  // Creation of the object
  oSinoModeCreation* sinomode = new oSinoModeCreation(scannerName, weighting, siddonDidier, seed, verbose);

  // Initialize input files
  if (sinomode->InitInputFiles(fileFloat,fileNetTrue,filePrompt,fileListMode,fileRandom,fileScatter,fileNormalization,fileAttenuation,attenuationMode,timeDelay,mode2D))
  {
    LogCerr ("***** Error while initializing input sinograms !" << endl);
    exit(1);
  }

  // Initialize conversion table
  if (sinomode->InitConversionTable())
  {
    LogCerr ("***** Error while initializing output sinogram !" << endl);
    exit(1);
  }

  // Compute total random amount
  if (fileRandom!="" && sinomode->ComputeTotalRandom())
  {
    LogCerr ("***** Error while computing total random counts !" << endl);
    exit(1);
  }

  // Compute total scatter amount
  if (fileScatter!="" && sinomode->ComputeTotalScatter())
  {
    LogCerr ("***** Error while computing total scatter counts !" << endl);
    exit(1);
  }

  // Initialize crystal map
  if (fileCrystalMap!="" && sinomode->InitCrystalMap(fileCrystalMap, castor))
  {
    LogCerr ("***** Error while initializing u-map !" << endl);
    exit(1);
  }

  // Process listmode file
  if (sinomode->ProcessListMode(fileBaseOut, fillEqualLORs, threads, castor, limited_angle_percent))
  {
    LogCerr ("***** Error while processing the list-mode file !" << endl);
    exit(1);
  }

  // ============================================================================================================
  //                                           THIS IS THE END MY FRIEND !
  // ============================================================================================================

  // Destroy the objects
  delete sinomode;
  delete output;

  // Ending
  cout << endl;
  return 0;
}

