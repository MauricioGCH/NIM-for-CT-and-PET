#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <ctime>
#include "oOutputManager.hh"
#include "oSimulator.hh"
using namespace std;

void showHelp(int code)
{
  cout << endl;
  cout << "Usage:  SMprojector  -m scannerName  -c crystal_map.ecm  -o output_root  -i image.i.hdr" << endl;
  cout << "                     [-s scatterFraction] [-r randomFraction] [-a mumap.i.hdr]" << endl;
  cout << "                     [-P nbCounts] [-E fixedECF] [-L fill]" << endl;
  cout << "                     [-M mashing] [-S span] [-R maxringdiff]" << endl;
  cout << "                     [-p psfTrans,psfAxial | -p psfIso]" << endl;
  cout << "                     [-f offX offY offZ]" << endl;
  cout << "                     [-v verbose] [-w] [-z seed]" << endl;
  cout << "                     [--noSino]" << endl;
  cout << "                     [-D]" << endl;
  #ifdef OMP_MODE
  cout << "                     [-T nbThread]" << endl;
  #endif
  cout << endl;
  cout << "   [Options]:" << endl;
  cout << "   -m  value: give the scanner name ('hrplus', 'inveon', 'mmr2d', 'biograph' or 'vision600' for the moment)" << endl;
  cout << "   -c  value: give the file name of the crystal map corresponding to the scanner model" << endl;
  cout << "   -i  value: give the file name of the header of the image to be projected" << endl;
  cout << "   -f  value: give the offset of the image in the field-of-view in mm (default: image centered 0. 0. 0.)" << endl;
  cout << "   -o  value: give the root output file name" << endl;
  cout << "   -p  value: give the FWHM in mm of the 3D gaussian PSF (can give 2 values separated by a comma for transaxial and axial FWHM)" << endl;
  cout << "   -s  value: include a scatter component from the given scatter fraction (scat/net_true)" << endl;
  cout << "   -r  value: include a random component from the given random fraction (rand/prompt)" << endl;
  cout << "   -l  value: ratio of LSO flat randoms in the total random component (default: 0.1)" << endl;
  cout << "   -a  value: include the attenuation effect from the given mumap in cm-1 (except for inveon in mm-1)" << endl;
  cout << "   -P  value: include Poisson noise for the given total number of count (or use -E)" << endl;
  cout << "   -E  value: give a fixed ECF, so that the number of counts is deduced from the raw prompts projection divided by the ECF (or use -P)" << endl;
  cout << "   -M  value: give the mashing power of the output sinogram (default: 1)" << endl;
  cout << "   -S  value: give the span level of the output sinogram (default: 1)" << endl;
  cout << "   -R  value: give the maximum allowed ring difference (default: 1)" << endl;
  cout << "   -L  value: simulate a list-mode acquisition (but also save the prompt sinogram) (have to give a number of counts)" << endl;
  cout << "              the first given value gives the mode for LOR filling (0: do not fill, 1: fill)" << endl;
  cout << "              the second given value is the number of replicates to simulate" << endl;
  #ifdef OMP_MODE
  cout << "   -T  value: give the number of threads for computation (default: 1)" << endl;
  #endif
  cout << "   -D       : use Didier's implementation of the exact Siddon projection algorithm" << endl;
  cout << "   -w       : also save the forward image" << endl;
  cout << "   --noSino : flag to say that we do not want to save any sinogram (when using -L for list-mode)" << endl;
  cout << "   -v  value: give the verbose level" << endl;
  cout << "   -z  value: give the seed of the random generator" << endl;
  cout << "   -h       : print this help page and quit" << endl;
  cout << endl;
  cout << "  This program is used to create a projected sinogram from an image." << endl;
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
  string fileCrystalMap = "";
  string fileImage = "";
  string fileAttenuation = "";
  float psfTransFWHM = -1.;
  float psfAxialFWHM = -1.;
  int span = 1;
  int mash = 1;
  int maxringdiff = 1;
  int nbThreads = 1;
  int verbose = 0;
  float scatterFraction = 0.;
  float randomFraction = 0.;
  float randomFractionLSO = 0.1;
  long int nbCounts = -1;
  float ECF = -1.;
  bool saveImage = false;
  bool listMode = false;
  bool fillEqualLORs = true;
  int seed = -1;
  int nbReplicates = 1;
  float offsetX = 0.;
  float offsetY = 0.;
  float offsetZ = 0.;
  bool saveSino = true;
  int projector = 0;

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
    else if (option=="-i")
    {
      if (i>=argc-1)
      {
        cerr << "***** Not enaugh arguments after option " << option << " !" << endl;
        showHelp(1);
      }
      fileImage = (string)argv[i+1];
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
      i++;
    }
    else if (option=="-s")
    {
      if (i>=argc-1)
      {
        cerr << "***** Not enaugh arguments after option " << option << " !" << endl;
        showHelp(1);
      }
      scatterFraction = atof(argv[i+1]);
      i++;
    }
    else if (option=="-r")
    {
      if (i>=argc-1)
      {
        cerr << "***** Not enaugh arguments after option " << option << " !" << endl;
        showHelp(1);
      }
      randomFraction = atof(argv[i+1]);
      i++;
    }
    else if (option=="-l")
    {
      if (i>=argc-1)
      {
        cerr << "***** Not enaugh arguments after option " << option << " !" << endl;
        showHelp(1);
      }
      randomFractionLSO = atof(argv[i+1]);
      i++;
    }
    else if (option=="-P")
    {
      if (i>=argc-1)
      {
        cerr << "***** Not enaugh arguments after option " << option << " !" << endl;
        showHelp(1);
      }
      nbCounts = (long int)atol(argv[i+1]);
      i++;
    }
    else if (option=="-E")
    {
      if (i>=argc-1)
      {
        cerr << "***** Not enaugh arguments after option " << option << " !" << endl;
        showHelp(1);
      }
      ECF = atof(argv[i+1]);
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
    else if (option=="-z")
    {
      if (i>=argc-1)
      {
        cerr << "***** Not enaugh arguments after option " << option << " !" << endl;
        showHelp(1);
      }
      seed = atoi(argv[i+1]);
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
      nbThreads = atoi(argv[i+1]);
      i++;
    }
    #endif
    else if (option=="-M")
    {
      if (i>=argc-1)
      {
        cerr << "***** Not enaugh arguments after option " << option << " !" << endl;
        showHelp(1);
      }
      mash = atoi(argv[i+1]);
      i++;
    }
    else if (option=="-S")
    {
      if (i>=argc-1)
      {
        cerr << "***** Not enaugh arguments after option " << option << " !" << endl;
        showHelp(1);
      }
      span = atoi(argv[i+1]);
      i++;
    }
    else if (option=="-R")
    {
      if (i>=argc-1)
      {
        cerr << "***** Not enaugh arguments after option " << option << " !" << endl;
        showHelp(1);
      }
      maxringdiff = atoi(argv[i+1]);
      i++;
    }
    else if (option=="-p")
    {
      if (i>=argc-1)
      {
        cerr << "***** Not enaugh arguments after option " << option << " !" << endl;
        showHelp(1);
      }
      string buf = (string)argv[i+1];
      int pos_comma=buf.find_first_of(",");
      if (pos_comma==string::npos)
      {
        psfTransFWHM = atof(buf.c_str());
        psfAxialFWHM = psfTransFWHM;
      }
      else
      {
        psfTransFWHM = atof( buf.substr(0,pos_comma).c_str() );
        psfAxialFWHM = atof( buf.substr(pos_comma+1).c_str() );
      }
      i++;
    }
    else if (option=="-L")
    {
      if (i>=argc-2)
      {
        cerr << "***** Not enaugh arguments after option " << option << " !" << endl;
        showHelp(1);
      }
      fillEqualLORs = (atoi(argv[i+1])!=0);
      nbReplicates  = atoi(argv[i+2]);
      listMode = true;
      i+=2;
    }
    else if (option=="-f")
    {
      if (i>=argc-3)
      {
        cerr << "***** Not enaugh arguments after option " << option << " !" << endl;
        showHelp(1);
      }
      offsetX = atof(argv[i+1]);
      offsetY = atof(argv[i+2]);
      offsetZ = atof(argv[i+3]);
      i+=3;
    }
    else if (option=="-w") saveImage = true;
    else if (option=="--noSino") saveSino = false;
    else if (option=="-D") projector = 1;
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
  if (fileImage=="")
  {
    cerr << "***** Please provide an input image !" << endl;
    showHelp(1);
  }
  if (fileCrystalMap=="")
  {
    cerr << "***** Please provide an esteban crystal map !" << endl;
    showHelp(1);
  }
  if (listMode && nbCounts==-1 && ECF==-1.)
  {
    cerr << "***** Please provide a number of counts or ECF when using list mode !" << endl;
    showHelp(1);
  }
  if (nbReplicates<=0)
  {
    cerr << "***** Please provide an accurate number of replicates !" << endl;
    showHelp(1);
  }

  // ============================================================================================================
  //                                           CREATION OF THE GENERIC LIST MODE
  // ============================================================================================================

  // Create the output manager and log command
  oOutputManager* output = new oOutputManager(fileBaseOut, verbose);
  output->LogCommandLine(argc,argv);

  // Creation of the object
  oSimulator* simulator = new oSimulator(nbThreads, seed, verbose);

  // Initialize all scanner related stuffs
  if (simulator->InitScannerStuff(scannerName, fileCrystalMap, mash, span, maxringdiff))
  {
    LogCerr ("***** Error while initializing all scanner related stuff !" << endl);
    exit(1);
  }

  // Initialize input image
  if (simulator->InitInputImages(fileImage,fileAttenuation,offsetX,offsetY,offsetZ,projector))
  {
    LogCerr ("***** Error while initializing input image !" << endl);
    exit(1);
  }

  // Initialize PSF stuff
  if (simulator->InitPSF(psfTransFWHM, psfAxialFWHM))
  {
    LogCerr ("***** Error while initializing the PSF !" << endl);
    exit(1);
  }

  // Project the image
  if (simulator->Project(saveImage,projector))
  {
    LogCerr ("***** Error while projecting input image for making output sinogram !" << endl);
    exit(1);
  }

  // Apply counts
  if (simulator->ApplyCounts(scatterFraction,randomFraction,randomFractionLSO,nbCounts,ECF,listMode,fillEqualLORs,nbReplicates))
  {
    LogCerr ("***** Error while initializing corrections !" << endl);
    exit(1);
  }

  // Save sinograms
  if (saveSino && simulator->SaveSinograms())
  {
    LogCerr ("***** Error while saving output sinograms !" << endl);
    exit(1);
  }

  // ============================================================================================================
  //                                           THIS IS THE END MY FRIEND !
  // ============================================================================================================

  // Destroy the objects
  delete simulator;
  delete output;

  // Ending
  cout << endl;
  return 0;
}

