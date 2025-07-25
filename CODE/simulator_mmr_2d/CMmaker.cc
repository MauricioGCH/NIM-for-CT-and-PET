#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <ctime>
#include "oScanner.hh"
#include "oOutputManager.hh"
using namespace std;

void showHelp(int code)
{
  cout << endl;
  cout << "Usage:  CMmaker  -m scannerName  {-r || -u}  -o outName" << endl;
  cout << "                 [-w] [-v verbose] [-b head,blockA,blockT]" << endl;
  cout << endl;
  cout << "  -m: scanner names could be 'hrrt', 'hrplus', 'biograph', 'inveon', 'mmr2d','vision600' or 'focus' for the moment." << endl;
  cout << "  -r: random deficiencies of crystal sensitivities, e.g. give 10 means 10% around 1 randomly shot." << endl;
  cout << "  -u: set uniform efficiency" << endl;
  cout << "  -o: gives the output root file name." << endl;
  cout << "  -w: write also the crystal efficiency map alone." << endl;
  cout << "  -b: declare a dead block by giving the head index, the axial block index and the transaxial block index (indices begin at 0)." << endl;
  cout << "  -castor: write the crystal map as a castor LUT file" << endl;
  cout << "  -v: be verbose." << endl;
  cout << "  -h: print this help page." << endl;
  cout << endl;
  cout << "  This programs creates a look-up table containing for each crystal, the 3D cartesian" << endl;
  cout << "  coordinates of 2 of its corners, and its associated crystal efficiency." << endl;
  cout << "  This kind of crystal map will be used with other programs for creating and reconstructing list-mode files." << endl;
  cout << endl;
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

  // File parameters
  string scannerName = "";
  float random_sensitivity_deficiency = 0.;
  bool uniform = false;
  bool random = false;
  string fileBaseOut = "";
  bool castor = false;
  int verbose = 0;
  bool writeCrystEffAlone = false;

  // Dead block parameters
  int nbDeadBlocks = 0;
  int* deadBlockHeadIndex = NULL;
  int* deadBlockAxialIndex = NULL;
  int* deadBlockTransIndex = NULL;

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
    else if (option=="-r")
    {
      if (i>=argc-1)
      {
        cerr << "***** Not enaugh arguments after option " << option << " !" << endl;
        showHelp(1);
      }
      random_sensitivity_deficiency = atof(argv[i+1]);
      random = true;
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
    else if (option=="-b")
    {
      if (i>=argc-1)
      {
        cerr << "***** Not enaugh arguments after option " << option << " !" << endl;
        showHelp(1);
      }
      int head_index, axial_block_index, trans_block_index;
      if (sscanf(argv[i+1],"%d,%d,%d",&head_index,&axial_block_index,&trans_block_index) != 3)
      {
        cerr << "***** Not enaugh arguments after option " << option << " !" << endl;
        showHelp(1);
      }
      i++;
      nbDeadBlocks++;
      deadBlockHeadIndex = (int*)realloc(deadBlockHeadIndex,nbDeadBlocks*sizeof(int));
      deadBlockAxialIndex = (int*)realloc(deadBlockAxialIndex,nbDeadBlocks*sizeof(int));
      deadBlockTransIndex = (int*)realloc(deadBlockTransIndex,nbDeadBlocks*sizeof(int));
      deadBlockHeadIndex[nbDeadBlocks-1] = head_index;
      deadBlockAxialIndex[nbDeadBlocks-1] = axial_block_index;
      deadBlockTransIndex[nbDeadBlocks-1] = trans_block_index;
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
    else if (option=="-w") writeCrystEffAlone = true;
    else if (option=="-u") uniform = true;
    else if (option=="-castor") castor = true;
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
  if (!uniform && !random)
  {
    cerr << "***** Please choose a crystal efficiency uniform or random !" << endl;
    showHelp(1);
  }
  if (fileBaseOut=="")
  {
    cerr << "***** Please provide an output root file name !" << endl;
    showHelp(1);
  }

  // ============================================================================================================
  //                                           CREATION OF THE GENERIC CRYSTAL MAP
  // ============================================================================================================

  // Create the output manager and log command
  oOutputManager* output = new oOutputManager(fileBaseOut, verbose);
  output->LogCommandLine(argc,argv);

  // Creation of the scanner
  oScanner* scanner = new oScanner(scannerName,verbose);

  // Read the crystal efficiency file
  if (scanner->ProcessCrystalEfficiency(uniform, random, random_sensitivity_deficiency))
  {
    LogCerr ("***** A problem occured while processing crystal efficiencies !" << endl);
    exit(1);
  }

  // Process the dead blocks
  if (scanner->ProcessDeadBlocks(nbDeadBlocks,deadBlockHeadIndex,deadBlockAxialIndex,deadBlockTransIndex))
  {
    LogCerr ("***** Failed to process dead blocks in crystal efficiency map !" << endl);
    exit(1);
  }

  // Write the crystal efficiency
  if (writeCrystEffAlone && scanner->WriteCrystalEfficiency(fileBaseOut))
  {
    LogCerr ("***** Failed to write crystal efficiency file alone !" << endl);
    exit(1);
  }

  // Calculate and write the crystal map
  if (scanner->WriteCrystalMap(fileBaseOut, castor))
  {
    LogCerr ("***** Failed to calculate and write the crystal map !" << endl);
    exit(1);
  }

  // ============================================================================================================
  //                                                     END
  // ============================================================================================================

  // Delete objects
  delete scanner;
  delete output;

  // Ending
  cout << endl;
  return 0;
}

