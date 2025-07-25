#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <ctime>
#include <fstream>
#include "oSinoModeCreation.hh"
#include "oTableBiograph.hh"
#include "oOutputManager.hh"
#include "oMiscellaneous.hh"
using namespace std;


// ==========================================================================================================================================
// int InitInputFile()
//    --> This function deals with the input files to switch on true, prompt or list mode.
//        It then calls the appropriate functions.
// ==========================================================================================================================================
int oSinoModeCreation::InitInputFiles( const string& f_FileFloat, const string& f_FileTrue, const string& f_FilePrompt, const string& f_FileListMode,
                                       const string& f_FileRandom, const string& f_FileScatter,
                                       const string& f_FileNormalization, const string& f_FileAttenuation, int f_AttenuationMode,
                                       PRECISION f_TimeDelay, bool f_Mode2D )
{
  if (m_Verbose>=1) LogCout ("oSinoModeCreation::InitInputFile() -> Initialize input file(s)" << endl);  

  // -------------------------------------------------------------------------------------------------------------------
  // Checks for multiple inputs or no input
  if ( f_FileTrue!="" && (f_FilePrompt!="" || f_FileListMode!="" || f_FileFloat!="") ||
       f_FilePrompt!="" && (f_FileTrue!="" || f_FileListMode!="" || f_FileFloat!="") ||
       f_FileListMode!="" && (f_FileTrue!="" || f_FilePrompt!="" || f_FileFloat!="") ||
       f_FileFloat!="" && (f_FileTrue!="" || f_FilePrompt!="" || f_FileListMode!="") )
  {
    LogCerr ("***** oSinoModeCreation::InitInputFiles() -> Please provide only one sort of input file !" << endl);
    return 1;
  }
  if (f_FileTrue=="" && f_FilePrompt=="" && f_FileListMode=="" && f_FileFloat=="")
  {
    LogCerr ("***** oSinoModeCreation::InitInputFiles() -> Please provide an input file !" << endl);
    return 1;
  }

  // 2D mode
  m_Mode2D = f_Mode2D;

  // -------------------------------------------------------------------------------------------------------------------
  // Time delay post-acquisition, if any, in milliseconds
  m_TimeDelay = (int)(f_TimeDelay*1000.);

  // -------------------------------------------------------------------------------------------------------------------
  // Preinit the number of replicates, if any (can be given in list-mode header when it was build on replicates)
  m_NbReplicates = 1;

  // -------------------------------------------------------------------------------------------------------------------
  // Switch on different input possibilities
  if (f_FileTrue!="")
  {
    // Initialize net trues
    if (InitNetTrue(f_FileTrue))
    {
      LogCerr ("***** oSinoModeCreation::InitInputFiles() -> An error occured while initalizing net-trues !" << endl);
      return 1;
    }
  }
  else if (f_FilePrompt!="")
  {
    // Initialize prompts
    if (InitPrompt(f_FilePrompt))
    {
      LogCerr ("***** oSinoModeCreation::InitInputFiles() -> An error occured while initalizing prompts !" << endl);
      return 1;
    }
  }
  else if (f_FileListMode!="")
  {
    // Initialize list-mode
    if (InitListMode(f_FileListMode))
    {
      LogCerr ("***** oSinoModeCreation::InitInputFiles() -> An error occured while initalizing list-mode !" << endl);
      return 1;
    }
  }
  else if (f_FileFloat!="")
  {
    // Check if a scatter or random sinogram is given, throw an error
    if (f_FileScatter!="" || f_FileRandom!="")
    {
      LogCerr ("!!!!! oSinoModeCreation::InitInputFiles() -> Scatter and random corrections are ignored when providing a float sinogram !" << endl);
//      f_FileScatter = "";
//      f_FileRandom  = "";
//      f_FileAttenuation = "";
    }
    // Initialize float sinogram
    if (InitFloat(f_FileFloat))
    {
      LogCerr ("***** oSinoModeCreation::InitInputFiles() -> An error occured while initalizing float sinogram !" << endl);
      return 1;
    }
  }
  m_InitBool = true;

  // -------------------------------------------------------------------------------------------------------------------
  // Initialize random if any
  if (f_FileRandom!="" && f_FileFloat=="")
  {
    if (InitRandom(f_FileRandom))
    {
      LogCerr ("***** oSinoModeCreation::InitInputFiles() -> An error occured while initalizing randoms !" << endl);
      return 1;
    }
  }

  // -------------------------------------------------------------------------------------------------------------------
  // Initialize scatter if any
  if (f_FileScatter!="" && f_FileFloat=="")
  {
    if (InitScatter(f_FileScatter))
    {
      LogCerr ("***** oSinoModeCreation::InitInputFiles() -> An error occured while initalizing scatters !" << endl);
      return 1;
    }
  }

  // -------------------------------------------------------------------------------------------------------------------
  // Initialize normalization if any
  if (f_FileNormalization!="")
  {
    if (InitNormalization(f_FileNormalization))
    {
      LogCerr ("***** oSinoModeCreation::InitInputFiles() -> An error occured while initalizing normalization !" << endl);
      return 1;
    }
  }

  // -------------------------------------------------------------------------------------------------------------------
  // Initialize attenuation if any
  if (f_FileAttenuation!="")
  {
    m_AttenuationMode = f_AttenuationMode;
    if (InitAttenuation(f_FileAttenuation))
    {
      LogCerr ("***** oSinoModeCreation::InitInputFiles() -> An error occured while initalizing attenuation !" << endl);
      return 1;
    }
  }

  // -------------------------------------------------------------------------------------------------------------------
  // Apply scatter fraction
  if (f_FileScatter!="" && f_FileFloat=="")
  {
    if (ApplyScatterFraction())
    {
      LogCerr ("***** oSinoModeCreation::InitInputFiles() -> An error occured while applying scatter fraction !" << endl);
      return 1;
    }
  }

  // -------------------------------------------------------------------------------------------------------------------
  // End
  if (m_Verbose>=1) LogCout ("oSinoModeCreation::InitInputFile() -> Finish initializing input file(s)" << endl);  
  return 0;
}

// ==========================================================================================================================================
// int InitNetTrue()
//    --> This function allocates and reads the net true sinogram.
// ==========================================================================================================================================
int oSinoModeCreation::InitNetTrue(const string& f_FileTrue)
{
  if (m_Verbose>=1) LogCout ("oSinoModeCreation::InitNetTrue() -> Initialize net true sinogram from file '" << f_FileTrue << "'" << endl);

  // -------------------------------------------------------------------------------------------------------------------
  // Read the corresponding file name and open it
  string header_name = f_FileTrue;
  ifstream fhead(header_name.c_str());
  if (!fhead)
  {
    LogCerr ("***** oSinoModeCreation::InitNetTrue() -> Input net true sinogram header '" << header_name << "' is missing or corrupted !" << endl);
    return 1;
  }
  // Verbose
  if (m_Verbose>=2) LogCout ("  --> Reading header file '" << header_name << "' ..." << endl);

  // -------------------------------------------------------------------------------------------------------------------
  // Now we have to find the data file name, the number of prompts and delays, span, maxring, elem, view and sino
  string data_name = "";
  int data_found = 0;
  char line[1024];
  // This is the mandotory number of informations to be found in the header file
  int nb_data_to_find = 12;
  // Special case for inveon where the number of sino is not found, so it will be computed from the span and max ring diff later
  // However, the dead-time correction factor is only present in the header for the inveon scanner, so we need it.
  // It is the same for the focus scanner
//  if (mp_Scanner->GetScannerModel()==SCANNER_INVEON || mp_Scanner->GetScannerModel()==SCANNER_FOCUS) nb_data_to_find = 12;
  if (mp_Scanner->GetScannerModel()==SCANNER_INVEON || mp_Scanner->GetScannerModel()==SCANNER_FOCUS) nb_data_to_find = 11;

  fhead.getline(line,1024);
//  while (!fhead.eof() && data_found<nb_data_to_find)
  while (!fhead.eof())
  {
    size_t found;
    // Transform it to string to benefit from useful c++ functions
    string test = (string)line;
    // Jump to next line if found the # character as the first one
    if ((found=test.find("#"))==0)
    {
      // Read a new line before continuing
      fhead.getline(line,1024);
      continue;
    }
    // ==========================================================================================
    // Here we apply all different tests to find the different infos in all possible header files
    // ==========================================================================================
    // Name of data file for biograph, hrplus and hrrt
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
    // Name of data file for inveon or focus
    found = test.find("file_name");
    if (found!=string::npos)
    {
      if (test.find("acquisition_file_name")==string::npos) // another keyword with file_name that we must ignore
      {
        // Get the file name
        data_name = test.substr(found+9);
        // Remove blancks and cariage return
        while ( (found=data_name.find(" ")) != string::npos) data_name.replace(found,1,"");
        while ( (found=data_name.find("\r")) != string::npos) data_name.replace(found,1,"");
        // Restrict to what is after the last backslash character
        found = data_name.find_last_of("\\");
        data_name = data_name.substr(found+1);
        // Increment the number of data found
        data_found++;
      }
    }
    // Frame number
    found = test.find("frame ");
    if (found!=string::npos)
    {
      // Get the number as a string
      string number = test.substr(found+6);
      // Remove blancks
      while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
      // Get the value
      m_Frame = atoi(number.c_str());
    }
    // Isotope name for biograph, hrplus and hrrt
    if (test.find("ose type")!=string::npos || test.find("sotope name")!=string::npos)
    {
      // Get the file name
      found = test.find("=");
      m_Isotope = test.substr(found+1);
      // Remove blancks and cariage return
      while ( (found=m_Isotope.find(" ")) != string::npos) m_Isotope.replace(found,1,"");
      while ( (found=m_Isotope.find("\r")) != string::npos) m_Isotope.replace(found,1,"");
      // Increment the number of data found
      data_found++;
    }
    // Isotope name for inveon or focus
    found = test.find("sotope ");
    if (mp_Scanner->GetScannerModel()==SCANNER_INVEON && found!=string::npos)
    {
      // Get the file name
      m_Isotope = test.substr(found+7);
      // Remove blancks and cariage return
      while ( (found=m_Isotope.find(" ")) != string::npos) m_Isotope.replace(found,1,"");
      while ( (found=m_Isotope.find("\r")) != string::npos) m_Isotope.replace(found,1,"");
      // Increment the number of data found
      data_found++;
    }
    // Sinogram dimensions for biograph, hrplus and hrrt
    found = test.find("matrix size");
    if (found!=string::npos)
    {
      found = test.find("[1]");
      if (found!=string::npos)
      {
        // Get the number as a string
        found = test.find("=");
        string number = test.substr(found+1);
        // Remove blancks
        while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
        // Get the value
        m_ReadNbElem = atoi(number.c_str());
        // Increment the number of data found
        data_found++;
      }
      found = test.find("[2]");
      if (found!=string::npos)
      {
        // Get the number as a string
        found = test.find("=");
        string number = test.substr(found+1);
        // Remove blancks
        while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
        // Get the value
        m_ReadNbView = atoi(number.c_str());
        // Increment the number of data found
        data_found++;
      }
      found = test.find("[3]");
      if (found!=string::npos)
      {
        // Get the number as a string
        found = test.find("=");
        string number = test.substr(found+1);
        // Remove blancks
        while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
        // Get the value
        m_ReadNbSino = atoi(number.c_str());
        // Increment the number of data found
        data_found++;
      }
    }
    // Sinogram dimensions for inveon or focus: nbElem
    found = test.find("x_dimension");
    if (found!=string::npos)
    {
      // Get the number as a string
      string number = test.substr(found+12);
      // Remove blancks
      while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
      // Get the value
      m_ReadNbElem = atoi(number.c_str());
      // Increment the number of data found
      data_found++;
    }
    // Sinogram dimensions for inveon or focus: nbView
    found = test.find("y_dimension");
    if (found!=string::npos)
    {
      // Get the number as a string
      string number = test.substr(found+12);
      // Remove blancks
      while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
      // Get the value
      m_ReadNbView = atoi(number.c_str());
      // Increment the number of data found
      data_found++;
    }
    // Axial compression for biograph, hrplus and hrrt
    found = test.find("axial compression");
    if (found!=string::npos)
    {
      // Get the number as a string
      found = test.find("=");
      string number = test.substr(found+1);
      // Remove blancks
      while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
      // Get the value
      m_ReadSpan = atoi(number.c_str());
      // Increment the number of data found
      data_found++;
    }
    // Axial compression for inveon or focus
    found = test.find("span");
    if (found!=string::npos)
    {
      // Get the number as a string
      string number = test.substr(found+5);
      // Remove blancks
      while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
      // Get the value
      m_ReadSpan = atoi(number.c_str());
      // Increment the number of data found
      data_found++;
    }
    // Maximum ring difference for biograph, hrplus and hrrt
    found = test.find("maximum ring difference");
    if (found!=string::npos)
    {
      // Get the number as a string
      found = test.find("=");
      string number = test.substr(found+1);
      // Remove blancks
      while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
      // Get the value
      m_ReadMaxRingDiff = atoi(number.c_str());
      // Increment the number of data found
      data_found++;
    }
    // Maximum ring difference for inveon or focus
    found = test.find("ring_difference");
    if (found!=string::npos)
    {
      // Get the number as a string
      string number = test.substr(found+16);
      // Remove blancks
      while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
      // Get the value
      m_ReadMaxRingDiff = atoi(number.c_str());
      // Increment the number of data found
      data_found++;
    }
    // Total number of prompts, randoms and net trues for biograph, hrplus and hrrt
    found = test.find("otal");
    if (found!=string::npos)
    {
      found = test.find("rompts");
      if (found!=string::npos)
      {
        // Get the number as a string
        found = test.find("=");
        string number = test.substr(found+1);
        // Remove blancks
        while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
        // Get the value
        m_NbPrompts = atol(number.c_str());
        // Increment the number of data found
        data_found++;
      }
      found = test.find("andoms");
      if (found!=string::npos)
      {
        // Get the number as a string
        found = test.find("=");
        string number = test.substr(found+1);
        // Remove blancks
        while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
        // Get the value
        m_NbDelays = atol(number.c_str());
        // Increment the number of data found
        data_found++;
      }
      found = test.find("rues");
      if (found!=string::npos)
      {
        // Get the number as a string
        found = test.find("=");
        string number = test.substr(found+1);
        // Remove blancks
        while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
        // Get the value
        m_NbNetTrues = atol(number.c_str());
        // Increment the number of data found
        data_found++;
      }
    }
    // Total number of prompts for inveon or focus (there are 3 successive numbers and we want the first one)
    found = test.find("prompts ");
    if (found!=string::npos)
    {
      // Get the 3 numbers as a string
      string number = test.substr(found+8);
      // Get the first value (atol does the job of sparsing of us !)
      m_NbPrompts = atol(number.c_str());
      // Increment the number of data found
      data_found++;
    }
    // Total number of delays for inveon or focus (there are 3 successive numbers and we want the first one)
    found = test.find("delays ");
    if (found!=string::npos)
    {
      // Get the 3 numbers as a string
      string number = test.substr(found+7);
      // Get the first value (atol does the job of sparsing of us !)
      m_NbDelays = atol(number.c_str());
      // Increment the number of data found
      data_found++;
    }
    // Total number of trues for inveon or focus (there are 3 successive numbers and we want the first one)
    found = test.find("trues ");
    if (found!=string::npos)
    {
      // Get the 3 numbers as a string
      string number = test.substr(found+6);
      // Get the first value (atol does the job of sparsing of us !)
      m_NbNetTrues = atol(number.c_str());
      // Increment the number of data found
      data_found++;
    }
    // Scan duration for biograph, hrplus and hrrt
    found = test.find("image duration");
    if (found!=string::npos)
    {
      // Get the number as a string
      found = test.find("=");
      string number = test.substr(found+1);
      // Remove blancks
      while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
      // Get the value
      m_Duration = (unsigned int)atoi(number.c_str());
      // Increment the number of data found
      data_found++;
    }
    // Scan duration for inveon or focus
    found = test.find("frame_duration ");
    if (found!=string::npos)
    {
      // Get the number as a string
      string number = test.substr(found+15);
      // Remove blancks
      while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
      // Get the value (in float in the file)
      m_Duration = (unsigned int)(atof(number.c_str()));
      // Increment the number of data found
      data_found++;
    }
    // Scan start time for biograph, hrplus and hrrt
    if (test.find("image start time")!=string::npos || test.find("image relative start time")!=string::npos)
    {
      // Get the number as a string
      found = test.find("=");
      string number = test.substr(found+1);
      // Remove blancks
      while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
      // Get the value
      m_StartTime = (unsigned int)atoi(number.c_str());
      // Increment the number of data found
      data_found++;
    }
    // Scan start time for inveon or focus
    found = test.find("frame_start ");
    if (found!=string::npos)
    {
      // Get the number as a string
      string number = test.substr(found+12);
      // Remove blancks
      while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
      // Get the value (in float in the file)
      m_StartTime = (unsigned int)(atof(number.c_str()));
      // Increment the number of data found
      data_found++;
    }
    // Dead time correction factor for inveon or focus
    found = test.find("deadtime_correction ");
    if (found!=string::npos)
    {
      // Get the number as a string
      string number = test.substr(found+20);
      // Remove blancks
      while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
      // Get the value (in float in the file)
      m_DeadTimeCorrectionFactor = atof(number.c_str());
      // Increment the number of data found
      data_found++;
    }
    // Read a new line
    fhead.getline(line,1024);
  }
  // Close header file
  fhead.close();
  // Check if all data were found
  if (data_found<nb_data_to_find)
  {
    LogCerr ("***** oSinoModeCreation::InitNetTrue() -> Failed to get all data (" << nb_data_to_find << ") in header file '" << header_name << "' (actually found: " << data_found << ") !" << endl);
    return 1;
  }

  // Calculate the sino size
  m_ReadSinoSize = m_ReadNbElem * m_ReadNbView;

  // For the inveon and focus scanners, we do not yet have the number of sinograms, so we call the scanner to compute it from the span and maxringdiff
  if (mp_Scanner->GetScannerModel()==SCANNER_INVEON || mp_Scanner->GetScannerModel()==SCANNER_FOCUS)
  {
    m_ReadNbSino = mp_Scanner->ComputeNbSinoFromSpanAndMaxRingDiff(m_ReadSpan,m_ReadMaxRingDiff);
    if (m_ReadNbSino<0)
    {
      LogCerr ("***** oSinoModeCreation::InitNetTrue() -> Inconsistencies found between span and maximum ring difference for inveon or focus scanner !" << endl);
      return 1;
    }
  }

  // -------------------------------------------------------------------------------------------------------------------
  // Allocate the corresponding sinogram
  mp_TrueSino = (short int***)malloc(m_ReadNbSino*sizeof(short int**));
  for (int s=0; s<m_ReadNbSino; s++)
  {
    mp_TrueSino[s] = (short int**)malloc(m_ReadNbView*sizeof(short int*));
    for (int v=0; v<m_ReadNbView; v++) mp_TrueSino[s][v] = (short int*)malloc(m_ReadNbElem*sizeof(short int));
  }

  // -------------------------------------------------------------------------------------------------------------------
  // Add the relative path to the data file name
  int pos; if ((pos=header_name.find_last_of("/"))!=string::npos && header_name.find_first_of("/")!=0) data_name = header_name.substr(0,pos) + "/" + data_name;
  // Open data file
  FILE* fsino = fopen(data_name.c_str(),"rb");
  if (fsino==NULL)
  {
    LogCerr ("***** oSinoModeCreation::InitNetTrue() -> Input net true sinogram '" << data_name << "' is missing or corrupted !" << endl);
    return 1;
  }
  // Verbose
  if (m_Verbose>=2) LogCout ("  --> Reading data file '" << data_name << "' ..." << endl);
  // Read it now
  int nb_data_read = 0;
  int nb_data_to_be_read = nb_data_to_be_read = m_ReadNbSino*m_ReadNbView*m_ReadNbElem;
  for (int s=0; s<m_ReadNbSino; s++) for (int v=0; v<m_ReadNbView; v++)
    nb_data_read += fread(&mp_TrueSino[s][v][0],sizeof(short int),m_ReadNbElem,fsino);
  // Close the file
  fclose(fsino);
  // Check the number of data read
  if (nb_data_read!=nb_data_to_be_read)
  {
    LogCerr ("***** oSinoModeCreation::InitNetTrue() -> Failed to read all data (" << nb_data_to_be_read << ") in file '" << data_name << "' (" << nb_data_read << " read) !" << endl);
    return 1;
  }

  // -------------------------------------------------------------------------------------------------------------------
  // Update start time if delay, and calculate relative timing (time delay was previously converted to ms)
  m_StartTime += m_TimeDelay/1000;
  m_RelativeStartTime = m_StartTime * 1000;
  m_RelativeStopTime  = (m_StartTime + m_Duration) * 1000;

  // -------------------------------------------------------------------------------------------------------------------
  // Compute number of prompts (thus net trues) from sinogram directly
  long int nb_nettrues = 0;
  for (int s=0; s<m_ReadNbSino; s++) for (int v=0; v<m_ReadNbView; v++) for (int e=0; e<m_ReadNbElem; e++)
    nb_nettrues += ((long int)(mp_TrueSino[s][v][e]));
  // Affect
//  m_NbPrompts = nb_nettrues;
  m_NbNetTrues = nb_nettrues;
//  m_NbDelays = 0;

  // Verbose
  if (m_Verbose>=2)
  {
    LogCout ("  --> Radial elements: " << m_ReadNbElem << " | Views: " << m_ReadNbView << " | Sinograms: " << m_ReadNbSino << " | Span: " << m_ReadSpan << " | Maximum ring difference: " << m_ReadMaxRingDiff << endl);
    LogCout ("  --> Number of prompts: " << m_NbPrompts << " | Delays: " << m_NbDelays << " | Net trues: " << m_NbNetTrues << endl);
    LogCout ("  --> Frame: " << m_Frame << " | Start acquisition: " << m_StartTime << " s | Duration: " << m_Duration << " s" << endl);
  }

  // -------------------------------------------------------------------------------------------------------------------
  // Ending
  m_TrueBool = true;
  return 0;
}

// ==========================================================================================================================================
// int InitPrompt()
//    --> This function allocates and reads the prompt sinogram.
// ==========================================================================================================================================
int oSinoModeCreation::InitPrompt(const string& f_FilePrompt)
{
  if (m_Verbose>=1) LogCout ("oSinoModeCreation::InitPrompt() -> Initialize prompt sinogram from file '" << f_FilePrompt << "'" << endl);

  // -------------------------------------------------------------------------------------------------------------------
  // Read the corresponding file name and open it
  string header_name = f_FilePrompt;
  ifstream fhead(header_name.c_str());
  if (!fhead)
  {
    LogCerr ("***** oSinoModeCreation::InitPrompt() -> Input prompt sinogram header '" << header_name << "' is missing or corrupted !" << endl);
    return 1;
  }
  // Verbose
  if (m_Verbose>=2) LogCout ("  --> Reading header file '" << header_name << "' ..." << endl);

  // -------------------------------------------------------------------------------------------------------------------
  // Now we have to find the data file name, the number of prompts and delays, span, maxring, elem, view and sino
  int nb_data_to_find = 10;
  string data_name = "";
  int data_found = 0;
  char line[1024];
  fhead.getline(line,1024);
  while (!fhead.eof() && data_found<nb_data_to_find)
  {
    size_t found;
    // Transform it to string to benefit from useful c++ functions
    string test = (string)line;
    // Apply tests
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
    if (test.find("ose type")!=string::npos || test.find("sotope name")!=string::npos)
    {
      // Get the file name
      found = test.find("=");
      m_Isotope = test.substr(found+1);
      // Remove blancks and cariage return
      while ( (found=m_Isotope.find(" ")) != string::npos) m_Isotope.replace(found,1,"");
      while ( (found=m_Isotope.find("\r")) != string::npos) m_Isotope.replace(found,1,"");
      // Do not increment the number of data found (not mandatory information)
      // data_found++;
    }
    found = test.find("matrix size");
    if (found!=string::npos)
    {
      found = test.find("[1]");
      if (found!=string::npos)
      {
        // Get the number as a string
        found = test.find("=");
        string number = test.substr(found+1);
        // Remove blancks
        while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
        // Get the value
        m_ReadNbElem = atoi(number.c_str());
        // Increment the number of data found
        data_found++;
      }
      found = test.find("[2]");
      if (found!=string::npos)
      {
        // Get the number as a string
        found = test.find("=");
        string number = test.substr(found+1);
        // Remove blancks
        while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
        // Get the value
        m_ReadNbView = atoi(number.c_str());
        // Increment the number of data found
        data_found++;
      }
      found = test.find("[3]");
      if (found!=string::npos)
      {
        // Get the number as a string
        found = test.find("=");
        string number = test.substr(found+1);
        // Remove blancks
        while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
        // Get the value
        m_ReadNbSino = atoi(number.c_str());
        // Increment the number of data found
        data_found++;
      }
    }
    found = test.find("axial compression");
    if (found!=string::npos)
    {
      // Get the number as a string
      found = test.find("=");
      string number = test.substr(found+1);
      // Remove blancks
      while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
      // Get the value
      m_ReadSpan = atoi(number.c_str());
      // Increment the number of data found
      data_found++;
    }
    found = test.find("maximum ring difference");
    if (found!=string::npos)
    {
      // Get the number as a string
      found = test.find("=");
      string number = test.substr(found+1);
      // Remove blancks
      while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
      // Get the value
      m_ReadMaxRingDiff = atoi(number.c_str());
      // Increment the number of data found
      data_found++;
    }
    // Total number of prompts, randoms and net trues
    found = test.find("otal");
    if (found!=string::npos)
    {
      found = test.find("rompts");
      if (found!=string::npos)
      {
        // Get the number as a string
        found = test.find("=");
        string number = test.substr(found+1);
        // Remove blancks
        while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
        // Get the value
        m_NbPrompts = atol(number.c_str());
        // Increment the number of data found
        data_found++;
      }
      found = test.find("andoms");
      if (found!=string::npos)
      {
        // Get the number as a string
        found = test.find("=");
        string number = test.substr(found+1);
        // Remove blancks
        while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
        // Get the value
        m_NbDelays = atol(number.c_str());
        // Increment the number of data found
        data_found++;
      }
    }
    found = test.find("image duration");
    if (found!=string::npos)
    {
      // Get the number as a string
      found = test.find("=");
      string number = test.substr(found+1);
      // Remove blancks
      while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
      // Get the value
      m_Duration = (unsigned int)atoi(number.c_str());
      // Increment the number of data found
      data_found++;
    }
    if (test.find("image start time")!=string::npos || test.find("image relative start time")!=string::npos)
    {
      // Get the number as a string
      found = test.find("=");
      string number = test.substr(found+1);
      // Remove blancks
      while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
      // Get the value
      m_StartTime = (unsigned int)atoi(number.c_str());
      // Increment the number of data found
      data_found++;
    }
    // Read a new line
    fhead.getline(line,1024);
  }
  // Close header file
  fhead.close();
  // Check if all data were found
  if (data_found!=nb_data_to_find)
  {
    LogCerr ("***** oSinoModeCreation::InitPrompt() -> Failed to get all data (" << nb_data_to_find << ") in header file '" << header_name << "' (actually found: " << data_found << ") !" << endl);
    return 1;
  }

  // Warning if no isotope found
  if (m_Isotope=="INF")
  {
    LogCerr ("!!!!! oSinoModeCreation::InitPrompt() -> No real isotope definition found, will not correct for decay nor branching ratio !" << endl);
  }

  // Calculate net trues
  m_NbNetTrues = m_NbPrompts - m_NbDelays;

  // Calculate the sino size
  m_ReadSinoSize = m_ReadNbElem * m_ReadNbView;

  // -------------------------------------------------------------------------------------------------------------------
  // Allocate the corresponding sinogram
  mp_PromptSino = (short int***)malloc(m_ReadNbSino*sizeof(short int**));
  for (int s=0; s<m_ReadNbSino; s++)
  {
    mp_PromptSino[s] = (short int**)malloc(m_ReadNbView*sizeof(short int*));
    for (int v=0; v<m_ReadNbView; v++) mp_PromptSino[s][v] = (short int*)malloc(m_ReadNbElem*sizeof(short int));
  }

  // -------------------------------------------------------------------------------------------------------------------
  // Add the relative path to the data file name
  int pos; if ((pos=header_name.find_last_of("/"))!=string::npos && header_name.find_first_of("/")!=0) data_name = header_name.substr(0,pos) + "/" + data_name;
  // Open data file
  FILE* fsino = fopen(data_name.c_str(),"rb");
  if (fsino==NULL)
  {
    LogCerr ("***** oSinoModeCreation::InitPrompt() -> Input prompt sinogram '" << data_name << "' is missing or corrupted !" << endl);
    return 1;
  }
  // Verbose
  if (m_Verbose>=2) LogCout ("  --> Reading data file '" << data_name << "' ..." << endl);
  // Read it now
  int nb_data_read = 0;
  int nb_data_to_be_read = nb_data_to_be_read = m_ReadNbSino*m_ReadNbView*m_ReadNbElem;
  for (int s=0; s<m_ReadNbSino; s++) for (int v=0; v<m_ReadNbView; v++)
    nb_data_read += fread(&mp_PromptSino[s][v][0],sizeof(short int),m_ReadNbElem,fsino);
  // Close the file
  fclose(fsino);
  // Check the number of data read
  if (nb_data_read!=nb_data_to_be_read)
  {
    LogCerr ("***** oSinoModeCreation::InitPrompt() -> Failed to read all data (" << nb_data_to_be_read << ") in file '" << data_name << "' (" << nb_data_read << " read) !" << endl);
    return 1;
  }

  // -------------------------------------------------------------------------------------------------------------------
  // Update start time if delay, and calculate relative timing
  m_StartTime += m_TimeDelay/1000;
  m_RelativeStartTime = m_StartTime * 1000;
  m_RelativeStopTime  = (m_StartTime + m_Duration) * 1000;

  // Verbose
  if (m_Verbose>=2)
  {
    LogCout ("  --> Radial elements: " << m_ReadNbElem << " | Views: " << m_ReadNbView << " | Sinograms: " << m_ReadNbSino << " | Span: " << m_ReadSpan << " | Maximum ring difference: " << m_ReadMaxRingDiff << endl);
    LogCout ("  --> Number of prompts: " << m_NbPrompts << " | Delays: " << m_NbDelays << " | Net trues: " << m_NbNetTrues << endl);
    LogCout ("  --> Start acquisition: " << m_StartTime << " s | Duration: " << m_Duration << " s" << endl);
  }

  // -------------------------------------------------------------------------------------------------------------------
  // Ending
  m_PromptBool = true;
  return 0;
}

// ==========================================================================================================================================
// int InitFloat()
//    --> This function allocates and reads the float sinogram.
// ==========================================================================================================================================
int oSinoModeCreation::InitFloat(const string& f_FileFloat)
{
  if (m_Verbose>=1) LogCout ("oSinoModeCreation::InitFloat() -> Initialize float sinogram from file '" << f_FileFloat << "'" << endl);

  // -------------------------------------------------------------------------------------------------------------------
  // Read the corresponding file name and open it
  string header_name = f_FileFloat;
  ifstream fhead(header_name.c_str());
  if (!fhead)
  {
    LogCerr ("***** oSinoModeCreation::InitFloat() -> Input float sinogram header '" << header_name << "' is missing or corrupted !" << endl);
    return 1;
  }
  // Verbose
  if (m_Verbose>=2) LogCout ("  --> Reading header file '" << header_name << "' ..." << endl);

  // Assume duration of 1s
  m_Duration  = 1.;
  m_StartTime = 0.;

  // -------------------------------------------------------------------------------------------------------------------
  // Now we have to find the data file name, span, maxring, elem, view and sino
  string data_name = "";
  int data_found = 0;
  char line[1024];
  fhead.getline(line,1024);
  int nb_data_to_find = 6;
  // We do not have the number of sinograms with inveon or focus
  if (mp_Scanner->GetScannerModel()==SCANNER_INVEON || mp_Scanner->GetScannerModel()==SCANNER_FOCUS) nb_data_to_find = 5; 
//  while (!fhead.eof() && data_found<nb_data_to_find)
  while (!fhead.eof())
  {
    size_t found;
    // Transform it to string to benefit from useful c++ functions
    string test = (string)line;
    // Jump to next line if found the # character as the first one
    if ((found=test.find("#"))==0)
    {
      // Read a new line before continuing
      fhead.getline(line,1024);
      continue;
    }
    // Name of data file for biograph and hrplus
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
    // Name of data file for inveon or focus
    found = test.find("file_name");
    if (found!=string::npos)
    {
      if (test.find("acquisition_file_name")==string::npos) // another keyword with file_name that we must ignore
      {
        // Get the file name
        data_name = test.substr(found+9);
        // Remove blancks and cariage return
        while ( (found=data_name.find(" ")) != string::npos) data_name.replace(found,1,"");
        while ( (found=data_name.find("\r")) != string::npos) data_name.replace(found,1,"");
        // Restrict to what is after the last backslash character
        found = data_name.find_last_of("\\");
        data_name = data_name.substr(found+1);
        // Increment the number of data found
        data_found++;
      }
    }
    // Sinogram dimensions for biograph and hrplus
    found = test.find("matrix size");
    if (found!=string::npos)
    {
      found = test.find("[1]");
      if (found!=string::npos)
      {
        // Get the number as a string
        found = test.find("=");
        string number = test.substr(found+1);
        // Remove blancks
        while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
        // Get the value
        m_ReadNbElem = atoi(number.c_str());
        // Increment the number of data found
        data_found++;
      }
      found = test.find("[2]");
      if (found!=string::npos)
      {
        // Get the number as a string
        found = test.find("=");
        string number = test.substr(found+1);
        // Remove blancks
        while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
        // Get the value
        m_ReadNbView = atoi(number.c_str());
        // Increment the number of data found
        data_found++;
      }
      found = test.find("[3]");
      if (found!=string::npos)
      {
        // Get the number as a string
        found = test.find("=");
        string number = test.substr(found+1);
        // Remove blancks
        while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
        // Get the value
        m_ReadNbSino = atoi(number.c_str());
        // Increment the number of data found
        data_found++;
      }
    }
    // Sinogram dimensions for inveon or focus: nbElem
    found = test.find("x_dimension");
    if (found!=string::npos)
    {
      // Get the number as a string
      string number = test.substr(found+12);
      // Remove blancks
      while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
      // Get the value
      m_ReadNbElem = atoi(number.c_str());
      // Increment the number of data found
      data_found++;
    }
    // Sinogram dimensions for inveon or focus: nbView
    found = test.find("y_dimension");
    if (found!=string::npos)
    {
      // Get the number as a string
      string number = test.substr(found+12);
      // Remove blancks
      while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
      // Get the value
      m_ReadNbView = atoi(number.c_str());
      // Increment the number of data found
      data_found++;
    }
    // Axial compression for biograph and hrplus
    found = test.find("axial compression");
    if (found!=string::npos)
    {
      // Get the number as a string
      found = test.find("=");
      string number = test.substr(found+1);
      // Remove blancks
      while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
      // Get the value
      m_ReadSpan = atoi(number.c_str());
      // Increment the number of data found
      data_found++;
    }
    // Axial compression for inveon or focus
    found = test.find("span"); // Problem if span is in the name of another line ...
    if (found!=string::npos)
    {
      // Test if span if before the sign =
      size_t tmp = test.find("=");
      if (tmp>found)
      {
        // Get the number as a string
        string number = test.substr(found+5);
        // Remove blancks
        while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
        // Get the value
        m_ReadSpan = atoi(number.c_str());
        // Increment the number of data found
        data_found++;
      }
    }
    // Maximum ring difference for biograph and hrplus
    found = test.find("maximum ring difference");
    if (found!=string::npos)
    {
      // Get the number as a string
      found = test.find("=");
      string number = test.substr(found+1);
      // Remove blancks
      while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
      // Get the value
      m_ReadMaxRingDiff = atoi(number.c_str());
      // Increment the number of data found
      data_found++;
    }
    // Maximum ring difference for inveon or focus
    found = test.find("ring_difference");
    if (found!=string::npos)
    {
      // Get the number as a string
      string number = test.substr(found+16);
      // Remove blancks
      while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
      // Get the value
      m_ReadMaxRingDiff = atoi(number.c_str());
      // Increment the number of data found
      data_found++;
    }
    // Read a new line
    fhead.getline(line,1024);
  }
  // Close header file
  fhead.close();
  // Check if all data were found
  if (data_found<nb_data_to_find)
  {
    LogCerr ("***** oSinoModeCreation::InitFloat() -> Failed to get all data (" << nb_data_to_find << ") in header file '" << header_name << "' (actually found: " << data_found << ") !" << endl);
    return 1;
  }

  // Overload some data if 2D mode
  if (m_Mode2D)
  {
    m_ReadMaxRingDiff = 1;
    m_ReadSpan = 3;
    m_ReadNbSino = -1;
  }

  // Cannot know the number of net trues nor the others
  m_NbNetTrues = m_NbPrompts = m_NbDelays = 0;

  // Calculate the sino size
  m_ReadSinoSize = m_ReadNbElem * m_ReadNbView;

  // For the inveon and focus scanners, we do not yet have the number of sinograms, so we call the scanner to compute it from the span and maxringdiff
  if (mp_Scanner->GetScannerModel()==SCANNER_INVEON || mp_Scanner->GetScannerModel()==SCANNER_FOCUS)
  {
    m_ReadNbSino = mp_Scanner->ComputeNbSinoFromSpanAndMaxRingDiff(m_ReadSpan,m_ReadMaxRingDiff);
    if (m_ReadNbSino<0)
    {
      LogCerr ("***** oSinoModeCreation::InitNetTrue() -> Inconsistencies found between span and maximum ring difference for inveon or focus scanner !" << endl);
      return 1;
    }
  }

  // -------------------------------------------------------------------------------------------------------------------
  // Allocate the corresponding sinogram
  mp_FloatSino = (float***)malloc(m_ReadNbSino*sizeof(float**));
  for (int s=0; s<m_ReadNbSino; s++)
  {
    mp_FloatSino[s] = (float**)malloc(m_ReadNbView*sizeof(float*));
    for (int v=0; v<m_ReadNbView; v++) mp_FloatSino[s][v] = (float*)malloc(m_ReadNbElem*sizeof(float));
  }

  // -------------------------------------------------------------------------------------------------------------------
  // Add the relative path to the data file name
  int pos; if ((pos=header_name.find_last_of("/"))!=string::npos && header_name.find_first_of("/")!=0) data_name = header_name.substr(0,pos) + "/" + data_name;
  // Open data file
  FILE* fsino = fopen(data_name.c_str(),"rb");
  if (fsino==NULL)
  {
    LogCerr ("***** oSinoModeCreation::InitFloat() -> Input float sinogram '" << data_name << "' is missing or corrupted !" << endl);
    return 1;
  }
  // Verbose
  if (m_Verbose>=2) LogCout ("  --> Reading data file '" << data_name << "' ..." << endl);
  // Read it now
  int nb_data_read = 0;
  int nb_data_to_be_read = nb_data_to_be_read = m_ReadNbSino*m_ReadNbView*m_ReadNbElem;
  for (int s=0; s<m_ReadNbSino; s++) for (int v=0; v<m_ReadNbView; v++)
    nb_data_read += fread(&mp_FloatSino[s][v][0],sizeof(float),m_ReadNbElem,fsino);
  // Close the file
  fclose(fsino);
  // Check the number of data read
  if (nb_data_read!=nb_data_to_be_read)
  {
    LogCerr ("***** oSinoModeCreation::InitFloat() -> Failed to read all data (" << nb_data_to_be_read << ") in file '" << data_name << "' (" << nb_data_read << " read) !" << endl);
    return 1;
  }

  // -------------------------------------------------------------------------------------------------------------------
  // Update start time if delay, and calculate relative timing
  m_StartTime += m_TimeDelay/1000;
  m_RelativeStartTime = m_StartTime * 1000;
  m_RelativeStopTime  = (m_StartTime + m_Duration) * 1000;

  // Verbose
  if (m_Verbose>=2)
  {
    LogCout ("  --> Radial elements: " << m_ReadNbElem << " | Views: " << m_ReadNbView << " | Sinograms: " << m_ReadNbSino << " | Span: " << m_ReadSpan << " | Maximum ring difference: " << m_ReadMaxRingDiff << endl);
    LogCout ("  --> Number of prompts: " << m_NbPrompts << " | Delays: " << m_NbDelays << " | Net trues: " << m_NbNetTrues << endl);
    LogCout ("  --> Start acquisition: " << m_StartTime << " s | Duration: " << m_Duration << " s" << endl);
  }

  // -------------------------------------------------------------------------------------------------------------------
  // Ending
  m_FloatBool = true;
  return 0;
}

// ==========================================================================================================================================
// int InitListMode()
//    --> This function read the list-mode to create a prompt sinogram.
// ==========================================================================================================================================
int oSinoModeCreation::InitListMode(const string& f_FileListMode)
{
  if (m_Verbose>=1) LogCout ("oSinoModeCreation::InitListMode() -> Initialize prompt sinogram from list-mode file '" << f_FileListMode << "'" << endl);

  // ======================================================================================================================
  // Check scanner model (only valid for biograph for the moment)
  if (mp_Scanner->GetScannerModel()!=SCANNER_BIOGRAPH)
  {
    LogCerr ("***** oSinoModeCreation::InitListMode() -> Giving a list-mode as input is only valid for Biograph !" << endl);
    return 1;
  }

  // ======================================================================================================================
  // Fix sinogram dimensions
  m_ReadNbElem = mp_Scanner->GetNbElem();
  m_ReadNbView = mp_Scanner->GetNbView();
  m_ReadNbSino = 559;
  m_ReadSpan = 11;
  m_ReadMaxRingDiff = 38;
  m_ReadSinoSize = m_ReadNbElem * m_ReadNbView;

  // ======================================================================================================================
  // Open and read list-mode input header file (duration, start-time, isotope, data-name)
  ifstream fhead(f_FileListMode.c_str());
  if (!fhead)
  {
    LogCerr ("***** oListModeCreation::ProcessListMode() -> Input header file '" << f_FileListMode << "' is missing or corrupted !" << endl);
    return 1;
  }
  // Verbose
  if (m_Verbose>=1) LogCout ("  --> Read input list-mode header file" << endl);
  // Read
  int nb_data_to_find = 4;
  string data_name = "";
  int data_found = 0;
  char line[1024];
  fhead.getline(line,1024);
  while (!fhead.eof())
  {
    size_t found;
    // Transform it to string to benefit from useful c++ functions
    string test = (string)line;
    // Test 1: data file name
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
    // test 2: isotope
    if (test.find("ose type")!=string::npos || test.find("sotope name")!=string::npos)
    {
      // Get the file name
      found = test.find("=");
      m_Isotope = test.substr(found+1);
      // Remove blancks and cariage return
      while ( (found=m_Isotope.find(" ")) != string::npos) m_Isotope.replace(found,1,"");
      while ( (found=m_Isotope.find("\r")) != string::npos) m_Isotope.replace(found,1,"");
      // Increment the number of data found
      data_found++;
    }
    found = test.find("image duration");
    if (found!=string::npos)
    {
      // Get the number as a string
      found = test.find("=");
      string number = test.substr(found+1);
      // Remove blancks
      while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
      // Get the value
      m_Duration = (unsigned int)atoi(number.c_str());
      // Increment the number of data found
      data_found++;
    }
    if (test.find("relative time of tracer injection")!=string::npos)
    {
      // Get the number as a string
      found = test.find("=");
      string number = test.substr(found+1);
      // Remove blancks
      while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
      // Get the value
      m_StartTime = (unsigned int)atoi(number.c_str());
      // Increment the number of data found
      data_found++;
    }
    found = test.find("replicates");
    if (found!=string::npos)
    {
      // Get the number as a string
      found = test.find("=");
      string number = test.substr(found+1);
      // Remove blancks
      while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
      // Get the value
      m_NbReplicates = atoi(number.c_str());
      // Do not increment the number of data found (not mandatory)
    }
    // Read a new line
    fhead.getline(line,1024);
  }
  // Close
  fhead.close();
  // Check if all data were found
  if (data_found!=nb_data_to_find)
  {
    LogCerr ("***** oListModeCreation::ProcessListMode() -> Failed to get all data in header file '" << f_FileListMode << "' !" << endl);
    return 1;
  }

  // -------------------------------------------------------------------------------------------------------------------
  // Allocate the corresponding sinogram
  mp_PromptSino = (short int***)malloc(m_ReadNbSino*sizeof(short int**));
  for (int s=0; s<m_ReadNbSino; s++)
  {
    mp_PromptSino[s] = (short int**)malloc(m_ReadNbView*sizeof(short int*));
    for (int v=0; v<m_ReadNbView; v++) mp_PromptSino[s][v] = (short int*)calloc(m_ReadNbElem,sizeof(short int));
  }

  // -------------------------------------------------------------------------------------------------------------------
  // Add the relative path to the data file name and open it
  int pos; if ((pos=f_FileListMode.find_last_of("/"))!=string::npos && f_FileListMode.find_first_of("/")!=0) data_name = f_FileListMode.substr(0,pos) + "/" + data_name;
  // Open data file
  FILE* flm = fopen(data_name.c_str(),"rb");
  if (flm==NULL)
  {
    LogCerr ("***** oSinoModeCreation::InitPrompt() -> Input prompt sinogram '" << data_name << "' is missing or corrupted !" << endl);
    return 1;
  }
  // Verbose
  if (m_Verbose>=2) LogCout ("  --> Reading data file '" << data_name << "' ..." << endl);

  // -------------------------------------------------------------------------------------------------------------------
  // Read list-mode (compute prompts and delays)
  int word_data = 0;
  int word_code = 0;
  while (1)
  {
    // 1. Read a word
    word_code = ReadWord32bits(flm,&word_data);
    // 2. Break if EOF
    if (word_code==WORD_EOF) break;
    // 3. Analyze the word if this is an event
    if (word_code==WORD_EVENT)
    {
      // Read the event type (1 is prompt, 0 is delay)
      unsigned short int event_type = (unsigned short int) (((word_data)>>30) & 1);
      if (event_type==0)
      {
        m_NbDelays++;
        continue;
      }
      else m_NbPrompts++;
      // Read the bin address
      int bin_adress = (word_data) & 268435455; // Set the first 4 bits to 0
      // Determine sinogram coordinates
      int bin_sino = bin_adress / m_ReadSinoSize;
      int bin_view = (bin_adress % m_ReadSinoSize) / m_ReadNbElem;
      int bin_elem = bin_adress % m_ReadNbElem;
      // Increment sinogram
      mp_PromptSino[bin_sino][bin_view][bin_elem]++;
    }
  }
  // Close file
  fclose(flm);
  // Calculate net trues
  m_NbNetTrues = m_NbPrompts - m_NbDelays;

  // -------------------------------------------------------------------------------------------------------------------
  // Update start time if delay, and calculate relative timing
  m_StartTime += m_TimeDelay/1000;
  m_RelativeStartTime = m_StartTime * 1000;
  m_RelativeStopTime  = (m_StartTime + m_Duration) * 1000;

  // Verbose
  if (m_Verbose>=2)
  {
    LogCout ("  --> Radial elements: " << m_ReadNbElem << " | Views: " << m_ReadNbView << " | Sinograms: " << m_ReadNbSino << " | Span: " << m_ReadSpan << " | Maximum ring difference: " << m_ReadMaxRingDiff << endl);
    LogCout ("  --> Number of prompts: " << m_NbPrompts << " | Delays: " << m_NbDelays << " | Net trues: " << m_NbNetTrues << endl);
    LogCout ("  --> Start acquisition: " << m_StartTime << " s | Duration: " << m_Duration << " s" << endl);
  }

  // -------------------------------------------------------------------------------------------------------------------
  // Ending
  m_PromptBool = true;
  return 0;
}
// ==========================================================================================================================================
// int ReadWord32bits()
//    --> This function reads a 32bits list-mode word.
// ==========================================================================================================================================
int oSinoModeCreation::ReadWord32bits(FILE* fin, int* f_Word)
{
  if (m_Verbose>=6) LogCout ("oListModeCreation::ReadWord32bits() -> Read a word" << endl);

  // Test if we are at the EOF (should never happen)
  if (feof(fin)) return WORD_EOF;

  // Otherwise first dumbly read the word
  int nb_data_read = fread(f_Word,sizeof(int),1,fin);

  // Test if we reach the end of file
  if (nb_data_read!=1) return WORD_EOF;

  // Then we begin to analyze the word type
  int test = ((*f_Word)>>31) & 1;
  if (test==0) return WORD_EVENT;
  else return WORD_OTHER;
}


