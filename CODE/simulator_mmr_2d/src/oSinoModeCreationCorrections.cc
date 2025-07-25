#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <ctime>
#include <fstream>
#include <vector>
#include <algorithm>
#include "oSinoModeCreation.hh"
#include "oTableBiograph.hh"
#include "oTableVision600.hh"
#include "oOutputManager.hh"
#include "oMiscellaneous.hh"
using namespace std;

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//                                                         S C A T T E R
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// ==========================================================================================================================================
// int InitScatter()
//    --> This function allocates and reads the scatter sinogram.
// ==========================================================================================================================================
int oSinoModeCreation::InitScatter(string f_FileScat)
{
  if (m_Verbose>=1) LogCout ("oSinoModeCreation::InitScatter() -> Initialize scatter sinogram from file '" << f_FileScat << "'" << endl);

  // Security: if input file has not been read, stop
  if (!m_InitBool)
  {
    LogCerr ("***** oSinoModeCreation::InitScatter() -> Cannot initialize scatter before input counts !" << endl);
    return 1;
  }

  // First allocate the corresponding sinogram
  if (mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH)
  {
    mp_ScatSino = (float***)malloc(mp_Scanner->GetNbPlanes()*sizeof(float**));
    for (int s=0; s<mp_Scanner->GetNbPlanes(); s++)
    {
      mp_ScatSino[s] = (float**)malloc(m_ReadNbView*sizeof(float*));
      for (int v=0; v<m_ReadNbView; v++) mp_ScatSino[s][v] = (float*)malloc(m_ReadNbElem*sizeof(float));
    }
  }
  else if (mp_Scanner->GetScannerModel()==SCANNER_CASTOR)
  {
    mp_ScatSino = (float***)malloc(mp_Scanner->GetNbPlanes()*sizeof(float**));
    for (int s=0; s<mp_Scanner->GetNbPlanes(); s++)
    {
      mp_ScatSino[s] = (float**)malloc(m_ReadNbView*sizeof(float*));
      for (int v=0; v<m_ReadNbView; v++) mp_ScatSino[s][v] = (float*)malloc(m_ReadNbElem*sizeof(float));
    }
  }
  else if (mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH2D || mp_Scanner->GetScannerModel()==SCANNER_MMR2D)
  {
    mp_ScatSino = (float***)malloc(mp_Scanner->GetNbPlanes()*sizeof(float**));
    for (int s=0; s<mp_Scanner->GetNbPlanes(); s++)
    {
      mp_ScatSino[s] = (float**)malloc(m_ReadNbView*sizeof(float*));
      for (int v=0; v<m_ReadNbView; v++) mp_ScatSino[s][v] = (float*)malloc(m_ReadNbElem*sizeof(float));
    }
  }
  else if (mp_Scanner->GetScannerModel()==SCANNER_HRPLUS || mp_Scanner->GetScannerModel()==SCANNER_INVEON || mp_Scanner->GetScannerModel()==SCANNER_FOCUS)
  {
    mp_ScatSino = (float***)malloc(m_ReadNbSino*sizeof(float**));
    for (int s=0; s<m_ReadNbSino; s++)
    {
      mp_ScatSino[s] = (float**)malloc(m_ReadNbView*sizeof(float*));
      for (int v=0; v<m_ReadNbView; v++) mp_ScatSino[s][v] = (float*)malloc(m_ReadNbElem*sizeof(float));
    }
  }
  else if (mp_Scanner->GetScannerModel()==SCANNER_HRRT)
  {
    mp_ScatSino = (float***)malloc(m_ReadNbSino*sizeof(float**));
    for (int s=0; s<m_ReadNbSino; s++)
    {
      mp_ScatSino[s] = (float**)malloc(m_ReadNbView*sizeof(float*));
      for (int v=0; v<m_ReadNbView; v++) mp_ScatSino[s][v] = (float*)malloc(m_ReadNbElem*sizeof(float));
    }
  }
  else if (mp_Scanner->GetScannerModel()==SCANNER_VISION600)
  {
    mp_ScatSino = (float***)malloc(m_ReadNbSino*sizeof(float**));
    for (int s=0; s<m_ReadNbSino; s++)
    {
      mp_ScatSino[s] = (float**)malloc(m_ReadNbView*sizeof(float*));
      for (int v=0; v<m_ReadNbView; v++) mp_ScatSino[s][v] = (float*)malloc(m_ReadNbElem*sizeof(float));
    }
  }
  else
  {
    LogCerr ("***** oSinoModeCreation::InitScatter() -> Unknown scanner model !" << endl);
    return 1;
  }

  // Special case for the FOCUS scanner for which we can provide the log file of the e7_sino program to get the actual scatter fraction
  // The two files are separated by a comma
  string ecat_sino_log_file = "";
  if (mp_Scanner->GetScannerModel()==SCANNER_FOCUS)
  {
    // Search for a comma
    size_t comma = f_FileScat.find(",");
    if (comma!=string::npos)
    {
      ecat_sino_log_file = f_FileScat.substr(comma+1);
      f_FileScat = f_FileScat.substr(0,comma);
      if (m_Verbose>=1) LogCout("  --> Ecat_sino log file provided as '" << ecat_sino_log_file << "'" << endl);
    }
    // Get the scatter fraction
    if (ecat_sino_log_file=="")
    {
      m_ScatterFraction = -1.;
    }
    else
    {
      // Open file
      ifstream flog(ecat_sino_log_file.c_str());
      if (!flog)
      {
        LogCerr ("***** oSinoModeCreation::InitScatter() -> Input ecat_sino log file '" << ecat_sino_log_file << "' is missing or corrupted !" << endl);
        return 1;
      }
      // Search inside for the last line indicating the scatter fraction
      char line[1024];
      flog.getline(line,1024);
      bool not_found = true;
      while (!flog.eof() && not_found)
      {
        string test = (string)line;
        size_t found = test.find("scatter fraction");
        // We found it
        if (found!=string::npos)
        {
          // Get the number
          found = test.find("=");
          string number = test.substr(found+1);
          // Remove blancks
          while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
          // Remove percentage sign
          while ( (found=number.find("%")) != string::npos) number.replace(found,1,"");
          // Convert the number
          m_ScatterFraction = ((double)(atof(number.c_str())/100.));
          // Say we found it to exit the while loop
          not_found = false;
        }
        // Get a new line
        flog.getline(line,1024);
      }
      // Close file
      flog.close();
      // If found it, say it
      if (!not_found && m_Verbose>=1) LogCout("  --> Scatter fraction: " << m_ScatterFraction << endl);
    }
  }

  // ???????????????????????????????????????????????????????????????????????????????????????????????????????????????
  // ???????????????????????????????????????????????????????????????????????????????????????????????????????????????
  //   Should also read dimensions in the header and check if everything is consistent with the input sinogram.
  // ???????????????????????????????????????????????????????????????????????????????????????????????????????????????
  // ???????????????????????????????????????????????????????????????????????????????????????????????????????????????

  // Then read the corresponding file name and open it
  string header_name = f_FileScat;
  ifstream fhead(header_name.c_str());
  if (!fhead)
  {
    LogCerr ("***** oSinoModeCreation::InitScatter() -> Input scatter sinogram header '" << header_name << "' is missing or corrupted !" << endl);
    return 1;
  }
  // Verbose
  if (m_Verbose>=2) LogCout ("  --> Reading header file '" << header_name << "' ..." << endl);
  // Now we have to find the data file name
  string data_name = "";
  int data_found = 0;
  char line[1024];
  fhead.getline(line,1024);
  while (!fhead.eof() && data_found<1)
  {
    size_t found;
    // Transform it to string to benefit from useful c++ functions
    string test = (string)line;
    // Data file name for biograph, hrplus and hrrt
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
    // Data file name for inveon or focus
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
    // Read a new line
    fhead.getline(line,1024);
  }
  // Close header file
  fhead.close();
  // Check if all data were found
  if (data_found!=1)
  {
    LogCerr ("***** oSinoModeCreation::InitScatter() -> Failed to get all data in header file '" << header_name << "' !" << endl);
    return 1;
  }

  // Add the relative path to the data file name
  int pos; if ((pos=header_name.find_last_of("/"))!=string::npos && header_name.find_first_of("/")!=0) data_name = header_name.substr(0,pos) + "/" + data_name;
  // Open data file
  FILE* fsino = fopen(data_name.c_str(),"rb");
  if (fsino==NULL)
  {
    LogCerr ("***** oSinoModeCreation::InitScatter() -> Input scatter sinogram '" << data_name << "' is missing or corrupted !" << endl);
    return 1;
  }
  // Verbose
  if (m_Verbose>=2) LogCout ("  --> Reading data file '" << data_name << "' ..." << endl);
  // Read it now
  int nb_data_read = 0;
  int nb_data_to_be_read = 0;
  if (mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH)
  {
    nb_data_to_be_read = mp_Scanner->GetNbPlanes()*m_ReadSinoSize;
    for (int s=0; s<mp_Scanner->GetNbPlanes(); s++) for (int v=0; v<m_ReadNbView; v++)
      nb_data_read += fread(&mp_ScatSino[s][v][0],sizeof(float),m_ReadNbElem,fsino);
  }
  else if (mp_Scanner->GetScannerModel()==SCANNER_CASTOR)
  {
    nb_data_to_be_read = mp_Scanner->GetNbPlanes()*m_ReadSinoSize;
    for (int s=0; s<mp_Scanner->GetNbPlanes(); s++) for (int v=0; v<m_ReadNbView; v++)
      nb_data_read += fread(&mp_ScatSino[s][v][0],sizeof(float),m_ReadNbElem,fsino);
  }
  else if (mp_Scanner->GetScannerModel()==SCANNER_VISION600)
  {
    nb_data_to_be_read = mp_Scanner->GetNbPlanes()*m_ReadSinoSize;
    for (int s=0; s<mp_Scanner->GetNbPlanes(); s++) for (int v=0; v<m_ReadNbView; v++)
      nb_data_read += fread(&mp_ScatSino[s][v][0],sizeof(float),m_ReadNbElem,fsino);
  }
  else if (mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH2D || mp_Scanner->GetScannerModel()==SCANNER_MMR2D)
  {
    nb_data_to_be_read = mp_Scanner->GetNbPlanes()*m_ReadSinoSize;
    for (int s=0; s<mp_Scanner->GetNbPlanes(); s++) for (int v=0; v<m_ReadNbView; v++)
      nb_data_read += fread(&mp_ScatSino[s][v][0],sizeof(float),m_ReadNbElem,fsino);
  }
  else if (mp_Scanner->GetScannerModel()==SCANNER_HRPLUS || mp_Scanner->GetScannerModel()==SCANNER_INVEON || mp_Scanner->GetScannerModel()==SCANNER_FOCUS)
  {
    nb_data_to_be_read = m_ReadNbSino*m_ReadSinoSize;
    for (int s=0; s<m_ReadNbSino; s++) for (int v=0; v<m_ReadNbView; v++)
      nb_data_read += fread(&mp_ScatSino[s][v][0],sizeof(float),m_ReadNbElem,fsino);
  }
  else if (mp_Scanner->GetScannerModel()==SCANNER_HRRT)
  {
    nb_data_to_be_read = m_ReadNbSino*m_ReadSinoSize;
    for (int s=0; s<m_ReadNbSino; s++) for (int v=0; v<m_ReadNbView; v++)
      nb_data_read += fread(&mp_ScatSino[s][v][0],sizeof(float),m_ReadNbElem,fsino);
  }
  else
  {
    LogCerr ("***** oSinoModeCreation::InitScatter() -> Unknown scanner model !" << endl);
    return 1;
  }
  // Close the file
  fclose(fsino);
  // Check the number of data read
  if (nb_data_read!=nb_data_to_be_read)
  {
    LogCerr ("***** oSinoModeCreation::InitScatter() -> Failed to read all data (" << nb_data_to_be_read << ") in file '" << data_name << "' (" << nb_data_read << " read) !" << endl);
    return 1;
  }

  // Set the scatter correction to true
  m_ScatCorr = true;

  // Ending
  return 0;
}

// ==========================================================================================================================================
// int ApplyScatterFraction()
//    --> This function uses the scatter fraction to scale the scatter sinogram.
// ==========================================================================================================================================
int oSinoModeCreation::ApplyScatterFraction()
{
  // If negative scatter fraction, just quit
  if (m_ScatterFraction<0.) return 0;

  // Verbose
  if (m_Verbose>=1) LogCout("oSinoModeCreation::ApplyScatterFraction() -> Apply scatter fraction " << m_ScatterFraction << endl);

  // Compute net trues: two cases
  double total_net_trues = 0.;
  if (m_TrueBool)
  {
    long int total = 0;
    for (int s=0; s<m_ReadNbSino; s++) for (int v=0; v<m_ReadNbView; v++) for (int e=0; e<m_ReadNbElem; e++)
      total += ((long int)(mp_TrueSino[s][v][e]));
    total_net_trues = ((double)(total));
  }
  else if (m_PromptBool && m_RandCorr)
  {
    long int total_prompt = 0;
    for (int s=0; s<m_ReadNbSino; s++) for (int v=0; v<m_ReadNbView; v++) for (int e=0; e<m_ReadNbElem; e++)
      total_prompt += ((long int)(mp_PromptSino[s][v][e]));
    total_net_trues = ((double)(total_prompt));
    for (int s=0; s<m_ReadNbSino; s++) for (int v=0; v<m_ReadNbView; v++) for (int e=0; e<m_ReadNbElem; e++)
      total_net_trues -= mp_RandSino[s][v][e];
  }
  else
  {
    return 0;
  }

  // Verbose
  if (m_Verbose>=1) LogCout("  --> Total net-trues: " << total_net_trues << endl);

  // Compute current total scatter, un-normalized
  double total_scatter = 0.;
  if (m_NormCorr)
  {
    for (int s=0; s<m_ReadNbSino; s++) for (int v=0; v<m_ReadNbView; v++) for (int e=0; e<m_ReadNbElem; e++)
      if (mp_NormSino[s][v][e]>0.) total_scatter += ((double)(mp_ScatSino[s][v][e])) / ((double)(mp_NormSino[s][v][e]));
  }
  else
  {
    for (int s=0; s<m_ReadNbSino; s++) for (int v=0; v<m_ReadNbView; v++) for (int e=0; e<m_ReadNbElem; e++)
      total_scatter += mp_ScatSino[s][v][e];
  }

  // Verbose
  if (m_Verbose>=1) LogCout("  --> Total un-normalized scatters: " << total_scatter << endl);

  // Compute scale factor
  double scale_factor = total_net_trues * m_ScatterFraction / total_scatter;

  // Verbose
  if (m_Verbose>=1) LogCout("  --> Scale factor: " << scale_factor << endl);

  // Apply it
  for (int s=0; s<m_ReadNbSino; s++) for (int v=0; v<m_ReadNbView; v++) for (int e=0; e<m_ReadNbElem; e++)
    mp_ScatSino[s][v][e] *= scale_factor;

  // Normal end
  return 0;
}

// ==========================================================================================================================================
// int ComputeScatterRate()
//    --> This function computes the scatter rate based on the input scatter sinogram.
// ==========================================================================================================================================
float oSinoModeCreation::ComputeScatterRate(int f_Elem, int f_View, int f_Sino)
{
  // Correction or not
  if (!m_ScatCorr) return 0.;

  // Verbose
  if (m_Verbose>=5) LogCout ("oSinoModeCreation::ComputeScatterRate() -> Compute rate for bin [" << f_Elem << ";" << f_View << ";" << f_Sino << "]" << endl);

  // The result
  float scatter = 0.;

  // Switch on scanner model
  if (mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH)
  {
    // Get the first ring pair for this sinogram  (the sum of rings is always the same)
    unsigned int ring1 = 0, ring2 = 0;
    mp_BiographTable->GetRingPairBySinoIndex(f_Sino, 0, &ring1, &ring2);
    // Get the scatter rate
    scatter = mp_ScatSino[ring1+ring2][f_View][f_Elem] / ((float)m_Duration);
  }
  else if (mp_Scanner->GetScannerModel()==SCANNER_CASTOR)
  {
    // Get the first ring pair for this sinogram  (the sum of rings is always the same)
    unsigned int ring1 = 0, ring2 = 0;
    mp_CastorTable->GetRingPairBySinoIndex(f_Sino, 0, &ring1, &ring2);
    // Get the scatter rate
    scatter = mp_ScatSino[ring1+ring2][f_View][f_Elem] / ((float)m_Duration);
  }
  else if (mp_Scanner->GetScannerModel()==SCANNER_VISION600)
  {
    // Get the first ring pair for this sinogram  (the sum of rings is always the same)
    unsigned int ring1 = 0, ring2 = 0;
    mp_Vision600Table->GetRingPairBySinoIndex(f_Sino, 0, &ring1, &ring2);
    // Get the scatter rate
    scatter = mp_ScatSino[ring1+ring2][f_View][f_Elem] / ((float)m_Duration);
  }
  else if (mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH2D || mp_Scanner->GetScannerModel()==SCANNER_MMR2D)
  {
    // Get the scatter rate
    scatter = mp_ScatSino[0][f_View][f_Elem] / ((float)m_Duration);
  }
  else if (mp_Scanner->GetScannerModel()==SCANNER_HRPLUS || mp_Scanner->GetScannerModel()==SCANNER_INVEON || mp_Scanner->GetScannerModel()==SCANNER_FOCUS)
  {
    scatter = mp_ScatSino[f_Sino][f_View][f_Elem] / ((float)m_Duration);
  }

  // Return rate
  return scatter;
}

// ==========================================================================================================================================
// int ComputeTotalScatter()
//    --> This function computes the total amount of un-normalized scatters.
// ==========================================================================================================================================
int oSinoModeCreation::ComputeTotalScatter()
{
  if (m_Verbose>=1) LogCout ("oSinoModeCreation::ComputeTotalScatter() -> Compute total scatter counts" << endl);

  // Check
  if (!m_ScatCorr || !m_NormCorr)
  {
    LogCerr ("***** oSinoModeCreation::ComputeTotalScatter() -> Cannot compute total scatter before loading scatters and normalization !" << endl);
    return 1;
  }
  
  // Initialize
  m_NbScatters = 0.;

  // Compute
  for (int s=0; s<m_ReadNbSino; s++)
  {
    // The sinogram index
    int sino_index = -1;
    // Switch on scanner model
    if (mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH)
    {
      // Get the first ring pair for this sinogram  (the sum of rings is always the same)
      unsigned int ring1 = 0, ring2 = 0;
      mp_BiographTable->GetRingPairBySinoIndex(s, 0, &ring1, &ring2);
      sino_index = ring1+ring2;
    }
    else if (mp_Scanner->GetScannerModel()==SCANNER_CASTOR)
    {
      // Get the first ring pair for this sinogram  (the sum of rings is always the same)
      unsigned int ring1 = 0, ring2 = 0;
      mp_CastorTable->GetRingPairBySinoIndex(s, 0, &ring1, &ring2);
      sino_index = ring1+ring2;
    }
    else if (mp_Scanner->GetScannerModel()==SCANNER_VISION600)
    {

      // Get the first ring pair for this sinogram  (the sum of rings is always the same)
      unsigned int ring1 = 0, ring2 = 0;
      //mp_Vision600Table->GetRingPairBySinoIndex(s, 0, &ring1, &ring2);
      
      //sino_index = ring1+ring2;
      sino_index = 0;
    }
    else if (mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH2D || mp_Scanner->GetScannerModel()==SCANNER_MMR2D)
    {
      sino_index = 0;
    }
    else if (mp_Scanner->GetScannerModel()==SCANNER_HRPLUS || mp_Scanner->GetScannerModel()==SCANNER_INVEON || mp_Scanner->GetScannerModel()==SCANNER_FOCUS)
    {
      sino_index = s;
    }
    else if (mp_Scanner->GetScannerModel()==SCANNER_HRRT)
    {
      sino_index = s;
    }
    else
    {
      LogCerr ("***** oSinoModeCreation::ComputeTotalScatter() -> Unknown scanner model !" << endl);
      return 1;
    }

    for (int v=0; v<m_ReadNbView; v++) for (int e=0; e<m_ReadNbElem; e++)
    {
      
      if (mp_NormSino[s][v][e]>0.)
        m_NbScatters += ((double)(mp_ScatSino[sino_index][v][e])) / ((double)(mp_NormSino[s][v][e]));;
    }
  }

  // Verbose
  if (m_Verbose>=1)
  {
    LogCout ("  --> Total un-normalized scatter counts: " << m_NbScatters << endl);
    if (m_NbRandoms==0.)
    {
      LogCout ("  --> Random proportion: " << ((double)m_NbDelays)*100./((double)m_NbPrompts) << " %" << endl);
      LogCout ("  --> Scatter proportion: " << ((double)m_NbScatters)*100./((double)m_NbPrompts) << " %" << endl);
      LogCout ("  --> Scatter fraction: " << ((double)m_NbScatters)*100./(((double)m_NbPrompts)-((double)m_NbDelays)) << " %" << endl);
      LogCout ("  --> Trues / Prompts: " << (((double)m_NbPrompts)-((double)m_NbDelays)-m_NbScatters)*100./((double)m_NbPrompts) << " %" << endl);
    }
    else
    {
      LogCout ("  --> Random proportion: " << ((double)m_NbRandoms)*100./((double)m_NbPrompts) << " %" << endl);
      LogCout ("  --> Scatter proportion: " << ((double)m_NbScatters)*100./((double)m_NbPrompts) << " %" << endl);
      LogCout ("  --> Scatter fraction: " << ((double)m_NbScatters)*100./(((double)m_NbPrompts)-((double)m_NbRandoms)) << " %" << endl);
      LogCout ("  --> Trues / Prompts: " << (((double)m_NbPrompts)-((double)m_NbRandoms)-m_NbScatters)*100./((double)m_NbPrompts) << " %" << endl);
    }
  }

  // Ending
  return 0;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//                                                         R A N D O M
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// ==========================================================================================================================================
// int InitRandom()
//    --> This function allocates and reads the random sinogram.
// ==========================================================================================================================================
int oSinoModeCreation::InitRandom(const string& f_FileRand)
{
  if (m_Verbose>=1) LogCout ("oSinoModeCreation::InitRandom() -> Initialize random sinogram from file '" << f_FileRand << "'" << endl);

  // Security: if input file has not been read, stop
  if (!m_InitBool)
  {
    LogCerr ("***** oSinoModeCreation::InitRandom() -> Cannot initialize random before input counts !" << endl);
    return 1;
  }

  // First allocate the corresponding sinogram
  mp_RandSino = (float***)malloc(m_ReadNbSino*sizeof(float**));
  for (int s=0; s<m_ReadNbSino; s++)
  {
    mp_RandSino[s] = (float**)malloc(m_ReadNbView*sizeof(float*));
    for (int v=0; v<m_ReadNbView; v++) mp_RandSino[s][v] = (float*)malloc(m_ReadNbElem*sizeof(float));
  }

  // ???????????????????????????????????????????????????????????????????????????????????????????????????????????????
  // ???????????????????????????????????????????????????????????????????????????????????????????????????????????????
  //   Should also read dimensions in the header and check if everything is consistent with the input sinogram.
  // ???????????????????????????????????????????????????????????????????????????????????????????????????????????????
  // ???????????????????????????????????????????????????????????????????????????????????????????????????????????????

  // Then read the corresponding file name and open it
  string header_name = f_FileRand;
  ifstream fhead(header_name.c_str());
  if (!fhead)
  {
    LogCerr ("***** oSinoModeCreation::InitRandom() -> Input random sinogram header '" << header_name << "' is missing or corrupted !" << endl);
    return 1;
  }
  // Verbose
  if (m_Verbose>=2) LogCout ("  --> Reading header file '" << header_name << "' ..." << endl);
  // Now we have to find the data file name
  string data_name = "";
  int data_found = 0;
  char line[1024];
  fhead.getline(line,1024);
  while (!fhead.eof() && data_found<1)
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
    // Read a new line
    fhead.getline(line,1024);
  }
  // Close header file
  fhead.close();
  // Check if all data were found
  if (data_found!=1)
  {
    LogCerr ("***** oSinoModeCreation::InitRandom() -> Failed to get all data in header file '" << header_name << "' !" << endl);
    return 1;
  }

  // Add the relative path to the data file name
  int pos; if ((pos=header_name.find_last_of("/"))!=string::npos && header_name.find_first_of("/")!=0) data_name = header_name.substr(0,pos) + "/" + data_name;
  // Open data file
  FILE* fsino = fopen(data_name.c_str(),"rb");
  if (fsino==NULL)
  {
    LogCerr ("***** oSinoModeCreation::InitRandom() -> Input random sinogram '" << data_name << "' is missing or corrupted !" << endl);
    return 1;
  }
  // Verbose
  if (m_Verbose>=2) LogCout ("  --> Reading data file '" << data_name << "' ..." << endl);
  // Read it now
  int nb_data_read = 0;
  int nb_data_to_be_read = 0;
  nb_data_to_be_read = m_ReadNbSino*m_ReadSinoSize;
  for (int s=0; s<m_ReadNbSino; s++) for (int v=0; v<m_ReadNbView; v++)
    nb_data_read += fread(&mp_RandSino[s][v][0],sizeof(float),m_ReadNbElem,fsino);
  // Close the file
  fclose(fsino);
  // Check the number of data read
  if (nb_data_read!=nb_data_to_be_read)
  {
    LogCerr ("***** oSinoModeCreation::InitRandom() -> Failed to read all data (" << nb_data_to_be_read << ") in file '" << data_name << "' (" << nb_data_read << " read) !" << endl);
    return 1;
  }

  // Set the random correction to true
  m_RandCorr = true;

  // Ending
  return 0;
}

// ==========================================================================================================================================
// int ComputeRandomRate()
//    --> This function computes the random rate based on the input random sinogram.
// ==========================================================================================================================================
float oSinoModeCreation::ComputeRandomRate(int f_Elem, int f_View, int f_Sino)
{
  // Correction or not
  if (!m_RandCorr) return 0.;

  // Verbose
  if (m_Verbose>=5) LogCout ("oSinoModeCreation::ComputeRandomRate() -> Compute rate for bin [" << f_Elem << ";" << f_View << ";" << f_Sino << "]" << endl);

  // The result
  float random = mp_RandSino[f_Sino][f_View][f_Elem] / ((float)m_Duration);

  // Return rate
  return random;
}

// ==========================================================================================================================================
// int ComputeTotalRandom()
//    --> This function computes the total amount of un-normalized random.
// ==========================================================================================================================================
int oSinoModeCreation::ComputeTotalRandom()
{
  if (m_Verbose>=1) LogCout ("oSinoModeCreation::ComputeTotalRandom() -> Compute total random counts" << endl);

  // Check
  if (!m_RandCorr)
  {
    LogCerr ("***** oSinoModeCreation::ComputeTotalRandom() -> Cannot compute total random before loading randoms !" << endl);
    return 1;
  }

  // Initialize
  m_NbRandoms = 0.;

  // Compute
  for (int s=0; s<m_ReadNbSino; s++)
  {
    for (int v=0; v<m_ReadNbView; v++) for (int e=0; e<m_ReadNbElem; e++)
    {
      m_NbRandoms += mp_RandSino[s][v][e];
    }
  }

  // Verbose
  if (m_Verbose>=1)
  {
    LogCout ("  --> Total un-normalized counts: " << m_NbRandoms << endl);
  }

  // Ending
  return 0;
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//                                                     N O R M A L I Z A T I O N
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// ==========================================================================================================================================
// int InitNormalization()
//    --> This function allocates the normalization sinogram and read it from the given file.
// ==========================================================================================================================================
int oSinoModeCreation::InitNormalization(const string& f_FileNorm)
{
  if (m_Verbose>=1) LogCout ("oSinoModeCreation::InitNormalization() -> Initialize normalization sinogram from file '" << f_FileNorm << "'" << endl);

  // Security: if net trues have not been read, stop
  if (!m_InitBool)
  {
    LogCerr ("***** oSinoModeCreation::InitNormalization() -> Cannot initialize normalization before net trues !" << endl);
    return 1;
  }

  // First allocate the corresponding sinogram
  mp_NormSino = (float***)malloc(m_ReadNbSino*sizeof(float**));
  for (int s=0; s<m_ReadNbSino; s++)
  {
    mp_NormSino[s] = (float**)malloc(m_ReadNbView*sizeof(float*));
    for (int v=0; v<m_ReadNbView; v++) mp_NormSino[s][v] = (float*)malloc(m_ReadNbElem*sizeof(float));
  }

  // ???????????????????????????????????????????????????????????????????????????????????????????????????????????????
  // ???????????????????????????????????????????????????????????????????????????????????????????????????????????????
  //   Should also read dimensions in the header and check if everything is consistent with the input sinogram.
  // ???????????????????????????????????????????????????????????????????????????????????????????????????????????????
  // ???????????????????????????????????????????????????????????????????????????????????????????????????????????????

  // Then read the corresponding file name and open it
  string header_name = f_FileNorm;
  ifstream fhead(header_name.c_str());
  if (!fhead)
  {
    LogCerr ("***** oSinoModeCreation::InitNormalization() -> Input normalization sinogram header '" << header_name << "' is missing or corrupted !" << endl);
    return 1;
  }
  // Verbose
  if (m_Verbose>=2) LogCout ("  --> Reading header file '" << header_name << "' ..." << endl);
  // Now we have to find the data file name
  string data_name = "";
  int data_found = 0;
  char line[1024];
  // Number of data info to find
  int nb_data_to_find = 1;
  if (mp_Scanner->GetScannerModel()==SCANNER_INVEON || mp_Scanner->GetScannerModel()==SCANNER_FOCUS) nb_data_to_find = 2; // Have to find the ECF too
  fhead.getline(line,1024);
//  while (!fhead.eof() && data_found<1)
  while (!fhead.eof())
  {
    size_t found;
    // Transform it to string to benefit from useful c++ functions
    string test = (string)line;
    // Data file name for biograph, hrplus and hrrt
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
    // Data file name for inveon or focus
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
    // ECF for inveon or focus
    found = test.find("calibration_factor ");
    if (found!=string::npos)
    {
      // Get the number as a string
      string number = test.substr(found+19);
      // Remove blancks
      while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
      // Get the value (in float in the file)
      m_ECF = (double)atof(number.c_str());
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
    LogCerr ("***** oSinoModeCreation::InitNormalization() -> Failed to get all data in header file '" << header_name << "' !" << endl);
    return 1;
  }

  // Add the relative path to the data file name
  int pos; if ((pos=header_name.find_last_of("/"))!=string::npos && header_name.find_first_of("/")!=0) data_name = header_name.substr(0,pos) + "/" + data_name;

  // Then read the corresponding file name and open it
  FILE* fsino = fopen(data_name.c_str(),"rb");
  if (fsino==NULL)
  {
    LogCerr ("***** oSinoModeCreation::InitNormalization() -> Input normalization sinogram '" << data_name << "' is missing or corrupted !" << endl);
    return 1;
  }
  if (m_Verbose>=2) LogCout ("  --> Reading file '" << data_name << "' ..." << endl);

  // Read it now
  int nb_data_read = 0;
  int nb_data_to_be_read = m_ReadNbSino*m_ReadSinoSize;
  for (int s=0; s<m_ReadNbSino; s++) for (int v=0; v<m_ReadNbView; v++)
    nb_data_read += fread(&mp_NormSino[s][v][0],sizeof(float),m_ReadNbElem,fsino);

  // Close the file
  fclose(fsino);

  // Check the number of data read
  if (nb_data_read!=nb_data_to_be_read)
  {
    LogCerr ("***** oSinoModeCreation::InitNormalization() -> Failed to read all data (" << nb_data_to_be_read << ") in file '" << data_name << "' (" << nb_data_read << " read) !" << endl);
    return 1;
  }

  // Set the normalization correction to true
  m_NormCorr = true;

  // Get the ECF for the HRPLUS (same for Biograph)
  if (mp_Scanner->GetScannerModel()==SCANNER_HRPLUS || mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH ||
      mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH2D || mp_Scanner->GetScannerModel()==SCANNER_CASTOR || mp_Scanner->GetScannerModel()==SCANNER_MMR2D)
  {
    // Get the calibration file name and open it
    string calib_name = data_name.substr(0,data_name.find_last_of("."))+".calib.hdr";
    ifstream fcal(calib_name.c_str());
    if (!fcal)
    {
      LogCerr ("***** oSinoModeCreation::InitNormalization() -> Input normalization calibration header '" << calib_name << "' is missing or corrupted !" << endl);
      return 1;
    }
    // Verbose
    if (m_Verbose>=2) LogCout ("  --> Reading calibration header file '" << calib_name << "' ..." << endl);
    // Now we have to find the ecf
    int data_found = 0;
    char line[1024];
    fcal.getline(line,1024);
    while (!fcal.eof() && data_found<1)
    {
      size_t found;
      // Transform it to string to benefit from useful c++ functions
      string test = (string)line;
      // Test 1: for the data file name
      if (test.find("calibration factor")!=string::npos || test.find("scanner quantification factor")!=string::npos)
      {
        // Get the ECF as a string
        found = test.find("=");
        string number = test.substr(found+1);
        // Remove blancks and cariage return
        while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
        // Get the ECF
        m_ECF = (double)atof(number.c_str());
        // Increment the number of data found
        data_found++;
      }
      // Read a new line
      fcal.getline(line,1024);
    }
    // Close header file
    fcal.close();
    // Check if all data were found
    if (data_found!=1)
    {
      LogCerr ("***** oSinoModeCreation::InitNormalization() -> Failed to get all data in calibration header file '" << calib_name << "' !" << endl);
      return 1;
    }
    // Verbose
    if (m_Verbose>=2) LogCout ("      Calibration factor: " << m_ECF << endl);
  }
  else if (mp_Scanner->GetScannerModel()==SCANNER_HRRT)
  {
    LogCerr ("!!!!! oSinoModeCreation::InitNormalization() -> No calibration yet !" << endl);
  }
  else if (mp_Scanner->GetScannerModel()==SCANNER_INVEON || mp_Scanner->GetScannerModel()==SCANNER_FOCUS)
  {
    // The ECF was already found into the standard header
    if (m_Verbose>=2) LogCout ("      Calibration factor: " << m_ECF << endl);
  }
  else if (mp_Scanner->GetScannerModel()==SCANNER_VISION600)
  {
    // Include sinogram bin size into the ECF
    m_ECF *=1;
    // Include LOR-DOI factor into the ECF
    m_ECF *=1;
    // Must introduce a factor of 4 here to be consistent with quantification (comes from the scanner crystal efficiencies)
    m_ECF /= 1.;
    if (m_Verbose>=2) LogCout ("      Calibration factor: " << m_ECF << " (tabke radial bin size, LOR-DOI factor and a correction on crystal efficiencies into account)" << endl);
  }
  else
  {
    LogCerr ("***** oSinoModeCreation::InitNormalization() -> Unknown scanner model !" << endl);
    return 1;
  }

  // For inveon and focus scanners, the dead-time correction factor is global and found in the header of the trues sinogram,
  // so we apply it now on the normalization sinogram
  if (mp_Scanner->GetScannerModel()==SCANNER_INVEON || mp_Scanner->GetScannerModel()==SCANNER_FOCUS)
  {
    // Verbose
    if (m_Verbose>=2) LogCout ("  --> Applying dead-time correction factor (" << m_DeadTimeCorrectionFactor << ") to normalization sinogram ..." << endl);
    for (int s=0; s<m_ReadNbSino; s++) for (int v=0; v<m_ReadNbView; v++) for (int e=0; e<m_ReadNbElem; e++)
      mp_NormSino[s][v][e] *= m_DeadTimeCorrectionFactor;
  }

  // Ending
  return 0;
}

// ==========================================================================================================================================
// int ComputeNormalizationFactor()
//    --> This function computes the normalization from the sinogram.
// ==========================================================================================================================================
float oSinoModeCreation::ComputeNormalizationFactor(int f_Elem, int f_View, int f_Sino)
{
  // Correction or not
  if (!m_NormCorr) return 1.;

  // Verbose
  if (m_Verbose>=5) LogCout ("oSinoModeCreation::ComputeNormalizationFactor() -> Compute factor for bin [" << f_Elem << ";" << f_View << ";" << f_Sino << "]" << endl);

  // Get the normalization factor
  float norm_factor = mp_NormSino[f_Sino][f_View][f_Elem];

  // Return rate
  return norm_factor;
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//                                                     A T T E N U A T I O N
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// ==========================================================================================================================================
// int InitAttenuation()
//    --> This function allocates the u-map and read it from the given file.
// ==========================================================================================================================================
int oSinoModeCreation::InitAttenuation(const string& f_FileAttn)
{
  // Particular case if attenuation sinogram is provided
  if (m_AttenuationMode==ATTENUATION_FROM_SINO)
  {
    if (m_Verbose>=1) LogCout ("oSinoModeCreation::InitAttenuation() -> Initialize attenuation sinogram from file '" << f_FileAttn << "'" << endl);
    // Allocate sinogram
    mp_AttnSino = (float***)malloc(m_ReadNbSino*sizeof(float**));
    for (int s=0; s<m_ReadNbSino; s++)
    {
      mp_AttnSino[s] = (float**)malloc(m_ReadNbView*sizeof(float*));
      for (int v=0; v<m_ReadNbView; v++) mp_AttnSino[s][v] = (float*)malloc(m_ReadNbElem*sizeof(float));
    }
    // Open file
    FILE* fdata = fopen(f_FileAttn.c_str(),"rb");
    if (fdata==NULL)
    {
      LogCerr ("***** oSinoModeCreation::InitAttenuation() -> Input sinogram data ACF file '" << f_FileAttn << "' is missing or corrupted !" << endl);
      return 1;
    }
    // Read data
    int nb_data_read = 0;
    for (int s=0; s<m_ReadNbSino; s++) for (int v=0; v<m_ReadNbView; v++)
      nb_data_read += fread(&mp_AttnSino[s][v][0],sizeof(float),m_ReadNbElem,fdata);
    // Close data file
    fclose(fdata);
    // Check
    if (nb_data_read!=m_ReadNbSino*m_ReadNbView*m_ReadNbElem)
    {
      LogCerr ("***** oSinoModeCreation::InitAttenuation() -> Failed to read all data (" << m_ReadNbSino*m_ReadNbView*m_ReadNbElem << ") in ACF sinogram (" << nb_data_read << " read) !" << endl);
      return 1;
    }
    // Set the attenuation correction to true
    m_AttnCorr = true;
    // And directly exit the function
    return 0;
  }

  // ==========  From there, we assume that an attenuation map has been provided and not a sinogram
  if (m_Verbose>=1) LogCout ("oSinoModeCreation::InitAttenuation() -> Initialize u-map from file '" << f_FileAttn << "'" << endl);

  // Open header file
  ifstream fhead(f_FileAttn.c_str());
  if (!fhead)
  {
    LogCerr ("***** oSinoModeCreation::InitAttenuation() -> Input u-map header '" << f_FileAttn << "' is missing or corrupted !" << endl);
    return 1;
  }

  // Now we have to find the data file name, the dimensions and voxel sizes
  int nb_data_to_find = 7;
  string data_name = "";
  int data_found = 0;
  char line[1024];
  fhead.getline(line,1024);

  while (!fhead.eof()) // Read the whole file
  {
    size_t found;
    // Transform it to string to benefit from useful c++ functions
    string test = (string)line;
    // Data file name for biograph, hrplus and hrrt
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
    // Data file name for inveon or focus
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
    // Image dimensions for biograph, hrplus and hrrt
    found = test.find("matrix size");
    if (found!=string::npos)
    {
      // Test 2: for the dimX
      found = test.find("[1]");
      if (found!=string::npos)
      {
        // Get the number as a string
        found = test.find("=");
        string number = test.substr(found+1);
        // Remove blancks
        while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
        // Get the value
        m_AttnDimX = atoi(number.c_str());
        // Increment the number of data found
        data_found++;
      }
      // Test 3: for the dimY
      found = test.find("[2]");
      if (found!=string::npos)
      {
        // Get the number as a string
        found = test.find("=");
        string number = test.substr(found+1);
        // Remove blancks
        while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
        // Get the value
        m_AttnDimY = atoi(number.c_str());
        // Increment the number of data found
        data_found++;
      }
      // Test 4: for the dimZ
      found = test.find("[3]");
      if (found!=string::npos)
      {
        // Get the number as a string
        found = test.find("=");
        string number = test.substr(found+1);
        // Remove blancks
        while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
        // Get the value
        m_AttnDimZ = atoi(number.c_str());
        // Increment the number of data found
        data_found++;
      }
    }
    // Image dimensions for inveon or focus: x
    found = test.find("x_dimension");
    if (found!=string::npos)
    {
      // Get the number as a string
      string number = test.substr(found+12);
      // Remove blancks
      while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
      // Get the value
      m_AttnDimX = atoi(number.c_str());
      // Increment the number of data found
      data_found++;
    }
    // Image dimensions for inveon or focus: y
    found = test.find("y_dimension");
    if (found!=string::npos)
    {
      // Get the number as a string
      string number = test.substr(found+12);
      // Remove blancks
      while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
      // Get the value
      m_AttnDimY = atoi(number.c_str());
      // Increment the number of data found
      data_found++;
    }
    // Image dimensions for inveon or focus: z
    found = test.find("z_dimension");
    if (found!=string::npos)
    {
      // Get the number as a string
      string number = test.substr(found+12);
      // Remove blancks
      while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
      // Get the value
      m_AttnDimZ = atoi(number.c_str());
      // Increment the number of data found
      data_found++;
    }
    // Voxel size for biograph, hrplus and hrrt
    if (test.find("scale factor")!=string::npos || test.find("scaling factor")!=string::npos)
    {
      // Test 5: for the voxel size X
      found = test.find("[1]");
      if (found!=string::npos)
      {
        // Get the number as a string
        found = test.find("=");
        string number = test.substr(found+1);
        // Remove blancks
        while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
        // Get the value
        m_AttnVoxSizeX = (double)atof(number.c_str());
        // Increment the number of data found
        data_found++;
      }
      // Test 6: for the voxel size Y
      found = test.find("[2]");
      if (found!=string::npos)
      {
        // Get the number as a string
        found = test.find("=");
        string number = test.substr(found+1);
        // Remove blancks
        while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
        // Get the value
        m_AttnVoxSizeY = (double)atof(number.c_str());
        // Increment the number of data found
        data_found++;
      }
      // Test 7: for the voxel size Z
      found = test.find("[3]");
      if (found!=string::npos)
      {
        // Get the number as a string
        found = test.find("=");
        string number = test.substr(found+1);
        // Remove blancks
        while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
        // Get the value
        m_AttnVoxSizeZ = (double)atof(number.c_str());
        // Increment the number of data found
        data_found++;
      }
    }
    // Voxel dimensions for inveon or focus: x
    found = test.find("pixel_size_x");
    if (found!=string::npos)
    {
      // Get the number as a string
      string number = test.substr(found+13);
      // Remove blancks
      while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
      // Get the value
      m_AttnVoxSizeX = (PRECISION)atof(number.c_str());
      // Increment the number of data found
      data_found++;
    }
    // Voxel dimensions for inveon or focus: y
    found = test.find("pixel_size_y");
    if (found!=string::npos)
    {
      // Get the number as a string
      string number = test.substr(found+13);
      // Remove blancks
      while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
      // Get the value
      m_AttnVoxSizeY = (PRECISION)atof(number.c_str());
      // Increment the number of data found
      data_found++;
    }
    // Voxel dimensions for inveon or focus: z
    found = test.find("pixel_size_z");
    if (found!=string::npos)
    {
      // Get the number as a string
      string number = test.substr(found+13);
      // Remove blancks
      while ( (found=number.find(" ")) != string::npos) number.replace(found,1,"");
      // Get the value
      m_AttnVoxSizeZ = (PRECISION)atof(number.c_str());
      // Increment the number of data found
      data_found++;
    }
    // Test 8: for the patient orientation
    found = test.find("atient orientation");
    if (found!=string::npos)
    {
      // Get the file name
      found = test.find("=");
      m_Orientation = test.substr(found+1);
      // Remove blancks and cariage return
      while ( (found=m_Orientation.find(" ")) != string::npos)  m_Orientation.replace(found,1,"");
      while ( (found=m_Orientation.find("\r")) != string::npos) m_Orientation.replace(found,1,"");
      // Do not increment the number of data found because this option is not mendatory
    }
    // Test 9: half voxel shift (if come from a simulation then no shift)
    found = test.find("alf voxel shift");
    if (found!=string::npos)
    {
      // Get the answer
      found = test.find("=");
      string answer = test.substr(found+1);
      // Remove blancks and cariage return
      while ( (found=answer.find(" ")) != string::npos)  answer.replace(found,1,"");
      while ( (found=answer.find("\r")) != string::npos) answer.replace(found,1,"");
      if (answer=="no") m_AttnHalfVoxelShift = false;
      // Do not increment the number of data found because this option is not mendatory
    }
    // Read a new line
    fhead.getline(line,1024);
  }

  // Close header file
  fhead.close();

  // Check if all data were found
  if (data_found<nb_data_to_find)
  {
    LogCerr ("***** oSinoModeCreation::InitAttenuation() -> Failed to get all data (" << nb_data_to_find << ") in header file '" << f_FileAttn << "' (actually found: " << data_found << ") !" << endl);
    return 1;
  }

  // Verbose
  if (!m_AttnHalfVoxelShift && m_Verbose>=1) LogCout ("  --> No half voxel shift applied here" << endl);

  // Allocate u-map (except for inveon or focus, we allocate the attenuation sinogram)
  m_AttnDimXY  = m_AttnDimX * m_AttnDimY;
  m_AttnDimTot = m_AttnDimX * m_AttnDimY * m_AttnDimZ;
  mp_AttnUMap = (float*)malloc(m_AttnDimTot*sizeof(float));

  // Add the relative path to the data file name
  int pos; if ((pos=f_FileAttn.find_last_of("/"))!=string::npos && f_FileAttn.find_first_of("/")!=0) data_name = f_FileAttn.substr(0,pos) + "/" + data_name;

  // Open data file
  FILE* fdata = fopen(data_name.c_str(),"rb");
  if (fdata==NULL)
  {
    LogCerr ("***** oSinoModeCreation::InitAttenuation() -> Input data file '" << data_name << "' is missing or corrupted !" << endl);
    return 1;
  }

  // Verbose
  if (m_Verbose>=1) LogCout ("  --> Reading file '" << data_name << "' ..." << endl);

  // Read u-map (if not inveon)
  int nb_data_read = fread(&mp_AttnUMap[0],sizeof(float),m_AttnDimTot,fdata);

  // Close data file
  fclose(fdata);

  // Check
  if (nb_data_read!=m_AttnDimTot)
  {
    LogCerr ("***** oSinoModeCreation::InitAttenuation() -> Failed to read all data (" << m_AttnDimTot << ") in u-map (" << nb_data_read << " read) !" << endl);
    return 1;
  }

  // Deal with patient orientation
  if ( mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH || mp_Scanner->GetScannerModel()==SCANNER_HRPLUS ||mp_Scanner->GetScannerModel()==SCANNER_VISION600 ||
       mp_Scanner->GetScannerModel()==SCANNER_BIOGRAPH2D || mp_Scanner->GetScannerModel()==SCANNER_INVEON ||
       mp_Scanner->GetScannerModel()==SCANNER_FOCUS || mp_Scanner->GetScannerModel()==SCANNER_CASTOR || mp_Scanner->GetScannerModel()==SCANNER_MMR2D )
  {
    // Switch on different possible orientations
    if (m_Orientation=="HFS" || m_Orientation=="UNKNOWN" || m_Orientation=="Y") // Unknown is the default constructor value
    {
      // Flip mu-map along y axis
      Misc_FlipY(mp_AttnUMap, m_AttnDimX, m_AttnDimY, m_AttnDimZ);
      // Set axis orientation for projections
      m_AxisOrientationX = 1.;
      m_AxisOrientationY = -1.;
    }
    else if (m_Orientation=="FFS")
    {
      // Flip mu-map along y and z axis
      Misc_FlipY(mp_AttnUMap, m_AttnDimX, m_AttnDimY, m_AttnDimZ);
      Misc_FlipZ(mp_AttnUMap, m_AttnDimX, m_AttnDimY, m_AttnDimZ);
      // Set axis orientation for projections
      m_AxisOrientationX = 1.;
      m_AxisOrientationY = -1.;
    }
    else if (m_Orientation=="HFP" || m_Orientation=="ORIGINAL")
    {
      // Set axis orientation for projections
      m_AxisOrientationX = 1.;
      m_AxisOrientationY = 1.;
    }
    else if (m_Orientation=="FFP")
    {
      // Flip mu-map along z axis
      Misc_FlipZ(mp_AttnUMap, m_AttnDimX, m_AttnDimY, m_AttnDimZ);
      // Set axis orientation for projections
      m_AxisOrientationX = 1.;
      m_AxisOrientationY = 1.;
    }
  }
  else
  {
    // Set axis orientation for projections
    m_AxisOrientationX = 1.;
    m_AxisOrientationY = 1.;
  }

  // Verbose
  if (m_Verbose>=1)
  {
    LogCout ("  --> Dimensions [" << m_AttnDimX << ";" << m_AttnDimY << ";" << m_AttnDimZ << "]" << endl);
    LogCout ("  --> Voxel size [" << m_AttnVoxSizeX << ";" << m_AttnVoxSizeY << ";" << m_AttnVoxSizeZ << "] mm" << endl);
    LogCout ("  --> FOV size [" << m_AttnVoxSizeX*((double)m_AttnDimX) << ";" << m_AttnVoxSizeY*((double)m_AttnDimY) << ";" << m_AttnVoxSizeZ*((double)m_AttnDimZ) << "] mm" << endl);
    LogCout ("  --> Patient orientation: " << m_Orientation << endl);
  }

  // Special case for inveon and focus where the umap is in mm-1, so we change it to cm-1
  if (mp_Scanner->GetScannerModel()==SCANNER_INVEON || mp_Scanner->GetScannerModel()==SCANNER_FOCUS) for (int v=0; v<m_AttnDimTot; v++) mp_AttnUMap[v] *= 10.;

  // Set the attenuation correction to true
  m_AttnCorr = true;

  // Ending
  return 0;
}

// ==========================================================================================================================================
// int ComputeACF()
//    --> This function computes the attenuation correction factor by forward-projection in the u-map using SIDDON.
// ==========================================================================================================================================
float oSinoModeCreation::ComputeACF(int f_CrystalID1, int f_CrystalID2)
{
  // Correction or not
  if (!m_AttnCorr) return 1.;

  // Verbose
  if (m_Verbose>=5) LogCout ("oSinoModeCreation::ComputeACF() -> Compute acf for crystal IDs [" << f_CrystalID1 << ";" << f_CrystalID2 << "]" << endl);

  // Get coordinates
  double x1 = (mp_Scanner->GetCornerX1(f_CrystalID1)+mp_Scanner->GetCornerX2(f_CrystalID1))/2.;
  double y1 = (mp_Scanner->GetCornerY1(f_CrystalID1)+mp_Scanner->GetCornerY2(f_CrystalID1))/2.;
  double z1 = (mp_Scanner->GetCornerZ1(f_CrystalID1)+mp_Scanner->GetCornerZ2(f_CrystalID1))/2.;
  double x2 = (mp_Scanner->GetCornerX1(f_CrystalID2)+mp_Scanner->GetCornerX2(f_CrystalID2))/2.;
  double y2 = (mp_Scanner->GetCornerY1(f_CrystalID2)+mp_Scanner->GetCornerY2(f_CrystalID2))/2.;
  double z2 = (mp_Scanner->GetCornerZ1(f_CrystalID2)+mp_Scanner->GetCornerZ2(f_CrystalID2))/2.;

  // Project the ACF
  double acf = SiddonAttenuationProjection(x1,y1,z1,x2,y2,z2);

  // Return ACF
  return ((float)acf);
}

// ==========================================================================================================================================
// double SiddonAttenuationProjection()
//    --> This function project a line through the u-map and return the exponential of the calculated value (in mm-1).
// ==========================================================================================================================================
// The Siddon function
double oSinoModeCreation::SiddonAttenuationProjection( double x1, double y1, double z1, double x2, double y2, double z2 )
{

  // Shift in z half the detector size
  z1 -= mp_Scanner->GetAxialScannerSize()/2. - mp_Scanner->GetAxialCrystalSize()/2.;
  z2 -= mp_Scanner->GetAxialScannerSize()/2. - mp_Scanner->GetAxialCrystalSize()/2.;

  // Half voxel shifting in X and Y directions
  if (m_AttnHalfVoxelShift)
  {
    x1 += m_AxisOrientationX*m_AttnVoxSizeX/2.;
    x2 += m_AxisOrientationX*m_AttnVoxSizeX/2.;
    y1 += m_AxisOrientationY*m_AttnVoxSizeY/2.;
    y2 += m_AxisOrientationY*m_AttnVoxSizeY/2.;
  }

  if (m_SiddonDidier)
  {

  PRECISION transmission = 0.;

  // Didier's stuff
  long int nPlane_[3];
  nPlane_[0] = m_AttnDimX+1;
  nPlane_[1] = m_AttnDimY+1;
  nPlane_[2] = m_AttnDimZ+1;
  PRECISION fovSizeX = ((PRECISION)m_AttnDimX)*m_AttnVoxSizeX;
  PRECISION fovSizeY = ((PRECISION)m_AttnDimY)*m_AttnVoxSizeY;
  PRECISION fovSizeZ = ((PRECISION)m_AttnDimZ)*m_AttnVoxSizeZ;
  PRECISION xPlane_[2]; xPlane_[1] = fovSizeX / 2.; xPlane_[0] = -xPlane_[1];
  PRECISION yPlane_[2]; yPlane_[1] = fovSizeY / 2.; yPlane_[0] = -yPlane_[1];
  PRECISION zPlane_[2]; zPlane_[1] = fovSizeZ / 2.; zPlane_[0] = -zPlane_[1];

	// Computing the distance between axis
	PRECISION const dX21 = x2 - x1;
	PRECISION const dY21 = y2 - y1;
	PRECISION const dZ21 = z2 - z1;

	// Computing distance between the two points
	PRECISION const d21 = sqrt( dX21 * dX21 + dY21 * dY21 + dZ21 * dZ21 );

	// Computing the parametric value alphaMin and alphaMax
	PRECISION alphaMin = 0.0, alphaMax = 1.0;

	// First step is to compute the extrem alpha for each axis
	// For the X-axis
	if( dX21 != 0.0 )
	{
		PRECISION const alphaX_1  = - ( fovSizeX/2. + x1 ) / dX21;
		PRECISION const alphaX_NX =   ( fovSizeX/2. - x1 ) / dX21;
		alphaMin = max( alphaMin, min( alphaX_1, alphaX_NX ) );
		alphaMax = min( alphaMax, max( alphaX_1, alphaX_NX ) );
	}

	// For the Y-axis
	if( dY21 != 0.0 )
	{
		PRECISION const alphaY_1  = - ( fovSizeY/2. + y1 ) / dY21;
		PRECISION const alphaY_NY =   ( fovSizeY/2. - y1 ) / dY21;
		alphaMin = max( alphaMin, min( alphaY_1, alphaY_NY ) );
		alphaMax = min( alphaMax, max( alphaY_1, alphaY_NY ) );
	}

	// For the Z-axis
	if( dZ21 != 0.0 )
	{
		PRECISION const alphaZ_1  = - ( fovSizeZ/2. + z1 ) / dZ21;
		PRECISION const alphaZ_NZ =   ( fovSizeZ/2. - z1 ) / dZ21;
		alphaMin = max( alphaMin, min( alphaZ_1, alphaZ_NZ ) );
		alphaMax = min( alphaMax, max( alphaZ_1, alphaZ_NZ ) );
	}

	// if alphaMax is less than or equal to alphaMin no intersection
	// and return an empty buffer
	if( alphaMax <= alphaMin )
  {
    // No voxels crossed
    return 1.;
  }

	// Now we have to find the indices of the particular plane
	// (iMin,iMax), (jMin,jMax), (kMin,kMax)
	long int iMin = 0, iMax = 0;
	long int jMin = 0, jMax = 0;
	long int kMin = 0, kMax = 0;

	// For the X-axis
	if( dX21 > 0.0 )
	{
		iMin = ceil( nPlane_[ 0 ] - ( xPlane_[ 1 ] - alphaMin * dX21 - x1 ) / m_AttnVoxSizeX );
		iMax = floor( 1 + ( x1 + alphaMax * dX21 - xPlane_[ 0 ] ) / m_AttnVoxSizeX );
	}
	else if( dX21 < 0.0 )
	{
		iMin = ceil( nPlane_[ 0 ] - ( xPlane_[ 1 ] - alphaMax * dX21 - x1 ) / m_AttnVoxSizeX );
		iMax = floor( 1 + ( x1 + alphaMin * dX21 - xPlane_[ 0 ] ) / m_AttnVoxSizeX );
	}
  else
	{
		iMin = 1, iMax = 0;
	}

	// For the Y-axis
	if( dY21 > 0.0 )
	{
		jMin = ceil( nPlane_[ 1 ] - ( yPlane_[ 1 ] - alphaMin * dY21 - y1 ) / m_AttnVoxSizeY );
		jMax = floor( 1 + ( y1 + alphaMax * dY21 - yPlane_[ 0 ] ) / m_AttnVoxSizeY );
	}
	else if( dY21 < 0.0 )
	{
		jMin = ceil( nPlane_[ 1 ] - ( yPlane_[ 1 ] - alphaMax * dY21 - y1 ) / m_AttnVoxSizeY );
		jMax = floor( 1 + ( y1 + alphaMin * dY21 - yPlane_[ 0 ] ) / m_AttnVoxSizeY );
	}
  else
	{
		jMin = 1, jMax = 0;
	}

	// For the Z-axis
	if( dZ21 > 0.0 )
	{
		kMin = ceil( nPlane_[ 2 ] - ( zPlane_[ 1 ] - alphaMin * dZ21 - z1 ) / m_AttnVoxSizeZ );
		kMax = floor( 1 + ( z1 + alphaMax * dZ21 - zPlane_[ 0 ] ) / m_AttnVoxSizeZ );
	}
	else if( dZ21 < 0.0 )
	{
		kMin = ceil( nPlane_[ 2 ] - ( zPlane_[ 1 ] - alphaMax * dZ21 - z1 ) / m_AttnVoxSizeZ );
		kMax = floor( 1 + ( z1 + alphaMin * dZ21 - zPlane_[ 0 ] ) / m_AttnVoxSizeZ );
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
			*idx++ = ( ( xPlane_[ 0 ] + ( i - 1 ) * m_AttnVoxSizeX ) - x1 ) / dX21;
		}
		else if( dX21 < 0 )
		{
			for( long int i = iMax; i >= iMin; --i )
				*idx++ = ( ( xPlane_[ 0 ] + ( i - 1 ) * m_AttnVoxSizeX ) - x1 ) / dX21;
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
				*idx++ = ( ( yPlane_[ 0 ] + ( j - 1 ) * m_AttnVoxSizeY ) - y1 ) / dY21;
		}
		else if( dY21 < 0 )
		{
			for( long int j = jMax; j >= jMin; --j )
				*idx++ = ( ( yPlane_[ 0 ] + ( j - 1 ) * m_AttnVoxSizeY ) - y1 ) / dY21;
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
				*idx++ = ( ( zPlane_[ 0 ] + ( k - 1 ) * m_AttnVoxSizeZ ) - z1 ) / dZ21;
		}
		else if( dZ21 < 0 )
		{
			for( long int k = kMax; k >= kMin; --k )
				*idx++ = ( ( zPlane_[ 0 ] + ( k - 1 ) * m_AttnVoxSizeZ ) - z1 ) / dZ21;
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
		i = 1 + ( x1 + alphaMid * dX21 - xPlane_[ 0 ] ) / m_AttnVoxSizeX;
		j = 1 + ( y1 + alphaMid * dY21 - yPlane_[ 0 ] ) / m_AttnVoxSizeY;
		k = 1 + ( z1 + alphaMid * dZ21 - zPlane_[ 0 ] ) / m_AttnVoxSizeZ;

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
                        index = (i-1) + (j-1)*m_AttnDimX + (k-1)*m_AttnDimXY;

                        // Update transmission values
                        transmission += length * mp_AttnUMap[index];
	}

            // Return the ACF (the division by 10 is here to convert cm-1 to mm-1)
            return max(exp(transmission/10.),1.);
  }
  // Standard siddon (including plane errors)
  else
  {
  // Parameters and variables
  int k1, ck;
  double xp, dl;
  double r1xyz[3], dxyz[3], pn_xyz[3];
  double alpha_1[3]={0.0, 0.0, 0.0}, alpha_n[3]={1.0, 1.0, 1.0}, alpha_min, alpha_max, oldvalue, almin[3], almax[3];
  double d_inc[3], d_xyz[3],i_min[3], i_max[3];
  int step_ijk[3];
  long i_c[3];
  int mxyz [] = {m_AttnDimX, m_AttnDimY, m_AttnDimZ};
  double acf = 0.;

/* ALREADY DONE ABOVE
  // Shift in z half the detector size
  z1 -= mp_Scanner->GetAxialScannerSize()/2.;
  z2 -= mp_Scanner->GetAxialScannerSize()/2.;

  // Half voxel shifting in X and Y directions
  if (m_AttnHalfVoxelShift)
  {
    x1 += m_AxisOrientationX*m_AttnVoxSizeX/2.;
    x2 += m_AxisOrientationX*m_AttnVoxSizeX/2.;
    y1 += m_AxisOrientationY*m_AttnVoxSizeY/2.;
    y2 += m_AxisOrientationY*m_AttnVoxSizeY/2.;
  }
*/
  // Voxel sizes
  double vox[3];
  vox[0] = m_AttnVoxSizeX;
  vox[1] = m_AttnVoxSizeY;
  vox[2] = m_AttnVoxSizeZ;

  dxyz[0] = (double)(x1-x2);
  dxyz[1] = (double)(y1-y2);
  dxyz[2] = (double)(z1-z2);
  r1xyz[0] = (double) x2;
  r1xyz[1] = (double) y2;
  r1xyz[2] = ((double) z2);
  dl=sqrt((dxyz[0]*dxyz[0])+(dxyz[1]*dxyz[1])+(dxyz[2]*dxyz[2]));
  pn_xyz[0] =((double) m_AttnDimX)/2.*m_AttnVoxSizeX;
  pn_xyz[1] =((double) m_AttnDimY)/2.*m_AttnVoxSizeY;
  pn_xyz[2] =((double) m_AttnDimZ)/2.*m_AttnVoxSizeZ;

  for (k1=0;k1<3;k1++)
  {
    if (fabs(dxyz[k1]) > 1.0e-5)
    {
      alpha_1[k1] = ( -pn_xyz[k1] - r1xyz[k1])/dxyz[k1];
      alpha_n[k1] = ( pn_xyz[k1] - r1xyz[k1])/dxyz[k1];
    }
  }

  alpha_min=0.;
  alpha_max=1.;
  for(k1=0;k1<3;k1++)
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
  if (alpha_min >= alpha_max) return 1.;

  for(k1 = 0; k1<3; k1++)
  {
    if (fabs(dxyz[k1]) >1.0e-5)
    {
      if (dxyz[k1] >0.)
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
          if (fabs(d_xyz[k1] - alpha_min > 1.0e-6))
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
      d_xyz[k1] = 2.*dl;
      i_c[k1] = (int) ((pn_xyz[k1] + r1xyz[k1] + 1.0e-4)/vox[k1]);
    }
  }		  

  ck = 0;
  for (oldvalue = alpha_min * dl; (i_c[0] < m_AttnDimX && i_c[1] < m_AttnDimY && i_c[2] < m_AttnDimZ && i_c[0]>-1 && i_c[1]>-1 && i_c[2]>-1);)
  {
    ck = (d_xyz[0] < d_xyz[1] ) ? 0 : 1;
    ck = (d_xyz[ck] < d_xyz[2] ) ? ck : 2;
    {
      long int index = i_c[2]*m_AttnDimXY + i_c[1]*m_AttnDimX + i_c[0];
      acf += (d_xyz[ck]-oldvalue) * mp_AttnUMap[index];
    }
    oldvalue = d_xyz[ck];
    if (fabs(d_xyz[0] - oldvalue) < 1.0e-6) {d_xyz[0] += d_inc[0];i_c[0] += step_ijk[0];}
    if (fabs(d_xyz[1] - oldvalue) < 1.0e-6) {d_xyz[1] += d_inc[1];i_c[1] += step_ijk[1];}
    if (fabs(d_xyz[2] - oldvalue) < 1.0e-6) {d_xyz[2] += d_inc[2];i_c[2] += step_ijk[2];}
  }

  // Return the ACF (the division by 10 is here to convert cm-1 to mm-1)
  return max(exp(acf/10.),1.);

  }
}

