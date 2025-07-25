// This file is not a class but contains functions usefull for other classes
#include "oMiscellaneous.hh"

// ========================================================================================================================
// Function that initialize the random generator
unsigned int Misc_InitRandomGenerator(int seed)
{
  // The final seed
  unsigned int the_seed = 0;
  // Switch on method
  if (seed<0) the_seed = ((unsigned int)(time(NULL)*getpid()));
  else the_seed = ((unsigned int)seed);
  // Initialize the generator
  srand(the_seed);
  // Return the seed
  return the_seed;
}

// ========================================================================================================================
// Big to little endian swap functions
float Misc_EndianSwap(const float aFloat)
{
  // The result
  float result;
  // Convert float* to char*
  char *inputFloat = (char*) &aFloat;
  char *outputFloat = (char*) &result;
  // Swap the bytes one per one
  outputFloat[0] = inputFloat[3];
  outputFloat[1] = inputFloat[2];
  outputFloat[2] = inputFloat[1];
  outputFloat[3] = inputFloat[0];
  // Return
  return result;
}

// ========================================================================================================================
// Functions to get isotope half-life and branching ratio (from ecat7 documentation)
// It returns 0 if the isotope is in the library, 1 if not.
int Misc_GetIsotopeCharacteristics(const string& isotope, PRECISION* half_life, PRECISION* branching_ratio)
{
  // Note: Half-lives are given in seconds
  // Note: First line is a hard-coded name when the isotope is not given
  if (isotope=="")           { *half_life = 0.;        *branching_ratio = 1.;    }
  else if (isotope=="Br-75") { *half_life = 5802.;     *branching_ratio = 0.725; }
  else if (isotope=="C-11")  { *half_life = 1223.;     *branching_ratio = 1.;    }
  else if (isotope=="Cu-62") { *half_life = 584.4;     *branching_ratio = 0.98;  }
  else if (isotope=="Cu-64") { *half_life = 45720.;    *branching_ratio = 0.179; }
  else if (isotope=="F-18")  { *half_life = 6586.2;    *branching_ratio = 0.967; }
  else if (isotope=="Mn-52") { *half_life = 483062.;   *branching_ratio = 0.296; }
  else if (isotope=="Ga-68") { *half_life = 4057.7;    *branching_ratio = 0.89;  }
  else if (isotope=="Ge-68") { *half_life = 23399000.; *branching_ratio = 0.891; }
  else if (isotope=="N-13")  { *half_life = 597.9;     *branching_ratio = 1.;    }
  else if (isotope=="O-14")  { *half_life = 70.606;    *branching_ratio = 1.;    }
  else if (isotope=="O-15")  { *half_life = 122.24;    *branching_ratio = 1.;    }
  else if (isotope=="Rb-82") { *half_life = 76.38;     *branching_ratio = 0.955; }
  else if (isotope=="Na-22") { *half_life = 82110000.; *branching_ratio = 0.9;   }
  else if (isotope=="Zn-62") { *half_life = 33070.;    *branching_ratio = 1.064; }
  else if (isotope=="Br-76") { *half_life = 58320.;    *branching_ratio = 0.57;  }
  else if (isotope=="K-38")  { *half_life = 458.2;     *branching_ratio = 0.995; }
  else if (isotope=="INF")   { *half_life = 0.;        *branching_ratio = 1.;    }
  else                       { *half_life = 0.;        *branching_ratio = 1.;      return 1; }
  return 0;
}

// ========================================================================================================================
// Functions to flip an image
// ------------------------------>  Overlay function
void Misc_FlipImage(const string& flip, float* image, int dimX, int dimY, int dimZ)
{
  if (flip=="NONE") return;
  else if (flip=="X")   Misc_FlipX   (image,dimX,dimY,dimZ);
  else if (flip=="Y")   Misc_FlipY   (image,dimX,dimY,dimZ);
  else if (flip=="Z")   Misc_FlipZ   (image,dimX,dimY,dimZ);
  else if (flip=="XY")  Misc_FlipXY  (image,dimX,dimY,dimZ);
  else if (flip=="XZ")  Misc_FlipXZ  (image,dimX,dimY,dimZ);
  else if (flip=="YZ")  Misc_FlipYZ  (image,dimX,dimY,dimZ);
  else if (flip=="XYZ") Misc_FlipXYZ (image,dimX,dimY,dimZ);
}
// ------------------------------>  Flip X
void Misc_FlipX  (float* image, int dimX, int dimY, int dimZ)
{
  // Dimensions
  int dimXY = dimX * dimY;
  int dimTot = dimXY * dimZ;

  // Allocate buffer
  float* buffer = (float*)malloc(dimTot*sizeof(float));

  // Flip
  for (int z=0; z<dimZ; z++) for (int y=0; y<dimY; y++) for (int x=0; x<dimX; x++)
  {
    int indexBuffer = z*dimXY + y*dimX + (dimX-1-x);
    int indexImage  = z*dimXY + y*dimX + x;
    buffer[indexBuffer] = image[indexImage];
  }

  // Copy buffer into image
  for (int v=0; v<dimTot; v++) image[v] = buffer[v];

  // Free buffer
  free(buffer);
}
// ------------------------------>  Flip Y
void Misc_FlipY  (float* image, int dimX, int dimY, int dimZ)
{
  // Dimensions
  int dimXY = dimX * dimY;
  int dimTot = dimXY * dimZ;

  // Allocate buffer
  float* buffer = (float*)malloc(dimTot*sizeof(float));

  // Flip
  for (int z=0; z<dimZ; z++) for (int y=0; y<dimY; y++) for (int x=0; x<dimX; x++)
  {
    int indexBuffer = z*dimXY + (dimY-1-y)*dimX + x;
    int indexImage  = z*dimXY + y*dimX + x;
    buffer[indexBuffer] = image[indexImage];
  }

  // Copy buffer into image
  for (int v=0; v<dimTot; v++) image[v] = buffer[v];

  // Free buffer
  free(buffer);
}
// ------------------------------>  Flip Z
void Misc_FlipZ  (float* image, int dimX, int dimY, int dimZ)
{
  // Dimensions
  int dimXY = dimX * dimY;
  int dimTot = dimXY * dimZ;

  // Allocate buffer
  float* buffer = (float*)malloc(dimTot*sizeof(float));

  // Flip
  for (int z=0; z<dimZ; z++) for (int y=0; y<dimY; y++) for (int x=0; x<dimX; x++)
  {
    int indexBuffer = (dimZ-1-z)*dimXY + y*dimX + x;
    int indexImage  = z*dimXY + y*dimX + x;
    buffer[indexBuffer] = image[indexImage];
  }

  // Copy buffer into image
  for (int v=0; v<dimTot; v++) image[v] = buffer[v];

  // Free buffer
  free(buffer);
}
// ------------------------------>  Flip XY
void Misc_FlipXY (float* image, int dimX, int dimY, int dimZ)
{
  // Dimensions
  int dimXY = dimX * dimY;
  int dimTot = dimXY * dimZ;

  // Allocate buffer
  float* buffer = (float*)malloc(dimTot*sizeof(float));

  // Flip
  for (int z=0; z<dimZ; z++) for (int y=0; y<dimY; y++) for (int x=0; x<dimX; x++)
  {
    int indexBuffer = z*dimXY + (dimY-1-y)*dimX + (dimX-1-x);
    int indexImage  = z*dimXY + y*dimX + x;
    buffer[indexBuffer] = image[indexImage];
  }

  // Copy buffer into image
  for (int v=0; v<dimTot; v++) image[v] = buffer[v];

  // Free buffer
  free(buffer);
}
// ------------------------------>  Flip XZ
void Misc_FlipXZ (float* image, int dimX, int dimY, int dimZ)
{
  // Dimensions
  int dimXY = dimX * dimY;
  int dimTot = dimXY * dimZ;

  // Allocate buffer
  float* buffer = (float*)malloc(dimTot*sizeof(float));

  // Flip
  for (int z=0; z<dimZ; z++) for (int y=0; y<dimY; y++) for (int x=0; x<dimX; x++)
  {
    int indexBuffer = (dimZ-1-z)*dimXY + y*dimX + (dimX-1-x);
    int indexImage  = z*dimXY + y*dimX + x;
    buffer[indexBuffer] = image[indexImage];
  }

  // Copy buffer into image
  for (int v=0; v<dimTot; v++) image[v] = buffer[v];

  // Free buffer
  free(buffer);
}
// ------------------------------>  Flip YZ
void Misc_FlipYZ (float* image, int dimX, int dimY, int dimZ)
{
  // Dimensions
  int dimXY = dimX * dimY;
  int dimTot = dimXY * dimZ;

  // Allocate buffer
  float* buffer = (float*)malloc(dimTot*sizeof(float));

  // Flip
  for (int z=0; z<dimZ; z++) for (int y=0; y<dimY; y++) for (int x=0; x<dimX; x++)
  {
    int indexBuffer = z*dimXY + (dimY-1-y)*dimX + (dimX-1-x);
    int indexImage  = z*dimXY + y*dimX + x;
    buffer[indexBuffer] = image[indexImage];
  }

  // Copy buffer into image
  for (int v=0; v<dimTot; v++) image[v] = buffer[v];

  // Free buffer
  free(buffer);
}
// ------------------------------>  Flip XYZ
void Misc_FlipXYZ (float* image, int dimX, int dimY, int dimZ)
{
  // Dimensions
  int dimXY = dimX * dimY;
  int dimTot = dimXY * dimZ;

  // Allocate buffer
  float* buffer = (float*)malloc(dimTot*sizeof(float));

  // Flip
  for (int z=0; z<dimZ; z++) for (int y=0; y<dimY; y++) for (int x=0; x<dimX; x++)
  {
    int indexBuffer = (dimZ-1-z)*dimXY + (dimY-1-y)*dimX + (dimX-1-x);
    int indexImage  = z*dimXY + y*dimX + x;
    buffer[indexBuffer] = image[indexImage];
  }

  // Copy buffer into image
  for (int v=0; v<dimTot; v++) image[v] = buffer[v];

  // Free buffer
  free(buffer);
}
// ------------------------------>  Flip X
void Misc_FlipX  (double* image, int dimX, int dimY, int dimZ)
{
  // Dimensions
  int dimXY = dimX * dimY;
  int dimTot = dimXY * dimZ;

  // Allocate buffer
  double* buffer = (double*)malloc(dimTot*sizeof(double));

  // Flip
  for (int z=0; z<dimZ; z++) for (int y=0; y<dimY; y++) for (int x=0; x<dimX; x++)
  {
    int indexBuffer = z*dimXY + y*dimX + (dimX-1-x);
    int indexImage  = z*dimXY + y*dimX + x;
    buffer[indexBuffer] = image[indexImage];
  }

  // Copy buffer into image
  for (int v=0; v<dimTot; v++) image[v] = buffer[v];

  // Free buffer
  free(buffer);
}
// ------------------------------>  Flip Y
void Misc_FlipY  (double* image, int dimX, int dimY, int dimZ)
{
  // Dimensions
  int dimXY = dimX * dimY;
  int dimTot = dimXY * dimZ;

  // Allocate buffer
  double* buffer = (double*)malloc(dimTot*sizeof(double));

  // Flip
  for (int z=0; z<dimZ; z++) for (int y=0; y<dimY; y++) for (int x=0; x<dimX; x++)
  {
    int indexBuffer = z*dimXY + (dimY-1-y)*dimX + x;
    int indexImage  = z*dimXY + y*dimX + x;
    buffer[indexBuffer] = image[indexImage];
  }

  // Copy buffer into image
  for (int v=0; v<dimTot; v++) image[v] = buffer[v];

  // Free buffer
  free(buffer);
}
// ------------------------------>  Flip Z
void Misc_FlipZ  (double* image, int dimX, int dimY, int dimZ)
{
  // Dimensions
  int dimXY = dimX * dimY;
  int dimTot = dimXY * dimZ;

  // Allocate buffer
  double* buffer = (double*)malloc(dimTot*sizeof(double));

  // Flip
  for (int z=0; z<dimZ; z++) for (int y=0; y<dimY; y++) for (int x=0; x<dimX; x++)
  {
    int indexBuffer = (dimZ-1-z)*dimXY + y*dimX + x;
    int indexImage  = z*dimXY + y*dimX + x;
    buffer[indexBuffer] = image[indexImage];
  }

  // Copy buffer into image
  for (int v=0; v<dimTot; v++) image[v] = buffer[v];

  // Free buffer
  free(buffer);
}
// ------------------------------>  Flip XY
void Misc_FlipXY (double* image, int dimX, int dimY, int dimZ)
{
  // Dimensions
  int dimXY = dimX * dimY;
  int dimTot = dimXY * dimZ;

  // Allocate buffer
  double* buffer = (double*)malloc(dimTot*sizeof(double));

  // Flip
  for (int z=0; z<dimZ; z++) for (int y=0; y<dimY; y++) for (int x=0; x<dimX; x++)
  {
    int indexBuffer = z*dimXY + (dimY-1-y)*dimX + (dimX-1-x);
    int indexImage  = z*dimXY + y*dimX + x;
    buffer[indexBuffer] = image[indexImage];
  }

  // Copy buffer into image
  for (int v=0; v<dimTot; v++) image[v] = buffer[v];

  // Free buffer
  free(buffer);
}
// ------------------------------>  Flip XZ
void Misc_FlipXZ (double* image, int dimX, int dimY, int dimZ)
{
  // Dimensions
  int dimXY = dimX * dimY;
  int dimTot = dimXY * dimZ;

  // Allocate buffer
  double* buffer = (double*)malloc(dimTot*sizeof(double));

  // Flip
  for (int z=0; z<dimZ; z++) for (int y=0; y<dimY; y++) for (int x=0; x<dimX; x++)
  {
    int indexBuffer = (dimZ-1-z)*dimXY + y*dimX + (dimX-1-x);
    int indexImage  = z*dimXY + y*dimX + x;
    buffer[indexBuffer] = image[indexImage];
  }

  // Copy buffer into image
  for (int v=0; v<dimTot; v++) image[v] = buffer[v];

  // Free buffer
  free(buffer);
}
// ------------------------------>  Flip YZ
void Misc_FlipYZ (double* image, int dimX, int dimY, int dimZ)
{
  // Dimensions
  int dimXY = dimX * dimY;
  int dimTot = dimXY * dimZ;

  // Allocate buffer
  double* buffer = (double*)malloc(dimTot*sizeof(double));

  // Flip
  for (int z=0; z<dimZ; z++) for (int y=0; y<dimY; y++) for (int x=0; x<dimX; x++)
  {
    int indexBuffer = z*dimXY + (dimY-1-y)*dimX + (dimX-1-x);
    int indexImage  = z*dimXY + y*dimX + x;
    buffer[indexBuffer] = image[indexImage];
  }

  // Copy buffer into image
  for (int v=0; v<dimTot; v++) image[v] = buffer[v];

  // Free buffer
  free(buffer);
}
// ------------------------------>  Flip XYZ
void Misc_FlipXYZ (double* image, int dimX, int dimY, int dimZ)
{
  // Dimensions
  int dimXY = dimX * dimY;
  int dimTot = dimXY * dimZ;

  // Allocate buffer
  double* buffer = (double*)malloc(dimTot*sizeof(double));

  // Flip
  for (int z=0; z<dimZ; z++) for (int y=0; y<dimY; y++) for (int x=0; x<dimX; x++)
  {
    int indexBuffer = (dimZ-1-z)*dimXY + (dimY-1-y)*dimX + (dimX-1-x);
    int indexImage  = z*dimXY + y*dimX + x;
    buffer[indexBuffer] = image[indexImage];
  }

  // Copy buffer into image
  for (int v=0; v<dimTot; v++) image[v] = buffer[v];

  // Free buffer
  free(buffer);
}

// ========================================================================================================================
// Function to round a PRECISION to the nearest integer value
int Round_Me(PRECISION x)
{
  int integer_part = floor(x);
  PRECISION decimal_part = x - ((PRECISION)integer_part);
  if (decimal_part<=((PRECISION)0.5)) return integer_part;
  else return integer_part+1;
}

// ========================================================================================================================
// 


