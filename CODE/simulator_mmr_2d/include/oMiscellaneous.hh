// This file is not a class but contains functions usefull for other classes
#ifndef OMISCELLANEOUS
#define OMISCELLANEOUS 1

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <ctime>
#include <sys/types.h>
#include <unistd.h>
using namespace std;

// ========================================================================================================================
// Some defined variables for precision
#define CONVOLUTION_SIGMAS_TRUNCATION 3.5
#define TWO_SQRT_TWO_LN_2 2.354820045
#define PRECISION float
#define MPI_PRECISION MPI_FLOAT

// ========================================================================================================================
// Function that initialize the random generator
unsigned int Misc_InitRandomGenerator(int seed);

// ========================================================================================================================
// Big to little endian swap functions
float Misc_EndianSwap(const float aFloat);

// ========================================================================================================================
// Functions to get isotope half-life and branching ratio (from ecat7 documentation)
// It returns 0 if the isotope is in the library, 1 if not.
int Misc_GetIsotopeCharacteristics(const string& isotope, PRECISION* half_life, PRECISION* branching_ratio);

// ========================================================================================================================
// Functions to flip an image
void Misc_FlipImage(const string& flip, float* image, int dimX, int dimY, int dimZ);
void Misc_FlipX  (float* image, int dimX, int dimY, int dimZ);
void Misc_FlipY  (float* image, int dimX, int dimY, int dimZ);
void Misc_FlipZ  (float* image, int dimX, int dimY, int dimZ);
void Misc_FlipXY (float* image, int dimX, int dimY, int dimZ);
void Misc_FlipXZ (float* image, int dimX, int dimY, int dimZ);
void Misc_FlipYZ (float* image, int dimX, int dimY, int dimZ);
void Misc_FlipXYZ(float* image, int dimX, int dimY, int dimZ);
void Misc_FlipX  (double* image, int dimX, int dimY, int dimZ);
void Misc_FlipY  (double* image, int dimX, int dimY, int dimZ);
void Misc_FlipZ  (double* image, int dimX, int dimY, int dimZ);
void Misc_FlipXY (double* image, int dimX, int dimY, int dimZ);
void Misc_FlipXZ (double* image, int dimX, int dimY, int dimZ);
void Misc_FlipYZ (double* image, int dimX, int dimY, int dimZ);
void Misc_FlipXYZ(double* image, int dimX, int dimY, int dimZ);

// ========================================================================================================================
// Function to round
int Round_Me(PRECISION x);

// ========================================================================================================================

#endif

