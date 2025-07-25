#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <fstream>
#include <omp.h>
using namespace std;

void showHelp()
{
  cout << endl;
  cout << "Usage: create_phantom  -o  baseName  -d  voxX  voxY  voxZ  -v  sizeX  sizeY  sizeZ  [Options]" << endl;
  cout << endl;
  cout << "  This program creates an interfile image in float 32bits format named 'baseName.img'. The image is" << endl;
  cout << "  composed of 'voxX', 'voxY', 'voxZ' voxels of size 'sizeX', 'sizeY', 'sizeZ' given in mm." << endl;
  cout << "  In this image, cylinders, spheres, ellipsoids, boxes, etc, can be inserted with different" << endl;
  cout << "  voxel values. Every distances and positions are given in mm. The last object declared wins" << endl;
  cout << "  the voxel value." << endl;
  cout << endl;
  cout << "  -o  baseName  : is the base name of the output image" << endl;
  cout << "  -d  vox[XYZ]  : is the number of voxels in the 3 dimensions" << endl;
  cout << "  -v  size[XYZ] : is the voxel dimensions (in mm)" << endl;
  cout << endl;
  cout << "  [Options]" << endl;
  cout << endl;
  cout << "  -s  pos[XYZ] rad val            : to insert a sphere of radius 'rad' and centered in position [posX,posY,posZ]," << endl;
  cout << "                                    given in mm. The point [0,0,0] is the center of the image. 'val' will be the" << endl;
  cout << "                                    value inside the sphere." << endl;
  cout << "  -e  pos[XYZ] radX radY leng val : to insert an object of length 'len' and of eliptical section (with 'radX' and" << endl;
  cout << "                                    'radY' as radius), and centered in position [posX,posY,posZ], all given in mm." << endl;
  cout << "                                    The point [0,0,0] is the center of the image. 'val' will be the value inside" << endl;
  cout << "                                    the object." << endl;
  cout << "  -c  pos[XYZ] rad len val        : to insert a cylinder of radius 'rad' and length 'len' and centered in" << endl;
  cout << "                                    position [posX,posY,posZ], given in mm. The point [0,0,0] is the center" << endl;
  cout << "                                    of the image. 'val' will be the value inside the cylinder." << endl;
  cout << "  -b  pos[XYZ] lenX lenY lenZ val : to insert a box of length 'lenX', 'lenY' and 'lenZ', and centered in position" << endl;
  cout << "                                    [posX,posY,posZ], given in mm. The point [0,0,0] is the center of the image." << endl;
  cout << "                                    'val' will be the value inside the box." << endl;
  cout << "  -p  pos[XYZ] val                : to insert a point in position [posX,posY,posZ], given in pixel. 'val' will be" << endl;
  cout << "                                    be the value of the point (i.e. the voxel)." << endl;
  cout << "  -bg  val                        : to fix the background value in the image (default: 0.)" << endl;
  cout << "  -i   file                       : use the given file as initialization (assume the same dimensions)" << endl;
  cout << "  -f   val                        : to multiply the whole image by a given value (default: 1.)" << endl;
  cout << "  -x  nb[XYZ]                     : perform a surpixelation to get smoother volume (default: 1 1 1)" << endl;	
  cout << "  -xx                             : when considering a surpixelation, any voxel containing a valid sub-voxel will be" << endl;
  cout << "                                    consider as valid too" << endl;
  cout << "  -yy                             : when considering a surpixelation, any voxel NOT containing all sub-voxels will be" << endl;
  cout << "                                    excluded" << endl;
  cout << "  -t  nb                          : to use multiple threads and give the number of threads" << endl;
  cout << "  -h                              : this help page" << endl;
  cout << endl;
  cout << "  Return codes:" << endl;
  cout << endl;
  cout << "  0: normal termination" << endl;
  cout << "  1: error while reading parameters" << endl;
  cout << "  2: error while allocating memory" << endl;
  cout << "  3: error while writing on disk" << endl;
  cout << endl;
  exit(0);
}

void errorMessage(const string& message, int code)
{
  cerr << endl;
  cerr << "***** Error message:" << endl;
  cerr << " => " << message << endl;
  cerr << "Type 'create_attenuation_map -h' to get help." << endl;
  cerr << endl;
  exit(code);
}

void warningMessage(const string& message)
{
  cerr << endl;
  cerr << "***** Warning message:" << endl;
  cerr << " => " << message << endl;
  cerr << endl;
}
/**************************************************************************************************************
 *
 *       ADDCYLINDER
 *
 **************************************************************************************************************/
void addCylinder( float*** image,
                  int dimX, int dimY, int dimZ,
                  int pixX, int pixY, int pixZ,
                  float sizeX, float sizeY, float sizeZ,
                  float posX, float posY, float posZ,
                  float radius, float length, float value,
                  bool subValid, bool subReject, int nbThread )
{
  // Calculate cylinder bounds
  float boxSizeX = radius*2./sizeX;
  float boxSizeY = radius*2./sizeY;
  float boxSizeZ = length/sizeZ;
  float boxPosX  = (posX+((float)dimX*sizeX/2.))/sizeX;
  float boxPosY  = (posY+((float)dimY*sizeY/2.))/sizeY;
  float boxPosZ  = (posZ+((float)dimZ*sizeZ/2.))/sizeZ;
  int boundXStart = ((int)(boxPosX-boxSizeX/2.)) - 1;
  int boundXStop  = ((int)(boxPosX+boxSizeX/2.)) + 1;
  int boundYStart = ((int)(boxPosY-boxSizeY/2.)) - 1;
  int boundYStop  = ((int)(boxPosY+boxSizeY/2.)) + 1;
  int boundZStart = ((int)(boxPosZ-boxSizeZ/2.)) - 1;
  int boundZStop  = ((int)(boxPosZ+boxSizeZ/2.)) + 1;
  if (boundXStart<0) boundXStart = 0;
  if (boundYStart<0) boundYStart = 0;
  if (boundZStart<0) boundZStart = 0;
  if (boundXStop>dimX-1) boundXStop = dimX-1;
  if (boundYStop>dimY-1) boundYStop = dimY-1;
  if (boundZStop>dimZ-1) boundZStop = dimZ-1;
  // Sur-pixelisation
  int pixDim = pixX*pixY*pixZ;
  double sizePixX = sizeX/((double)pixX);
  double sizePixY = sizeY/((double)pixY);
  double sizePixZ = sizeZ/((double)pixZ);
  // OPENMP
  int z; omp_set_num_threads(nbThread);
  #pragma omp parallel for private(z) schedule(static,4)
  for (z=boundZStart; z<=boundZStop; z++)
  {
    // For optimization
    bool isInY = false;
    for (int y=boundYStart; y<=boundYStop; y++)
    {
      // For optimization
      bool isInX = false;
      for (int x=boundXStart; x<=boundXStop; x++)
      {
        bool subValidBreak = false;
        int subValidPixels = 0;
        double mean = 0.;
        // BEGIN: Sub-pixel loops
        for (int xx=0; xx<pixX; xx++)
        {
          for (int yy=0; yy<pixY; yy++)
          {
            for (int zz=0; zz<pixZ; zz++)
            {
              double coordX = ( ((double)x)-((double)dimX)/2. )*sizeX + ((double)xx)*sizePixX + sizePixX/2. - posX; // position du centre
              double coordY = ( ((double)y)-((double)dimY)/2. )*sizeY + ((double)yy)*sizePixY + sizePixY/2. - posY; // du survoxel dans le
              double coordZ = ( ((double)z)-((double)dimZ)/2. )*sizeZ + ((double)zz)*sizePixZ + sizePixZ/2. - posZ; // referentiel cylindre
              if ( ((coordX*coordX + coordY*coordY) <= radius*radius) && (coordZ>=-length/2.) && (coordZ<=length/2.) )
              {
                mean += value;
                subValidPixels++;
                if (subValid)
                {
                  subValidBreak = true;
                  break;
                }
                // For optimization
                isInX = true;
              }
              else mean += image[x][y][z];
            }
            if (subValidBreak) break;
          }
          if (subValidBreak) break;
        }
        // END: Sub-pixel loops
        // Test for optimization
        isInY = isInX;
        if (isInX && mean == image[x][y][z]/((double)pixDim)) break;
        // Update image value
        if (subValid && subValidBreak) image[x][y][z] = value;
        else if (subReject && subValidPixels!=pixDim) continue;
        else image[x][y][z] = mean/((double)pixDim);
      }
      // Test for optimization
      if (isInY && !isInX) break;
    }
  }
}
/**************************************************************************************************************
 *
 *       ADDSPHERE
 *
 **************************************************************************************************************/
void addSphere( float*** image,
                int dimX, int dimY, int dimZ,
                int pixX, int pixY, int pixZ,
                float sizeX, float sizeY, float sizeZ,
                float posX, float posY, float posZ,
                float radius, float value,
                bool subValid, bool subReject, int nbThread )
{
  // Calculate sphere bounds
  float boxSizeX = radius*2./sizeX;
  float boxSizeY = radius*2./sizeY;
  float boxSizeZ = radius*2./sizeZ;
  float boxPosX  = (posX+((float)dimX*sizeX/2.))/sizeX;
  float boxPosY  = (posY+((float)dimY*sizeY/2.))/sizeY;
  float boxPosZ  = (posZ+((float)dimZ*sizeZ/2.))/sizeZ;
  int boundXStart = ((int)(boxPosX-boxSizeX/2.)) - 1;
  int boundXStop  = ((int)(boxPosX+boxSizeX/2.)) + 1;
  int boundYStart = ((int)(boxPosY-boxSizeY/2.)) - 1;
  int boundYStop  = ((int)(boxPosY+boxSizeY/2.)) + 1;
  int boundZStart = ((int)(boxPosZ-boxSizeZ/2.)) - 1;
  int boundZStop  = ((int)(boxPosZ+boxSizeZ/2.)) + 1;
  if (boundXStart<0) boundXStart = 0;
  if (boundYStart<0) boundYStart = 0;
  if (boundZStart<0) boundZStart = 0;
  if (boundXStop>dimX-1) boundXStop = dimX-1;
  if (boundYStop>dimY-1) boundYStop = dimY-1;
  if (boundZStop>dimZ-1) boundZStop = dimZ-1;
  // Sur-pixelization
  int pixDim = pixX*pixY*pixZ;
  double sizePixX = sizeX/((double)pixX);
  double sizePixY = sizeY/((double)pixY);
  double sizePixZ = sizeZ/((double)pixZ);
  // OPENMP
  int x; omp_set_num_threads(nbThread);
  #pragma omp parallel for private(x) schedule(static,4)
  for (x=boundXStart; x<=boundXStop; x++)
  {
    // For optimization
    bool isInY = false;
    for (int y=boundYStart; y<=boundYStop; y++)
    {
      // For optimization
      bool isInZ = false;
      for (int z=boundZStart; z<=boundZStop; z++)
      {
        bool subValidBreak = false;
        int subValidPixels = 0;
        double mean = 0.;
        // BEGIN: Sub-pixel loops
        for (int xx=0; xx<pixX; xx++)
        {
          for (int yy=0; yy<pixY; yy++)
          {
            for (int zz=0; zz<pixZ; zz++)
            {
              double coordX = ( ((double)x)-((double)dimX)/2. )*sizeX + ((double)xx)*sizePixX + sizePixX/2. - posX; // position du centre
              double coordY = ( ((double)y)-((double)dimY)/2. )*sizeY + ((double)yy)*sizePixY + sizePixY/2. - posY; // du survoxel dans le
              double coordZ = ( ((double)z)-((double)dimZ)/2. )*sizeZ + ((double)zz)*sizePixZ + sizePixZ/2. - posZ; // referentiel sphere
              if ( ((coordX*coordX + coordY*coordY + coordZ*coordZ ) <= radius*radius) )
              {
                mean += value;
                subValidPixels++;
                if (subValid)
                {
                  subValidBreak = true;
                  break;
                }
                // For optimization
                isInZ = true;
              }
              else mean += image[x][y][z];
            }
            if (subValidBreak) break;
          }
          if (subValidBreak) break;
        }
        // END: Sub-pixel loops
        // Test for optimization
        isInY = isInZ;
        if (isInZ && mean == image[x][y][z]/((double)pixDim)) break;
        // Update image value
        if (subValid && subValidBreak) image[x][y][z] = value;
        else if (subReject && subValidPixels!=pixDim) continue;
        else image[x][y][z] = mean/((double)pixDim);
      }
      // Test for optimization
      if (isInY && !isInZ) break;
    }
  }
}
/**************************************************************************************************************
 *
 *       ADDELIPTICAL
 *
 **************************************************************************************************************/
void addEliptical( float*** image,
                   int dimX, int dimY, int dimZ,
                   int pixX, int pixY, int pixZ,
                   float sizeX, float sizeY, float sizeZ,
                   float posX, float posY, float posZ,
                   float radiusX, float radiusY, float length, float value,
                   bool subValid, bool subReject, int nbThread )
{
  int pixDim = pixX*pixY*pixZ;
  double sizePixX = sizeX/((double)pixX);
  double sizePixY = sizeY/((double)pixY);
  double sizePixZ = sizeZ/((double)pixZ);
  int x; omp_set_num_threads(nbThread);
  #pragma omp parallel for private(x) schedule(static,4)
  for (x=0; x<dimX; x++)
  {
    for (int y=0; y<dimY; y++)
    {
      for (int z=0; z<dimZ; z++)
      {
        bool subValidBreak = false;
        int subValidPixels = 0;
        double mean = 0.;
        for (int xx=0; xx<pixX; xx++)
        {
          for (int yy=0; yy<pixY; yy++)
          {
            for (int zz=0; zz<pixZ; zz++)
            {
              double coordX = ( ((double)x)-((double)dimX)/2. )*sizeX + ((double)xx)*sizePixX + sizePixX/2. - posX; // position du centre
              double coordY = ( ((double)y)-((double)dimY)/2. )*sizeY + ((double)yy)*sizePixY + sizePixY/2. - posY; // du survoxel dans le
              double coordZ = ( ((double)z)-((double)dimZ)/2. )*sizeZ + ((double)zz)*sizePixZ + sizePixZ/2. - posZ; // referentiel elipse
              if ( ((coordX*coordX/(radiusX*radiusX) + coordY*coordY/(radiusY*radiusY)) <= 1.) && (coordZ>=-length/2.) && (coordZ<=length/2.) )
              {
                mean += value;
                subValidPixels++;
                if (subValid)
                {
                  subValidBreak = true;
                  break;
                }
              }
              else mean += image[x][y][z];
            }
            if (subValidBreak) break;
          }
          if (subValidBreak) break;
        }
        if (subValid && subValidBreak) image[x][y][z] = value;
        else if (subReject && subValidPixels!=pixDim) continue;
        else image[x][y][z] = mean/((double)pixDim);
      }
    }
  }
}
/**************************************************************************************************************
 *
 *       ADDBOX
 *
 **************************************************************************************************************/
void addBox( float*** image,
             int dimX, int dimY, int dimZ,
             int pixX, int pixY, int pixZ,
             float sizeX, float sizeY, float sizeZ,
             float posX, float posY, float posZ,
             float lengthX, float lengthY, float lengthZ, float value,
             bool subValid, bool subReject, int nbThread )
{
  int pixDim = pixX*pixY*pixZ;
  double sizePixX = sizeX/((double)pixX);
  double sizePixY = sizeY/((double)pixY);
  double sizePixZ = sizeZ/((double)pixZ);
  int x; omp_set_num_threads(nbThread);
  #pragma omp parallel for private(x) schedule(static,4)
  for (x=0; x<dimX; x++)
  {
    for (int y=0; y<dimY; y++)
    {
      for (int z=0; z<dimZ; z++)
      {
        bool subValidBreak = false;
        int subValidPixels = 0;
        double mean = 0.;
        for (int xx=0; xx<pixX; xx++)
        {
          for (int yy=0; yy<pixY; yy++)
          {
            for (int zz=0; zz<pixZ; zz++)
            {
              double coordX = ( ((double)x)-((double)dimX)/2. )*sizeX + ((double)xx)*sizePixX + sizePixX/2. - posX; // position du centre
              double coordY = ( ((double)y)-((double)dimY)/2. )*sizeY + ((double)yy)*sizePixY + sizePixY/2. - posY; // du survoxel dans le
              double coordZ = ( ((double)z)-((double)dimZ)/2. )*sizeZ + ((double)zz)*sizePixZ + sizePixZ/2. - posZ; // referentiel elipse
              if ( (coordX>=-lengthX/2.) && (coordX<=lengthX/2.) && (coordY>=-lengthY/2.) && (coordY<=lengthY/2.) && (coordZ>=-lengthZ/2.) && (coordZ<=lengthZ/2.) )
              {
                mean += value;
                subValidPixels++;
                if (subValid)
                {
                  subValidBreak = true;
                  break;
                }
              }
              else mean += image[x][y][z];
            }
            if (subValidBreak) break;
          }
          if (subValidBreak) break;
        }
        if (subValid && subValidBreak) image[x][y][z] = value;
        else if (subReject && subValidPixels!=pixDim) continue;
        else image[x][y][z] = mean/((double)pixDim);
      }
    }
  }
}
/**************************************************************************************************************
 *
 *       ADDPOINT
 *
 **************************************************************************************************************/
void addPoint(float*** image, int dimX, int dimY, int dimZ, int posX, int posY, int posZ, float value)
{
  if (posX>=dimX)
  {
    cerr << "***** Unable to add a point in image since the X position of the voxel (" << posX << ") is greater than the X dimension (" << dimX << ") !" << endl;
    exit(1);
  }
  if (posY>=dimY)
  {
    cerr << "***** Unable to add a point in image since the Y position of the voxel (" << posY << ") is greater than the Y dimension (" << dimY << ") !" << endl;
    exit(1);
  }
  if (posZ>=dimZ)
  {
    cerr << "***** Unable to add a point in image since the Z position of the voxel (" << posZ << ") is greater than the Z dimension (" << dimZ << ") !" << endl;
    exit(1);
  }
  image[posX][posY][posZ] = value;
}
/**************************************************************************************************************
 *
 *       INIT FROM IMAGE
 *
 **************************************************************************************************************/
int initFromImage(float*** image, int dimX, int dimY, int dimZ, const string& file)
{
  // Open file
  FILE* fin = fopen(file.c_str(),"rb");
  if (fin==NULL)
  {
    cerr << "***** Input file '" << file << "' is missing or corrupted !" << endl;
    return 1;
  }
  // Init image
  int nb_read_voxels = 0;
  for (int z=0; z<dimZ; z++) for (int y=0; y<dimY; y++) for (int x=0; x<dimX; x++)
  {
    nb_read_voxels += fread(&image[x][y][z],sizeof(float),1,fin);
  }
  // Close file
  fclose(fin);
  // Check
  if (nb_read_voxels!=dimX*dimY*dimZ)
  {
    cerr << "***** Input image '" << file << "' does not have the given dimensions !" << endl;
    return 1;
  }
  // End
  return 0;
}
/**************************************************************************************************************
 *
 *       CHECK NON-BKGD VOXELS
 *
 **************************************************************************************************************/
int computeNonBackgroundVoxels(float*** image, int dimX, int dimY, int dimZ, float background, int nbThread)
{
  int x; omp_set_num_threads(nbThread);
  int* nb_voxels = (int*)calloc(nbThread,sizeof(int));
  #pragma omp parallel for private(x) schedule(static,4)
  for (x=0; x<dimX; x++)
  {
    int th = omp_get_thread_num();
    for (int y=0; y<dimY; y++)
    {
      for (int z=0; z<dimZ; z++)
      {
        if (image[x][y][z]!=background) nb_voxels[th]++;
      }
    }
  }

  int nb_total_voxels = 0;
  for (int th=0; th<nbThread; th++) nb_total_voxels += nb_voxels[th];
  return nb_total_voxels;
}
/**************************************************************************************************************
 *
 *       MAIN
 *
 **************************************************************************************************************/
int main (int argc, char** argv)
{
  // Number of threads
  int nbThread = 1;
  // On fixe ici le nombre d'objets differents existants
  int nbObjectTypes = 5;
  // Voici la liste des objets differents existants
  enum {CYLINDER=1, SPHERE=2, ELIPSOID=3, BOX=4, POINT=5};
  // Pas le temps de faire du dynamique donc on fixe un nombre d'objets max par type d'objet
  int nbMaxObjectsPerType = 1000;
  // Ceci sera la liste qui contiendra uniquement des valeurs indiquant l'ordre des objets a ajouter
  int listOfObjects[nbObjectTypes*nbMaxObjectsPerType];
  for (int i=0; i<nbObjectTypes*nbMaxObjectsPerType; i++) listOfObjects[i]=-1;

  // Compteurs pour chaque type d'objet
  int nbSph = 0; int cntSph = 0;
  int nbCyl = 0; int cntCyl = 0;
  int nbEli = 0; int cntEli = 0;
  int nbBox = 0; int cntBox = 0;
  int nbPnt = 0; int cntPnt = 0;
  // Compteur pour le nombre total d'objet
  int nbObjects = 0;

  // Structure contenant les infos pour les cylindres:
  //   - le premier [] designe le numero du cylindre,
  //   - le deuxieme [] designe tous les parametres du cylindre (0: posX, 1: posY, 2: posZ, 3: radius, 4: length, 5: value)
  float** cylinders = (float**)malloc(nbMaxObjectsPerType*sizeof(float*));
  for (int c=0; c<nbMaxObjectsPerType; c++) cylinders[c] = (float*)malloc(6*sizeof(float));

  // Structure contenant les infos pour les spheres:
  //   - le premier [] designe le numero de la sphere,
  //   - le deuxieme [] designe tous les parametres de la sphere (0: posX, 1: posY, 2: posZ, 3: radius, 4: value)
  float** spheres = (float**)malloc(nbMaxObjectsPerType*sizeof(float*));
  for (int c=0; c<nbMaxObjectsPerType; c++) spheres[c] = (float*)malloc(5*sizeof(float));

  // Structure contenant les infos pour les elipses:
  //   - le premier [] designe le numero de l'elipse,
  //   - le deuxieme [] designe tous les parametres de l'elipse (0: posX, 1: posY, 2: posZ, 3: radiusX, 4: radiusY, 5: length, 6: value)
  float** elipses = (float**)malloc(nbMaxObjectsPerType*sizeof(float*));
  for (int c=0; c<nbMaxObjectsPerType; c++) elipses[c] = (float*)malloc(7*sizeof(float));

  // Structure contenant les infos pour les boites:
  //   - le premier [] designe le numero de la boite,
  //   - le deuxieme [] designe tous les parametres de la boite (0: posX, 1: posY, 2: posZ, 3: lengthX, 4: lengthY, 5: lengthZ, 6: value)
  float** boxes = (float**)malloc(nbMaxObjectsPerType*sizeof(float*));
  for (int c=0; c<nbMaxObjectsPerType; c++) boxes[c] = (float*)malloc(7*sizeof(float));

  // Structure contenant les infos pour les voxels:
  //   - le premier [] designe le numero du voxel,
  //   - le deuxieme [] designe tous les parametres du voxel (0: posX, 1: posY, 2: posZ, 3: value)
  float** points = (float**)malloc(nbMaxObjectsPerType*sizeof(float*));
  for (int c=0; c<nbMaxObjectsPerType; c++) points[c] = (float*)malloc(4*sizeof(float));

  //===============================================================================================================//
  //  LECTURE DES PARAMETRES
  //===============================================================================================================//

  int argument=0;
  string option="";

  int dimX=0, dimY=0, dimZ=0;
  int pixX=1, pixY=1, pixZ=1;
  float sizeX=0., sizeY=0., sizeZ=0.;
  float background = 0.;
  float scaleFactor = 1.;
  string initFile = "";
  bool subValid  = false;
  bool subReject = false;
  string baseName="";

  if (argc<9) showHelp();
  for (int i=1; i<argc; i++)
  {
    option = (string)argv[i];
    if (argument==0)
    {
      if (option=="-o") argument=10;        // nom de l'image
      else if (option=="-d")  argument=20;  // dimensions de l'image
      else if (option=="-v")  argument=30;  // tailles du voxel
      else if (option=="-bg") argument=40;  // changement du background
      else if (option=="-s")  argument=50;  // creation d'une sphere
      else if (option=="-c")  argument=60;  // creation d'un cylindre
      else if (option=="-e")  argument=70;  // creation d'une elipse
      else if (option=="-b")  argument=80;  // creation d'une boite
      else if (option=="-p")  argument=90;  // creation d'un point
      else if (option=="-x")  argument=100; // surpixelisation
      else if (option=="-t")  argument=110; // nombre de threads
      else if (option=="-f")  argument=120; // facteur d'echelle
      else if (option=="-i")  argument=130; // image d'initialisation
      else if (option=="-xx") subValid = true;
      else if (option=="-yy") subReject = true;
      else if (option=="-h") showHelp();
      else errorMessage("Unknown option '"+option+"' !",1);
    }
    else if (argument==10) // nom de l'image
    {
      baseName = option;
      argument = 0;
    }
    else if (argument==20) // dimensions de l'image (n°1: dimX)
    {
      dimX = (int)atoi(option.c_str());
      argument=21;
    }
    else if (argument==21) // dimensions de l'image (n°2: dimY)
    {
      dimY = (int)atoi(option.c_str());
      argument=22;
    }
    else if (argument==22) // dimensions de l'image (n°3: dimZ)
    {
      dimZ = (int)atoi(option.c_str());
      argument=0;
    }
    else if (argument==100) // surpixelisation de l'image (n°1: pixX)
    {
      pixX = (int)atoi(option.c_str());
      argument=101;
    }
    else if (argument==101) // surpixelisation de l'image (n°2: pixY)
    {
      pixY = (int)atoi(option.c_str());
      argument=102;
    }
    else if (argument==102) // surpixelisation de l'image (n°3: pixZ)
    {
      pixZ = (int)atoi(option.c_str());
      argument=0;
    }
    else if (argument==110) // nombre de threads
    {
      nbThread = atoi(option.c_str());
      argument=0;
    }
    else if (argument==120) // facteur d'echelle
    {
      scaleFactor = atof(option.c_str());
      argument=0;
    }
    else if (argument==130) // image d'initialisation
    {
      initFile = option;
      argument=0;
    }
    else if (argument==30) // tailles du voxel (n°1: sizeX)
    {
      sizeX = atof(option.c_str());
      argument=31;
    }
    else if (argument==31) // tailles du voxel (n°2: sizeY)
    {
      sizeY = atof(option.c_str());
      argument=32;
    }
    else if (argument==32) // tailles du voxel (n°3: sizeZ)
    {
      sizeZ = atof(option.c_str());
      argument=0;
    }
    else if (argument==40) // changement de background
    {
      background = atof(option.c_str());
      argument=0;
    }
    else if (argument==50) // creation d'une sphere (n°1: posX)
    {
      if (nbSph==nbMaxObjectsPerType) {
        char tmp[10];
        sprintf(tmp,"%d",nbMaxObjectsPerType);
        errorMessage("Maximum number of spheres reached ("+((string)tmp)+") !",1);
      }
      spheres[nbSph][0] = atof(option.c_str());
      argument=51;
    }
    else if (argument==51) // creation d'une sphere (n°2: posY)
    {
      spheres[nbSph][1] = atof(option.c_str());
      argument=52;
    }
    else if (argument==52) // creation d'une sphere (n°3: posZ)
    {
      spheres[nbSph][2] = atof(option.c_str());
      argument=53;
    }
    else if (argument==53) // creation d'une sphere (n°4: radius)
    {
      spheres[nbSph][3] = atof(option.c_str());
      argument=54;
    }
    else if (argument==54) // creation d'une sphere (n°5: value)
    {
      spheres[nbSph][4] = atof(option.c_str());
      nbSph++;
      listOfObjects[nbObjects]=SPHERE;
      nbObjects++;
      argument=0;
    }
    else if (argument==60) // creation d'un cylindre (n°1: posX)
    {
      if (nbCyl==nbMaxObjectsPerType) {
        char tmp[10];
        sprintf(tmp,"%d",nbMaxObjectsPerType);
        errorMessage("Maximum number of cylinders reached ("+((string)tmp)+") !",1);
      }
      cylinders[nbCyl][0] = atof(option.c_str());
      argument=61;
    }
    else if (argument==61) // creation d'un cylindre (n°2: posY)
    {
      cylinders[nbCyl][1] = atof(option.c_str());
      argument=62;
    }
    else if (argument==62) // creation d'un cylindre (n°3: posZ)
    {
      cylinders[nbCyl][2] = atof(option.c_str());
      argument=63;
    }
    else if (argument==63) // creation d'un cylindre (n°4: radius)
    {
      cylinders[nbCyl][3] = atof(option.c_str());
      argument=64;
    }
    else if (argument==64) // creation d'un cylindre (n°5: length)
    {
      cylinders[nbCyl][4] = atof(option.c_str());
      argument=65;
    }
    else if (argument==65) // creation d'un cylindre (n°6: value)
    {
      cylinders[nbCyl][5] = atof(option.c_str());
      nbCyl++;
      listOfObjects[nbObjects]=CYLINDER;
      nbObjects++;
      argument=0;
    }
    else if (argument==70) // creation d'une elipse (n°1: posX)
    {
      if (nbEli==nbMaxObjectsPerType) {
        char tmp[10];
        sprintf(tmp,"%d",nbMaxObjectsPerType);
        errorMessage("Maximum number of eliptical reached ("+((string)tmp)+") !",1);
      }
      elipses[nbEli][0] = atof(option.c_str());
      argument=71;
    }
    else if (argument==71) // creation d'une elipse (n°2: posY)
    {
      elipses[nbEli][1] = atof(option.c_str());
      argument=72;
    }
    else if (argument==72) // creation d'une elipse (n°3: posZ)
    {
      elipses[nbEli][2] = atof(option.c_str());
      argument=73;
    }
    else if (argument==73) // creation d'une elipse (n°4: radiusX)
    {
      elipses[nbEli][3] = atof(option.c_str());
      argument=74;
    }
    else if (argument==74) // creation d'une elipse (n°5: radiusY)
    {
      elipses[nbEli][4] = atof(option.c_str());
      argument=75;
    }
    else if (argument==75) // creation d'une elipse (n°6: length)
    {
      elipses[nbEli][5] = atof(option.c_str());
      argument=76;
    }
    else if (argument==76) // creation d'une elipse (n°7: value)
    {
      elipses[nbEli][6] = atof(option.c_str());
      nbEli++;
      listOfObjects[nbObjects]=ELIPSOID;
      nbObjects++;
      argument=0;
    }
    else if (argument==80) // creation d'une boite (n°1: posX)
    {
      if (nbBox==nbMaxObjectsPerType) {
        char tmp[10];
        sprintf(tmp,"%d",nbMaxObjectsPerType);
        errorMessage("Maximum number of boxes reached ("+((string)tmp)+") !",1);
      }
      boxes[nbBox][0] = atof(option.c_str());
      argument=81;
    }
    else if (argument==81) // creation d'une boite (n°2: posY)
    {
      boxes[nbBox][1] = atof(option.c_str());
      argument=82;
    }
    else if (argument==82) // creation d'une boite (n°3: posZ)
    {
      boxes[nbBox][2] = atof(option.c_str());
      argument=83;
    }
    else if (argument==83) // creation d'une boite (n°4: lengthX)
    {
      boxes[nbBox][3] = atof(option.c_str());
      argument=84;
    }
    else if (argument==84) // creation d'une boite (n°5: lengthY)
    {
      boxes[nbBox][4] = atof(option.c_str());
      argument=85;
    }
    else if (argument==85) // creation d'une boite (n°6: lengthZ)
    {
      boxes[nbBox][5] = atof(option.c_str());
      argument=86;
    }
    else if (argument==86) // creation d'une boite (n°7: value)
    {
      boxes[nbBox][6] = atof(option.c_str());
      nbBox++;
      listOfObjects[nbObjects]=BOX;
      nbObjects++;
      argument=0;
    }
    else if (argument==90) // creation d'un point (n°1: posX)
    {
      if (nbPnt==nbMaxObjectsPerType) {
        char tmp[10];
        sprintf(tmp,"%d",nbMaxObjectsPerType);
        errorMessage("Maximum number of points reached ("+((string)tmp)+") !",1);
      }
      points[nbPnt][0] = atof(option.c_str());
      argument=91;
    }
    else if (argument==91) // creation d'un point (n°2: posY)
    {
      points[nbPnt][1] = atof(option.c_str());
      argument=92;
    }
    else if (argument==92) // creation d'un point (n°3: posZ)
    {
      points[nbPnt][2] = atof(option.c_str());
      argument=93;
    }
    else if (argument==93) // creation d'un point (n°4: value)
    {
      points[nbPnt][3] = atof(option.c_str());
      nbPnt++;
      listOfObjects[nbObjects]=POINT;
      nbObjects++;
      argument=0;
    }
    if (i==argc-1 && argument!=0) errorMessage("Some parameters are missing !",1);
  }

  // Check des parametres
  if (dimX==0 || dimY==0 || dimZ==0)
  {
    errorMessage("One of the dimensions is 0 !",1);
  }
  if (sizeX<=0. || sizeY<=0. || sizeZ<=0.)
  {
    errorMessage("Voxel dimensions must be strickly positive !",1);
  }
  if (baseName=="")
  {
    errorMessage("Maybe giving an output file name will be better !",1);
  }
  if (nbObjects==0)
  {
    //warningMessage("No object will be created, resulting in an empty image !");
  }
  if (nbThread<1)
  {
    errorMessage("Thread number is not valid !",1);
  }
  if (subReject && subValid)
  {
    errorMessage("Cannot use -xx and -yy options together !",1);
  }

  //===============================================================================================================//
  //  CREATION DE L'IMAGE
  //===============================================================================================================//  

  // Verbose
  cout << "  --> Creating image (background: " << background << ")" << endl;

  float *** image = (float***)malloc(dimX*sizeof(float**));
  for (int x=0; x<dimX; x++)
  {
    image[x] = (float**)malloc(dimY*sizeof(float*));
    for (int y=0; y<dimY; y++)
    {
      image[x][y] = (float*)malloc(dimZ*sizeof(float));
      for (int z=0; z<dimZ; z++) image[x][y][z] = background;
    }
  }

  //===============================================================================================================//
  //  INITIALIZATION DE L'IMAGE
  //===============================================================================================================//  

  if (initFile!="")
  {
    cout << "  --> Initializing from file '" << initFile << "'" << endl;
    if (initFromImage(image, dimX, dimY, dimZ, initFile))
    {
      errorMessage("Problem while initializing image from input file !",1);
    }
  }

  //===============================================================================================================//
  //  INSERTION DES OBJETS
  //===============================================================================================================//  

  for (int n=0; n<nbObjects; n++)
  {
    if (listOfObjects[n]==CYLINDER)
    {
      // Verbose
      cout << "  --> Inserting cylinder at [" << cylinders[cntCyl][0] << ";" << cylinders[cntCyl][1] << ";" << cylinders[cntCyl][2] << "] of radius "
           << cylinders[cntCyl][3] << " and length " << cylinders[cntCyl][4] << " with value " << cylinders[cntCyl][5] << endl << flush;
      // (0: posX, 1: posY, 2: posZ, 3: radius, 4: length, 5: value)
      addCylinder(image,dimX,dimY,dimZ,pixX,pixY,pixZ,sizeX,sizeY,sizeZ,cylinders[cntCyl][0],cylinders[cntCyl][1],cylinders[cntCyl][2],cylinders[cntCyl][3],cylinders[cntCyl][4],cylinders[cntCyl][5],subValid,subReject,nbThread);
      cntCyl++;
    }
    else if (listOfObjects[n]==SPHERE)
    {
      // Verbose
      cout << "  --> Inserting sphere at [" << spheres[cntSph][0] << ";" << spheres[cntSph][1] << ";" << spheres[cntSph][2] << "] of radius "
           << spheres[cntSph][3] << " with value " << spheres[cntSph][4] << endl << flush;
      // (0: posX, 1: posY, 2: posZ, 3: radius, 4: value)
      addSphere(image,dimX,dimY,dimZ,pixX,pixY,pixZ,sizeX,sizeY,sizeZ,spheres[cntSph][0],spheres[cntSph][1],spheres[cntSph][2],spheres[cntSph][3],spheres[cntSph][4],subValid,subReject,nbThread);
      cntSph++;
    }
    else if (listOfObjects[n]==ELIPSOID)
    {
      // Verbose
      cout << "  --> Inserting elipsoid at [" << elipses[cntEli][0] << ";" << elipses[cntEli][1] << ";" << elipses[cntEli][2] << "] of radii [ "
           << elipses[cntEli][3] << ";" << elipses[cntEli][5] << "] and length " << elipses[cntEli][5] << " with value " << elipses[cntEli][6] << endl << flush;
      // (0: posX, 1: posY, 2: posZ, 3: radiusX, 4: radiusY, 5: length, 6: value)
      addEliptical(image,dimX,dimY,dimZ,pixX,pixY,pixZ,sizeX,sizeY,sizeZ,elipses[cntEli][0],elipses[cntEli][1],elipses[cntEli][2],elipses[cntEli][3],elipses[cntEli][4],elipses[cntEli][5],elipses[cntEli][6],subValid,subReject,nbThread);
      cntEli++;
    }
    else if (listOfObjects[n]==BOX)
    {
      // Verbose
      cout << "  --> Inserting box at [" << boxes[cntBox][0] << ";" << boxes[cntBox][1] << ";" << boxes[cntBox][2] << "] of dimensions ["
           << boxes[cntBox][3] << ";" << boxes[cntBox][4] << ";" << boxes[cntBox][5] << "] with value " << boxes[cntBox][6] << endl << flush;
      // (0: posX, 1: posY, 2: posZ, 3: lengthX, 4: lengthY, 5: lengthZ, 6: value)
      addBox(image,dimX,dimY,dimZ,pixX,pixY,pixZ,sizeX,sizeY,sizeZ,boxes[cntBox][0],boxes[cntBox][1],boxes[cntBox][2],boxes[cntBox][3],boxes[cntBox][4],boxes[cntBox][5],boxes[cntBox][6],subValid,subReject,nbThread);
      cntBox++;
    }
    else if (listOfObjects[n]==POINT)
    {
      // (0: posX, 1: posY, 2: posZ, 3: value)
      addPoint(image,dimX,dimY,dimZ,((int)points[cntPnt][0]),((int)points[cntPnt][1]),((int)points[cntPnt][2]),points[cntPnt][3]);
      cntPnt++;
    }
    // Compute number of non-background voxels
    int nb_voxels = computeNonBackgroundVoxels(image, dimX, dimY, dimZ, background, nbThread);
    cout << "      Number of non-background voxels: " << nb_voxels << endl << flush;
  }

  //===============================================================================================================//
  //  APPLICATION DU FACTEUR D'ECHELLE
  //===============================================================================================================//

  for (int z=0; z<dimZ; z++) for (int y=0; y<dimY; y++) for (int x=0; x<dimX; x++) image[x][y][z] *= scaleFactor;

  //===============================================================================================================//
  //  ECRITURE DE L'IMAGE
  //===============================================================================================================//

  // Binary image
  string fileImg = baseName + ".img";
  FILE* fout = fopen(fileImg.c_str(),"wb");
  if (fout==NULL) errorMessage("Cannot create output file '"+fileImg+"' !",3);
  for (int z=0; z<dimZ; z++)
  {
    for (int y=0; y<dimY; y++)
    {
      for (int x=0; x<dimX; x++)
      {
        float buffer = image[x][y][z];
        fwrite(&buffer,1,4,fout);
      }
    }
  }
  fclose(fout);

  // Interfile header
  string fileHdr = baseName + ".hdr";
  ofstream fhdr(fileHdr.c_str());
  if (!fhdr) errorMessage("Cannot create output file '"+fileHdr+"' !",3);
  fhdr << "!INTERFILE := " << endl;
  fhdr << "!imaging modality := phantom" << endl;
  fhdr << "!version of keys := CASToRv3.1" << endl;
  fhdr << "CASToR version := 3.1" << endl;
  fhdr << endl;
  fhdr << "!GENERAL DATA := " << endl;
  fhdr << "!originating system := create_phantom" << endl;
  fhdr << "!data offset in bytes := 0" << endl;
  fhdr << "!name of data file := " << fileImg << endl;
  fhdr << "patient name := " << baseName << endl;
  fhdr << endl;
  fhdr << "!GENERAL IMAGE DATA " << endl;
  fhdr << "!type of data := Static" << endl;
  fhdr << "!total number of images := " << dimZ << endl;
  fhdr << "imagedata byte order := LITTLEENDIAN" << endl;
  fhdr << "!study duration (sec) := 1" << endl;
  fhdr << endl;
  fhdr << "!STATIC STUDY (General) :=" << endl;
  fhdr << "number of dimensions := 3" << endl;
  fhdr << "!matrix size [1] := " << dimX << endl;
  fhdr << "!matrix size [2] := " << dimY << endl;
  fhdr << "!matrix size [3] := " << dimZ << endl;
  fhdr << "!number format := short float" << endl;
  fhdr << "!number of bytes per pixel := 4" << endl;
  fhdr << "scaling factor (mm/pixel) [1] := " << sizeX << endl;
  fhdr << "scaling factor (mm/pixel) [2] := " << sizeY << endl;
  fhdr << "scaling factor (mm/pixel) [3] := " << sizeZ << endl;
  fhdr << "first pixel offset (mm) [1] := 0" << endl;
  fhdr << "first pixel offset (mm) [2] := 0" << endl;
  fhdr << "first pixel offset (mm) [3] := 0" << endl;
  fhdr << "data rescale offset := 0" << endl;
  fhdr << "data rescale slope := 1" << endl;
  fhdr << "quantification units := 1" << endl;
  fhdr << "!image duration (sec) := 1" << endl;
  fhdr << "!image start time (sec) := 0" << endl;
  fhdr << "!END OF INTERFILE :=" << endl;
  fhdr.close();

  //===============================================================================================================//
  //  DESTRUCTION DES STRUCTURES
  //===============================================================================================================//

  for (int i=0; i<nbMaxObjectsPerType; i++)
  {
    free(cylinders[i]);
    free(spheres[i]);
    free(elipses[i]);
    free(boxes[i]);
    free(points[i]);
  }
  free(cylinders);
  free(spheres);
  free(elipses);
  free(boxes);
  free(points);
  for (int x=0; x<dimX; x++)
  {
    for (int y=0; y<dimY; y++)
    {
      free(image[x][y]);
    }
    free(image[x]);
  }
  free(image);

  return 0;
}

