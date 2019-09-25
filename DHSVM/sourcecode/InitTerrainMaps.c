/*
 * SUMMARY:      InitTerrainMaps() - Initialize terrain coverages
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Initialize terrain coverages
 * DESCRIP-END.
 * FUNCTIONS:    InitTerrainMaps()
 *               InitTopoMap()
 *               InitSoilMap()
 *               InitVegMap()
 * COMMENTS:
 * $Id: InitTerrainMaps.c,v 3.1 2013/2/3 00:08:33 Ning Exp $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "fileio.h"
#include "functions.h"
#include "constants.h"
#include "getinit.h"
#include "sizeofnt.h"
#include "slopeaspect.h"
#include "varid.h"

 /*****************************************************************************
   InitTerrainMaps()
 *****************************************************************************/
void InitTerrainMaps(LISTPTR Input, OPTIONSTRUCT *Options, MAPSIZE *Map,
  LAYER *Soil, LAYER *Veg, TOPOPIX ***TopoMap, SOILTABLE *SType, SOILPIX ***SoilMap,
  VEGTABLE *VType, VEGPIX ***VegMap)

{
  printf("\nInitializing terrain maps\n");

  InitTopoMap(Input, Options, Map, TopoMap);
  InitSoilMap(Input, Options, Map, Soil, *TopoMap, SoilMap);
  InitVegMap(Options, Input, Map, VegMap);
  if (Options->CanopyGapping)
    InitCanopyGapMap(Options, Input, Map, Soil, Veg, VType, VegMap, SType, SoilMap);
  if (Options->CanopyTiling)
    InitTileMap(Options, Input, Map, Soil, Veg, VType, VegMap, SType, SoilMap);
}

/*****************************************************************************
  InitTopoMap()
*****************************************************************************/
void InitTopoMap(LISTPTR Input, OPTIONSTRUCT * Options, MAPSIZE * Map,
  TOPOPIX *** TopoMap)
{
  const char *Routine = "InitTopoMap";
  char VarName[BUFSIZE + 1];	/* Variable name */
  int i;			/* Counter */
  int x;			/* Counter */
  int y;			/* Counter */
  int flag;         /* either or not reverse the matrix */
  int NumberType;		/* Number type of data set */
  unsigned char *Mask = NULL;	/* Basin mask */
  float *Elev;			/* Surface elevation */
  STRINIENTRY StrEnv[] = {
    {"TERRAIN", "DEM FILE", "", ""},
    {"TERRAIN", "BASIN MASK FILE", "", ""},
    {NULL, NULL, "", NULL}
  };

  /* Process the [TERRAIN] section in the input file */
  if (!(*TopoMap = (TOPOPIX **)calloc(Map->NY, sizeof(TOPOPIX *))))
    ReportError((char *)Routine, 1);
  for (y = 0; y < Map->NY; y++) {
    if (!((*TopoMap)[y] = (TOPOPIX *)calloc(Map->NX, sizeof(TOPOPIX))))
      ReportError((char *)Routine, 1);
  }

  /* Read the key-entry pairs from the input file */
  for (i = 0; StrEnv[i].SectionName; i++) {
    GetInitString(StrEnv[i].SectionName, StrEnv[i].KeyName, StrEnv[i].Default,
      StrEnv[i].VarStr, (unsigned long)BUFSIZE, Input);
    if (IsEmptyStr(StrEnv[i].VarStr))
      ReportError(StrEnv[i].KeyName, 51);
  }

  /* Read the elevation data from the DEM dataset */
  GetVarName(001, 0, VarName);
  GetVarNumberType(001, &NumberType);
  if (!(Elev = (float *)calloc(Map->NX * Map->NY,
    SizeOfNumberType(NumberType))))
    ReportError((char *)Routine, 1);

  flag = Read2DMatrix(StrEnv[demfile].VarStr, Elev, NumberType, Map, 0,
    VarName, 0);

  /* Assign the attributes to the map pixel */
  /* Reverse the matrix is flag = 1 & netcdf option is selected */
  if ((Options->FileFormat == NETCDF && flag == 0) || (Options->FileFormat == BIN)) {
    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*TopoMap)[y][x].Dem = Elev[i];
      }
    }
  }
  else if (Options->FileFormat == NETCDF && flag == 1) {
    for (y = Map->NY - 1, i = 0; y >= 0; y--) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*TopoMap)[y][x].Dem = Elev[i];
      }
    }
  }
  else ReportError((char *)Routine, 57);
  free(Elev);

  /* Read the mask */
  GetVarName(002, 0, VarName);
  GetVarNumberType(002, &NumberType);
  if (!(Mask = (unsigned char *)calloc(Map->NX * Map->NY,
    SizeOfNumberType(NumberType))))
    ReportError((char *)Routine, 1);
  flag = Read2DMatrix(StrEnv[maskfile].VarStr, Mask, NumberType, Map, 0,
    VarName, 0);

  if ((Options->FileFormat == NETCDF && flag == 0)
    || (Options->FileFormat == BIN))
  {
    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*TopoMap)[y][x].Mask = Mask[i];
      }
    }
  }
  else if (Options->FileFormat == NETCDF && flag == 1) {
    for (y = Map->NY - 1, i = 0; y >= 0; y--) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*TopoMap)[y][x].Mask = Mask[i];
      }
    }
  }
  else ReportError((char *)Routine, 57);
  free(Mask);


  /* Calculate slope, aspect, magnitude of subsurface flow gradient, and
     fraction of flow flowing in each direction based on the land surface
     slope. */
  ElevationSlopeAspect(Map, *TopoMap);

  /* After calculating the slopes and aspects for all the points, reset the
     mask if the model is to be run in point mode */
  if (Options->Extent == POINT) {
    for (y = 0; y < Map->NY; y++)
      for (x = 0; x < Map->NX; x++)
        (*TopoMap)[y][x].Mask = OUTSIDEBASIN;
    (*TopoMap)[Options->PointY][Options->PointX].Mask = (1 != OUTSIDEBASIN);
  }
  /* find out the minimum grid elevation of the basin */
  MINELEV = 9999;
  for (y = 0, i = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++, i++) {
      if (INBASIN((*TopoMap)[y][x].Mask)) {
      if ((*TopoMap)[y][x].Dem < MINELEV) {
        MINELEV = (*TopoMap)[y][x].Dem;
      }
      }
    }
  }
}

/*****************************************************************************
  InitSoilMap()
*****************************************************************************/
void InitSoilMap(LISTPTR Input, OPTIONSTRUCT * Options, MAPSIZE * Map,
  LAYER * Soil, TOPOPIX ** TopoMap, SOILPIX *** SoilMap)
{
  const char *Routine = "InitSoilMap";
  char VarName[BUFSIZE + 1];	/* Variable name */
  int i;			/* counter */
  int x;			/* counter */
  int y;			/* counter */
  int NumberType;		/* number type */
  unsigned char *Type;		/* Soil type */
  float *Depth;			/* Soil depth */
  int flag;
  STRINIENTRY StrEnv[] = {
    {"SOILS", "SOIL MAP FILE", "", ""},
    {"SOILS", "SOIL DEPTH FILE", "", ""},
    {NULL, NULL, "", NULL}
  };

  /* Process the filenames in the [SOILS] section in the input file */
  /* Assign the attributes to the correct map pixel */
  if (!(*SoilMap = (SOILPIX **)calloc(Map->NY, sizeof(SOILPIX *))))
    ReportError((char *)Routine, 1);
  for (y = 0; y < Map->NY; y++) {
    if (!((*SoilMap)[y] = (SOILPIX *)calloc(Map->NX, sizeof(SOILPIX))))
      ReportError((char *)Routine, 1);
  }

  /* Read the key-entry pairs from the input file */
  for (i = 0; StrEnv[i].SectionName; i++) {
    GetInitString(StrEnv[i].SectionName, StrEnv[i].KeyName, StrEnv[i].Default,
      StrEnv[i].VarStr, (unsigned long)BUFSIZE, Input);
    if (IsEmptyStr(StrEnv[i].VarStr))
      ReportError(StrEnv[i].KeyName, 51);
  }

  /* Read the soil type */
  GetVarName(003, 0, VarName);
  GetVarNumberType(003, &NumberType);
  if (!(Type = (unsigned char *)calloc(Map->NX * Map->NY,
    SizeOfNumberType(NumberType))))
    ReportError((char *)Routine, 1);
  flag = Read2DMatrix(StrEnv[soiltype_file].VarStr, Type, NumberType,
	Map, 0, VarName, 0);

  if ((Options->FileFormat == NETCDF && flag == 0)
    || (Options->FileFormat == BIN))
  {
    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        if (((int)Type[i]) > Soil->NTypes)
          ReportError(StrEnv[soiltype_file].VarStr, 32);
        (*SoilMap)[y][x].Soil = Type[i];
      }
    }
  }
  else if (Options->FileFormat == NETCDF && flag == 1) {
    for (y = Map->NY - 1, i = 0; y >= 0; y--) {
      for (x = 0; x < Map->NX; x++, i++) {
        if (((int)Type[i]) > Soil->NTypes)
          ReportError(StrEnv[soiltype_file].VarStr, 32);
        (*SoilMap)[y][x].Soil = Type[i];
      }
    }
  }
  else ReportError((char *)Routine, 57);

  /* Read the total soil depth  */
  GetVarName(004, 0, VarName);
  GetVarNumberType(004, &NumberType);
  if (!(Depth = (float *)calloc(Map->NX * Map->NY,
    SizeOfNumberType(NumberType))))
    ReportError((char *)Routine, 1);
  flag = Read2DMatrix(StrEnv[soildepth_file].VarStr, Depth, NumberType,
	Map, 0, VarName, 0);

  /* Assign the attributes to the correct map pixel */
  if ((Options->FileFormat == NETCDF && flag == 0)
    || (Options->FileFormat == BIN))
  {
    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*SoilMap)[y][x].Depth = Depth[i];
      }
    }
  }
  else if (Options->FileFormat == NETCDF && flag == 1) {
    for (y = Map->NY - 1, i = 0; y >= 0; y--) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*SoilMap)[y][x].Depth = Depth[i];
      }
    }
  }
  else ReportError((char *)Routine, 57);

  for (y = 0, i = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++, i++) {
      if (Options->Infiltration == DYNAMIC)
        (*SoilMap)[y][x].InfiltAcc = 0.;
      (*SoilMap)[y][x].MoistInit = 0.;

      /* allocate memory for the number of root layers, plus an additional
       layer below the deepest root layer */
      if (INBASIN(TopoMap[y][x].Mask)) {
        if (!((*SoilMap)[y][x].Moist =
          (float *)calloc((Soil->NLayers[Type[i] - 1] + 1), sizeof(float))))
          ReportError((char *)Routine, 1);
        if (!((*SoilMap)[y][x].Perc =
          (float *)calloc(Soil->NLayers[Type[i] - 1], sizeof(float))))
          ReportError((char *)Routine, 1);
        if (!((*SoilMap)[y][x].Temp =
          (float *)calloc(Soil->NLayers[Type[i] - 1], sizeof(float))))
          ReportError((char *)Routine, 1);
      }
      else {
        (*SoilMap)[y][x].Moist = NULL;
        (*SoilMap)[y][x].Perc = NULL;
        (*SoilMap)[y][x].Temp = NULL;
      }
    }
  }
  free(Type);
  free(Depth);
}

/*****************************************************************************
  InitVegMap()
*****************************************************************************/
void InitVegMap(OPTIONSTRUCT * Options, LISTPTR Input, MAPSIZE * Map, VEGPIX *** VegMap)
{
  const char *Routine = "InitVegMap";
  char VarName[BUFSIZE + 1];
  char VegMapFileName[BUFSIZE + 1];
  int i;			/* counter */
  int x;			/* counter */
  int y;			/* counter */
  int flag;
  int NumberType;		/* number type */
  unsigned char *Type;		/* Vegetation type */

  /* Get the map filename from the [VEGETATION] section */
  GetInitString("VEGETATION", "VEGETATION MAP FILE", "", VegMapFileName,
    (unsigned long)BUFSIZE, Input);
  if (!VegMapFileName)
    ReportError("VEGETATION MAP FILE", 51);

  /* Read the vegetation type */
  GetVarName(005, 0, VarName);
  GetVarNumberType(005, &NumberType);
  if (!(Type = (unsigned char *)calloc(Map->NX * Map->NY,
    SizeOfNumberType(NumberType))))
    ReportError((char *)Routine, 1);
  flag = Read2DMatrix(VegMapFileName, Type, NumberType, Map, 0, VarName, 0);

  /* Assign the attributes to the correct map pixel */
  if (!(*VegMap = (VEGPIX **)calloc(Map->NY, sizeof(VEGPIX *))))
    ReportError((char *)Routine, 1);
  for (y = 0; y < Map->NY; y++) {
    if (!((*VegMap)[y] = (VEGPIX *)calloc(Map->NX, sizeof(VEGPIX))))
      ReportError((char *)Routine, 1);
  }

  if ((Options->FileFormat == NETCDF && flag == 0)
    || (Options->FileFormat == BIN))
  {
    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*VegMap)[y][x].Veg = Type[i];
        (*VegMap)[y][x].Tcanopy = 0.0;
      }
    }
  }
  else if (Options->FileFormat == NETCDF && flag == 1) {
    for (y = Map->NY - 1, i = 0; y >= 0; y--) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*VegMap)[y][x].Veg = Type[i];
        (*VegMap)[y][x].Tcanopy = 0.0;
      }
    }
  }
  else ReportError((char *)Routine, 57);

  free(Type);
}


/*****************************************************************************
InitCanopyGapMap()
*****************************************************************************/
void InitCanopyGapMap(OPTIONSTRUCT *Options, LISTPTR Input, MAPSIZE *Map,
  LAYER *Soil, LAYER *Veg, VEGTABLE *VType, VEGPIX ***VegMap,
  SOILTABLE *SType, SOILPIX ***SoilMap)
{
  const char *Routine = "InitCanopyGapMap";
  char VarName[BUFSIZE + 1];
  char CanopyMapFileName[BUFSIZE + 1];
  int i, j;			/* counter */
  int x;			/* counter */
  int y;			/* counter */
  int flag;
  int NVeg;
  int NSoil;
  int NumberType;		/* number type */
  float *Gap;		/* gap diameter */

  /* Get the canopy gap map filename from the [VEGETATION] section */
  GetInitString("VEGETATION", "CANOPY GAP MAP FILE", "", CanopyMapFileName,
    (unsigned long)BUFSIZE, Input);
  if (!CanopyMapFileName)
    ReportError("CANOPY GAP MAP FILE", 51);

  /* Read the vegetation type */
  GetVarName(007, 0, VarName);
  GetVarNumberType(007, &NumberType);
  if (!(Gap = (float *)calloc(Map->NX * Map->NY,
    SizeOfNumberType(NumberType))))
    ReportError((char *)Routine, 1);
  flag = Read2DMatrix(CanopyMapFileName, Gap, NumberType, Map, 0, VarName, 0);

  /* if NetCDF, may need to reverse the matrix */
  if ((Options->FileFormat == NETCDF && flag == 0)
    || (Options->FileFormat == BIN))
  {
    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*VegMap)[y][x].Gapping = Gap[i];
        /* set gapping to false for cells with no overstory */
        if (VType[(*VegMap)[y][x].Veg - 1].OverStory == FALSE)
          (*VegMap)[y][x].Gapping = 0.0;
      }
    }
  }
  else if (Options->FileFormat == NETCDF && flag == 1) {
    for (y = Map->NY - 1, i = 0; y >= 0; y--) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*VegMap)[y][x].Gapping = Gap[i];
        /* set gapping to false for cells with no overstory */
        if (VType[(*VegMap)[y][x].Veg - 1].OverStory == FALSE)
          (*VegMap)[y][x].Gapping = 0.0;
        /* set gapping to false given glacier cell */
        if (VType[(*VegMap)[y][x].Veg - 1].Index == GLACIER)
          (*VegMap)[y][x].Gapping = 0.0;
      }
    }
  }
  else ReportError((char *)Routine, 57);

  for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
          printf("Gap Size %f \n",(*VegMap)[y][x].Gapping);
      NVeg = Veg->MaxLayers;
      NSoil = Soil->MaxLayers;
      if (Options->CanopyGapping) {
        if (!((*VegMap)[y][x].Type = (CanopyGapStruct *)calloc(2, sizeof(CanopyGapStruct))))
          ReportError((char *)Routine, 1);
        for (i = 0; i < CELL_PARTITION; i++) {
          if (!((*VegMap)[y][x].Type[i].IntRain = (float *)calloc(NVeg, sizeof(float))))
            ReportError((char *)Routine, 1);
          if (!((*VegMap)[y][x].Type[i].IntSnow = (float *)calloc(NVeg, sizeof(float))))
            ReportError((char *)Routine, 1);
          if (!((*VegMap)[y][x].Type[i].Moist = (float *)calloc(NSoil+1, sizeof(float))))
            ReportError((char *)Routine, 1);
          if (!((*VegMap)[y][x].Type[i].EPot = (float *)calloc(NVeg+1, sizeof(float))))
            ReportError((char *)Routine, 1);
          if (!((*VegMap)[y][x].Type[i].EAct = (float *)calloc(NVeg+1, sizeof(float))))
            ReportError((char *)Routine, 1);
          if (!((*VegMap)[y][x].Type[i].EInt = (float *)calloc(NVeg, sizeof(float))))
            ReportError((char *)Routine, 1);
          if (!((*VegMap)[y][x].Type[i].ESoil = (float **)calloc(NVeg, sizeof(float *))))
            ReportError((char *)Routine, 1);

          for (j = 0; j < NVeg; j++) {
            if (!((*VegMap)[y][x].Type[i].ESoil[j] = (float *)calloc(NSoil, sizeof(float))))
              ReportError((char *)Routine, 1);
          }
        }
      }
    }
  }
  free(Gap);
}




/*****************************************************************************
InitTileMap()
*****************************************************************************/
void InitTileMap(OPTIONSTRUCT *Options, LISTPTR Input, MAPSIZE *Map,
  LAYER *Soil, LAYER *Veg, VEGTABLE *VType, VEGPIX ***VegMap,
  SOILTABLE *SType, SOILPIX ***SoilMap)
{
  const char *Routine = "InitTileMap";
  char VarName[BUFSIZE + 1];
  char NorthFacingFileName[BUFSIZE + 1];
  char SouthFacingFileName[BUFSIZE + 1];
  char ExposedFileName[BUFSIZE + 1];
  char ForestTileFileName[BUFSIZE + 1];
  int i, j;			/* counter */
  int x;			/* counter */
  int y;			/* counter */
  int flag;
  int NVeg;
  int NSoil;
  int NumberType;		/* number type */
  float *NFfrac;		/* fraction of grid cell that is North facing forest edge */
  float *SFfrac;		/* fraction of grid cell that is South facing forest edge */
  float *EXPfrac;		/* fraction of grid cell that is Exposed  */
  float *FORfrac;		/* fraction of grid cell that is Forest */

  /* Get the canopy gap map filename from the [VEGETATION] section */
  GetInitString("VEGETATION", "NORTH FACING EDGE MAP FILE", "", NorthFacingFileName,
    (unsigned long)BUFSIZE, Input);
  GetInitString("VEGETATION", "SOUTH FACING EDGE MAP FILE", "", SouthFacingFileName,
    (unsigned long)BUFSIZE, Input);
  GetInitString("VEGETATION", "EXPOSED TILE MAP FILE", "", ExposedFileName,
    (unsigned long)BUFSIZE, Input);
  GetInitString("VEGETATION", "FOREST TILE MAP FILE", "", ForestTileFileName,
    (unsigned long)BUFSIZE, Input);

  if (!NorthFacingFileName)
    ReportError("NORTH FACING EDGE MAP FILE", 51);
  if (!SouthFacingFileName)
    ReportError("SOUTH FACING EDGE MAP FILE", 51);
  if (!ExposedFileName)
    ReportError("EXPOSED TILE MAP FILE", 51);
  if (!ForestTileFileName)
    ReportError("FOREST TILE MAP FILE", 51);

  /* Read the vegetation type */
  GetVarName(010, 0, VarName);
  GetVarNumberType(010, &NumberType);

  if (!(NFfrac = (float *)calloc(Map->NX * Map->NY,
    SizeOfNumberType(NumberType))))
    ReportError((char *)Routine, 1);
  flag = Read2DMatrix(NorthFacingFileName, NFfrac, NumberType, Map, 0, VarName, 0);
  printf("Read in North Facing Map named %s \n", VarName);

  GetVarName(011, 0, VarName);
  GetVarNumberType(011, &NumberType);

  if (!(SFfrac = (float *)calloc(Map->NX * Map->NY,
    SizeOfNumberType(NumberType))))
    ReportError((char *)Routine, 1);
  flag = Read2DMatrix(SouthFacingFileName, SFfrac, NumberType, Map, 0, VarName, 0);
  printf("Read in South Facing Map named %s \n", VarName);

  GetVarName(012, 0, VarName);
  GetVarNumberType(012, &NumberType);

  if (!(EXPfrac = (float *)calloc(Map->NX * Map->NY,
    SizeOfNumberType(NumberType))))
    ReportError((char *)Routine, 1);
  flag = Read2DMatrix(ExposedFileName, EXPfrac, NumberType, Map, 0, VarName, 0);
  printf("Read in Exposed Map named %s \n", VarName);

  GetVarName(013, 0, VarName);
  GetVarNumberType(013, &NumberType);

  if (!(FORfrac = (float *)calloc(Map->NX * Map->NY,
    SizeOfNumberType(NumberType))))
    ReportError((char *)Routine, 1);
  flag = Read2DMatrix(ForestTileFileName, FORfrac, NumberType, Map, 0, VarName, 0);
  printf("Read in Forest Map named %s \n", VarName);

  /* if NetCDF, may need to reverse the matrix */
  if ((Options->FileFormat == NETCDF && flag == 0)
    || (Options->FileFormat == BIN))
  {
    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*VegMap)[y][x].NFfrac   = NFfrac[i];
        (*VegMap)[y][x].SFfrac   = SFfrac[i];
        (*VegMap)[y][x].EXPfrac  = EXPfrac[i];
        (*VegMap)[y][x].FORfrac  = FORfrac[i];
        /* set All Tile Fractions to false for cells with no overstory */
      }
    }
  }
  else if (Options->FileFormat == NETCDF && flag == 1) {
    for (y = Map->NY - 1, i = 0; y >= 0; y--) {
      for (x = 0; x < Map->NX; x++, i++) {
        (*VegMap)[y][x].NFfrac   = NFfrac[i];
        (*VegMap)[y][x].SFfrac   = SFfrac[i];
        (*VegMap)[y][x].EXPfrac  = EXPfrac[i];
        (*VegMap)[y][x].FORfrac  = FORfrac[i];
        /* set gapping to false for cells with no overstory */
        if (VType[(*VegMap)[y][x].Veg - 1].OverStory == FALSE)
          (*VegMap)[y][x].NFfrac  = 0.0;
          (*VegMap)[y][x].SFfrac  = 0.0;
          (*VegMap)[y][x].EXPfrac = 0.0;
          (*VegMap)[y][x].FORfrac = 0.0;
        /* set gapping to false given glacier cell */
        if (VType[(*VegMap)[y][x].Veg - 1].Index == GLACIER)
          (*VegMap)[y][x].NFfrac  = 0.0;
          (*VegMap)[y][x].SFfrac  = 0.0;
          (*VegMap)[y][x].EXPfrac = 0.0;
          (*VegMap)[y][x].FORfrac = 0.0;

      }
    }
  }
  else ReportError((char *)Routine, 57);

    for (y = 0, i = 0; y < Map->NY; y++) {
      for (x = 0; x < Map->NX; x++, i++) {
        printf("NFfrac %f , SFfrac %f , EXPfrac %f , FORfrac %f \n",(*VegMap)[y][x].NFfrac,(*VegMap)[y][x].SFfrac,(*VegMap)[y][x].EXPfrac,(*VegMap)[y][x].FORfrac);
      }
    }
   for (y = 0; y < Map->NY; y++) {
    for (x = 0; x < Map->NX; x++) {
      NVeg = Veg->MaxLayers;
      NSoil = Soil->MaxLayers;
      if (Options->CanopyTiling) {
        if (!((*VegMap)[y][x].Tile = (TileStruct *)calloc(4, sizeof(TileStruct))))
          ReportError((char *)Routine, 1);
        for (i = 0; i < TILE_PARTITION; i++) {
          if (!((*VegMap)[y][x].Tile[i].IntRain = (float *)calloc(NVeg, sizeof(float))))
            ReportError((char *)Routine, 1);
          if (!((*VegMap)[y][x].Tile[i].IntSnow = (float *)calloc(NVeg, sizeof(float))))
            ReportError((char *)Routine, 1);
          if (!((*VegMap)[y][x].Tile[i].Moist = (float *)calloc(NSoil+1, sizeof(float))))
            ReportError((char *)Routine, 1);
          if (!((*VegMap)[y][x].Tile[i].EPot = (float *)calloc(NVeg+1, sizeof(float))))
            ReportError((char *)Routine, 1);
          if (!((*VegMap)[y][x].Tile[i].EAct = (float *)calloc(NVeg+1, sizeof(float))))
            ReportError((char *)Routine, 1);
          if (!((*VegMap)[y][x].Tile[i].EInt = (float *)calloc(NVeg, sizeof(float))))
            ReportError((char *)Routine, 1);
          if (!((*VegMap)[y][x].Tile[i].ESoil = (float **)calloc(NVeg, sizeof(float *))))
            ReportError((char *)Routine, 1);

          for (j = 0; j < NVeg; j++) {
            if (!((*VegMap)[y][x].Tile[i].ESoil[j] = (float *)calloc(NSoil, sizeof(float))))
              ReportError((char *)Routine, 1);
          }
        }
      }
    }
  }
  free(NFfrac);
  free(SFfrac);
  free(EXPfrac);
  free(FORfrac);
}



