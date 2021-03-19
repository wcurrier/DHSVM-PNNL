/*
* SUMMARY:      CanopyGapRadiation.c
* USAGE:        Part of DHSVM
*
* AUTHOR:       Wiliam Ryan Currier
* E-MAIL:       currierw@uw.edu
* ORIG-DATE:    Jul-19
* DESCRIPTION:  Calculate radiation balance for different tiles in relation to the forest
                Note: heavily adapted from CanopyGapRadiation.c and RadiationBalance.c
* DESCRIP-END.
* FUNCTIONS:    TileRadiation()
* Reference:


 * Reference:

   Wigmosta, M. S., L. W. Vail, and D. P. Lettenmaier, A distributed
   hydrology-vegetation model for complex terrain, Water Resour. Res.,
   30(6), 1665-1679, 1994.

   Nijssen and Lettenmaier, A simplified approach for predicting shortwave
   radiation transfer through boreal forest canopies, JGR, 1999.

   **references of improved radiation scheme **

   Reifsnyder, W. E., and H. W. Lull (1965), Radiant energy in relation to forests,
   USDA For. Serv. Tech. Bull., 1344, 111.

   Thyer M. et al. (2004), Diagnosing a distributed hydrologic model for two
   high-elevation forested catchments based on detailed stand- and basin-scale data.
   Water Resources Research. DOI:10.1029/2003WR002414.

    C.R. Ellis, J.W. Pomeroy, and T.E. Link, Modeling increases in snowmelt
    yield and desynchronization resulting from forest gap-thinning treatments
    in a northern mountain headwater basin, Water Resour. Res., 49, 936-949, 2013.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "functions.h"
#include "data.h"
#include "DHSVMerror.h"
#include "massenergy.h"
#include "constants.h"

#define MAXIT 100

/********************************************************************************
Function Name: TileRadiation()

Purpose      :

Required     :

     int HeatFluxOption - TRUE if surface temperature is being calculated
     float Rs           - Incoming shortwave radiation adjusted by topo shading
                          if shading option is on (W/m2)
     float Rsb          - Incoming direct shortwave radiation separated from Rs
     float Rsd          - Incoming diffuse shortwave radiation separated from Rs
     float Ld           - Incoming longwave radiation (W/m2)
     float Tair         - Ambient air temperature (C)
     float Tcanopy      - Canopy temperature from the previous timestep (C)
     float Tsoil        - soil surface temperature from the previous timestep (C)
     VEGTABLE VType     - Information about number of veg layers
     TileStruct **Tile  - Information about veg/met/snow conditions at current pixel tile
     Extn/Tau           - light extinction coeff.
     h                  - canopy height

Returns      :
Comments     :
********************************************************************************/
void TileRadiation(TileStruct **Tile, OPTIONSTRUCT *Options, float SineSolarAltitude, float Rs,
  float Rsb, float Rsd, int CanopyRadAttOption, VEGTABLE *VType, int HeatFluxOption,
  float Ld, float Tcanopy, float Tair, float Tsoil, float SoilAlbedo)

{
  /* ==== Setup the canopy conditions and some caveats for whether  ====
     ==== the longwave is enhanced by forest of forest is shaded.   ====
     ==== This happens in SouthFacing and NorthFacing flags         ====*/

  (*Tile)[NorthFacing].OverStory      = VType->OverStory;
  (*Tile)[NorthFacing].UnderStory     = VType->UnderStory;
  (*Tile)[NorthFacing].NorthFacingInt = TRUE;
  (*Tile)[NorthFacing].SouthFacingInt = FALSE;
  (*Tile)[NorthFacing].ExposedInt     = FALSE;
  (*Tile)[NorthFacing].ForestInt      = FALSE;

  (*Tile)[SouthFacing].OverStory      = VType->OverStory;
  (*Tile)[SouthFacing].UnderStory     = VType->UnderStory;
  (*Tile)[SouthFacing].NorthFacingInt = FALSE;
  (*Tile)[SouthFacing].SouthFacingInt = TRUE;
  (*Tile)[SouthFacing].ExposedInt     = FALSE;
  (*Tile)[SouthFacing].ForestInt      = FALSE;

  (*Tile)[Exposed].OverStory          = VType->OverStory;
  (*Tile)[Exposed].UnderStory         = VType->UnderStory;
  (*Tile)[Exposed].NorthFacingInt     = FALSE;
  (*Tile)[Exposed].SouthFacingInt     = FALSE;
  (*Tile)[Exposed].ExposedInt         = TRUE;
  (*Tile)[Exposed].ForestInt          = FALSE;

  (*Tile)[ForestTile].OverStory       = VType->OverStory;
  (*Tile)[ForestTile].UnderStory      = VType->UnderStory;
  (*Tile)[ForestTile].NorthFacingInt  = FALSE;
  (*Tile)[ForestTile].SouthFacingInt  = FALSE;
  (*Tile)[ForestTile].ExposedInt      = FALSE;
  (*Tile)[ForestTile].ForestInt       = TRUE;

}

/********************************************************************************
Function Name: TileShortLongRadiation()
********************************************************************************/
void TileShortRadiation(VEGTABLE *VType, TileStruct *Tile, OPTIONSTRUCT *Options,
                        float SineSolarAltitude, float Rs, float Rsb, float Rsd,
                        int CanopyRadAttOption, float SoilAlbedo)
{

  float F;			        /* Fraction of pixel covered by top canopy layer [0-1] */
  float h;                  /* Canopy height (m) */
  float Albedo[2];		    /* Albedo of each layer */
  float Tau;			    /* Transmittance for overstory vegetation layer */
  float Taub, Taud;         /* Transmittance for overstory vegetation layer for
                           direct and diffuse radiation, respectively */
  float TaudNF, TaudSF;     /* Diffused transmittance for north facing and south
                               south-facing edges*/
  float TauMult;            /* North-facing edge direct beam radiation multiplier [-]*/
  float TaubNF;				/* Final Radiaiton transmittance for north-facing edge [-]*/

  unsigned char UnderStory;
  unsigned char OverStory;

  OverStory  = Tile->OverStory;
  UnderStory = Tile->UnderStory;

  Taud   = VType->Taud;
  TaudNF = VType->TaudNF;
  TaudSF = VType->TaudSF;
  TauMult =  VType->TauMult;

  F = VType->Fract[0];
  h = VType->Height[0];
  
  /* Determine Albedo */
  Albedo[0] = VType->Albedo[0];
  /* With snow, understory canopy albedo is set equal to snow albedo */
  if (Tile->HasSnow == TRUE)
    Albedo[1] = Tile->Albedo;
  else if (VType->UnderStory == TRUE)
    Albedo[1] = VType->Albedo[1];
  else
    Albedo[1] = SoilAlbedo;

  /* Improved radiation scheme taking into account solar position */
  if (OverStory == TRUE) {
    if (Options->ImprovRadiation) {
      if (SineSolarAltitude > 0. && Rs > 0.) {
        Tau = exp(-VType->ExtnCoeff * h / SineSolarAltitude);
      }
      else
        Tau = 0.;
    }
    else if (CanopyRadAttOption == FIXED) { /* conventional radiation scheme */
      Tau = exp(-VType->Atten * VType->LAI[0]);
    }
    /* Nijssen's simplified radiation scheme as in Nijssen and Lettenmaier, 1999 */
    else if (CanopyRadAttOption == VARIABLE) {
      /* Calculate transmittance of overstory canopy for direct radiation:
         1) LAI * ClumpingFactor = Effective LAI
         2) Formulation is typically based on the cos of the solar zenith angle,
         which is the sin of the solar altitude (SA = 90 - SZA) */
      Taub = exp(-VType->LAI[0] / VType->ClumpingFactor *
        (VType->LeafAngleA / SineSolarAltitude + VType->LeafAngleB));

      /* transmittance for diffuse radiation (cacluated in CheckOut.c as a function of
         LeafAngleA and LeafAngleB and solar altitude) */
      Taud = VType->Taud;

      /* cacluate the total canopy transimittance for shortwave radiation (adjusted to
         scattering and multiple reflection */
      if (Rs > 0.0) {
        Tau = Taub * Rsb / Rs + Taud * Rsd / Rs;
        /* adjust Tau to scaterring parameter */
        Tau = pow(Tau, (VType->Scat));
        /* adjust Tau to over- and under- story reflection */
        Tau = Tau / (1 - Albedo[0] * Albedo[1]);
      }
      else
        Tau = 0.;
    }
  }


  /* North Facing Shortwave Radiation */
  if (Tile->NorthFacingInt == TRUE) {
	// Tile->NetShort[0] = Rs * F * ((1-Albedo[0])-Tau*(1-Albedo[1]));          /* Overstory = 0 */
	Tile->NetShort[0] = 0.;          /* Overstory = 0 */

    if (Options->ImprovRadiation == TRUE) {
      TaubNF = Tau + (TauMult*Tau);
      if (TaubNF > 1) {
		TaubNF=1;
	  }
      printf("TauMult = %f \n", TauMult);
      printf("TaubNF = %f \n", TaubNF);
	  Tile->NetShort[1] = (1-Albedo[1])*(Rs*(1-F) + F*(Rsb*TaubNF + Rsd*TaudNF)); /* Understory = 1 */
    }
    else {
      Tile->NetShort[1] = Rs * (1-Albedo[1]) * ((1-F)+(Tau*F));
    }
  }
  /* South Facing Shortwave Radiation */
  if (Tile->SouthFacingInt == TRUE) {
      Tile->NetShort[1] = (Rsb + Rsd*TaudSF) * (1-Albedo[1]); /* Understory = 1  */
      Tile->NetShort[0] = 0.;                 /* Overstory = 0 */
  }
  /* Exposed Shortwave Radiation */
  if (Tile->ExposedInt == TRUE) {
      Tile->NetShort[1] = Rs * (1-Albedo[1]); /* Understory = 1  */
      Tile->NetShort[0] = 0.;                 /* Overstory = 0 */
  }
  /* Forest Radiation */
  if (Tile->ForestInt == TRUE) {
	Tile->NetShort[0] = Rs * F * ((1-Albedo[0])-Tau*(1-Albedo[1]));          /* Overstory = 0 */

    if (Options->ImprovRadiation == TRUE) {
	  Tile->NetShort[1] = (1-Albedo[1])*(Rs*(1-F) + F*(Rsb*Tau + Rsd*Taud)); /* Understory = 1 */
    }
    else {
      Tile->NetShort[1] = Rs * (1-Albedo[1]) * ((1-F)+(Tau*F));
    }
  }

}

/********************************************************************************
Function Name: TileLongRadiation()
Purpose      : Calculate net long radiation incident on canopy gap
********************************************************************************/
void TileLongRadiation(VEGTABLE *VType,TileStruct *Tile, OPTIONSTRUCT *Options,
                        int CanopyRadAttOption,float Ld, int HeatFluxOption,
                        float Tair, float Tcanopy, float Tsoil, float SoilAlbedo)
{
  float F;			        /* Fraction of pixel covered by top canopy layer [0-1] */
  float Fnf;                /* Fraction of the north-facing edge within the pixel that is exposed to the forest */
  float Fsf;                /* Fraction of the south-facing edge within the pixel that is exposed to the forest */
  float h;                  /* Canopy height (m) */
  float Tsurf;			    /* Surface temperature (C) */
  float Vf;
  double Tmp;               /* used in longwave calulatnions */
  unsigned char OverStory;

  F = VType->Fract[0];  /* Fraction of pixel covered by top canopy layer [0-1] */
  h = VType->Height[0]; /* Canopy height (m) */

  Fnf = VType->FractNF[0];
  Fsf = VType->FractSF[0];

  if (CanopyRadAttOption == VARIABLE) {
    F = VType->HemiFract[0];
  }

  OverStory  = Tile->OverStory;

  Vf=VType->Vf;

  if (Tile->HasSnow == TRUE)
    Tsurf = Tile->TSurf;
  else if (HeatFluxOption == TRUE)
    Tsurf = Tsoil;
  else
    Tsurf = Tair;


  /* North Facing Longwave Radiation */
  if (Tile->NorthFacingInt == TRUE) {
      Tmp = Tsurf + 273.15;
      Tile->LongOut[1] = STEFAN * (Tmp * Tmp * Tmp * Tmp); /* Understory = 1  */
      Tmp = Tcanopy + 273.15;
      Tile->LongOut[0] = STEFAN * (Tmp * Tmp * Tmp * Tmp); /* Overstory = 0 */
	  // Tile->LongOut[0] = 0.;

      Tile->LongIn[1] = Ld * (1-Fnf) + Tile->LongOut[0] * Fnf; /* Understory = 1 - could use some work  */
      // Tile->LongIn[0] = (Ld + Tile->LongOut[1]) * Fnf;       /* Overstory = 0 */
	  Tile->LongIn[0] = 0.;
  }
  /* South Facing Longwave Radiation */
  if (Tile->SouthFacingInt == TRUE) {
    Tmp = Tsurf + 273.15;
    Tile->LongOut[1] = STEFAN * (Tmp * Tmp * Tmp * Tmp); /* Understory = 1  */
    Tmp = Tcanopy + 273.15;
    Tile->LongOut[0] = STEFAN * (Tmp * Tmp * Tmp * Tmp); /* Overstory = 0 */
	// Tile->LongOut[0] = 0.;

    Tile->LongIn[1] = Ld * (1-Fsf) + Tile->LongOut[0] * Fsf; /* Understory = 1 - could use some work  */
    // Tile->LongIn[0] = (Ld + Tile->LongOut[1]) * Fsf;       /* Overstory = 0 */
	Tile->LongIn[0] = 0.;
  }
  /* Exposed Longwave Radiation */
  if (Tile->ExposedInt == TRUE) {
      Tmp = Tsurf + 273.15;
      Tile->LongOut[1] = STEFAN * (Tmp * Tmp * Tmp * Tmp); /* Understory = 1  */
      Tile->LongOut[0] = 0.;                               /* Overstory = 0 */

      Tile->LongIn[1]  = Ld;                               /* Understory = 1  */
      Tile->LongIn[0]  = 0.;                               /* Overstory = 0 */
  }
  /* Forest Radiation */
  if (Tile->ForestInt == TRUE) {
    Tmp = Tsurf + 273.15;
    Tile->LongOut[1] = STEFAN * (Tmp * Tmp * Tmp * Tmp); /* Understory = 1  */
    Tmp = Tcanopy + 273.15;
    Tile->LongOut[0] = STEFAN * (Tmp * Tmp * Tmp * Tmp); /* Overstory = 0 */

    if (Options->ImprovRadiation == TRUE) {
        Tile->LongIn[1] = Ld * (1-Vf) + Tile->LongOut[0] * Vf; /* Understory = 1  */
        Tile->LongIn[0] = (Ld + Tile->LongOut[1]) * Vf;       /* Overstory = 0 */
    }
    else {
        Tile->LongIn[1] = Ld * (1-F) + Tile->LongOut[0] * F; /* Understory = 1  */
        Tile->LongIn[0] = (Ld + Tile->LongOut[1]) * F;       /* Overstory = 0 */
    }
  }

}

