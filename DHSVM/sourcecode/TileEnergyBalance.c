/*
* SUMMARY:      TileEnergyBalance.c
* USAGE:        Part of DHSVM
*
* AUTHOR:       William Ryan Currier
* E-MAIL:       currierw@uw.edu
* ORIG-DATE:    Aug-19
* DESCRIPTION:  Calculate snow balance for NorthFacing/SouthFacing/Exposed/Forest
* DESCRIP-END.
* FUNCTIONS:
* Reference:
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "data.h"
#include "DHSVMerror.h"
#include "functions.h"
#include "massenergy.h"
#include "constants.h"
#include "snow.h"

/*****************************************************************************
Function name:CalcNoOverStoryAerodynamic()

Purpose      : Calculate the aerodynamic resistance for each vegetation
layer, and the wind 2m above the layer boundary.
*****************************************************************************/
void CalcNoOverStoryAerodynamic(TileStruct **Tile, int NVegLayers,
  float *Height)

{
  float K2;
  float Z0_Lower;    /* roughness length for understory (m) */
  float d_Lower;     /* displacement height (m) */

  K2 = VON_KARMAN * VON_KARMAN;


  /* Exposed */

  /* bare soil */
  if ((*Tile)[Exposed].UnderStory == FALSE) {
    Z0_Lower = Z0_GROUND;
    d_Lower = 0;
  }
  /* understory present*/
  else {
    Z0_Lower = Z0_MULTIPLIER * Height[1]; /* Understory Height */
    d_Lower = D0_MULTIPLIER * Height[1];  /* Understory Height */
  }

  /* No snow: get wind speed & aerodynamic resistence value */
  (*Tile)[Exposed].U[1] =
    log((2.+Z0_Lower)/Z0_Lower) / log((Zref-d_Lower)/Z0_Lower);
  (*Tile)[Exposed].Ra[1] =
    log((2.+Z0_Lower)/Z0_Lower) * log((Zref-d_Lower)/Z0_Lower)/K2;

  /* get the wind speed and aerodynamic resistence value for snow surface*/
  (*Tile)[Exposed].USnow = log((2.+Z0_SNOW)/Z0_SNOW) / log(Zref/Z0_SNOW);
  (*Tile)[Exposed].RaSnow = log((2.+Z0_SNOW)/Z0_SNOW) * log(Zref/Z0_SNOW)/K2;


  /* North Facing */

  /* bare soil */
  if ((*Tile)[NorthFacing].UnderStory == FALSE) {
    Z0_Lower = Z0_GROUND;
    d_Lower = 0;
  }
  /* understory present*/
  else {
    Z0_Lower = Z0_MULTIPLIER * Height[1]; /* Understory Height */
    d_Lower = D0_MULTIPLIER * Height[1];  /* Understory Height */
  }
  /* No snow: get wind speed & aerodynamic resistence value */
  (*Tile)[NorthFacing].U[1] =
    log((2.+Z0_Lower)/Z0_Lower) / log((Zref-d_Lower)/Z0_Lower);
  (*Tile)[NorthFacing].Ra[1] =
    log((2.+Z0_Lower)/Z0_Lower) * log((Zref-d_Lower)/Z0_Lower)/K2;

  /* get the wind speed and aerodynamic resistence value for snow surface*/
  (*Tile)[NorthFacing].USnow = log((2.+Z0_SNOW)/Z0_SNOW) / log(Zref/Z0_SNOW);
  (*Tile)[NorthFacing].RaSnow = log((2.+Z0_SNOW)/Z0_SNOW) * log(Zref/Z0_SNOW)/K2;


  /* South Facing */

  /* bare soil */
  if ((*Tile)[SouthFacing].UnderStory == FALSE) {
    Z0_Lower = Z0_GROUND;
    d_Lower = 0;
  }
  /* understory present*/
  else {
    Z0_Lower = Z0_MULTIPLIER * Height[1];  /* Understory Height */
    d_Lower = D0_MULTIPLIER * Height[1];   /* Understory Height */
  }

  /* No snow: get wind speed & aerodynamic resistence value */
  (*Tile)[SouthFacing].U[1] =
    log((2.+Z0_Lower)/Z0_Lower) / log((Zref-d_Lower)/Z0_Lower);
  (*Tile)[SouthFacing].Ra[1] =
    log((2.+Z0_Lower)/Z0_Lower) * log((Zref-d_Lower)/Z0_Lower)/K2;

  /* get the wind speed and aerodynamic resistence value for snow surface*/
  (*Tile)[SouthFacing].USnow = log((2.+Z0_SNOW)/Z0_SNOW) / log(Zref/Z0_SNOW);
  (*Tile)[SouthFacing].RaSnow = log((2.+Z0_SNOW)/Z0_SNOW) * log(Zref/Z0_SNOW)/K2;

}

/********************************************************************************
Function Name: TileNoOverStoryInterception()

Purpose      : Calculate snow/rain interception
Returns      :
Comments     :
********************************************************************************/
void TileNoOverStoryInterception(OPTIONSTRUCT *Options, TileStruct **Tile,
  int HeatFluxOption, int y, int x, int Dt, int NVegLActual,
  float DX, float DY, float UpperRa, float UpperWind, VEGTABLE *VType,
  SOILPIX *LocalSoil, VEGPIX *LocalVeg, SNOWPIX *LocalSnow,
  PRECIPPIX *LocalPrecip, PIXRAD *LocalRad, PIXMET *LocalMet) {

  float Tsurf;
  float SnowLongIn;			/* Incoming longwave radiation at snow surface (W/m2) */
  float SnowNetShort;		/* Net amount of short wave radiation at the snow surface (W/m2) */
  float SnowRa;				/* Aerodynamic resistance for snow */
  float SnowWind;		    /* Wind 2 m above snow */

  if ((*Tile)[Exposed].UnderStory == TRUE) {
    (*Tile)[Exposed].Tcanopy = LocalMet->Tair;
    (*Tile)[Exposed].CanopyVaporMassFlux = 0.0;
    (*Tile)[Exposed].TempIntStorage = 0.0;

    /* calculate rain interception in the exposed area*/
    CanopyGapInterceptionStorage((*Tile)[Exposed].NVegLActual, VType->MaxInt,
      VType->Fract, (*Tile)[Exposed].IntRain, &((*Tile)[Exposed].RainFall));
  }

  if ((*Tile)[NorthFacing].UnderStory == TRUE) {
    (*Tile)[NorthFacing].Tcanopy = LocalMet->Tair;
    (*Tile)[NorthFacing].CanopyVaporMassFlux = 0.0;
    (*Tile)[NorthFacing].TempIntStorage = 0.0;

    /* calculate rain interception in the north facing area*/
    CanopyGapInterceptionStorage((*Tile)[NorthFacing].NVegLActual, VType->MaxInt,
      VType->Fract, (*Tile)[NorthFacing].IntRain, &((*Tile)[NorthFacing].RainFall));
  }

  if ((*Tile)[SouthFacing].UnderStory == TRUE) {
    (*Tile)[SouthFacing].Tcanopy = LocalMet->Tair;
    (*Tile)[SouthFacing].CanopyVaporMassFlux = 0.0;
    (*Tile)[SouthFacing].TempIntStorage = 0.0;

    /* calculate rain interception in the south facing area*/
    CanopyGapInterceptionStorage((*Tile)[SouthFacing].NVegLActual, VType->MaxInt,
      VType->Fract, (*Tile)[SouthFacing].IntRain, &((*Tile)[SouthFacing].RainFall));

  }

}

/********************************************************************************
Function Name: NoOverstorySnowMelt()

Purpose      : Calculate snow melt
Returns      :
Comments     :
********************************************************************************/
void NoOverStorySnowMelt(OPTIONSTRUCT *Options, int y, int x, int Dt,
  TileStruct **Tile, float DX, float DY, VEGTABLE *VType, SOILTABLE *SType,
  VEGPIX *LocalVeg, SNOWPIX *LocalSnow, PRECIPPIX *LocalPrecip, PIXRAD *LocalRad,
  PIXMET *LocalMet, SOILPIX *LocalSoil, int HeatFluxOption, int CanopyRadAttOption)
{
  float SnowLongIn;			/* Incoming longwave radiation at snow surface (W/m2) */
  float SnowNetShort;		/* Net amount of short wave radiation at the snow surface (W/m2) */
  float SnowRa;				/* Aerodynamic resistance for snow */
  float SnowWind;		    /* Wind 2 m above snow */
  float Tsurf;              /* Surface temperature */
  float OldTSurf;           /* Effective surface temperature at the end of the last timestep (C)*/
  float Tmean;              /* Average snow surface temperature*/
  float tmp;                /* temporary variable */
  float Tmp;                /* temporary variable */
  float Ls;			        /* Latent heat of sublimation (J/kg) */
  float GAPWIND_FACTOR_Tile; /* Local GAPWIND Factor [-] */

  GAPWIND_FACTOR_Tile = 1 - LocalVeg->FORfrac;

  /********************* calculate NorthFacing snow melt *********************/

  /* NorthFacing */

  if ((*Tile)[NorthFacing].HasSnow || (*Tile)[NorthFacing].SnowFall > 0.0) {
    SnowLongIn = (*Tile)[NorthFacing].LongIn[1];
    SnowNetShort = (*Tile)[NorthFacing].NetShort[1];

//	printf("(*Tile)[NorthFacing].USnow %f \n", (*Tile)[NorthFacing].USnow);
//	printf("(*Tile)[SouthFacing].USnow %f \n", (*Tile)[SouthFacing].USnow);
//	printf("(*Tile)[Exposed].USnow %f \n", (*Tile)[Exposed].USnow);
//	printf("VType->USnow %f \n", VType->USnow);
//	printf("VType->USnow %f \n", VType->USnow);


	SnowWind = (*Tile)[NorthFacing].USnow * LocalMet->Wind;
    SnowRa = (*Tile)[NorthFacing].RaSnow / LocalMet->Wind;

//	printf("North Facing SnowWind Before Correction: %f \n", SnowWind);
//	printf("North Facing SnowRa Before Correction: %f \n", SnowRa);

    /* adjust the wind and Ra values for gap so they fall between open
    and forested values */

    tmp = VType->USnow*LocalMet->Wind; //forested snow wind
	SnowWind = tmp + (SnowWind- tmp)*GAPWIND_FACTOR_Tile;

	tmp = VType->RaSnow/LocalMet->Wind;
    SnowRa = tmp - (tmp-SnowRa)*GAPWIND_FACTOR_Tile;

//	printf("Exposed Facing SnowWind Before Correction After GapWind: %f \n", SnowWind);
//	printf("Exposed Facing SnowRa Before Correction After GapWind: %f \n", SnowRa);

    OldTSurf = (*Tile)[NorthFacing].TSurf;
    (*Tile)[NorthFacing].SnowPackOutflow =
      SnowMelt(y, x, Dt, 2.+Z0_SNOW, 0.f, Z0_SNOW, SnowRa, LocalMet->AirDens,
        LocalMet->Eact, LocalMet->Lv, SnowNetShort, SnowLongIn,
        LocalMet->Press, (*Tile)[NorthFacing].RainFall, (*Tile)[NorthFacing].SnowFall,
        LocalMet->Tair, LocalMet->Vpd, SnowWind,
        &((*Tile)[NorthFacing].PackWater), &((*Tile)[NorthFacing].SurfWater),
        &((*Tile)[NorthFacing].Swq), &((*Tile)[NorthFacing].VaporMassFlux),
        &((*Tile)[NorthFacing].TPack), &((*Tile)[NorthFacing].TSurf), &((*Tile)[NorthFacing].MeltEnergy));

    /* Calculate the terms of the snow energy balance.  This is similar to the
    code in SnowPackEnergyBalance.c */
    Tmean = 0.5 * (OldTSurf + (*Tile)[NorthFacing].TSurf);

    /* Apply the stability correction to the aerodynamic resistance */
    if (SnowWind > 0.0)
      SnowRa /= StabilityCorrection(2.0f, 0.f, Tmean, LocalMet->Tair, SnowWind, Z0_SNOW);
    else
      SnowRa = DHSVM_HUGE;

	/* Calculate the saturated vapor pressure in the snow pack,
     (Equation 3.32, Bras 1990) */
	/* Write this in so they can be written out for model evaluation */
	(*Tile)[NorthFacing].EsSnow_tile = SatVaporPressure(Tmean);
    (*Tile)[NorthFacing].Eact_tile   = LocalMet->Eact;
    (*Tile)[NorthFacing].Ra_tile     = SnowRa;

//	printf("North Facing SnowRa After Correction: %f \n", SnowRa);

    /* convert snow surface temperature from C to K */
    Tmp = Tmean + 273.15;
    /* net shortwave radiation */
    (*Tile)[NorthFacing].Qsw = SnowNetShort;
    /* net longwave radiation */
	(*Tile)[NorthFacing].Qlin = SnowLongIn;
    (*Tile)[NorthFacing].Qlw = SnowLongIn - STEFAN * (Tmp * Tmp * Tmp * Tmp);
    /* sensible heat */
    (*Tile)[NorthFacing].Qs = LocalMet->AirDens * CP * (LocalMet->Tair - Tmean) / SnowRa;

    /* Calculate latent heat flux */
    if (Tmean >= 0.0) {
      /* Melt conditions: use latent heat of vaporization */
      (*Tile)[NorthFacing].Qe = LocalMet->Lv * (*Tile)[NorthFacing].VaporMassFlux * WATER_DENSITY;
    }
    else {
      /* Accumulation: use latent heat of sublimation (Eq. 3.19, Bras 1990 */
      Ls = (677. - 0.07 * Tmean) * JOULESPCAL * GRAMSPKG;
      (*Tile)[NorthFacing].Qe = Ls * (*Tile)[NorthFacing].VaporMassFlux * WATER_DENSITY;
    }
    (*Tile)[NorthFacing].Qe /= Dt;

    /* Calculate advected heat flux from rain */
    (*Tile)[NorthFacing].Qp = (CH_WATER * LocalMet->Tair * (*Tile)[NorthFacing].RainFall) / Dt;

    /* end of snow energy tems */

    /* Rainfall was added to SurfWater of the snow pack and has to be set to zero */
    (*Tile)[NorthFacing].RainFall = 0.0;
    (*Tile)[NorthFacing].MoistureFlux -= (*Tile)[NorthFacing].VaporMassFlux;

    /* Because we now have a new estimate of the snow surface temperature we
    can recalculate the longwave balance */
    Tsurf = Tile[NorthFacing]->TSurf;


    TileLongRadiation(VType, &(LocalVeg->Tile[NorthFacing]), Options,
                    CanopyRadAttOption, LocalMet->Lin, HeatFluxOption, LocalMet->Tair,
                    LocalVeg->Tcanopy, LocalSoil->TSurf, SType->Albedo);
  }
  else {
    (*Tile)[NorthFacing].SnowPackOutflow = 0.0;
    (*Tile)[NorthFacing].VaporMassFlux = 0.0;

	/* snow is all melt */
	(*Tile)[NorthFacing].Qs = 0.;
	(*Tile)[NorthFacing].Qe = 0.;
	(*Tile)[NorthFacing].Qp = 0.;
	(*Tile)[NorthFacing].Qsw = 0;
	(*Tile)[NorthFacing].Qlin = 0;
	(*Tile)[NorthFacing].Qlw = 0; /* not used in calculation */
	(*Tile)[NorthFacing].MeltEnergy = 0;
  }

  if ((*Tile)[NorthFacing].Swq > 0.0)
    (*Tile)[NorthFacing].HasSnow = TRUE;
  else
    (*Tile)[NorthFacing].HasSnow = FALSE;


  /********************* calculate SouthFacing snow melt *********************/

  /* SouthFacing */

  if ((*Tile)[SouthFacing].HasSnow || (*Tile)[SouthFacing].SnowFall > 0.0) {
    SnowLongIn = (*Tile)[SouthFacing].LongIn[1];
    SnowNetShort = (*Tile)[SouthFacing].NetShort[1];

    SnowWind = (*Tile)[SouthFacing].USnow * LocalMet->Wind;
    SnowRa = (*Tile)[SouthFacing].RaSnow / LocalMet->Wind;

//	printf("South Facing SnowWind Before Correction: %f \n", SnowWind);
//	printf("South Facing SnowRa Before Correction: %f \n", SnowRa);

    /* adjust the wind and Ra values for gap so they fall betwene open
    and forested values */

    tmp = VType->USnow*LocalMet->Wind; //forested snow wind
    SnowWind = tmp + (SnowWind- tmp)*GAPWIND_FACTOR_Tile;

    tmp = VType->RaSnow/LocalMet->Wind;
    SnowRa = tmp - (tmp-SnowRa)*GAPWIND_FACTOR_Tile;

//	printf("Exposed Facing SnowWind Before Correction After GapWind: %f \n", SnowWind);
//	printf("Exposed Facing SnowRa Before Correction After GapWind: %f \n", SnowRa);

    OldTSurf = (*Tile)[SouthFacing].TSurf;
    (*Tile)[SouthFacing].SnowPackOutflow =
      SnowMelt(y, x, Dt, 2.+Z0_SNOW, 0.f, Z0_SNOW, SnowRa, LocalMet->AirDens,
        LocalMet->Eact, LocalMet->Lv, SnowNetShort, SnowLongIn,
        LocalMet->Press, (*Tile)[SouthFacing].RainFall, (*Tile)[SouthFacing].SnowFall,
        LocalMet->Tair, LocalMet->Vpd, SnowWind,
        &((*Tile)[SouthFacing].PackWater), &((*Tile)[SouthFacing].SurfWater),
        &((*Tile)[SouthFacing].Swq), &((*Tile)[SouthFacing].VaporMassFlux),
        &((*Tile)[SouthFacing].TPack), &((*Tile)[SouthFacing].TSurf), &((*Tile)[SouthFacing].MeltEnergy));

    /* Calculate the terms of the snow energy balance.  This is similar to the
    code in SnowPackEnergyBalance.c */
    Tmean = 0.5 * (OldTSurf + (*Tile)[SouthFacing].TSurf);

    /* Apply the stability correction to the aerodynamic resistance */
    if (SnowWind > 0.0)
      SnowRa /= StabilityCorrection(2.0f, 0.f, Tmean, LocalMet->Tair, SnowWind, Z0_SNOW);
    else
      SnowRa = DHSVM_HUGE;

	/* Calculate the saturated vapor pressure in the snow pack,
     (Equation 3.32, Bras 1990) */
	/* Write this in so they can be written out for model evaluation */
	(*Tile)[SouthFacing].EsSnow_tile = SatVaporPressure(Tmean);
    (*Tile)[SouthFacing].Eact_tile   = LocalMet->Eact;
    (*Tile)[SouthFacing].Ra_tile     = SnowRa;
//	printf("South Facing SnowRa After Correction: %f \n", SnowRa);

    /* convert snow surface temperature from C to K */
    Tmp = Tmean + 273.15;
    /* net shortwave radiation */
    (*Tile)[SouthFacing].Qsw = SnowNetShort;
    /* net longwave radiation */
	(*Tile)[SouthFacing].Qlin = SnowLongIn;
    (*Tile)[SouthFacing].Qlw = SnowLongIn - STEFAN * (Tmp * Tmp * Tmp * Tmp);
    /* sensible heat */
    (*Tile)[SouthFacing].Qs = LocalMet->AirDens * CP * (LocalMet->Tair - Tmean) / SnowRa;

    /* Calculate latent heat flux */
    if (Tmean >= 0.0) {
      /* Melt conditions: use latent heat of vaporization */
      (*Tile)[SouthFacing].Qe = LocalMet->Lv * (*Tile)[SouthFacing].VaporMassFlux * WATER_DENSITY;
    }
    else {
      /* Accumulation: use latent heat of sublimation (Eq. 3.19, Bras 1990 */
      Ls = (677. - 0.07 * Tmean) * JOULESPCAL * GRAMSPKG;
      (*Tile)[SouthFacing].Qe = Ls * (*Tile)[SouthFacing].VaporMassFlux * WATER_DENSITY;
    }
    (*Tile)[SouthFacing].Qe /= Dt;

    /* Calculate advected heat flux from rain */
    (*Tile)[SouthFacing].Qp = (CH_WATER * LocalMet->Tair * (*Tile)[SouthFacing].RainFall) / Dt;

    /* end of snow energy tems */

    /* Rainfall was added to SurfWater of the snow pack and has to be set to zero */
    (*Tile)[SouthFacing].RainFall = 0.0;
    (*Tile)[SouthFacing].MoistureFlux -= (*Tile)[SouthFacing].VaporMassFlux;

    /* Because we now have a new estimate of the snow surface temperature we
    can recalculate the longwave balance */
    Tsurf = Tile[SouthFacing]->TSurf;


    TileLongRadiation(VType, &(LocalVeg->Tile[SouthFacing]), Options,
                    CanopyRadAttOption, LocalMet->Lin, HeatFluxOption, LocalMet->Tair,
                    LocalVeg->Tcanopy, LocalSoil->TSurf, SType->Albedo);
  }
  else {
    (*Tile)[SouthFacing].SnowPackOutflow = 0.0;
    (*Tile)[SouthFacing].VaporMassFlux = 0.0;

	/* snow is all melt */
	(*Tile)[SouthFacing].Qs = 0.;
	(*Tile)[SouthFacing].Qe = 0.;
	(*Tile)[SouthFacing].Qp = 0.;
	(*Tile)[SouthFacing].Qsw = 0;
	(*Tile)[SouthFacing].Qlin = 0;
	(*Tile)[SouthFacing].Qlw = 0; /* not used in calculation */
	(*Tile)[SouthFacing].MeltEnergy = 0;
  }

  if ((*Tile)[SouthFacing].Swq > 0.0)
    (*Tile)[SouthFacing].HasSnow = TRUE;
  else
    (*Tile)[SouthFacing].HasSnow = FALSE;

  /********************* calculate Exposed snow melt *********************/


  /* Exposed */

  if ((*Tile)[Exposed].HasSnow || (*Tile)[Exposed].SnowFall > 0.0) {
    SnowLongIn = (*Tile)[Exposed].LongIn[1];
    SnowNetShort = (*Tile)[Exposed].NetShort[1];

    SnowWind = (*Tile)[Exposed].USnow * LocalMet->Wind;
    SnowRa = (*Tile)[Exposed].RaSnow / LocalMet->Wind;

//	printf("Exposed Facing SnowWind Before Correction: %f \n", SnowWind);
//	printf("Exposed Facing SnowRa Before Correction: %f \n", SnowRa);

    /* adjust the wind and Ra values for gap so they fall betwene open
    and forested values */

    tmp = VType->USnow*LocalMet->Wind; //forested snow wind
    SnowWind = tmp + (SnowWind- tmp)*GAPWIND_FACTOR_Tile;

    tmp = VType->RaSnow/LocalMet->Wind;
    SnowRa = tmp - (tmp-SnowRa)*GAPWIND_FACTOR_Tile;

//	printf("Exposed Facing SnowWind Before Correction After GapWind: %f \n", SnowWind);
//	printf("Exposed Facing SnowRa Before Correction After GapWind: %f \n", SnowRa);

    OldTSurf = (*Tile)[Exposed].TSurf;
    (*Tile)[Exposed].SnowPackOutflow =
      SnowMelt(y, x, Dt, 2.+Z0_SNOW, 0.f, Z0_SNOW, SnowRa, LocalMet->AirDens,
        LocalMet->Eact, LocalMet->Lv, SnowNetShort, SnowLongIn,
        LocalMet->Press, (*Tile)[Exposed].RainFall, (*Tile)[Exposed].SnowFall,
        LocalMet->Tair, LocalMet->Vpd, SnowWind,
        &((*Tile)[Exposed].PackWater), &((*Tile)[Exposed].SurfWater),
        &((*Tile)[Exposed].Swq), &((*Tile)[Exposed].VaporMassFlux),
        &((*Tile)[Exposed].TPack), &((*Tile)[Exposed].TSurf), &((*Tile)[Exposed].MeltEnergy));

    /* Calculate the terms of the snow energy balance.  This is similar to the
    code in SnowPackEnergyBalance.c */
    Tmean = 0.5 * (OldTSurf + (*Tile)[Exposed].TSurf);

    /* Apply the stability correction to the aerodynamic resistance */
    if (SnowWind > 0.0)
      SnowRa /= StabilityCorrection(2.0f, 0.f, Tmean, LocalMet->Tair, SnowWind, Z0_SNOW);
    else
      SnowRa = DHSVM_HUGE;

	/* Calculate the saturated vapor pressure in the snow pack,
     (Equation 3.32, Bras 1990) */
	/* Write this in so they can be written out for model evaluation */
	(*Tile)[Exposed].EsSnow_tile = SatVaporPressure(Tmean);
    (*Tile)[Exposed].Eact_tile   = LocalMet->Eact;
    (*Tile)[Exposed].Ra_tile     = SnowRa;
//	printf("Exposed Facing SnowRa After Correction: %f \n", SnowRa);

    /* convert snow surface temperature from C to K */
    Tmp = Tmean + 273.15;
    /* net shortwave radiation */
    (*Tile)[Exposed].Qsw = SnowNetShort;
    /* net longwave radiation */
	(*Tile)[Exposed].Qlin = SnowLongIn;
    (*Tile)[Exposed].Qlw = SnowLongIn - STEFAN * (Tmp * Tmp * Tmp * Tmp);
    /* sensible heat */
    (*Tile)[Exposed].Qs = LocalMet->AirDens * CP * (LocalMet->Tair - Tmean) / SnowRa;

    /* Calculate latent heat flux */
    if (Tmean >= 0.0) {
      /* Melt conditions: use latent heat of vaporization */
      (*Tile)[Exposed].Qe = LocalMet->Lv * (*Tile)[Exposed].VaporMassFlux * WATER_DENSITY;
    }
    else {
      /* Accumulation: use latent heat of sublimation (Eq. 3.19, Bras 1990 */
      Ls = (677. - 0.07 * Tmean) * JOULESPCAL * GRAMSPKG;
      (*Tile)[Exposed].Qe = Ls * (*Tile)[Exposed].VaporMassFlux * WATER_DENSITY;
    }
    (*Tile)[Exposed].Qe /= Dt;

    /* Calculate advected heat flux from rain */
    (*Tile)[Exposed].Qp = (CH_WATER * LocalMet->Tair * (*Tile)[Exposed].RainFall) / Dt;

    /* end of snow energy tems */

    /* Rainfall was added to SurfWater of the snow pack and has to be set to zero */
    (*Tile)[Exposed].RainFall = 0.0;
    (*Tile)[Exposed].MoistureFlux -= (*Tile)[Exposed].VaporMassFlux;

    /* Because we now have a new estimate of the snow surface temperature we
    can recalculate the longwave balance */
    Tsurf = Tile[Exposed]->TSurf;


    TileLongRadiation(VType, &(LocalVeg->Tile[Exposed]), Options,
                    CanopyRadAttOption, LocalMet->Lin, HeatFluxOption, LocalMet->Tair,
                    LocalVeg->Tcanopy, LocalSoil->TSurf, SType->Albedo);
  }
  else {
    (*Tile)[Exposed].SnowPackOutflow = 0.0;
    (*Tile)[Exposed].VaporMassFlux = 0.0;

	/* snow is all melt */
	(*Tile)[Exposed].Qs = 0.;
	(*Tile)[Exposed].Qe = 0.;
	(*Tile)[Exposed].Qp = 0.;
	(*Tile)[Exposed].Qsw = 0;
	(*Tile)[Exposed].Qlin = 0;
	(*Tile)[Exposed].Qlw = 0; /* not used in calculation */
	(*Tile)[Exposed].MeltEnergy = 0;
  }

  if ((*Tile)[Exposed].Swq > 0.0)
    (*Tile)[Exposed].HasSnow = TRUE;
  else
    (*Tile)[Exposed].HasSnow = FALSE;
}

/*****************************************************************************
Function name: NoOverStoryET()

Purpose      : Calculate the ET
*****************************************************************************/
void NoOverStoryET(TileStruct **Tile, int NSoil, VEGTABLE *VType,
  VEGPIX *LocalVeg, SOILTABLE *SType, SOILPIX *LocalSoil, PIXMET *LocalMet,
  EVAPPIX *LocalEvap, ROADSTRUCT *LocalNetwork, int Dt, float UpperRa,
  float LowerRa)
{
  float NetRadiation;		/* Total Net long- and shortwave radiation (W/m2) */
  float Rp;					/* radiation flux in visible part of the spectrum (W/m^2) */

  /********************** for NorthFacing **********************/

  if ((*Tile)[NorthFacing].HasSnow != TRUE && VType->UnderStory == TRUE) {

    Rp = VISFRACT * (*Tile)[NorthFacing].NetShort[1];
    NetRadiation =
      (*Tile)[NorthFacing].NetShort[1] +
      (*Tile)[NorthFacing].LongIn[1] - (*Tile)[NorthFacing].LongOut[1];

    EvapoTranspiration(1, 1, Dt, 1, LocalMet, NetRadiation,
      Rp, VType, SType, (*Tile)[NorthFacing].MoistureFlux, (*Tile)[NorthFacing].Moist,
      LocalSoil->Temp, &((*Tile)[NorthFacing].IntRain[1]),
      (*Tile)[NorthFacing].EPot, (*Tile)[NorthFacing].EInt, (*Tile)[NorthFacing].ESoil,
      (*Tile)[NorthFacing].EAct, &((*Tile)[NorthFacing].ETot), LocalNetwork->Adjust, LowerRa,
	  (*Tile)[NorthFacing].NorthFacingInt,LocalVeg->NFfrac);

    (*Tile)[NorthFacing].MoistureFlux += (*Tile)[NorthFacing].EAct[1] + (*Tile)[NorthFacing].EInt[1];

    (*Tile)[NorthFacing].NetRadiation[0] = 0.;
    (*Tile)[NorthFacing].NetRadiation[1] = NetRadiation;
  }
  else if (VType->UnderStory == TRUE) {
    (*Tile)[NorthFacing].EAct[1] = 0.;
    (*Tile)[NorthFacing].EInt[1] = 0.;
    (*Tile)[NorthFacing].NetRadiation[1] = 0.;
    (*Tile)[NorthFacing].NetRadiation[0] = 0.;
  }

  /* Calculate soil evaporation from the upper soil layer if no snow is
  present and there is no understory */
  if ((*Tile)[NorthFacing].HasSnow != TRUE && VType->UnderStory != TRUE) {
    NetRadiation =
      (*Tile)[NorthFacing].NetShort[1] + (*Tile)[NorthFacing].LongIn[1] - (*Tile)[NorthFacing].LongOut[1];
    (*Tile)[NorthFacing].NetRadiation[1] = NetRadiation;
    (*Tile)[NorthFacing].NetRadiation[0] = 0.;
    (*Tile)[NorthFacing].EvapSoil =
    SoilEvaporation(Dt, LocalMet->Tair, LocalMet->Slope, LocalMet->Gamma,
        LocalMet->Lv, LocalMet->AirDens, LocalMet->Vpd,
        NetRadiation, UpperRa, (*Tile)[NorthFacing].MoistureFlux, SType->Porosity[1],
        SType->FCap[1], SType->Ks[1], SType->Press[1], SType->PoreDist[1],
        VType->RootDepth[1], &((*Tile)[NorthFacing].Moist[1]),
        LocalNetwork->Adjust[1]);
  }
  else
    (*Tile)[NorthFacing].EvapSoil = 0.0;

  (*Tile)[NorthFacing].MoistureFlux += (*Tile)[NorthFacing].EvapSoil;
  (*Tile)[NorthFacing].ETot += (*Tile)[NorthFacing].EvapSoil;


  /********************** for SouthFacing **********************/

  if ((*Tile)[SouthFacing].HasSnow != TRUE && VType->UnderStory == TRUE) {
    Rp = VISFRACT * (*Tile)[SouthFacing].NetShort[1];
    NetRadiation =
      (*Tile)[SouthFacing].NetShort[1] +
      (*Tile)[SouthFacing].LongIn[1] - (*Tile)[SouthFacing].LongOut[1];

    EvapoTranspiration(1, 1, Dt, 1, LocalMet, NetRadiation,
      Rp, VType, SType, (*Tile)[SouthFacing].MoistureFlux, (*Tile)[SouthFacing].Moist,
      LocalSoil->Temp, &((*Tile)[SouthFacing].IntRain[1]),
      (*Tile)[SouthFacing].EPot, (*Tile)[SouthFacing].EInt, (*Tile)[SouthFacing].ESoil,
      (*Tile)[SouthFacing].EAct, &((*Tile)[SouthFacing].ETot), LocalNetwork->Adjust, LowerRa,
      (*Tile)[SouthFacing].SouthFacingInt,LocalVeg->SFfrac);

    (*Tile)[SouthFacing].MoistureFlux += (*Tile)[SouthFacing].EAct[1] + (*Tile)[SouthFacing].EInt[1];

    (*Tile)[SouthFacing].NetRadiation[1] = NetRadiation;
    (*Tile)[SouthFacing].NetRadiation[0] = 0.;
  }
  else if (VType->UnderStory == TRUE) {
    (*Tile)[SouthFacing].EAct[1] = 0.;
    (*Tile)[SouthFacing].EInt[1] = 0.;
    (*Tile)[SouthFacing].NetRadiation[1] = 0.;
    (*Tile)[SouthFacing].NetRadiation[0] = 0.;
  }

  /* Calculate soil evaporation from the upper soil layer if no snow is
  present and there is no understory */
  if ((*Tile)[SouthFacing].HasSnow != TRUE && VType->UnderStory != TRUE) {
    NetRadiation =
      (*Tile)[SouthFacing].NetShort[1] + (*Tile)[SouthFacing].LongIn[1] - (*Tile)[SouthFacing].LongOut[1];
    (*Tile)[SouthFacing].NetRadiation[1] = NetRadiation;
    (*Tile)[SouthFacing].NetRadiation[0] = 0.;
    (*Tile)[SouthFacing].EvapSoil =
    SoilEvaporation(Dt, LocalMet->Tair, LocalMet->Slope, LocalMet->Gamma,
        LocalMet->Lv, LocalMet->AirDens, LocalMet->Vpd,
        NetRadiation, UpperRa, (*Tile)[SouthFacing].MoistureFlux, SType->Porosity[1],
        SType->FCap[1], SType->Ks[1], SType->Press[1], SType->PoreDist[1],
        VType->RootDepth[1], &((*Tile)[SouthFacing].Moist[1]),
        LocalNetwork->Adjust[1]);
  }
  else
    (*Tile)[SouthFacing].EvapSoil = 0.0;

  (*Tile)[SouthFacing].MoistureFlux += (*Tile)[SouthFacing].EvapSoil;
  (*Tile)[SouthFacing].ETot += (*Tile)[SouthFacing].EvapSoil;



    /********************** for Exposed **********************/

  if ((*Tile)[Exposed].HasSnow != TRUE && VType->UnderStory == TRUE) {
    Rp = VISFRACT * (*Tile)[Exposed].NetShort[1];
    NetRadiation =
      (*Tile)[Exposed].NetShort[1] +
      (*Tile)[Exposed].LongIn[1] - (*Tile)[Exposed].LongOut[1];

    EvapoTranspiration(1, 1, Dt, 1, LocalMet, NetRadiation,
      Rp, VType, SType, (*Tile)[Exposed].MoistureFlux, (*Tile)[Exposed].Moist,
      LocalSoil->Temp, &((*Tile)[Exposed].IntRain[1]),
      (*Tile)[Exposed].EPot, (*Tile)[Exposed].EInt, (*Tile)[Exposed].ESoil,
      (*Tile)[Exposed].EAct, &((*Tile)[Exposed].ETot), LocalNetwork->Adjust, LowerRa,
      (*Tile)[Exposed].ExposedInt,LocalVeg->EXPfrac);

    (*Tile)[Exposed].MoistureFlux += (*Tile)[Exposed].EAct[1] + (*Tile)[Exposed].EInt[1];

    (*Tile)[Exposed].NetRadiation[1] = NetRadiation;
    (*Tile)[Exposed].NetRadiation[0] = 0.;
  }
  else if (VType->UnderStory == TRUE) {
    (*Tile)[Exposed].EAct[1] = 0.;
    (*Tile)[Exposed].EInt[1] = 0.;
    (*Tile)[Exposed].NetRadiation[1] = 0.;
    (*Tile)[Exposed].NetRadiation[0] = 0.;
  }

  /* Calculate soil evaporation from the upper soil layer if no snow is
  present and there is no understory */
  if ((*Tile)[Exposed].HasSnow != TRUE && VType->UnderStory != TRUE) {
    NetRadiation =
      (*Tile)[Exposed].NetShort[1] + (*Tile)[Exposed].LongIn[1] - (*Tile)[Exposed].LongOut[1];
    (*Tile)[Exposed].NetRadiation[1] = NetRadiation;
    (*Tile)[Exposed].NetRadiation[0] = 0.;
    (*Tile)[Exposed].EvapSoil =
    SoilEvaporation(Dt, LocalMet->Tair, LocalMet->Slope, LocalMet->Gamma,
        LocalMet->Lv, LocalMet->AirDens, LocalMet->Vpd,
        NetRadiation, UpperRa, (*Tile)[Exposed].MoistureFlux, SType->Porosity[1],
        SType->FCap[1], SType->Ks[1], SType->Press[1], SType->PoreDist[1],
        VType->RootDepth[1], &((*Tile)[Exposed].Moist[1]),
        LocalNetwork->Adjust[1]);
  }
  else
    (*Tile)[Exposed].EvapSoil = 0.0;

  (*Tile)[Exposed].MoistureFlux += (*Tile)[Exposed].EvapSoil;
  (*Tile)[Exposed].ETot += (*Tile)[Exposed].EvapSoil;

}

/*****************************************************************************
Function name: OverStoryInterceptSnowMelt()

Purpose      : .
*****************************************************************************/

void OverStoryInterceptSnowMelt(OPTIONSTRUCT *Options, int HeatFluxOption,
  int y, int x, int Dt, int NVegLActual, TileStruct **Tile, VEGTABLE *VType,
  PIXRAD *LocalRad, PIXMET *LocalMet, float UpperRa, float UpperWind, VEGPIX *LocalVeg,
  SOILPIX *LocalSoil, int CanopyRadAttOption, SOILTABLE *SType)
{
  float Tsurf;
  float SnowLongIn;			/* Incoming longwave radiation at snow surface (W/m2) */
  float SnowNetShort;		/* Net amount of short wave radiation at the snow surface (W/m2) */
  float SnowRa;				/* Aerodynamic resistance for snow */
  float SnowWind;		    /* Wind 2 m above snow */
  float OldTSurf;
  float Tmean;


  if (((*Tile)[ForestTile].IntSnow[0] || (*Tile)[ForestTile].SnowFall > 0.0)) {
    SnowInterception(Options, y, x, Dt, VType->Fract[0], VType->Vf,
      VType->LAI[0], VType->MaxInt[0], VType->MaxSnowInt, VType->MDRatio,
      VType->SnowIntEff, UpperRa, LocalMet->AirDens,
      LocalMet->Eact, LocalMet->Lv, LocalRad, LocalMet->Press,
      LocalMet->Tair, LocalMet->Vpd, UpperWind,
      &((*Tile)[ForestTile].RainFall), &((*Tile)[ForestTile].SnowFall),
      &((*Tile)[ForestTile].IntRain[0]), &((*Tile)[ForestTile].IntSnow[0]),
      &((*Tile)[ForestTile].TempIntStorage), &((*Tile)[ForestTile].CanopyVaporMassFlux),
      &((*Tile)[ForestTile].Tcanopy), &((*Tile)[ForestTile].MeltEnergy), VType->Height,
      VType->UnderStory);
    (*Tile)[ForestTile].MoistureFlux -= (*Tile)[ForestTile].CanopyVaporMassFlux;

    /* Because we now have a new estimate of the canopy temperature we can
    recalculate the longwave balance */
    /* update longwave radiation */
    TileLongRadiation(VType, &(LocalVeg->Tile[ForestTile]), Options,
                    CanopyRadAttOption, LocalMet->Lin, HeatFluxOption, LocalMet->Tair,
                    LocalVeg->Tcanopy, LocalSoil->TSurf, SType->Albedo);
  }
  /* if no snow */
  else if (VType->NVegLayers > 0) {
    (*Tile)[ForestTile].Tcanopy = LocalMet->Tair;
    (*Tile)[ForestTile].CanopyVaporMassFlux = 0.0;
    (*Tile)[ForestTile].TempIntStorage = 0.0;
    InterceptionStorage(NVegLActual, VType->MaxInt, VType->Fract, (*Tile)[ForestTile].IntRain,
      &((*Tile)[ForestTile].RainFall));
  }

  /* if snow is present, simulate the snow pack dynamics */
  if ((*Tile)[ForestTile].HasSnow || (*Tile)[ForestTile].SnowFall > 0.0) {

    /* SnowLongIn = LocalRad->LongIn[1];
    SnowNetShort = LocalRad->NetShort[1]; */
    SnowLongIn = (*Tile)[ForestTile].LongIn[1];
    SnowNetShort = (*Tile)[ForestTile].NetShort[1];
    SnowWind = VType->USnow * LocalMet->Wind;
    SnowRa = VType->RaSnow / LocalMet->Wind;

    OldTSurf = (*Tile)[ForestTile].TSurf;
    (*Tile)[ForestTile].SnowPackOutflow =
      SnowMelt(y, x, Dt, 2.+Z0_SNOW, 0.f, Z0_SNOW, SnowRa, LocalMet->AirDens,
        LocalMet->Eact, LocalMet->Lv, SnowNetShort, SnowLongIn,
        LocalMet->Press, (*Tile)[ForestTile].RainFall, (*Tile)[ForestTile].SnowFall,
        LocalMet->Tair, LocalMet->Vpd, SnowWind,
        &((*Tile)[ForestTile].PackWater), &((*Tile)[ForestTile].SurfWater),
        &((*Tile)[ForestTile].Swq), &((*Tile)[ForestTile].VaporMassFlux),
        &((*Tile)[ForestTile].TPack), &((*Tile)[ForestTile].TSurf), &((*Tile)[ForestTile].MeltEnergy));

    /* Rainfall was added to SurfWater of the snow pack and has to be set to zero */
    (*Tile)[ForestTile].RainFall = 0.0;
    (*Tile)[ForestTile].MoistureFlux -= (*Tile)[ForestTile].VaporMassFlux;

    /* Because we now have a new estimate of the snow surface temperature we
    can recalculate the longwave balance */
    /* Because we now have a new estimate of the canopy temperature we can
    recalculate the longwave balance */
    /* update longwave radiation */
    TileLongRadiation(VType, &(LocalVeg->Tile[ForestTile]), Options,
                    CanopyRadAttOption, LocalMet->Lin, HeatFluxOption, LocalMet->Tair,
                    LocalVeg->Tcanopy, LocalSoil->TSurf, SType->Albedo);

    /* Calculate the terms of the snow energy balance.  This is similar to the
    code in SnowPackEnergyBalance.c */
    Tmean = 0.5 * (OldTSurf + (*Tile)[ForestTile].TSurf);

    /* Apply the stability correction to the aerodynamic resistance */
    if (SnowWind > 0.0)
      SnowRa /= StabilityCorrection(2.0f, 0.f, Tmean, LocalMet->Tair, SnowWind, Z0_SNOW);
    else
      SnowRa = DHSVM_HUGE;

	/* Write this in so they can be written out for model evaluation */
	(*Tile)[ForestTile].EsSnow_tile = SatVaporPressure(Tmean);
    (*Tile)[ForestTile].Eact_tile   = LocalMet->Eact;
    (*Tile)[ForestTile].Ra_tile     = SnowRa;

  }
  else {
    (*Tile)[ForestTile].SnowPackOutflow = 0.0;
    (*Tile)[ForestTile].VaporMassFlux = 0.0;
  }

  /* Determine whether a snow pack is still present, or whether everything
  has melted */
  if ((*Tile)[ForestTile].Swq > 0.0)
    (*Tile)[ForestTile].HasSnow = TRUE;
  else
    (*Tile)[ForestTile].HasSnow = FALSE;
}


/*****************************************************************************
Function name: OverStoryET()

Purpose      :
*****************************************************************************/
void OverStoryET(int Dt, TileStruct **Tile,
  SOILTABLE *SType, VEGTABLE *VType, PIXRAD *LocalRad, PIXMET *LocalMet,
  SOILPIX *LocalSoil, ROADSTRUCT *LocalNetwork, float UpperRa, float LowerRa, VEGPIX *LocalVeg)
{
  float Rp;
  float NetRadiation;

  if (VType->OverStory == TRUE) {
    Rp = VISFRACT * (*Tile)[ForestTile].NetShort[0];
    NetRadiation = (*Tile)[ForestTile].NetShort[0] +
      (*Tile)[ForestTile].LongIn[0] - 2 * VType->Vf *(*Tile)[ForestTile].LongOut[0];
    (*Tile)[ForestTile].NetRadiation[0] = NetRadiation;

    EvapoTranspiration(0, 1, Dt, VType->Fract[0], LocalMet, NetRadiation,
      Rp, VType, SType, (*Tile)[ForestTile].MoistureFlux, (*Tile)[ForestTile].Moist, LocalSoil->Temp,
      &((*Tile)[ForestTile].IntRain[0]), (*Tile)[ForestTile].EPot, (*Tile)[ForestTile].EInt, (*Tile)[ForestTile].ESoil,
      (*Tile)[ForestTile].EAct, &((*Tile)[ForestTile].ETot), LocalNetwork->Adjust, UpperRa,
      (*Tile)[ForestTile].ForestInt,LocalVeg->FORfrac);
    (*Tile)[ForestTile].MoistureFlux += (*Tile)[ForestTile].EAct[0] + (*Tile)[ForestTile].EInt[0];

    if ((*Tile)[ForestTile].HasSnow != TRUE && VType->UnderStory == TRUE) {
      Rp = VISFRACT * LocalRad->NetShort[1];
      NetRadiation =
        LocalRad->NetShort[1] +
        LocalRad->LongIn[1] - VType->Fract[1] * LocalRad->LongOut[1];
      LocalRad->NetRadiation[1] = NetRadiation;
      EvapoTranspiration(1, 1, Dt, VType->Fract[1], LocalMet, NetRadiation,
        Rp, VType, SType, (*Tile)[ForestTile].MoistureFlux, (*Tile)[ForestTile].Moist, LocalSoil->Temp,
        &((*Tile)[ForestTile].IntRain[1]), (*Tile)[ForestTile].EPot, (*Tile)[ForestTile].EInt, (*Tile)[ForestTile].ESoil,
        (*Tile)[ForestTile].EAct, &((*Tile)[ForestTile].ETot), LocalNetwork->Adjust, LowerRa,
        (*Tile)[ForestTile].ForestInt,LocalVeg->FORfrac);
      (*Tile)[ForestTile].MoistureFlux += (*Tile)[ForestTile].EAct[1] + (*Tile)[ForestTile].EInt[1];
    }
    else if (VType->UnderStory == TRUE) {
      (*Tile)[ForestTile].EAct[1] = 0.;
      (*Tile)[ForestTile].EInt[1] = 0.;
      (*Tile)[ForestTile].NetRadiation[1] = 0.;
    }
  }

  /* Calculate soil evaporation from the upper soil layer if no snow is
  present and there is no understory */
  if ((*Tile)[ForestTile].HasSnow != TRUE && VType->UnderStory != TRUE) {
    if (VType->OverStory == TRUE) {
      NetRadiation =
        (*Tile)[ForestTile].NetShort[1] + (*Tile)[ForestTile].LongIn[1] - (*Tile)[ForestTile].LongOut[1];
      (*Tile)[ForestTile].NetRadiation[1] = NetRadiation;
    }
    (*Tile)[ForestTile].EvapSoil =
    SoilEvaporation(Dt, LocalMet->Tair, LocalMet->Slope, LocalMet->Gamma,
        LocalMet->Lv, LocalMet->AirDens, LocalMet->Vpd,
        NetRadiation, LowerRa, (*Tile)[ForestTile].MoistureFlux, SType->Porosity[0],
        SType->FCap[0], SType->Ks[0], SType->Press[0], SType->PoreDist[0],
        VType->RootDepth[0], &((*Tile)[ForestTile].Moist[0]),
        LocalNetwork->Adjust[0]);
  }
  else
    (*Tile)[ForestTile].EvapSoil = 0.0;

  (*Tile)[ForestTile].MoistureFlux += (*Tile)[ForestTile].EvapSoil;
  (*Tile)[ForestTile].ETot += (*Tile)[ForestTile].EvapSoil;
}


/*****************************************************************************
Function name: AggregateTile()

Purpose      : Aggregate the gap and non-gap mass balance variables based
on area weight.
*****************************************************************************/
void AggregateTile(TileStruct **Tile, VEGPIX *LocalVeg,
  SOILPIX *LocalSoil, SNOWPIX *LocalSnow, EVAPPIX *LocalEvap,
  PRECIPPIX *LocalPrecip, PIXRAD *LocalRad, int NSoil, int NVeg)
{
  int i, j;
    /* printf("Exposed fraction 1 = %f \n", (1 - (LocalVeg->FORfrac + LocalVeg->NFfrac + LocalVeg->SFfrac)));
    printf("FORfrac fraction 2 = %f \n", LocalVeg->FORfrac);
    printf("NFfrac fraction 2 = %f \n", LocalVeg->NFfrac);
    printf("SFfrac fraction 2 = %f \n", LocalVeg->SFfrac);
    printf("EXPfrac fraction 2 = %f \n", LocalVeg->EXPfrac);
    printf("Sum fraction 2 = %f \n", (LocalVeg->FORfrac + LocalVeg->NFfrac + LocalVeg->SFfrac + LocalVeg->EXPfrac)); */

  LocalPrecip->RainFall =
    (1 - (LocalVeg->FORfrac + LocalVeg->NFfrac + LocalVeg->SFfrac))*(*Tile)[Exposed].RainFall + LocalVeg->FORfrac*(*Tile)[ForestTile].RainFall + LocalVeg->NFfrac*(*Tile)[NorthFacing].RainFall + LocalVeg->SFfrac*(*Tile)[SouthFacing].RainFall;
  LocalPrecip->SnowFall =
    (1 - (LocalVeg->FORfrac + LocalVeg->NFfrac + LocalVeg->SFfrac))*(*Tile)[Exposed].SnowFall + LocalVeg->FORfrac*(*Tile)[ForestTile].SnowFall + LocalVeg->NFfrac*(*Tile)[NorthFacing].SnowFall + LocalVeg->SFfrac*(*Tile)[SouthFacing].SnowFall;
  LocalPrecip->Precip =
    (1 - (LocalVeg->FORfrac + LocalVeg->NFfrac + LocalVeg->SFfrac))*(*Tile)[Exposed].Precip + LocalVeg->FORfrac*(*Tile)[ForestTile].Precip + LocalVeg->NFfrac*(*Tile)[NorthFacing].Precip + LocalVeg->SFfrac*(*Tile)[SouthFacing].Precip;

  LocalSnow->Outflow =
    (1 - (LocalVeg->FORfrac + LocalVeg->NFfrac + LocalVeg->SFfrac))*(*Tile)[Exposed].SnowPackOutflow + LocalVeg->FORfrac*(*Tile)[ForestTile].SnowPackOutflow + LocalVeg->NFfrac*(*Tile)[NorthFacing].SnowPackOutflow + LocalVeg->SFfrac*(*Tile)[SouthFacing].SnowPackOutflow;
  LocalSnow->CanopyVaporMassFlux =
    (1 - (LocalVeg->FORfrac + LocalVeg->NFfrac + LocalVeg->SFfrac))*(*Tile)[Exposed].CanopyVaporMassFlux + LocalVeg->FORfrac*(*Tile)[ForestTile].CanopyVaporMassFlux + LocalVeg->NFfrac*(*Tile)[NorthFacing].CanopyVaporMassFlux + LocalVeg->SFfrac*(*Tile)[SouthFacing].CanopyVaporMassFlux;
  LocalSnow->VaporMassFlux =
    (1 - (LocalVeg->FORfrac + LocalVeg->NFfrac + LocalVeg->SFfrac))*(*Tile)[Exposed].VaporMassFlux + LocalVeg->FORfrac*(*Tile)[ForestTile].VaporMassFlux + LocalVeg->NFfrac*(*Tile)[NorthFacing].VaporMassFlux + LocalVeg->SFfrac*(*Tile)[SouthFacing].VaporMassFlux;

  for (i = 0; i < 2; i++) {
    LocalRad->NetShort[i] =
      (1 - (LocalVeg->FORfrac + LocalVeg->NFfrac + LocalVeg->SFfrac))*(*Tile)[Exposed].NetShort[i] + LocalVeg->FORfrac*(*Tile)[ForestTile].NetShort[i] +  LocalVeg->NFfrac*(*Tile)[NorthFacing].NetShort[i] + LocalVeg->SFfrac*(*Tile)[SouthFacing].NetShort[i];
    LocalRad->LongIn[i] =
      (1 - (LocalVeg->FORfrac + LocalVeg->NFfrac + LocalVeg->SFfrac))*(*Tile)[Exposed].LongIn[i] + LocalVeg->FORfrac*(*Tile)[ForestTile].LongIn[i] + LocalVeg->NFfrac*(*Tile)[NorthFacing].LongIn[i] + LocalVeg->SFfrac*(*Tile)[SouthFacing].LongIn[i];
    LocalRad->LongOut[i] =
      (1 - (LocalVeg->FORfrac + LocalVeg->NFfrac + LocalVeg->SFfrac))*(*Tile)[Exposed].LongOut[i] + LocalVeg->FORfrac*(*Tile)[ForestTile].LongOut[i] + LocalVeg->NFfrac*(*Tile)[NorthFacing].LongOut[i] + LocalVeg->SFfrac*(*Tile)[SouthFacing].LongOut[i];
  }

  LocalSnow->Swq =
    (1 - (LocalVeg->FORfrac + LocalVeg->NFfrac + LocalVeg->SFfrac))*(*Tile)[Exposed].Swq + LocalVeg->FORfrac*(*Tile)[ForestTile].Swq + LocalVeg->NFfrac*(*Tile)[NorthFacing].Swq + LocalVeg->SFfrac*(*Tile)[SouthFacing].Swq;
  LocalSnow->TPack =
    (1 - (LocalVeg->FORfrac + LocalVeg->NFfrac + LocalVeg->SFfrac))*(*Tile)[Exposed].TPack + LocalVeg->FORfrac*(*Tile)[ForestTile].TPack + LocalVeg->NFfrac*(*Tile)[NorthFacing].TPack + LocalVeg->SFfrac*(*Tile)[SouthFacing].TPack;
  LocalSnow->PackWater =
    (1 - (LocalVeg->FORfrac + LocalVeg->NFfrac + LocalVeg->SFfrac))*(*Tile)[Exposed].PackWater + LocalVeg->FORfrac*(*Tile)[ForestTile].PackWater + LocalVeg->NFfrac*(*Tile)[NorthFacing].PackWater + LocalVeg->SFfrac*(*Tile)[SouthFacing].PackWater;
  LocalSnow->SurfWater =
    (1 - (LocalVeg->FORfrac + LocalVeg->NFfrac + LocalVeg->SFfrac))*(*Tile)[Exposed].SurfWater + LocalVeg->FORfrac*(*Tile)[ForestTile].SurfWater + LocalVeg->NFfrac*(*Tile)[NorthFacing].SurfWater + LocalVeg->SFfrac*(*Tile)[SouthFacing].SurfWater;

  for (j = 0; j <= NSoil; j++) {
    LocalSoil->Moist[j] = (1 - (LocalVeg->FORfrac + LocalVeg->NFfrac + LocalVeg->SFfrac))*(*Tile)[Exposed].Moist[j] +
      LocalVeg->FORfrac*(*Tile)[ForestTile].Moist[j] + LocalVeg->NFfrac*(*Tile)[NorthFacing].Moist[j] + LocalVeg->SFfrac*(*Tile)[SouthFacing].Moist[j];
  }

  LocalVeg->MoistureFlux =
    (1 - (LocalVeg->FORfrac + LocalVeg->NFfrac + LocalVeg->SFfrac))*(*Tile)[Exposed].MoistureFlux + LocalVeg->FORfrac*(*Tile)[ForestTile].MoistureFlux + LocalVeg->NFfrac*(*Tile)[NorthFacing].MoistureFlux + LocalVeg->SFfrac*(*Tile)[SouthFacing].MoistureFlux;

  LocalVeg->MeltEnergy =
    (1 - (LocalVeg->FORfrac + LocalVeg->NFfrac + LocalVeg->SFfrac))*(*Tile)[Exposed].MeltEnergy + LocalVeg->FORfrac*(*Tile)[ForestTile].MeltEnergy + LocalVeg->NFfrac*(*Tile)[NorthFacing].MeltEnergy + LocalVeg->SFfrac*(*Tile)[SouthFacing].MeltEnergy;

  /* Intercepted rain/snow */
  for (i = 0; i < NVeg; i++) {
    LocalPrecip->IntRain[i] =
      (1 - (LocalVeg->FORfrac + LocalVeg->NFfrac + LocalVeg->SFfrac))*(*Tile)[Exposed].IntRain[i] + LocalVeg->FORfrac*(*Tile)[ForestTile].IntRain[i] + LocalVeg->NFfrac*(*Tile)[NorthFacing].IntRain[i] + LocalVeg->SFfrac*(*Tile)[SouthFacing].IntRain[i];
    LocalPrecip->IntSnow[i] =
      (1 - (LocalVeg->FORfrac + LocalVeg->NFfrac + LocalVeg->SFfrac))*(*Tile)[Exposed].IntSnow[i] + LocalVeg->FORfrac*(*Tile)[ForestTile].IntSnow[i] + LocalVeg->NFfrac*(*Tile)[NorthFacing].IntSnow[i] + LocalVeg->SFfrac*(*Tile)[SouthFacing].IntSnow[i];
  }
  /* ET */
  for (i = 0; i <= NVeg; i++) {
    LocalEvap->EPot[i] =
      (1 - (LocalVeg->FORfrac + LocalVeg->NFfrac + LocalVeg->SFfrac))*(*Tile)[Exposed].EPot[i] + LocalVeg->FORfrac*(*Tile)[ForestTile].EPot[i] + LocalVeg->NFfrac*(*Tile)[NorthFacing].EPot[i] + LocalVeg->SFfrac*(*Tile)[SouthFacing].EPot[i];
    LocalEvap->EAct[i] =
      (1 - (LocalVeg->FORfrac + LocalVeg->NFfrac + LocalVeg->SFfrac))*(*Tile)[Exposed].EAct[i] + LocalVeg->FORfrac*(*Tile)[ForestTile].EAct[i] + LocalVeg->NFfrac*(*Tile)[NorthFacing].EAct[i] + LocalVeg->SFfrac*(*Tile)[SouthFacing].EAct[i];
  }

  for (i = 0; i < NVeg; i++) {
    LocalEvap->EInt[i] =
      (1 - (LocalVeg->FORfrac + LocalVeg->NFfrac + LocalVeg->SFfrac))*(*Tile)[Exposed].EInt[i] + LocalVeg->FORfrac*(*Tile)[ForestTile].EInt[i] + LocalVeg->NFfrac*(*Tile)[NorthFacing].EInt[i] + LocalVeg->SFfrac*(*Tile)[SouthFacing].EInt[i];
  }

  for (i = 0; i < NVeg; i++) {
    for (j = 0; j < NSoil; j++) {
      LocalEvap->ESoil[i][j] =
        (1 - (LocalVeg->FORfrac + LocalVeg->NFfrac + LocalVeg->SFfrac))*(*Tile)[Exposed].ESoil[i][j] + LocalVeg->FORfrac*(*Tile)[ForestTile].ESoil[i][j] + LocalVeg->NFfrac*(*Tile)[NorthFacing].ESoil[i][j] + LocalVeg->SFfrac*(*Tile)[SouthFacing].ESoil[i][j];
    }
  }



  LocalEvap->ETot =
    (1 - (LocalVeg->FORfrac + LocalVeg->NFfrac + LocalVeg->SFfrac))*(*Tile)[Exposed].ETot + LocalVeg->FORfrac*(*Tile)[ForestTile].ETot + LocalVeg->NFfrac*(*Tile)[NorthFacing].ETot + LocalVeg->SFfrac*(*Tile)[SouthFacing].ETot;

}
