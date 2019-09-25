/*
 * SUMMARY:      myNewScript.c - This does awesome stuff
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Rhinoceros
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Calculate the snow albedo as a function of snow age
 * DESCRIP-END.
 * FUNCTIONS:    CalcSnowAlbedo()
 * COMMENTS:
 * $Id: CalcSnowAlbedo.c,v 1.4 2003/07/01 21:26:10 olivier Exp $     
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "settings.h"
#include "constants.h"
#include "data.h"
#include "Calendar.h"
#include "functions.h"

/*****************************************************************************
  myNewScript()

  Source:
  Laramie, R. L., and J. C. Schaake, Jr., Simulation of the continuous
  snowmelt process, Ralph M. Parsons Laboratory, Mass. Inst. of Technol.,
  1972

  Snow albedo is calculated as a function of the number of days since the
  last observed snow fall. There are separete albedo curves for the freeze
  and thaw conditions.
*****************************************************************************/
float myNewScript(float TSurf, unsigned short Last, SNOWPIX *LocalSnow, 
  int StepsPerDay)
{
/*THIS DOES EVERYTHING
*/
}
