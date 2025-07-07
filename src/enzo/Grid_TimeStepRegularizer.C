/***********************************************************************
/
/  GRID CLASS (REGULARIZE TIME STEP)
/
/  written by: Greg Bryan
/  date:       May 2025
/  modified1: 
/
/  PURPOSE:
/
/  RETURNS:
/    
/
************************************************************************/
 
// Compute the timestep from all the constrains for this grid.
//

 
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "RadiativeTransferParameters.h"
#include "hydro_rk/EOS.h"
#include "hydro_rk/tools.h"
#include "phys_constants.h"
 
/* function prototypes */
 
int CosmologyComputeExpansionTimestep(FLOAT time, float *dtExpansion);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
 
int grid::TimeStepRegularizer(double dtMinimum)
{
 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber || NumberOfBaryonFields == 0)
    return SUCCESS;
 
  this->DebugCheck("RegularizeTimeStep");
 
  /* initialize */
 
  float TemperatureUnits, DensityUnits, LengthUnits, 
    VelocityUnits, TimeUnits, aUnits = 1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
  }

  /* Compute the field size. */
 
  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* If using comoving coordinates, compute the expansion factor a.  Otherwise,
     set it to one. */
 
  FLOAT a = 1, dadt;
  if (ComovingCoordinates)
    CosmologyComputeExpansionFactor(Time, &a, &dadt);
 
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "ComputeTimeStep: IdentifyPhysicalQuantities error.\n");
    exit(FAIL);
  }

  FLOAT dx = 1.0/CellWidth[0][0]/a;
  double rho, vx, vy, vz, v2, eint, etot, cs, dt_cs, dt_vel, rho_new;
  double rho_max, vx_max, vy_max, vz_max, v2_max, eint_max, etot_max;
  double delta_rho, eint_new, vel_new, vel, vel_max;
  int dim, i, j, k, n, imax, ioffset;
  
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      n = GridStartIndex[0] + (j + k*GridDimension[1])*GridDimension[0];
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
	rho = BaryonField[DensNum][n];
	vx  = BaryonField[Vel1Num][n];
	vy  = BaryonField[Vel2Num][n];
	vz  = BaryonField[Vel3Num][n];
	v2 = vx*vx + vy*vy + vz*vz;
	
	if (DualEnergyFormalism) {
	  eint = BaryonField[GENum][n];
	}
	else {
	  etot = BaryonField[TENum][n];
	  eint = etot - 0.5*v2;
	}

	/* Compute thermal signal speed and resulting timestep constraints,
	    split by thermal and advective. */
	
	cs = sqrt(max(eint, tiny_number));
	dt_cs = dx/cs;
	dt_vel = dx/max(sqrt(v2), tiny_number);
	
	/* Check for sound speed failure. */

	  if (min(dt_cs, dt_vel) < dtMinimum) {

	  /* Identify cell in nearest neighbours with largest density.
	     missing: this should respect grid edges but currently does not. */
	  
	  int ioffset = 1, imax = 0;
	  for (dim = 0; dim < GridRank; dim++) {
	    if (BaryonField[DensNum][n-ioffset] > BaryonField[DensNum][n+imax])
	      imax = -ioffset;
	    if (BaryonField[DensNum][n+ioffset] > BaryonField[DensNum][n+imax])
	      imax = +ioffset;
	    ioffset *= GridDimension[dim];
	  }
	  rho_max = BaryonField[DensNum][n+imax];
	  vx_max  = BaryonField[Vel1Num][n+imax];
	  vy_max  = BaryonField[Vel2Num][n+imax];
	  vz_max  = BaryonField[Vel3Num][n+imax];
	  v2_max = vx_max*vx_max + vy_max*vy_max + vz_max*vz_max;
	  if (DualEnergyFormalism) {
	    eint_max = BaryonField[GENum][n+imax];
	  }
	  else {
	    etot_max = BaryonField[TENum][n+imax];
	    eint_max = etot_max - 0.5*v2_max;
	  }

	  /* Calculate require new density (assuming thermal energy conservation)
	     required to decrease timestep by required factor. */

	  delta_rho = 0.0;
	  if (dt_cs < dtMinimum) {
	    eint_new = eint * pow(dt_cs/dtMinimum, 2);
	    if (eint_new > eint_max) {
	      delta_rho = rho * (eint - eint_new)/(eint_new - eint_max);
	    }
	  }	    

	  if (dt_vel < dtMinimum) {
	    vel_new = vel * (dt_vel/dtMinimum);
	    if (vel_new > vel_max) {
	      delta_rho = max(rho * (vel - vel_new)/(vel_new - vel_max), delta_rho);
	    } 
	  }

	  /* Check that we are moving at most a fixed fraction of the mass. */
	  const double REGULARIZER_MAX_FACTOR = 0.5;
	  
	  if (delta_rho/rho_max > REGULARIZER_MAX_FACTOR || delta_rho == 0.0) {
	    fprintf(stderr, "DT regularizer failure, rho = %g, rho_new = %g, delta_rho = %g, imax = %d\n", rho, rho_new, delta_rho, imax);
	  } else {

	    fprintf(stderr, "DT regularizer success, rho = %g, rho_new = %g, delta_rho = %g, imax = %d ratio = %g\n", rho, rho_new, delta_rho, imax, delta_rho/rho_max);

	    /* Set conserved quantities for the two cells */

	    BaryonField[Vel1Num][n] *= rho;
	    BaryonField[Vel2Num][n] *= rho;
	    BaryonField[Vel3Num][n] *= rho;
	    eint *= rho;
	    BaryonField[Vel1Num][n+imax] *= rho_max;
	    BaryonField[Vel2Num][n+imax] *= rho_max;
	    BaryonField[Vel3Num][n+imax] *= rho_max;
	    eint_max *= rho_max;

	    /* missing: colour variables! */
	    
	    /* Move mass from old (max) cell to current cell. */
	    
	    BaryonField[DensNum][n]      += (rho_new - rho);
	    BaryonField[DensNum][n+imax] -= (rho_new - rho);

	    /* Revert to primitives. */

	    rho = BaryonField[DensNum][n];
	    rho_max = BaryonField[DensNum][n+imax];
	    
	    BaryonField[Vel1Num][n] /= rho;
	    BaryonField[Vel2Num][n] /= rho;
	    BaryonField[Vel3Num][n] /= rho;
	    eint /= rho;
	    BaryonField[Vel1Num][n+imax] /= rho_max;
	    BaryonField[Vel2Num][n+imax] /= rho_max;
	    BaryonField[Vel3Num][n+imax] /= rho_max;
	    eint_max /= rho_max;

	    if (DualEnergyFormalism) {
	      BaryonField[GENum][n] = eint;
	      BaryonField[GENum][n+imax] = eint_max;
	    } else {
	      BaryonField[TENum][n] = eint;
	      BaryonField[TENum][n+imax] = eint_max;
	      for (dim = 0; dim < GridRank; dim++) {
		BaryonField[TENum][n] += 0.5*pow(BaryonField[Vel1Num][n], 2.0);
		BaryonField[TENum][n+imax] += 0.5*pow(BaryonField[Vel1Num][n+imax], 2.0);
	      }
	    }

	    /* missing: colour variables! */

	    
	  } /* end: check for regularizer density OK */

	} /* end: check for dt too small. */
	
      }
    }
  }
  
  return SUCCESS;
}
