/***********************************************************************
/
/  GRID CLASS (COMPUTE TIME STEP)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1: 2010 Tom Abel, added MHD part 
/
/  PURPOSE:
/
/  RETURNS:
/    dt   - timestep
/
************************************************************************/
 
// Compute the timestep from all the constrains for this grid.
//
// Somebody fix the error handling in this routine! please.
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

#include <float.h>   /* for DBL_MAX */

/* function prototypes */
void find_min_dt_cell(
  long long *rank,
  long long *nx,  long long *ny,  long long *nz,
  long long *i1,  long long *i2,
  long long *j1,  long long *j2,
  long long *k1,  long long *k2,
  long long *hydroMethod,
  double     C2,
  double    *dx,
  double    *dy,
  double    *dz,
  double    *vgx,
  double    *vgy,
  double    *vgz,
  double     gamma,
  long long *ipfree,
  double     aye,
  double    *dens,
  double    *pres,
  double    *u,
  double    *v,
  double    *w,
  double    *out_rho,
  double    *out_cs,
  double    *out_dx,
  double    *out_dt
);

int CosmologyComputeExpansionTimestep(FLOAT time, float *dtExpansion);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
extern "C" void PFORTRAN_NAME(calc_dt)(
                  int *rank, int *idim, int *jdim, int *kdim,
                  int *i1, int *i2, int *j1, int *j2, int *k1, int *k2,
			     hydro_method *ihydro, float *C2,
                  FLOAT *dx, FLOAT *dy, FLOAT *dz, float *vgx, float *vgy,
                             float *vgz, float *gamma, int *ipfree, float *aye,
                  float *d, float *p, float *u, float *v, float *w,
			     float *dt, float *dtviscous);
 
extern "C" void
FORTRAN_NAME(mhd_dt)(float *bxc, float *byc, float *bzc,
                     float *vx, float *vy, float *vz,
                     float *d, float *p, float *gamma, float *dt,
                     FLOAT *dx, FLOAT *dy, FLOAT *dz,
                     int *idim, int *jdim, int *kdim, int * rank,
                     int *i1, int *i2,
                     int *j1, int *j2,
                     int *k1, int *k2, float* eng);
 
float grid::ComputeTimeStep()
{
 
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return huge_number;
 
  this->DebugCheck("ComputeTimeStep");
 
  /* initialize */
 
  float dt, dtTemp;
  float dtBaryons      = huge_number;
  float dtViscous      = huge_number;
  float dtParticles    = huge_number;
  float dtExpansion    = huge_number;
  float dtAcceleration = huge_number;
  float dtMHD          = huge_number;
  float dtConduction   = huge_number;
  float dtCR           = huge_number;
  float dtGasDrag      = huge_number;
  float dtCooling      = huge_number;
  int dim, i, j, k, index, result;
  float dtQuantum        = huge_number;  //FDM

  float TemperatureUnits, DensityUnits, LengthUnits, 
    VelocityUnits, TimeUnits, aUnits = 1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
  }

  /* Compute the field size. */
 
  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* If using comoving coordinates, compute the expansion factor a.  Otherwise,
     set it to one. */
 
  FLOAT a = 1, dadt;
  if (ComovingCoordinates)
    CosmologyComputeExpansionFactor(Time, &a, &dadt);
  float afloat = float(a);
 
  /* 1) Compute Courant condition for baryons. */
 
  if (NumberOfBaryonFields > 0 && (HydroMethod != HD_RK) && (HydroMethod != MHD_RK)) {
 
    /* Find fields: density, total energy, velocity1-3. */

    int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, 
        B1Num, B2Num, B3Num, PhiNum;
    this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, 
                                     TENum, B1Num, B2Num, B3Num, PhiNum); 

    /* For one-zone free-fall test, just compute free-fall time. */
    if (ProblemType == 63) {
      float *force_factor = new float[size];
      if (this->ComputeOneZoneCollapseFactor(force_factor) == FAIL) {
	ENZO_FAIL("Error in ComputeOneZoneCollapseFactor.\n");
      }

      dt = huge_number;
      for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) { // nothing
	for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) { // metallicity
	  for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) { // energy

	    index = i + j*GridDimension[0] + k*GridDimension[0]*GridDimension[1];

	    dt = min(dt, POW(((3 * pi) / 
			      (32 * GravitationalConstant * 
			       BaryonField[DensNum][index] *
			       (1 - force_factor[index]))), 0.5));
	  }
	}
      }

      delete [] force_factor;
      dt *= TestProblemData.OneZoneFreefallTimestepFraction;
      return dt;
    }
 
    /* Compute the pressure. */
 
    float *pressure_field = new float[size];
    this->ComputePressure(Time, pressure_field,0,1); // Note: Force use of CRs to get sound speed correct
 
#ifdef UNUSED
    int Zero[3] = {0,0,0}, TempInt[3] = {0,0,0};
    for (dim = 0; dim < GridRank; dim++)
      TempInt[dim] = GridDimension[dim]-1;
#endif /* UNUSED */
 
    /* Call fortran routine to do calculation. */
 
    if( HydroMethod != MHD_Li)
      PFORTRAN_NAME(calc_dt)(&GridRank, GridDimension, GridDimension+1,
			     GridDimension+2,
			     GridStartIndex, GridEndIndex,
			     GridStartIndex+1, GridEndIndex+1,
			     GridStartIndex+2, GridEndIndex+2,
			     &HydroMethod, &ZEUSQuadraticArtificialViscosity,
			     CellWidth[0], CellWidth[1], CellWidth[2],
			     GridVelocity, GridVelocity+1, GridVelocity+2,
			     &Gamma, &PressureFree, &afloat,
			     BaryonField[DensNum], pressure_field,
			     BaryonField[Vel1Num], BaryonField[Vel2Num],
			     BaryonField[Vel3Num], &dtBaryons, &dtViscous);
 
    if ((HydroMethod != MHD_Li) & (dtBaryons*CourantSafetyNumber*TimeUnits/yr_s < 100)){
            double rho_min, cs_min, dx_min, dt_min;
            find_min_dt_cell(&GridRank, GridDimension, GridDimension+1,
                              GridDimension+2,
                              GridStartIndex, GridEndIndex,
                              GridStartIndex+1, GridEndIndex+1,
                              GridStartIndex+2, GridEndIndex+2,
                              &HydroMethod, ZEUSQuadraticArtificialViscosity,
                              CellWidth[0], CellWidth[1], CellWidth[2],
                              GridVelocity, GridVelocity+1, GridVelocity+2,
                              Gamma, &PressureFree, afloat,
                              BaryonField[DensNum], pressure_field,
                              BaryonField[Vel1Num], BaryonField[Vel2Num],
                              BaryonField[Vel3Num], 
                               &rho_min, &cs_min, &dx_min, &dt_min);
            printf("Cell with min(dt): dt = %e, rho = %g, cs = %g, dx = %g\n",
                   dt_min, rho_min, cs_min, dx_min);
          }

    if(HydroMethod == MHD_Li){
      /* 1.5) Calculate minimum dt due to MHD: Maximum Fast MagnetoSonic Shock Speed */
      
      //Cosmos nees this, for some reason.
      if(GridRank < 3 ){
	if( CellWidth[2] == NULL ) CellWidth[2] = new FLOAT;
	CellWidth[2][0] = 1.0;
	if( GridRank < 2 ){
	  if( CellWidth[1] == NULL ) CellWidth[1] = new FLOAT;
	  CellWidth[1][0] = 1.0;
	}
      }
      int Rank_Hack = 3; //MHD needs a 3d timestep always.
      FORTRAN_NAME(mhd_dt)(BaryonField[B1Num], BaryonField[B2Num], BaryonField[B3Num],
			   BaryonField[Vel1Num], BaryonField[Vel2Num], BaryonField[Vel3Num],
			   BaryonField[DensNum], pressure_field, &Gamma, &dtMHD, 
			   CellWidth[0], CellWidth[1], CellWidth[2],
			   GridDimension, GridDimension + 1, GridDimension +2, &Rank_Hack,
			   GridStartIndex, GridEndIndex,
			   GridStartIndex+1, GridEndIndex+1,
			   GridStartIndex+2, GridEndIndex+2, BaryonField[TENum]);
      
      dtMHD *= CourantSafetyNumber;
      dtMHD *= afloat;  
    }//if HydroMethod== MHD_Li

    /* Clean up */
 
    delete [] pressure_field;
 
    /* Multiply resulting dt by CourantSafetyNumber (for extra safety!). */
 
    dtBaryons *= CourantSafetyNumber;
 
  }


  if (NumberOfBaryonFields > 0 && 
      (HydroMethod == HD_RK)) {

    int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
    if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
					 Vel3Num, TENum) == FAIL) {
      fprintf(stderr, "ComputeTimeStep: IdentifyPhysicalQuantities error.\n");
      exit(FAIL);
    }

    FLOAT dxinv = 1.0 / CellWidth[0][0]/a;
    FLOAT dyinv = (GridRank > 1) ? 1.0 / CellWidth[1][0]/a : 0.0;
    FLOAT dzinv = (GridRank > 2) ? 1.0 / CellWidth[2][0]/a : 0.0;
    float dt_temp = 1.e-20, dt_ltemp, dt_x, dt_y, dt_z;
    float rho, p, vx, vy, vz, v2, eint, etot, h, cs, dpdrho, dpde,
      v_signal_x, v_signal_y, v_signal_z;
    int n = 0;
    for (k = 0; k < GridDimension[2]; k++) {
      for (j = 0; j < GridDimension[1]; j++) {
	for (i = 0; i < GridDimension[0]; i++, n++) {
	  rho = BaryonField[DensNum][n];
	  vx  = BaryonField[Vel1Num][n];
	  vy  = BaryonField[Vel2Num][n];
	  vz  = BaryonField[Vel3Num][n];

	  if (DualEnergyFormalism) {
	    eint = BaryonField[GENum][n];
	  }
	  else {
	    etot = BaryonField[TENum][n];
	    v2 = vx*vx + vy*vy + vz*vz;
	    eint = etot - 0.5*v2;
	  }

	  EOS(p, rho, eint, h, cs, dpdrho, dpde, EOSType, 2);

	  v_signal_y = v_signal_z = 0;

	  v_signal_x = (cs + fabs(vx));
	  if (GridRank > 1) v_signal_y = (cs + fabs(vy));
	  if (GridRank > 2) v_signal_z = (cs + fabs(vz));

	  dt_x = v_signal_x * dxinv;
	  dt_y = v_signal_y * dyinv;
	  dt_z = v_signal_z * dzinv;

	  dt_ltemp = my_MAX(dt_x, dt_y, dt_z);

	  if (dt_ltemp > dt_temp) {
	    dt_temp = dt_ltemp;
	  }
        }
      }
    }

    dtBaryons = CourantSafetyNumber / dt_temp;

  }


  // MHD
  if (NumberOfBaryonFields > 0 && HydroMethod == MHD_RK) {

    int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, 
      B1Num, B2Num, B3Num, PhiNum, CRNum;
    if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
					 Vel3Num, TENum, B1Num, B2Num, B3Num, PhiNum) == FAIL)
      ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
    if (CRModel){
      if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
					   Vel3Num, TENum, CRNum) == FAIL) {
	ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
      }
    }
    FLOAT dxinv = 1.0 / CellWidth[0][0]/a;
    FLOAT dyinv = (GridRank > 1) ? 1.0 / CellWidth[1][0]/a : 0.0;
    FLOAT dzinv = (GridRank > 2) ? 1.0 / CellWidth[2][0]/a : 0.0;
    float vxm, vym, vzm, Bm, rhom;
    float dt_temp = 1.e-20, dt_ltemp, dt_x, dt_y, dt_z;
    float rho, p, vx, vy, vz, v2, eint, etot, h, cs, cs2, dpdrho, dpde,
      v_signal_x, v_signal_y, v_signal_z, cf, cf2, temp1, Bx, By, Bz, B2, ca2, Pcr;
    int n = 0;
    float rho_dt, B_dt, v_dt;
    for (k = 0; k < GridDimension[2]; k++) {
      for (j = 0; j < GridDimension[1]; j++) {
	for (i = 0; i < GridDimension[0]; i++, n++) {
	  rho = BaryonField[DensNum][n];
	  vx  = BaryonField[Vel1Num][n];
	  vy  = BaryonField[Vel2Num][n];
	  vz  = BaryonField[Vel3Num][n];
	  Bx  = BaryonField[B1Num][n];
	  By  = BaryonField[B2Num][n];
	  Bz  = BaryonField[B3Num][n];
          if (CRModel)
            Pcr = (CRgamma - 1.0) * BaryonField[CRNum][n];

          B2 = Bx*Bx + By*By + Bz*Bz;
	  if (DualEnergyFormalism) {
	    eint = BaryonField[GENum][n];
	  }
	  else {
	    etot = BaryonField[TENum][n];
	    v2 = vx*vx + vy*vy + vz*vz;
	    eint = etot - 0.5*v2 - 0.5*B2/rho;
	  }

	  v_signal_y = v_signal_z = 0;

	  EOS(p, rho, eint, h, cs, dpdrho, dpde, EOSType, 2);
	  cs2 = cs*cs;
          if (CRModel){
            cs2 += CRgamma*Pcr/rho;
          }

	  temp1 = cs2 + B2/rho;

	  ca2 = Bx*Bx/rho;
	  cf2 = 0.5 * (temp1 + sqrt(temp1*temp1 - 4.0*cs2*ca2));
	  cf = sqrt(cf2);
	  v_signal_x = (cf + fabs(vx));

	  if (GridRank > 1) {
	    ca2 = By*By/rho;
	    cf2 = 0.5 * (temp1 + sqrt(temp1*temp1 - 4.0*cs2*ca2));
	    cf = sqrt(cf2);
	    v_signal_y = (cf + fabs(vy));
	  }

	  if (GridRank > 2) {
	    ca2 = Bz*Bz/rho;
	    cf2 = 0.5 * (temp1 + sqrt(temp1*temp1 - 4.0*cs2*ca2));
	    cf = sqrt(cf2);
	    v_signal_z = (cf + fabs(vz));
	  }

	  dt_x = v_signal_x * dxinv;
	  dt_y = v_signal_y * dyinv;
	  dt_z = v_signal_z * dzinv;

	  dt_ltemp = my_MAX(dt_x, dt_y, dt_z);

	  if (dt_ltemp > dt_temp) {
	    dt_temp = dt_ltemp;
	    rho_dt = rho;
	    B_dt = sqrt(Bx*Bx+By*By+Bz*Bz);
	    v_dt = max(fabs(vx), fabs(vy));
	    v_dt = max(v_dt, fabs(vz));
	  }

        }
      }
    }
    dtMHD = CourantSafetyNumber / dt_temp;
    //    fprintf(stderr, "ok %g %g %g\n", dt_x,dt_y,dt_z);
    //    if (dtMHD*TimeUnits/yr < 5) {
    //float ca = B_dt/sqrt(rho_dt)*VelocityUnits;
    //printf("dt=%g, rho=%g, B=%g\n, v=%g, ca=%g, dt=%g", dtMHD*TimeUnits/yr, rho_dt*DensityUnits, B_dt*MagneticUnits, 
    //    v_dt*VelocityUnits/1e5, ca/1e5, LengthUnits/(dxinv*ca)/yr*CourantSafetyNumber);
    // }

  } // HydroMethod = MHD_RK

 
  /* 2) Calculate dt from particles. */
 
  if (NumberOfParticles > 0 || NumberOfActiveParticles > 0) {
 
    /* Compute dt constraint from particle velocities. */
 
    for (dim = 0; dim < GridRank; dim++) {
      float dCell = CellWidth[dim][0]*a;
      for (i = 0; i < NumberOfParticles; i++) {
        dtTemp = dCell/max(fabs(ParticleVelocity[dim][i]), tiny_number);
	dtParticles = min(dtParticles, dtTemp);
      }
      for (i = 0; i < NumberOfActiveParticles; i++) {
    dtTemp = dCell/max(fabs(ActiveParticles[i]->ReturnVelocity()[dim]), tiny_number);
    dtParticles = min(dtParticles, dtTemp);
      }
    }
 
    /* Multiply resulting dt by ParticleCourantSafetyNumber. */
 
    dtParticles *= ParticleCourantSafetyNumber;
 
  }
 

  /* 3) Find dt from expansion. */
 
  if (ComovingCoordinates)
    if (CosmologyComputeExpansionTimestep(Time, &dtExpansion) == FAIL) {
      fprintf(stderr, "nudt: Error in ComputeExpansionTimestep.\n");
      exit(FAIL);
    }
 
  /* 4) Calculate minimum dt due to acceleration field (if present). */
 
  if (SelfGravity) {
    for (dim = 0; dim < GridRank; dim++)
      if (AccelerationField[dim])
	for (i = 0; i < size; i++) {
	  dtTemp = sqrt(CellWidth[dim][0]/
			fabs(AccelerationField[dim][i])+tiny_number);
	  dtAcceleration = min(dtAcceleration, dtTemp);
	}
    if (dtAcceleration != huge_number)
      dtAcceleration *= 0.5;
  }

  /* 5) Calculate minimum dt due to thermal conduction. */

  if(IsotropicConduction || AnisotropicConduction){
    if (this->ComputeConductionTimeStep(dtConduction) == FAIL) 
      ENZO_FAIL("Error in ComputeConductionTimeStep.\n");

    dtConduction *= float(NumberOfGhostZones);     // for subcycling 
  }
  
  /* 6) Calculate minimum dt due to CR diffusion */

  if(CRModel){
    if(CRDiffusion){
      if( this->ComputeCRDiffusionTimeStep(dtCR) == FAIL) {
	fprintf(stderr, "Error in ComputeCRDiffusionTimeStep.\n");
        return FAIL;
      }
    }
    if(CRStreaming){
      if( this->ComputeCRStreamingTimeStep(dtCR) == FAIL) {
        fprintf(stderr, "Error in ComputeCRStreamingTimeStep.\n");
        return FAIL;
      }
    }
    dtCR *= CRCourantSafetyNumber;
    if (CRDiffusion == 1)
      dtCR *= float(NumberOfGhostZones); // for subcycling
  }

  /* 7) GasDrag time step */
  if (UseGasDrag && GasDragCoefficient != 0.) {
    dtGasDrag = 0.5/GasDragCoefficient;
  }


  /* Cooling time */
  
  if (UseCoolingTimestep == TRUE) {
    float *cooling_time = new float[size];
    if (this->ComputeCoolingTime(cooling_time, TRUE) == FAIL) {
      ENZO_FAIL("Error in grid->ComputeCoolingTime.\n");
    }

    for (k = GridStartIndex[2]; k < GridEndIndex[2]; k++) {
      for (j = GridStartIndex[1]; j < GridEndIndex[1]; j++) {
	index = GRIDINDEX_NOGHOST(GridStartIndex[0], j, k);
	for (i = GridStartIndex[0]; i < GridEndIndex[0]; i++, index++) {
	  dtCooling = min(dtCooling, cooling_time[index]);
	}
      }
    }
    dtCooling *= CoolingTimestepSafetyFactor;
 
    delete [] cooling_time;
  }

   /* FDM: Calculate minimum dt due to quantum pressure. */

  if(QuantumPressure){

    double hmcoef = 5.9157166856e27*TimeUnits/POW(LengthUnits/afloat,2)/FDMMass; // 5.916e27 is hbar/m with m=1e-22eV, FDMMass is in unit of 1e-22eV.

    FLOAT dx = CellWidth[0][0]*afloat;

    if (GridRank>1)
      dx = min( dx, CellWidth[1][0]*afloat);
  
    if (GridRank>2)
      dx = min( dx, CellWidth[2][0]*afloat);

    dtQuantum = POW(dx,2)/hmcoef/2.;

    dtQuantum *= CourantSafetyNumber;

    if (SelfGravity && (PotentialField != NULL)){
      int gsize = GravitatingMassFieldDimension[0]*GravitatingMassFieldDimension[1]*GravitatingMassFieldDimension[2];

      for (int i=0; i<gsize; ++i){
	dtQuantum = min(dtQuantum, fabs(hmcoef*1./(PotentialField[i])));
      }
    }

  }

  /* 8) calculate minimum timestep */
 
  dt = min(dtBaryons, dtParticles);
  dt = min(dt, dtMHD);
  dt = min(dt, dtViscous);
  dt = min(dt, dtAcceleration);
  dt = min(dt, dtExpansion);
  dt = min(dt, dtConduction);
  dt = min(dt, dtCR);
  dt = min(dt, dtGasDrag);
  dt = min(dt, dtCooling);
  dt = min(dt, dtQuantum); //FDM
  
#ifdef TRANSFER

  /* 8) If using radiation pressure, calculate minimum dt */

  float dtRadPressure = huge_number;
  float absVel, absAccel;

  if (RadiationPressure && RadiativeTransfer) {

    int RPresNum1, RPresNum2, RPresNum3;
    if (IdentifyRadiationPressureFields(RPresNum1, RPresNum2, RPresNum3) 
	== FAIL) {
      ENZO_FAIL("Error in IdentifyRadiationPressureFields.");
    }

    for (i = 0; i < size; i++)
      for (dim = 0; dim < GridRank; dim++) {
	dtTemp = sqrt(CellWidth[dim][0] / (fabs(BaryonField[RPresNum1+dim][i])+
					   tiny_number));
	dtRadPressure = min(dtRadPressure, dtTemp);
      }
    
    if (dtRadPressure < huge_number)
      dtRadPressure *= 0.5;

    dt = min(dt, dtRadPressure);

  } /* ENDIF RadiationPressure */

  /* 9) Safety Velocity to limit timesteps */

  float dtSafetyVelocity = huge_number;
  if (TimestepSafetyVelocity > 0)
    dtSafetyVelocity = a*CellWidth[0][0] / 
      (TimestepSafetyVelocity*1e5 / VelocityUnits);    // parameter in km/s

  dt = min(dt, dtSafetyVelocity);


  /* 10) FLD Radiative Transfer timestep limitation */
  if (RadiativeTransferFLD)
    dt = min(dt, MaxRadiationDt);
  
#endif /* TRANSFER */
 
  /* Debugging info. */
  
  if (debug1) {
    printf("ComputeTimeStep = %"ESYM" (", dt);
    if (HydroMethod != MHD_RK && HydroMethod != MHD_Li && NumberOfBaryonFields > 0)
      printf("Bar = %"ESYM" ", dtBaryons);
    if (HydroMethod == MHD_RK || HydroMethod == MHD_Li)
      printf("dtMHD = %"ESYM" ", dtMHD);
    if (CRModel)
      printf("dtCR = %"ESYM" ", dtCR);
    if (HydroMethod == Zeus_Hydro)
      printf("Vis = %"ESYM" ", dtViscous);
    if (ComovingCoordinates)
      printf("Exp = %"ESYM" ", dtExpansion);
    if (dtAcceleration != huge_number)
      printf("Acc = %"ESYM" ", dtAcceleration);
    if (NumberOfParticles)
      printf("Part = %"ESYM" ", dtParticles);
    if (UseCoolingTimestep)
      printf("Cool = %"ESYM" ", dtCooling);
    if (IsotropicConduction || AnisotropicConduction)
      printf("Cond = %"ESYM" ",(dtConduction));
    if (UseGasDrag)
      printf("Drag = %"ESYM" ",(dtGasDrag));
    if (QuantumPressure)
      printf("Quantum = %"ESYM" ",(dtQuantum));//FDM
    printf(")\n");
  }
 
  return dt;
}



void find_min_dt_cell(
  long long *rank,
  long long *nx,  long long *ny,  long long *nz,
  long long *i1,  long long *i2,
  long long *j1,  long long *j2,
  long long *k1,  long long *k2,
  long long *hydroMethod,
  double     C2,
  double    *dx,
  double    *dy,
  double    *dz,
  double    *vgx,
  double    *vgy,
  double    *vgz,
  double     gamma,
  long long *ipfree,
  double     aye,
  double    *dens,
  double    *pres,
  double    *u,
  double    *v,
  double    *w,
  double    *out_rho,
  double    *out_cs,
  double    *out_dx,
  double    *out_dt
) {
  const double tiny = 1e-20;
  double best_dt = DBL_MAX;
  long long best_i = *i1, best_j = *j1, best_k = *k1;

  // Bulk‐velocity values (dereference pointers):
  double vgx_val = (vgx != NULL ? *vgx : 0.0);
  double vgy_val = (vgy != NULL ? *vgy : 0.0);
  double vgz_val = (vgz != NULL ? *vgz : 0.0);

  // Loop over active subregion [i1..i2] × [j1..j2] × [k1..k2]
  for (long long kk = *k1;  kk <= *k2;  ++kk) {
      for (long long jj = *j1;  jj <= *j2;  ++jj) {
          for (long long ii = *i1;  ii <= *i2;  ++ii) {
              // Compute flattened index
              long long idx = ii
                            + jj * (*nx)
                            + kk * ((*nx) * (*ny));
              double rho = dens[idx];
              double P   = pres[idx];
              // Skip cells with non-positive density or pressure
              if (!(rho > 0.0) || !(P > 0.0)) {
                  continue;
              }
              // Sound speed
              double cs_local = sqrt(gamma * P / rho);
              if (*ipfree == 1) {
                  cs_local = tiny;
              }
              // Velocity differences
              double du = fabs(u[idx] - vgx_val);
              double dv = fabs(v[idx] - vgy_val);
              double dw = 0.0;
              if (*rank == 3 && w != NULL) {
                  dw = fabs(w[idx] - vgz_val);
              }
              // Build the denominator = (cs+|du|)/dx[ii] + (cs+|dv|)/dy[jj] + (if 3D) (cs+|dw|)/dz[kk]
              double denom = (cs_local + du) / dx[ii]
                           + (cs_local + dv) / dy[jj];
              if (*rank == 3 && w != NULL) {
                  denom += (cs_local + dw) / dz[kk];
              }
              if (denom <= tiny) {
                  continue;
              }
              double dt_cell = aye / denom;

              // If you want the ZEUS‐viscous term when *hydroMethod == 2, you can insert it here.
              // For pure Godunov‐style dt, we leave it out.

              if (dt_cell < best_dt) {
                  best_dt = dt_cell;
                  best_i  = ii;
                  best_j  = jj;
                  best_k  = kk;
              }
          }
      }
  }

  // If we found a valid “best” cell, fill outputs
  if (best_dt < DBL_MAX) {
      long long idxm = best_i
                     + best_j * (*nx)
                     + best_k * ((*nx) * (*ny));
      double rho_m = dens[idxm];
      double P_m   = pres[idxm];
      double cs_m  = sqrt(gamma * P_m / rho_m);
      double T_m   = P_m / rho_m;  // “temperature” proxy

      *out_rho = rho_m;
      *out_cs   = cs_m;
      *out_dx   = dx[best_i];
      *out_dt   = best_dt;
  } else {
      // No valid cell → zero out results
      *out_rho = 0.0;
      *out_cs   = 0.0;
      *out_dx   = 0.0;
      *out_dt   = 0.0;
  }
}