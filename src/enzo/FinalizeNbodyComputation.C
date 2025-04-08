/***********************************************************************
	/
	/  FIND ALL NBODY PARTICLES OVER ALL PROCESSORS
	/
	/  written by: Yongseok Jo
	/  date:       November, 2022
	/
	/  PURPOSE: First synchronizes particle information in the normal and 
	/           nbody particles.  Then we make a global particle list, which
	/           simplifies Nbody calculations.
	/
 ************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
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
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "CommunicationUtilities.h"
#include "NbodyRoutines.h"  //added  
#include "phys_constants.h"


void InitializeNbodyArrays(int);
int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		HierarchyEntry **Grids[]);
int GetUnits(double *DensityUnits, double *LengthUnits,
		double *TemperatureUnits, double *TimeUnits,
		double *VelocityUnits, double *MassUnits, double Time);

#ifdef NBODY
int FinalizeNbodyComputation(LevelHierarchyEntry *LevelArray[], int level)
{


	if (LevelArray[level+1] == NULL) {
		int i, GridNum, LocalNumberOfNbodyParticles=0, NewLocalNumberOfNbodyParticles=0;
		LevelHierarchyEntry *Temp;
		int start_index, start_index_new;

		FindTotalNumberOfNbodyParticles(LevelArray, &LocalNumberOfNbodyParticles, &NewLocalNumberOfNbodyParticles);

		if (NumberOfNbodyParticles == 0 && NumberOfNewNbodyParticles == 0)
			return SUCCESS;

		double *NbodyParticlePositionTemp[MAX_DIMENSION];
		double *NbodyParticleVelocityTemp[MAX_DIMENSION];
#ifdef SEVN
		double *NbodyParticleInitialMassTemp;
		double *NbodyParticleWindEjectedMassTemp;
		double *NbodyParticleSNEjectedMassTemp;
		double *NbodyParticleTemperatureTemp;

		double *NbodyParticleMassTemp;
#endif

		double *NewNbodyParticlePositionTemp[MAX_DIMENSION];
		double *NewNbodyParticleVelocityTemp[MAX_DIMENSION];

#ifdef SEVN
		double *NewNbodyParticleInitialMassTemp;
		double *NewNbodyParticleWindEjectedMassTemp;
		double *NewNbodyParticleSNEjectedMassTemp;
		double *NewNbodyParticleTemperatureTemp;

		double *NewNbodyParticleMassTemp;
#endif

		for (int dim=0; dim<MAX_DIMENSION; dim++) {
			NbodyParticlePositionTemp[dim]            = new double[LocalNumberOfNbodyParticles];
			NbodyParticleVelocityTemp[dim]            = new double[LocalNumberOfNbodyParticles];
		}
#ifdef SEVN
		NbodyParticleInitialMassTemp			= new double[LocalNumberOfNbodyParticles];
		NbodyParticleWindEjectedMassTemp		= new double[LocalNumberOfNbodyParticles];
		NbodyParticleSNEjectedMassTemp			= new double[LocalNumberOfNbodyParticles];
		NbodyParticleTemperatureTemp			= new double[LocalNumberOfNbodyParticles];

		NbodyParticleMassTemp					= new double[LocalNumberOfNbodyParticles];
#endif

		for (int dim=0; dim<MAX_DIMENSION; dim++) {
			NewNbodyParticlePositionTemp[dim]            = new double[NewLocalNumberOfNbodyParticles];
			NewNbodyParticleVelocityTemp[dim]            = new double[NewLocalNumberOfNbodyParticles];
		}
#ifdef SEVN
		NewNbodyParticleInitialMassTemp			= new double[NewLocalNumberOfNbodyParticles];
		NewNbodyParticleWindEjectedMassTemp		= new double[NewLocalNumberOfNbodyParticles];
		NewNbodyParticleSNEjectedMassTemp		= new double[NewLocalNumberOfNbodyParticles];
		NewNbodyParticleTemperatureTemp			= new double[NewLocalNumberOfNbodyParticles];

		NewNbodyParticleMassTemp				= new double[NewLocalNumberOfNbodyParticles];
#endif


		/* Find the index of the array */
		start_index = FindStartIndex(&LocalNumberOfNbodyParticles);
		start_index_new = FindStartIndex(&NewLocalNumberOfNbodyParticles);
		//fprintf(stderr,"NumberOfParticles=%d in FINAL\n",NumberOfNbodyParticles);
		//fprintf(stderr,"NewNumberOfParticles=%d in FINAL\n",NumberOfNewNbodyParticles);


#ifdef USE_MPI
		if (MyProcessorNumber == ROOT_PROCESSOR) {

			double *NewNbodyParticleVelocity[MAX_DIMENSION]; // feedback can affect velocity
			double *NewNbodyParticlePosition[MAX_DIMENSION]; // feedback can affect velocity
#ifdef SEVN
			double *NewNbodyParticleInitialMass;
			double *NewNbodyParticleWindEjectedMass;
			double *NewNbodyParticleSNEjectedMass;
			double *NewNbodyParticleTemperature;

			double *NewNbodyParticleMass;
#endif

			for (int dim=0; dim<MAX_DIMENSION; dim++) {
				NewNbodyParticlePosition[dim]            = new double[NumberOfNewNbodyParticles];
				NewNbodyParticleVelocity[dim]            = new double[NumberOfNewNbodyParticles];
			}
#ifdef SEVN
			NewNbodyParticleInitialMass				= new double[NumberOfNewNbodyParticles];
			NewNbodyParticleWindEjectedMass			= new double[NumberOfNewNbodyParticles];
			NewNbodyParticleSNEjectedMass			= new double[NumberOfNewNbodyParticles];
			NewNbodyParticleTemperature				= new double[NumberOfNewNbodyParticles];

			NewNbodyParticleMass					= new double[NumberOfNewNbodyParticles];
#endif

			/* Receiving Index, NumberOfParticles, NbodyArrays from other processs */
			int* start_index_all;
			int* LocalNumberAll;
			start_index_all = new int[NumberOfProcessors];
			LocalNumberAll = new int[NumberOfProcessors];

			int* start_index_all_new;
			int* NewLocalNumberAll;
			start_index_all_new = new int[NumberOfProcessors];
			NewLocalNumberAll = new int[NumberOfProcessors];

			MPI_Request request;
			MPI_Status status;
			int ierr;
			int errclass,resultlen;
			char err_buffer[MPI_MAX_ERROR_STRING];

			if (NumberOfNbodyParticles != 0) {
				MPI_Gather(&LocalNumberOfNbodyParticles, 1, IntDataType, LocalNumberAll, 1, IntDataType, ROOT_PROCESSOR, enzo_comm);
				MPI_Gather(&start_index, 1, IntDataType, start_index_all, 1, IntDataType, ROOT_PROCESSOR, enzo_comm);
			}

			if (NumberOfNewNbodyParticles != 0) {
				MPI_Gather(&NewLocalNumberOfNbodyParticles, 1, IntDataType, NewLocalNumberAll, 1, IntDataType, ROOT_PROCESSOR, enzo_comm);
				MPI_Gather(&start_index_new, 1, IntDataType, start_index_all_new, 1, IntDataType, ROOT_PROCESSOR, enzo_comm);
			}



			/*-----------------------------------------------*/
			/******** Recv Arrays from Nbody+    *****/
			/*-----------------------------------------------*/
			//fprintf(stderr,"NumberOfParticles=%d in Final\n",NumberOfNbodyParticles);


			fprintf(stdout, "ENZO: Waiting for NBODY+ to receive data \n");
			//fprintf(stderr, "ENZO: NNB    = %d \n", NumberOfNbodyParticles);
			//fprintf(stderr, "ENZO: newNNB = %d \n", NumberOfNewNbodyParticles);
			InitializeNbodyArrays(1);
			CommunicationInterBarrier();
			//fprintf(stderr,"NumberOfParticles=%d\n",NumberOfNbodyParticles);
			if (NumberOfNbodyParticles != 0)
			{
				for (int dim = 0; dim < MAX_DIMENSION; dim++)
				{
					ierr = MPI_Recv(NbodyParticlePosition[dim], NumberOfNbodyParticles, MPI_DOUBLE, 1, 300, inter_comm, &status);
					ierr = MPI_Recv(NbodyParticleVelocity[dim], NumberOfNbodyParticles, MPI_DOUBLE, 1, 400, inter_comm, &status);
				}
#ifdef SEVN
				ierr = MPI_Recv(NbodyParticleInitialMass,		NumberOfNbodyParticles, MPI_DOUBLE, 1, 800,		inter_comm, &status);
				ierr = MPI_Recv(NbodyParticleWindEjectedMass,	NumberOfNbodyParticles, MPI_DOUBLE, 1, 900,		inter_comm, &status);
				ierr = MPI_Recv(NbodyParticleSNEjectedMass,		NumberOfNbodyParticles, MPI_DOUBLE, 1, 1000,	inter_comm, &status);
				ierr = MPI_Recv(NbodyParticleTemperature,		NumberOfNbodyParticles, MPI_DOUBLE, 1, 1100,	inter_comm, &status);

				ierr = MPI_Recv(NbodyParticleMass,				NumberOfNbodyParticles, MPI_DOUBLE, 1, 1600,	inter_comm, &status);
#endif
			}
			//fprintf(stderr,"NewNumberOfParticles=%d\n",NumberOfNewNbodyParticles);
			if (NumberOfNewNbodyParticles > 0) {
				for (int dim=0; dim<MAX_DIMENSION; dim++) {
					ierr = MPI_Recv(NewNbodyParticlePosition[dim], NumberOfNewNbodyParticles, MPI_DOUBLE, 1, 500, inter_comm, &status);
					ierr = MPI_Recv(NewNbodyParticleVelocity[dim], NumberOfNewNbodyParticles, MPI_DOUBLE, 1, 600, inter_comm, &status);
				}
#ifdef SEVN
				ierr = MPI_Recv(NewNbodyParticleInitialMass,		NumberOfNewNbodyParticles, MPI_DOUBLE, 1, 1200,	inter_comm, &status);
				ierr = MPI_Recv(NewNbodyParticleWindEjectedMass,	NumberOfNewNbodyParticles, MPI_DOUBLE, 1, 1300,	inter_comm, &status);
				ierr = MPI_Recv(NewNbodyParticleSNEjectedMass,		NumberOfNewNbodyParticles, MPI_DOUBLE, 1, 1400,	inter_comm, &status);
				ierr = MPI_Recv(NewNbodyParticleTemperature,		NumberOfNewNbodyParticles, MPI_DOUBLE, 1, 1500,	inter_comm, &status);

				ierr = MPI_Recv(NewNbodyParticleMass,				NumberOfNewNbodyParticles, MPI_DOUBLE, 1, 1700,	inter_comm, &status);
#endif
			}

			if ((NumberOfNbodyParticles+NumberOfNewNbodyParticles)!=0 && isNbodyParticleIdentification && isIdentificationOnTheFly) {
				ierr = MPI_Recv(NbodyClusterPosition, 3, MPI_DOUBLE, 1, 700, inter_comm, &status);
				//fprintf(stdout, "In Final, NbodyClusterPosition = (%e, %e, %e)\n", NbodyClusterPosition[0], NbodyClusterPosition[1], NbodyClusterPosition[2]);
				//fprintf(stderr, "In Final, NbodyClusterPosition = (%e, %e, %e)\n", NbodyClusterPosition[0], NbodyClusterPosition[1], NbodyClusterPosition[2]);
			}
			CommunicationInterBarrier();
			fprintf(stdout, "ENZO: Data received.\n");

			//fprintf(stderr,"NumberOfParticles after NBODY=%d\n",NumberOfNbodyParticles);
			//fprintf(stderr,"enzo: X=%e, V=%e\n ",NbodyParticlePosition[0][0], NbodyParticleVelocity[0][0]);


			/* Sending Index, NumberOfParticles, NbodyArrays to other processs */
			/*
				 MPI_Iscatterv(NbodyParticleID, LocalNumberAll, start_index_all, IntDataType,
				 NbodyParticleIDTemp, LocalNumberOfNbodyParticles, IntDataType, ROOT_PROCESSOR, enzo_comm,&request);
				 MPI_Wait(&request, &status);
				 */

			//CommunicationBarrier();
			if (NumberOfNbodyParticles != 0)
			{
				for (int dim = 0; dim < MAX_DIMENSION; dim++)
				{
					MPI_Iscatterv(NbodyParticlePosition[dim], LocalNumberAll, start_index_all, MPI_DOUBLE,
								  NbodyParticlePositionTemp[dim], LocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
								  &request);
					ierr = MPI_Wait(&request, &status);
					MPI_Iscatterv(NbodyParticleVelocity[dim], LocalNumberAll, start_index_all, MPI_DOUBLE,
								  NbodyParticleVelocityTemp[dim], LocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
								  &request);
					ierr = MPI_Wait(&request, &status);
				}
#ifdef SEVN
				MPI_Iscatterv(NbodyParticleInitialMass, LocalNumberAll, start_index_all, MPI_DOUBLE,
							  NbodyParticleInitialMassTemp, LocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
								&request);
				ierr = MPI_Wait(&request, &status);
				MPI_Iscatterv(NbodyParticleWindEjectedMass, LocalNumberAll, start_index_all, MPI_DOUBLE,
							  NbodyParticleWindEjectedMassTemp, LocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
							  &request);
				ierr = MPI_Wait(&request, &status);
				MPI_Iscatterv(NbodyParticleSNEjectedMass, LocalNumberAll, start_index_all, MPI_DOUBLE,
							  NbodyParticleSNEjectedMassTemp, LocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
							  &request);
				ierr = MPI_Wait(&request, &status);
				MPI_Iscatterv(NbodyParticleTemperature, LocalNumberAll, start_index_all, MPI_DOUBLE,
							  NbodyParticleTemperatureTemp, LocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
							  &request);
				ierr = MPI_Wait(&request, &status);

				MPI_Iscatterv(NbodyParticleMass, LocalNumberAll, start_index_all, MPI_DOUBLE,
							  NbodyParticleMassTemp, LocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
							  &request);
				ierr = MPI_Wait(&request, &status);
#endif
			}
			if (NumberOfNewNbodyParticles > 0) {
				for (int dim=0; dim<MAX_DIMENSION; dim++) {
					MPI_Iscatterv(NewNbodyParticlePosition[dim], NewLocalNumberAll, start_index_all_new, MPI_DOUBLE,
							NewNbodyParticlePositionTemp[dim], NewLocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
							&request);
					ierr = MPI_Wait(&request, &status);
					MPI_Iscatterv(NewNbodyParticleVelocity[dim], NewLocalNumberAll, start_index_all_new, MPI_DOUBLE,
							NewNbodyParticleVelocityTemp[dim], NewLocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
							&request);
					ierr  = MPI_Wait(&request, &status);
				}
#ifdef SEVN
				MPI_Iscatterv(NewNbodyParticleInitialMass, NewLocalNumberAll, start_index_all_new, MPI_DOUBLE,
							  NewNbodyParticleInitialMassTemp, NewLocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
							  &request);
				ierr = MPI_Wait(&request, &status);
				MPI_Iscatterv(NewNbodyParticleWindEjectedMass, NewLocalNumberAll, start_index_all_new, MPI_DOUBLE,
							  NewNbodyParticleWindEjectedMassTemp, NewLocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
							  &request);
				ierr = MPI_Wait(&request, &status);
				MPI_Iscatterv(NewNbodyParticleSNEjectedMass, NewLocalNumberAll, start_index_all_new, MPI_DOUBLE,
							  NewNbodyParticleSNEjectedMassTemp, NewLocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
							  &request);
				ierr = MPI_Wait(&request, &status);
				MPI_Iscatterv(NewNbodyParticleTemperature, NewLocalNumberAll, start_index_all_new, MPI_DOUBLE,
							  NewNbodyParticleTemperatureTemp, NewLocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
							  &request);
				ierr = MPI_Wait(&request, &status);

				MPI_Iscatterv(NewNbodyParticleMass, NewLocalNumberAll, start_index_all_new, MPI_DOUBLE,
							  NewNbodyParticleMassTemp, NewLocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
							  &request);
				ierr = MPI_Wait(&request, &status);
#endif
			}



			//fprintf(stderr,"Root:Done?2-5\n");

			fprintf(stderr, "FNC root... 1 start!\n");
			fflush(stderr);

			if (start_index_all != NULL)
				delete [] start_index_all;
			start_index_all = NULL;
			if (start_index_all != NULL)
				delete [] start_index_all_new;
			start_index_all_new = NULL;
			if (LocalNumberAll != NULL)
				delete [] LocalNumberAll;
			LocalNumberAll = NULL;
			if (NewLocalNumberAll != NULL)
				delete [] NewLocalNumberAll;
			NewLocalNumberAll = NULL;


			for (int dim=0; dim<MAX_DIMENSION; dim++) {
				if (NewNbodyParticlePosition[dim] != NULL)
					delete [] NewNbodyParticlePosition[dim];
				NewNbodyParticlePosition[dim] = NULL;

				if (NewNbodyParticleVelocity[dim] != NULL)
					delete [] NewNbodyParticleVelocity[dim];
				NewNbodyParticleVelocity[dim] = NULL;
			}

			fprintf(stderr, "FNC root... 1 done!\n");
			fflush(stderr);
#ifdef SEVN
			fprintf(stderr, "FNC root... 2 start!\n");
			fflush(stderr);
			if (NewNbodyParticleInitialMass != NULL)
				delete [] NewNbodyParticleInitialMass;
			NewNbodyParticleInitialMass = NULL;
			if (NewNbodyParticleWindEjectedMass != NULL)
				delete [] NewNbodyParticleWindEjectedMass;
			NewNbodyParticleWindEjectedMass = NULL;
			if (NewNbodyParticleSNEjectedMass != NULL)
				delete [] NewNbodyParticleSNEjectedMass;
			NewNbodyParticleSNEjectedMass = NULL;
			if (NewNbodyParticleTemperature != NULL)
				delete [] NewNbodyParticleTemperature;
			NewNbodyParticleTemperature = NULL;

			if (NewNbodyParticleMass != NULL)
				delete [] NewNbodyParticleMass;
			NewNbodyParticleMass = NULL;
			fprintf(stderr, "FNC root... 2 done!\n");
			fflush(stderr);
#endif


			DeleteNbodyArrays();

		}  // ENDIF: Root processor
		else {

			if (NumberOfNbodyParticles != 0)
			{
				MPI_Gather(&LocalNumberOfNbodyParticles, 1, IntDataType, NULL, NULL, IntDataType, ROOT_PROCESSOR, enzo_comm);
				MPI_Gather(&start_index, 1, IntDataType, NULL, NULL, IntDataType, ROOT_PROCESSOR, enzo_comm);
			}

			if (NumberOfNewNbodyParticles > 0)  {
				MPI_Gather(&NewLocalNumberOfNbodyParticles, 1, IntDataType, NULL, NULL, IntDataType, ROOT_PROCESSOR, enzo_comm);
				MPI_Gather(&start_index_new, 1, IntDataType, NULL, NULL, IntDataType, ROOT_PROCESSOR, enzo_comm);
			}

			MPI_Request request;
			MPI_Status status;
			int ierr;


			/* Receiving Index, NumberOfParticles, NbodyArrays from the root processs */
			if (NumberOfNbodyParticles != 0)
			{
				for (int dim = 0; dim < MAX_DIMENSION; dim++)
				{
					MPI_Iscatterv(NULL, NULL, NULL, MPI_DOUBLE,
								  NbodyParticlePositionTemp[dim], LocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
								  &request);
					ierr = MPI_Wait(&request, &status);
					MPI_Iscatterv(NULL, NULL, NULL, MPI_DOUBLE,
								  NbodyParticleVelocityTemp[dim], LocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
								  &request);
					ierr = MPI_Wait(&request, &status);
				}
#ifdef SEVN
				MPI_Iscatterv(NULL, NULL, NULL, MPI_DOUBLE,
							  NbodyParticleInitialMassTemp, LocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
							  &request);
				ierr = MPI_Wait(&request, &status);
				MPI_Iscatterv(NULL, NULL, NULL, MPI_DOUBLE,
							  NbodyParticleWindEjectedMassTemp, LocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
							  &request);
				ierr = MPI_Wait(&request, &status);
				MPI_Iscatterv(NULL, NULL, NULL, MPI_DOUBLE,
							  NbodyParticleSNEjectedMassTemp, LocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
							  &request);
				ierr = MPI_Wait(&request, &status);
				MPI_Iscatterv(NULL, NULL, NULL, MPI_DOUBLE,
							  NbodyParticleTemperatureTemp, LocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
							  &request);
				ierr = MPI_Wait(&request, &status);

				MPI_Iscatterv(NULL, NULL, NULL, MPI_DOUBLE,
							  NbodyParticleMassTemp, LocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
							  &request);
				ierr = MPI_Wait(&request, &status);
#endif
			}
			if (NumberOfNewNbodyParticles > 0) {
				for (int dim=0; dim<MAX_DIMENSION; dim++) {
					MPI_Iscatterv(NULL, NULL, NULL, MPI_DOUBLE,
							NewNbodyParticlePositionTemp[dim], NewLocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
							&request);
					ierr = MPI_Wait(&request, &status);
					MPI_Iscatterv(NULL, NULL, NULL, MPI_DOUBLE,
							NewNbodyParticleVelocityTemp[dim], NewLocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
							&request);
					ierr = MPI_Wait(&request, &status);
				}
#ifdef SEVN
				MPI_Iscatterv(NULL, NULL, NULL, MPI_DOUBLE,
							  NewNbodyParticleInitialMassTemp, NewLocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
							  &request);
				ierr = MPI_Wait(&request, &status);
				MPI_Iscatterv(NULL, NULL, NULL, MPI_DOUBLE,
							  NewNbodyParticleWindEjectedMassTemp, NewLocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
							  &request);
				ierr = MPI_Wait(&request, &status);
				MPI_Iscatterv(NULL, NULL, NULL, MPI_DOUBLE,
								NewNbodyParticleSNEjectedMassTemp, NewLocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
								&request);
				ierr = MPI_Wait(&request, &status);
				MPI_Iscatterv(NULL, NULL, NULL, MPI_DOUBLE,
								NewNbodyParticleTemperatureTemp, NewLocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
								&request);
				ierr = MPI_Wait(&request, &status);

				MPI_Iscatterv(NULL, NULL, NULL, MPI_DOUBLE,
							  NewNbodyParticleMassTemp, NewLocalNumberOfNbodyParticles, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm,
							  &request);
				ierr = MPI_Wait(&request, &status);
#endif
			}

		} // end else
#endif


		if (isNbodyParticleIdentification && isIdentificationOnTheFly)  {
			MPI_Bcast(NbodyClusterPosition, 3, MPI_DOUBLE, ROOT_PROCESSOR, enzo_comm);
			//fprintf(stderr, "In Final all, NbodyClusterPosition = (%e, %e, %e)\n", NbodyClusterPosition[0], NbodyClusterPosition[1], NbodyClusterPosition[2]);
		}

		/* Update Particle Velocity and Position Back to Grids */
		int count = 0;
		for (int level1=0; level1<MAX_DEPTH_OF_HIERARCHY-1;level1++)
			for (Temp = LevelArray[level1]; Temp; Temp = Temp->NextGridThisLevel)
#ifdef SEVN
				if (Temp->GridData->UpdateNbodyParticles(&count, 
							LocalNumberOfNbodyParticles, NbodyParticleIDTemp, 
							NbodyParticlePositionTemp, NbodyParticleVelocityTemp,
							NbodyParticleInitialMassTemp, NbodyParticleWindEjectedMassTemp,
							NbodyParticleSNEjectedMassTemp, NbodyParticleTemperatureTemp,
							NbodyParticleMassTemp,
							NewLocalNumberOfNbodyParticles, NewNbodyParticleIDTemp, 
							NewNbodyParticlePositionTemp, NewNbodyParticleVelocityTemp,
							NewNbodyParticleInitialMassTemp, NewNbodyParticleWindEjectedMassTemp,
							NewNbodyParticleSNEjectedMassTemp, NewNbodyParticleTemperatureTemp,
							NewNbodyParticleMassTemp
							) == FAIL) {
					ENZO_FAIL("Error in grid::CopyNbodyParticles.");
				}
#else
				if (Temp->GridData->UpdateNbodyParticles(&count, 
							LocalNumberOfNbodyParticles, NbodyParticleIDTemp, 
							NbodyParticlePositionTemp, NbodyParticleVelocityTemp,
							NewLocalNumberOfNbodyParticles, NewNbodyParticleIDTemp, 
							NewNbodyParticlePositionTemp, NewNbodyParticleVelocityTemp
							) == FAIL) {
					ENZO_FAIL("Error in grid::CopyNbodyParticles.");
				}
#endif

		/*
		fprintf(stderr,"Proc %d, # of Nbody = %d, count = %d\n", 
				MyProcessorNumber, LocalNumberOfNbodyParticles, count);
				*/

		fprintf(stderr, "FNC root... 3 start!\n");
		fflush(stderr);

		/* Destruct Arrays*/
		if (NbodyParticleIDTemp != NULL)
			delete [] NbodyParticleIDTemp;
		NbodyParticleIDTemp = NULL;

		if (NewNbodyParticleIDTemp != NULL)
			delete [] NewNbodyParticleIDTemp;
		NewNbodyParticleIDTemp = NULL;

		for (int dim=0; dim<MAX_DIMENSION; dim++) {
			if (NbodyParticlePositionTemp[dim] != NULL)
				delete [] NbodyParticlePositionTemp[dim];
			NbodyParticlePositionTemp[dim] = NULL;

			if (NbodyParticleVelocityTemp[dim] != NULL)
				delete [] NbodyParticleVelocityTemp[dim];
			NbodyParticleVelocityTemp[dim] = NULL;
		}
		fprintf(stderr, "FNC root... 3 done!\n");
		fflush(stderr);
#ifdef SEVN
		fprintf(stderr, "FNC root... 4 start!\n");
		fflush(stderr);
		if (NbodyParticleInitialMassTemp != NULL)
			delete [] NbodyParticleInitialMassTemp;
		NbodyParticleInitialMassTemp = NULL;
		if (NbodyParticleWindEjectedMassTemp != NULL)
			delete [] NbodyParticleWindEjectedMassTemp;
		NbodyParticleWindEjectedMassTemp = NULL;
		if (NbodyParticleSNEjectedMassTemp != NULL)
			delete [] NbodyParticleSNEjectedMassTemp;
		NbodyParticleSNEjectedMassTemp = NULL;
		if (NbodyParticleTemperatureTemp != NULL)
			delete [] NbodyParticleTemperatureTemp;
		NbodyParticleTemperatureTemp = NULL;

		if (NbodyParticleMassTemp != NULL)
			delete [] NbodyParticleMassTemp;
		NbodyParticleMassTemp = NULL;
		fprintf(stderr, "FNC root... 4 done!\n");
		fflush(stderr);
#endif

		fprintf(stderr, "FNC root... 5 start!\n");
		fflush(stderr);

		for (int dim=0; dim<MAX_DIMENSION; dim++) {
			if (NewNbodyParticlePositionTemp[dim] != NULL)
				delete [] NewNbodyParticlePositionTemp[dim];
			NewNbodyParticlePositionTemp[dim] = NULL;

			if (NewNbodyParticleVelocityTemp[dim] != NULL)
				delete [] NewNbodyParticleVelocityTemp[dim];
			NewNbodyParticleVelocityTemp[dim] = NULL;
		}
		fprintf(stderr, "FNC root... 5 done!\n");
		fflush(stderr);
#ifdef SEVN
		fprintf(stderr, "FNC root... 6 start!\n");
		fflush(stderr);
		if (NewNbodyParticleInitialMassTemp != NULL)
			delete [] NewNbodyParticleInitialMassTemp;
		NewNbodyParticleInitialMassTemp = NULL;
		if (NewNbodyParticleWindEjectedMassTemp != NULL)
			delete [] NewNbodyParticleWindEjectedMassTemp;
		NewNbodyParticleWindEjectedMassTemp = NULL;
		if (NewNbodyParticleSNEjectedMassTemp != NULL)
			delete [] NewNbodyParticleSNEjectedMassTemp;
		NewNbodyParticleSNEjectedMassTemp = NULL;
		if (NewNbodyParticleTemperatureTemp != NULL)
			delete [] NewNbodyParticleTemperatureTemp;
		NewNbodyParticleTemperatureTemp = NULL;

		if (NewNbodyParticleMassTemp != NULL)
			delete [] NewNbodyParticleMassTemp;
		NewNbodyParticleMassTemp = NULL;
		fprintf(stderr, "FNC root... 6 done!\n");
		fflush(stderr);
#endif
		NumberOfNewNbodyParticles = 0;

		} // ENDIF level
		return SUCCESS;
	}
#endif

