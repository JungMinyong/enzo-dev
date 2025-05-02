/***********************************************************************
	/
	/  FIND ALL STAR PARTICLES OVER ALL PROCESSORS
	/
	/  written by: Yongseok Jo
	/  date:       November, 2022
	/
 ************************************************************************/

#ifndef __NBODYROUTINES_H
#define __NBODYROUTINES_H
#ifdef ENZO_ONLY
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

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
#include "phys_constants.h"
#include "Star.h"
#else
#define MAX_DIMENSION 3
#endif

#ifndef INDIVIDUALSTAR
#ifdef ENZO_ONLY
#include "macros_and_parameters.h"
#include "global_data.h"
#include "typedefs.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"


int FindTotalNumberOfNbodyParticles(LevelHierarchyEntry *LevelArray[],int *LocalNumberOfNbodyParticles, bool prepareNbodyComputation);
int FindTotalNumberOfNbodyParticles(LevelHierarchyEntry *LevelArray[],int *LocalNumberOfNbodyParticles, int *NewLocalNumberOfNbodyParticles, bool prepareNbodyComputation);
int FindStartIndex(int* LocalNumberOfNbodyParticles);
void InitializeNbodyArrays(void);
void CopyNbodyArrayToOld(void);
void MatchAccelerationWithIndex(void);
void DeleteNbodyArrays(void);
#ifdef USE_MPI
void Scan(int *in, int *inout, int *len, MPI_Datatype *dptr);
#endif
void DeleteNbodyArrays(void);
#endif
#else


struct ParticleDataType{
	int ID;
	double Position[MAX_DIMENSION];
	double Velocity[MAX_DIMENSION];
	double BackgroundAcceleration[MAX_DIMENSION];
	double Mass;
	double CreationTime;
	double DynamicalTime;
	double Metallicity;

#ifdef ENZO_ONLY
	void copyFrom(Star *ptcl) {
		for (int dim=0; dim<MAX_DIMENSION; dim++) {
			this->Position[dim]              = ptcl->ReturnPosition()[dim];
			this->Velocity[dim]              = ptcl->ReturnVelocity()[dim];
			this->BackgroundAcceleration[dim] = ptcl->ReturnBackgroundAcceleration()[dim];
		}
		this->ID = ptcl->ReturnID();
		this->Mass          = ptcl->ReturnMass();
		this->CreationTime  = ptcl->ReturnBirthTime();
		this->DynamicalTime = ptcl->ReturnLifetime();
		this->Metallicity   = ptcl->ReturnMetallicity();
	};
#endif
};

struct ParticleSendDataType{
	int ID;
	double BackgroundAcceleration[MAX_DIMENSION];
	//double Mass;

#ifdef ENZO_ONLY
	void copyFrom(Star *ptcl) {
		for (int dim=0; dim<MAX_DIMENSION; dim++) {
			this->BackgroundAcceleration[dim] = ptcl->ReturnBackgroundAcceleration()[dim];
		}
		this->ID = ptcl->ReturnID();
		//this->Mass          = ptcl->Mass;
	};
#endif
};

struct ParticleReceiveDataType{
	int	ID;
	double Position[MAX_DIMENSION];
	double Velocity[MAX_DIMENSION];
	//double Mass;
	//double CreationTime;
	//double DynamicalTime;
	//double Metallicity;

#ifdef SEVN
	double InitialMass;
	double WindEjectedMass;
	double SNEjectedMass;
	double Temperature;
#endif

#ifdef ENZO_ONLY
	void copyTo(Star *ptcl) {
		// I might generate a map to boost this process (ID matching needed), if so the maps should go into Star.
		for (int dim=0; dim<MAX_DIMENSION; dim++) {
			//ptcl->pos[dim]              = this->Position[dim];
			//ptcl->vel[dim]              = this->Velocity[dim];
		}
		//ptcl->ID = this->identifier;
		//this->Mass          = ptcl->Mass;
		//this->CreationTime  = ptcl->BirthTime;
		//this->DynamicalTime = ptcl->LifeTime;
		//this->Metallicity   = ptcl->Metallicity;
	};
#endif
};
#endif
#endif