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
#include <stdio.h>
#ifdef NBODY
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
		ptcl->GetBackgroundAcceleration();
		for (int dim=0; dim<MAX_DIMENSION; dim++) {
			this->Position[dim]              = ptcl->ReturnPosition()[dim];
			this->Velocity[dim]              = ptcl->ReturnVelocity()[dim];
			this->BackgroundAcceleration[dim] = ptcl->ReturnBackgroundAcceleration()[dim];
		}
		this->ID 						= ptcl->ReturnID();
		this->Mass          = ptcl->ReturnMass();
		this->CreationTime  = ptcl->ReturnBirthTime();
		this->DynamicalTime = ptcl->ReturnLifetime();
		this->Metallicity   = ptcl->ReturnMetallicity();
	};
#endif
#ifdef ABYSS_ONLY
#ifdef COM_EVOLUTION
	void reposition(const double *pos, const double *vel, const double *acc) {
      for (int dim = 0; dim < MAX_MAX_DIMENSION; dim++) {
				Velocity[dim]	-= vel[dim];
        Position[dim] -= pos[dim];
        BackgroundAcceleration[dim] -= acc[dim];
			}
	};
#else
	void reposition(const double *pos) {
      for (int dim = 0; dim < MAX_DIMENSION; dim++) {
        Position[dim] -= pos[dim];
			}
	};
#endif
#endif
#define DEBUG
#ifdef DEBUG
	void print(const double &umass, const double &upos, const double &uvel){
          fprintf(stderr,
                  "PID: %d. Mass: %e Msun\n"
									"x: %e pc, y: %e pc, z: %e\nvx: %e "
                  "km/s, vy: %e km/s, vz: %e km/s\n",
                  ID, Mass*umass, 
									Position[0]*upos, Position[1]*upos,Position[2]*upos,
									Velocity[0]*uvel, Velocity[1]*uvel,Velocity[2]*uvel);
	};
#endif
};

struct ParticleSendDataType{
	int ID;
	double BackgroundAcceleration[MAX_DIMENSION];
	//double Mass;

#ifdef ENZO_ONLY
	void copyFrom(Star *ptcl) {
		ptcl->GetBackgroundAcceleration();
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
	double Mass;
	double WindEjectedMass;
	double SNEjectedMass;
	double Temperature;
#endif


#ifdef ENZO_ONLY
	void copyTo(Star *ptcl) {
		ptcl->SetPosition(this->Position);
		ptcl->SetVelocity(this->Velocity);
		ptcl->UpdateToGridParticle(this->Position, this->Velocity);
		ptcl->DeleteBackgroundAcceleration();
#ifdef SEVN
		//ptcl->ID = this->identifier;
		//this->Mass          = ptcl->Mass;
		//this->CreationTime  = ptcl->BirthTime;
		//this->DynamicalTime = ptcl->LifeTime;
		//this->Metallicity   = ptcl->Metallicity;
#endif
	};
#endif
};
#endif
#endif
#endif