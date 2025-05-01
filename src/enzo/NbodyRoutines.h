/***********************************************************************
	/
	/  FIND ALL STAR PARTICLES OVER ALL PROCESSORS
	/
	/  written by: Yongseok Jo
	/  date:       November, 2022
	/
 ************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
#include "macros_and_parameters.h"
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
void Scan(int *in, int *inout, int *len, MPI_Datatype *dptr);
void DeleteNbodyArrays(void);



struct ParticleDataType{
	PINT ID;
	double Position[MAX_DIMENSION];
	double Velocity[MAX_DIMENSION];
	double BackgrounAcceleration[MAX_DIMENSION];
	double Mass;
	double CreationTime;
	double DynamicalTime;
	double Metallicity;

	void copyFrom(Star *ptcl) {
		for (int dim=0; dim<MAX_DIMENSION; dim++) {
			this->Position[dim]              = ptcl->ReturnPosition()[dim];
			this->Velocity[dim]              = ptcl->ReturnVelocity()[dim];
			this->BackgrounAcceleration[dim] = ptcl->ReturnBackgroundAcceleration()[dim];
		}
		this->ID = ptcl->ReturnID();
		this->Mass          = ptcl->ReturnMass();
		this->CreationTime  = ptcl->ReturnBirthTime();
		this->DynamicalTime = ptcl->ReturnLifeTime();
		this->Metallicity   = ptcl->ReturnMetallicity();
	};
};

struct ParticleSendDataType{
	PINT ID;
	double BackgrounAcceleration[MAX_DIMENSION];
	//double Mass;

	void copyFrom(Star *ptcl) {
		for (int dim=0; dim<MAX_DIMENSION; dim++) {
			this->BackgrounAcceleration[dim] = ptcl->ReturnBackgroundAcceleration()[dim];
		}
		this->ID = ptcl->ReturnID();
		//this->Mass          = ptcl->Mass;
	};
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

/*
	void copyTo(Star *ptcl) {
		// I might generate a map to boost this process (ID matching needed), if so the maps should go into Star.
		for (int dim=0; dim<MAX_DIMENSION; dim++) {
			ptcl->Position[dim]              = this->pos[dim];
			ptcl->Velocity[dim]              = this->vel[dim];
		}
		ptcl->ID = this->identifier;
		//this->Mass          = ptcl->Mass;
		//this->CreationTime  = ptcl->BirthTime;
		//this->DynamicalTime = ptcl->LifeTime;
		//this->Metallicity   = ptcl->Metallicity;
	};*/
};