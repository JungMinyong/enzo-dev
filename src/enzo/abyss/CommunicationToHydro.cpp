#include <mpi.h>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <cmath>
#include "global.h"
#include "def.h"
#define FLOAT double
#define MAX_NUMBER_OF_OUTPUT_REDSHIFTS    500
#include "../CosmologyParameters.h"
#include "../phys_constants.h"


#define No_COM_EVOLUTION



extern int StarParticleFeedback;
extern double StarMassEjectionFraction;
extern int ComovingCoordinates;
double EnzoCurrentTime, ClusterRadius2;
double ClusterAcceleration[Dim], ClusterPosition[Dim], ClusterVelocity[Dim], EnzoClusterPosition[Dim+1];
double eta_tmp;
int FixNumNeighbor0, IdentifyNbodyParticles;
int BinaryRegularization, IdentifyOnTheFly;

//double KSTime;
//double KSDistance;


void deleteParticle(int &PID, int &index);
//int writeParticle(std::vector<Particle*> &particle, double MinRegTime, int outputNum);
//void InitializeParticle(std::vector<Particle*> &particle);
//void InitializeNewParticle(std::vector<Particle*> &particle, int offset, int newSize);
void GetCenterOfMass(double *mass, double *x[Dim], double *v[Dim], double x_com[], double v_com[], int N);
void GetNewCenterOfMass(int *PID, double *mass2, double *x2[Dim], double *v2[Dim], int n2, double x_X[], double v_X[]);
//void UpdateNextRegTime(std::vector<Particle*> &particle);
int CommunicationInterBarrier();
void broadcastFromRoot(int &data);
void broadcastFromRoot(double &data);
void broadcastFromRoot(ULL &data);
//void CalculateAllAccelerationOnGPU(std::vector<Particle*> &particle);
//void KSTermination(Particle* ptclCM, std::vector<Particle*> &particle, double current_time, ULL current_block);
void InitializationAfterCommunication();


using namespace std;
const int width = 18;

//int InitialCommunication(std::vector<Particle*> &particle) {
int InitialCommunication() {

	int *PID;
	double *Mass, *Position[Dim], *Velocity[Dim], *BackgroundAcceleration[Dim];
	double *CreationTime, *DynamicalTime;
	double TimeStep, TimeUnits, LengthUnits, VelocityUnits, DensityUnits;//MassUnits;

	MPI_Request request;
	MPI_Status status;

	fprintf(nbpout, "NBODY+: First Waiting for Enzo to receive data...\n");
	fflush(nbpout);

	CommunicationInterBarrier();
	fprintf(nbpout, "NBODY+: Frist Receiving data from Enzo...\n");
	MPI_Recv(&NumberOfSingleParticle, 1, MPI_INT, 0, 100, inter_comm, &status);
	fprintf(nbpout, "NBODY+: NumberOfSingleParticle=%d\n", NumberOfSingleParticle);
	if (NumberOfSingleParticle != 0) {
		PID           = new int[NumberOfSingleParticle];
		Mass          = new double[NumberOfSingleParticle];
		CreationTime  = new double[NumberOfSingleParticle];
		DynamicalTime = new double[NumberOfSingleParticle];
		MPI_Recv(PID          , NumberOfSingleParticle, MPI_INT   , 0, 200, inter_comm, &status);
		MPI_Recv(Mass         , NumberOfSingleParticle, MPI_DOUBLE, 0, 201, inter_comm, &status);
		MPI_Recv(CreationTime , NumberOfSingleParticle, MPI_DOUBLE, 0, 202, inter_comm, &status);
		MPI_Recv(DynamicalTime, NumberOfSingleParticle, MPI_DOUBLE, 0, 203, inter_comm, &status);

		for (int dim=0; dim<Dim; dim++) {
			Position[dim] = new double[NumberOfSingleParticle];
			Velocity[dim] = new double[NumberOfSingleParticle];
			MPI_Recv(Position[dim], NumberOfSingleParticle, MPI_DOUBLE, 0, 300, inter_comm, &status);
			MPI_Recv(Velocity[dim], NumberOfSingleParticle, MPI_DOUBLE, 0, 400, inter_comm, &status);
		}

		for (int dim=0; dim<Dim; dim++) {
			BackgroundAcceleration[dim] = new double[NumberOfSingleParticle];
			MPI_Recv(BackgroundAcceleration[dim], NumberOfSingleParticle, MPI_DOUBLE, 0, 500, inter_comm, &status);
		}
	}
	MPI_Recv(&TimeStep,                 1, MPI_DOUBLE, 0,  600, inter_comm, &status);
	MPI_Recv(&TimeUnits,                1, MPI_DOUBLE, 0,  700, inter_comm, &status);
	MPI_Recv(&LengthUnits,              1, MPI_DOUBLE, 0,  800, inter_comm, &status);
	MPI_Recv(&DensityUnits,             1, MPI_DOUBLE, 0,  900, inter_comm, &status);
	MPI_Recv(&VelocityUnits,            1, MPI_DOUBLE, 0, 1000, inter_comm, &status);
	MPI_Recv(&StarParticleFeedback    , 1, MPI_INT   , 0, 1100, inter_comm, &status);
	MPI_Recv(&StarMassEjectionFraction, 1, MPI_DOUBLE, 0, 1200, inter_comm, &status);
	MPI_Recv(&EnzoCurrentTime         , 1, MPI_DOUBLE, 0, 1300, inter_comm, &status);
	MPI_Recv(&EPS2                    , 1, MPI_DOUBLE, 0, 1400, inter_comm, &status);
	MPI_Recv(&eta_tmp                 , 1, MPI_DOUBLE, 0, 1500, inter_comm, &status);
	MPI_Recv(&InitialNeighborRadius   , 1, MPI_DOUBLE, 0, 1600, inter_comm, &status);
	MPI_Recv(EnzoClusterPosition      , 4, MPI_DOUBLE, 0, 1700, inter_comm, &status);
	MPI_Recv(&IdentifyNbodyParticles  , 1, MPI_INT   , 0, 1750, inter_comm, &status);
	MPI_Recv(&IdentifyOnTheFly        , 1, MPI_INT   , 0, 1775, inter_comm, &status);
	MPI_Recv(&FixNumNeighbor          , 1, MPI_INT   , 0, 1800, inter_comm, &status);
	MPI_Recv(&MaxNumNeighbor          , 1, MPI_INT   , 0, 1850, inter_comm, &status);
	MPI_Recv(&BinaryRegularization    , 1, MPI_INT   , 0, 1900, inter_comm, &status);
	//MPI_Recv(&KSDistance              , 1, MPI_DOUBLE, 0, 2000, inter_comm, &status);
	//MPI_Recv(&KSTime                  , 1, MPI_DOUBLE, 0, 2100, inter_comm, &status);
	//MPI_Recv(&HydroMethod         , 1, MPI_INT   , 0, 1200, inter_comm, &status);
	fprintf(nbpout, "data receiving!\n");
	fflush(nbpout);

	MPI_Recv(&ComovingCoordinates        , 1, MPI_INT   , 0, 3000, inter_comm, &status);
	fprintf(nbpout, "ComovingCoordinates=%d\n",ComovingCoordinates);
	fflush(nbpout);
	if (ComovingCoordinates) {
		fprintf(nbpout, "Cosmo data receiving!\n");
		MPI_Recv(&HubbleConstantNow        , 1, MPI_DOUBLE, 0, 3010, inter_comm, &status);
		MPI_Recv(&OmegaMatterNow           , 1, MPI_DOUBLE, 0, 3020, inter_comm, &status);
		MPI_Recv(&OmegaDarkMatterNow       , 1, MPI_DOUBLE, 0, 3030, inter_comm, &status);
		MPI_Recv(&OmegaLambdaNow           , 1, MPI_DOUBLE, 0, 3040, inter_comm, &status);
		MPI_Recv(&OmegaRadiationNow        , 1, MPI_DOUBLE, 0, 3050, inter_comm, &status);
		MPI_Recv(&ComovingBoxSize          , 1, MPI_DOUBLE, 0, 3060, inter_comm, &status);
		MPI_Recv(&MaxExpansionRate         , 1, MPI_DOUBLE, 0, 3070, inter_comm, &status);
		MPI_Recv(&InitialTimeInCodeUnits   , 1, MPI_DOUBLE, 0, 3080, inter_comm, &status);
		MPI_Recv(&InitialRedshift          , 1, MPI_DOUBLE, 0, 3090, inter_comm, &status);
		MPI_Recv(&FinalRedshift            , 1, MPI_DOUBLE, 0, 3100, inter_comm, &status);
		MPI_Recv(&CosmologyTableNumberOfBins, 1, MPI_INT  , 0, 3110, inter_comm, &status);
		MPI_Recv(&CosmologyTableLogtIndex   , 1, MPI_INT  , 0, 3120, inter_comm, &status);
		MPI_Recv(&CosmologyTableLogaInitial , 1, MPI_DOUBLE, 0, 3130, inter_comm, &status);
		MPI_Recv(&CosmologyTableLogaFinal   , 1, MPI_DOUBLE, 0, 3140, inter_comm, &status);
		CosmologyTableLoga = new double[CosmologyTableNumberOfBins];
		CosmologyTableLogt = new double[CosmologyTableNumberOfBins];
		MPI_Recv(CosmologyTableLoga, CosmologyTableNumberOfBins, MPI_DOUBLE, 0, 3150, inter_comm, &status);
		MPI_Recv(CosmologyTableLogt, CosmologyTableNumberOfBins, MPI_DOUBLE, 0, 3160, inter_comm, &status);
	}
	CommunicationInterBarrier();
	fprintf(nbpout, "Data received!\n");


	ClusterRadius2 = EnzoClusterPosition[3];

	// Enzo to Nbody unit convertors
	//EnzoMass         = MassUnits/Msun/mass_unit;
	EnzoMass         = DensityUnits*pow(LengthUnits,3.)/Msun/mass_unit;
	EnzoLength       = LengthUnits/pc/position_unit;
	EnzoVelocity     = VelocityUnits/pc*yr/velocity_unit;
	EnzoTime         = TimeUnits/yr/time_unit;
	//EnzoAcceleration = LengthUnits/TimeUnits/TimeUnits/pc*yr*yr/position_unit*time_unit*time_unit;
	EnzoAcceleration = EnzoLength/EnzoTime/EnzoTime;

	// Unit conversion
	EnzoTimeStep       = TimeStep*EnzoTime;
	EnzoCurrentTime   *= EnzoTime;
	if (EPS2 < 0)
		EPS2 = -1;
	else {
		EPS2 *= EnzoLength;
		EPS2 *= EPS2;
	}

	InitialNeighborRadius *= EnzoLength;
	FixNumNeighbor0    = FixNumNeighbor;

	fprintf(nbpout, "Enzo Time                = %lf\n", TimeStep);
	fprintf(nbpout, "Nbody Time               = %lf\n", EnzoTimeStep);
	fprintf(nbpout, "EPS2                     = %lf pc**2\n", EPS2*position_unit*position_unit);
	fprintf(nbpout, "InitialNeighborRadius        = %.2e pc\n", InitialNeighborRadius*position_unit);
	fprintf(nbpout, "eta                      = %lf\n", eta);
	fprintf(nbpout, "ClusterRadius2           = %.2e pc**2\n", ClusterRadius2*EnzoLength*EnzoLength*position_unit*position_unit);
	fprintf(nbpout, "StarMassEjectionFraction = %lf\n", StarMassEjectionFraction);
	fprintf(nbpout, "StarParticleFeedback     = %d\n", StarParticleFeedback);
	fprintf(nbpout, "FixNumNeighbor           = %d\n", FixNumNeighbor);
	fprintf(nbpout, "BinaryRegularization     = %d\n", BinaryRegularization);
	//fprintf(nbpout, "KSTime                   = %lf\n", KSTime);
	//fprintf(nbpout, "KSDistance               = %lf\n", KSDistance);
	fprintf(nbpout, "IdentifyNbodyParticles   = %d\n", IdentifyNbodyParticles);
	fprintf(nbpout, "IdentifyOnTheFly         = %d\n\n", IdentifyOnTheFly);


	for (int dim=0; dim<Dim; dim++) {
		ClusterAcceleration[dim] = 0;
		ClusterPosition[dim]     = 0;
		ClusterVelocity[dim]     = 0;
	}
	if (NumberOfSingleParticle != 0) {
		// set COM for background acc
		GetCenterOfMass(Mass, Position, Velocity, ClusterPosition, ClusterVelocity, NumberOfSingleParticle);

#ifdef COM_EVOLUTION
		ClusterAcceleration[0] = 0.;
		ClusterAcceleration[1] = 0.;
		ClusterAcceleration[2] = 0.;
		double total_mass=0.;
		for (int i=0; i<NumberOfSingleParticle; i++) {
			for (int dim=0; dim<Dim; dim++) {
				ClusterAcceleration[dim] += Mass[i]*BackgroundAcceleration[dim][i];
			}
			total_mass += Mass[i];
		}

		for (int dim=0; dim<Dim; dim++)
			ClusterAcceleration[dim] /= total_mass;
#endif

		for (int i=0; i<NumberOfSingleParticle; i++) {
			for (int dim=0; dim<Dim; dim++) {
#ifdef COM_EVOLUTION
				BackgroundAcceleration[dim][i] -= ClusterAcceleration[dim];
				Velocity[dim][i]               -= ClusterVelocity[dim];
#endif
				Position[dim][i]               -= ClusterPosition[dim];
			}
			EnzoPIDs[i] = PID[i];
			particles[i].set(PID, Mass, CreationTime, DynamicalTime, Position, Velocity,
							 BackgroundAcceleration, i);
			PIDtoIndexMap.insert({PID[i], i});
			particles[i].ParticleIndex = i;
		}

		
		std::cerr << "enzo Pos     :" << Position[0][0] << ", " << Position[0][1] << std::endl;
		std::cerr << "nbody Pos    :" << Position[0][0]*EnzoLength << ", " << Position[0][1]*EnzoLength << std::endl;
		std::cerr << "enzo Vel     :" << Velocity[1][0] << ", " << Velocity[1][1] << std::endl;
		std::cerr << "nbody Vel    :" << Velocity[1][0]*EnzoVelocity << ", " << Velocity[1][1]*EnzoVelocity << std::endl;
		std::cerr << "km/s Vel     :" << Velocity[1][0]*VelocityUnits/1e5 << ", " << Velocity[1][1]*VelocityUnits/1e5 << std::endl;
		std::cerr << "enzo  Mass   :" << Mass[0] << std::endl;
		std::cerr << "nbody Mass   :" << Mass[0]*EnzoMass << std::endl;
		std::cerr << "CreationTime :" << CreationTime[0] << std::endl;
		std::cerr << "DynamicalTime:" << DynamicalTime[0] << std::endl;

		/*
	for (Particle* ptcl:particle) {
		fprintf(nbpout, "NBODY0: PID=%d\n", ptcl->PID);
		fprintf(nbpout, "NBODY0: Mass Of NewNbodyParticles=%.3e\n", ptcl->Mass*mass_unit);
		fprintf(nbpout, "NBODY0: Vel  Of news=(%.3e, %.3e, %.3e)\n", 
				ptcl->Velocity[0]*velocity_unit/yr*pc/1e5, ptcl->Velocity[1]*velocity_unit/yr*pc/1e5, ptcl->Velocity[2]*velocity_unit/yr*pc/1e5);
		fprintf(nbpout, "NBODY0: Pos  Of news=(%.3e, %.3e, %.3e)\n",
				ptcl->Position[0]*position_unit, ptcl->Position[1]*position_unit, ptcl->Position[2]*position_unit);
		fprintf(nbpout, "NBODY0: Acc  Of regs=(%.3e, %.3e, %.3e)\n",
				ptcl->a_reg[0][0], ptcl->a_reg[1][0], ptcl->a_reg[2][0]);
		fprintf(nbpout, "NBODY0: Acc  Of irrs=(%.3e, %.3e, %.3e)\n",
				ptcl->a_irr[0][0], ptcl->a_irr[1][0], ptcl->a_irr[2][0]);
		fprintf(nbpout, "NBODY0: Back Acc  Of news=(%.3e, %.3e, %.3e)\n",
				ptcl->BackgroundAcceleration[0], ptcl->BackgroundAcceleration[1], ptcl->BackgroundAcceleration[2]);
	}
	*/

		delete [] PID;
		delete [] Mass;
		delete [] CreationTime;
		delete [] DynamicalTime;
		for (int dim=0; dim<Dim; dim++) {
			delete[] Position[dim];
			delete[] Velocity[dim];
			delete[] BackgroundAcceleration[dim];
		}
	}

	NumberOfParticle = NumberOfSingleParticle;
	LastParticleIndex = NumberOfSingleParticle-1;


	FixNumNeighbor = std::min((int) std::floor(NumberOfSingleParticle/2), FixNumNeighbor0);



	fprintf(nbpout, "NBODY+: %d particles loaded!\n", NumberOfSingleParticle);
	fflush(stdout);
	fflush(stderr);
	fflush(nbpout);
	//fflush(gpuout);
	//fflush(binout);
	return true;
}


int ReceiveFromEnzo() {

	int i;
	int *PID, *newPID;
	double *BackgroundAcceleration[Dim];
	double *Mass, *newMass, *newPosition[Dim], *newVelocity[Dim], *newBackgroundAcceleration[Dim];
	double *newCreationTime, *newDynamicalTime;
	double TimeStep;

	MPI_Request request;
	MPI_Status status;
	fprintf(nbpout, "NBODY+: Waiting for Enzo to receive data...\n");


	fprintf(stdout, "NBODY+: Waiting for Enzo to receive data...\n");
	CommunicationInterBarrier();
	fprintf(nbpout, "NBODY+: Receiving data from Enzo...\n");

	// Existing Particle Information
	
	int EnzoNumberOfSingleParticle;
	MPI_Recv(&EnzoNumberOfSingleParticle, 1, MPI_INT, 0, 10, inter_comm, &status);
	if (EnzoNumberOfSingleParticle != NumberOfSingleParticle) {
		fprintf(stderr, "The numbers of nbody and hydro do not match: NumberOfSingleParticle=%d, EnzoNumberOfSingleParticle=%d\n",
			 	NumberOfSingleParticle, EnzoNumberOfSingleParticle);
		fflush(stderr);
		throw runtime_error("");
	}
	if (NumberOfSingleParticle != 0) {
		PID = new int[NumberOfSingleParticle];
		Mass = new double[NumberOfSingleParticle];
		MPI_Recv(PID , NumberOfSingleParticle, MPI_INT   , 0, 25, inter_comm, &status);
		MPI_Recv(Mass, NumberOfSingleParticle, MPI_DOUBLE, 0, 40, inter_comm, &status);
		for (int dim=0; dim<Dim; dim++) {
			BackgroundAcceleration[dim] = new double[NumberOfSingleParticle];
			MPI_Recv(BackgroundAcceleration[dim], NumberOfSingleParticle, MPI_DOUBLE, 0, 50, inter_comm, &status);
		}
	}


	// New Particle Information
	MPI_Recv(&newNumberOfSingleParticle, 1, MPI_INT, 0, 100, inter_comm, &status);
	if (newNumberOfSingleParticle > 0)
	{
		newPID           = new int[newNumberOfSingleParticle];
		newMass          = new double[newNumberOfSingleParticle];
		newCreationTime  = new double[newNumberOfSingleParticle];
		newDynamicalTime = new double[newNumberOfSingleParticle];
		MPI_Recv(newPID          , newNumberOfSingleParticle, MPI_INT   , 0, 200, inter_comm, &status);
		MPI_Recv(newMass         , newNumberOfSingleParticle, MPI_DOUBLE, 0, 201, inter_comm, &status);
		MPI_Recv(newCreationTime , newNumberOfSingleParticle, MPI_DOUBLE, 0, 202, inter_comm, &status);
		MPI_Recv(newDynamicalTime, newNumberOfSingleParticle, MPI_DOUBLE, 0, 203, inter_comm, &status);

		for (int dim=0; dim<Dim; dim++) {
			newPosition[dim]               = new double[newNumberOfSingleParticle];
			newVelocity[dim]               = new double[newNumberOfSingleParticle];
			newBackgroundAcceleration[dim] = new double[newNumberOfSingleParticle];
			MPI_Recv(newPosition[dim]              , newNumberOfSingleParticle, MPI_DOUBLE, 0, 300, inter_comm, &status);
			MPI_Recv(newVelocity[dim]              , newNumberOfSingleParticle, MPI_DOUBLE, 0, 400, inter_comm, &status);
			MPI_Recv(newBackgroundAcceleration[dim], newNumberOfSingleParticle, MPI_DOUBLE, 0, 500, inter_comm, &status);
		}
	}

	/***************************************************
	 * at some point we have to reconstruct *particle* vector because there is redundance due to particle removal.
	***************************************************/

	// Timestep
	MPI_Recv(&TimeStep       , 1, MPI_DOUBLE, 0, 600, inter_comm, &status);
	double OldEnzoCurrentTime = EnzoCurrentTime;
	MPI_Recv(&EnzoCurrentTime, 1, MPI_DOUBLE, 0, 700, inter_comm, &status);
	CommunicationInterBarrier();

	std::cout << "Enzo  Time    :" << EnzoCurrentTime << std::endl;
	//std::cout << "Nbody Time    :" << OldEnzoCurrentTime+particle[0]->CurrentTimeReg*EnzoTimeStep << std::endl;
	EnzoTimeStep = TimeStep*EnzoTime;

	std::cout << "NBODY+: Data trnsferred!" << std::endl;
	fprintf(stdout, "NBODY+: Data trnsferred!\n");
	EnzoCurrentTime = EnzoCurrentTime*EnzoTime;
	std::cout << "Enzo Time    :" << EnzoCurrentTime*1e4 << " Myr" << std::endl;
	std::cout << "Enzo TimeStep:" << TimeStep  << std::endl;
	std::cout << "EnzoTimeStep :" << EnzoTimeStep*1e4 << " Myr" << std::endl;
	//std::cout << "Nbody Mass    :" << particle[0]->Mass << std::endl;
	//std::cout << "Enzo  Mass    :" << Mass[0]*EnzoMass << std::endl;
	//std::cout << "enzo Time :" << TimeStep << std::endl;
	//std::cout << "nbody Time:" << EnzoTimeStep << std::endl;
	FixNumNeighbor = std::min((int) std::floor((NumberOfSingleParticle+newNumberOfSingleParticle-1)/2.0), FixNumNeighbor0);



	// COM conversion
	// later on we might need to take mass weight into account 
	// 1. F=ma; 2. F -> F_com; 3. F_com -> a_com
	ClusterAcceleration[0] = 0.0;
	ClusterAcceleration[1] = 0.0;
	ClusterAcceleration[2] = 0.0;

#ifdef COM_EVOLUTION
	double total_mass = 0.;
	if (NumberOfSingleParticle != 0) {
		for (int i=0; i<NumberOfSingleParticle; i++) {
			for (int dim=0; dim<Dim; dim++) {
				ClusterAcceleration[dim] += Mass[i]*BackgroundAcceleration[dim][i];
			}
			total_mass += Mass[i];
		}
		std::cout << 
			std::setw(width)  << "NBODY : ClusterAcceleration = " << 
			std::setw(width)  << ClusterAcceleration[0] << 
			std::setw(width)  << ClusterAcceleration[1] << 
			std::setw(width)  << ClusterAcceleration[2] << 
			std::endl;

		std::cout << 
			std::setw(width)  << "NBODY : ClusterPosition = " << 
			std::setw(width)  << ClusterPosition[0] << 
			std::setw(width)  << ClusterPosition[1] << 
			std::setw(width)  << ClusterPosition[2] << 
			std::endl;	 	
		std::cout << 
			std::setw(width)  << "NBODY : ClusterVelocity = " << 
			std::setw(width)  << ClusterVelocity[0] << 
			std::setw(width)  << ClusterVelocity[1] << 
			std::setw(width)  << ClusterVelocity[2] << 
			std::endl;
	}

	if (newNumberOfSingleParticle != 0) {
		// we need to make adjustment to COM
		GetNewCenterOfMass(PID, newMass, newPosition, newVelocity, newNumberOfSingleParticle, 
				ClusterPosition, ClusterVelocity);

		for (int i=0; i<newNumberOfSingleParticle; i++) {
			for (int dim=0; dim<Dim; dim++) {
				ClusterAcceleration[dim] += newMass[i]*newBackgroundAcceleration[dim][i];
			}
			total_mass += newMass[i];
		}
	}
	if (total_mass != 0.) {
		for (int dim=0; dim<Dim; dim++) {
			ClusterAcceleration[dim] /= total_mass;
		}
	}
#endif



	// Update Existing Particles
	// need to update if the ids match between Enzo and Nbody
	if (NumberOfSingleParticle != 0) {
		std::cerr << "In ReceiveFromEnzo, NumberOfSingleParticle = "<< NumberOfSingleParticle<< std::endl;
		// loop for PID, going backwards to update the NextParticle
		for (int i=NumberOfSingleParticle-1; i>=0; i--) {
#ifdef COM_EVOLUTION
			for (int dim=0; dim<Dim; dim++) {
				BackgroundAcceleration[dim][i] -= ClusterAcceleration[dim];
			}
#endif
			EnzoPIDs[i] = PID[i];
			particles[PIDtoIndexMap[PID[i]]].update(Mass, BackgroundAcceleration, i);
		} //endfor i
	} //endif nnb




#define no_star_formation_location_test
#ifdef star_formation_location_test
	std::vector<Particle*> new_particle; 
#endif
	//Particle *newPtcl = new Particle[newNumberOfSingleParticle]; // we shouldn't delete it, maybe make it to vector
	
	if (newNumberOfSingleParticle > 0) {
		std::cerr << "In ReceiveFromEnzo, newNumberOfSingleParticle = "<< newNumberOfSingleParticle << std::endl;
		/*
		for (int i = 0; i < newNumberOfSingleParticle; i++) {
			fprintf(stderr, "NBODY: Mass Of NewNbodyParticles=%.3e\n", newMass[i]);
			fprintf(stderr, "NBODY: Vel  Of news=(%.3e, %.3e, %.3e)\n", 
					newVelocity[0][i], newVelocity[1][i], newVelocity[2][i]);
			fprintf(stderr, "NBODY: Pos  Of news=(%.3e, %.3e, %.3e)\n",
					newPosition[0][i], newPosition[1][i], newPosition[2][i]);
			fprintf(stderr, "NBODY: Acc  Of news=(%.3e, %.3e, %.3e)\n",
					newBackgroundAcceleration[0][i], newBackgroundAcceleration[1][i], newBackgroundAcceleration[2][i]);
		}
		*/

		// Update New Particles
		int index = -1;
		for (int i=0; i<newNumberOfSingleParticle; i++) {
			for (int dim=0; dim<Dim; dim++) {
#ifdef COM_EVOLUTION
				newBackgroundAcceleration[dim][i] -= ClusterAcceleration[dim];
				newVelocity[dim][i]               -= ClusterVelocity[dim];
#endif
				newPosition[dim][i]               -= ClusterPosition[dim];
			}

			if (NumberOfAvailableIndices != 0) {
				index = AvailableIndices[NumberOfAvailableIndices-1];
				AvailableIndices[NumberOfAvailableIndices-1] = -1;
                NumberOfAvailableIndices--;
			}
			else {
				index = LastParticleIndex + 1;
				LastParticleIndex++;
			}

			particles[index].set(newPID, newMass, newCreationTime, newDynamicalTime, newPosition, newVelocity,
					//newBackgroundAcceleration, NormalStar+SingleParticle+NewParticle, i);
					newBackgroundAcceleration, NormalStar, i);
			particles[index].ParticleIndex = index;
			PIDtoIndexMap.insert({newPID[i], index});
			EnzoPIDs[NumberOfSingleParticle+i] = newPID[i];

#ifdef star_formation_location_test
			new_particle.push_back(particles[NumberOfSingleParticle+i]);
#endif
		}

#ifdef star_formation_location_test
		writeParticle(new_particle, EnzoCurrentTime, outNum++);
#endif

		/*
		std::cout << "enzo Pos  :" << newPosition[0][0] << ", " << newPosition[0][1] << std::endl;
		std::cout << "nbody Pos :" << newPosition[0][0]*EnzoLength << ", " << newPosition[0][1]*EnzoLength << std::endl;
		std::cout << "enzo Vel  :" << newVelocity[1][0] << ", " << newVelocity[1][1] << std::endl;
		std::cout << "nbody Vel :" << newVelocity[1][0]*EnzoVelocity << ", " << newVelocity[1][1]*EnzoVelocity << std::endl;
		std::cout << "km/s Vel  :" << newVelocity[1][0]*velocity_unit/1e5/yr*pc << ", " << newVelocity[1][1]*velocity_unit/1e5/yr*pc << std::endl;
		std::cout << "enzo  Mass:" << newMass[0] << std::endl;
		std::cout << "nbody Mass:" << newMass[0]*EnzoMass << std::endl;
		*/
		std::cout << "NBODY+    : "  << newNumberOfSingleParticle << " new particles loaded!" << std::endl;

		// This includes modification of regular force and irregular force
	} // endif newNumberOfSingleParticle != 0

	/*
	// update particle order 
	if (particle.size() > 1) {
		i=0;
		for (Particle* ptcl:particle) {
			ptcl->ParticleOrder = i;
			i++;
		}
	}
	*/

	/*
	if (newNumberOfSingleParticle > 0) {
		for (Particle* ptcl:particle) {
			fprintf(nbpout, "before init, %3d neighbors of %3d : ", ptcl->NumberOfAC, ptcl->PID);
			fprintf(stderr, "before init, %3d neighbors of %3d : ", ptcl->NumberOfAC, ptcl->PID);
			for (Particle* nn:ptcl->ACList) {
				fprintf(nbpout, "%3d, ", nn->PID);
				fprintf(stderr, "%3d, ", nn->PID);
			}
			fprintf(nbpout, "\n");
			fprintf(stderr, "\n");
		}
	}
	*/

	/* Initialize New Particles */
	/*
	if (newNumberOfSingleParticle > 0) {
		if (NumberOfSingleParticle < 2)
			InitializeParticle(particle);
		else
			InitializeNewParticle(particle, NumberOfSingleParticle, newNumberOfSingleParticle);

		fprintf(stderr, "NBODY:(New) PID=\n");
		for (int i=NumberOfSingleParticle; i<NumberOfSingleParticle+newNumberOfSingleParticle; i++) {
			fprintf(stderr, "%d, ", particles[i].PID);
		}
		fprintf(stderr, "\n ");
	*/
	/*
	for (Particle* ptcl:particle) {
		fprintf(nbpout, "after init, %3d neighbors of %3d : ", ptcl->NumberOfAC, ptcl->PID);
		fprintf(stderr, "after init, %3d neighbors of %3d : ", ptcl->NumberOfAC, ptcl->PID);
		for (Particle* nn:ptcl->ACList) {
			fprintf(nbpout, "%3d, ", nn->PID);
			fprintf(stderr, "%3d, ", nn->PID);
			if (nn->PID == ptcl->PID) {
				fflush(nbpout);
				fflush(stderr);
				throw std::runtime_error("Fatal error in the neighbor list: CommunicationToHydro.\n");
			}
		}
		fprintf(nbpout, "\n");
		fprintf(stderr, "\n");
	}
	Particle* ptcl;
	fprintf(stderr, "New particles PID = ");
	for (int i=0; i<newNumberOfSingleParticle; i++) {
		ptcl = particle[NumberOfSingleParticle+i];
		fprintf(stderr, "%d ", ptcl->PID);
		fprintf(nbpout, "NBODY2: Mass Of %d =%.3e Msun\n", ptcl->PID, ptcl->Mass*mass_unit);
		fprintf(nbpout, "NBODY2: Vel  Of news=(%.3e, %.3e, %.3e) km/s\n",
				ptcl->Velocity[0]*velocity_unit/yr*pc/1e5, ptcl->Velocity[1]*velocity_unit/yr*pc/1e5, ptcl->Velocity[2]*velocity_unit/yr*pc/1e5);
		fprintf(nbpout, "NBODY2: Pos  Of news=(%.3e, %.3e, %.3e) pc\n",
				ptcl->Position[0]*position_unit, ptcl->Position[1]*position_unit, ptcl->Position[2]*position_unit);
		fprintf(nbpout, "NBODY2: Acc  Of regs=(%.3e, %.3e, %.3e)\n",
				ptcl->a_reg[0][0], ptcl->a_reg[1][0], ptcl->a_reg[2][0]);
		fprintf(nbpout, "NBODY2: Acc  Of irrs=(%.3e, %.3e, %.3e)\n",
				ptcl->a_irr[0][0], ptcl->a_irr[1][0], ptcl->a_irr[2][0]);
		fprintf(nbpout, "NBODY2: Back Acc  Of news=(%.3e, %.3e, %.3e)\n",
				ptcl->BackgroundAcceleration[0], ptcl->BackgroundAcceleration[1], ptcl->BackgroundAcceleration[2]);
		fprintf(nbpout, "NBODY2: Timestep (%.3e, %.3e), nn=%d\n",
				ptcl->TimeStepReg, ptcl->TimeStepIrr, ptcl->NumberOfAC);
		fprintf(nbpout, "NBODY2: Timestep (%llu, %llu), (%d, %d)\n",
				ptcl->TimeBlockReg, ptcl->TimeBlockIrr, ptcl->TimeLevelReg, ptcl->TimeLevelIrr);
	}
	fprintf(stderr, "\n");
	*/
	//}

	//RegularList.clear();

	NumberOfSingleParticle += newNumberOfSingleParticle;
	NumberOfParticle 	   += newNumberOfSingleParticle;


	if (NumberOfSingleParticle != 0) {
		delete[] PID;
		delete[] Mass;
		for (int dim=0; dim<Dim; dim++) {
			delete[] BackgroundAcceleration[dim];
		}
	}
	if (newNumberOfSingleParticle != 0) {
		delete[] newPID;
		delete[] newMass;
		delete[] newDynamicalTime;
		delete[] newCreationTime;
		for (int dim=0; dim<Dim; dim++) {
			delete[] newBackgroundAcceleration[dim];
			delete[] newPosition[dim];
			delete[] newVelocity[dim];
		}
	}


	//  (Query) Do I need this?
	fprintf(nbpout, "NBODY+    : Acceleration for particles on GPU.\n");
	/*
	if (NumberOfSingleParticle > 1) {
		CalculateAllAccelerationOnGPU(particle);
	}*/

	//UpdateNextRegTime(particle);
	fprintf(nbpout, "NBODY+    : Acceleration and neighbors are updated.\n");

	fprintf(nbpout, "NBODY+    : In ReceiveFromEzno (after new particle might be added): \n");
	fprintf(nbpout, "NBODY+    : original NumberOfSingleParticle      = %d (+%d)\n", NumberOfSingleParticle-newNumberOfSingleParticle, newNumberOfSingleParticle);
	fprintf(nbpout, "NBODY+    : newly updated NumberOfSingleParticle = %d\n", NumberOfSingleParticle, newNumberOfSingleParticle);
	//fprintf(nbpout, "NBODY+    : Particle size     = %d\n", particle.size());
	fprintf(nbpout, "NBODY+    : NextRegTimeStep   = %.3e\n", NextRegTimeBlock*time_step);
	fprintf(nbpout, "NBODY+    : NextRegTimeBlock  = %d\n", NextRegTimeBlock);
	//fprintf(nbpout, "NBODY+    : RegularList size  = %d\n", RegularList.size());
	fprintf(nbpout, "NBODY+    : FixNumNeighbor    = %d\n", FixNumNeighbor);


	fprintf(stderr, "NBODY+    : In ReceiveFromEzno (after new particle might be added): \n");
	fprintf(stderr, "NBODY+    : original NumberOfSingleParticle      = %d (+%d)\n", NumberOfSingleParticle-newNumberOfSingleParticle, newNumberOfSingleParticle);
	fprintf(stderr, "NBODY+    : newly updated NumberOfSingleParticle = %d\n", NumberOfSingleParticle, newNumberOfSingleParticle);
	//fprintf(stderr, "NBODY+    : Particle size     = %d\n", particle.size());
	//fprintf(stderr, "NBODY+    : RegularList size = %d\n", RegularList.size());
	fflush(stderr);
	fflush(nbpout);
	//fflush(gpuout);
	//fflush(binout);


	return true;
}



int SendToEnzo() {

	std::cout << "NBODY+: Entering SendToEnzo..." << std::endl;
	if (NumberOfSingleParticle == 0 && newNumberOfSingleParticle == 0) {
		std::cout << "NBODY+: Skipping SendToEnzo..." << std::endl;
		fflush(stdout);
		fflush(stderr);
		return SUCCESS;
	}
	MPI_Request request;
	MPI_Status status;

	double *Position[Dim], *Velocity[Dim], *newPosition[Dim], *newVelocity[Dim];
	int index;
	Particle *ptcl;

	for (int dim=0; dim<Dim; dim++) {
		if (NumberOfSingleParticle-newNumberOfSingleParticle != 0) {
			Position[dim]    = new double[NumberOfSingleParticle-newNumberOfSingleParticle];
			Velocity[dim]    = new double[NumberOfSingleParticle-newNumberOfSingleParticle];
		}

		if (newNumberOfSingleParticle > 0) {
			newPosition[dim] = new double[newNumberOfSingleParticle];
			newVelocity[dim] = new double[newNumberOfSingleParticle];
		}
	}

	/*
	std::cout << "NBODY+: NumberOfSingleParticle=" << NumberOfSingleParticle << ", newNumberOfSingleParticle=" << newNumberOfSingleParticle << std::endl;
	std::cout << "NBODY+: size=" << particle.size() << std::endl;
	std::cout << "NBODY+: FirstParticleInEnzo PID=" << \
		FirstParticleInEnzo->PID << " in SendToEzno" << std::endl;
		*/


	//for (Particle* ptcl:particle) {


	// (Query to EW) in this part, all the sdar objects should return to individual particles.
#ifdef FEWBODY
	if (BinaryList.size() > 0) {

		int offset = particle.size() - BinaryList.size();

		for (Binary* bin:BinaryList) {
			if (bin->ptclCM->Position[0] != bin->ptclCM->Position[0])
				fprintf(stderr, "COM=%d, %e\n", bin->ptclCM->PID, bin->ptclCM->Position[0]);
			bin->ptclCM->convertBinaryCoordinatesToCartesian();
			bin->ptclCM->isErase = true;
			particle.push_back(bin->ptclCM->BinaryParticleI);
			particle.push_back(bin->ptclCM->BinaryParticleJ);
			if (bin->ptclCM->BinaryParticleI->Position[0] != bin->ptclCM->BinaryParticleI->Position[0]) {
				fprintf(stderr, "%d, %e\n", bin->ptclCM->BinaryParticleI->PID, bin->ptclCM->BinaryParticleI->Position[0]);
			}
			std::cerr << bin->ptclCM->BinaryParticleI->PID << std::endl;
			std::cerr << bin->ptclCM->BinaryParticleJ->PID << std::endl;
			delete bin;
		}
		fflush(stderr);
		particle.erase(
				std::remove_if(particle.begin(), particle.end(),
					[](Particle* p) {
					bool to_remove = p->isErase;
					if (to_remove) delete p;
					return to_remove;
					}),
				particle.end());
		std::cerr << "CM Particle dissociation failed!\n" << std::endl;
		if (offset+BinaryList.size()*2 != particle.size()) {
			std::cerr << "CM Particle dissociation failed!\n" << std::endl;
			fprintf(stderr, "offset = %d, BinaryList=%d, ptclsize=%d\n", offset, BinaryList.size(), particle.size());
			fflush(stderr);
			throw runtime_error("");
		}

		/*
		int i=0;
		for (Particle* ptcl:particle) {
			ptcl->ParticleOrder = i;
			i++;
		}
		*/
		InitializeNewParticle(particle, offset, BinaryList.size()*2);
	}
	BinaryList.clear();
	std::cerr << "Binary Done!\n" << std::endl;
#endif // FewBody

	int NumberOfEscapeParticle = 0;
	double r2;
	double TimeStep=EnzoTimeStep/EnzoTime;



#ifdef COM_EVOLUTION
	// COM evolution	
	for (int dim=0; dim<Dim; dim++) {	
		ClusterPosition[dim] += ClusterVelocity[dim]*TimeStep;
		ClusterPosition[dim] += ClusterAcceleration[dim]*TimeStep*TimeStep/2;
		ClusterVelocity[dim] += ClusterAcceleration[dim]*TimeStep;
	}
#else
	double NbodyCOM[Dim] = {0,0,0};
	double mass = 0;
	for (int i = 0; i < NumberOfSingleParticle; i++)
	{
		index = PIDtoIndexMap[EnzoPIDs[i]];
		ptcl = &particles[index];

		for (int dim = 0; dim < Dim; dim++)
		{
			NbodyCOM[dim] += ptcl->Mass * ptcl->Position[dim];
		}
		mass += ptcl->Mass;
	}

	for (int dim = 0; dim < Dim; dim++)
	{
		NbodyCOM[dim] /= mass;
		NbodyCOM[dim] /= EnzoLength;
		NbodyCOM[dim] += ClusterPosition[dim];
	}
#endif

	if (NumberOfSingleParticle - newNumberOfSingleParticle > 0) {
		//fprintf(stderr, "Sending PID order= ");
		for (int i=0; i<NumberOfSingleParticle-newNumberOfSingleParticle; i++) {
			//fprintf(stderr, "%d, ",ptcl->PID);
			index = PIDtoIndexMap[EnzoPIDs[i]];
			ptcl = &particles[index];
			r2 = 0;
			for (int dim=0; dim<Dim; dim++) {
				Position[dim][i]  = ptcl->Position[dim]/EnzoLength;
				Velocity[dim][i]  = ptcl->Velocity[dim]/EnzoVelocity;

#ifdef COM_EVOLUTION
				if (IdentifyNbodyParticles && IdentifyOnTheFly)
					r2 += Position[dim][i]*Position[dim][i];
				Velocity[dim][i] += ClusterVelocity[dim];

				// COM correction
				Position[dim][i] += ClusterPosition[dim];
#else
				// COM correction
				Position[dim][i] += ClusterPosition[dim];
				if (IdentifyNbodyParticles && IdentifyOnTheFly)
					r2 += (Position[dim][i]-NbodyCOM[dim])*(Position[dim][i]-NbodyCOM[dim]);
#endif

				if (IdentifyNbodyParticles && !IdentifyOnTheFly)
					r2 += (Position[dim][i]-EnzoClusterPosition[dim])*(Position[dim][i]-EnzoClusterPosition[dim]);
			}

			if (IdentifyNbodyParticles && ClusterRadius2 > 0 && r2 > ClusterRadius2) { // in Enzo Unit
				Position[0][i] -= 20;
				deleteParticle(EnzoPIDs[i],index);
				NumberOfEscapeParticle++;
			}
			//fprintf(stdout, "NBODY+: pid= %d, x=%e\n",ptcl->PID,Position[0][i]);

			if ((ptcl == nullptr) && (i != NumberOfSingleParticle-newNumberOfSingleParticle-1))
			{
				std::cout << "NBODY+: Warning! ParticleChain for Communication has been broken!" << std::endl;
				std::cerr << "NBODY+: Warning! ParticleChain for Communication has been broken!" << std::endl;
				fprintf(stderr, "%d-th particle, NumberOfSingleParticle=%d, newNumberOfSingleParticle=%d\n",i, NumberOfSingleParticle, newNumberOfSingleParticle);
				fprintf(nbpout, "%d-th particle, NumberOfSingleParticle=%d, newNumberOfSingleParticle=%d\n",i, NumberOfSingleParticle, newNumberOfSingleParticle);
				fflush(stderr);
				fflush(stdout);
				fflush(nbpout);
				throw std::runtime_error("CommunicationToHydro.cpp:757");
			}
		}
		//fprintf(stderr, "\n");
	}

	if (newNumberOfSingleParticle > 0) {
		int offset = NumberOfSingleParticle-newNumberOfSingleParticle;
		for (int i=0; i<newNumberOfSingleParticle; i++) {
			index = PIDtoIndexMap[EnzoPIDs[i+offset]];
			ptcl = &particles[index];
			r2 = 0;
			for (int dim=0; dim<Dim; dim++) {
				newPosition[dim][i]  = ptcl->Position[dim]/EnzoLength;
				newVelocity[dim][i]  = ptcl->Velocity[dim]/EnzoVelocity;
				newVelocity[dim][i] += ClusterVelocity[dim];

#ifdef COM_EVOLUTION
				if (IdentifyNbodyParticles && IdentifyOnTheFly)
					r2 += newPosition[dim][i]*newPosition[dim][i];

				// COM correction
				newPosition[dim][i] += ClusterPosition[dim];
#else
				// COM correction
				newPosition[dim][i] += ClusterPosition[dim];
				if (IdentifyNbodyParticles && IdentifyOnTheFly)
					r2 += (newPosition[dim][i]-NbodyCOM[dim])*(newPosition[dim][i]-NbodyCOM[dim]);
#endif

				if (IdentifyNbodyParticles && !IdentifyOnTheFly)
					r2 += (newPosition[dim][i]-EnzoClusterPosition[dim])*(newPosition[dim][i]-EnzoClusterPosition[dim]);
			}
			if (IdentifyNbodyParticles && ClusterRadius2 > 0 && r2 > ClusterRadius2) {
				newPosition[0][i] -= 20;
				deleteParticle(EnzoPIDs[i+offset],index);
				NumberOfEscapeParticle++;
				// (Query) binary termination?
			}
		}
	}

	//std::cerr << "NBODY+: Waiting for Enzo to sent data..." << std::endl;
	fprintf(nbpout, "NBODY+: Waiting for Enzo to sent data...\n");
	CommunicationInterBarrier();
	if (NumberOfSingleParticle-newNumberOfSingleParticle != 0)
	{
		for (int dim = 0; dim < Dim; dim++)
		{
			MPI_Send(Position[dim], NumberOfSingleParticle - newNumberOfSingleParticle, MPI_DOUBLE, 0, 300, inter_comm);
			MPI_Send(Velocity[dim], NumberOfSingleParticle - newNumberOfSingleParticle, MPI_DOUBLE, 0, 400, inter_comm);
		}
	}
	//std::cerr << "NBODY+: Escape particles=" << EscapeParticleNum << std::endl;

	//fprintf(stderr,"NewNumberOfSingleParticles=%d\n",NumberOfNewNbodyParticles);
	if (newNumberOfSingleParticle > 0) {
		for (int dim=0; dim<Dim; dim++) {
			MPI_Send(newPosition[dim], newNumberOfSingleParticle, MPI_DOUBLE, 0, 500, inter_comm);
			MPI_Send(newVelocity[dim], newNumberOfSingleParticle, MPI_DOUBLE, 0, 600, inter_comm);
		}
	}


	if (NumberOfSingleParticle != 0 && IdentifyOnTheFly) {
		//fprintf(stderr, "NBODY: ClusterPosition =(%lf, %lf, %lf)\n",
				//ClusterPosition[0]-0.5, ClusterPosition[1]-0.5, ClusterPosition[2]-0.5);
		/*
		fprintf(stderr, "NBODY: ClusterPosition =(%lf, %lf, %lf)\n",
				(ClusterPosition[0]-0.5)*EnzoLength*position_unit,
				(ClusterPosition[1]-0.5)*EnzoLength*position_unit,
				(ClusterPosition[2]-0.5)*EnzoLength*position_unit);*/
#ifdef COM_EVOLUTION
		MPI_Send(ClusterPosition, 3, MPI_DOUBLE, 0, 700, inter_comm);
#else
		MPI_Send(NbodyCOM, 3, MPI_DOUBLE, 0, 700, inter_comm);
#endif
	}

	CommunicationInterBarrier();
	fprintf(stdout, "NBODY+: Data sent!\n");




	/*
	fprintf(stderr, "NBODY:(particle before) PID=");
	for (Particle* ptcl:particle) {
		fprintf(stderr, "%d, ", ptcl->PID);
	}
	fprintf(stderr, "\n");
	*/



	// erase from other particles' neighbor
	// and correct force, but not too much worry about a_reg for now

	/*
	int size_tmp = 0;
	for (Particle* ptcl: particle) {
		size_tmp = ptcl->NumberOfAC;
		auto it = ptcl->ACList.erase(
				std::remove_if(ptcl->ACList.begin(), ptcl->ACList.end(),
					[ptcl](Particle* p) {
					bool to_remove = p->isErase;
					if (to_remove) {
					//fprintf(stderr, "earsing %d of %d\n", p->PID, ptcl->PID);
						double a[Dim], adot[Dim];
						ptcl->ComputeAcceleration(p,a,adot);
						for (int dim=0; dim<Dim; dim++) {
							ptcl->a_irr[dim][0] -= a[dim];
							ptcl->a_irr[dim][1] -= adot[dim];
						}
					}
					return to_remove;
					}),
				ptcl->ACList.end());

		if (size_tmp != ptcl->ACList.size())  {
			//fprintf(stderr, "Acceleration Correction PID=%d\n", ptcl->PID);
			for (int dim=0;dim<Dim;dim++) {
				ptcl->a_tot[dim][0] = ptcl->a_reg[dim][0] + ptcl->a_irr[dim][0];
				ptcl->a_tot[dim][1] = ptcl->a_reg[dim][1] + ptcl->a_irr[dim][1];
			}
			ptcl->NumberOfAC = ptcl->ACList.size();
			//fprintf(stderr, "NN=%d, a_irr=%e\n", ptcl->NumberOfAC, ptcl->a_irr[0][0]);
		}
	}
	*/


	/*
	fprintf(stderr, "NBODY:(particle after) PID=\n");
	for (Particle* ptcl:particle) {
		fprintf(stderr, "%d, ", ptcl->PID);
		}
	fprintf(stderr, "\n ");
	*/



	for (int dim = 0; dim < Dim; dim++)
	{
		if (NumberOfSingleParticle-newNumberOfSingleParticle != 0) {
			delete[] Position[dim];
			delete[] Velocity[dim];
		}
		if (newNumberOfSingleParticle > 0)
		{
			delete[] newPosition[dim];
			delete[] newVelocity[dim];
		}
	}

	NumberOfSingleParticle -= NumberOfEscapeParticle;
	NumberOfParticle -= NumberOfEscapeParticle;
	std::cout << "NBODY+: Sending data finished." << std::endl;


	fprintf(stderr, "NBODY+    : In SendToEnzo (after particle might be escaped): \n");
	fprintf(stderr, "NBODY+    : original NumberOfSingleParticle      = %d (-%d)\n", NumberOfSingleParticle+NumberOfEscapeParticle, NumberOfEscapeParticle);
	fprintf(stderr, "NBODY+    : newly updated NumberOfSingleParticle = %d\n", NumberOfSingleParticle);

	fprintf(nbpout, "NBODY+    : In SendToEnzo (after particle might be escaped): \n");
	fprintf(nbpout, "NBODY+    : original NumberOfSingleParticle      = %d (-%d)\n", NumberOfSingleParticle+NumberOfEscapeParticle, NumberOfEscapeParticle);
	fprintf(nbpout, "NBODY+    : newly updated NumberOfSingleParticle = %d\n", NumberOfSingleParticle);

	fflush(stderr);
	fflush(nbpout);
	//fflush(gpuout);
	fflush(binout);
	return true;
}

/* Adjust COM according to new particles
X = (a1+b1+c1)/M
Y = (a2+b2)/N
Z = (a1+b1+c1+a2+b2)/(M+N) = (X*M+Y*N)/(M+N) -> New COM
Z-X = N*(Y-X)/(M+N) -> This should be applied to particles
*/
// this should be improved by using iterative loop
void GetNewCenterOfMass(int* PID, double *mass2, double *x2[Dim], double *v2[Dim], int n2, double x_X[], double v_X[]) {

	double M=0., N=0.;
	double x_Y[Dim], v_Y[Dim], x_Z[Dim], v_Z[Dim];

	for (int dim=0; dim<Dim; dim++) {
		x_Y[dim] = 0;
		x_Z[dim] = 0;
		v_Y[dim] = 0;
		v_Z[dim] = 0;
	}

	for (int i=0; i<NumberOfSingleParticle; i++) {
		M += particles[PIDtoIndexMap[PID[i]]].Mass/EnzoMass;
	}


	for (int i=0; i<n2; i++) {
		for (int dim=0; dim<Dim; dim++) {
			x_Y[dim] += mass2[i]*x2[dim][i];
			v_Y[dim] += mass2[i]*v2[dim][i];
		}
		N += mass2[i];
	}
	for (int dim=0; dim<Dim; dim++) {
		x_Y[dim] = x_Y[dim]/N;
		v_Y[dim] = v_Y[dim]/N;
		x_Z[dim] = (x_X[dim]*M+x_Y[dim]*N)/(M+N);
		v_Z[dim] = (v_X[dim]*M+v_Y[dim]*N)/(M+N);
	}


	// Adjustment to particles
	for (int i=0; i<NumberOfSingleParticle; i++) {
		for (int dim=0; dim<Dim;dim++) {
			particles[PIDtoIndexMap[PID[i]]].Position[dim] - (x_Z[dim] - x_X[dim]) * EnzoLength;
			particles[PIDtoIndexMap[PID[i]]].Velocity[dim] - (v_Z[dim] - v_X[dim]) * EnzoVelocity;
		}
	}

	for (int dim=0; dim<Dim; dim++) {
		x_X[dim] = x_Z[dim];
		v_X[dim] = v_Z[dim];
	}
}

void deleteParticle(int &PID, int &index) {
	particles[index].isActive = false;
	PIDtoIndexMap.erase(PID);
	AvailableIndices[NumberOfAvailableIndices] = index;
	NumberOfAvailableIndices++;
}

