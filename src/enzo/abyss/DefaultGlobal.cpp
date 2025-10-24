#include <iostream>
#include <stdio.h>
#include "global.h"
#include <unordered_map>
#ifdef SEVN
#include <map>
#endif

int NumberOfWorker;

// Custom MPI data types
MPI_Datatype QueueType;
MPI_Datatype IparticleType;
MPI_Datatype JparticleType;

int *AvailableIndices; 
int NumberOfAvailableIndices;
std::unordered_map<int,int> PIDtoIndexMap; // (Query) EW: How about CM particles?
#ifndef INDIVIDUALSTAR
int *EnzoPIDs;
#else
int *displs = NULL;
#endif
int newNumberOfSingleParticle;
int NumberOfParticle; // The number of active particles (single + CM ptcl)
int NumberOfSingleParticle; // The number of single particles (only single, not CM ptcl)
int NewCMPID;
int LastParticleIndex; // The last index of particle array; for few-body case by EW 2025.3.10

// Enzo to Nbody
Particle* FirstEnzoParticle;
double EnzoLength, EnzoMass, EnzoVelocity, EnzoTime, EnzoForce, EnzoAcceleration;
double EnzoCurrentTime;
double AbyssCenter[3];

// Few-Body
std::unordered_map<int, int> CMPtclWorker;	   // by EW 2025.1.4 // unordered_map by EW 2025.1.11
std::unordered_map<int, int> PrevCMPtclWorker; // by EW 2025.1.4 // unordered_map by EW 2025.1.11

// Time
double global_time;
double global_time_irr;
ULL NextRegTimeBlock;
double outputTimeStep;
double endTime;

// i/o
char* fname;
double inputTime;
bool restart;
char* foutput;
bool IsOutput;
double outputTime;
int outNum;

FILE* binout;
FILE* nbpout;
FILE* mergerout;
#ifdef SEVN
FILE* SEVNout;
IO* sevnio = nullptr;
std::multimap<double, int> SEVNList; // This constains the time of next SEVN evolution time and the particle index by EW 2025.3.27

int NumberOfEnzoSEVNParticle;		// This is the number of SEVN particles in Enzo, not in Abyss by EW 2025.4.27
int newNumberOfEnzoSEVNParticle;	// This is the number of SEVN particles newly added in Enzo, not in Abyss by EW 2025.4.27
std::unordered_map<int,int> PIDtoIndexMap_SEVN;
int *EnzoPIDs_SEVN;

std::vector<StarSEVN*> SEVNList_Enzo;
std::vector<double> creation_time_Enzo;
std::vector<double> world_time_Enzo;
#endif
FILE* workerout;

#ifdef PerformanceTrace
Performance performance;
#endif



double InitialNeighborRadius2;
double EPS2;
int FixNumNeighbor;

void DefaultGlobal() {

	/* Timesteps */
	endTime = 1;
	outputTimeStep = outputTimeStep/endTime; // endTime should be Myr

	global_variable->time_block = -30;
	global_variable->block_max = static_cast<ULL>(pow(2, -global_variable->time_block));
	global_variable->time_step = std::pow(2,global_variable->time_block);

	inputTime = 0.0;
	endTime = 0.0;
	outputTimeStep = 0.;

	global_time = 0.;
	outputTime = 0.;

	if (AbyssProcessorNumber == 0)  {
#ifndef INDIVIDUALSTAR
		EnzoPIDs         = new int[MaxNumParticle];
#else
		displs = new int[NumberOfAbyssProcessors+1];
		for (int i=0; i<NumberOfAbyssProcessors+1; i++) {
			displs[i] = 0;
		}
#endif
		AvailableIndices = new int[MaxNumParticle];
		NumberOfAvailableIndices = 0;
		for (int i=0; i<MaxNumParticle; i++) {
			AvailableIndices[i] = -1;
		}
#ifdef SEVN
		EnzoPIDs_SEVN    = new int[MaxNumParticle];
#endif
	}

	NumberOfWorker = NumberOfAbyssProcessors-1;

	NewCMPID = -1;

#ifdef SEVN
	std::vector<std::string> args = {"empty", // Not used
		// "-myself", "/data/vinicius/NbodyPlus/SEVN",
		"-tables", "/home/vinicius/install/sevn/tables/SEVNtracks_parsec_ov04_AGB", 
		//  "-tables", "/data/vinicius/NbodyPlus/SEVN/tables/SEVNtracks_MIST_AGB",
		// "-tables_HE", "/data/vinicius/NbodyPlus/SEVN/tables/SEVNtracks_parsec_pureHe36",
		// "-turn_WR_to_pureHe", "false",
		"-xspinmode", "geneva",
		"-hardmode", "disabled",
		"-tmode", "disabled",
		// "-snmode", "delayed",
		// "-Z", "0.0002",
		// "-spin", "0.0",
		// "-tini", "zams",
		// "-tf", "end",
		// "-dtout", "events",
		};
	std::vector<char*> c_args;
	for (auto& arg : args) {
		c_args.push_back(&arg[0]);
	}

	sevnio = new IO; // (SEVN Query) We should initialize sevnio only once. Here in Abyss, and somewhere else in Enzo by EW 2025.3.27
	sevnio->load(c_args.size(), c_args.data());
#endif

}
