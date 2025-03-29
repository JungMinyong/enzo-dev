#include <iostream>
#include <stdio.h>
#include "global.h"
#include <unordered_map>
#ifdef SEVN
#include <map>
#endif

std::unordered_set<int> RegularList;
int NumberOfWorker;
int NumberOfCommunication;
int *AvailableIndices; 
int NumberOfAvailableIndices;
std::unordered_map<int,int> PIDtoIndexMap; // (Query) EW: How about CM particles?
int *EnzoPIDs; 
int newNumberOfSingleParticle;
int NumberOfParticle; // The number of active particles (single + CM ptcl)
int NumberOfSingleParticle; // The number of single particles (only single, not CM ptcl)
int NewPID;
int LastParticleIndex; // The last index of particle array; for few-body case by EW 2025.3.10

// Few-Body
std::unordered_map<int, int> CMPtclWorker;	   // by EW 2025.1.4 // unordered_map by EW 2025.1.11
std::unordered_map<int, int> PrevCMPtclWorker; // by EW 2025.1.4 // unordered_map by EW 2025.1.11

// Task
int Task[NumberOfTask];

// Time
double global_time;
double global_time_irr;
ULL NextRegTimeBlock;
double outputTimeStep;
double endTime;

double binary_time;
double binary_time_prev;
ULL binary_block;

// Enzo to Nbody
Particle* FirstEnzoParticle;
double EnzoLength, EnzoMass, EnzoVelocity, EnzoTime, EnzoForce, EnzoAcceleration;


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
double EnzoElapsedTime; // (SEVN Query) [in Myr unit] I introduced this value for SEVN stellar evolution by EW 2025.3.27
#endif
FILE* workerout;

#ifdef PerformanceTrace
Performance performance;
#endif



double InitialNeighborRadius2;
double EPS2;
int FixNumNeighbor, MaxNumNeighbor;

void DefaultGlobal() {


	/* Task initialization */
	//int Task[NumberOfTask];
	for (int i=0;i<NumberOfTask; i++) {
		Task[i] = i;
	}

	NumberOfCommunication = 0;

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
		EnzoPIDs         = new int[MaxNumberOfParticle];
		AvailableIndices = new int[MaxNumberOfParticle];
		NumberOfAvailableIndices = 0;
		for (int i=0; i<MaxNumberOfParticle; i++) {
			AvailableIndices[i] = -1;
		}
	}

	NumberOfWorker = NumberOfAbyssProcessors-1;

	NewPID = -1;

#ifdef SEVN
	std::vector<std::string> args = {"empty", // Not used
		// "-myself", "/data/vinicius/NbodyPlus/SEVN",
		"-tables", "/data/vinicius/sevn/tables/SEVNtracks_parsec_ov04_AGB", 
		//  "-tables", "/data/vinicius/NbodyPlus/SEVN/tables/SEVNtracks_MIST_AGB",
		// "-tables_HE", "/data/vinicius/NbodyPlus/SEVN/tables/SEVNtracks_parsec_pureHe36",
		// "-turn_WR_to_pureHe", "false",
		"-snmode", "delayed",
		"-Z", "0.0002",
		"-spin", "0.0",
		"-tini", "zams",
		"-tf", "end",
		"-dtout", "events",
		"-xspinmode", "geneva"};
	std::vector<char*> c_args;
	for (auto& arg : args) {
		c_args.push_back(&arg[0]);
	}

	sevnio = new IO; // (SEVN Query) We should initialize sevnio only once. Here in Abyss, and somewhere else in Enzo by EW 2025.3.27
	sevnio->load(c_args.size(), c_args.data());
#endif

}
