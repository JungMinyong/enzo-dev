#ifndef GLOBAL_H
#define GLOBAL_H
#include "def.h"
#include "particle.h"
#include "GlobalVariable.h"
#include "./FewBody/Group.h"
#include "performance.h"
#include <mpi.h>
#include <unordered_map>

#ifdef SEVN
#include "IO.h"
#include <map>
#endif

#include <unordered_set>

const int ROOT = 0;

// Task
const int TASK_TAG = 1;
const int PTCL_TAG = 2;
const int TIME_TAG = 3;
const int QUEUE_TAG = 4;
const int ANY_TAG = 100;
const int TERMINATE_TAG = 666;

/* Communicators */
extern int AbyssProcessorNumber;
extern int NumberOfAbyssProcessors;
extern int NumberOfWorker;

extern MPI_Comm abyss_comm;
extern MPI_Comm inter_comm;

extern MPI_Win win;
extern Particle *particles;

extern MPI_Win win2;
extern GlobalVariable *global_variable;

extern MPI_Win win3;
extern int* Neighbors;

extern MPI_Win win4;
extern int* NewNeighbors;

extern MPI_Datatype QueueType;
extern MPI_Datatype IparticleType;
extern MPI_Datatype JparticleType;

// Particle array
extern int LastParticleIndex;
extern int NumberOfParticle;
extern int NewCMPID;

extern int FixNumNeighbor;

// Parameters for Enzo-Abyss
extern int NumberOfSingleParticle;
extern int *AvailableIndices; // stores available indices due to inactive particles (< LastParticleIndex).
extern int NumberOfAvailableIndices;
extern std::unordered_map<int,int> PIDtoIndexMap;

// Enzo to Abyss
extern double EnzoLength, EnzoMass, EnzoVelocity, EnzoTime, EnzoForce, EnzoAcceleration;
extern double EnzoCurrentTime;
extern double ClusterRadius2;
extern double ClusterAcceleration[Dim];
extern double ClusterPosition[Dim];
extern double ClusterVelocity[Dim];
extern double EnzoClusterPosition[Dim+1];
extern int BinaryRegularization;
extern int IdentifyOnTheFly;
#ifndef INDIVIDUALSTAR
extern int *EnzoPIDs; // stores the order of pids from enzo.
#else
extern int *displs;
#endif
extern int newNumberOfSingleParticle;
extern int StarParticleFeedback;
extern double StarMassEjectionFraction;
extern int ComovingCoordinates;
extern double eta_tmp;
extern double InitialNeighborRadius2;
extern double EPS2;
extern bool OnlyIrregularRoutine;


// Few-Body
extern std::unordered_map<int, int> CMPtclWorker;	   // by EW 2025.1.4 // unordered_map by EW 2025.1.11
extern std::unordered_map<int, int> PrevCMPtclWorker; // by EW 2025.1.4 // unordered_map by EW 2025.1.11



// Time
extern double global_time;
extern double global_time_irr;
extern ULL NextRegTimeBlock;

extern double endTime;




// i/o
extern char* fname;
extern bool restart;
extern char* foutput;
extern double outputTime;
extern int outNum;
extern double outputTimeStep;

extern double AbyssCenter[3];

extern FILE* nbpout;
extern FILE* binout;
extern FILE* mergerout;
#ifdef SEVN
extern FILE* SEVNout;
extern IO* sevnio;
extern std::multimap<double, int> SEVNList;

extern int NumberOfEnzoSEVNParticle;		// This is the number of SEVN particles in Enzo, not in Abyss by EW 2025.4.27
extern int newNumberOfEnzoSEVNParticle;	// This is the number of SEVN particles newly added in Enzo, not in Abyss by EW 2025.4.27
extern std::unordered_map<int,int> PIDtoIndexMap_SEVN;
extern int *EnzoPIDs_SEVN; // stores the order of pids from enzo.

extern std::vector<StarSEVN*> SEVNList_Enzo;
extern std::vector<double> creation_time_Enzo;
extern std::vector<double> world_time_Enzo;
#endif
extern FILE* workerout;

#ifdef PerformanceTrace
// Performance trace
extern Performance performance;
#endif

#endif
