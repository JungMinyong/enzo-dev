#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <cassert>
#include <mpi.h>
#include "global.h"
#include "SkipList.h"
#include "Worker.h"
#include "QueueScheduler.h"

#ifdef NSIGHT
#include <nvToolsExt.h>
#endif

#define noDEBUG

void InitializationRoutines(QueueScheduler &queue_scheduler, Worker *workers);
bool IrregularRoutines(QueueScheduler &queue_scheduler, Worker *workers);
void RegularRoutines(QueueScheduler &queue_scheduler, Worker *workers);
void InitialAssignmentOfTasks(std::vector<int> &data, double next_time, int NumTask, int TAG);
void InitialAssignmentOfTasks(std::vector<int> &data, int NumTask, int TAG);
void InitialAssignmentOfTasks(int data, int NumTask, int TAG);
void InitialAssignmentOfTasks(int *data, int NumTask, int TAG);
void broadcastFromRoot(double &data);
void broadcastFromRoot(ULL &data);
void broadcastFromRoot(int &data);
void ParticleSynchronization();
void updateNextRegTime(std::unordered_set<int> &RegularList);
bool createSkipList(SkipList *skiplist);
bool updateSkipList(SkipList *skiplist, int ptcl_id);
int writeParticle(double current_time, int outputNum);

int SendToEnzo(Worker *workers);
int ReceiveFromEnzo();
void InitializationAfterCommunication(QueueScheduler &queue_scheduler, Worker *workers);
#ifdef SEVN
void StellarEvolution();
#endif

Worker *workers;

void RootRoutines()
{

	std::cout << "Root processor is ready." << std::endl;
	fprintf(nbpout, "Abyss Processor %d is ready.", AbyssProcessorNumber);

	Particle *ptcl;
	// int worker_rank;
	TaskName task;
	int total_tasks;
	int remaining_tasks = 0, completed_tasks = 0, completed_rank;

	std::vector<int> EmptyIndex;				   // by EW 2025.1.7  empty slots in particles e.g., due to mergers
	// unordered_set? by EW 2025.1.11
	// merged particles & PISN will be contained here
	// new single Particle formed in Enzo can be formed in ParticleIndex of these ptcls
	// if empty, LastParticleIndex++

	// MPI_Request requests[NumberOfProcessor];  // Pointer to the request handle
	// MPI_Status statuses[NumberOfProcessor];    // Pointer to the status object
	MPI_Request request; // Pointer to the request handle
	MPI_Status status;	 // Pointer to the status object

	// int sender_rank, sender_tag;
	int ptcl_id;

	workers = new Worker[NumberOfWorker + 1];

	for (int i = 0; i <= NumberOfWorker; i++)
	{
		workers[i].initialize(i);
	}

	QueueScheduler queue_scheduler;

#ifdef PerformanceTrace
	std::chrono::high_resolution_clock::time_point start_point;
	std::chrono::high_resolution_clock::time_point end_point;
#endif

	/* Particle loading Check */
	/*
	{
		//, NextRegTime= %.3e Myr(%llu),
		for (int i=0; i<=LastParticleIndex; i++) {
			ptcl = &particles[i];
			fprintf(stdout, "PID=%d, pos=(%lf, %lf, %lf), vel=(%lf, %lf, %lf)\n",
					ptcl->PID,
					ptcl->Position[0],
					ptcl->Position[1],
					ptcl->Position[2],
					ptcl->Velocity[0],
					ptcl->Velocity[1],
					ptcl->Velocity[2]
					);
		}
		fflush(stdout);
	}*/

	InitializationRoutines(queue_scheduler, workers);

	/* Main Loop */
	while (1)
	{
		// create output at appropriate time intervals
		/*
		if (global_time >= outputTime) {
			writeParticle(global_time, outNum++);
			outputTime += outputTimeStep;
		}
		*/

#ifdef NSIGHT
		nvtxRangePushA("updateNextRegTime");
#endif

		updateNextRegTime(RegularList);

#ifdef NSIGHT
		nvtxRangePop();
#endif
		/*
		std::cout << "NextRegTimeBlock=" << NextRegTimeBlock << std::endl;
		std::cout << "PID= ";
		for (int i : RegularList)
			std::cout << i<< ", ";
		std::cout << std::endl;
		std::cout << "size of regularlist= " << RegularList.size() << std::endl;
		*/	

		if (!IrregularRoutines(queue_scheduler, workers))
			continue;

		RegularRoutines(queue_scheduler, workers);

		global_time = NextRegTimeBlock * global_variable->time_step;
#ifdef SEVN // (Query) EW: PISN should be deleted in PIDtoIndexMap, EnzoPID, ...
		StellarEvolution(); // How about evolving particles inside RegularList only? by EW 2025.1.19
							// Currently, evolving all the particles upto global_time
#endif
		// Time to communicate with enzo
		if (global_time >= 1)
		{
			SendToEnzo(workers);
			ReceiveFromEnzo();
			InitializationAfterCommunication(queue_scheduler, workers);
			NextRegTimeBlock = 0;
			global_time = 0;
		}

	} //  Main loop ends here.
}

void updateNextRegTime(std::unordered_set<int> &RegularList)
{

	ULL time_tmp = 0, time = global_variable->block_max;
	Particle *ptcl = nullptr;

	RegularList.clear();

	for (int i = 0; i <= global_variable->LastParticleIndex; i++)
	{
		// std::cout << i << std::endl;
		ptcl = &particles[i];
		if (!ptcl->isActive)
			continue;
		// Next regular time step
		time_tmp = ptcl->CurrentBlockReg + ptcl->TimeBlockReg;

		// Find the minum regular time step
		if (time_tmp <= time)
		{
			// fprintf(stderr, "PID=%d, time_tme=%llu\n", ptcl->PID, time_tmp);
			if (time_tmp < time)
			{
				RegularList.clear();
				time = time_tmp;
			}
			// RegularList.push_back(ptcl->ParticleIndex);

			// RegularList.insert(PIDtoIndexMap[ptcl->PID]); // If PID is duplicated in ENZO, this causes an error by EW 2025.3.12
			RegularList.insert(ptcl->ParticleIndex);
		}
	}
	NextRegTimeBlock = time;
	global_variable->NextRegTimeBlock = NextRegTimeBlock;
}
