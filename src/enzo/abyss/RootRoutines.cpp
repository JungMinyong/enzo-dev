#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <cassert>
#include <mpi.h>
#include <chrono>
#include "global.h"
#include "QueueScheduler.h"

#ifdef NSIGHT
#include <nvToolsExt.h>
#endif

#define noDEBUG

void InitializationRoutines(QueueScheduler &queue_scheduler, Worker *workers);
bool IrregularRoutines(QueueScheduler &queue_scheduler, Worker *workers, std::unordered_set<int>& RegularList);
void RegularRoutines(QueueScheduler &queue_scheduler, Worker *workers, std::unordered_set<int>& RegularList);
void InitialAssignmentOfTasks(std::vector<int> &data, double next_time, int NumTask, int TAG);
void InitialAssignmentOfTasks(std::vector<int> &data, int NumTask, int TAG);
void InitialAssignmentOfTasks(int data, int NumTask, int TAG);
void InitialAssignmentOfTasks(int *data, int NumTask, int TAG);
void broadcastFromRoot(double &data);
void broadcastFromRoot(ULL &data);
void broadcastFromRoot(int &data);
void ParticleSynchronization();
void updateNextRegTime(std::unordered_set<int> &RegularList);
int writeParticle(double current_time, int outputNum);
int updateParticleBackground();

#ifdef INDIVIDUALSTAR
int SendParticleToEnzo(Worker *workers);
int ReceiveParticleFromEnzo();
#else
int SendToEnzo(Worker *workers);
int ReceiveFromEnzo();
#endif





void InitializationAfterCommunication(QueueScheduler &queue_scheduler, Worker *workers);
#ifdef SEVN
void StellarEvolution();
#endif

Worker *workers;

void RootRoutines()
{
	outNum = 0;
	int countSave = 0; // for StoreTimeStep
	// int MinParticles = 5; // Minimum number of particles to start the nbody routine

	std::cout << "Root processor is ready." << std::endl;
	fprintf(nbpout, "Abyss Processor %d is ready.", AbyssProcessorNumber);

	Particle *ptcl;
	TaskName task;
	int total_tasks;
	int completed_tasks = 0;

	std::chrono::high_resolution_clock::time_point start_point_routine;
	std::chrono::high_resolution_clock::time_point end_point_routine;
	long nbody_durationtime = 0;

	std::unordered_set<int> RegularList;

	QueueScheduler queue_scheduler;
	Queue queue;
	workers = new Worker[NumberOfWorker + 1];

	for (int i = 0; i <= NumberOfWorker; i++) {
		workers[i].initialize(i);
	}

	

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

	if (NumberOfParticle >= MinParticles)
		InitializationRoutines(queue_scheduler, workers);

	/* Main Loop */
	while (1)
	{


		if (NumberOfParticle >= MinParticles) {

			start_point_routine = std::chrono::high_resolution_clock::now();

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
#ifdef DEBUG_ABYSS
			fprintf(nbpout, "Before IrregularRoutines...\n");
			fflush(nbpout);
#endif
			if (!IrregularRoutines(queue_scheduler, workers, RegularList)) {
#ifdef DEBUG_ABYSS
				fprintf(nbpout, "return false in IrregularRoutines...\n");
				fflush(nbpout);
#endif
				end_point_routine = std::chrono::high_resolution_clock::now();
				nbody_durationtime += std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
				continue;
			}
#ifdef DEBUG_ABYSS
			fprintf(nbpout, "After IrregularRoutines...\n");
			fflush(nbpout);
#endif

#ifdef DEBUG_ABYSS
			fprintf(nbpout, "Before RegularRoutines...\n");
			fflush(nbpout);
#endif
			RegularRoutines(queue_scheduler, workers, RegularList); // If OnlyIrregularRoutine is true, this function returns immediately by EW 2025.9.17
#ifdef DEBUG_ABYSS
			fprintf(nbpout, "After RegularRoutines...\n");
			fflush(nbpout);
#endif

			global_time = NextRegTimeBlock * global_variable->time_step;

#ifdef SEVN // (Query) EW: PISN should be deleted in PIDtoIndexMap, EnzoPID, ...
#ifdef DEBUG_ABYSS
			fprintf(nbpout, "Before StellarEvolution...\n");
			fflush(nbpout);
#endif
			if (!SEVNList.empty() && SEVNList.begin()->first <= global_time*global_variable->EnzoTimeStep*1e4 + global_variable->EnzoCurrentTime)
				StellarEvolution(); // Currently, evolving all the particles upto global_time
#ifdef DEBUG_ABYSS
			fprintf(nbpout, "After StellarEvolution...\n");
			fflush(nbpout);
#endif
#endif

			end_point_routine = std::chrono::high_resolution_clock::now();
			nbody_durationtime += std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();

		}
		else  {		// NumberOfParticle < 2 case
					// update position and velocity from background acceleration
					// assume there is no CM particles
				fprintf(stderr, "NumberOfParticle < MinParticles (%d). Only background potential is applied.\n", MinParticles);
#ifdef DEBUG_ABYSS
				fprintf(nbpout, "Before updateParticleBackground...\n");
				fflush(nbpout);
#endif
				updateParticleBackground();
				global_time = 1;
#ifdef DEBUG_ABYSS
				fprintf(nbpout, "After updateParticleBackground...\n");
				fflush(nbpout);
#endif
		}

		// Time to communicate with enzo
		if (global_time >= 1)
		{	
			if (StoreTimeStep > 0){
				countSave++;
				if (countSave >= StoreTimeStep) {
					countSave = 0;
					writeParticle(global_time, outNum++);
				}
			}
		
			fprintf(stderr, "NbodyRoutine: %e (s)\n", nbody_durationtime*1e-9);
			nbody_durationtime = 0;

#ifdef INDIVIDUALSTAR
			start_point_routine = std::chrono::high_resolution_clock::now();
			SendParticleToEnzo(workers);
			end_point_routine = std::chrono::high_resolution_clock::now();
			fprintf(stderr, "SendToEnzo: %e (s)\n", std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count()*1e-9);
			fprintf(nbpout, "SendToEnzo: %e (s)\n", std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count()*1e-9);
			
			start_point_routine = std::chrono::high_resolution_clock::now();
			ReceiveParticleFromEnzo();
			end_point_routine = std::chrono::high_resolution_clock::now();
			fprintf(stderr, "ReceiveFromEnzo: %e (s)\n", std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count()*1e-9);
			fprintf(nbpout, "ReceiveFromEnzo: %e (s)\n", std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count()*1e-9);
#else
			start_point_routine = std::chrono::high_resolution_clock::now();
			SendToEnzo(workers);
			end_point_routine = std::chrono::high_resolution_clock::now();
			fprintf(stderr, "SendToEnzo: %e (s)\n", std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count()*1e-9);

			start_point_routine = std::chrono::high_resolution_clock::now();
			ReceiveFromEnzo();
			end_point_routine = std::chrono::high_resolution_clock::now();
			fprintf(stderr, "ReceiveFromEnzo: %e (s)\n", std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count()*1e-9);
#endif

			start_point_routine = std::chrono::high_resolution_clock::now();
			// we don't need to initialize for N > MinParticles. similar routine has been already implemented for N < 2 within the funciton though
			if (NumberOfParticle >= MinParticles) InitializationAfterCommunication(queue_scheduler, workers);
			end_point_routine = std::chrono::high_resolution_clock::now();
			fprintf(stderr, "InitializationAfterCommunication: %e (s)\n", std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count()*1e-9);
#ifdef DEBUG_ABYSS
			fflush(nbpout);
			fclose(nbpout);
			nbpout = fopen("abyss_output.txt", "w");
			fprintf(nbpout, "Abyss Output Starts!\n");
#endif

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
		if (OnlyIrregularRoutine)
			time_tmp = ptcl->CurrentBlockIrr + ptcl->TimeBlockIrr;
		else
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

