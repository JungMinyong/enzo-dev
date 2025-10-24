#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>
#include <algorithm>
#include "../global.h"
#include "../QueueScheduler.h"
#include "cuda_functions.h"

#include <random>

#ifdef NSIGHT
#include <nvToolsExt.h>
#endif

#define noDEBUG

void sendAllParticlesToGPU(double new_time, const int& RegularListSize, std::vector<int>& RegularListIndices);

void InitializationOnGPU(QueueScheduler &queue_scheduler, Worker *workers);
void sendAllParticlesToGPU_init(std::vector<int>& RegularList_init);

/*
 *  Purporse: calculate acceleration and neighbors of regular particles by sending them to GPU
 *
 *  Date    : 2024.01.18  by Seoyoung Kim
 *
 */
void calculateRegAccelerationOnGPU(std::unordered_set<int>& RegularList, QueueScheduler &queue_scheduler){



	int ListSize = RegularList.size();
	double new_time = NextRegTimeBlock*global_variable->time_step;  // next regular time

	Particle* ptcl;

#ifdef DEBUG_ABYSS
	fprintf(nbpout, "sendAllParticlesToGPU starts\n");
	fflush(nbpout);
#endif

#ifdef NSIGHT
	nvtxRangePushA("sendAllParticlesToGPU");
#endif
	int RegularListSize = RegularList.size();
	std::vector<int> RegularListIndices;
	sendAllParticlesToGPU(new_time, RegularListSize, RegularListIndices);
#ifdef NSIGHT
	nvtxRangePop();
#endif

#ifdef DEBUG_ABYSS
	fprintf(nbpout, "sendAllParticlesToGPU ended\n");
	fflush(nbpout);
#endif
	
	

#ifdef DEBUG_ABYSS
	fprintf(nbpout, "CalculateAccelerationOnDevice starts\n");
	fflush(nbpout);
#endif
  
#ifdef NSIGHT
	nvtxRangePushA("CalculateAccelerationOnDevice");
#endif
  
#ifdef CUDA_FLOAT
	CalculateAccelerationOnDevice(&ListSize, RegularListIndices);
#else
	CalculateAccelerationOnDevice(&ListSize, IndexList, AccRegReceive, AccRegDotReceive, NumNeighborReceive, ACListReceive);
#endif
  
#ifdef NSIGHT
	nvtxRangePop();
#endif
  
#ifdef DEBUG_ABYSS
	fprintf(nbpout, "CalculateAccelerationOnDevice ended\n");
	fflush(nbpout);
#endif
	

#ifdef DEBUG_ABYSS
	fprintf(nbpout, "Adjust Regular Gravity starts\n");
	fflush(nbpout);
#endif

#ifdef NSIGHT
	nvtxRangePushA("RegCuda");
#endif

	queue_scheduler.initialize(RegCuda);
	queue_scheduler.takeQueueRegularList(RegularList);
	do
	{
		queue_scheduler.assignQueueAutoRegularList();
		queue_scheduler.runQueueAuto();
		queue_scheduler.waitQueue(0); // blocking wait
	} while (queue_scheduler.isComplete());

#ifdef NSIGHT
	nvtxRangePop();
#endif

#ifdef DEBUG_ABYSS
	fprintf(nbpout, "Adjust Regular Gravity ended\n");
	fflush(nbpout);
#endif
	//CloseDevice();
} // calculate 0th, 1st derivative of force + neighbors on GPU ends





void sendAllParticlesToGPU(double new_time, const int& RegularListSize, std::vector<int>& RegularListIndices) {

#ifdef COMOVE
	double a = global_variable->a_i + (global_variable->a_f-global_variable->a_i)*(new_time);
#endif

/*
#ifdef PERFORMANCETRACE
	start_point_routine = std::chrono::high_resolution_clock::now();
#endif
*/

	Queue queue = {PrepareGPUCalc, -1, new_time};
	MPI_Request requests[NumberOfWorker];
	for (int i = 0; i < NumberOfWorker; i++)
		MPI_Isend(&queue, 1, QueueType, i+1, QUEUE_TAG, abyss_comm, &requests[i]);

	std::vector<Jparticle> Jparticles;
	Jparticles.resize(NumberOfParticle);

	std::vector<Iparticle> Iparticles;
	Iparticles.resize(RegularListSize);

	RegularListIndices.resize(RegularListSize);

	std::vector<int> counts;
	counts.resize(NumberOfAbyssProcessors * 2);
	int send_buf[2] = {0, 0};

	std::vector<int> Jcounts;
	Jcounts.resize(NumberOfAbyssProcessors);
	Jcounts[0] = 0;

	std::vector<int> Icounts;
	Icounts.resize(NumberOfAbyssProcessors);
	Icounts[0] = 0;

	std::vector<int> Jdispls;
	Jdispls.resize(NumberOfAbyssProcessors);
	Jdispls[0] = 0;

	std::vector<int> Idispls;
	Idispls.resize(NumberOfAbyssProcessors);
	Idispls[0] = 0;

	MPI_Waitall(NumberOfWorker, requests, MPI_STATUSES_IGNORE);

/*
#ifdef PERFORMANCETRACE
	end_point_routine = std::chrono::high_resolution_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end_point_routine - start_point_routine);
	fprintf(stdout, "Send job took %lld microseconds\n", static_cast<long long>(elapsed.count()));
#endif
*/
/*
#ifdef PERFORMANCETRACE
	start_point_routine = std::chrono::high_resolution_clock::now();
#endif
*/

	MPI_Gather(send_buf, 2, MPI_INT, counts.data(), 2, MPI_INT, ROOT, abyss_comm);

	for (int rank = 1; rank < NumberOfAbyssProcessors; rank++) {
	
		Jcounts[rank] = counts[rank * 2 + 0];
		Icounts[rank] = counts[rank * 2 + 1];
	
		Jdispls[rank] = Jdispls[rank - 1] + Jcounts[rank - 1];
		Idispls[rank] = Idispls[rank - 1] + Icounts[rank - 1];
	}
	assert(Jdispls[NumberOfAbyssProcessors - 1] + Jcounts[NumberOfAbyssProcessors - 1] == NumberOfParticle);
	assert(Idispls[NumberOfAbyssProcessors - 1] + Icounts[NumberOfAbyssProcessors - 1] == RegularListSize);

	MPI_Gatherv(nullptr, 0, JparticleType, 
				Jparticles.data(), Jcounts.data(), Jdispls.data(), JparticleType, ROOT, abyss_comm);
	MPI_Gatherv(nullptr, 0, IparticleType, 
				Iparticles.data(), Icounts.data(), Idispls.data(), IparticleType, ROOT, abyss_comm);
	MPI_Gatherv(nullptr, 0, MPI_INT,
				RegularListIndices.data(), Icounts.data(), Idispls.data(), MPI_INT, ROOT, abyss_comm);

/*
#ifdef PERFORMANCETRACE
	end_point_routine = std::chrono::high_resolution_clock::now();
	elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end_point_routine - start_point_routine);
	fprintf(stdout, "Gather calculation took %lld microseconds\n", static_cast<long long>(elapsed.count()));
#endif
*/
	// send the arrays to GPU
/*
#ifdef PERFORMANCETRACE
	start_point_routine = std::chrono::high_resolution_clock::now();
#endif
*/
	SendToDevice(Jparticles, Iparticles);
/*
#ifdef PERFORMANCETRACE
	end_point_routine = std::chrono::high_resolution_clock::now();
	elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end_point_routine - start_point_routine);
	fprintf(stdout, "SendToDevice took %lld microseconds\n", static_cast<long long>(elapsed.count()));
#endif
*/

}
	
void sendAllParticlesToGPU_Worker(double new_time) {

	Particle* ptcl;
	std::vector<Jparticle> Jparticles;
	std::vector<Iparticle> Iparticles;
	std::vector<int> LocalRegularList;

	int J_start = (AbyssProcessorNumber - 1) * (global_variable->LastParticleIndex + 1) / NumberOfWorker;
	int J_end   = AbyssProcessorNumber * (global_variable->LastParticleIndex + 1) / NumberOfWorker;

	// Create a vector of indices from 0 to LastParticleIndex
	std::vector<int> indices(global_variable->LastParticleIndex + 1);
	std::iota(indices.begin(), indices.end(), 0);

	// Shuffle the indices randomly
	// std::random_device rd;
	// std::mt19937 g(rd());
	unsigned int seed = 714;
	std::mt19937 g(seed);
	std::shuffle(indices.begin(), indices.end(), g);

	for (int j = J_start; j < J_end; j++) {

		ptcl = &particles[indices[j]];

		if (!ptcl->isActive)
			continue;

		if (ptcl->NumberOfNeighbor == 0)
			ptcl->predictParticleSecondOrder(new_time-ptcl->CurrentTimeReg, Jparticles, Iparticles, LocalRegularList);
		else
			ptcl->predictParticleSecondOrder(new_time-ptcl->CurrentTimeIrr, Jparticles, Iparticles, LocalRegularList);

	}

	int sizes[2] = {Jparticles.size(), Iparticles.size()};
	MPI_Gather(sizes, 2, MPI_INT, nullptr, 0, MPI_INT, ROOT, abyss_comm);

	MPI_Gatherv(Jparticles.data(), sizes[0], JparticleType,	
				nullptr, nullptr, nullptr, JparticleType, ROOT, abyss_comm);
	MPI_Gatherv(Iparticles.data(), sizes[1], IparticleType,	
				nullptr, nullptr, nullptr, IparticleType, ROOT, abyss_comm);
	MPI_Gatherv(LocalRegularList.data(), sizes[1], MPI_INT,
				nullptr, nullptr, nullptr, MPI_INT, ROOT, abyss_comm);

}


/*
 *  Purporse: calculate acceleration and neighbors of regular particles by sending them to GPU after Enzo-abyss communication
 *
 *  Date    : 2025.03.18  by Eunwoo Chung
 *
 */
void InitializationOnGPU(QueueScheduler &queue_scheduler, Worker *workers) {

	std::vector<int> RegularList_init;
	assert(RegularList_init.empty());

	int ListSize = NumberOfParticle;

#ifdef DEBUG_ABYSS
	fprintf(nbpout, "sendAllParticlesToGPU starts\n");
	fflush(nbpout);
#endif

#ifdef NSIGHT
	nvtxRangePushA("sendAllParticlesToGPU");
#endif
	sendAllParticlesToGPU_init(RegularList_init);  // needs to be updated
	assert(RegularList_init.size() == NumberOfParticle);
#ifdef NSIGHT
	nvtxRangePop();
#endif

#ifdef DEBUG_ABYSS
	fprintf(nbpout, "sendAllParticlesToGPU ended\n");
	fflush(nbpout);
#endif

#ifdef DEBUG_ABYSS
	fprintf(nbpout, "InitializationOnDevice starts\n");
	fflush(nbpout);
#endif
  
#ifdef NSIGHT
	nvtxRangePushA("InitializationOnDevice");
#endif
  
#ifdef CUDA_FLOAT
	// (EW to MY): We should get AccIrrReceive_f, AccIrrDotReceive_f too
	// (EW to MY): Then, we should calculate acc_tot[0], [1] ( acc_irr[0] + acc_reg[0] & acc_irr[1] + acc_reg[1] )
	// (EW to MY): Then, send acc_tot[0] & acc_tot[1] to GPU and calculate acc_irr[2], acc_reg[2], acc_irr[3], acc_reg[3]
	// (EW to MY): Reference: Particle/Initialize.cpp CalculateAcceleration23 function
	// CalculateAccelerationOnDevice(&ListSize, IndexList, AccRegReceive_f, AccRegDotReceive_f, NumNeighborReceive, ACListReceive);
	InitializationOnDevice(&ListSize, RegularList_init);
#endif
  
#ifdef NSIGHT
	nvtxRangePop();
#endif
  
#ifdef DEBUG_ABYSS
	fprintf(nbpout, "InitializationOnDevice ended\n");
	fflush(nbpout);
#endif

	for (const auto& pair: CMPtclWorker) {
		Queue queue;
		int rank = pair.second;
		queue.task = ResetSDARTime;
		queue.pid = pair.first;
		workers[rank].addQueue(queue);
		workers[rank].runQueue();
		workers[rank].callback();
	}
	

#ifdef DEBUG_ABYSS
	fprintf(nbpout, "Adjust Regular Gravity starts\n");
	fflush(nbpout);
#endif

#ifdef NSIGHT
	nvtxRangePushA("RegCuda");
#endif

	// Adjust Regular Gravity
	queue_scheduler.initialize(InitOnGPU);
	queue_scheduler.takeQueue(RegularList_init);
	do
	{
		queue_scheduler.assignQueueAuto();
		queue_scheduler.runQueueAuto();
		queue_scheduler.waitQueue(0); // blocking wait
	} while (queue_scheduler.isComplete());

#ifdef NSIGHT
	nvtxRangePop();
#endif

#ifdef DEBUG_ABYSS
	fprintf(nbpout, "Adjust Regular Gravity ended\n");
	fflush(nbpout);
#endif

	//CloseDevice();
} // calculate 0th, 1st derivative of force + neighbors on GPU ends


void sendAllParticlesToGPU_init(std::vector<int>& RegularList_init) {

/*
#ifdef PERFORMANCETRACE
	start_point_routine = std::chrono::high_resolution_clock::now();
#endif
*/

	Queue queue = {PrepareGPUCalc_init, -1, -1};
	MPI_Request requests[NumberOfWorker];
	for (int i = 0; i < NumberOfWorker; i++)
		MPI_Isend(&queue, 1, QueueType, i+1, QUEUE_TAG, abyss_comm, &requests[i]);

	std::vector<Jparticle> Jparticles;
	Jparticles.resize(NumberOfParticle);

	std::vector<Iparticle> Iparticles;
	Iparticles.resize(NumberOfParticle);

	RegularList_init.resize(NumberOfParticle);

	std::vector<int> counts;
	counts.resize(NumberOfAbyssProcessors);
	counts[0] = 0;
	int send_buf = 0;

	std::vector<int> displs;
	displs.resize(NumberOfAbyssProcessors);
	displs[0] = 0;

	MPI_Waitall(NumberOfWorker, requests, MPI_STATUSES_IGNORE);

/*
#ifdef PERFORMANCETRACE
	end_point_routine = std::chrono::high_resolution_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end_point_routine - start_point_routine);
	fprintf(stdout, "Send job took %lld microseconds\n", static_cast<long long>(elapsed.count()));
#endif
*/
/*
#ifdef PERFORMANCETRACE
	start_point_routine = std::chrono::high_resolution_clock::now();
#endif
*/

	MPI_Gather(&send_buf, 1, MPI_INT, counts.data(), 1, MPI_INT, ROOT, abyss_comm);

	for (int rank = 1; rank < NumberOfAbyssProcessors; rank++)
		displs[rank] = displs[rank - 1] + counts[rank - 1];

	assert(displs[NumberOfAbyssProcessors - 1] + counts[NumberOfAbyssProcessors - 1] == NumberOfParticle);

	MPI_Gatherv(nullptr, 0, JparticleType, 
				Jparticles.data(), counts.data(), displs.data(), JparticleType, ROOT, abyss_comm);
	MPI_Gatherv(nullptr, 0, IparticleType, 
				Iparticles.data(), counts.data(), displs.data(), IparticleType, ROOT, abyss_comm);
	MPI_Gatherv(nullptr, 0, MPI_INT,
				RegularList_init.data(), counts.data(), displs.data(), MPI_INT, ROOT, abyss_comm);

/*
#ifdef PERFORMANCETRACE
	end_point_routine = std::chrono::high_resolution_clock::now();
	elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end_point_routine - start_point_routine);
	fprintf(stdout, "Gather calculation took %lld microseconds\n", static_cast<long long>(elapsed.count()));
#endif
*/
	// send the arrays to GPU
/*
#ifdef PERFORMANCETRACE
	start_point_routine = std::chrono::high_resolution_clock::now();
#endif
*/
	SendToDevice(Jparticles, Iparticles);
/*
#ifdef PERFORMANCETRACE
	end_point_routine = std::chrono::high_resolution_clock::now();
	elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end_point_routine - start_point_routine);
	fprintf(stdout, "SendToDevice took %lld microseconds\n", static_cast<long long>(elapsed.count()));
#endif
*/
}
		
void sendAllParticlesToGPU_init_Worker() {

	Particle* ptcl;
	std::vector<Jparticle> Jparticles;
	std::vector<Iparticle> Iparticles;
	std::vector<int> LocalRegularList;

	int J_start = (AbyssProcessorNumber - 1) * (global_variable->LastParticleIndex + 1) / NumberOfWorker;
	int J_end   = AbyssProcessorNumber * (global_variable->LastParticleIndex + 1) / NumberOfWorker;

	// Create a vector of indices from 0 to LastParticleIndex
	std::vector<int> indices(global_variable->LastParticleIndex + 1);
	std::iota(indices.begin(), indices.end(), 0);

	// Shuffle the indices randomly
	// std::random_device rd;
	// std::mt19937 g(rd());
	unsigned int seed = 714;
	std::mt19937 g(seed);
	std::shuffle(indices.begin(), indices.end(), g);

	for (int j = J_start; j < J_end; j++) {

		ptcl = &particles[indices[j]];
		Iparticle iptcl;
		Jparticle jptcl;

		if (ptcl->TimeStepIrr != 0)
			ptcl->setNewTimeStepWithNewEnzoTimeStep(global_variable->OldEnzoTimeStep, global_variable->EnzoTimeStep);

		ptcl->CurrentTimeIrr = 0.;
		ptcl->CurrentBlockIrr = 0;
		ptcl->CurrentTimeReg = 0.;
		ptcl->CurrentBlockReg = 0;
		ptcl->NewCurrentBlockIrr = 0;

		if (!ptcl->isActive)
			continue;

		LocalRegularList.push_back(ptcl->ParticleIndex);
		jptcl.posx		= static_cast<CUDA_REAL>(ptcl->Position[0]);
		jptcl.posy		= static_cast<CUDA_REAL>(ptcl->Position[1]);
		jptcl.posz		= static_cast<CUDA_REAL>(ptcl->Position[2]);
		jptcl.velx		= static_cast<CUDA_REAL>(ptcl->Velocity[0]);
		jptcl.vely		= static_cast<CUDA_REAL>(ptcl->Velocity[1]);
		jptcl.velz		= static_cast<CUDA_REAL>(ptcl->Velocity[2]);
		jptcl.mass		= static_cast<CUDA_REAL>(ptcl->Mass);
		jptcl.index		= ptcl->ParticleIndex;
		Jparticles.push_back(jptcl);

		iptcl.posx		= jptcl.posx;
		iptcl.posy		= jptcl.posy;
		iptcl.posz		= jptcl.posz;
		iptcl.velx		= jptcl.velx;
		iptcl.vely		= jptcl.vely;
		iptcl.velz		= jptcl.velz;
		iptcl.r2		= static_cast<CUDA_REAL>(ptcl->RadiusOfNeighbor);
		iptcl.dtr		= static_cast<CUDA_REAL>(ptcl->TimeBlockReg*global_variable->time_step*global_variable->EnzoTimeStep);
		Iparticles.push_back(iptcl);
	}

	int size = LocalRegularList.size();

	MPI_Gather(&size, 1, MPI_INT, nullptr, 0, MPI_INT, ROOT, abyss_comm);

	MPI_Gatherv(Jparticles.data(), size, JparticleType,	
				nullptr, nullptr, nullptr, JparticleType, ROOT, abyss_comm);
	MPI_Gatherv(Iparticles.data(), size, IparticleType,	
				nullptr, nullptr, nullptr, IparticleType, ROOT, abyss_comm);
	MPI_Gatherv(LocalRegularList.data(), size, MPI_INT,
				nullptr, nullptr, nullptr, MPI_INT, ROOT, abyss_comm);

}