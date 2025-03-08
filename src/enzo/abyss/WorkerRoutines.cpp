#include <iostream>
#include <vector>
#include <errno.h>
#include "global.h"
#include "Queue.h"
#include <cassert>

void ComputeAcceleration(int ptcl_index, double next_time);
void broadcastFromRoot(double &data);
void broadcastFromRoot(ULL &data);
void broadcastFromRoot(int &data);
void CalculateAcceleration01(Particle* ptcl1);
void CalculateAcceleration23(Particle* ptcl1);
void makePrimordialGroup(Particle* ptclCM);
void NewFBInitialization(Particle* ptclCM);
void deleteGroup(Particle* ptclCM);
void NewFBInitialization3(Group* group);

void WorkerRoutines() {

	std::cout << "Processor " << AbyssProcessorNumber << " is ready." << std::endl;
	fprintf(nbpout, "Abyss Processor %d is ready.", AbyssProcessorNumber);

	TaskName task = Error;
	MPI_Status status;
	MPI_Request request;
	int ptcl_index;
	double next_time;
	int NewNumberOfNeighbor;
	int NewNeighbors[NumNeighborMax];
	int size=0;
	double new_a[Dim];
	double new_adot[Dim];
	Particle *ptcl;
	std::chrono::high_resolution_clock::time_point start_point;
	std::chrono::high_resolution_clock::time_point end_point;

	while (true) {
		MPI_Recv(&task, 1, MPI_INT, ROOT, TASK_TAG, abyss_comm, &status);
		//MPI_Irecv(&task, 1, MPI_INT, ROOT, TASK_TAG, abyss_comm, &request);
		//MPI_Wait(&request, &status);
		//if (status.MPI_TAG == TERMINATE_TAG) break;
		//std::cerr << "Processor " << AbyssProcessorNumber << " received task " << task << std::endl;

		switch (task) {
			case IrrForce: // Irregular Acceleration
				MPI_Recv(&ptcl_index,   1, MPI_INT   , ROOT, PTCL_TAG, abyss_comm, &status);
				//std::cout << "(IRR_FORCE) Processor " << AbyssProcessorNumber<< ": PID= "<<ptcl_index << std::endl;
				MPI_Recv(&next_time, 1, MPI_DOUBLE, ROOT, TIME_TAG, abyss_comm, &status); // (Query to myself) it seems like it's not needed.
#ifdef PerformanceTrace
				ptcl = &particles[ptcl_index];
				start_point = std::chrono::high_resolution_clock::now();
				ptcl->computeAccelerationIrr();
				end_point = std::chrono::high_resolution_clock::now();
				performance.IrregularForce +=
					std::chrono::duration_cast<std::chrono::nanoseconds>(end_point - start_point).count();
#else
				ptcl->computeAccelerationIrr();
#endif

				ptcl->NewCurrentBlockIrr = ptcl->CurrentBlockIrr + ptcl->TimeBlockIrr; // of this particle
				ptcl->calculateTimeStepIrr();
				ptcl->NextBlockIrr = ptcl->NewCurrentBlockIrr + ptcl->TimeBlockIrr; // of this particle
				ptcl->isUpdateToDate = true;
				//std::cout << "IrrCal done " << AbyssProcessorNumber << std::endl;
				break;

			case RegForce: // Regular Acceleration
				//std::cout << "RegCal start " << AbyssProcessorNumber << std::endl;
				MPI_Recv(&ptcl_index,   1, MPI_INT,    ROOT, PTCL_TAG, abyss_comm, &status);
				MPI_Recv(&next_time, 1, MPI_DOUBLE, ROOT, TIME_TAG, abyss_comm, &status);

				particles[ptcl_index].computeAccelerationReg();
				//ComputeAcceleration(ptcl_index, next_time);
				//std::cout << "RegCal end" << AbyssProcessorNumber << std::endl;
				break;

			case IrrUpdate: // Irregular Update Particle
				//std::cout << "IrrUp Processor " << AbyssProcessorNumber << std::endl;
				MPI_Recv(&ptcl_index  , 1, MPI_INT   , ROOT, PTCL_TAG, abyss_comm, &status);
				ptcl = &particles[ptcl_index];

				if (ptcl->NumberOfNeighbor != 0) // IAR modified
					ptcl->updateParticle();
				ptcl->CurrentBlockIrr = ptcl->NewCurrentBlockIrr;
				ptcl->CurrentTimeIrr  = ptcl->CurrentBlockIrr*global_variable->time_step;
				//std::cout << "pid=" << ptcl_index << ", CurrentBlockIrr=" << particles[ptcl_index].CurrentBlockIrr << std::endl;
				//std::cout << "IrrUp end " << AbyssProcessorNumber << std::endl;
				break;

			case RegUpdate: // Regular Update Particle
				//std::cout << "RegUp start " << AbyssProcessorNumber << std::endl;
				MPI_Recv(&ptcl_index, 1, MPI_INT, ROOT, PTCL_TAG, abyss_comm, &status);
				//std::cout << "ptcl " << ptcl_index << std::endl;

				ptcl = &particles[ptcl_index];
				ptcl->updateParticle();

				for (int i=0; i<ptcl->NewNumberOfNeighbor; i++)
					ptcl->Neighbors[i] = ptcl->NewNeighbors[i];
				ptcl->NumberOfNeighbor = ptcl->NewNumberOfNeighbor;

				ptcl->CurrentBlockReg += ptcl->TimeBlockReg;
				ptcl->CurrentTimeReg   = ptcl->CurrentBlockReg*global_variable->time_step;
				ptcl->calculateTimeStepReg();
				ptcl->NewCurrentBlockIrr = ptcl->CurrentBlockReg;
				ptcl->calculateTimeStepIrr();
				ptcl->updateRadius();
				if (ptcl->NumberOfNeighbor == 0) {
					ptcl->CurrentBlockIrr = ptcl->CurrentBlockReg;
					ptcl->CurrentTimeIrr = ptcl->CurrentBlockReg*global_variable->time_step;
				}
				ptcl->NextBlockIrr = ptcl->CurrentBlockIrr + ptcl->TimeBlockIrr; // of this particle
				break;

			case RegCuda: // Update Regular Particle CUDA
				//std::cout << "(REG_CUDA) Processor " << std::endl;
				MPI_Recv(&ptcl_index, 1, MPI_INT, ROOT, PTCL_TAG, abyss_comm, MPI_STATUS_IGNORE);
				//std::cout << "(REG_CUDA) Processor " << AbyssProcessorNumber<< ": Particle Index= "<<ptcl_index << std::endl;
				MPI_Recv(&NewNumberOfNeighbor, 1, MPI_INT, ROOT, 10, abyss_comm, &status);
				MPI_Recv(NewNeighbors, NewNumberOfNeighbor, MPI_INT, ROOT, 11, abyss_comm, &status);
				MPI_Recv(new_a, 3, MPI_DOUBLE, ROOT, 12, abyss_comm, &status);
				MPI_Recv(new_adot, 3, MPI_DOUBLE, ROOT, 13, abyss_comm, &status);
				particles[ptcl_index].updateRegularParticleCuda(NewNeighbors, NewNumberOfNeighbor, new_a, new_adot);
				break;

			case RegCudaUpdate: // Update Regular Particle CUDA II
				MPI_Recv(&ptcl_index, 1, MPI_INT, ROOT, PTCL_TAG, abyss_comm, MPI_STATUS_IGNORE);
				//std::cout << "(REG_UPDATE) Processor " << AbyssProcessorNumber<< ": PID= "<<ptcl_index << std::endl;
				ptcl = &particles[ptcl_index];

				for (int j = 0; j < ptcl->NewNumberOfNeighbor; j++)
					ptcl->Neighbors[j] = ptcl->NewNeighbors[j];
				ptcl->NumberOfNeighbor = ptcl->NewNumberOfNeighbor;

				ptcl->updateParticle();
				ptcl->CurrentBlockReg = ptcl->CurrentBlockReg + ptcl->TimeBlockReg;
				ptcl->CurrentTimeReg = ptcl->CurrentBlockReg * global_variable->time_step;
				ptcl->calculateTimeStepReg();
				ptcl->calculateTimeStepIrr();
				if (ptcl->NumberOfNeighbor == 0) {
					//ptcl->CurrentBlockIrr = ptcl->CurrentBlockReg;
					//ptcl->CurrentTimeIrr = ptcl->CurrentBlockReg*time_step;
					if (ptcl->CurrentBlockIrr != ptcl->CurrentBlockReg || ptcl->CurrentTimeIrr != ptcl->CurrentBlockReg*global_variable->time_step) {
						fprintf(stderr, "PID: %d\n", ptcl->PID);
						fprintf(stderr, "TimeBlockIrr: %llu, TimeBlockReg: %llu\n", ptcl->TimeBlockIrr, ptcl->TimeBlockReg);
						fprintf(stderr, "CurrentBlockIrr: %llu, CurrentBlockReg: %llu\n", ptcl->CurrentBlockIrr, ptcl->CurrentBlockReg);
						fprintf(stderr, "CurrentTimeIrr: %e, CurrentTimeReg: %e\n", ptcl->CurrentTimeIrr, ptcl->CurrentTimeReg);
						fprintf(stderr, "NextRegTimeBlock: %llu\n", global_variable->NextRegTimeBlock);
						fflush(stderr);
						assert(ptcl->CurrentBlockIrr == ptcl->CurrentBlockReg);
						assert(ptcl->CurrentTimeIrr == ptcl->CurrentBlockReg*global_variable->time_step);
					}
				}
				ptcl->updateRadius();
				ptcl->NextBlockIrr = ptcl->CurrentBlockIrr + ptcl->TimeBlockIrr; // of ptcl particle
				break;

			case InitAcc1: // Initialize Acceleration(01)
				//std::cout << "Processor " << AbyssProcessorNumber<< " initialization starts." << std::endl;
				MPI_Recv(&ptcl_index, 1, MPI_INT, ROOT, PTCL_TAG, abyss_comm, &status);
				//std::cout << "Processor " << AbyssProcessorNumber<< ": PID= "<<ptcl_index << std::endl;
				ptcl = &particles[ptcl_index];
				//std::cerr << ptcl_index+i << std::endl;
				CalculateAcceleration01(ptcl);
				//std::cout << "Processor " << AbyssProcessorNumber<< " done." << std::endl;
				break;

			case InitAcc2: // Initialize Acceleration(23)
				MPI_Recv(&ptcl_index, 1, MPI_INT, ROOT, PTCL_TAG, abyss_comm, &status);
				//std::cout << "Processor " << AbyssProcessorNumber<< ": PID= "<<ptcl_index << std::endl;
				ptcl = &particles[ptcl_index];
				CalculateAcceleration23(ptcl);
				break;

			case InitTime: // Initialize Time Step
				MPI_Recv(&ptcl_index, 1, MPI_INT, ROOT, PTCL_TAG, abyss_comm, &status);
				//std::cout << "Processor " << AbyssProcessorNumber<< ": PID= "<<ptcl_index << std::endl;
				ptcl = &particles[ptcl_index];
				//if (ptcl->isActive)
				ptcl->initializeTimeStep();
				break;

			case TimeSync: // Initialize Timestep variables
				fprintf(stderr, "time sync (%d)\n", AbyssProcessorNumber);
				//broadcastFromRoot(time_block);
				//broadcastFromRoot(block_max);
				//broadcastFromRoot(time_step);
				//MPI_Win_sync(win);  // Synchronize memory
				//MPI_Barrier(abyss_comm);
				//MPI_Win_fence(0, win);
				fprintf(stderr, "(%d) nbody+:time_block = %d, EnzoTimeStep=%e\n", AbyssProcessorNumber, global_variable->time_block, global_variable->EnzoTimeStep);
				fflush(stderr);
				break;


#ifdef FEWBODY
			case SearchPrimordialGroup: // Primordial binary search
				MPI_Recv(&ptcl_index, 1, MPI_INT, ROOT, PTCL_TAG, abyss_comm, &status);
				ptcl = &particles[ptcl_index];

				ptcl->NewNumberOfNeighbor = 0;
				ptcl->checkNewGroup2();

				break;

			case SearchGroup: // Few-body group search
				MPI_Recv(&ptcl_index  , 1, MPI_INT   , ROOT, PTCL_TAG, abyss_comm, &status);
				ptcl = &particles[ptcl_index];
				// std::cerr << "FB search of particle  " << ptcl_index << " is initiated on rank " << AbyssProcessorNumber << "." <<std::endl;

				if (ptcl->getBinaryInterruptState()==BinaryInterruptState::threebody) {
					ptcl->setBinaryInterruptState(BinaryInterruptState::none);
				}
				else if (ptcl->getBinaryInterruptState()==BinaryInterruptState::manybody) {
					ptcl->NewNumberOfNeighbor = 0;
					ptcl->checkNewGroup2();
					ptcl->setBinaryInterruptState(BinaryInterruptState::none);
				}
				else {
					ptcl->NewNumberOfNeighbor = 0;
					if (ptcl->TimeStepIrr*global_variable->EnzoTimeStep*1e4 < TSEARCH)
						ptcl->checkNewGroup();
				}
				/*
				if (ptcl->getBinaryInterruptState()==BinaryInterruptState::manybody) {
					ptcl->setBinaryInterruptState(BinaryInterruptState::none);
					std::cout << "ptcl PID: " << ptcl->PID << ", ptcl NewNumberOfNeighbor: " << ptcl->NewNumberOfNeighbor << std::endl;
				}
				else {
					ptcl->NewNumberOfNeighbor = 0;
					if (ptcl->TimeStepIrr*EnzoTimeStep*1e4 < TSEARCH)
						ptcl->checkNewGroup();
				}
				*/
				// std::cerr << "FB search of particle  " << ptcl_index << " is successfully finished on rank " << AbyssProcessorNumber << "." <<std::endl;

				break;

			case MakePrimordialGroup: // Make a primordial group
				MPI_Recv(&ptcl_index  , 1, MPI_INT   , ROOT, PTCL_TAG, abyss_comm, &status);
				ptcl = &particles[ptcl_index];

				makePrimordialGroup(ptcl);

				break;

			case MakeGroup: // Make a group
				MPI_Recv(&ptcl_index  , 1, MPI_INT   , ROOT, PTCL_TAG, abyss_comm, &status);
				ptcl = &particles[ptcl_index];

				NewFBInitialization(ptcl);
#ifdef DEBUG
				std::cout << "FewBody object of particle " << ptcl->PID
						  << " is successfully initialized on rank " << AbyssProcessorNumber << "." <<std::endl;
#endif
				break;

			case DeleteGroup: // Delete a Group struct
				MPI_Recv(&ptcl_index  , 1, MPI_INT   , ROOT, PTCL_TAG, abyss_comm, &status);
				ptcl = &particles[ptcl_index];
				deleteGroup(ptcl);
				break;

			case ARIntegration: // SDAR for few body encounters
				MPI_Recv(&ptcl_index,   1, MPI_INT   , ROOT, PTCL_TAG, abyss_comm, &status);
				MPI_Recv(&next_time, 1, MPI_DOUBLE, ROOT, TIME_TAG, abyss_comm, &status);
				
				ptcl = &particles[ptcl_index];
				// std::cout << "(SDAR) Processor " << AbyssProcessorNumber<< ": PID= "<<ptcl->PID << std::endl;

				/* (Query) this will be done already. 
				ptcl->computeAccelerationIrr();
				ptcl->NewCurrentBlockIrr = ptcl->CurrentBlockIrr + ptcl->TimeBlockIrr; // of this particle
				ptcl->calculateTimeStepIrr();
				ptcl->NextBlockIrr = ptcl->NewCurrentBlockIrr + ptcl->TimeBlockIrr; // of this particle
				*/

				if (!ptcl->isCMptcl || ptcl->GroupInfo == nullptr) {
					fprintf(stderr, "Something is wrong. ptcl->isCMptcl=%d ptcl->GroupInfo=%p\n", ptcl->isCMptcl, ptcl->GroupInfo);
					exit(EXIT_FAILURE);
				}
				
				ptcl->GroupInfo->ARIntegration(next_time);
				if (!ptcl->GroupInfo->isMerger && !ptcl->GroupInfo->isTerminate)
					ptcl->GroupInfo->isTerminate = ptcl->GroupInfo->CheckBreak();

				if (ptcl->GroupInfo->isTerminate) {
					if (ptcl->getBinaryInterruptState() == BinaryInterruptState::none)
						ptcl->setBinaryInterruptState(BinaryInterruptState::terminated);

					delete ptcl->GroupInfo;
#ifdef DEBUG
					std::cout << "(SDAR) Processor " << AbyssProcessorNumber<< ": PID= "<<ptcl->PID << " deleted!" <<std::endl;
#endif
				}
#ifdef DEBUG
				else
					std::cout << "(SDAR) Processor " << AbyssProcessorNumber<< ": PID= "<<ptcl->PID << " done!" <<std::endl;
#endif
				break;
			
			case MergeManyBody: // Merger insided many-body (>2) group

				MPI_Recv(&ptcl_index,   1, MPI_INT   , ROOT, PTCL_TAG, abyss_comm, &status);
				
				ptcl = &particles[ptcl_index];
				std::cout << "(SDAR) Processor " << AbyssProcessorNumber<< ": PID= "<<ptcl->PID << std::endl;

				if (!ptcl->isCMptcl || ptcl->GroupInfo == nullptr) {
					fprintf(stderr, "Something is wrong. ptcl->isCMptcl=%d ptcl->GroupInfo=%p\n", ptcl->isCMptcl, ptcl->GroupInfo);
					exit(EXIT_FAILURE);
				}

				NewFBInitialization3(ptcl->GroupInfo);

				ptcl->GroupInfo->isMerger = false;
				ptcl->setBinaryInterruptState(BinaryInterruptState::none);

				std::cout << "(SDAR) Processor " << AbyssProcessorNumber<< ": PID= "<<ptcl->PID << " NewFBInitialization3 done!" <<std::endl;
				break;
#endif 

			case Synchronize: // Synchronize
				MPI_Win_sync(win);  // Synchronize memory
				MPI_Barrier(abyss_comm);
				break;

			case Ends: // Simualtion ends
				std::cout << "Processor " << AbyssProcessorNumber<< " returns." << std::endl;
				return;
				break;

			case Error:
				perror("Error task assignments");
				exit(EXIT_FAILURE);
				break;
			default:
				break;
		}

		// return that it's over
		//task = -1;
		if (task == IrrForce || task == RegForce || task == IrrUpdate || task == RegUpdate)
			MPI_Isend(&ptcl_index, 1, MPI_INT, ROOT, TERMINATE_TAG, abyss_comm,&request);
		else
			MPI_Isend(&task, 1, MPI_INT, ROOT, TERMINATE_TAG, abyss_comm,&request);

		MPI_Wait(&request, &status);
		//std::cerr << "Processor " << AbyssProcessorNumber << " done." << std::endl;
	}
}


void ComputeAcceleration(int ptcl_index, double next_time) {
	Particle *ptcl = &particles[ptcl_index];
	Particle *neighbor;
	double neighbor_pos[Dim], neighbor_vel[Dim];
	for (int i=0; i<ptcl->NumberOfNeighbor; i++) {
		neighbor = &particles[ptcl->Neighbors[i]];
		//neighbor->predictParticle(next_time, neighbor_pos, neighbor_vel);
		//ptcl->getAcceleration(neighbor_pos, neighbor_vel);
	}
}

