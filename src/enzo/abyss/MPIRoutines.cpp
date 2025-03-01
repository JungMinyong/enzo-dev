#include <iostream>
#include <vector>
#include <unistd.h>
#include "def.h"
#include "particle.h"
#include "GlobalVariable.h"
#include "global.h"
#include <mpi.h>


template <typename T>
void InitialAssignmentOfTasks(T data, int NumTask);
void InitialAssignmentOfTasks(std::vector<int>& data, int NumTask, int TAG) {
	for (int i=0; i<NumberOfWorker; i++) {
		if (i >= NumTask) break;
		//MPI_Isend(&data[i],   1, MPI_INT,    i+1, TAG, abyss_comm, &requests[NumberOfCommunication++]);
		MPI_Send(&data[i],   1, MPI_INT,    i+1, TAG, abyss_comm);
	}
}


void InitialAssignmentOfTasks(std::vector<int>& data, double next_time, int NumTask, int TAG) {
	for (int i=0; i<NumberOfWorker; i++) {
		if (i >= NumTask) break;
		//MPI_Isend(&data[i],   1, MPI_INT,    i+1, PTCL_TAG, abyss_comm, &requests[NumberOfCommunication++]);
		//MPI_Isend(&next_time, 1, MPI_DOUBLE, i+1, TIME_TAG, abyss_comm, &requests[NumberOfCommunication++]);

		MPI_Send(&data[i],   1, MPI_INT,    i+1, PTCL_TAG, abyss_comm);
		MPI_Send(&next_time, 1, MPI_DOUBLE, i+1, TIME_TAG, abyss_comm);
	}
}


void InitialAssignmentOfTasks(int* data, int NumTask, int TAG) {
	for (int i=0; i<NumberOfWorker; i++) {
		if (i >= NumTask) break;
		//MPI_Isend(&data[i], 1, MPI_INT, i+1, TAG, abyss_comm, &requests[NumberOfCommunication++]);
		MPI_Send(&data[i], 1, MPI_INT, i+1, TAG, abyss_comm);
	}
}

void InitialAssignmentOfTasks(int data, int NumTask, int TAG) {
	for (int i=0; i<NumberOfWorker; i++) {
		//std::cerr << "InitialAssignmentOfTasks out of" << NumTask<< ": " << i << std::endl;
		if (i >= NumTask)
			break;
		//MPI_Isend(&data, 1, MPI_INT, i+1, TAG, abyss_comm, &requests[NumberOfCommunication++]);
		MPI_Send(&data, 1, MPI_INT, i+1, TAG, abyss_comm);
	}
	//fprintf(stderr, "Number of tasks assigned = %d\n", i);
	//fflush(stderr);
}




void broadcastFromRoot(int &data) {
	MPI_Bcast(&data, 1, MPI_INT, ROOT, abyss_comm);
}
void broadcastFromRoot(double &data) {
	MPI_Bcast(&data, 1, MPI_DOUBLE, ROOT, abyss_comm);
}
void broadcastFromRoot(ULL &data) {
	MPI_Bcast(&data, 1, MPI_UNSIGNED_LONG_LONG, ROOT, abyss_comm);
}

void ParticleSynchronization() {
		//std::cout << "Particle synchronization." << std::endl;
		MPI_Request requests[NumberOfAbyssProcessors];  // Pointer to the request handle
		MPI_Status statuses[NumberOfAbyssProcessors];    // Pointer to the status object
		int task=100;
		int completed_tasks = 0;
		//int ptcl_id=104;
		InitialAssignmentOfTasks(task, NumberOfWorker, TASK_TAG);
		//std::cerr << "before, Rank=" << MyRank <<" pid=" << ptcl_id << ", current_time=" << particles[ptcl_id].CurrentTimeIrr << std::endl;
		MPI_Win_sync(win);  // Synchronize memory
		MPI_Barrier(abyss_comm);
		//std::cerr << "after, Rank=" << MyRank <<" pid=" << ptcl_id << ", current_time=" << particles[ptcl_id].CurrentTimeIrr << std::endl;
		for (int i=0; i<NumberOfWorker; i++) {

			MPI_Irecv(&task, 1, MPI_INT, MPI_ANY_SOURCE, TERMINATE_TAG, abyss_comm, &requests[i]);
		}
		MPI_Waitall(NumberOfWorker, requests, statuses);
		//std::cerr << "Done" << std::endl;
}

void FurtherAssignmentOfTasks() {

}

