#ifdef SEVN
#include "sevn.h"
#endif

#include <iostream>
#include <fstream>
#include <unistd.h>
#include "def.h"
#include "particle.h"
#include "GlobalVariable.h"
#include "global.h"
#include <mpi.h>
#ifdef CUDA
#include <cuda_runtime.h>
#include "cuda/cuda_functions.h"
#endif

#include <sys/stat.h>
#include <sys/types.h>
#include <cerrno>
#include <cstring>


void broadcastFromRoot(int &data);
void broadcastFromRoot(double &data);
void DefaultGlobal();
void WorkerRoutines();
void RootRoutines();
int InitialCommunication();
bool directoryExists(const std::string &path);

MPI_Datatype createQueueType();
MPI_Datatype createIparticleType();
MPI_Datatype createJparticleType();

int ABYSS() {

	fprintf(stdout, "Abyss Starts!\n");

	/* Initialize global variables */
	DefaultGlobal();


	/*********************************************************************
	 *  Configuration for outputing log files
	 *********************************************************************/
	nbpout = fopen("abyss_output.txt", "w");
	//gpuout = fopen("cuda_output.txt", "w");
	fprintf(nbpout, "Abyss Output Starts!\n");
	//fprintf(gpuout, "CUDA Output Starts!\n");

	binout = fopen("binary_output.txt", "w"); // (Query) EW: how to open output files?
	fprintf(binout, "Starting nbody - Binary OUTPUT\n");
	fflush(binout);
	mergerout = fopen("merger_output.txt", "w");
	fprintf(mergerout, "Starting nbody - Merger OUTPUT\n");
	fflush(mergerout);
#ifdef SEVN
	SEVNout = fopen("SEVN_output.txt", "w");
	fprintf(SEVNout, "Starting nbody - SEVN OUTPUT\n");
	fflush(SEVNout);
#endif

	QueueType		= createQueueType();
	IparticleType	= createIparticleType();
	JparticleType	= createJparticleType();

#ifdef CUDA
	if (AbyssProcessorNumber == ROOT) {
		OpenDevice();
		cudaDeviceSynchronize();
	}
#endif

	if (AbyssProcessorNumber == ROOT)
		InitialCommunication();

	// things that should be synchronized. this can be moved to GlobalVariable
	MPI_Barrier(abyss_comm);
	broadcastFromRoot(EnzoMass);
	broadcastFromRoot(EnzoLength);
	broadcastFromRoot(EnzoVelocity);
	broadcastFromRoot(EnzoTime);
	broadcastFromRoot(EnzoAcceleration);
	broadcastFromRoot(EPS2);
	broadcastFromRoot(InitialNeighborRadius2);
	broadcastFromRoot(FixNumNeighbor);

	if (AbyssProcessorNumber == ROOT) {
		RootRoutines();
	} else {
		// /* // by EW 2025.1.27
		std::string filename = "log/worker/worker_output_" + std::to_string(AbyssProcessorNumber) + ".txt";
		std::string dir_name;

		dir_name = "log";
		if (!directoryExists(dir_name))
		{
			// Create directory with permission 0755
			if (mkdir(dir_name.c_str(), 0755) == 0)
			{
				std::cout << "Directory created successfully." << std::endl;
			}
			else
			{
				std::cerr << "Failed to create directory: "
						<< std::strerror(errno) << std::endl;
			}
		}
		dir_name = "log/worker";
		if (!directoryExists(dir_name))
		{
			// Create directory with permission 0755
			if (mkdir(dir_name.c_str(), 0755) == 0)
			{
				std::cout << "Directory created successfully." << std::endl;
			}
			else
			{
				std::cerr << "Failed to create directory: "
						<< std::strerror(errno) << std::endl;
			}
		}

		workerout = fopen(filename.c_str(), "w");
		fprintf(workerout, "Starting nbody - WORKER OUTPUT\n");
		fflush(workerout);
		// */
		WorkerRoutines();
	}

	return 0;
}


bool directoryExists(const std::string &path) {
    struct stat info;
    // stat returns 0 if the path exists
    if (stat(path.c_str(), &info) != 0) {
        return false;
    }
    // Check if it's a directory
    return (info.st_mode & S_IFDIR) != 0;
}