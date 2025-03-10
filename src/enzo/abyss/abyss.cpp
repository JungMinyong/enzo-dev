#include <sys/stat.h>
#include <sys/types.h>
#include <cerrno>
#include <cstring>
#include <iostream>
#include <fstream>
#include <iostream>
#include <vector>
#include <unistd.h>
#include "def.h"
#include "particle.h"
#include "GlobalVariable.h"
#include "global.h"
#include <mpi.h>
#include <unistd.h>
#include <cuda_runtime.h>
#include "cuda/cuda_functions.h"

#ifdef SEVN
#include "sevn.h"
#endif

void broadcastFromRoot(int &data);
void broadcastFromRoot(double &data);
void DefaultGlobal();
void WorkerRoutines();
void RootRoutines();
int InitialCommunication();
bool directoryExists(const std::string &path);

int ABYSS() {

	fprintf(stdout, "Abyss Starts!\n");

	/* Initialize global variables */
	DefaultGlobal();


	/*********************************************************************
	 *  Configuration for outputing log files
	 *********************************************************************/
	binout = fopen("binary_output.txt", "w");
	nbpout = fopen("abyss_output.txt", "w");
	//gpuout = fopen("cuda_output.txt", "w");
	fprintf(nbpout, "Abyss Output Starts!\n");
	//fprintf(binout, "Binary Output Starts!\n");
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


#ifdef CUDA
	int root_proc = 0;
	//if (AbyssProcessorNumber == ROOT)
	OpenDevice(&root_proc);
	cudaDeviceSynchronize(); 
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
	broadcastFromRoot(InitialNeighborRadius);
	broadcastFromRoot(FixNumNeighbor);

	/*
	// Insert this function definition at the top of your code after the include directives.
	char hostname[256];
	gethostname(hostname, sizeof(hostname));

	// Insert this code right after the  MPI initialization routines (though not a mandatory requirement 
	// to add there only). Please make a judgement based on your code.
	// Retrieve process ID and hostname
	pid_t pid = getpid();

	volatile int i = 0;
	while (0 == i)
	{
		std::cout << "My rank = " << AbyssProcessorNumber << " PID = " << pid << " running on Host = " << hostname << " in sleep " << std::endl;
		sleep(5);
	}
	*/

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

#ifdef PerformanceTrace
	//MPI_Reduce(&performance.IrregularForce, &IrrPerformance, 1, MPI_INT, MPI_SUM, 0, abyss_comm);
	MPI_Allreduce(MPI_IN_PLACE, &performance.IrregularForce, 1, MPI_LONG, MPI_SUM, abyss_comm);
	if (AbyssProcessorNumber == ROOT) {
		std::cerr << "Irr Time (ns) = " << performance.IrregularForce << std::endl;

		// Open the file in append mode
		std::ofstream outFile("performance", std::ios::app);

		// Check if the file opened successfully
		if (outFile.is_open()) {
			// Write the variable to the file
			outFile << performance.IrregularForce << std::endl;
			outFile << "Irregular Routine " << performance.IrregularRoutine << std::endl;
			outFile << "Regular Routine " << performance.RegularRoutine << std::endl;

			// Close the file
			outFile.close();
		} else {
			std::cerr << "Error opening file!" << std::endl;
		}
	}
#endif

	// Finalize the window and MPI environment
	//MPI_Win_free(&win);
	//MPI_Finalize();
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