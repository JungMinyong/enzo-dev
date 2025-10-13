/***********************************************************************
/
/  COMMUNICATION ROUTINE: INITIALIZE
/
/  written by: Greg Bryan
/  date:       December, 1997
/  modified1: Yongseok Jo, 2025
/
/  PURPOSE:
/
************************************************************************/
 
#ifdef USE_MPI
#include "mpi.h"
#include <stdlib.h>
#endif /* USE_MPI */

#include <stdio.h>
#include <string.h>
#include <vector>
#include <algorithm>
#include <cstdlib> // For getenv

#define no_TEST_CUDA // (To be deleted)
#ifdef TEST_CUDA
#include <cuda.h> 
#include <cuda_runtime.h>
#endif

typedef int MPI_Arg;
//extern int local_rank, local_size;
//extern int MyProcessorNumber;
//extern int NumberOfProcessors;
//extern int TotalNumberOfProcessors;


//extern MPI_Errhandler CommunicationErrorHandler;
//extern float CommunicationTime;
//extern int CommunicationDirection;
//#include "ErrorExceptions.h"
//#include "macros_and_parameters.h"
//#include "typedefs.h"
//#include "global_data.h"

#ifdef NBODY
#undef NormalStar
#undef BlackHole
#undef max
#include "abyss/particle.h"
#include "abyss/global.h"
#include "abyss/def.h"
#undef NormalStar
#undef BlackHole
#endif

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "communication.h"
//#define ENZO_ONLY
#if defined (NBODY) && defined (INDIVIDUALSTAR)
#include "NbodyRoutines.h"
#endif

/*
extern int WorldProcessorNumber;
extern int AbyssProcessorNumber;
extern int NumberOfAbyssProcessors;
extern MPI_Datatype MPI_ENZO_PTCL;
extern MPI_Datatype MPI_ENZO_PTCL_SEND;
extern MPI_Datatype MPI_ENZO_PTCL_RECV;
extern MPI_Comm enzo_comm;
extern MPI_Comm abyss_comm;
extern MPI_Comm inter_comm;
extern MPI_Comm local_comm;
*/

#ifdef NBODY
Particle *particles_original;
Particle *particles;
int *ActiveIndexToOriginalIndex;
int *ActiveIndexToOriginalIndex_orginal;

MPI_Win win;
MPI_Win win2;
MPI_Win win3;

GlobalVariable *global_variable;
GlobalVariable *global_variable_original;

MPI_Comm abyss_comm;
MPI_Comm inter_comm;
MPI_Comm enzo_comm;
MPI_Comm local_comm;

#ifdef INDIVIDUALSTAR
MPI_Datatype MPI_ENZO_PTCL = MPI_DATATYPE_NULL;
MPI_Datatype MPI_ENZO_PTCL_SEND = MPI_DATATYPE_NULL;
MPI_Datatype MPI_ENZO_PTCL_RECV = MPI_DATATYPE_NULL;
#endif
#endif

 
/* function prototypes */
void my_exit(int exit_status);

#ifdef USE_MPI
void CommunicationErrorHandlerFn(MPI_Comm *comm, MPI_Arg *err, ...);
#ifdef NBODY
//#define NumberOfNbodyProcessors 1
	int local_rank, local_size;
	//int NumberOfAbyssProcessors;
	//int WorldProcessorNumber;
	//int AbyssProcessorNumber;
	//int TotalNumberOfProcessors;
#endif
#endif



int CommunicationInitialize(int &argc, char *argv[])
{



#ifdef USE_MPI
 
  /* Initialize MPI and get info. */

  MPI_Arg world_rank;
  MPI_Arg world_size;
  MPI_Comm comm = MPI_COMM_WORLD;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_create_errhandler(CommunicationErrorHandlerFn, &CommunicationErrorHandler);
  MPI_Comm_set_errhandler(comm, CommunicationErrorHandler);

#ifdef NBODY
	// by YS, N body function will be activated only if world_size > 2
	if (world_size >= 2) {

		std::cout << "Enzo-Abyss Communication Initialization Begins." << std::endl;

		TotalNumberOfProcessors = world_size;
		WorldProcessorNumber = world_rank;

		int NumberOfEnzoNodes = 0;
		int NumberOfEnzoProcessors = 0;
		char *env;

		env = std::getenv("ENZO_PROCESSOR_NUM");
		NumberOfEnzoProcessors = std::atoi(env);
		env = std::getenv("ENZO_NODE_NUM");
		NumberOfEnzoNodes = std::atoi(env);

		NumberOfProcessors = NumberOfEnzoProcessors;
		NumberOfAbyssProcessors = TotalNumberOfProcessors - NumberOfEnzoProcessors;


		/***********************************
		 *     Communicator Setting        *
		 ***********************************/
		// by YS, create the group of processes in MPI_COMM_WORLD
		MPI_Group world_group;
		MPI_Comm_group(MPI_COMM_WORLD, &world_group);

	#define MUST
	#ifdef MUST

		// 1) Get and normalize hostname (helps when some ranks see FQDNs, others short names)
		char raw_name[MPI_MAX_PROCESSOR_NAME]; int len = 0;
		MPI_Get_processor_name(raw_name, &len);

		// Make a normalized, short, lowercase name up to the first '.'
		char name[MPI_MAX_PROCESSOR_NAME];
		int n = 0;
		for (int i = 0; i < len && n < MPI_MAX_PROCESSOR_NAME-1; ++i) {
				char c = raw_name[i];
				if (c == '.') break;                    // drop domain suffixes
				if (c >= 'A' && c <= 'Z') c += 32;      // tolower ASCII
				name[n++] = c;
		}
		name[n] = '\0';

		// 2) Make a stable, non-negative color from hostname (FNV-1a masked to 31 bits)
		uint32_t h = 2166136261u;
		for (int i = 0; i < n; ++i) { h ^= (uint8_t)name[i]; h *= 16777619u; }
		int color = (int)(h & 0x7fffffff);          // ensure color >= 0

		// (Optional) Debug: ensure color is valid and show mapping
		#ifndef NDEBUG
		if (color < 0) { fprintf(stderr, "[%d] bad color=%d\n", world_rank, color); MPI_Abort(MPI_COMM_WORLD, 1); }
		fprintf(stderr, "[%d] host=%s color=%d\n", world_rank, name, color);
		#endif

		// 3) Split by color; 'key' controls the ordering inside each node communicator
		MPI_Comm local_comm = MPI_COMM_NULL;
		MPI_Comm_split(MPI_COMM_WORLD, color, /*key=*/world_rank, &local_comm);

		// 4) Now you have contiguous local ranks 0..local_size-1 per node
		int local_rank = -1, local_size = 0;
		MPI_Comm_rank(local_comm, &local_rank);
		MPI_Comm_size(local_comm, &local_size);
#else
		// Create a shared memory communicator to identify node-local processes
		//MPI_Comm local_comm = MPI_COMM_NULL;
		//MPI_Comm_split(MPI_COMM_WORLD, (int)color, world_rank, &local_comm);
		MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, world_rank, MPI_INFO_NULL, &local_comm);
		int rc = MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, world_rank, MPI_INFO_NULL, &local_comm);
		assert(rc == MPI_SUCCESS && local_comm != MPI_COMM_NULL);

		// Get node-local rank
		MPI_Comm_rank(local_comm, &local_rank);
		MPI_Comm_size(local_comm, &local_size);
#endif

		int local_root = world_rank;
		MPI_Allreduce(&world_rank, &local_root, 1, MPI_INT, MPI_MIN, local_comm);

		// Gather all local roots to determine node ID
		std::vector<int> node_roots(world_size, -1);
		MPI_Allgather(&local_root, 1, MPI_INT, node_roots.data(), 1, MPI_INT, MPI_COMM_WORLD);

		// Assign a unique node_id based on the unique sorted node roots
		std::vector<int> unique_roots = node_roots;
		std::sort(unique_roots.begin(), unique_roots.end());
		unique_roots.erase(std::unique(unique_roots.begin(), unique_roots.end()), unique_roots.end());

		int node_id = std::distance(unique_roots.begin(),
									std::find(unique_roots.begin(), unique_roots.end(), local_root));
		int ordered_rank =  node_id * local_size + local_rank;
		std::cout << "world rank " << world_rank << " is on node id of " << node_id
				  << " (local rank: " << local_rank << ", total on node: " << local_size << ")  "
				  << "ordered rank: " << ordered_rank << " " << std::endl;


		// This is for inter_comm
		/*
		int ranks_inter[2];
		if (ordered_rank == 0) {
			ranks_inter[0] = world_rank;
		}
		if (ordered_rank == NumberOfEnzoProcessors) {
			ranks_inter[1] = world_rank;
		}*/

		// Create a Enzo communicator 
		if (ordered_rank < NumberOfEnzoProcessors)
		{
			MPI_Comm_split(MPI_COMM_WORLD, 1, world_rank, &enzo_comm);
		}
		else
		{
			MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, world_rank, &enzo_comm);
		}

		// Create a Abyss communicator
		if (ordered_rank >= NumberOfEnzoProcessors)
		{
			MPI_Comm_split(MPI_COMM_WORLD, 1, world_rank, &abyss_comm);
		}
		else
		{
			MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, world_rank, &abyss_comm);
		}

		// Create a Inter communicator
#ifdef INDIVIDUAL
		if (ordered_rank <= NumberOfEnzoProcessors)
#else
		if (ordered_rank == 0 || ordered_rank == NumberOfEnzoProcessors)
#endif
		{
			MPI_Comm_split(MPI_COMM_WORLD, 1, world_rank, &inter_comm);
		}
		else
		{
			MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, world_rank, &inter_comm);
		}




		// Rank for each communicator 
		if (ordered_rank < NumberOfEnzoProcessors)
		{
			MPI_Comm_rank(enzo_comm, &MyProcessorNumber);
			std::cout << "World rank " << world_rank << " (local rank " << local_rank << " of node " << node_id << ") is in enzo_comm"
					  << std::endl;
		}
		if (ordered_rank >= NumberOfEnzoProcessors)
		{
			MPI_Comm_rank(abyss_comm, &AbyssProcessorNumber);
			std::cout << "World rank " << world_rank << " (local rank " << local_rank << " of node " << node_id << ") is in abyss_comm"
					  << std::endl;
		}
		int inter_rank;
#ifdef INDIVIDUAL
		if (ordered_rank <= NumberOfEnzoProcessors)
#else
		if (ordered_rank == 0 || ordered_rank == NumberOfEnzoProcessors)
#endif
		{
			MPI_Comm_rank(inter_comm, &inter_rank);
			std::cout << "World rank " << world_rank << " (local rank " << local_rank << " of node " << node_id << ") is in inter_comm"
					  << std::endl;
		}


		if (MPI_COMM_NULL != abyss_comm) {
			fprintf(stderr,"nbody: (%d, %d)\n", world_rank, AbyssProcessorNumber);
		}
		if (MPI_COMM_NULL != enzo_comm) {
			fprintf(stderr,"enzo: (%d, %d)\n", world_rank, MyProcessorNumber);
		}
		if (MPI_COMM_NULL != inter_comm) {
			fprintf(stderr,"inter: (%d, %d)\n", world_rank, inter_rank);
		}

#define COMM_TEST // (To be deleted)
#ifdef COMM_TEST
 		if (MyProcessorNumber == ROOT_PROCESSOR) {
			// get its Fortran handle too
			MPI_Fint id_enzo = MPI_Comm_c2f(enzo_comm);
			MPI_Fint id_inter = MPI_Comm_c2f(inter_comm);
			MPI_Fint id_abyss = MPI_Comm_c2f(abyss_comm);
			fprintf(stderr, "COMM ID (%d): enzo_comm=%d, inter_comm=%d, abyss_comm=%d\n", world_rank, (int) id_enzo, (int) id_inter, (int) id_abyss);
		}
#endif

#ifdef TEST_CUDA // (To be deleted)
	int deviceCount, devid;
	cudaGetDeviceCount(&deviceCount);
	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, devid);
	//  char *hostname = getenv("HOSTNAME");
	fprintf(stderr, "# GPU initialization - rank: %d; NGPU %d; device: %d %s\n", world_rank, deviceCount, devid, prop.name);
#endif

#ifdef INDIVIDUALSTAR
		/***********************************
		 *     Struct MPI Data Type        *
		 ***********************************/
		{
			MPI_Datatype MPI_ENZO_PTCL_RAW;

			ParticleDataType dummy;

			int block_lengths[5] = {
				1,                      // ID
				MAX_DIMENSION,          // Position
				MAX_DIMENSION,          // Velocity
				MAX_DIMENSION,          // BackgroundAcceleration
				4                       // Mass, CreationTime, DynamicalTime, Metallicity
			};

			MPI_Aint displacements[5];
			MPI_Datatype types[5] = {
				MPI_INT,
				MPI_DOUBLE,
				MPI_DOUBLE,
				MPI_DOUBLE,
				MPI_DOUBLE
			};

			MPI_Aint base;
			MPI_Get_address(&dummy, &base);
			MPI_Get_address(&dummy.ID, &displacements[0]);
			MPI_Get_address(&dummy.Position, &displacements[1]);
			MPI_Get_address(&dummy.Velocity, &displacements[2]);
			MPI_Get_address(&dummy.BackgroundAcceleration, &displacements[3]);
			MPI_Get_address(&dummy.Mass, &displacements[4]);

			for (int i = 0; i < 5; ++i)
				displacements[i] -= base;

			MPI_Type_create_struct(5, block_lengths, displacements, types, &MPI_ENZO_PTCL_RAW);
			MPI_Type_commit(&MPI_ENZO_PTCL_RAW);

			MPI_Aint lb=0, extent=sizeof(ParticleDataType);
			//MPI_Type_get_extent(MPI_ENZO_PTCL, &lb, &extent);
			MPI_Type_create_resized(MPI_ENZO_PTCL_RAW, lb, extent, &MPI_ENZO_PTCL);
			MPI_Type_commit(&MPI_ENZO_PTCL);
			MPI_Type_free(&MPI_ENZO_PTCL_RAW);
			fprintf(stderr,"Extent = %ld, sizeof = %zu\n", (long)extent, sizeof(ParticleDataType));
		}

		{
			MPI_Datatype MPI_ENZO_PTCL_SEND_RAW;
			ParticleSendDataType dummy;

			int block_lengths[2] = {1, MAX_DIMENSION};
			MPI_Aint displacements[2];
			MPI_Datatype types[2] = {MPI_INT, MPI_DOUBLE};

			MPI_Aint base;
			MPI_Get_address(&dummy, &base);
			MPI_Get_address(&dummy.ID, &displacements[0]);
			MPI_Get_address(&dummy.BackgroundAcceleration, &displacements[1]);
			//MPI_Get_address(&dummy.isNew, &displacements[2]);

			for (int i = 0; i < 2; ++i)
				displacements[i] -= base;

			MPI_Type_create_struct(2, block_lengths, displacements, types, &MPI_ENZO_PTCL_SEND_RAW);
			MPI_Type_commit(&MPI_ENZO_PTCL_SEND_RAW);

			MPI_Aint lb=0, extent=sizeof(ParticleSendDataType);
			MPI_Type_create_resized(MPI_ENZO_PTCL_SEND_RAW, lb, extent, &MPI_ENZO_PTCL_SEND);
			MPI_Type_commit(&MPI_ENZO_PTCL_SEND);
			MPI_Type_free(&MPI_ENZO_PTCL_SEND_RAW);
		}

		{
			MPI_Datatype MPI_ENZO_PTCL_RECV_RAW;
			ParticleReceiveDataType dummy;

			int block_lengths[3] = {1, MAX_DIMENSION, MAX_DIMENSION};
			MPI_Aint displacements[3];
			MPI_Datatype types[3] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE};

			MPI_Aint base;
			MPI_Get_address(&dummy, &base);
			MPI_Get_address(&dummy.ID, &displacements[0]);
			MPI_Get_address(&dummy.Position, &displacements[1]);
			MPI_Get_address(&dummy.Velocity, &displacements[2]);

			displacements[0] -= base;
			displacements[1] -= base;
			displacements[2] -= base;

			MPI_Type_create_struct(3, block_lengths, displacements, types, &MPI_ENZO_PTCL_RECV_RAW);
			MPI_Type_commit(&MPI_ENZO_PTCL_RECV_RAW);

			MPI_Aint lb=0, extent=sizeof(ParticleReceiveDataType);
			MPI_Type_create_resized(MPI_ENZO_PTCL_RECV_RAW, lb, extent, &MPI_ENZO_PTCL_RECV);
			MPI_Type_commit(&MPI_ENZO_PTCL_RECV);
			MPI_Type_free(&MPI_ENZO_PTCL_RECV_RAW);
		}
#endif

		/***********************************
		 *     Shared Memeory Setting      *
		 ***********************************/

		// Allocate shared memory
		//if (abyss_comm != MPI_COMM_NULL)
		//{
		if (AbyssProcessorNumber == 0)
		{
			// MPI_Win_allocate_shared(sizeof(int), sizeof(int), MPI_INFO_NULL, local_comm, &shared_mem, &win);
			MPI_Win_allocate_shared(sizeof(Particle) * MaxNumberOfParticle, sizeof(Particle),
									MPI_INFO_NULL, local_comm, &particles_original, &win);
			MPI_Win_allocate_shared(sizeof(GlobalVariable), sizeof(GlobalVariable),
									MPI_INFO_NULL, local_comm, &global_variable_original, &win2);
			MPI_Win_allocate_shared(sizeof(int) * MaxNumberOfParticle, sizeof(int),
									MPI_INFO_NULL, local_comm, &ActiveIndexToOriginalIndex_orginal, &win3);
		}
		else
		{
			MPI_Win_allocate_shared(0, sizeof(Particle), MPI_INFO_NULL, local_comm, &particles_original, &win);
			MPI_Win_allocate_shared(0, sizeof(GlobalVariable), MPI_INFO_NULL, local_comm, &global_variable_original, &win2);
			MPI_Win_allocate_shared(0, sizeof(int), MPI_INFO_NULL, local_comm, &ActiveIndexToOriginalIndex_orginal, &win3);
		}
		// Query shared memory of rank 0

		MPI_Aint size_bytes;
		int disp_unit;

		MPI_Win_shared_query(win, 0, &size_bytes, &disp_unit, &particles);
		MPI_Win_shared_query(win2, 0, &size_bytes, &disp_unit, &global_variable);
		MPI_Win_shared_query(win3, 0, &size_bytes, &disp_unit, &ActiveIndexToOriginalIndex);
		//}
	}
	else
	{
		std::cout << "The number of cores ought to be at least two for Enzo-Abyss!!!" << std::endl;
		std::cerr << "The number of cores ought to be at least two for Enzo-Abyss!!!" << std::endl;
		exit(EXIT_FAILURE);
	}
#else //  ABYSS
  NumberOfProcessors = world_size;
  MyProcessorNumber = world_rank;
  if (MyProcessorNumber == ROOT_PROCESSOR)
  {
	  printf("MPI_Init: NumberOfProcessors = %" ISYM "\n", NumberOfProcessors);
	  printf("MPI_Init: TotalNumberOfProcessors = %" ISYM "\n", TotalNumberOfProcessors);
  }
#endif // ABYSS
#else /* USE_MPI */
 
  //MyProcessorNumber  = 0;
  //NumberOfProcessors = 1;
 
#endif /* USE_MPI */
 
  CommunicationTime = 0;


  CommunicationDirection = 0; //COMMUNICATION_SEND_RECEIVE
 

//   return SUCCESS;
  return 1; // SUCCESS -> 1 by EW 2025.3.11
  return 1; //SUCCESS;
}
 
#ifdef USE_MPI
void CommunicationErrorHandlerFn(MPI_Comm *comm, MPI_Arg *err, ...)
{
  char error_string[1024];
  MPI_Arg length, error_class;
  if (*err != MPI_ERR_OTHER) {
      MPI_Error_class(*err, &error_class);
      MPI_Error_string(error_class, error_string, &length);
      //fprintf(stderr, "P%"ISYM": %s\n", MyProcessorNumber, error_string);
      MPI_Error_string(*err, error_string, &length);
      //fprintf(stderr, "P%"ISYM": %s\n", MyProcessorNumber, error_string);
	  throw(EnzoFatalException("MPI communication error.", __FILE__, __LINE__));
  } // ENDIF MPI_ERROR
  return;
}
#endif /* USE_MPI */
 
int CommunicationFinalize()
{
 
#ifdef USE_MPI
#ifdef NBODY
	MPI_Comm_free(&enzo_comm);
	MPI_Comm_free(&abyss_comm);
	MPI_Comm_free(&inter_comm);
	MPI_Type_free(&MPI_ENZO_PTCL);
	MPI_Type_free(&MPI_ENZO_PTCL_RECV);
	MPI_Type_free(&MPI_ENZO_PTCL_SEND);
#endif
  MPI_Errhandler_free(&CommunicationErrorHandler);
  MPI_Finalize();
#endif /* USE_MPI */
 
//   return SUCCESS;
  return 1; // SUCCESS -> 1 by EW 2025.3.11
  return 1; //SUCCESS;
}

void CommunicationAbort(int status)
{

#ifdef USE_MPI
  //MPI_Abort(MPI_COMM_WORLD,status);
  MPI_Abort(enzo_comm,status);
#else
  //  my_exit(status);
#endif

  return;
}

