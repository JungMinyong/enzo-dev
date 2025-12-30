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

typedef int MPI_Arg;

#ifdef NBODY
#include "abyss/particle.h"
#include "abyss/global.h"
#include "abyss/def.h"
#include "abyss/Queue.h"
#include "abyss/cuda/cuda_defs.h"
#endif

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "communication.h"

#if defined (NBODY) && defined (INDIVIDUALSTAR)
#include "NbodyRoutines.h"
#endif

MPI_Comm enzo_comm;
#ifdef NBODY
MPI_Win win;
MPI_Win win2;
MPI_Win win3;
MPI_Win win4;

Particle *particles;
GlobalVariable *global_variable;
int* Neighbors;
int* NewNeighbors;

MPI_Comm abyss_comm;
MPI_Comm inter_comm;
MPI_Comm local_comm;

#ifdef INDIVIDUALSTAR
MPI_Datatype MPI_ENZO_PTCL = MPI_DATATYPE_NULL;
MPI_Datatype MPI_ENZO_PTCL_SEND = MPI_DATATYPE_NULL;
MPI_Datatype MPI_ENZO_PTCL_RECV = MPI_DATATYPE_NULL;
#endif
#endif
int err; //just for debug 

/* function prototypes */
void my_exit(int exit_status);

#ifdef USE_MPI
void CommunicationErrorHandlerFn(MPI_Comm *comm, MPI_Arg *err, ...);
#ifdef NBODY
	int local_rank, local_size;
#endif
#endif


int CommunicationInitialize(Eint32 *argc, char **argv[])
{
 
#ifdef USE_MPI
 
  /* Initialize MPI and get info. */

  MPI_Arg world_rank;
  MPI_Arg world_size;
  MPI_Comm comm = MPI_COMM_WORLD;

  MPI_Init(argc, argv);
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

		int NumberOfCpuNodes = 0;
		int NumberOfGpuNodes = 0;
		int NumberOfEnzoProcessors = 0;
		char *env;

		env = std::getenv("ENZO_PROCESSOR_NUM");
		if (env == NULL) {
			fprintf(stderr, "Missing ENZO_PROCESSOR_NUM\n");
			exit(EXIT_FAILURE);
		}
		NumberOfEnzoProcessors = std::atoi(env);
		env = std::getenv("ENZO_CPU_NODE_NUM");
		if (env == NULL) {
			fprintf(stderr, "Missing ENZO_CPU_NODE_NUM\n");
			exit(EXIT_FAILURE);
		}
		NumberOfCpuNodes = std::atoi(env);
		env = std::getenv("ENZO_GPU_NODE_NUM");
		if (env == NULL) {
			fprintf(stderr, "Missing ENZO_GPU_NODE_NUM\n");
			exit(EXIT_FAILURE);
		}
		NumberOfGpuNodes = std::atoi(env);


		/***********************************
		 *     Communicator Setting        *
		 ***********************************/
		// by YS, create the group of processes in MPI_COMM_WORLD
		MPI_Group world_group;
		MPI_Comm_group(MPI_COMM_WORLD, &world_group);

		// Create a node-local communicator for shared-memory windows
		local_comm = MPI_COMM_NULL;
		MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, world_rank,
							MPI_INFO_NULL, &local_comm);

		local_rank = -1;
		local_size = 0;
		MPI_Comm_rank(local_comm, &local_rank);
		MPI_Comm_size(local_comm, &local_size);

		char raw_name[MPI_MAX_PROCESSOR_NAME];
		int raw_len = 0;
		MPI_Get_processor_name(raw_name, &raw_len);

		char short_name[MPI_MAX_PROCESSOR_NAME];
		int sn = 0;
		for (int i = 0; i < raw_len && sn < MPI_MAX_PROCESSOR_NAME - 1; ++i) {
			char c = raw_name[i];
			if (c == '.') break;
			if (c >= 'A' && c <= 'Z') c += 32;
			short_name[sn++] = c;
		}
		short_name[sn] = '\0';


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

		int total_nodes = NumberOfCpuNodes + NumberOfGpuNodes;
		int detected_nodes = static_cast<int>(unique_roots.size());
		if (detected_nodes != total_nodes) {
			if (world_rank == 0) {
				fprintf(stderr,
						"ENZO_CPU_NODE_NUM (%d) + ENZO_GPU_NODE_NUM (%d) must match detected nodes (%d)\n",
						NumberOfCpuNodes, NumberOfGpuNodes, detected_nodes);
			}
			exit(EXIT_FAILURE);
		}

		bool is_cpu_node = (node_id < NumberOfCpuNodes);
		bool is_gpu_node = (node_id >= NumberOfCpuNodes && node_id < total_nodes);

		if (local_rank == 0) {
			std::cout << "node " << node_id << " (" << short_name << ") classified as "
					  << (is_gpu_node ? "gpu" : "cpu") << std::endl;
		}

		std::cout << "world rank " << world_rank << " is on node id of " << node_id
				  << " (" << (is_gpu_node ? "gpu" : "cpu") << " node, local rank: " << local_rank
				  << ", total on node: " << local_size << ") " << std::endl;

		int cpu_rank_flag = is_cpu_node ? 1 : 0;
		int cpu_rank_count = 0;
		MPI_Allreduce(&cpu_rank_flag, &cpu_rank_count, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		std::vector<int> local_sizes(world_size, 0);
		MPI_Allgather(&local_size, 1, MPI_INT, local_sizes.data(), 1, MPI_INT, MPI_COMM_WORLD);

		std::vector<int> node_sizes(total_nodes, 0);
		for (int r = 0; r < world_size; ++r) {
			int root = node_roots[r];
			std::vector<int>::iterator it = std::lower_bound(unique_roots.begin(), unique_roots.end(), root);
			if (it == unique_roots.end() || *it != root) {
				continue;
			}
			int nid = std::distance(unique_roots.begin(), it);
			node_sizes[nid] = local_sizes[r];
		}

		int gpu_rank_count = 0;
		for (int nid = NumberOfCpuNodes; nid < total_nodes; ++nid) {
			gpu_rank_count += node_sizes[nid];
		}

		if (NumberOfEnzoProcessors < cpu_rank_count) {
			if (world_rank == 0) {
				fprintf(stderr,
						"ENZO_PROCESSOR_NUM (%d) is less than CPU-node ranks (%d); using %d\n",
						NumberOfEnzoProcessors, cpu_rank_count, cpu_rank_count);
			}
			NumberOfEnzoProcessors = cpu_rank_count;
		}
		if (NumberOfEnzoProcessors > TotalNumberOfProcessors) {
			if (world_rank == 0) {
				fprintf(stderr,
						"ENZO_PROCESSOR_NUM (%d) exceeds total ranks (%d); using %d\n",
						NumberOfEnzoProcessors, TotalNumberOfProcessors, TotalNumberOfProcessors);
			}
			NumberOfEnzoProcessors = TotalNumberOfProcessors;
		}

		int max_enzo = cpu_rank_count + gpu_rank_count;
		if (NumberOfEnzoProcessors > max_enzo) {
			if (world_rank == 0) {
				fprintf(stderr,
						"ENZO_PROCESSOR_NUM (%d) exceeds CPU+GPU ranks (%d); using %d\n",
						NumberOfEnzoProcessors, max_enzo, max_enzo);
			}
			NumberOfEnzoProcessors = max_enzo;
		}

		int remaining_enzo = NumberOfEnzoProcessors - cpu_rank_count;
		if (remaining_enzo < 0) {
			remaining_enzo = 0;
		}
		if (remaining_enzo > gpu_rank_count) {
			remaining_enzo = gpu_rank_count;
		}

		int gpu_order = -1;
		if (is_gpu_node) {
			int gpu_rank_offset = 0;
			for (int nid = NumberOfCpuNodes; nid < node_id; ++nid) {
				gpu_rank_offset += node_sizes[nid];
			}
			gpu_order = gpu_rank_offset + local_rank;
		}

		bool is_enzo_rank = is_cpu_node;
		bool is_abyss_rank = false;
		if (is_gpu_node) {
			if (gpu_order >= 0 && gpu_order < remaining_enzo) {
				is_enzo_rank = true;
			} else {
				is_abyss_rank = true;
			}
		}

		int local_has_abyss = 0;
		int abyss_flag = is_abyss_rank ? 1 : 0;
		MPI_Allreduce(&abyss_flag, &local_has_abyss, 1, MPI_INT, MPI_MAX, local_comm);

		NumberOfProcessors = NumberOfEnzoProcessors;
		NumberOfAbyssProcessors = TotalNumberOfProcessors - NumberOfEnzoProcessors;

		// Create a Enzo communicator
		if (is_enzo_rank)
		{
			MPI_Comm_split(MPI_COMM_WORLD, 1, world_rank, &enzo_comm);
		}
		else
		{
			MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, world_rank, &enzo_comm);
		}

		// Create a Abyss communicator
		if (is_abyss_rank)
		{
			MPI_Comm_split(MPI_COMM_WORLD, 1, world_rank, &abyss_comm);
		}
		else
		{
			MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, world_rank, &abyss_comm);
		}

		// Rank for each communicator
		if (MPI_COMM_NULL != enzo_comm)
		{
			MPI_Comm_rank(enzo_comm, &MyProcessorNumber);
			std::cout << "World rank " << world_rank << " (local rank " << local_rank << " of node " << node_id << ") is in enzo_comm"
					  << std::endl;
		}
		if (MPI_COMM_NULL != abyss_comm)
		{
			MPI_Comm_rank(abyss_comm, &AbyssProcessorNumber);
			std::cout << "World rank " << world_rank << " (local rank " << local_rank << " of node " << node_id << ") is in abyss_comm"
					  << std::endl;
		}

		// Create a Inter communicator
		int in_inter = 0;
#ifdef INDIVIDUAL
		if (MPI_COMM_NULL != enzo_comm) {
			in_inter = 1;
		}
		if (MPI_COMM_NULL != abyss_comm && AbyssProcessorNumber == 0) {
			in_inter = 1;
		}
#else
		if (MPI_COMM_NULL != enzo_comm && MyProcessorNumber == 0) {
			in_inter = 1;
		}
		if (MPI_COMM_NULL != abyss_comm && AbyssProcessorNumber == 0) {
			in_inter = 1;
		}
#endif

		if (in_inter)
		{
			MPI_Comm_split(MPI_COMM_WORLD, 1, world_rank, &inter_comm);
		}
		else
		{
			MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, world_rank, &inter_comm);
		}

		int inter_rank = -1;
		if (MPI_COMM_NULL != inter_comm)
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
#ifdef SEVN
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
#else
			int block_lengths[3] = {1, MAX_DIMENSION, 1};
			MPI_Aint displacements[3];
			MPI_Datatype types[3] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE};

			MPI_Aint base;
			MPI_Get_address(&dummy, &base);
			MPI_Get_address(&dummy.ID, &displacements[0]);
			MPI_Get_address(&dummy.BackgroundAcceleration, &displacements[1]);
			MPI_Get_address(&dummy.Mass, &displacements[2]);
			//MPI_Get_address(&dummy.isNew, &displacements[2]);

			for (int i = 0; i < 3; ++i)
					displacements[i] -= base;

			MPI_Type_create_struct(3, block_lengths, displacements, types, &MPI_ENZO_PTCL_SEND_RAW);
#endif
			MPI_Type_commit(&MPI_ENZO_PTCL_SEND_RAW);

			MPI_Aint lb=0, extent=sizeof(ParticleSendDataType);
			MPI_Type_create_resized(MPI_ENZO_PTCL_SEND_RAW, lb, extent, &MPI_ENZO_PTCL_SEND);
			MPI_Type_commit(&MPI_ENZO_PTCL_SEND);
			MPI_Type_free(&MPI_ENZO_PTCL_SEND_RAW);
		}

		{
			MPI_Datatype MPI_ENZO_PTCL_RECV_RAW;
			ParticleReceiveDataType dummy;
#ifdef SEVN
			int block_lengths[4] = {1, MAX_DIMENSION, MAX_DIMENSION, 5};
			MPI_Aint displacements[4];
			MPI_Datatype types[4] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};

			MPI_Aint base;
			MPI_Get_address(&dummy, &base);
			MPI_Get_address(&dummy.ID, &displacements[0]);
			MPI_Get_address(&dummy.Position, &displacements[1]);
			MPI_Get_address(&dummy.Velocity, &displacements[2]);
			MPI_Get_address(&dummy.InitialMass, &displacements[3]);

			for (int i = 0; i < 4; ++i)
				displacements[i] -= base;

			MPI_Type_create_struct(4, block_lengths, displacements, types, &MPI_ENZO_PTCL_RECV);
#else
			int block_lengths[3] = {1, MAX_DIMENSION, MAX_DIMENSION};
			MPI_Aint displacements[3];
			MPI_Datatype types[3] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE};

			MPI_Aint base;
			MPI_Get_address(&dummy, &base);
			MPI_Get_address(&dummy.ID, &displacements[0]);
			MPI_Get_address(&dummy.Position, &displacements[1]);
			MPI_Get_address(&dummy.Velocity, &displacements[2]);

			for (int i = 0; i < 3; ++i)
				displacements[i] -= base;

			MPI_Type_create_struct(3, block_lengths, displacements, types, &MPI_ENZO_PTCL_RECV_RAW);
 #endif
      
			MPI_Type_commit(&MPI_ENZO_PTCL_RECV_RAW);

			MPI_Aint lb=0, extent=sizeof(ParticleReceiveDataType);
			MPI_Type_create_resized(MPI_ENZO_PTCL_RECV_RAW, lb, extent, &MPI_ENZO_PTCL_RECV);
			MPI_Type_commit(&MPI_ENZO_PTCL_RECV);
			MPI_Type_free(&MPI_ENZO_PTCL_RECV_RAW);
		}
		#endif // INDIVIDUALSTAR


		/***********************************
		 *     Shared Memeory Setting      *
		 ***********************************/

		// Allocate shared memory on nodes that have Abyss ranks
		if (local_has_abyss)
		{
			if (local_rank == 0)
			{
				MPI_Win_allocate_shared(sizeof(Particle) * MaxNumParticle, sizeof(Particle), MPI_INFO_NULL, local_comm, &particles, &win);
				MPI_Win_allocate_shared(sizeof(GlobalVariable), sizeof(GlobalVariable), MPI_INFO_NULL, local_comm, &global_variable, &win2);
				MPI_Win_allocate_shared(sizeof(int) * MaxNumParticle * MaxNumNeighbor, sizeof(int), MPI_INFO_NULL, local_comm, &Neighbors, &win3);
				MPI_Win_allocate_shared(sizeof(int) * MaxNumParticle * MaxNumNeighbor, sizeof(int), MPI_INFO_NULL, local_comm, &NewNeighbors, &win4);
			}
			else
			{
				MPI_Win_allocate_shared(0, sizeof(Particle), MPI_INFO_NULL, local_comm, &particles, &win);
				MPI_Win_allocate_shared(0, sizeof(GlobalVariable), MPI_INFO_NULL, local_comm, &global_variable, &win2);
				MPI_Win_allocate_shared(0, sizeof(int), MPI_INFO_NULL, local_comm, &Neighbors, &win3);
				MPI_Win_allocate_shared(0, sizeof(int), MPI_INFO_NULL, local_comm, &NewNeighbors, &win4);
			}
			// Query shared memory of rank 0

			MPI_Aint size_bytes;
			int disp_unit;

			MPI_Win_shared_query(win, 0, &size_bytes, &disp_unit, &particles);
			MPI_Win_shared_query(win2, 0, &size_bytes, &disp_unit, &global_variable);
			MPI_Win_shared_query(win3, 0, &size_bytes, &disp_unit, &Neighbors);
			MPI_Win_shared_query(win4, 0, &size_bytes, &disp_unit, &NewNeighbors);
		}
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
  enzo_comm = MPI_COMM_WORLD;
  if (MyProcessorNumber == ROOT_PROCESSOR)
  {
	  printf("MPI_Init: NumberOfProcessors = %" ISYM "\n", NumberOfProcessors);
	  printf("MPI_Init: TotalNumberOfProcessors = %" ISYM "\n", TotalNumberOfProcessors);
  }
#endif // ABYSS

#else /* USE_MPI */
  MyProcessorNumber  = 0;
  NumberOfProcessors = 1;
 
#endif /* USE_MPI */
 
  CommunicationTime = 0;
 

  CommunicationDirection = 0; //COMMUNICATION_SEND_RECEIVE
 

  return SUCCESS;
  // return 1; // SUCCESS -> 1 by EW 2025.3.11
}

#ifdef USE_MPI
#ifdef NBODY
MPI_Datatype createQueueType() {
    MPI_Datatype QueueType;
    int blocklen[3] = {1, 1, 1};
    MPI_Datatype types[3] = {MPI_INT8_T, MPI_INT, MPI_DOUBLE};
	MPI_Aint disp[3], base;

	Queue sample;
    MPI_Get_address(&sample,           &base);
    MPI_Get_address(&sample.task,      &disp[0]);
    MPI_Get_address(&sample.pid,       &disp[1]);
    MPI_Get_address(&sample.next_time, &disp[2]);

	for (int i = 0; i < 3; ++i) disp[i] -= base;

    MPI_Type_create_struct(3, blocklen, disp, types, &QueueType);
    MPI_Type_commit(&QueueType);

    return QueueType;
}

MPI_Datatype createIparticleType() {
	MPI_Datatype IparticleType;
	int blocklen[8] = {1,1,1,1,1,1,1,1};
#ifdef CUDA_FLOAT
	MPI_Datatype types[8] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT,
							MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT};
#else
	MPI_Datatype types[8] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
							MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
#endif
	MPI_Aint disp[8], base;

	Iparticle sample;
	MPI_Get_address(&sample,		&base);
	MPI_Get_address(&sample.posx,	&disp[0]);
	MPI_Get_address(&sample.posy,	&disp[1]);
	MPI_Get_address(&sample.posz,	&disp[2]);
	MPI_Get_address(&sample.r2,		&disp[3]);
	MPI_Get_address(&sample.velx,	&disp[4]);
	MPI_Get_address(&sample.vely,	&disp[5]);
	MPI_Get_address(&sample.velz,	&disp[6]);
	MPI_Get_address(&sample.dtr,	&disp[7]);

	for (int i = 0; i < 8; ++i) disp[i] -= base;

	MPI_Type_create_struct(8, blocklen, disp, types, &IparticleType);
	MPI_Type_commit(&IparticleType);

	return IparticleType;
}

MPI_Datatype createJparticleType() {
	MPI_Datatype JparticleType;
	int blocklen[8] = {1,1,1,1,1,1,1,1};
#ifdef CUDA_FLOAT
	MPI_Datatype types[8] = {MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT,
							MPI_FLOAT, MPI_FLOAT, MPI_FLOAT, MPI_INT};
#else
	MPI_Datatype types[8] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,
							MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_LONG_LONG};
#endif
	MPI_Aint disp[8], base;

	Jparticle sample;
	MPI_Get_address(&sample,		&base);
	MPI_Get_address(&sample.posx,	&disp[0]);
	MPI_Get_address(&sample.posy,	&disp[1]);
	MPI_Get_address(&sample.posz,	&disp[2]);
	MPI_Get_address(&sample.mass,	&disp[3]);
	MPI_Get_address(&sample.velx,	&disp[4]);
	MPI_Get_address(&sample.vely,	&disp[5]);
	MPI_Get_address(&sample.velz,	&disp[6]);
	MPI_Get_address(&sample.index,	&disp[7]);

	for (int i = 0; i < 8; ++i) disp[i] -= base;

	MPI_Type_create_struct(8, blocklen, disp, types, &JparticleType);
	err = MPI_Type_commit(&JparticleType);
	if (err != MPI_SUCCESS) fprintf(stderr,"Error: MPI_Type_commit JparticleType failed!\n");

	return JparticleType;
}
#endif
#endif
 
#ifdef USE_MPI
void CommunicationErrorHandlerFn(MPI_Comm *comm, MPI_Arg *err, ...)
{
  char error_string[1024];
  MPI_Arg length, error_class;
  if (*err != MPI_ERR_OTHER) {
      MPI_Error_class(*err, &error_class);
      MPI_Error_string(error_class, error_string, &length);
      fprintf(stderr, "P%" ISYM": %s\n", MyProcessorNumber, error_string);
      MPI_Error_string(*err, error_string, &length);
      fprintf(stderr, "P%" ISYM": %s\n", MyProcessorNumber, error_string);
      ENZO_FAIL("MPI communication error.");
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
 
  return SUCCESS;
}

void CommunicationAbort(int status)
{

#ifdef USE_MPI
  MPI_Abort(enzo_comm,status);
#else
  //  my_exit(status);
#endif

  return;
}
