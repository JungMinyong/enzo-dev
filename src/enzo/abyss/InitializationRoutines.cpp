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

void InitialAssignmentOfTasks(std::vector<int>& data, double next_time, int NumTask, int TAG);
void InitialAssignmentOfTasks(std::vector<int>& data, int NumTask, int TAG);
void InitialAssignmentOfTasks(int data, int NumTask, int TAG);
void InitialAssignmentOfTasks(int* data, int NumTask, int TAG);
void broadcastFromRoot(double &data);
void broadcastFromRoot(ULL &data);
void broadcastFromRoot(int &data);

/* Initialization */
void InitializationRoutines(QueueScheduler &queue_scheduler, Worker *workers)
{

    Particle* ptcl;
	int min_time_level=0;
	TaskName task;
	int total_tasks;
	MPI_Request request;  // Pointer to the request handle
	MPI_Status status;    // Pointer to the status object
	int completed_tasks=0;

    std::vector<int> ParticleIndices;
    ParticleIndices.reserve(NumberOfParticle);
    ParticleIndices.resize(NumberOfParticle);

    for (int i = 0; i <= LastParticleIndex; i++)
    {
        ParticleIndices[i] = i;
    }

    std::cout << "Initialization of particles starts." << std::endl;
    queue_scheduler.initialize(InitAcc1);
    queue_scheduler.takeQueue(ParticleIndices);
    do
    {
        // queue_scheduler.printFreeWorker();
        // queue_scheduler.printWorkerToGo();
        queue_scheduler.assignQueueAuto();
        // queue_scheduler.printFreeWorker();
        // queue_scheduler.printWorkerToGo();
        queue_scheduler.runQueueAuto();
        queue_scheduler.waitQueue(0); // blocking wait
    } while (queue_scheduler.isComplete());

    std::cout << "Init 01 done" << std::endl;
    queue_scheduler.initialize(InitAcc2);
    queue_scheduler.takeQueue(ParticleIndices);
    do
    {
        queue_scheduler.assignQueueAuto();
        queue_scheduler.runQueueAuto();
        queue_scheduler.waitQueue(0); // blocking wait
    } while (queue_scheduler.isComplete());

    std::cout << "Init 02 done" << std::endl;

#ifdef FEWBODY
    // Primordial binary search

    queue_scheduler.initialize(SearchPrimordialGroup);
    queue_scheduler.takeQueue(ParticleIndices);
    do
    {
        queue_scheduler.assignQueueAuto();
        queue_scheduler.runQueueAuto();
        queue_scheduler.waitQueue(0); // blocking wait
    } while (queue_scheduler.isComplete());

    std::cout << "Primordial binary search done" << std::endl;

    // example code by EW 2025.1.7
    Queue queue;
    int rank;
    int OriginalLastParticleIndex = LastParticleIndex;
    formPrimordialBinaries(OriginalLastParticleIndex);
    assert(OriginalLastParticleIndex <= LastParticleIndex); // for debugging by EW 2025.1.4
    assert(CMPtclWorker.empty());                           // for debugging by EW 2025.1.4
    if (OriginalLastParticleIndex != LastParticleIndex)
    {
        std::cout << "In total, " << LastParticleIndex - OriginalLastParticleIndex
                  << " primordial binaries are created." << std::endl;
        queue_scheduler.initialize(MakePrimordialGroup);
        for (int i = OriginalLastParticleIndex + 1; i <= LastParticleIndex; i++)
        {
            std::cout << "New Primordial Binary of PID="
                      << i << " is created with being assigned to a worker of rank "
                      << rank << "." << std::endl;
            ptcl = &particles[i];
            CMPtclWorker.insert({ptcl->ParticleIndex, CMPtclWorker.size() % NumberOfWorker + 1});
            ParticleIndices.push_back(ptcl->ParticleIndex);
            rank = CMPtclWorker[ptcl->ParticleIndex];

            queue.task = MakePrimordialGroup;
            queue.pid = ptcl->ParticleIndex;
            workers[rank].addQueue(queue);
            queue_scheduler.assignWorker(&workers[rank]);
        }
        queue_scheduler.setTotalQueue(CMPtclWorker.size());
        do
        {
            queue_scheduler.runQueueAuto();
            queue_scheduler.waitQueue(0);
        } while (queue_scheduler.isComplete());
    }
    else
    {
        std::cout << "There is no primordial binary." << std::endl;
    }
    fprintf(stdout, "PrimordialBinariesRoutine has ended...\n"
                    "The total number of particles is  %d\n",
            NumberOfParticle);
    fflush(stdout);
#endif

    // Initialize Time Step
    queue_scheduler.initialize(InitTime);
    queue_scheduler.takeQueue(ParticleIndices);
    do
    {
        queue_scheduler.assignQueueAuto();
        queue_scheduler.runQueueAuto();
        queue_scheduler.waitQueue(0); // blocking wait
    } while (queue_scheduler.isComplete());

    /*
    for (int i=0; i<=LastParticleIndex; i++) {
        ptcl = &particles[i];
        if (ptcl->isActive)
            fprintf(stdout, "PID=%d, CurrentTime (Irr, Reg) = (%.3e(%llu), %.3e(%llu)) Myr\n"
                            "dtIrr = %.4e Myr, dtReg = %.4e Myr, blockIrr=%llu (%d), blockReg=%llu (%d)\n"
                            "NumNeighbor= %d\n",
                    ptcl->PID,
                    ptcl->CurrentTimeIrr * EnzoTimeStep * 1e10 / 1e6,
                    ptcl->CurrentBlockIrr,
                    ptcl->CurrentTimeReg * EnzoTimeStep * 1e10 / 1e6,
                    ptcl->CurrentBlockReg,
                    // NextRegTimeBlock*time_step*EnzoTimeStep*1e10/1e6,
                    // NextRegTimeBlock,
                    ptcl->TimeStepIrr * EnzoTimeStep * 1e10 / 1e6,
                    ptcl->TimeStepReg * EnzoTimeStep * 1e10 / 1e6,
                    ptcl->TimeBlockIrr,
                    ptcl->TimeLevelIrr,
                    ptcl->TimeBlockReg,
                    ptcl->TimeLevelReg,
                    ptcl->NumberOfNeighbor);
    }
    */

    /* synchronization */
    // ParticleSynchronization();

    /* timestep correction */
    {
        std::cout << "Time Step correction." << std::endl;
        for (int i = 0; i <= LastParticleIndex; i++)
        {
            ptcl = &particles[i];

            if (!ptcl->isActive)
                continue;

            if (ptcl->NumberOfNeighbor != 0)
            {
                while (ptcl->TimeLevelIrr >= ptcl->TimeLevelReg)
                {
                    ptcl->TimeStepIrr *= 0.5;
                    ptcl->TimeBlockIrr *= 0.5;
                    ptcl->TimeLevelIrr--;
                }
            }
            if (ptcl->TimeLevelIrr < min_time_level)
            {
                min_time_level = ptcl->TimeLevelIrr;
            }
        }

        // resetting time_block based on the system
        time_block = std::max(-60, min_time_level - MIN_LEVEL_BUFFER);
        block_max = static_cast<ULL>(pow(2, -time_block));
        time_step = pow(2, time_block);

        for (int i = 0; i <= LastParticleIndex; i++)
        {
            ptcl = &particles[i];

            if (!ptcl->isActive)
                continue;

            ptcl->TimeBlockIrr = static_cast<ULL>(pow(2, ptcl->TimeLevelIrr - time_block));
            ptcl->TimeBlockReg = static_cast<ULL>(pow(2, ptcl->TimeLevelReg - time_block));
#ifdef IRR_TEST
            ptcl->TimeStepReg = 1;
            ptcl->TimeLevelReg = 0;
            ptcl->TimeBlockReg = block_max;
#endif
            ptcl->NextBlockIrr = ptcl->CurrentBlockIrr + ptcl->TimeBlockIrr; // of this particle
        }
        std::cout << "Time Step done." << std::endl;
    }

    /*
        for (int i=0; i<=LastParticleIndex; i++) {
            ptcl = &particles[i];
            if (ptcl->isActive)
                fprintf(stdout, "PID=%d, CurrentTime (Irr, Reg) = (%.3e(%llu), %.3e(%llu)) Myr\n"
                                "dtIrr = %.4e Myr, dtReg = %.4e Myr, blockIrr=%llu (%d), blockReg=%llu (%d)\n"
                                "NumNeighbor= %d\n",
                        ptcl->PID,
                        ptcl->CurrentTimeIrr * EnzoTimeStep * 1e10 / 1e6,
                        ptcl->CurrentBlockIrr,
                        ptcl->CurrentTimeReg * EnzoTimeStep * 1e10 / 1e6,
                        ptcl->CurrentBlockReg,
                        // NextRegTimeBlock*time_step*EnzoTimeStep*1e10/1e6,
                        // NextRegTimeBlock,
                        ptcl->TimeStepIrr * EnzoTimeStep * 1e10 / 1e6,
                        ptcl->TimeStepReg * EnzoTimeStep * 1e10 / 1e6,
                        ptcl->TimeBlockIrr,
                        ptcl->TimeLevelIrr,
                        ptcl->TimeBlockReg,
                        ptcl->TimeLevelReg,
                        ptcl->NumberOfNeighbor);
        }
        */

    /* timestep variable synchronization */
    {
        std::cout << "Time Step synchronization." << std::endl;
        task = TimeSync;
        completed_tasks = 0;
        total_tasks = NumberOfWorker;
        InitialAssignmentOfTasks(task, NumberOfWorker, TASK_TAG);
        // MPI_Waitall(NumberOfCommunication, requests, statuses);
        // NumberOfCommunication = 0;
        broadcastFromRoot(time_block);
        broadcastFromRoot(block_max);
        broadcastFromRoot(time_step);
        // MPI_Win_sync(win);  // Synchronize memory
        // MPI_Barrier(shared_comm);
        while (completed_tasks < total_tasks)
        {
            MPI_Irecv(&task, 1, MPI_INT, MPI_ANY_SOURCE, TERMINATE_TAG, abyss_comm, &request);
            MPI_Wait(&request, &status);
            completed_tasks++;
        }
        fprintf(nbpout, "nbody+:time_block = %d, EnzoTimeStep=%e\n", time_block, EnzoTimeStep);
        //fflush(stderr);
    }

    /* Particle Initialization Check */
    {
        //, NextRegTime= %.3e Myr(%llu),
        for (int i=0; i<=LastParticleIndex; i++) {
            ptcl = &particles[i];
            fprintf(nbpout, "PID=%d,PI=%d, (%d)=",ptcl->PID,ptcl->ParticleIndex,ptcl->NumberOfNeighbor);
            for (int j=0;j<ptcl->NumberOfNeighbor;j++) {
                fprintf(nbpout, "%d, ",ptcl->Neighbors[j]);
            }
            fprintf(nbpout, "\n");
        }
    }
        for (int i=0; i<=LastParticleIndex; i++) {
            ptcl = &particles[i];
            fprintf(nbpout, "PID=%d(%d), CurrentTime (Irr, Reg) = (%.3e(%llu), %.3e(%llu)) Myr\n"\
                    "dtIrr = %.4e Myr, dtReg = %.4e Myr, blockIrr=%llu (%d), blockReg=%llu (%d)\n"\
                    "NumNeighbor= %d\n",
                    ptcl->PID,
                    ptcl->ParticleIndex,
                    ptcl->CurrentTimeIrr*EnzoTimeStep*1e10/1e6,
                    ptcl->CurrentBlockIrr,
                    ptcl->CurrentTimeReg*EnzoTimeStep*1e10/1e6,
                    ptcl->CurrentBlockReg,
                    //NextRegTimeBlock*time_step*EnzoTimeStep*1e10/1e6,
                    //NextRegTimeBlock,
                    ptcl->TimeStepIrr*EnzoTimeStep*1e10/1e6,
                    ptcl->TimeStepReg*EnzoTimeStep*1e10/1e6,
                    ptcl->TimeBlockIrr,
                    ptcl->TimeLevelIrr,
                    ptcl->TimeBlockReg,
                    ptcl->TimeLevelReg,
                    ptcl->NumberOfNeighbor
                    );

            fprintf(nbpout, " a_tot = (%.4e,%.4e,%.4e), a_reg = (%.4e,%.4e,%.4e), a_irr = (%.4e,%.4e,%.4e), n_n=%d, R=%.3e\n\
                    a1_reg = (%.4e,%.4e,%.4e), a2_reg = (%.4e,%.4e,%.4e), a3_reg = (%.4e,%.4e,%.4e)\n\
                    a1_irr = (%.4e,%.4e,%.4e), a2_irr = (%.4e,%.4e,%.4e), a3_irr = (%.4e,%.4e,%.4e)\n",
                    ptcl->a_tot[0][0],
                    ptcl->a_tot[1][0],
                    ptcl->a_tot[2][0],
                    ptcl->a_reg[0][0],
                    ptcl->a_reg[1][0],
                    ptcl->a_reg[2][0],
                    ptcl->a_irr[0][0],
                    ptcl->a_irr[1][0],
                    ptcl->a_irr[2][0],
                    ptcl->NumberOfNeighbor,
                    ptcl->RadiusOfNeighbor,
                    ptcl->a_reg[0][1],
                    ptcl->a_reg[1][1],
                    ptcl->a_reg[2][1],
                    ptcl->a_reg[0][2],
                    ptcl->a_reg[1][2],
                    ptcl->a_reg[2][2],
                    ptcl->a_reg[0][3],
                    ptcl->a_reg[1][3],
                    ptcl->a_reg[2][3],
                    ptcl->a_irr[0][1],
                    ptcl->a_irr[1][1],
                    ptcl->a_irr[2][1],
                    ptcl->a_irr[0][2],
                    ptcl->a_irr[1][2],
                    ptcl->a_irr[2][2],
                    ptcl->a_irr[0][3],
                    ptcl->a_irr[1][3],
                    ptcl->a_irr[2][3]
                        );

        }
        fflush(nbpout);
    /* Particle Initialization Check */
    // /*
        //, NextRegTime= %.3e Myr(%llu),
        /*
    {
        for (int i = 0; i <= LastParticleIndex; i++)
        {
            ptcl = &particles[i];
            if (ptcl->isActive)
                fprintf(nbpout, "PID=%d, CurrentTime (Irr, Reg) = (%.3e(%llu), %.3e(%llu)) Myr\n"
                                "dtIrr = %.4e Myr, dtReg = %.4e Myr, blockIrr=%llu (%d), blockReg=%llu (%d)\n"
                                "NumNeighbor= %d\n",
                        ptcl->PID,
                        ptcl->CurrentTimeIrr * EnzoTimeStep * 1e10 / 1e6,
                        ptcl->CurrentBlockIrr,
                        ptcl->CurrentTimeReg * EnzoTimeStep * 1e10 / 1e6,
                        ptcl->CurrentBlockReg,
                        // NextRegTimeBlock*time_step*EnzoTimeStep*1e10/1e6,
                        // NextRegTimeBlock,
                        ptcl->TimeStepIrr * EnzoTimeStep * 1e10 / 1e6,
                        ptcl->TimeStepReg * EnzoTimeStep * 1e10 / 1e6,
                        ptcl->TimeBlockIrr,
                        ptcl->TimeLevelIrr,
                        ptcl->TimeBlockReg,
                        ptcl->TimeLevelReg,
                        ptcl->NumberOfNeighbor);
        }
        fflush(nbpout);
    }
        */
}





void InitializationAfterCommunication() {

    /* Initialize New Particle */ 
    /*  Neighbor inclusion might be needed (to be updated) */


    /* Initialize Particle Attributes */
    /* Since we're not doing full-initialization, we have to do more work on time steps
    e.g., if enzo time can be smaller than regualr time steps. we gotta re-normalize it.
    but this part is not complete yet. */
    Particle *ptcl;
    for (int i = 0; i <= LastParticleIndex; i++) {
        ptcl = &particles[i];
        ptcl->CurrentTimeIrr = 0.;
        ptcl->CurrentBlockIrr = 0;
        ptcl->CurrentTimeReg = 0.;
        ptcl->CurrentBlockReg = 0;
        ptcl->NewCurrentBlockIrr = 0;
        ptcl->NextBlockIrr = ptcl->CurrentBlockIrr + ptcl->TimeBlockIrr; // of this particle
        if (ptcl->Position[0] != ptcl->Position[0])
        {
            fprintf(stderr, "%d, %e, %e\n", ptcl->PID, ptcl->Position[0]);
            throw std::runtime_error("InitializationAfterCommunication\n");
        }
    }
}