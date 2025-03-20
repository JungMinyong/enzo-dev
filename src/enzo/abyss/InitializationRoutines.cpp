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
void initializeTime(QueueScheduler &queue_scheduler, Worker *workers, std::vector<int> ParticleIndices);

void formPrimordialBinaries(int beforeLastParticleIndex);

void InitializationOnGPU(QueueScheduler &queue_scheduler, Worker *workers);

/* Initialization */
void InitializationRoutines(QueueScheduler &queue_scheduler, Worker *workers)
{

    Particle* ptcl;
	TaskName task;
	int total_tasks;
	MPI_Request request;  // Pointer to the request handle
	MPI_Status status;    // Pointer to the status object
	int completed_tasks=0;

    std::vector<int> ParticleIndices;
    ParticleIndices.reserve(NumberOfParticle);
    ParticleIndices.resize(NumberOfParticle);

    for (int i = 0; i <= global_variable->LastParticleIndex; i++)
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
    int OriginalLastParticleIndex = global_variable->LastParticleIndex;
    LastParticleIndex = global_variable->LastParticleIndex; // for formPrimordialBinaires function by EW 2025.3.11
    formPrimordialBinaries(OriginalLastParticleIndex);
    assert(OriginalLastParticleIndex <= global_variable->LastParticleIndex); // for debugging by EW 2025.1.4
    assert(CMPtclWorker.empty());                           // for debugging by EW 2025.1.4
    if (OriginalLastParticleIndex != global_variable->LastParticleIndex)
    {
        std::cout << "In total, " << global_variable->LastParticleIndex - OriginalLastParticleIndex
                  << " primordial binaries are created." << std::endl;
        queue_scheduler.initialize(MakePrimordialGroup);
        for (int i = OriginalLastParticleIndex + 1; i <= global_variable->LastParticleIndex; i++)
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

        // Erase member particles from ParticleIndices by EW 2025.3.11
        ParticleIndices.erase(
            std::remove_if(ParticleIndices.begin(), ParticleIndices.end(),
                [](int i) {
                return !particles[i].isActive;
                }
            ),
            ParticleIndices.end()
        );

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


    initializeTime(queue_scheduler, workers, ParticleIndices);

    /* Particle Initialization Check */
    {
        //, NextRegTime= %.3e Myr(%llu),
        for (int i=0; i<=global_variable->LastParticleIndex; i++) {
            ptcl = &particles[i];
            fprintf(nbpout, "PID=%d,PI=%d, (%d)=",ptcl->PID,ptcl->ParticleIndex,ptcl->NumberOfNeighbor);
            for (int j=0;j<ptcl->NumberOfNeighbor;j++) {
                fprintf(nbpout, "%d, ",ptcl->Neighbors[j]);
            }
            fprintf(nbpout, "\n");
        }
    }
    for (int i=0; i<=global_variable->LastParticleIndex; i++) {
        ptcl = &particles[i];
        fprintf(nbpout, "PID=%d(%d), CurrentTime (Irr, Reg) = (%.3e(%llu), %.3e(%llu)) Myr\n"\
                "dtIrr = %.4e Myr, dtReg = %.4e Myr, blockIrr=%llu (%d), blockReg=%llu (%d)\n"\
                "NumNeighbor= %d\n",
                ptcl->PID,
                ptcl->ParticleIndex,
                ptcl->CurrentTimeIrr*global_variable->EnzoTimeStep*1e10/1e6,
                ptcl->CurrentBlockIrr,
                ptcl->CurrentTimeReg*global_variable->EnzoTimeStep*1e10/1e6,
                ptcl->CurrentBlockReg,
                //NextRegTimeBlock*time_step*EnzoTimeStep*1e10/1e6,
                //NextRegTimeBlock,
                ptcl->TimeStepIrr*global_variable->EnzoTimeStep*1e10/1e6,
                ptcl->TimeStepReg*global_variable->EnzoTimeStep*1e10/1e6,
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







void initializeTime(QueueScheduler &queue_scheduler, Worker *workers, std::vector<int> ParticleIndices) {
    Particle *ptcl;
    int min_time_level=0;

    // Initialize Time Step
    queue_scheduler.initialize(InitTime);
    queue_scheduler.takeQueue(ParticleIndices);
    do
    {
        queue_scheduler.assignQueueAuto();
        queue_scheduler.runQueueAuto();
        queue_scheduler.waitQueue(0); // blocking wait
    } while (queue_scheduler.isComplete());

    /* timestep correction */
    {
        std::cout << "Time Step correction." << std::endl;
        for (int i = 0; i <= global_variable->LastParticleIndex; i++)
        {
            ptcl = &particles[i];

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
        global_variable->time_block = std::max(-60, min_time_level - MIN_LEVEL_BUFFER);
        global_variable->block_max = static_cast<ULL>(pow(2, -global_variable->time_block));
        global_variable->time_step = pow(2, global_variable->time_block);

        for (int i = 0; i <= global_variable->LastParticleIndex; i++)
        {
            ptcl = &particles[i];

            ptcl->TimeBlockIrr = static_cast<ULL>(pow(2, ptcl->TimeLevelIrr - global_variable->time_block));
            ptcl->TimeBlockReg = static_cast<ULL>(pow(2, ptcl->TimeLevelReg - global_variable->time_block));
#ifdef IRR_TEST
            ptcl->TimeStepReg = 1;
            ptcl->TimeLevelReg = 0;
            ptcl->TimeBlockReg = global_variable->block_max;
#endif
            ptcl->NextBlockIrr = ptcl->CurrentBlockIrr + ptcl->TimeBlockIrr; // of this particle
        }
        std::cout << "Time Step done." << std::endl;
    }

    /* timestep variable synchronization */
    {
        fprintf(stderr, "nbody+:time_block = %d, EnzoTimeStep=%e\n", global_variable->time_block, global_variable->EnzoTimeStep);
        fprintf(nbpout, "nbody+:time_block = %d, EnzoTimeStep=%e\n", global_variable->time_block, global_variable->EnzoTimeStep);
        //fflush(stderr);
    }

    /*
    {
        for (int i = 0; i <= global_variable->LastParticleIndex; i++)
        {
            ptcl = &particles[i];
            fprintf(nbpout, "PID=%d, CurrentTime (Irr, Reg) = (%.3e(%llu), %.3e(%llu)) Myr\n"
                            "dtIrr = %.4e Myr, dtReg = %.4e Myr, blockIrr=%llu (%d), blockReg=%llu (%d)\n"
                            "NumNeighbor= %d\n",
                    ptcl->PID,
                    ptcl->CurrentTimeIrr * global_variable->EnzoTimeStep * 1e10 / 1e6,
                    ptcl->CurrentBlockIrr,
                    ptcl->CurrentTimeReg * global_variable->EnzoTimeStep * 1e10 / 1e6,
                    ptcl->CurrentBlockReg,
                    // NextRegTimeBlock*time_step*EnzoTimeStep*1e10/1e6,
                    // NextRegTimeBlock,
                    ptcl->TimeStepIrr * global_variable->EnzoTimeStep * 1e10 / 1e6,
                    ptcl->TimeStepReg * global_variable->EnzoTimeStep * 1e10 / 1e6,
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


void InitializationAfterCommunication(QueueScheduler &queue_scheduler, Worker *workers) {
    
    /* Initialize New Particle */ 
    /*  Neighbor inclusion might be needed (to be updated) */

    /* Initialize Particle Attributes */
    /* Since we're not doing full-initialization, we have to do more work on time steps
    e.g., if enzo time can be smaller than regualr time steps. we gotta re-normalize it.
    but this part is not complete yet. */
    Particle* ptcl;

    // Example code by EW 2025.3.18
    if (newNumberOfSingleParticle > 0) {

        InitializationOnGPU(queue_scheduler, workers); // GPU Initialization code

    } else {

        // Particle *ptcl;

        for (int i = 0; i <= global_variable->LastParticleIndex; i++) {

            ptcl = &particles[i];

            ptcl->setNewTimeStepWithNewEnzoTimeStep(global_variable->OldEnzoTimeStep, global_variable->EnzoTimeStep);
            ptcl->CurrentTimeIrr = 0.;
            ptcl->CurrentBlockIrr = 0;
            ptcl->CurrentTimeReg = 0.;
            ptcl->CurrentBlockReg = 0;
            ptcl->NewCurrentBlockIrr = 0;
            ptcl->NextBlockIrr = ptcl->CurrentBlockIrr + ptcl->TimeBlockIrr; // of this particle
#ifdef FEWBODY
            if (ptcl->isCMptcl) { // We have to reset the SDAR clock;
                Queue queue;
                int rank = CMPtclWorker[ptcl->ParticleIndex];
                queue.task = ResetSDARTime;
                queue.pid = ptcl->ParticleIndex;
                workers[rank].addQueue(queue);
                workers[rank].runQueue();
                workers[rank].callback();
            }
#endif
        }
    }

    

    //initializeTime(queue_scheduler, workers, ParticleIndices);
    // ParticleIndices.clear(); // commented out by EW 2025.3.11
    // ParticleIndices.shrink_to_fit(); // commented out by EW 2025.3.11

    for (int i = 0; i <= global_variable->LastParticleIndex; i++)
    {
        ptcl = &particles[i];
        if (!ptcl->isActive)
            continue;
        fprintf(nbpout, "PID=%d, CurrentTime (Irr, Reg) = (%.3e(%llu), %.3e(%llu)) Myr\n"
                        "dtIrr = %.4e Myr, dtReg = %.4e Myr, blockIrr=%llu (%d), blockReg=%llu (%d)\n"
                        "NumNeighbor= %d\n",
                ptcl->PID,
                ptcl->CurrentTimeIrr * global_variable->EnzoTimeStep * 1e10 / 1e6,
                ptcl->CurrentBlockIrr,
                ptcl->CurrentTimeReg * global_variable->EnzoTimeStep * 1e10 / 1e6,
                ptcl->CurrentBlockReg,
                // NextRegTimeBlock*time_step*EnzoTimeStep*1e10/1e6,
                // NextRegTimeBlock,
                ptcl->TimeStepIrr * global_variable->EnzoTimeStep * 1e10 / 1e6,
                ptcl->TimeStepReg * global_variable->EnzoTimeStep * 1e10 / 1e6,
                ptcl->TimeBlockIrr,
                ptcl->TimeLevelIrr,
                ptcl->TimeBlockReg,
                ptcl->TimeLevelReg,
                ptcl->NumberOfNeighbor);
        /*
        Particle* members = ptcl;
        if (ptcl->PID == 2259065) {

            fprintf(nbpout, "PID: %d. Position (pc) - x:%e, y:%e, z:%e, \n", members->PID, members->Position[0]*position_unit, members->Position[1]*position_unit, members->Position[2]*position_unit);
            fprintf(nbpout, "PID: %d. Velocity (km/s) - vx:%e, vy:%e, vz:%e, \n", members->PID, members->Velocity[0]*velocity_unit/yr*pc/1e5, members->Velocity[1]*velocity_unit/yr*pc/1e5, members->Velocity[2]*velocity_unit/yr*pc/1e5);
            fprintf(nbpout, "PID: %d. Mass (Msol) - %e, \n", members->PID, members->Mass*mass_unit);
            fprintf(nbpout, "PID: %d. RadiusOfNeighbor - %e pc, \n", members->PID, sqrt(members->RadiusOfNeighbor)*position_unit);
            fprintf(nbpout, "PID: %d. Background Acceleration - ax%e, ay:%e, az:%e \n", members->PID, members->BackgroundAcceleration[0], members->BackgroundAcceleration[1], members->BackgroundAcceleration[2]);
            fprintf(nbpout, "PID: %d. Total Acceleration - ax:%e, ay:%e, az:%e \n", members->PID, members->a_tot[0][0], members->a_tot[1][0], members->a_tot[2][0]);
            fprintf(nbpout, "PID: %d. Total Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->PID, members->a_tot[0][1], members->a_tot[1][1], members->a_tot[2][1]);
            fprintf(nbpout, "PID: %d. Total Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->PID, members->a_tot[0][2], members->a_tot[1][2], members->a_tot[2][2]);
            fprintf(nbpout, "PID: %d. Total Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->PID, members->a_tot[0][3], members->a_tot[1][3], members->a_tot[2][3]);
            fprintf(nbpout, "PID: %d. Reg Acceleration - ax:%e, ay:%e, az:%e, \n", members->PID, members->a_reg[0][0], members->a_reg[1][0], members->a_reg[2][0]);
            fprintf(nbpout, "PID: %d. Reg Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->PID, members->a_reg[0][1], members->a_reg[1][1], members->a_reg[2][1]);
            fprintf(nbpout, "PID: %d. Reg Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->PID, members->a_reg[0][2], members->a_reg[1][2], members->a_reg[2][2]);
            fprintf(nbpout, "PID: %d. Reg Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->PID, members->a_reg[0][3], members->a_reg[1][3], members->a_reg[2][3]);
            fprintf(nbpout, "PID: %d. Irr Acceleration - ax:%e, ay:%e, az:%e, \n", members->PID, members->a_irr[0][0], members->a_irr[1][0], members->a_irr[2][0]);
            fprintf(nbpout, "PID: %d. Irr Acceleration - axdot:%e, aydot:%e, azdot:%e, \n", members->PID, members->a_irr[0][1], members->a_irr[1][1], members->a_irr[2][1]);
            fprintf(nbpout, "PID: %d. Irr Acceleration - ax2dot:%e, ay2dot:%e, az2dot:%e, \n", members->PID, members->a_irr[0][2], members->a_irr[1][2], members->a_irr[2][2]);
            fprintf(nbpout, "PID: %d. Irr Acceleration - ax3dot:%e, ay3dot:%e, az3dot:%e, \n", members->PID, members->a_irr[0][3], members->a_irr[1][3], members->a_irr[2][3]);
            fprintf(nbpout, "PID: %d. Time Steps (Myr) - irregular:%e, regular:%e \n", members->PID, members->TimeStepIrr*global_variable->EnzoTimeStep*1e4, members->TimeStepReg*global_variable->EnzoTimeStep*1e4);
            fprintf(nbpout, "PID: %d. Time Blocks - irregular:%llu, regular:%llu \n", members->PID, members->TimeBlockIrr, members->TimeBlockReg);
            fprintf(nbpout, "PID: %d. Time Level - irregular:%d, regular:%d\n", members->PID, members->TimeLevelIrr, members->TimeLevelReg);
            fprintf(nbpout, "PID: %d. Current Blocks - irregular: %llu, regular:%llu \n", members->PID, members->CurrentBlockIrr, members->CurrentBlockReg);
        }
        */

    }
    fflush(nbpout);
}