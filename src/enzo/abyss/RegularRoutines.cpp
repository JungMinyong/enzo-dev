
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

void calculateRegAccelerationOnGPU(std::unordered_set<int> RegularList, QueueScheduler &queue_scheduler);

#ifdef CUDA
void RegularRoutines(QueueScheduler &queue_scheduler, Worker *workers)
{
#ifdef PerformanceTrace
    std::chrono::high_resolution_clock::time_point start_point;
    std::chrono::high_resolution_clock::time_point end_point;
#endif
    double next_time = 0;

#ifdef PerformanceTrace
    start_point = std::chrono::high_resolution_clock::now();
#endif
    // total_tasks = RegularList.size();
    next_time = NextRegTimeBlock * global_variable->time_step;

    // fprintf(stdout, "Regular starts\n");
#ifdef NSIGHT
    nvtxRangePushA("calculateRegAccelerationOnGPU");
#endif

#ifdef DEBUG
    std::cout << "calculateRegAccelerationOnGPU starts" << std::endl;
    std::cout << "RegularList size: " << RegularList.size() << std::endl;
#endif

    calculateRegAccelerationOnGPU(RegularList, queue_scheduler);

#ifdef DEBUG
    std::cout << "calculateRegAccelerationOnGPU ended" << std::endl;
#endif

#ifdef NSIGHT
    nvtxRangePop();
#endif

#ifdef NSIGHT
    nvtxRangePushA("RegCudaUpdate");
#endif

#ifdef DEBUG
    std::cout << "update regular starts" << std::endl;
#endif

    // Update Regular
    queue_scheduler.initialize(RegCudaUpdate);
    queue_scheduler.takeQueueRegularList(RegularList);
    do
    {
        queue_scheduler.assignQueueRegularList();
        queue_scheduler.runQueueAuto();
        queue_scheduler.waitQueue(0); // blocking wait
    } while (queue_scheduler.isComplete());
#ifdef DEBUG
    std::cout << "update regular ended" << std::endl;
#endif

#ifdef NSIGHT
    nvtxRangePop();
#endif

#ifdef PerformanceTrace
    end_point = std::chrono::high_resolution_clock::now();
    performance.RegularRoutine +=
        std::chrono::duration_cast<std::chrono::nanoseconds>(end_point - start_point).count();
#endif

#ifdef DEBUG
    {
        Particle *ptcl;
        for (int index: RegularList)
        {
            ptcl = &particles[index];
            fprintf(stdout, "PID=%d, CurrentTime (Irr, Reg) = (%.3e(%llu), %.3e(%llu)) Myr, NextReg = %.3e (%llu)\n"
                            "dtIrr = %.4e Myr, dtReg = %.4e Myr, blockIrr=%llu (%d), blockReg=%llu (%d), NextBlockIrr= %.3e(%llu)\n"
                            "NumNeighbor= %d\n",
                    ptcl->PID,
                    ptcl->CurrentTimeIrr * global_variable->EnzoTimeStep * 1e10 / 1e6,
                    ptcl->CurrentBlockIrr,
                    ptcl->CurrentTimeReg * global_variable->EnzoTimeStep * 1e10 / 1e6,
                    ptcl->CurrentBlockReg,
                    NextRegTimeBlock * global_variable->time_step * global_variable->EnzoTimeStep * 1e10 / 1e6,
                    NextRegTimeBlock,
                    ptcl->TimeStepIrr * global_variable->EnzoTimeStep * 1e10 / 1e6,
                    ptcl->TimeStepReg * global_variable->EnzoTimeStep * 1e10 / 1e6,
                    ptcl->TimeBlockIrr,
                    ptcl->TimeLevelIrr,
                    ptcl->TimeBlockReg,
                    ptcl->TimeLevelReg,
                    ptcl->NextBlockIrr * global_variable->time_step * global_variable->EnzoTimeStep * 1e10 / 1e6,
                    ptcl->NextBlockIrr,
                    ptcl->NumberOfNeighbor);

            fprintf(stdout, " a_tot = (%.4e,%.4e,%.4e), a_reg = (%.4e,%.4e,%.4e), a_irr = (%.4e,%.4e,%.4e), n_n=%d, R=%.3e\n\
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
                    ptcl->a_irr[2][3]);
        }
        // fflush(stdout);
    }
#endif // endif debug
}
#else
void RegularRoutines()
{
    // Regular Gravity
    task = RegForce;
    completed_tasks = 0;
    total_tasks = RegularList.size();
    next_time = NextRegTimeBlock * global_variable->time_step;

    // std::cout << "TotalTask=" << total_tasks << std::endl;

    /*
    std::cout << "RegularList, PID= ";
    for (int i=0; i<total_tasks; i++) {
        std::cout << RegularList[i]<< ", ";
    }*/
    // std::cout << std::endl;

    InitialAssignmentOfTasks(task, total_tasks, TASK_TAG);
    InitialAssignmentOfTasks(RegularList, next_time, total_tasks, PTCL_TAG);
    MPI_Waitall(NumberOfCommunication, requests, statuses);
    NumberOfCommunication = 0;

    // further assignments
    remaining_tasks = total_tasks - NumberOfWorker;
    while (completed_tasks < total_tasks)
    {
        // Check which worker is done
        MPI_Irecv(&ptcl_id_return, 1, MPI_INT, MPI_ANY_SOURCE, TERMINATE_TAG, abyss_comm, &request);
        // Poll until send completes
        /*
             flag=0;
             while (!flag) {
             MPI_Test(&request, &flag, &status);
        // Perform other work while waiting
        }
        */
        MPI_Wait(&request, &status);
        completed_rank = status.MPI_SOURCE;
        // printf("Rank %d: Send operation completed (%d).\n",completed_rank, ptcl_id_return);

        if (remaining_tasks > 0)
        {
            ptcl_id = RegularList[NumberOfWorker + completed_tasks];
            MPI_Send(&task, 1, MPI_INT, completed_rank, TASK_TAG, abyss_comm);
            MPI_Send(&ptcl_id, 1, MPI_INT, completed_rank, PTCL_TAG, abyss_comm);
            MPI_Send(&next_time, 1, MPI_DOUBLE, completed_rank, TIME_TAG, abyss_comm);
            remaining_tasks--;
        }
        else
        {
            // printf("Rank %d: No more tasks to assign\n", completed_rank);
        }
        // updateSkipList(skiplist, ptcl_id_return);
        completed_tasks++;
    }

    // ParticleSynchronization();

    // Regular Update
    // std::cout<< "Reg Acc Done." <<std::endl;
    task = RegUpdate;
    completed_tasks = 0;

    InitialAssignmentOfTasks(task, total_tasks, TASK_TAG);
    InitialAssignmentOfTasks(RegularList, total_tasks, PTCL_TAG);
    MPI_Waitall(NumberOfCommunication, requests, statuses);
    NumberOfCommunication = 0;

    // further assignments
    remaining_tasks = total_tasks - NumberOfWorker;
    while (completed_tasks < total_tasks)
    {
        // Check which worker is done
        MPI_Irecv(&ptcl_id_return, 1, MPI_INT, MPI_ANY_SOURCE, TERMINATE_TAG, abyss_comm, &request);
        MPI_Wait(&request, &status);
        completed_rank = status.MPI_SOURCE;

        if (remaining_tasks > 0)
        {
            ptcl_id = RegularList[NumberOfWorker + completed_tasks];
            // MPI_Isend(&task,      1, MPI_INT, completed_rank, TASK_TAG, abyss_comm, &request);
            // MPI_Isend(&ptcl_id,   1, MPI_INT, completed_rank, PTCL_TAG, abyss_comm, &request);
            MPI_Send(&task, 1, MPI_INT, completed_rank, TASK_TAG, abyss_comm);
            MPI_Send(&ptcl_id, 1, MPI_INT, completed_rank, PTCL_TAG, abyss_comm);
            remaining_tasks--;
        }
        else
        {
            // printf("Rank %d: No more tasks to assign\n", completed_rank);
        }
        completed_tasks++;
    }

        //, NextRegTime= %.3e Myr(%llu),
        for (int i=0; i<total_tasks; i++) {
            ptcl = &particles[RegularList[i]];
            fprintf(stdout, "PID=%d, CurrentTime (Irr, Reg) = (%.3e(%llu), %.3e(%llu)) Myr, NextReg = %.3e (%llu)\n"\
                    "dtIrr = %.4e Myr, dtReg = %.4e Myr, blockIrr=%llu (%d), blockReg=%llu (%d), NextBlockIrr= %.3e(%llu)\n"\
                    "NumNeighbor= %d\n",
                    ptcl->PID,
                    ptcl->CurrentTimeIrr*EnzoTimeStep*1e10/1e6,
                    ptcl->CurrentBlockIrr,
                    ptcl->CurrentTimeReg*EnzoTimeStep*1e10/1e6,
                    ptcl->CurrentBlockReg,
                    NextRegTimeBlock*time_step*EnzoTimeStep*1e10/1e6,
                    NextRegTimeBlock,
                    ptcl->TimeStepIrr*EnzoTimeStep*1e10/1e6,
                    ptcl->TimeStepReg*EnzoTimeStep*1e10/1e6,
                    ptcl->TimeBlockIrr,
                    ptcl->TimeLevelIrr,
                    ptcl->TimeBlockReg,
                    ptcl->TimeLevelReg,
                    ptcl->NextBlockIrr*time_step*EnzoTimeStep*1e10/1e6,
                    ptcl->NextBlockIrr,
                    ptcl->NumberOfNeighbor
                    );

/*
            fprintf(stdout, " a_tot = (%.4e,%.4e,%.4e), a_reg = (%.4e,%.4e,%.4e), a_irr = (%.4e,%.4e,%.4e), n_n=%d, R=%.3e\n\
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
                    ptcl->NewNumberOfNeighbor,
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
        //fflush(stdout);
    }

    */
    // current_time_irr = particles[ThisLevelNode->ParticleList[0]].CurrentBlockIrr;
} // Regular Done.

#endif
