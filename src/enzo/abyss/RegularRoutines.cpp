
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

#define DEBUG_ABYSS

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

#ifdef DEBUG_ABYSS
    std::cout << "calculateRegAccelerationOnGPU starts" << std::endl;
    std::cout << "RegularList size: " << RegularList.size() << std::endl;
#endif

    calculateRegAccelerationOnGPU(RegularList, queue_scheduler);

#ifdef DEBUG_ABYSS
    std::cout << "calculateRegAccelerationOnGPU ended" << std::endl;
#endif

#ifdef NSIGHT
    nvtxRangePop();
#endif

#ifdef NSIGHT
    nvtxRangePushA("RegCudaUpdate");
#endif

#ifdef DEBUG_ABYSS
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
#ifdef DEBUG_ABYSS
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

#ifdef DEBUG_ABYSS
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
void RegularRoutines(QueueScheduler &queue_scheduler, Worker *workers)
{
    double next_time = NextRegTimeBlock * global_variable->time_step;

#ifdef PERFORMANCETRACE
				start_point_routine = std::chrono::high_resolution_clock::now();
#endif

#ifdef NSIGHT
				nvtxRangePushA("RegForce");
#endif

#ifdef DEBUG
				std::cout << "Regular force starts" << std::endl;
#endif
				// Regular force
				queue_scheduler.initialize(RegForce);
				queue_scheduler.takeQueueRegularList(RegularList);
				do
				{
					queue_scheduler.assignQueueAutoRegularList();
					queue_scheduler.runQueueAuto();
					queue_scheduler.waitQueue(0); // blocking wait
				} while (queue_scheduler.isComplete());
#ifdef DEBUG
				std::cout << "Regular force ended" << std::endl;
#endif

#ifdef NSIGHT
				nvtxRangePop();
#endif

#ifdef PERFORMANCETRACE
                end_point_routine = std::chrono::high_resolution_clock::now();
                performance.RegularForce +=
                    std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
#endif


#ifdef PERFORMANCETRACE
                start_point_routine = std::chrono::high_resolution_clock::now();
#endif

#ifdef NSIGHT
				nvtxRangePushA("RegUpdate");
#endif

#ifdef DEBUG
				std::cout << "update regular starts" << std::endl;
#endif
				// Update Regular
				queue_scheduler.initialize(RegUpdate);
				queue_scheduler.takeQueueRegularList(RegularList);
				do
				{
					queue_scheduler.assignQueueAutoRegularList();
					queue_scheduler.runQueueAuto();
					queue_scheduler.waitQueue(0); // blocking wait
				} while (queue_scheduler.isComplete());
#ifdef DEBUG
				std::cout << "update regular ended" << std::endl;
#endif

#ifdef NSIGHT
				nvtxRangePop();
#endif

#ifdef PERFORMANCETRACE
                end_point_routine = std::chrono::high_resolution_clock::now();
                performance.RegularUpdate +=
                    std::chrono::duration_cast<std::chrono::nanoseconds>(end_point_routine - start_point_routine).count();
#endif


#ifdef DEBUG_ABYSS
    //, NextRegTime= %.3e Myr(%llu),
    //for (int i=0; i<RegularList.size(); i++) {
    for (int index: RegularList){
        Particle *ptcl = &particles[index];
        fprintf(stdout, "PID=%d, CurrentTime (Irr, Reg) = (%.3e(%llu), %.3e(%llu)) Myr, NextReg = %.3e (%llu)\n"\
                "dtIrr = %.4e Myr, dtReg = %.4e Myr, blockIrr=%llu (%d), blockReg=%llu (%d), NextBlockIrr= %.3e(%llu)\n"\
                "NumNeighbor= %d\n",
                ptcl->PID,
                ptcl->CurrentTimeIrr* global_variable->time_step*global_variable->EnzoTimeStep*1e10/1e6,
                ptcl->CurrentBlockIrr,
                ptcl->CurrentTimeReg* global_variable->time_step*global_variable->EnzoTimeStep*1e10/1e6,
                ptcl->CurrentBlockReg,
                NextRegTimeBlock* global_variable->time_step*global_variable->EnzoTimeStep*1e10/1e6,
                NextRegTimeBlock,
                ptcl->TimeStepIrr* global_variable->time_step*global_variable->EnzoTimeStep*1e10/1e6,
                ptcl->TimeStepReg* global_variable->time_step*global_variable->EnzoTimeStep*1e10/1e6,
                ptcl->TimeBlockIrr,
                ptcl->TimeLevelIrr,
                ptcl->TimeBlockReg,
                ptcl->TimeLevelReg,
                ptcl->NextBlockIrr* global_variable->time_step*global_variable->EnzoTimeStep*1e10/1e6,
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
*/
#endif
    }

    // current_time_irr = particles[ThisLevelNode->ParticleList[0]].CurrentBlockIrr;
} // Regular Done.

#endif
