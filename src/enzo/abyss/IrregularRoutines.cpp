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

bool createSkipList(SkipList *skiplist);
bool updateSkipList(SkipList *skiplist, int ptcl_id);

void formBinaries(std::vector<int>& ParticleList, std::vector<int>& newCMptcls, std::unordered_map<int, int>& existing, std::unordered_map<int, int>& terminated);
void FBTermination(Particle* ptclCM);
void Merge(Particle* p1, Particle* p2);

#ifdef SEVN_BINARY
bool makeSEVNBinary(Particle* ptclCM);
void deleteSEVNBinary(Particle* ptclCM);
void BinaryEvolution(Particle* ptclCM);
#endif

bool IrregularRoutines(QueueScheduler &queue_scheduler, Worker *workers, std::unordered_set<int>& RegularList) {
#ifdef PerformanceTrace
    std::chrono::high_resolution_clock::time_point start_point;
    std::chrono::high_resolution_clock::time_point end_point;
#endif
    int max_level = 20;
    double prob = 0.5;
    SkipList *skiplist;
    Node *ThisLevelNode;
    Particle *ptcl;
    double current_time_irr = 0;
    double next_time = 0;
    Worker* worker;

#ifdef FEWBODY
    bool bin_termination = false;
	bool new_binaries = false;
    std::vector<int> newCMptcls; // by EW 2025.1.6 // unordered_set? by EW 2025.1.11
#endif

#ifdef NSIGHT
    nvtxRangePushA("createSkipList");
#endif

#ifdef DEBUG_ABYSS
    fprintf(nbpout, "createSkipList starts\n");
    fflush(nbpout);
#endif

    skiplist = new SkipList(max_level, prob);
    if (createSkipList(skiplist) == false)
        fprintf(stderr, "There are no irregular particles!\nBut is it really happening? check skiplist->display()\n");

#ifdef DEBUG_ABYSS
    fprintf(nbpout, "createSkipList ends\n");
    fflush(nbpout);
#endif

#ifdef NSIGHT
    nvtxRangePop();
#endif

    /* Irregular Loop starts */
    while (skiplist->getFirstNode() != nullptr)
    {
        // update_idx=0; // commented out by EW 2025.1.6
        ThisLevelNode = skiplist->getFirstNode();
#ifdef DEBUG_ABYSS
        fprintf(nbpout, "After getFirstNode... Irregular list size: %d\n", ThisLevelNode->ParticleList.size());
        fflush(nbpout);
        // ptcl = &particles[ThisLevelNode->ParticleList[0]];
        /*
        for (int i = 0; i < ThisLevelNode->ParticleList.size(); i++) {
            ptcl = &particles[ThisLevelNode->ParticleList[i]];
            // fprintf(nbpout, "PID: %d. RadiusOfNeighbor: %e\n", ptcl->PID, ptcl->RadiusOfNeighbor);
            if (ptcl->isActive)
                fprintf(stderr, "PID: %d. Mass: %e Msun\n", ptcl->PID, ptcl->Mass * mass_unit);
        }
        */
#endif
        ThisLevelNode->ParticleList.erase(
            std::remove_if(ThisLevelNode->ParticleList.begin(), ThisLevelNode->ParticleList.end(),
                           [](int i)
                           {
                               return !particles[i].isActive;
                           }),
            ThisLevelNode->ParticleList.end()
        );
#ifdef DEBUG_ABYSS
        fprintf(nbpout, "2. Irregular list size: %d\n", ThisLevelNode->ParticleList.size());
        fflush(nbpout);
#endif
        if (ThisLevelNode->ParticleList.size() == 0) {
            skiplist->deleteFirstNode();
            continue;
        }

        next_time = particles[ThisLevelNode->ParticleList[0]].CurrentTimeIrr + particles[ThisLevelNode->ParticleList[0]].TimeStepIrr;

#ifdef DEBUG_ABYSS
        // print out particlelist
        // fprintf(stdout, "(IRR_FORCE) next_time: %e Myr\n", next_time*global_variable->EnzoTimeStep*1e4);
        fprintf(nbpout, "(IRR_FORCE) next_time: %e Myr\n", next_time*global_variable->EnzoTimeStep*1e4);
        fprintf(nbpout, "Irregular list size: %d\n", ThisLevelNode->ParticleList.size());
        fflush(nbpout);
        /*
        fprintf(stdout, "PID: %d. CurrentTimeIrr: %e Myr, TimeStepIrr: %e Myr\n", 
                    particles[ThisLevelNode->ParticleList[0]].PID, 
                    particles[ThisLevelNode->ParticleList[0]].CurrentTimeIrr*global_variable->EnzoTimeStep*1e4, 
                    particles[ThisLevelNode->ParticleList[0]].TimeStepIrr*global_variable->EnzoTimeStep*1e4);

        // fprintf(stdout, "PID (%d) = ", ThisLevelNode->ParticleList.size());
        for (int i=0; i<ThisLevelNode->ParticleList.size(); i++) {
            ptcl = &particles[ThisLevelNode->ParticleList[i]];
            // fprintf(stdout, "%d, ", ptcl->PID);
            fprintf(stdout, "PID: %d. %e Myr, %e Myr\n", 
                    ptcl->PID,
                    ptcl->CurrentTimeIrr*global_variable->EnzoTimeStep*1e4,
                    ptcl->TimeStepIrr*global_variable->EnzoTimeStep*1e4);
        }
        fprintf(stdout, "\n");
        // fflush(stdout);
        */
#endif

        // Irregular Force
#ifdef FEWBODY
#ifdef NSIGHT
        nvtxRangePushA("IrregularForce");
#endif
#ifdef DEBUG_ABYSS
        fprintf(nbpout, "Irr force starts\n");
        fflush(nbpout);
#endif
        int cm_pid;
        Queue queue;
        queue_scheduler.initializeIrr(IrrForce, next_time, ThisLevelNode->ParticleList);
        auto iter = queue_scheduler.CMPtcls.begin();
#ifdef SEVN_BINARY
        std::vector<int> CMPtclsForSEVN;
#endif
        do
        {
            queue_scheduler.assignQueueAuto();
            queue_scheduler.runQueueAuto();
            // queue_scheduler.printStatus();
            do
            {
                worker = queue_scheduler.waitQueue(1); // non-blocking wait
                // if there's any CMPtcl
                if (queue_scheduler.CMPtcls.size() > 0)
                {
                    /* check if there's any CM ptcl ready to go for SDAR*/
                    if (iter == queue_scheduler.CMPtcls.end())
                        iter = queue_scheduler.CMPtcls.begin();
                    cm_pid = *(iter);
                    ptcl = &particles[cm_pid];
                    /* // commented out by EW 2025.3.10
                    for (int j = 0; j < ptcl->NumberOfNeighbor; j++)
                    {
                        // if (particles[ptcl->Neighbors[j]].isUpdateToDate == false) // original code
                        if (particles[ptcl->Neighbors[j]].isActive && !particles[ptcl->Neighbors[j]].isUpdateToDate) // modified by EW 2025.2.26
                        {
                            iter++;
                            goto skip_to_next;
                        }
                    }
                    */
                    // queue_scheduler.printFreeWorker();
                    // queue_scheduler.printWorkerToGo();
                    // std::cout << "before: The number of CM ptcl is " << queue_scheduler.CMPtcls.size() << std::endl;
                    queue.task = ARIntegration;
                    queue.pid = cm_pid;
                    queue.next_time = next_time;
                    workers[CMPtclWorker[cm_pid]].addQueue(queue);
                    queue_scheduler.assignWorker(&workers[CMPtclWorker[cm_pid]]);
                    iter = queue_scheduler.CMPtcls.erase(iter);
#ifdef SEVN_BINARY
                    if (ptcl->BinaryEvolution != nullptr)
                        CMPtclsForSEVN.push_back(cm_pid);
#endif
                    // std::cout << "after: The number of CM ptcl is " << queue_scheduler.CMPtcls.size() << std::endl;
                    // queue_scheduler.printFreeWorker();
                    // queue_scheduler.printWorkerToGo();
                // skip_to_next:;
                }
                if (worker != nullptr)
                {
                    // fprintf(stdout, "Worker rank: %d\n", worker->AbyssProcessorNumber);
                }
            } while (worker == nullptr);
            queue_scheduler.callback(worker);
            // queue_scheduler.printStatus();
        } while (queue_scheduler.isComplete());
#ifdef DEBUG_ABYSS
        fprintf(nbpout, "Irregular Force done\n");
        fflush(nbpout);
#endif
#ifdef NSIGHT
        nvtxRangePop();
#endif
#else
        queue_scheduler.initialize(IrrForce, next_time);
        queue_scheduler.takeQueue(ThisLevelNode->ParticleList);
        do
        {
            queue_scheduler.assignQueueAuto();
            queue_scheduler.runQueueAuto();
            queue_scheduler.waitQueue(0); // blocking wait
        } while (queue_scheduler.isComplete());
#endif

        // if remaining, further sorting
        /*
        if (update_idx < total_tasks) {
                updateSkipList(skiplist, ThisLevelNode->ParticleList[update_idx]);
                update_idx++;
            }				updateSkipList(skiplist, ptcl_id_return);
            */
#ifdef NSIGHT
        nvtxRangePushA("IrregularUpdate");
#endif
/*
        // Irregular Update
        queue_scheduler.initialize(IrrUpdate);
        queue_scheduler.takeQueue(ThisLevelNode->ParticleList);
        do
        {
            queue_scheduler.assignQueueAuto();
            queue_scheduler.runQueueAuto();
            queue_scheduler.waitQueue(0); // blocking wait
        } while (queue_scheduler.isComplete());
*/
        for (int ptcl_id : ThisLevelNode->ParticleList)
        {
            ptcl = &particles[ptcl_id];

            if (ptcl->NumberOfNeighbor != 0) // IAR modified
                ptcl->updateParticle();
            ptcl->CurrentBlockIrr = ptcl->NewCurrentBlockIrr;
            ptcl->CurrentTimeIrr = ptcl->CurrentBlockIrr * global_variable->time_step;
        }
#ifdef DEBUG_ABYSS
        // for (int i : ThisLevelNode->ParticleList)
        // {
        //     ptcl = &particles[i];
        //     if (ptcl->CurrentTimeIrr != next_time)
        //     {
        //         // fprintf(stdout, "Error! PID: %d, CurrentTimeIrr: %e Myr, next_time: %e Myr\n", ptcl->PID, ptcl->CurrentTimeIrr * global_variable->EnzoTimeStep * 1e4, next_time * global_variable->EnzoTimeStep * 1e4);
        //         fprintf(nbpout, "Error! PID: %d, CurrentTimeIrr: %e Myr, next_time: %e Myr\n", ptcl->PID, ptcl->CurrentTimeIrr * global_variable->EnzoTimeStep * 1e4, next_time * global_variable->EnzoTimeStep * 1e4);
        //         fflush(nbpout);
        //         assert(ptcl->CurrentTimeIrr == next_time);
        //     }
        // }
        // std::cout << "Irregular update done" << std::endl;
        fprintf(nbpout, "Irregular update done\n");
        fflush(nbpout);
#endif
#ifdef NSIGHT
        nvtxRangePop();
#endif

#ifdef FEWBODY

#ifdef NSIGHT
        nvtxRangePushA("FewBodyTermination");
#endif
        int OriginalParticleListSize;
        OriginalParticleListSize = ThisLevelNode->ParticleList.size();
        for (int i = 0; i < OriginalParticleListSize; i++)
        {
            ptcl = &particles[ThisLevelNode->ParticleList[i]];
            if (ptcl->getBinaryInterruptState() == BinaryInterruptState::merger ||
                ptcl->getBinaryInterruptState() == BinaryInterruptState::terminated)
            {

                assert(ptcl->isCMptcl); // for debugging by EW 2025.1.20

                if (ptcl->getBinaryInterruptState() == BinaryInterruptState::merger)
                {
                    if (ptcl->NumberOfMember == 2)
                    { // binary merger

                        Particle *donor = &particles[ptcl->Members[0]];
                        Particle *accretor = &particles[ptcl->Members[1]];

                        Merge(donor, accretor);
#ifdef SEVN
                        if (donor->StellarEvolution != nullptr && accretor->StellarEvolution != nullptr) {
                            fprintf(nbpout, "After Merge... Donor (PID: %d). Mass: %e Msun,  StellarEvolution->get_zams: %e Msun\n", donor->PID, donor->Mass*mass_unit, donor->StellarEvolution->get_zams());
                            fprintf(nbpout, "After Merge... Accretor (PID: %d). Mass: %e Msun, StellarEvolution->get_zams: %e Msun\n", accretor->PID, accretor->Mass*mass_unit, accretor->StellarEvolution->get_zams());
                            fprintf(nbpout, "Donor: amiempty(): %d\n", donor->StellarEvolution->amiempty());
                            fprintf(nbpout, "Accretor: amiempty(): %d\n", accretor->StellarEvolution->amiempty());
                            fflush(nbpout);
                        }
#endif
                    }
                    else
                    { // from NewFBInitialization3

                        assert(ptcl->NumberOfMember > 2); // for debugging by EW 2025.1.20

                        Particle *donor;
                        Particle *accretor;

                        for (int j = 0; j < ptcl->NumberOfMember; j++)
                        {
                            if (particles[ptcl->Members[j]].getBinaryInterruptState() == BinaryInterruptState::collision)
                            {
                                donor = &particles[ptcl->Members[j]];
                                accretor = &particles[donor->getBinaryPairID()];

                                assert(accretor->getBinaryInterruptState() == BinaryInterruptState::collision);
                                assert(accretor->getBinaryPairID() == donor->ParticleIndex);
                                break;
                            }
                        }

                        Merge(donor, accretor);
#ifdef SEVN
                        if (donor->StellarEvolution != nullptr && accretor->StellarEvolution != nullptr) {
                            fprintf(nbpout, "After Merge... Donor (PID: %d). Mass: %e Msun,  StellarEvolution->get_zams: %e Msun\n", donor->PID, donor->Mass*mass_unit, donor->StellarEvolution->get_zams());
                            fprintf(nbpout, "After Merge... Accretor (PID: %d). Mass: %e Msun, StellarEvolution->get_zams: %e Msun\n", accretor->PID, accretor->Mass*mass_unit, accretor->StellarEvolution->get_zams());
                            fprintf(nbpout, "Donor: amiempty(): %d\n", donor->StellarEvolution->amiempty());
                            fprintf(nbpout, "Accretor: amiempty(): %d\n", accretor->StellarEvolution->amiempty());
                            fflush(nbpout);
                        }
#endif
                        int rank = CMPtclWorker[ptcl->ParticleIndex];
                        queue.task = MergeManyBody;
                        queue.pid = ptcl->ParticleIndex;
                        workers[rank].addQueue(queue);
                        workers[rank].runQueue();
                        workers[rank].callback();

#ifdef SEVN // newly added by EW 2025.5.14 // not updated in Abyss code yet!!!
                        Particle* ptcl_erased = donor->Mass < 0.0 ? donor : accretor;
                        fprintf(nbpout, "ptcl_erased... PID: %d\n", ptcl_erased->PID);
                        if (ptcl_erased->StellarEvolution != nullptr) {

                            auto it = SEVNList.begin();
                            while (it != SEVNList.end()) {
                                if (it->second == ptcl_erased->ParticleIndex) {
                                    it = SEVNList.erase(it);
                                    fprintf(nbpout, "Merger induced zero mass particle (PID: %d) is deleted from SEVNList\n", ptcl_erased->PID);
                                    break;
                                }
                                else
                                    it++;
                            }

                            delete ptcl_erased->StellarEvolution;
                            ptcl_erased->StellarEvolution = nullptr;
                            fprintf(nbpout, "Merger induced zero mass particle (PID: %d) SEVN memory is free now\n", ptcl_erased->PID);
                        }
                        fflush(nbpout);
#endif
                        continue;
                    }
                }

                bin_termination = true;

                if (ptcl->ParticleIndex == global_variable->LastParticleIndex)
                {
                    global_variable->LastParticleIndex--;
                    //global_variable->LastParticleIndex == LastParticleIndex;
                }
                else
                    PrevCMPtclWorker.insert({ptcl->ParticleIndex, CMPtclWorker[ptcl->ParticleIndex]});
                CMPtclWorker.erase(ptcl->ParticleIndex);

                for (int j = 0; j < ptcl->NumberOfMember; j++)
                {
                    particles[ptcl->Members[j]].CMPtclIndex = -1;
                    if (particles[ptcl->Members[j]].Mass < 0.0) {
#ifdef SEVN
                        Particle* ptcl_erased = &particles[ptcl->Members[j]];
                        fprintf(nbpout, "ptcl_erased... PID: %d\n", ptcl_erased->PID);
                        if (ptcl_erased->StellarEvolution != nullptr) {

                            auto it = SEVNList.begin();
                            while (it != SEVNList.end()) {
                                if (it->second == ptcl_erased->ParticleIndex) {
                                    it = SEVNList.erase(it);
                                    fprintf(nbpout, "Merger induced zero mass particle (PID: %d) is deleted from SEVNList\n", ptcl_erased->PID);
                                    break;
                                }
                                else
                                    it++;
                            }

                            delete ptcl_erased->StellarEvolution;
                            ptcl_erased->StellarEvolution = nullptr;
                            fprintf(nbpout, "Merger induced zero mass particle (PID: %d) SEVN memory is free now\n", ptcl_erased->PID);
                        }
                        fflush(nbpout);
#endif
                        
                        /* // (Query) EW: Here? I don't think so
                        NumberOfSingleParticle--;
                        AvailableIndices[NumberOfAvailableIndices] = particles[ptcl->Members[j]].ParticleIndex;
                        NumberOfAvailableIndices++;
                        */

                        continue;
                    }
                    ThisLevelNode->ParticleList.push_back(ptcl->Members[j]);
#ifdef DEBUG_ABYSS
                    fprintf(nbpout, "PID: %d is added to ParticleList\n", particles[ptcl->Members[j]].PID);
                    fprintf(nbpout, "ParticleList size: %d\n", ThisLevelNode->ParticleList.size());
                    fflush(nbpout);
#endif
                }
#ifdef DEBUG_ABYSS
                fprintf(nbpout, "FBTermination (PID: %d)\n", ptcl->PID);
                fflush(nbpout);
#endif

#ifdef SEVN_BINARY
                if (ptcl->BinaryEvolution != nullptr) {
                    deleteSEVNBinary(ptcl);
                }
#endif
                RegularList.erase(ptcl->ParticleIndex);
                FBTermination(ptcl);
            }
        }

        if (bin_termination) {
            for (int i=OriginalParticleListSize; i<ThisLevelNode->ParticleList.size(); i++) {
                ptcl = &particles[ThisLevelNode->ParticleList[i]];

                if (ptcl->CurrentBlockReg + ptcl->TimeBlockReg == NextRegTimeBlock)
                    RegularList.insert(ptcl->ParticleIndex);
            }

            // Erase terminated CM particles by EW 2025.1.6
            ThisLevelNode->ParticleList.erase(
                std::remove_if(ThisLevelNode->ParticleList.begin(), ThisLevelNode->ParticleList.end(),
                    [](int i) {
                    return !particles[i].isActive;
                    }
                ),
                ThisLevelNode->ParticleList.end()
            );
#ifdef DEBUG_ABYSS
            fprintf(nbpout, "After bin_termination, ParticleList size: %d\n", ThisLevelNode->ParticleList.size());
            fflush(nbpout);
#endif
        }
#ifdef unused
#ifdef NSIGHT
        nvtxRangePop();
#endif

#ifdef NSIGHT
        nvtxRangePushA("FewBodySearch");
#endif

#ifdef DEBUG_ABYSS
        // std::cout << "FB search starts" << std::endl;
        fprintf(nbpout, "FB search starts\n");
        fflush(nbpout);
#endif
        /*
        // std::cerr << "FB search starts" << std::endl;
        // Few-body group search
        queue_scheduler.initialize(SearchGroup);
        queue_scheduler.takeQueue(ThisLevelNode->ParticleList);
        do
        {
            queue_scheduler.assignQueueAuto();
            queue_scheduler.runQueueAuto();
            // queue_scheduler.printStatus();
            queue_scheduler.waitQueue(0); // blocking wait
        } while (queue_scheduler.isComplete());

        // std::cerr << "FB search ended" << std::endl;
        */
        for (int ptcl_id : ThisLevelNode->ParticleList)
        {
            ptcl = &particles[ptcl_id];

            if (ptcl->getBinaryInterruptState() == BinaryInterruptState::threebody)
            {
                ptcl->setBinaryInterruptState(BinaryInterruptState::none);
            }
            else if (ptcl->getBinaryInterruptState() == BinaryInterruptState::manybody)
            {
                ptcl->NewNumberOfNeighbor = 0;
                ptcl->checkNewGroup2();
                ptcl->setBinaryInterruptState(BinaryInterruptState::none);
            }
            else if (ptcl->NewNumberOfNeighbor != 0)
            {	
                ptcl->checkNewGroup3();
            }
        }
#ifdef DEBUG_ABYSS
        fprintf(nbpout, "FB search ended\n");
        fflush(nbpout);
#endif

#ifdef NSIGHT
        nvtxRangePop();
#endif
#endif // unused

#ifdef NSIGHT
        nvtxRangePushA("FormBinaries");
#endif
        OriginalParticleListSize = ThisLevelNode->ParticleList.size();
        int rank_delete, rank_new;
#ifdef DEBUG_ABYSS
        fprintf(nbpout, "formBinaries starts\n");
        fflush(nbpout);
#endif
        LastParticleIndex = global_variable->LastParticleIndex; // for formBinareis function by EW 2025.3.11
        formBinaries(ThisLevelNode->ParticleList, newCMptcls, CMPtclWorker, PrevCMPtclWorker);
#ifdef DEBUG_ABYSS
        fprintf(nbpout, "formBinaries ended\n");
        fflush(nbpout);
#endif
        if (OriginalParticleListSize != ThisLevelNode->ParticleList.size())
        {
#ifdef DEBUG_ABYSS
            fprintf(nbpout, "New Binary!\n");
            fflush(nbpout);
#endif
            new_binaries = true;

            // new code by EW 2025.1.26
            Particle *ptclCM;
            Particle *mem_ptclCM;
            for (int i = 0; i < newCMptcls.size(); i++)
            {
                ptclCM = &particles[newCMptcls[i]]; // 2025.01.10 edited to newCMptcls[i] by YS

                for (int j = 0; j < ptclCM->NewNumberOfMember; j++)
                {
                    mem_ptclCM = &particles[ptclCM->NewMembers[j]];
                    if (mem_ptclCM->isCMptcl)
                    {
                        fprintf(stdout, "manybody group detected; PID %d should be deleted first\n", mem_ptclCM->PID);
                        rank_delete = CMPtclWorker[mem_ptclCM->ParticleIndex];
                        fprintf(stdout, "Rank of CM ptcl %d: %d\n", mem_ptclCM->PID, rank_delete);
                        queue.task = DeleteGroup;
                        queue.pid = mem_ptclCM->ParticleIndex;
                        workers[rank_delete].addQueue(queue);
                        workers[rank_delete].runQueue();
                        workers[rank_delete].callback();

                        if (mem_ptclCM->ParticleIndex == LastParticleIndex) {
                            LastParticleIndex--;
                            global_variable->LastParticleIndex = LastParticleIndex;
                        }
                        else
                            PrevCMPtclWorker.insert({mem_ptclCM->ParticleIndex, CMPtclWorker[mem_ptclCM->ParticleIndex]});
                        CMPtclWorker.erase(mem_ptclCM->ParticleIndex);
                    }
                }

                rank_new = CMPtclWorker[ptclCM->ParticleIndex];
#ifdef DEBUG_ABYSS
                fprintf(nbpout, "Rank of CM ptcl %d: %d\n", ptclCM->PID, rank_new);
                fflush(nbpout);
#endif

#ifdef SEVN_BINARY
                if (!makeSEVNBinary(ptclCM)) {
                    if (ptclCM->ParticleIndex == LastParticleIndex) {
                        LastParticleIndex--;
                        global_variable->LastParticleIndex = LastParticleIndex;
                    }
                    else
                        PrevCMPtclWorker.insert({ptclCM->ParticleIndex, CMPtclWorker[ptclCM->ParticleIndex]});
                    CMPtclWorker.erase(ptclCM->ParticleIndex);
                    continue;
                }
#endif
                queue.task = MakeGroup;
                queue.pid = ptclCM->ParticleIndex;
                workers[rank_new].addQueue(queue);
                workers[rank_new].runQueue();
                workers[rank_new].callback();

                for (int i = 0; i < ptclCM->NewNumberOfMember; i++) {
                    RegularList.erase(ptclCM->NewMembers[i]);
                    particles[ptclCM->NewMembers[i]].NewNumberOfMember = 0;
                }
                if (ptclCM->CurrentBlockReg + ptclCM->TimeBlockReg == NextRegTimeBlock)
                    RegularList.insert(ptclCM->ParticleIndex); // VERY IMPORTANT BUG FIXED by EW 2025.7.18
            }
#ifdef DEBUG_ABYSS
            fprintf(nbpout, "All new fewbody objects are initialized.\n");
            fflush(nbpout);
#endif

            ThisLevelNode->ParticleList.erase(
                std::remove_if(
                    ThisLevelNode->ParticleList.begin(),
                    ThisLevelNode->ParticleList.end(),
                    [](int i)
                    { return !particles[i].isActive; }),
                ThisLevelNode->ParticleList.end());
        }
        newCMptcls.clear();

        // std::cout << "erase success" << std::endl;
#ifdef NSIGHT
        nvtxRangePop();
#endif

#endif
        OriginalParticleListSize = ThisLevelNode->ParticleList.size();
        for (int i = 0; i < OriginalParticleListSize; i++)
            updateSkipList(skiplist, ThisLevelNode->ParticleList[i]);
#ifdef DEBUG_ABYSS
        fprintf(nbpout, "updateSkipList ended\n");
        fflush(nbpout);
#endif

        // std::cout << "update success" << std::endl;
        // skiplist->display();
        /*
        {
            //, NextRegTime= %.3e Myr(%llu),
            for (int i=0; i<ThisLevelNode->ParticleList.size(); i++) {
                ptcl = &particles[ThisLevelNode->ParticleList[i]];
                fprintf(stdout, "PID=%d, CurrentTime (Irr, Reg) = (%.3e(%llu), %.3e(%llu)) Myr, NextReg = %.3e (%llu)\n"\
                        "dtIrr = %.4e Myr, dtReg = %.4e Myr, blockIrr=%llu (%d), blockReg=%llu (%d), NextBlockIrr= %.3e(%llu)\n"\
                        "NumNeighbor= %d\n",
                        ptcl->PID,
                        ptcl->CurrentTimeIrr*global_variable->EnzoTimeStep*1e10/1e6,
                        ptcl->CurrentBlockIrr,
                        ptcl->CurrentTimeReg*global_variable->EnzoTimeStep*1e10/1e6,
                        ptcl->CurrentBlockReg,
                        NextRegTimeBlock*global_variable->time_step*global_variable->EnzoTimeStep*1e10/1e6,
                        NextRegTimeBlock,
                        ptcl->TimeStepIrr*global_variable->EnzoTimeStep*1e10/1e6,
                        ptcl->TimeStepReg*global_variable->EnzoTimeStep*1e10/1e6,
                        ptcl->TimeBlockIrr,
                        ptcl->TimeLevelIrr,
                        ptcl->TimeBlockReg,
                        ptcl->TimeLevelReg,
                        ptcl->NextBlockIrr*global_variable->time_step*global_variable->EnzoTimeStep*1e10/1e6,
                        ptcl->NextBlockIrr,
                        ptcl->NumberOfNeighbor
                        );

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
                        ptcl->a_irr[2][3]
                            );
            }
            fflush(stdout);
        }
        */
        current_time_irr = particles[ThisLevelNode->ParticleList[0]].CurrentBlockIrr * global_variable->time_step;
#ifdef DEBUG_ABYSS
        fprintf(nbpout, "skiplist->deleteFirstNode() stacrts\n");
        fflush(nbpout);
#endif
        skiplist->deleteFirstNode();
#ifdef DEBUG_ABYSS
        fprintf(nbpout, "skiplist->deleteFirstNode() ended\n");
        fflush(nbpout);
        fprintf(nbpout, "After deleteFirstNode...\n");
        /*
        for (int i = 0; i < ThisLevelNode->ParticleList.size(); i++) {
            ptcl = &particles[ThisLevelNode->ParticleList[i]];
            fprintf(nbpout, "PID: %d. RadiusOfNeighbor: %e\n", ptcl->PID, ptcl->RadiusOfNeighbor);
        }
        */
#endif

        // std::cout << "deleteFirstNode success" << std::endl;
        // ThisLevelNode = skiplist->getFirstNode();

        /*
        {
            //, NextRegTime= %.3e Myr(%llu),
            for (int i=0; i<ThisLevelNode->ParticleList.size(); i++) {
                ptcl = &particles[ThisLevelNode->ParticleList[i]];
                fprintf(stdout, "PID=%d, CurrentTime (Irr, Reg) = (%.3e(%llu), %.3e(%llu)) Myr, NextReg = %.3e (%llu)\n"\
                        "dtIrr = %.4e Myr, dtReg = %.4e Myr, blockIrr=%llu (%d), blockReg=%llu (%d)\n"\
                        "NumNeighbor= %d\n",
                        ptcl->PID,
                        ptcl->CurrentTimeIrr*global_variable->EnzoTimeStep*1e10/1e6,
                        ptcl->CurrentBlockIrr,
                        ptcl->CurrentTimeReg*global_variable->EnzoTimeStep*1e10/1e6,
                        ptcl->CurrentBlockReg,
                        NextRegTimeBlock*time_step*global_variable->EnzoTimeStep*1e10/1e6,
                        NextRegTimeBlock,
                        ptcl->TimeStepIrr*global_variable->EnzoTimeStep*1e10/1e6,
                        ptcl->TimeStepReg*global_variable->EnzoTimeStep*1e10/1e6,
                        ptcl->TimeBlockIrr,
                        ptcl->TimeLevelIrr,
                        ptcl->TimeBlockReg,
                        ptcl->TimeLevelReg,
                        ptcl->NumberOfNeighbor
                        );
            }
        }
        */
        // exit(SUCCESS);
    } // Irregular Loop ends

#ifdef PerformanceTrace
    start_point = std::chrono::high_resolution_clock::now();
#endif
#ifdef DEBUG_ABYSS
    fprintf(nbpout, "delete skiplist\n");
    fflush(nbpout);
#endif
    delete skiplist;
    skiplist = nullptr;
    // exit(SUCCESS);

#ifdef FEWBODY
    if (bin_termination || new_binaries) {

        if (RegularList.empty())
            return false;
    }
#endif

#ifdef PerformanceTrace
    end_point = std::chrono::high_resolution_clock::now();
    performance.IrregularRoutine +=
        std::chrono::duration_cast<std::chrono::nanoseconds>(end_point - start_point).count();
#endif
    return true;
}


bool createSkipList(SkipList *skiplist) {

	bool debug = false;
	//fprintf(stdout, "create level starts!\n");
	//fflush(stdout);

	if (debug) {
		fprintf(stdout, "create level starts!\n");
		fflush(stdout);
	}

	/*
#ifdef time_trace
	_time.irr_chain.markStart();
#endif
*/

	Particle* ptcl;

	for (int i=0; i<=global_variable->LastParticleIndex; i++) {
		ptcl =  &particles[i];

		// if ((ptcl->NumberOfNeighbor != 0) && (ptcl->NextBlockIrr <= NextRegTimeBlock)) { // IAR original
		if (ptcl->isActive && ptcl->NextBlockIrr <= NextRegTimeBlock) {	// IAR modified
			//fprintf(stdout, "PID=%d, NBI=%llu\n", ptcl->PID, ptcl->NextBlockIrr);
			if (!skiplist->search(ptcl->NextBlockIrr, ptcl->ParticleIndex))
				skiplist->insert(ptcl->NextBlockIrr, ptcl->ParticleIndex);
		}
	}

	/*
#ifdef time_trace
	_time.irr_chain.markEnd();
	_time.irr_chain.getDuration();
#endif
*/

	/*
	if (debug) {
		skiplist->display();
	}
	*/

	//fprintf(stdout, "create level ends!\n");
	//fflush(stdout);

	if (skiplist->getFirstNode() == nullptr)
		return false;
	else
		return true;
}



bool updateSkipList(SkipList *skiplist, int ptcl_id) {
	bool debug = false;
	if (debug) {
		fprintf(stdout, "update level starts!\n");
		fflush(stdout);
	}

	/*
#ifdef time_trace
	_time.irr_sort.markStart();
#endif
*/

	/* Update New Time Steps */
	//Node* ThisLevelNode = skiplist->getFirstNode();

    /*
         if (this->debug) {
         fprintf(stdout, "PID=%d, NBI=%llu, size=%lu\n", ptcl->PID, ptcl->NextBlockIrr, ThisLevelNode->particle_list.size());
         fprintf(stdout, "NextBlockIrr=%llu\n",ptcl->NextBlockIrr);
         fflush(stdout);
         }
         */

    Particle * ptcl = &particles[ptcl_id];

	//std::cout << "NextBlockIrr of "<< ptcl_id<<" = " << ptcl->NextBlockIrr << std::endl;
	if (ptcl->NextBlockIrr > NextRegTimeBlock)
		return true;

	if (!skiplist->search(ptcl->NextBlockIrr, ptcl->ParticleIndex))
		skiplist->insert(ptcl->NextBlockIrr, ptcl->ParticleIndex);

	if (debug) {
	}

	/*
#ifdef time_trace
	_time.irr_sort.markEnd();
	_time.irr_sort.getDuration();
#endif
*/

	if (debug) {
		//skiplist->display();
		//fprintf(stderr, "This is it.\n\n\n");
	}

	return true;
}

