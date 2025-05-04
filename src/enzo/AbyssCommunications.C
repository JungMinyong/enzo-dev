/***********************************************************************
        /
        /  FIND ALL NBODY PARTICLES OVER ALL PROCESSORS
        /
        /  written by: Yongseok Jo
        /  date:       April, 2025
        /
        /  PURPOSE: First synchronizes particle information in the normal and
        /           abyss particles.
        /
 ************************************************************************/


#if defined (NBODY) && defined (INDIVIDUALSTAR)
#include <unordered_map>
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "CosmologyParameters.h"
#include "phys_constants.h"
#define ENZO_ONLY
#include "NbodyRoutines.h" //added
#include "Star.h"
//#undef max
//#undef min
//#include "abyss/global.h"

int GetUnits(double *DensityUnits, double *LengthUnits,
             double *TemperatureUnits, double *TimeUnits, double *VelocityUnits,
             double *MassUnits, double Time);



#ifdef TEST



int CommunicationToAbyss(LevelHierarchyEntry *LevelArray[], int level,
                         Star *&AllStar,
                         std::unordered_map<int, Star *> LocalStarLookupMap);

#endif

int CommunicationToAbyssInitialize(
   // LevelHierarchyEntry *LevelArray[], int level, 
    Star *&AllStars,
    std::unordered_map<int, Star *> &LocalStarLookupMap);


int SendParticleToAbyss(//LevelHierarchyEntry *LevelArray[], int level,
                        Star *&AllStars,
                        std::unordered_map<int, Star *> &LocalStarLookupMap) {

    if (NbodyFirst) {
      if (CommunicationToAbyssInitialize(AllStars, LocalStarLookupMap))
        NbodyFirst = FALSE;
    } else {
      //CommunicationToAbyss(LevelArray, level, AllStars, LocalStarLookupMap);
    }
  return SUCCESS;
}

int CommunicationToAbyssInitialize(
    // LevelHierarchyEntry *LevelArray[], int level,
    Star *&AllStars, std::unordered_map<int, Star *> &LocalStarLookupMap) {

  fprintf(stdout, "ENZO: CommunicationToAbyssInitialize ...\n");
  fprintf(stderr, "ENZO: CommunicationToAbyssInitialize ...\n");

  /* Do direct calculation!*/
  double dt = 1e-3, scale_factor = 1.0;
  double DensityUnits = 1, LengthUnits = 1, VelocityUnits = 1, TimeUnits = 1,
         TemperatureUnits = 1;
  double MassUnits = 1;
  double Time, TimeStep;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits,
               &VelocityUnits, &MassUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
  }

  Star *ThisStar;

  /* Number of Star Particles */
  int NumberOfParticles = LocalStarLookupMap.size();
  fprintf(stderr, "ENZO: (%d) NumberOfParticles=%d\n", MyProcessorNumber,
          NumberOfParticles);



#ifdef USE_MPI
  MPI_Request requests[2];

  /*-------------------------------------------*/
  /********   Send Parameters to ABYSS  ********/
  /*-------------------------------------------*/
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    /* Prepare buffer array */
    int double_size = 13;
    int int_size = 6;
    if (ComovingCoordinates) {
      double_size += 14;
      double_size += 2 * CosmologyTableNumberOfBins;
      int_size += 1;
    }
    double *params_double = new double[double_size];
    int *params_int = new int[int_size];
    params_double[0] = TimeStep;
    params_double[1] = TimeUnits;
    params_double[2] = LengthUnits;
    params_double[3] = DensityUnits;
    params_double[4] = VelocityUnits;
    params_double[5] = StarMassEjectionFraction;
    params_double[6] = Time;
    params_double[7] = NbodySmoothingLength;
    params_double[8] = NbodyTimeStepConstant;
    params_double[9] = NbodyNeighborRadius;
    params_double[10] = NbodyClusterPosition[0];
    params_double[11] = NbodyClusterPosition[1];
    params_double[12] = NbodyClusterPosition[2];
    params_double[13] = NbodyClusterPosition[3];

    params_int[0] = StarParticleFeedback;
    params_int[1] = isNbodyParticleIdentification;
    params_int[2] = isIdentificationOnTheFly;
    params_int[3] = NbodyFixNumNeighbor;
    params_int[4] = NbodyMaxNumNeighbor;
    params_int[5] = NbodyBinaryRegularization;

    if (ComovingCoordinates) {
      params_double[14] = HubbleConstantNow;
      params_double[15] = OmegaMatterNow;
      params_double[16] = OmegaDarkMatterNow;
      params_double[17] = OmegaLambdaNow;
      params_double[18] = OmegaRadiationNow;
      params_double[19] = ComovingBoxSize;
      params_double[20] = MaxExpansionRate;
      params_double[21] = InitialTimeInCodeUnits;
      params_double[22] = InitialRedshift;
      params_double[23] = FinalRedshift;
      params_double[24] = CosmologyTableLogaInitial;
      params_double[25] = CosmologyTableLogaFinal;
      for (int i = 0; i < CosmologyTableNumberOfBins; i++) {
        params_double[26 + i] = CosmologyTableLoga[i];
      }
      for (int i = 0; i < CosmologyTableNumberOfBins; i++) {
        params_double[26 + CosmologyTableNumberOfBins + i] =
            CosmologyTableLogt[i];
      }

      params_int[7] = CosmologyTableLogtIndex;
    }

    /* Send Parameters First*/
    fprintf(stderr, "ENZO: ComovingCoordinates=%d\n", ComovingCoordinates);
    MPI_Send(&ComovingCoordinates, 1, MPI_INT, 1, 100, inter_comm);
    if (ComovingCoordinates)
      MPI_Send(&CosmologyTableNumberOfBins, 1, MPI_INT, 1, 150, inter_comm);
    MPI_Isend(params_double, double_size, MPI_DOUBLE, 1, 200, inter_comm,
              &requests[0]);
    MPI_Isend(params_int, int_size, MPI_INT, 1, 300, inter_comm, &requests[1]);
    if (MyProcessorNumber == ROOT_PROCESSOR)
      MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);
    delete[] params_double;
    delete[] params_int;
  }

  /*-------------------------------------------*/
  /********   Send Particles to ABYSS  *********/
  /*-------------------------------------------*/

  /* Prepare Send Buffer Packet */
  ParticleDataType *packet = new ParticleDataType[NumberOfParticles];

  /* Copy  Data to Send Buffer */
  int count = 0;
  for (auto &kv : LocalStarLookupMap) {
    Star *star = kv.second;
    packet[count++].copyFrom(star);
  }
  fprintf(stdout, "ENZO: Buffer Ready!\n");


  //MPI_Send(&NumberOfParticles, 1, MPI_INT, ReceiverRank, 100, MPI_COMM_WORLD);
  //fprintf(stderr, "ENZO: done %d\n", MyProcessorNumber);
  //MPI_Barrier(MPI_COMM_WORLD);
  //MPI_Send(packet, NumberOfParticles, MPI_ENZO_PTCL, ReceiverRank, 200, MPI_COMM_WORLD);

  MPI_Isend(&NumberOfParticles, 1, MPI_INT, 1, 100, inter_comm, &requests[0]);
  MPI_Isend(packet, NumberOfParticles, MPI_ENZO_PTCL, 1, 200, inter_comm, &requests[1]);

  fprintf(stdout, "ENZO: Waiting for ABYSS to send data (first) \n");
  fprintf(stderr, "ENZO: Waiting for ABYSS to send data (first) \n");
  MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);
  delete [] packet;
#endif

  return SUCCESS;
}

#ifdef TEST
int CommunicationToAbyss(Star *&AllStars,
                         std::unordered_map<int, Star *> &LocalStarLookupMap) {

  /* Do direct calculation!*/
  double dt = 1e-3, scale_factor = 1.0;
  double DensityUnits = 1, LengthUnits = 1, VelocityUnits = 1, TimeUnits = 1,
         TemperatureUnits = 1;
  double MassUnits = 1;
  double Time, TimeStep;
  Time = LevelArray[level]->GridData->ReturnTime();         // Not sure ?
  TimeStep = LevelArray[level]->GridData->ReturnTimeStep(); // Not sure ?
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits,
               &VelocityUnits, &MassUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
  }
  fprintf(stdout, "TimeStep: %e, %f\n", TimeStep, TimeUnits);

#ifdef USE_MPI
  /*---------------------------------------*/
  /******** Send PARTICLES to ABYSS    *****/
  /*---------------------------------------*/
  int OldNumberOfParticles = 0, NewNumberOfParticles = 0;

  /* Number of Star Particles */
  for (auto &kv : LocalStarLookupMap) {
    Star *star = kv.second;
    if (star->ReturnNewStarFlag())
      OldNumberOfParticles++;
    else
      NewNumberOfParticles++;
  }
  fprintf(stderr, "ENZO: (%d) OldNumberOfParticles=%d\n", MyProcessorNumber,
          OldNumberOfParticles);
  fprintf(stderr, "ENZO: (%d) NewNumberOfParticles=%d\n", MyProcessorNumber,
          NewNumberOfParticles);

  /* Prepare Send Buffer */
  ParticleSendDataType *ptcl_old =
      new ParticleSendDataType[OldNumberOfParticles];
  ParticleDataType *ptcl_new = new ParticleDataType[NewNumberOfParticles];

  /* Copy  Data to Send Buffer */
  int count_new = 0, count_old = 0;
  for (auto &kv : LocalStarLookupMap) {
    Star *star = kv.second;
    // if (ThisStar->isABYSS) this is not implemented yet.
    if (star->ReturnNewStarFlag())
      ptcl_new[count_new++].copyFrom(star);
    else
      ptcl_old[count_old++].copyFrom(star);
  }
  fprintf(stdout, "ENZO: Buffer Ready!\n");

  fprintf(stdout, "ENZO: Waiting for ABYSS to send data \n");
  fprintf(stderr, "ENZO: Waiting for ABYSS to send data \n");
  int ReceiverRank =
      NumberOfProcessors + MyProcessorNumber % (NumberOfAbyssProcessors);

  MPI_Request request;
  MPI_Send(&OldNumberOfParticles, 1, MPI_INT, ReceiverRank, 100, inter_comm);
  if (OldNumberOfParticles != 0)
    MPI_Isend(ptcl_old, OldNumberOfParticles, MPI_ENZO_PTCL_SEND, ReceiverRank,
              200, inter_comm, &request);
  MPI_Wait(&request, MPI_STATUSES_IGNORE);
  MPI_Send(&NewNumberOfParticles, 1, MPI_INT, ReceiverRank, 300, inter_comm);
  if (NewNumberOfParticles != 0)
    MPI_Isend(ptcl_new, NewNumberOfParticles, MPI_ENZO_PTCL, ReceiverRank, 400,
              inter_comm, &request);
  MPI_Wait(&request, MPI_STATUSES_IGNORE);
  fprintf(stderr, "ENZO: data sent! \n");

  delete[] ptcl_new;
  delete[] ptcl_old;

#endif

  return SUCCESS;
}

int ReceiveParticleFromAbyss(
    LevelHierarchyEntry *LevelArray[], int level, Star *&AllStars,
    std::unordered_map<int, Star *> &LocalStarLookupMap) {
  Star *ThisStar;
  if (LevelArray[level + 1] != NULL) {
    return SUCCESS;
  }

  /* Number of Star Particles */
  int NumberOfParticles = LocalStarLookupMap.size();


  fprintf(stderr, "ENZO: (%d) NumberOfParticles=%d\n", MyProcessorNumber,
          NumberOfParticles);

#ifdef USE_MPI
  MPI_Status statuses[2];
  MPI_Request requests[2];
  // let's use broadcast
  if (isNbodyParticleIdentification && isIdentificationOnTheFly) {
    MPI_Irecv(NbodyClusterPosition, 3, MPI_DOUBLE, 1, 200, inter_comm,
              &requests[0]);
  }
  fprintf(stdout, "In Final, NbodyClusterPosition = (%e, %e, %e)\n",
          NbodyClusterPosition[0], NbodyClusterPosition[1],
          NbodyClusterPosition[2]);
  /*
  if (isNbodyParticleIdentification && isIdentificationOnTheFly)  {
          MPI_Bcast(NbodyClusterPosition, 3, MPI_DOUBLE, ROOT_PROCESSOR,
  enzo_comm);
          //fprintf(stderr, "In Final all, NbodyClusterPosition = (%e, %e,
  %e)\n", NbodyClusterPosition[0], NbodyClusterPosition[1],
  NbodyClusterPosition[2]);
  }*/
#endif

  if (NumberOfParticles == 0)
    return SUCCESS;

#ifdef USE_MPI
  /*------------------------------------------*/
  /******** Receive PARTICLES to ABYSS    *****/
  /*------------------------------------------*/

  /* Prepare Send Buffer */
  ParticleReceiveDataType *ptcl =
      new ParticleReceiveDataType[NumberOfParticles];

  fprintf(stdout, "ENZO: Waiting for ABYSS to receive data \n");
  fprintf(stderr, "ENZO: Waiting for ABYSS to receive data \n");
  int ReceiverRank =
      NumberOfProcessors + MyProcessorNumber % (NumberOfAbyssProcessors);
  MPI_Irecv(ptcl, NumberOfParticles, MPI_ENZO_PTCL_RECV, ReceiverRank, 100,
            inter_comm, &requests[1]);

  MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);

  // think of broadcast

  fprintf(stderr, "ENZO: data sent! \n");

  /* Copy Data to Star */
  int count = 0;
	for (int i=0; i<NumberOfParticles; i++) {
		Star *star = LocalStarLookupMap.at(ptcl[i].ID);
		star->copyFrom(&ptcl[i]);
	}
  for (ThisStar = AllStars; ThisStar; ThisStar = ThisStar->NextStar) {
    ptcl[count++].copyTo(ThisStar);
  }
  fprintf(stdout, "ENZO: Copy to Star!\n");

  delete[] ptcl;
#endif
  return SUCCESS;
}
#endif

/*
void Star::copyFrom(ParticleReceiveDataType *ptcl) {
}*/


#endif
