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
int CommunicationToAbyss(LevelHierarchyEntry *LevelArray[], int level,
                         Star *&AllStar,
                         std::unordered_map<int, Star *> &LocalStarLookupMap);
int CommunicationToAbyssInitialize(
    LevelHierarchyEntry *LevelArray[], int level, 
    Star *&AllStars,
    std::unordered_map<int, Star *> &LocalStarLookupMap);

/*
*
*
*
*
*
*/

int SendParticleToAbyss(LevelHierarchyEntry *LevelArray[], int level,
                        Star *&AllStars,
                        std::unordered_map<int, Star *> &LocalStarLookupMap) {

    if (NbodyFirst) {
      if (CommunicationToAbyssInitialize(LevelArray, level, AllStars, LocalStarLookupMap))
        NbodyFirst = FALSE;
    } else {
      CommunicationToAbyss(LevelArray, level, AllStars, LocalStarLookupMap);
    }
  return SUCCESS;
}

/*
*
*
*
*
*
*/


int CommunicationToAbyssInitialize(
    LevelHierarchyEntry *LevelArray[], int level,
    Star *&AllStars, std::unordered_map<int, Star *> &LocalStarLookupMap) {



  fprintf(stdout, "ENZO: CommunicationToAbyssInitialize ...\n");
  fprintf(stderr, "ENZO: CommunicationToAbyssInitialize ...\n");

  /* Do direct calculation!*/
  double dt = 1e-3, scale_factor = 1.0;
  double DensityUnits = 1, LengthUnits = 1, VelocityUnits = 1, TimeUnits = 1,
         TemperatureUnits = 1;
  double MassUnits = 1;
  double Time, TimeStep;
	Time = LevelArray[level]->GridData->ReturnTime(); 
	TimeStep = LevelArray[level]->GridData->ReturnTimeStep(); 
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits,
               &VelocityUnits, &MassUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
  }

  fprintf(stderr, "LengthUnit                = %lf\n", LengthUnits);
  fprintf(stderr, "DensityUnit               = %lf\n", DensityUnits);
  fprintf(stderr, "TimeUnit                  = %lf\n", TimeUnits);
  fprintf(stderr, "VelocityUnit              = %lf\n", VelocityUnits);


  /* Number of Star Particles */
  int LocalNumberOfParticles = LocalStarLookupMap.size();
  fprintf(stderr, "ENZO: (%d) NumberOfParticles=%d\n", MyProcessorNumber,
          LocalNumberOfParticles);



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
    MPI_Send(&ComovingCoordinates, 1, MPI_INT, NumberOfProcessors, 100, inter_comm);
    if (ComovingCoordinates)
      MPI_Send(&CosmologyTableNumberOfBins, 1, MPI_INT, NumberOfProcessors, 150, inter_comm);
    MPI_Isend(params_double, double_size, MPI_DOUBLE, NumberOfProcessors, 200, inter_comm,
              &requests[0]);
    MPI_Isend(params_int, int_size, MPI_INT, NumberOfProcessors, 300, inter_comm, &requests[1]);
    MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);
    delete[] params_double;
    delete[] params_int;
  }


#define no_TEST1
#ifdef TEST1
  fprintf(stderr, "ENZO: in CTABI, ID (%d) = ", MyProcessorNumber);
  for (auto &kv : LocalStarLookupMap) {
    Star *star = kv.second;
    fprintf(stderr, "(1) %d,", star->ReturnID());
    star->GetBackgroundAcceleration();
    fprintf(stderr, "(2) %d,", star->ReturnID());
  }
  fprintf(stderr, "\n");
#endif

  /*-------------------------------------------*/
  /********   Send Particles to ABYSS  *********/
  /*-------------------------------------------*/

  // Step 1: Gather sizes //I can make this MPI_Igather for a slight speep-up.
  MPI_Gather(&LocalNumberOfParticles, 1, MPI_INT, NULL, 1, MPI_INT,
            NumberOfProcessors, inter_comm);

  /* Step 2: Prepare Send Buffer sendbuf */
  ParticleDataType *sendbuf = new ParticleDataType[LocalNumberOfParticles];

  /* Step 3: Copy  Data to Send Buffer */
  int count = 0;
  fprintf(stderr, "ENZO: ID (%d) = ", MyProcessorNumber);
  //fprintf(stderr, "ENZO: Pos of x (%d) = ", MyProcessorNumber);
  for (auto &kv : LocalStarLookupMap) {
    Star *star = kv.second;
    sendbuf[count].copyFrom(star);
    //fprintf(stderr, "(%d, %d, %e)", star->ReturnID(), star->ReturnType(), star->ReturnMass());
    fprintf(stderr, "(%d, ", star->ReturnID());
    //fprintf(stderr, "(%.5e,", star->ReturnPosition()[0]);
    //fprintf(stderr, "%.5e, ),", sendbuf[count].Position[0]);
    count++;
  }
  fprintf(stderr, ")\n");
  fprintf(stderr, "ENZO: Buffer Ready!\n");

  /* Step 4: Gatherv  */
  MPI_Gatherv(sendbuf, LocalNumberOfParticles, MPI_ENZO_PTCL, NULL, NULL, NULL,
              MPI_ENZO_PTCL, NumberOfProcessors, inter_comm);

  fprintf(stdout, "ENZO: data sent to ABYSS (first) \n");
  fprintf(stderr, "ENZO: data sent to ABYSS (first) \n");
  delete [] sendbuf;
#endif
  return SUCCESS;
}

/*
*
*
*
*
*
*/

int CommunicationToAbyss(LevelHierarchyEntry *LevelArray[], int level, Star *&AllStars,
                         std::unordered_map<int, Star *> &LocalStarLookupMap) {

  fprintf(stdout, "ENZO: CommunicationToAbyss ...\n");
  fprintf(stderr, "ENZO: CommunicationToAbyss ...\n");
  fflush(stderr);
  //MPI_Barrier(inter_comm); 


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


  /* Number of Star Particles */
  int LocalNumberOfParticlesOld = 0;
  int LocalNumberOfParticlesNew = 0;
  for (auto &kv : LocalStarLookupMap) {
    Star *star = kv.second;
    if (star->ReturnNewStarFlag())
      LocalNumberOfParticlesNew++;
    else
      LocalNumberOfParticlesOld++;
  }
  fprintf(stderr, "ENZO: (%d) NumberOfParticlesOld=%d/ NewStar=%d\n", MyProcessorNumber,
          LocalNumberOfParticlesOld, LocalNumberOfParticlesNew);



  /*-------------------------------------------*/
  /********   Send Parameters to ABYSS  ********/
  /*-------------------------------------------*/
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    double *params_double = new double[6];
    params_double[0] = TimeStep;
    params_double[1] = TimeUnits;
    params_double[2] = Time;
    params_double[3] = LengthUnits;
    params_double[4] = DensityUnits;
    params_double[5] = VelocityUnits;

    /* Send Parameters First*/
    MPI_Send(params_double, 6, MPI_DOUBLE, NumberOfProcessors, 100, inter_comm);

    delete[] params_double;
  }


#ifdef USE_MPI
  /*---------------------------------------*/
  /******** Send PARTICLES to ABYSS    *****/
  /*---------------------------------------*/

  // Step 1: Gather sizes //I can make this MPI_Igather for a slight speep-up.
  MPI_Gather(&LocalNumberOfParticlesOld, 1, MPI_INT, NULL, 1, MPI_INT,
            NumberOfProcessors, inter_comm);
  MPI_Gather(&LocalNumberOfParticlesNew, 1, MPI_INT, NULL, 1, MPI_INT,
            NumberOfProcessors, inter_comm);


  /* Step 2: Prepare Send Buffer sendbuf */
  ParticleSendDataType *sendbuf_old = new ParticleSendDataType[LocalNumberOfParticlesOld];
  ParticleDataType *sendbuf_new = new ParticleDataType[LocalNumberOfParticlesNew];

  /* Step 3: Copy  Data to Send Buffer */
  int count = 0;
  fprintf(stderr, "ENZO: ID (%d) = ", MyProcessorNumber);
  //fprintf(stderr, "ENZO: Pos of x (%d) = ", MyProcessorNumber);
  for (auto &kv : LocalStarLookupMap) {
    Star *star = kv.second;
    if (star->ReturnNewStarFlag())
      sendbuf_new[count].copyFrom(star);
    else
      sendbuf_old[count].copyFrom(star);
    //fprintf(stderr, "(%d, %d, %e)", star->ReturnID(), star->ReturnType(), star->ReturnMass());
    fprintf(stderr, "(%d, ", star->ReturnID());
    //fprintf(stderr, "(%.5e,", star->ReturnPosition()[0]);
    //fprintf(stderr, "%.5e, ),", sendbuf[count].Position[0]);
    count++;
  }
  fprintf(stderr, ")\n");
  fprintf(stderr, "ENZO: Buffer Ready!\n");

  /* Step 4: Gatherv  */
  MPI_Gatherv(sendbuf_old, LocalNumberOfParticlesOld, MPI_ENZO_PTCL_SEND, NULL, NULL, NULL,
              MPI_ENZO_PTCL_SEND, NumberOfProcessors, inter_comm);
  MPI_Gatherv(sendbuf_new, LocalNumberOfParticlesNew, MPI_ENZO_PTCL, NULL, NULL, NULL,
              MPI_ENZO_PTCL, NumberOfProcessors, inter_comm);

  fprintf(stdout, "ENZO: data sent to ABYSS \n");
  fprintf(stderr, "ENZO: data sent to ABYSS \n");
  delete [] sendbuf_old;
  delete [] sendbuf_new;
#endif

    
  //MPI_Barrier(inter_comm); 
  return SUCCESS;
}

/*
*
*
*
*
*
*/

int ReceiveParticleFromAbyss(
    Star *&AllStars, std::unordered_map<int, Star *> &LocalStarLookupMap) {

  fprintf(stderr, "ENZO: Starting ReceiveParticleFromAbyss...\n");
          
  /* Number of Star Particles */
  int LocalNumberOfParticles = LocalStarLookupMap.size();
  fprintf(stderr, "ENZO: (%d) NumberOfParticles=%d\n", MyProcessorNumber,
          LocalNumberOfParticles);

  //MPI_Barrier(inter_comm); 

#ifdef USE_MPI
  /*------------------------------------------*/
  /******** Receive PARTICLES to ABYSS    *****/
  /*------------------------------------------*/

  /* Prepare Send Buffer */
  ParticleReceiveDataType *recvbuf = new ParticleReceiveDataType[LocalNumberOfParticles];
      
  fprintf(stdout, "ENZO: Waiting for ABYSS, to receive data \n");
  fprintf(stderr, "ENZO: Waiting for ABYSS, to receive data \n");
  fflush(stderr);

  MPI_Scatterv(NULL, NULL, NULL, MPI_ENZO_PTCL_RECV,
          recvbuf, LocalNumberOfParticles, MPI_ENZO_PTCL_RECV, NumberOfProcessors, inter_comm);
  fprintf(stderr, "ENZO: data received! \n");



  // Synchronize Cluster Posiition
  MPI_Request request;
  if (isNbodyParticleIdentification && isIdentificationOnTheFly) {
    MPI_Ibcast(NbodyClusterPosition, 3, MPI_DOUBLE, NumberOfProcessors, inter_comm,
              &request);
  }

  fprintf(stderr, "ENZO: Data received from ABYSS\n");

  /* Copy Data to Star */
  int count = 0;
	for (int i=0; i<LocalNumberOfParticles; i++) {
		auto it  = LocalStarLookupMap.find(recvbuf[i].ID);
    if (it != LocalStarLookupMap.end()) {
      recvbuf[i].copyTo(it->second);
      fprintf(stderr, "ENZO: %d successfully copied on %d.\n",
      recvbuf[i].ID, MyProcessorNumber);
    }
    else {
      fprintf(stderr, "ENZO: LocalStarLookupMap doesn't have %d on %d!\n",
      recvbuf[i].ID, MyProcessorNumber);
      exit(0);
    }
	}
  fprintf(stdout, "ENZO: Copy-to-Star DONE!\n");
  fprintf(stderr, "ENZO: Copy-to-Star DONE %d!\n", MyProcessorNumber);
  delete[] recvbuf;

  if (isNbodyParticleIdentification && isIdentificationOnTheFly) {
    MPI_Wait(&request, MPI_STATUS_IGNORE);
    fprintf(stderr, "In Final, NbodyClusterPosition = (%e, %e, %e)\n",
            NbodyClusterPosition[0], NbodyClusterPosition[1],
            NbodyClusterPosition[2]);
  }
#endif

  //MPI_Barrier(inter_comm); 
  std::cerr << "ENZO: Receiving data done!" << std::endl;
  return SUCCESS;
}
#endif

/*
void Star::copyFrom(ParticleReceiveDataType *ptcl) {
}*/
