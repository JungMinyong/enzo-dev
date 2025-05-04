#include <cstdio>
#ifdef INDIVIDUALSTAR
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#define MAX_NUMBER_OF_OUTPUT_REDSHIFTS 500
#define FLOAT double
#include "../CosmologyParameters.h"
#undef FLOAT
//#include "../phys_constants.h"

#ifdef FEWBODY
#include <unordered_set>
#endif
#include "../NbodyRoutines.h"

#include "Worker.h"
#undef NormalStar
#include "def.h"
#include "global.h"
#define NO_COM_EVOLUTION // (Query) eventaully I think this should adopt COM
                         // evolution due to bulk motion.
// Bulk motion in irregular force will cause some distorts since a fraction of
// particle will advance due to bulk motion
//  and that will cause artificial tidal force.

//

extern int StarParticleFeedback;
extern double StarMassEjectionFraction;
extern int ComovingCoordinates;
extern MPI_Datatype MPI_ENZO_PTCL;
extern MPI_Datatype MPI_ENZO_PTCL_SEND;
extern MPI_Datatype MPI_ENZO_PTCL_RECV;
extern int NumberOfProcessors;
extern int TotalNumberOfProcessors;
double EnzoCurrentTime, ClusterRadius2;
extern int WorldProcessorNumber;
double ClusterAcceleration[Dim], ClusterPosition[Dim], ClusterVelocity[Dim],
    EnzoClusterPosition[Dim + 1];
double eta_tmp;
int FixNumNeighbor0, IdentifyNbodyParticles;
int BinaryRegularization, IdentifyOnTheFly;

// double KSTime;
// double KSDistance;

void deleteParticle(int &PID, int &index);
// int writeParticle(std::vector<Particle*> &particle, double MinRegTime, int
// outputNum); void InitializeParticle(std::vector<Particle*> &particle); void
// InitializeNewParticle(std::vector<Particle*> &particle, int offset, int
// newSize);
void GetCenterOfMass(double *mass, double *x[Dim], double *v[Dim],
                    double x_com[], double v_com[], int N);
void GetNewCenterOfMass(int *PID, double *mass2, double *x2[Dim],
                        double *v2[Dim], int n2, double x_X[], double v_X[]);
// void UpdateNextRegTime(std::vector<Particle*> &particle);
int CommunicationInterBarrier();
void broadcastFromRoot(int &data);
void broadcastFromRoot(double &data);
void broadcastFromRoot(ULL &data);
// void CalculateAllAccelerationOnGPU(std::vector<Particle*> &particle);
void InitializationAfterCommunication();

#ifdef SEVN
void initializeStellarEvolution();
void initializeStellarEvolution(int ParticleIndex);
#endif

using namespace std;
const int width = 18;

// int InitialCommunication(std::vector<Particle*> &particle) {
int InitialCommunication() {

  fprintf(nbpout, "ABYSS: First Waiting for Enzo to receive data...\n");
  fprintf(stderr, "ABYSS: First Waiting for Enzo to receive data...\n");

  /*-------------------------------------------*/
  /****   Receive Parameters from ENZO  ********/
  /*-------------------------------------------*/
  MPI_Request requests[2];
  MPI_Status statuses[2];
  MPI_Recv(&ComovingCoordinates, 1, MPI_INT, 0, 100, inter_comm, &statuses[0]);
  fprintf(stderr, "ABYSS: ComovingCoordinates=%d\n", ComovingCoordinates);
  if (ComovingCoordinates)
    MPI_Recv(&CosmologyTableNumberOfBins, 1, MPI_INT, 0, 150, inter_comm,
             &statuses[1]);

  /* Prepare buffer array */
  double *CosmologyTableLoga;
  double *CosmologyTableLogt;
  double TimeStep, TimeUnits, LengthUnits, VelocityUnits, DensityUnits;
  int double_size = 13;
  int int_size = 6;
  if (ComovingCoordinates) {
    double_size += 14;
    double_size += 2 * CosmologyTableNumberOfBins;
    int_size += 1;
    CosmologyTableLoga =
        new double[CosmologyTableNumberOfBins]; // (Query) EW: not
                                                // deleted later?
    CosmologyTableLogt =
        new double[CosmologyTableNumberOfBins]; // (Query) EW: not
                                                // deleted later?
  }
  double *params_double = new double[double_size];
  int *params_int = new int[int_size];

  MPI_Irecv(params_double, double_size, MPI_DOUBLE, 0, 200, inter_comm,
            &requests[0]);
  MPI_Irecv(params_int, int_size, MPI_INT, 0, 300, inter_comm, &requests[1]);
  MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);
  fprintf(nbpout, "Data received!\n");

  TimeStep = params_double[0];
  TimeUnits = params_double[1];
  LengthUnits = params_double[2];
  DensityUnits = params_double[3];
  VelocityUnits = params_double[4];
  StarMassEjectionFraction = params_double[5];
  EnzoCurrentTime = params_double[6];
  EPS2 = params_double[7];
  eta_tmp = params_double[8];
  InitialNeighborRadius2 = params_double[9];
  EnzoClusterPosition[0] = params_double[10];
  EnzoClusterPosition[1] = params_double[11];
  EnzoClusterPosition[2] = params_double[12];
  EnzoClusterPosition[3] = params_double[13];

  StarParticleFeedback = params_int[0];
  IdentifyNbodyParticles = params_int[1];
  IdentifyOnTheFly = params_int[2];
  FixNumNeighbor = params_int[3];
  MaxNumNeighbor = params_int[4];
  BinaryRegularization = params_int[5];

  if (ComovingCoordinates) {
    HubbleConstantNow = params_double[14];
    OmegaMatterNow = params_double[15];
    OmegaDarkMatterNow = params_double[16];
    OmegaLambdaNow = params_double[17];
    OmegaRadiationNow = params_double[18];
    ComovingBoxSize = params_double[19];
    MaxExpansionRate = params_double[20];
    InitialTimeInCodeUnits = params_double[21];
    InitialRedshift = params_double[22];
    FinalRedshift = params_double[23];
    CosmologyTableLogaInitial = params_double[24];
    CosmologyTableLogaFinal = params_double[25];
    for (int i = 0; i < CosmologyTableNumberOfBins; i++) {
      CosmologyTableLoga[i] = params_double[26 + i];
    }
    for (int i = 0; i < CosmologyTableNumberOfBins; i++) {
      CosmologyTableLogt[i] =
          params_double[26 + CosmologyTableNumberOfBins + i];
    }
    CosmologyTableLogtIndex = params_int[7];
  }
  ClusterRadius2 = EnzoClusterPosition[3]; // it's already squared

  // Enzo to Nbody unit convertors
  // EnzoMass         = MassUnits/Msun/mass_unit;
  EnzoMass = DensityUnits * pow(LengthUnits, 3.) / Msun / mass_unit;
  EnzoLength = LengthUnits / pc / position_unit;
  EnzoVelocity = VelocityUnits / pc * yr / velocity_unit;
  EnzoTime = TimeUnits / yr / time_unit;
  // EnzoAcceleration =
  // LengthUnits/TimeUnits/TimeUnits/pc*yr*yr/position_unit*time_unit*time_unit;
  EnzoAcceleration = EnzoLength / EnzoTime / EnzoTime;

  // Unit conversion
  global_variable->EnzoTimeStep = TimeStep * EnzoTime;

  EnzoCurrentTime *= EnzoTime;

  if (EPS2 < 0)
    EPS2 = -1;
  else {
    EPS2 *= EnzoLength;
    EPS2 *= EPS2;
  }

  InitialNeighborRadius2 *= EnzoLength;
  InitialNeighborRadius2 *= InitialNeighborRadius2;
  FixNumNeighbor0 = FixNumNeighbor;

  // need to fix units
  // fprintf(nbpout, "Enzo Time                = %lf\n", TimeStep);
  fprintf(nbpout, "Nbody Time               = %lf\n", EnzoCurrentTime);
  fprintf(nbpout, "Nbody TimeStep           = %lf\n",
          global_variable->EnzoTimeStep);
  fprintf(nbpout, "EPS2                     = %lf pc**2\n",
          EPS2 * position_unit * position_unit);
  fprintf(nbpout, "InitialNeighborRadius2        = %.2e pc**2\n",
          InitialNeighborRadius2 * position_unit * position_unit);
  fprintf(nbpout, "eta                      = %lf\n", eta);
  fprintf(nbpout, "ClusterRadius2           = %.2e pc**2\n",
          ClusterRadius2 * position_unit * position_unit);
  fprintf(nbpout, "StarMassEjectionFraction = %lf\n", StarMassEjectionFraction);
  fprintf(nbpout, "StarParticleFeedback     = %d\n", StarParticleFeedback);
  fprintf(nbpout, "FixNumNeighbor           = %d\n", FixNumNeighbor);
  fprintf(nbpout, "BinaryRegularization     = %d\n",
          BinaryRegularization); // (Query) EW: What is this?
  // fprintf(nbpout, "KSTime                   = %lf\n", KSTime);
  // fprintf(nbpout, "KSDistance               = %lf\n", KSDistance);
  fprintf(nbpout, "IdentifyNbodyParticles   = %d\n", IdentifyNbodyParticles);
  fprintf(nbpout, "IdentifyOnTheFly         = %d\n\n", IdentifyOnTheFly);
  fflush(nbpout);

  /*------------------===-------------------------*/
  /********   Receive Particles to ABYSS  *********/
  /*----------------===---------------------------*/
  // int LocalNumberOfParticles=0;

  MPI_Recv(&NumberOfSingleParticle, 1, MPI_INT, 0, 100, inter_comm,
           MPI_STATUS_IGNORE);
  ParticleDataType *packet = new ParticleDataType[NumberOfSingleParticle];
  MPI_Recv(packet, NumberOfSingleParticle, MPI_ENZO_PTCL, 0, 200, inter_comm,
           MPI_STATUS_IGNORE);

  /* Unpack the Packet */
  for (int dim = 0; dim < Dim; dim++) {
    ClusterAcceleration[dim] = 0;
    ClusterPosition[dim] = 0;
    ClusterVelocity[dim] = 0;
  }

  /*
  if (NumberOfSingleParticle != 0) {
    // set COM for background acc
    GetCenterOfMass(Mass, Position, Velocity, ClusterPosition, ClusterVelocity,
                    NumberOfSingleParticle);

#ifdef COM_EVOLUTION
    ClusterAcceleration[0] = 0.;
    ClusterAcceleration[1] = 0.;
    ClusterAcceleration[2] = 0.;
    double total_mass = 0.;
    for (int i = 0; i < NumberOfSingleParticle; i++) {
      for (int dim = 0; dim < Dim; dim++) {
        ClusterAcceleration[dim] += Mass[i] * BackgroundAcceleration[dim][i];
      }
      total_mass += Mass[i];
    }

    for (int dim = 0; dim < Dim; dim++)
      ClusterAcceleration[dim] /= total_mass;
#endif

    for (int i = 0; i < LocalNumberOfParticles; i++) {
      for (int dim = 0; dim < Dim; dim++) {
#ifdef COM_EVOLUTION
        BackgroundAcceleration[dim][i] -= ClusterAcceleration[dim];
        Velocity[dim][i] -= ClusterVelocity[dim];
#endif
        // Position[dim][i]               -= ClusterPosition[dim];
      }
      EnzoPIDs[i] = packet[i].ID;
      particles[i].set(&packet[i]);
      PIDtoIndexMap.insert({packet[i].ID, i});
      particles[i].ParticleIndex = i;
    }
    delete [] packet;

    global_variable->EnzoCurrentTime =
        EnzoCurrentTime *
        1e4; // This should be 0.0 // Unlike original value, this is Myr unit
    std::cerr << "EnzoCurrentTime:" << global_variable->EnzoCurrentTime
              << std::endl;
#ifdef SEVN
    initializeStellarEvolution();
#endif
  }
     */

  // NumberOfParticle = NumberOfSingleParticle;
  // global_variable->LastParticleIndex = NumberOfSingleParticle - 1;

  fprintf(nbpout, "ABYSS: %d particles loaded on %d!\n", NumberOfSingleParticle,
          AbyssProcessorNumber);
  fflush(nbpout);
  // exit(0);
  return true;
}

#ifdef TEST
int ReceiveParticleFromEnzo() {

  int i;
  int *PID, *newPID;
  double *BackgroundAcceleration[Dim];
  double *Mass, *newMass, *newPosition[Dim], *newVelocity[Dim],
      *newBackgroundAcceleration[Dim];
  double *newCreationTime, *newDynamicalTime, *newMetallicity;
  double TimeStep;

  MPI_Request request;
  MPI_Status status;
  fprintf(nbpout, "ABYSS: Waiting for Enzo to receive data...\n");

  fprintf(stdout, "ABYSS: Waiting for Enzo to receive data...\n");
  CommunicationInterBarrier();
  fprintf(nbpout, "ABYSS: Receiving data from Enzo...\n");
  fprintf(stdout, "ABYSS: Receiving data from Enzo...\n");

  // Existing Particle Information

  int EnzoNumberOfSingleParticle;
  MPI_Recv(&EnzoNumberOfSingleParticle, 1, MPI_INT, 0, 10, inter_comm, &status);
  if (EnzoNumberOfSingleParticle != NumberOfSingleParticle) {
    fprintf(stderr,
            "The numbers of nbody and hydro do not match: "
            "NumberOfSingleParticle=%d, EnzoNumberOfSingleParticle=%d\n",
            NumberOfSingleParticle, EnzoNumberOfSingleParticle);
    fflush(stderr);
    throw runtime_error("");
  }
  if (NumberOfSingleParticle != 0) {
    PID = new int[NumberOfSingleParticle];
    Mass = new double[NumberOfSingleParticle];
    MPI_Recv(PID, NumberOfSingleParticle, MPI_INT, 0, 25, inter_comm, &status);
    MPI_Recv(Mass, NumberOfSingleParticle, MPI_DOUBLE, 0, 40, inter_comm,
             &status);
    for (int dim = 0; dim < Dim; dim++) {
      BackgroundAcceleration[dim] = new double[NumberOfSingleParticle];
      MPI_Recv(BackgroundAcceleration[dim], NumberOfSingleParticle, MPI_DOUBLE,
               0, 50, inter_comm, &status);
    }
  }

  // New Particle Information
  MPI_Recv(&newNumberOfSingleParticle, 1, MPI_INT, 0, 100, inter_comm, &status);
  if (newNumberOfSingleParticle > 0) {
    newPID = new int[newNumberOfSingleParticle];
    newMass = new double[newNumberOfSingleParticle];
    newCreationTime = new double[newNumberOfSingleParticle];
    newDynamicalTime = new double[newNumberOfSingleParticle];
    newMetallicity = new double[newNumberOfSingleParticle];
    MPI_Recv(newPID, newNumberOfSingleParticle, MPI_INT, 0, 200, inter_comm,
             &status);
    MPI_Recv(newMass, newNumberOfSingleParticle, MPI_DOUBLE, 0, 201, inter_comm,
             &status);
    MPI_Recv(newCreationTime, newNumberOfSingleParticle, MPI_DOUBLE, 0, 202,
             inter_comm, &status);
    MPI_Recv(newDynamicalTime, newNumberOfSingleParticle, MPI_DOUBLE, 0, 203,
             inter_comm, &status);
    MPI_Recv(newMetallicity, newNumberOfSingleParticle, MPI_DOUBLE, 0, 204,
             inter_comm, &status);

    for (int dim = 0; dim < Dim; dim++) {
      newPosition[dim] = new double[newNumberOfSingleParticle];
      newVelocity[dim] = new double[newNumberOfSingleParticle];
      newBackgroundAcceleration[dim] = new double[newNumberOfSingleParticle];
      MPI_Recv(newPosition[dim], newNumberOfSingleParticle, MPI_DOUBLE, 0, 300,
               inter_comm, &status);
      MPI_Recv(newVelocity[dim], newNumberOfSingleParticle, MPI_DOUBLE, 0, 400,
               inter_comm, &status);
      MPI_Recv(newBackgroundAcceleration[dim], newNumberOfSingleParticle,
               MPI_DOUBLE, 0, 500, inter_comm, &status);
    }
  }

  // Timestep
  MPI_Recv(&TimeStep, 1, MPI_DOUBLE, 0, 600, inter_comm, &status);
  double OldEnzoCurrentTime = EnzoCurrentTime;
  MPI_Recv(&EnzoCurrentTime, 1, MPI_DOUBLE, 0, 700, inter_comm, &status);
  // for cosmological runs, we need to update the units
  double TimeUnits, LengthUnits, VelocityUnits, DensityUnits; // MassUnits;
  MPI_Recv(&TimeUnits, 1, MPI_DOUBLE, 0, 701, inter_comm, &status);
  MPI_Recv(&LengthUnits, 1, MPI_DOUBLE, 0, 800, inter_comm, &status);
  MPI_Recv(&DensityUnits, 1, MPI_DOUBLE, 0, 900, inter_comm, &status);
  MPI_Recv(&VelocityUnits, 1, MPI_DOUBLE, 0, 1000, inter_comm, &status);
  CommunicationInterBarrier();

  // for cosmological runs
  EnzoMass = DensityUnits * pow(LengthUnits, 3.) / Msun / mass_unit;
  EnzoLength = LengthUnits / pc / position_unit;
  EnzoVelocity = VelocityUnits / pc * yr / velocity_unit;
  EnzoTime = TimeUnits / yr / time_unit;
  // EnzoAcceleration =
  // LengthUnits/TimeUnits/TimeUnits/pc*yr*yr/position_unit*time_unit*time_unit;
  EnzoAcceleration = EnzoLength / EnzoTime / EnzoTime;

  // std::cout << "Enzo  Time    :" << EnzoCurrentTime << std::endl;
  // std::cout << "Nbody Time    :" <<
  // OldEnzoCurrentTime+particle[0]->CurrentTimeReg*EnzoTimeStep << std::endl;
  global_variable->OldEnzoTimeStep = global_variable->EnzoTimeStep;
  global_variable->EnzoTimeStep = TimeStep * EnzoTime;

  std::cout << "ABYSS: Data transferred!" << std::endl;
  fprintf(stdout, "ABYSS: Data transferred!\n");
  EnzoCurrentTime = EnzoCurrentTime * EnzoTime;
  global_variable->EnzoCurrentTime = EnzoCurrentTime * 1e4; // in Myr unit

  std::cout << "Enzo Time    :" << global_variable->EnzoCurrentTime << " Myr"
            << std::endl;
  std::cerr << "Enzo Time    :" << global_variable->EnzoCurrentTime << " Myr"
            << std::endl;

  std::cout << "EnzoTimeStep :" << global_variable->EnzoTimeStep * 1e4 << " Myr"
            << std::endl;
  std::cerr << "EnzoTimeStep :" << global_variable->EnzoTimeStep * 1e4 << " Myr"
            << std::endl;
  std::cerr << "Next EnzoTimeStep :" << global_variable->EnzoTimeStep * 1e4
            << " Myr" << std::endl;
  // std::cout << "Nbody Mass    :" << particle[0]->Mass << std::endl;
  // std::cout << "Enzo  Mass    :" << Mass[0]*EnzoMass << std::endl;
  // std::cout << "enzo Time :" << TimeStep << std::endl;
  // std::cout << "nbody Time:" << EnzoTimeStep << std::endl;
  //  FixNumNeighbor = std::min((int)
  //  std::floor((NumberOfSingleParticle+newNumberOfSingleParticle-1)/2.0),
  //  FixNumNeighbor0); // original by EW 2025.3.27

  /*
          if (NumberOfSingleParticle + newNumberOfSingleParticle < 1000)
                  FixNumNeighbor = 30;
          else
                  FixNumNeighbor =
     static_cast<int>(sqrt(NumberOfSingleParticle+newNumberOfSingleParticle));
     // test by EW 2025.3.27
  */

  // COM conversion
  // later on we might need to take mass weight into account
  // 1. F=ma; 2. F -> F_com; 3. F_com -> a_com
  ClusterAcceleration[0] = 0.0;
  ClusterAcceleration[1] = 0.0;
  ClusterAcceleration[2] = 0.0;

#ifdef COM_EVOLUTION
  double total_mass = 0.;
  if (NumberOfSingleParticle != 0) {
    for (int i = 0; i < NumberOfSingleParticle; i++) {
      for (int dim = 0; dim < Dim; dim++) {
        ClusterAcceleration[dim] += Mass[i] * BackgroundAcceleration[dim][i];
      }
      total_mass += Mass[i];
    }
    std::cout << std::setw(width)
              << "NBODY : ClusterAcceleration = " << std::setw(width)
              << ClusterAcceleration[0] << std::setw(width)
              << ClusterAcceleration[1] << std::setw(width)
              << ClusterAcceleration[2] << std::endl;

    std::cout << std::setw(width)
              << "NBODY : ClusterPosition = " << std::setw(width)
              << ClusterPosition[0] << std::setw(width) << ClusterPosition[1]
              << std::setw(width) << ClusterPosition[2] << std::endl;
    std::cout << std::setw(width)
              << "NBODY : ClusterVelocity = " << std::setw(width)
              << ClusterVelocity[0] << std::setw(width) << ClusterVelocity[1]
              << std::setw(width) << ClusterVelocity[2] << std::endl;
  }

  if (newNumberOfSingleParticle != 0) {
    // we need to make adjustment to COM
    GetNewCenterOfMass(PID, newMass, newPosition, newVelocity,
                       newNumberOfSingleParticle, ClusterPosition,
                       ClusterVelocity);

    for (int i = 0; i < newNumberOfSingleParticle; i++) {
      for (int dim = 0; dim < Dim; dim++) {
        ClusterAcceleration[dim] +=
            newMass[i] * newBackgroundAcceleration[dim][i];
      }
      total_mass += newMass[i];
    }
  }
  if (total_mass != 0.) {
    for (int dim = 0; dim < Dim; dim++) {
      ClusterAcceleration[dim] /= total_mass;
    }
  }
#endif

  // Update Existing Particles
  // need to update if the ids match between Enzo and Nbody
  if (NumberOfSingleParticle != 0) {
    std::cerr << "In ReceiveFromEnzo, NumberOfSingleParticle = "
              << NumberOfSingleParticle << std::endl;

#ifdef FEWBODY
    std::unordered_set<int> CMPtclsSet;
    CMPtclsSet.reserve(CMPtclWorker.size());
#endif
    Particle *ptcl;

    // loop for PID, going backwards to update the NextParticle
    for (int i = 0; i < NumberOfSingleParticle; i++) {

      ptcl = &particles[PIDtoIndexMap[PID[i]]];

#ifdef FEWBODY
      if (!ptcl->isActive && ptcl->CMPtclIndex != -1) {
        CMPtclsSet.insert(ptcl->CMPtclIndex);
      }
#endif

#ifdef COM_EVOLUTION
      for (int dim = 0; dim < Dim; dim++) {
        BackgroundAcceleration[dim][i] -= ClusterAcceleration[dim];
      }
#endif
      EnzoPIDs[i] = PID[i];
      ptcl->update(Mass, BackgroundAcceleration, i);
    } // endfor i
#ifdef FEWBODY
    for (int i :
         CMPtclsSet) { // Calculate background acceleration on the CM particles
      Particle *ptcl = &particles[i];
      double BackgroundAccelerationCM[Dim] = {0, 0, 0};

      for (int j = 0; j < ptcl->NumberOfMember; j++) {
        Particle *members = &particles[ptcl->Members[j]];

        for (int dim = 0; dim < Dim; dim++)
          BackgroundAccelerationCM[dim] +=
              members->Mass * members->BackgroundAcceleration[dim];
      }
      for (int dim = 0; dim < Dim; dim++)
        ptcl->BackgroundAcceleration[dim] =
            BackgroundAccelerationCM[dim] / ptcl->Mass;
    }
#endif
  } // endif nnb

#define no_star_formation_location_test
#ifdef star_formation_location_test
  std::vector<Particle *> new_particle;
#endif
  // Particle *newPtcl = new Particle[newNumberOfSingleParticle]; // we
  // shouldn't delete it, maybe make it to vector

  if (newNumberOfSingleParticle > 0) {
    std::cerr << "In ReceiveFromEnzo, newNumberOfSingleParticle = "
              << newNumberOfSingleParticle << std::endl;
    /*
    for (int i = 0; i < newNumberOfSingleParticle; i++) {
            fprintf(stderr, "NBODY: Mass Of NewNbodyParticles=%.3e\n",
    newMass[i]); fprintf(stderr, "NBODY: Vel  Of news=(%.3e, %.3e, %.3e)\n",
                            newVelocity[0][i], newVelocity[1][i],
    newVelocity[2][i]); fprintf(stderr, "NBODY: Pos  Of news=(%.3e, %.3e,
    %.3e)\n", newPosition[0][i], newPosition[1][i], newPosition[2][i]);
            fprintf(stderr, "NBODY: Acc  Of news=(%.3e, %.3e, %.3e)\n",
                            newBackgroundAcceleration[0][i],
    newBackgroundAcceleration[1][i], newBackgroundAcceleration[2][i]);
    }
    */

    // Update New Particles
    int index = -1;
    for (int i = 0; i < newNumberOfSingleParticle; i++) {
      for (int dim = 0; dim < Dim; dim++) {
#ifdef COM_EVOLUTION
        newBackgroundAcceleration[dim][i] -= ClusterAcceleration[dim];
        newVelocity[dim][i] -= ClusterVelocity[dim];
#endif
        newPosition[dim][i] -= ClusterPosition[dim];
      }

      if (NumberOfAvailableIndices != 0) {
        index = AvailableIndices[NumberOfAvailableIndices - 1];
        AvailableIndices[NumberOfAvailableIndices - 1] = -1;
        NumberOfAvailableIndices--;
      } else {
        index = global_variable->LastParticleIndex + 1;
        global_variable->LastParticleIndex++;
      }

      particles[index].set(newPID, newMass, newCreationTime, newDynamicalTime,
                           newMetallicity, newPosition, newVelocity,
                           newBackgroundAcceleration, i);
      particles[index].ParticleIndex = index;
      PIDtoIndexMap.insert({newPID[i], index});
      EnzoPIDs[NumberOfSingleParticle + i] = newPID[i];
      fprintf(stderr, "PID: %d is newly added from Enzo!\n", newPID[i]);
      fprintf(stderr, "PID: %d. Mass: %e Msun\n", newPID[i],
              particles[index].Mass * mass_unit);
      fprintf(stderr, "PID: %d. x: %e pc, y: %e pc, z: %e pc\n", newPID[i],
              particles[index].Position[0] * position_unit,
              particles[index].Position[1] * position_unit,
              particles[index].Position[2] * position_unit);
      fprintf(stderr, "PID: %d. vx: %e km/s, vy: %e km/s, vz: %e km/s\n",
              newPID[i],
              particles[index].Velocity[0] * velocity_unit / yr * pc / 1e5,
              particles[index].Velocity[1] * velocity_unit / yr * pc / 1e5,
              particles[index].Velocity[2] * velocity_unit / yr * pc / 1e5);

      // added for debugging in back acceleration
      fprintf(nbpout, "PID: %d is newly added from Enzo!\n", newPID[i]);
      fprintf(nbpout, "PID: %d. Mass: %e Msun\n", newPID[i],
              particles[index].Mass * mass_unit);
      fprintf(nbpout, "PID: %d. x: %e pc, y: %e pc, z: %e pc\n", newPID[i],
              particles[index].Position[0] * position_unit,
              particles[index].Position[1] * position_unit,
              particles[index].Position[2] * position_unit);
      fprintf(nbpout, "PID: %d. vx: %e km/s, vy: %e km/s, vz: %e km/s\n",
              newPID[i],
              particles[index].Velocity[0] * velocity_unit / yr * pc / 1e5,
              particles[index].Velocity[1] * velocity_unit / yr * pc / 1e5,
              particles[index].Velocity[2] * velocity_unit / yr * pc / 1e5);
      fprintf(nbpout, "Reg ax: %e, ay: %e, az: %e\n",
              particles[index].a_reg[0][0], particles[index].a_reg[1][0],
              particles[index].a_reg[2][0]);
      fprintf(nbpout, "Irr ax: %e, ay: %e, az: %e\n",
              particles[index].a_irr[0][0], particles[index].a_irr[1][0],
              particles[index].a_irr[2][0]);
      fprintf(nbpout, "Bgd ax: %e, ay: %e, az: %e\n",
              particles[index].BackgroundAcceleration[0],
              particles[index].BackgroundAcceleration[1],
              particles[index].BackgroundAcceleration[2]);

#ifdef star_formation_location_test
      new_particle.push_back(particles[NumberOfSingleParticle + i]);
#endif

#ifdef SEVN
      if (newCreationTime[i] > 0.0 && newMass[i] * EnzoMass * mass_unit > 2.2)
        initializeStellarEvolution(index);
#endif
    }

#ifdef star_formation_location_test
    writeParticle(new_particle, EnzoCurrentTime, outNum++);
#endif

    /*
    std::cout << "enzo Pos  :" << newPosition[0][0] << ", " << newPosition[0][1]
    << std::endl; std::cout << "nbody Pos :" << newPosition[0][0]*EnzoLength <<
    ", " << newPosition[0][1]*EnzoLength << std::endl; std::cout << "enzo Vel :"
    << newVelocity[1][0] << ", " << newVelocity[1][1] << std::endl; std::cout <<
    "nbody Vel :" << newVelocity[1][0]*EnzoVelocity << ", " <<
    newVelocity[1][1]*EnzoVelocity << std::endl; std::cout << "km/s Vel  :" <<
    newVelocity[1][0]*velocity_unit/1e5/yr*pc << ", " <<
    newVelocity[1][1]*velocity_unit/1e5/yr*pc << std::endl; std::cout << "enzo
    Mass:" << newMass[0] << std::endl; std::cout << "nbody Mass:" <<
    newMass[0]*EnzoMass << std::endl;
    */
    std::cout << "ABYSS    : " << newNumberOfSingleParticle
              << " new particles loaded!" << std::endl;

    // This includes modification of regular force and irregular force
  } // endif newNumberOfSingleParticle != 0

  /*
  // update particle order
  if (particle.size() > 1) {
          i=0;
          for (Particle* ptcl:particle) {
                  ptcl->ParticleOrder = i;
                  i++;
          }
  }
  */

  // (SEVN Query) After processing, wind & SN feedback,
  // 1. dm should be set to 0 - in update function
  // 2. If the particle is kicked, kicked velocity should be accounted - in
  // update function
  // 3. For PISN case, SEVN memory should be free - not yet

  if (NumberOfSingleParticle != 0) {
    delete[] PID;
    delete[] Mass;
    for (int dim = 0; dim < Dim; dim++) {
      delete[] BackgroundAcceleration[dim];
    }
  }

  if (newNumberOfSingleParticle != 0) {
    delete[] newPID;
    delete[] newMass;
    delete[] newCreationTime;
    delete[] newDynamicalTime;
    delete[] newMetallicity;
    for (int dim = 0; dim < Dim; dim++) {
      delete[] newBackgroundAcceleration[dim];
      delete[] newPosition[dim];
      delete[] newVelocity[dim];
    }
  }

  NumberOfSingleParticle += newNumberOfSingleParticle;
  NumberOfParticle += newNumberOfSingleParticle;

  //  (Query) Do I need this?
  // fprintf(nbpout, "ABYSS    : Acceleration for particles on GPU.\n");
  /*
  if (NumberOfSingleParticle > 1) {
          CalculateAllAccelerationOnGPU(particle);
  }*/

  // UpdateNextRegTime(particle);
  //  fprintf(nbpout, "ABYSS    : Acceleration and neighbors are updated.\n");

  fprintf(
      nbpout,
      "ABYSS    : In ReceiveFromEnzo (after new particle might be added): \n");
  fprintf(nbpout,
          "ABYSS    : original NumberOfSingleParticle      = %d (+%d)\n",
          NumberOfSingleParticle - newNumberOfSingleParticle,
          newNumberOfSingleParticle);
  fprintf(nbpout, "ABYSS    : newly updated NumberOfSingleParticle = %d\n",
          NumberOfSingleParticle, newNumberOfSingleParticle);
  // fprintf(nbpout, "ABYSS    : Particle size     = %d\n", particle.size());
  //  fprintf(nbpout, "ABYSS    : NextRegTimeStep   = %.3e\n",
  //  NextRegTimeBlock*global_variable->time_step); fprintf(nbpout, "ABYSS    :
  //  NextRegTimeBlock  = %d\n", NextRegTimeBlock);
  // fprintf(nbpout, "ABYSS    : RegularList size  = %d\n",
  // RegularList.size());
  fprintf(nbpout, "ABYSS    : FixNumNeighbor    = %d\n", FixNumNeighbor);

  fprintf(
      stderr,
      "ABYSS    : In ReceiveFromEnzo (after new particle might be added): \n");
  fprintf(stderr,
          "ABYSS    : original NumberOfSingleParticle      = %d (+%d)\n",
          NumberOfSingleParticle - newNumberOfSingleParticle,
          newNumberOfSingleParticle);
  fprintf(stderr, "ABYSS    : newly updated NumberOfSingleParticle = %d\n",
          NumberOfSingleParticle, newNumberOfSingleParticle);
  // fprintf(stderr, "ABYSS    : Particle size     = %d\n", particle.size());
  // fprintf(stderr, "ABYSS    : RegularList size = %d\n", RegularList.size());
  fflush(stderr);
  fflush(nbpout);
  // fflush(gpuout);
  // fflush(binout);

  return true;
}

// (Query) For binary mergers and SEVN, mass should be sent to Enzo as it was
// changed in Abyss by EW 2025.3.13 (Query) For SEVN, not only mass but also dm
// should be sent to Enzo, and it should be distributed into the grid by EW
// 2025.3.13
int SendParticleToEnzo(Worker *workers) {

  std::cout << "ABYSS: Entering SendToEnzo..." << std::endl;
  if (NumberOfSingleParticle == 0 && newNumberOfSingleParticle == 0) {
    std::cout << "ABYSS: Skipping SendToEnzo..." << std::endl;
    return 1; // SUCCESS -> 1 by EW 2025.3.11
  }
  MPI_Request request;
  MPI_Status status;

  double *Position[Dim], *Velocity[Dim], *newPosition[Dim], *newVelocity[Dim];
#ifdef SEVN
  double *InitialMass, *WindEjectedMass, *SNEjectedMass,
      *Temperature; // Msun & K unit
  double *newInitialMass, *newWindEjectedMass, *newSNEjectedMass,
      *newTemperature;    // Msun & K unit
  double *Mass, *newMass; // We need to update mass here
#endif
  int index;
  Particle *ptcl;

  for (int dim = 0; dim < Dim; dim++) {
    if (NumberOfSingleParticle - newNumberOfSingleParticle != 0) {
      Position[dim] =
          new double[NumberOfSingleParticle - newNumberOfSingleParticle];
      Velocity[dim] =
          new double[NumberOfSingleParticle - newNumberOfSingleParticle];
    }

    if (newNumberOfSingleParticle != 0) {
      newPosition[dim] = new double[newNumberOfSingleParticle];
      newVelocity[dim] = new double[newNumberOfSingleParticle];
    }
  }
#ifdef SEVN
  if (NumberOfSingleParticle - newNumberOfSingleParticle != 0) {
    InitialMass =
        new double[NumberOfSingleParticle - newNumberOfSingleParticle];
    WindEjectedMass =
        new double[NumberOfSingleParticle - newNumberOfSingleParticle];
    SNEjectedMass =
        new double[NumberOfSingleParticle - newNumberOfSingleParticle];
    Temperature =
        new double[NumberOfSingleParticle - newNumberOfSingleParticle];

    Mass = new double[NumberOfSingleParticle - newNumberOfSingleParticle];
  }

  if (newNumberOfSingleParticle != 0) {
    newInitialMass = new double[newNumberOfSingleParticle];
    newWindEjectedMass = new double[newNumberOfSingleParticle];
    newSNEjectedMass = new double[newNumberOfSingleParticle];
    newTemperature = new double[newNumberOfSingleParticle];

    newMass = new double[newNumberOfSingleParticle];
  }
#endif

  /*
  std::cout << "ABYSS: NumberOfSingleParticle=" << NumberOfSingleParticle << ",
  newNumberOfSingleParticle=" << newNumberOfSingleParticle << std::endl;
  std::cout << "ABYSS: size=" << particle.size() << std::endl;
  std::cout << "ABYSS: FirstParticleInEnzo PID=" << \
          FirstParticleInEnzo->PID << " in SendToEzno" << std::endl;
  */

  int NumberOfEscapeParticle = 0;
  double r2;
  double TimeStep = global_variable->EnzoTimeStep / EnzoTime;

#ifdef COM_EVOLUTION
  // COM evolution
  for (int dim = 0; dim < Dim; dim++) {
    ClusterPosition[dim] += ClusterVelocity[dim] * TimeStep;
    ClusterPosition[dim] += ClusterAcceleration[dim] * TimeStep * TimeStep / 2;
    ClusterVelocity[dim] += ClusterAcceleration[dim] * TimeStep;
  }
#else

#define UPDATE_InitialNeighborRadius2
#ifdef UPDATE_InitialNeighborRadius2
  std::vector<std::vector<double>> star(NumberOfSingleParticle,
                                        std::vector<double>(3, 0.0));
  std::vector<double> rarray(NumberOfSingleParticle);
  double NbodyCOM_nbodyunit[Dim] = {0, 0, 0};
#endif
  double NbodyCOM[Dim] = {0, 0, 0};
  double mass = 0;
  for (int i = 0; i < NumberOfSingleParticle; i++) {
    index = PIDtoIndexMap[EnzoPIDs[i]];
    ptcl = &particles[index];

    for (int dim = 0; dim < Dim; dim++) {
      NbodyCOM[dim] += ptcl->Mass * ptcl->Position[dim];
#ifdef UPDATE_InitialNeighborRadius2
      star[i][dim] = ptcl->Position[dim];
#endif
    }
    mass += ptcl->Mass;
  }

  for (int dim = 0; dim < Dim; dim++) {
    NbodyCOM[dim] /= mass;
#ifdef UPDATE_InitialNeighborRadius2
    NbodyCOM_nbodyunit[dim] = NbodyCOM[dim];
#endif
    NbodyCOM[dim] /= EnzoLength;
    NbodyCOM[dim] += ClusterPosition[dim];
  }
#endif

#ifdef FEWBODY
  std::unordered_set<int> CMPtclsSet;
  CMPtclsSet.reserve(CMPtclWorker.size());
  Particle *ptclCM;
#endif

  if (NumberOfSingleParticle - newNumberOfSingleParticle > 0) {
    // fprintf(stderr, "Sending PID order= ");

    for (int i = 0; i < NumberOfSingleParticle - newNumberOfSingleParticle;
         i++) {
      // fprintf(stderr, "%d, ",ptcl->PID);
      index = PIDtoIndexMap[EnzoPIDs[i]];
      ptcl = &particles[index];

#ifdef UPDATE_InitialNeighborRadius2
      for (int dim = 0; dim < Dim; dim++)
        star[i][dim] -= NbodyCOM_nbodyunit[dim];
      rarray[i] = sqrt(star[i][0] * star[i][0] + star[i][1] * star[i][1] +
                       star[i][2] * star[i][2]);
#endif

#ifdef FEWBODY
      if (!ptcl->isActive) {
        for (int dim = 0; dim < Dim; dim++) {
          Position[dim][i] =
              ptcl->Position[dim] / EnzoLength + ClusterPosition[dim];
          Velocity[dim][i] =
              ptcl->Velocity[dim] / EnzoVelocity + ClusterVelocity[dim];
        }
#ifdef SEVN
        InitialMass[i] = ptcl->InitialMass; // This is already in Msun unit!!!
        WindEjectedMass[i] = ptcl->dm * mass_unit;
        SNEjectedMass[i] = ptcl->SNEjectedMass * mass_unit;
        Temperature[i] = ptcl->T_eff;

        if (WindEjectedMass[i] > 0.0 || SNEjectedMass[i] > 0.0) {
          fprintf(stderr, "Feedback info send to Enzo...\n");
          fprintf(stderr,
                  "\tPID: %d. InitialMass: %e Msun, WindEjectedMass: %e Msun, "
                  "SNEjectedMass: %e Msun, T_eff: %e K\n",
                  ptcl->PID, InitialMass[i], WindEjectedMass[i],
                  SNEjectedMass[i], Temperature[i]);
        }

        Mass[i] = ptcl->Mass / EnzoMass;
#endif

        if (ptcl->CMPtclIndex != -1) {
          ptclCM = &particles[ptcl->CMPtclIndex];
          if (CMPtclsSet.find(ptcl->CMPtclIndex) == CMPtclsSet.end()) {

            CMPtclsSet.insert(ptcl->CMPtclIndex);
            ptclCM->NewNumberOfNeighbor = 0;
          }
          ptclCM->NewNeighbors[ptclCM->NewNumberOfNeighbor++] = i;
        } else if (ptcl->Mass <
                   0.0) { // merger induced zero-mass particles, PISN case

          fprintf(stdout,
                  "In CommunicationToHydro... PID: %d should be removed!\n",
                  ptcl->PID);

          // /* // This is temporarilly commented out by EW 2025.4.13
          deleteParticle(EnzoPIDs[i], index); // delete this particle in Abyss
          NumberOfEscapeParticle++;
          // */
        } else {
          fprintf(stderr, "What's wrong? PID: %d\n", ptcl->PID);
          throw std::runtime_error("Why inactive particle in EnzoPIDs?");
        }
        continue;
      }
#endif

      r2 = 0;
      for (int dim = 0; dim < Dim; dim++) {
        Position[dim][i] = ptcl->Position[dim] / EnzoLength;
        Velocity[dim][i] = ptcl->Velocity[dim] / EnzoVelocity;

#ifdef COM_EVOLUTION
        if (IdentifyNbodyParticles && IdentifyOnTheFly)
          r2 += Position[dim][i] * Position[dim][i];
        // COM correction
        Position[dim][i] += ClusterPosition[dim];
        Velocity[dim][i] += ClusterVelocity[dim];
#else

        // COM correction
        Position[dim][i] += ClusterPosition[dim];
        if (IdentifyNbodyParticles && IdentifyOnTheFly)
          r2 += (Position[dim][i] - NbodyCOM[dim]) *
                (Position[dim][i] - NbodyCOM[dim]);
#endif

        if (IdentifyNbodyParticles && !IdentifyOnTheFly)
          r2 += (Position[dim][i] - EnzoClusterPosition[dim]) *
                (Position[dim][i] - EnzoClusterPosition[dim]);
      }

      if (IdentifyNbodyParticles && ClusterRadius2 > 0 &&
          r2 > ClusterRadius2) { // in Enzo Unit
        Position[0][i] -= 20;    // original code
        // Position[0][i] -= 10*EnzoClusterPosition[0]; // 20; //fix this?
        deleteParticle(EnzoPIDs[i], index);
        fprintf(stderr, "In SendToEnzo... PID: %d is escaping!\n", ptcl->PID);
#ifdef FEWBODY
        NumberOfParticle--;
#ifdef SEVN
        if (ptcl->StellarEvolution != nullptr) {
          delete ptcl->StellarEvolution;
          ptcl->StellarEvolution = nullptr;
          fprintf(stderr, "In SendToEnzo... PID: %d. SEVN memory is freed!\n",
                  ptcl->PID);
          auto it = SEVNList.begin();
          while (it != SEVNList.end()) {
            if (it->second == ptcl->ParticleIndex) {
              it = SEVNList.erase(it);
              fprintf(stderr, "In SendToEnzo... PID: %d. SEVNList is erased!\n",
                      ptcl->PID);
              break;
            } else
              it++;
          }
        }
#endif
#endif
        NumberOfEscapeParticle++;
      }
      // fprintf(stdout, "ABYSS: pid= %d, x=%e\n",ptcl->PID,Position[0][i]);

#ifdef SEVN
      InitialMass[i] = ptcl->InitialMass; // This is already in Msun unit!!!
      WindEjectedMass[i] = ptcl->dm * mass_unit;
      SNEjectedMass[i] = ptcl->SNEjectedMass * mass_unit;
      Temperature[i] = ptcl->T_eff;

      if (WindEjectedMass[i] > 0.0 || SNEjectedMass[i] > 0.0) {
        fprintf(stderr, "Feedback info send to Enzo...\n");
        fprintf(stderr,
                "\tPID: %d. InitialMass: %e Msun, WindEjectedMass: %e Msun, "
                "SNEjectedMass: %e Msun, T_eff: %e K\n",
                ptcl->PID, InitialMass[i], WindEjectedMass[i], SNEjectedMass[i],
                Temperature[i]);
      }

      Mass[i] = ptcl->Mass / EnzoMass;
#endif
    }
  }

  int offset = NumberOfSingleParticle;

  if (newNumberOfSingleParticle > 0) {

    offset = NumberOfSingleParticle - newNumberOfSingleParticle;
    for (int i = 0; i < newNumberOfSingleParticle; i++) {
      index = PIDtoIndexMap[EnzoPIDs[i + offset]];
      ptcl = &particles[index];

#ifdef UPDATE_InitialNeighborRadius2
      for (int dim = 0; dim < Dim; dim++)
        star[i + offset][dim] -= NbodyCOM_nbodyunit[dim];
      rarray[i + offset] = sqrt(star[i + offset][0] * star[i + offset][0] +
                                star[i + offset][1] * star[i + offset][1] +
                                star[i + offset][2] * star[i + offset][2]);
#endif

#ifdef FEWBODY
      if (!ptcl->isActive) {
        for (int dim = 0; dim < Dim; dim++) {
          newPosition[dim][i] =
              ptcl->Position[dim] / EnzoLength + ClusterPosition[dim];
          newVelocity[dim][i] =
              ptcl->Velocity[dim] / EnzoVelocity + ClusterVelocity[dim];
        }
        // (SEVN Query) This particle should be deleted in Enzo too!!!
#ifdef SEVN
        newInitialMass[i] =
            ptcl->InitialMass; // This is already in Msun unit!!!
        newWindEjectedMass[i] = ptcl->dm * mass_unit;
        newSNEjectedMass[i] = ptcl->SNEjectedMass * mass_unit;
        newTemperature[i] = ptcl->T_eff;

        if (newWindEjectedMass[i] > 0.0 || newSNEjectedMass[i] > 0.0) {
          fprintf(stderr, "Feedback info send to Enzo...\n");
          fprintf(stderr,
                  "\tPID: %d. InitialMass: %e Msun, WindEjectedMass: %e Msun, "
                  "SNEjectedMass: %e Msun, T_eff: %e K\n",
                  ptcl->PID, newInitialMass[i], newWindEjectedMass[i],
                  newSNEjectedMass[i], newTemperature[i]);
        }

        newMass[i] = ptcl->Mass / EnzoMass;
#endif
        if (ptcl->CMPtclIndex != -1) {
          ptclCM = &particles[ptcl->CMPtclIndex];
          if (CMPtclsSet.find(ptcl->CMPtclIndex) == CMPtclsSet.end()) {

            CMPtclsSet.insert(ptcl->CMPtclIndex);
            ptclCM->NewNumberOfNeighbor = 0;
          }
          ptclCM->NewNeighbors[ptclCM->NewNumberOfNeighbor++] = i + offset;
        } else if (ptcl->Mass <
                   0.0) { // merger induced zero-mass particles, PISN case

          fprintf(stdout,
                  "In CommunicationToHydro... PID: %d should be removed!\n",
                  ptcl->PID);
          // /* // This is temporarilly commented out by EW 2025.4.13
          deleteParticle(EnzoPIDs[i + offset],
                         index); // delete this particle in Abyss
          NumberOfEscapeParticle++;
          // */
        } else {
          fprintf(stderr, "What's wrong? PID: %d\n", ptcl->PID);
          throw std::runtime_error("Why inactive particle in EnzoPIDs?");
        }
        continue;
      }
#endif

      r2 = 0;
      for (int dim = 0; dim < Dim; dim++) {
        newPosition[dim][i] = ptcl->Position[dim] / EnzoLength;
        newVelocity[dim][i] = ptcl->Velocity[dim] / EnzoVelocity;

#ifdef COM_EVOLUTION
        if (IdentifyNbodyParticles && IdentifyOnTheFly)
          r2 += newPosition[dim][i] * newPosition[dim][i];
        // COM correction
        newPosition[dim][i] += ClusterPosition[dim];
        newVelocity[dim][i] += ClusterVelocity[dim];
#else
        // COM correction
        newPosition[dim][i] += ClusterPosition[dim];
        if (IdentifyNbodyParticles && IdentifyOnTheFly)
          r2 += (newPosition[dim][i] - NbodyCOM[dim]) *
                (newPosition[dim][i] - NbodyCOM[dim]);
#endif

        if (IdentifyNbodyParticles && !IdentifyOnTheFly)
          r2 += (newPosition[dim][i] - EnzoClusterPosition[dim]) *
                (newPosition[dim][i] - EnzoClusterPosition[dim]);
      }
      if (IdentifyNbodyParticles && ClusterRadius2 > 0 && r2 > ClusterRadius2) {
        newPosition[0][i] -= 20; // original code
        // newPosition[0][i] -= 10*EnzoClusterPosition[0];
        deleteParticle(EnzoPIDs[i + offset], index);
#ifdef FEWBODY
        NumberOfParticle--;
#ifdef SEVN
        if (ptcl->StellarEvolution != nullptr) {
          delete ptcl->StellarEvolution;
          ptcl->StellarEvolution = nullptr;
          fprintf(stderr, "In SendToEnzo... PID: %d. SEVN memory is freed!\n",
                  ptcl->PID);
          auto it = SEVNList.begin();
          while (it != SEVNList.end()) {
            if (it->second == ptcl->ParticleIndex) {
              it = SEVNList.erase(it);
              fprintf(stderr, "In SendToEnzo... PID: %d. SEVNList is erased!\n",
                      ptcl->PID);
              break;
            } else
              it++;
          }
        }
#endif
#endif
        NumberOfEscapeParticle++;
        // (Query) binary termination?
      }

#ifdef SEVN
      newInitialMass[i] = ptcl->InitialMass; // This is already in Msun unit!!!
      newWindEjectedMass[i] = ptcl->dm * mass_unit;
      newSNEjectedMass[i] = ptcl->SNEjectedMass * mass_unit;
      newTemperature[i] = ptcl->T_eff;

      if (WindEjectedMass[i] > 0.0 || SNEjectedMass[i] > 0.0) {
        fprintf(stderr, "Feedback info send to Enzo...\n");
        fprintf(stderr,
                "\tPID: %d. InitialMass: %e Msun, WindEjectedMass: %e Msun, "
                "SNEjectedMass: %e Msun, T_eff: %e K\n",
                ptcl->PID, InitialMass[i], WindEjectedMass[i], SNEjectedMass[i],
                Temperature[i]);
      }

      newMass[i] = ptcl->Mass / EnzoMass;
#endif
    }
  }

#ifdef FEWBODY
  for (int i : CMPtclsSet) {
    ptcl = &particles[i];

    assert(ptcl->NewNumberOfNeighbor == ptcl->NumberOfMember);

    double memPosition[3];
    double memVelocity[3];
    bool memEscape = true;

    Particle *members;

    for (int j = 0; j < ptcl->NumberOfMember; j++) {
      members = &particles[ptcl->Members[j]];

      r2 = 0;
      for (int dim = 0; dim < Dim; dim++) {
        memPosition[dim] = members->Position[dim] / EnzoLength;
        memVelocity[dim] = members->Velocity[dim] / EnzoVelocity;

#ifdef COM_EVOLUTION
        if (IdentifyNbodyParticles && IdentifyOnTheFly)
          r2 += memPosition[dim] * memPosition[dim];
        // COM correction
        memPosition[dim] += ClusterPosition[dim];
        memVelocity[dim] += ClusterVelocity[dim];
#else

        // COM correction
        memPosition[dim] += ClusterPosition[dim];
        if (IdentifyNbodyParticles && IdentifyOnTheFly)
          r2 += (memPosition[dim] - NbodyCOM[dim]) *
                (memPosition[dim] - NbodyCOM[dim]);
#endif

        if (IdentifyNbodyParticles && !IdentifyOnTheFly)
          r2 += (memPosition[dim] - EnzoClusterPosition[dim]) *
                (memPosition[dim] - EnzoClusterPosition[dim]);
      }

      if (IdentifyNbodyParticles && ClusterRadius2 > 0 &&
          r2 <= ClusterRadius2) { // in Enzo Unit
        memEscape = false;
        break;
      }
      // fprintf(stdout, "ABYSS: pid= %d, x=%e\n",ptcl->PID,Position[0][i]);
    }
    if (memEscape) {
      fprintf(stderr, "Binary escape!\n");
      fprintf(stdout, "Binary escape!\n");

      ptcl->isActive = false;
      NumberOfParticle--; // CM particle should be removed from active particles

      fprintf(stderr,
              "Binary escape... CM PID: %d, NumberOfMember: %d, "
              "NewNumberOfNeighbors: %d\n",
              ptcl->PID, ptcl->NumberOfMember, ptcl->NewNumberOfNeighbor);
      for (int j = 0; j < ptcl->NumberOfMember; j++) {
        members = &particles[ptcl->Members[j]];
        if (ptcl->NewNeighbors[j] < offset) {
          Position[0][ptcl->NewNeighbors[j]] -= 20;
          fprintf(stderr, "Binary escape1... mem PID: %d, i: %d\n",
                  members->PID, ptcl->NewNeighbors[j]);
        } else {
          newPosition[0][ptcl->NewNeighbors[j] - offset] -= 20;
          fprintf(stderr, "Binary escape2... mem PID: %d, i: %d\n",
                  members->PID, ptcl->NewNeighbors[j]);
        }
        deleteParticle(members->PID, ptcl->Members[j]);
        NumberOfEscapeParticle++;
#ifdef SEVN
        if (members->StellarEvolution != nullptr) {
          delete members->StellarEvolution;
          members->StellarEvolution = nullptr;
          fprintf(stderr, "In SendToEnzo... PID: %d. SEVN memory is freed!\n",
                  members->PID);
          auto it = SEVNList.begin();
          while (it != SEVNList.end()) {
            if (it->second == members->ParticleIndex) {
              it = SEVNList.erase(it);
              fprintf(stderr, "In SendToEnzo... PID: %d. SEVNList is erased!\n",
                      members->PID);
              break;
            } else
              it++;
          }
        }
#endif
      }

      int rank_delete = CMPtclWorker[ptcl->ParticleIndex];
      fprintf(nbpout, "Rank of CM ptcl %d: %d\n", ptcl->PID, rank_delete);
      Queue queue;
      queue.task = DeleteGroup;
      queue.pid = ptcl->ParticleIndex;
      workers[rank_delete].addQueue(queue);
      workers[rank_delete].runQueue();
      workers[rank_delete].callback();

      if (ptcl->ParticleIndex == global_variable->LastParticleIndex) {
        global_variable->LastParticleIndex--;
        // global_variable->LastParticleIndex == LastParticleIndex;
      } else
        PrevCMPtclWorker.insert(
            {ptcl->ParticleIndex, CMPtclWorker[ptcl->ParticleIndex]});
      CMPtclWorker.erase(ptcl->ParticleIndex);
    }
  }
#endif

#ifdef UPDATE_InitialNeighborRadius2
  std::sort(rarray.begin(), rarray.end());
  int NNBMAX = static_cast<int>(std::sqrt(NumberOfSingleParticle));
  double RS0 = rarray[NNBMAX];
  double originalInitialNeighborRadius2 = InitialNeighborRadius2;
  fprintf(stderr, "Original InitialNeighborRadius2: %e\n",
          originalInitialNeighborRadius2);
  fprintf(stderr, "Newly calculated InitialNeighborRadius2: %e\n", RS0 * RS0);
  if (NumberOfSingleParticle + newNumberOfSingleParticle > 1000)
    InitialNeighborRadius2 =
        std::min(RS0 * RS0, originalInitialNeighborRadius2);
#endif

  // std::cerr << "ABYSS: Waiting for Enzo to send data..." << std::endl;
  fprintf(nbpout, "ABYSS: Waiting for Enzo to send data...\n");
  CommunicationInterBarrier();
  if (NumberOfSingleParticle - newNumberOfSingleParticle != 0) {
    for (int dim = 0; dim < Dim; dim++) {
      MPI_Send(Position[dim],
               NumberOfSingleParticle - newNumberOfSingleParticle, MPI_DOUBLE,
               0, 300, inter_comm);
      MPI_Send(Velocity[dim],
               NumberOfSingleParticle - newNumberOfSingleParticle, MPI_DOUBLE,
               0, 400, inter_comm);
    }
#ifdef SEVN
    MPI_Send(InitialMass, NumberOfSingleParticle - newNumberOfSingleParticle,
             MPI_DOUBLE, 0, 800, inter_comm);
    MPI_Send(WindEjectedMass,
             NumberOfSingleParticle - newNumberOfSingleParticle, MPI_DOUBLE, 0,
             900, inter_comm);
    MPI_Send(SNEjectedMass, NumberOfSingleParticle - newNumberOfSingleParticle,
             MPI_DOUBLE, 0, 1000, inter_comm);
    MPI_Send(Temperature, NumberOfSingleParticle - newNumberOfSingleParticle,
             MPI_DOUBLE, 0, 1100, inter_comm);

    MPI_Send(Mass, NumberOfSingleParticle - newNumberOfSingleParticle,
             MPI_DOUBLE, 0, 1600, inter_comm);
#endif
  }
  // std::cerr << "ABYSS: Escape particles=" << EscapeParticleNum << std::endl;

  // fprintf(stderr,"NewNumberOfSingleParticles=%d\n",NumberOfNewNbodyParticles);
  if (newNumberOfSingleParticle != 0) {
    for (int dim = 0; dim < Dim; dim++) {
      MPI_Send(newPosition[dim], newNumberOfSingleParticle, MPI_DOUBLE, 0, 500,
               inter_comm);
      MPI_Send(newVelocity[dim], newNumberOfSingleParticle, MPI_DOUBLE, 0, 600,
               inter_comm);
    }
#ifdef SEVN
    MPI_Send(newInitialMass, newNumberOfSingleParticle, MPI_DOUBLE, 0, 1200,
             inter_comm);
    MPI_Send(newWindEjectedMass, newNumberOfSingleParticle, MPI_DOUBLE, 0, 1300,
             inter_comm);
    MPI_Send(newSNEjectedMass, newNumberOfSingleParticle, MPI_DOUBLE, 0, 1400,
             inter_comm);
    MPI_Send(newTemperature, newNumberOfSingleParticle, MPI_DOUBLE, 0, 1500,
             inter_comm);

    MPI_Send(newMass, newNumberOfSingleParticle, MPI_DOUBLE, 0, 1700,
             inter_comm);
#endif
  }

  if (NumberOfSingleParticle != 0 && IdentifyOnTheFly) {
    // fprintf(stderr, "NBODY: ClusterPosition =(%lf, %lf, %lf)\n",
    // ClusterPosition[0]-0.5, ClusterPosition[1]-0.5, ClusterPosition[2]-0.5);
    /*
    fprintf(stderr, "NBODY: ClusterPosition =(%lf, %lf, %lf)\n",
                    (ClusterPosition[0]-0.5)*EnzoLength*position_unit,
                    (ClusterPosition[1]-0.5)*EnzoLength*position_unit,
                    (ClusterPosition[2]-0.5)*EnzoLength*position_unit);*/
#ifdef COM_EVOLUTION
    MPI_Send(ClusterPosition, 3, MPI_DOUBLE, 0, 700, inter_comm);
#else
    MPI_Send(NbodyCOM, 3, MPI_DOUBLE, 0, 700, inter_comm);
#endif
  }

  CommunicationInterBarrier();
  fprintf(stdout, "ABYSS: Data sent!\n");

  for (int dim = 0; dim < Dim; dim++) {
    if (NumberOfSingleParticle - newNumberOfSingleParticle != 0) {
      delete[] Position[dim];
      delete[] Velocity[dim];
    }
    if (newNumberOfSingleParticle != 0) {
      delete[] newPosition[dim];
      delete[] newVelocity[dim];
    }
  }
#ifdef SEVN
  if (NumberOfSingleParticle - newNumberOfSingleParticle != 0) {
    delete[] InitialMass;
    delete[] WindEjectedMass;
    delete[] SNEjectedMass;
    delete[] Temperature;

    delete[] Mass;
  }
  if (newNumberOfSingleParticle != 0) {
    delete[] newInitialMass;
    delete[] newWindEjectedMass;
    delete[] newSNEjectedMass;
    delete[] newTemperature;

    delete[] newMass;
  }
#endif

  NumberOfSingleParticle -= NumberOfEscapeParticle;
#ifdef FEWBODY
  // NumberOfParticle -= NumberOfEscapeParticle; // (Query) EW: We don't have to
  // do this.
#else
  NumberOfParticle -= NumberOfEscapeParticle;
#endif
  std::cout << "ABYSS: Sending data finished." << std::endl;

  fprintf(stderr,
          "ABYSS    : In SendToEnzo (after particle might be escaped): \n");
  fprintf(
      stderr, "ABYSS    : original NumberOfSingleParticle      = %d (-%d)\n",
      NumberOfSingleParticle + NumberOfEscapeParticle, NumberOfEscapeParticle);
  fprintf(stderr, "ABYSS    : newly updated NumberOfSingleParticle = %d\n",
          NumberOfSingleParticle);
  fprintf(stderr, "ABYSS    : newly updated NumberOfParticle = %d\n",
          NumberOfParticle);

  fprintf(nbpout,
          "ABYSS    : In SendToEnzo (after particle might be escaped): \n");
  fprintf(
      nbpout, "ABYSS    : original NumberOfSingleParticle      = %d (-%d)\n",
      NumberOfSingleParticle + NumberOfEscapeParticle, NumberOfEscapeParticle);
  fprintf(nbpout, "ABYSS    : newly updated NumberOfSingleParticle = %d\n",
          NumberOfSingleParticle);

  fflush(stderr);
  fflush(nbpout);
  // fflush(gpuout);
  //  fflush(binout);
  return true;
}

/* Adjust COM according to new particles
X = (a1+b1+c1)/M
Y = (a2+b2)/N
Z = (a1+b1+c1+a2+b2)/(M+N) = (X*M+Y*N)/(M+N) -> New COM
Z-X = N*(Y-X)/(M+N) -> This should be applied to particles
*/
// this should be improved by using iterative loop

void GetNewCenterOfMass(int *PID, double *mass2, double *x2[Dim],
                        double *v2[Dim], int n2, double x_X[], double v_X[]) {

  double M = 0., N = 0.;
  double x_Y[Dim], v_Y[Dim], x_Z[Dim], v_Z[Dim];

  for (int dim = 0; dim < Dim; dim++) {
    x_Y[dim] = 0;
    x_Z[dim] = 0;
    v_Y[dim] = 0;
    v_Z[dim] = 0;
  }

  for (int i = 0; i < NumberOfSingleParticle; i++) {
    M += particles[PIDtoIndexMap[PID[i]]].Mass / EnzoMass;
  }

  for (int i = 0; i < n2; i++) {
    for (int dim = 0; dim < Dim; dim++) {
      x_Y[dim] += mass2[i] * x2[dim][i];
      v_Y[dim] += mass2[i] * v2[dim][i];
    }
    N += mass2[i];
  }
  for (int dim = 0; dim < Dim; dim++) {
    x_Y[dim] = x_Y[dim] / N;
    v_Y[dim] = v_Y[dim] / N;
    x_Z[dim] = (x_X[dim] * M + x_Y[dim] * N) / (M + N);
    v_Z[dim] = (v_X[dim] * M + v_Y[dim] * N) / (M + N);
  }

  // Adjustment to particles
  for (int i = 0; i < NumberOfSingleParticle; i++) {
    for (int dim = 0; dim < Dim; dim++) {
      particles[PIDtoIndexMap[PID[i]]].Position[dim] -
          (x_Z[dim] - x_X[dim]) * EnzoLength;
      particles[PIDtoIndexMap[PID[i]]].Velocity[dim] -
          (v_Z[dim] - v_X[dim]) * EnzoVelocity;
    }
  }

  for (int dim = 0; dim < Dim; dim++) {
    x_X[dim] = x_Z[dim];
    v_X[dim] = v_Z[dim];
  }
}

void deleteParticle(int &PID, int &index) {
  particles[index].isActive = false;
  PIDtoIndexMap.erase(PID);
  AvailableIndices[NumberOfAvailableIndices] = index;
  NumberOfAvailableIndices++;
}

#endif
#endif