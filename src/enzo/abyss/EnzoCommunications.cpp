// #include "FewBody/ar_interaction.hpp"
#include <cstdio>
#ifdef INDIVIDUALSTAR
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>
#define MAX_NUMBER_OF_OUTPUT_REDSHIFTS 500
#define FLOAT double
#include "../CosmologyParameters.h"
#undef FLOAT
//#include "../phys_constants.h"

#ifdef FEWBODY
#include <unordered_set>
#endif
#define ABYSS_ONLY
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
double ClusterRadius2;
extern int WorldProcessorNumber;
double ClusterAcceleration[Dim], ClusterPosition[Dim], ClusterVelocity[Dim], EnzoClusterPosition[Dim + 1];
double eta_tmp;
double a_i, a_f, dadt_i, dadt_f; // RedshiftEnzoStart, RedshiftEnzoEnd;
int FixNumNeighbor0, IdentifyNbodyParticles;
int BinaryRegularization, IdentifyOnTheFly, StoreTimeStep;

// double KSTime;
// double KSDistance;

void deleteParticle(int &PID, int &index);
// int writeParticle(std::vector<Particle*> &particle, double MinRegTime, int
// outputNum); void InitializeParticle(std::vector<Particle*> &particle); void
// InitializeNewParticle(std::vector<Particle*> &particle, int offset, int
// newSize);

void GetCenterOfMass(ParticleDataType *ptcl, const int &N,
                    double x_com[], double v_com[]);
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

/*
*
*
*
*
*
*
*
*/
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
        MPI_Recv(&CosmologyTableNumberOfBins, 1, MPI_INT, 0, 150, inter_comm, &statuses[1]);

    /* Prepare buffer array */
    double *CosmologyTableLoga;
    double *CosmologyTableLogt;
    double TimeStep, TimeUnits, LengthUnits, VelocityUnits, DensityUnits;
    int double_size = 14;
    int int_size = 7;
    if (ComovingCoordinates) {
        double_size += 12;
        double_size += 2 * CosmologyTableNumberOfBins;
        double_size += 4; // a_i, dadt_i, a_f, dadt_f
        int_size += 1;
        CosmologyTableLoga = new double[CosmologyTableNumberOfBins]; // (Query) EW: not deleted later?
        CosmologyTableLogt = new double[CosmologyTableNumberOfBins]; // (Query) EW: not deleted later?
    }
    double  *params_double   = new double[double_size];
    int     *params_int      = new int[int_size];

    MPI_Irecv(params_double,    double_size,    MPI_DOUBLE, 0, 200, inter_comm, &requests[0]);
    MPI_Irecv(params_int,       int_size,       MPI_INT,    0, 300, inter_comm, &requests[1]);
    MPI_Waitall(2, requests, MPI_STATUSES_IGNORE);
    fprintf(nbpout, "Data received!\n");

    TimeStep                    = params_double[0];
    TimeUnits                   = params_double[1];
    LengthUnits                 = params_double[2];
    DensityUnits                = params_double[3];
    VelocityUnits               = params_double[4];
    StarMassEjectionFraction    = params_double[5];
    EnzoCurrentTime             = params_double[6];
    EPS2                        = params_double[7];
    eta_tmp                     = params_double[8];
    InitialNeighborRadius2      = params_double[9];
    EnzoClusterPosition[0]      = params_double[10];
    EnzoClusterPosition[1]      = params_double[11];
    EnzoClusterPosition[2]      = params_double[12];
    EnzoClusterPosition[3]      = params_double[13];

    StarParticleFeedback        = params_int[0];
    IdentifyNbodyParticles      = params_int[1];
    IdentifyOnTheFly            = params_int[2];
    FixNumNeighbor              = params_int[3];
    MaxNumNeighbor              = params_int[4]; // (Query) This is compile parameter
    BinaryRegularization        = params_int[5];
    StoreTimeStep               = params_int[6];

    if (ComovingCoordinates) {
        HubbleConstantNow           = params_double[14];
        OmegaMatterNow              = params_double[15];
        OmegaDarkMatterNow          = params_double[16];
        OmegaLambdaNow              = params_double[17];
        OmegaRadiationNow           = params_double[18];
        ComovingBoxSize             = params_double[19];
        MaxExpansionRate            = params_double[20];
        InitialTimeInCodeUnits      = params_double[21];
        InitialRedshift             = params_double[22];
        FinalRedshift               = params_double[23];
        CosmologyTableLogaInitial   = params_double[24];
        CosmologyTableLogaFinal     = params_double[25];
        for (int i = 0; i < CosmologyTableNumberOfBins; i++) {
            CosmologyTableLoga[i]   = params_double[26 + i];
        }
        for (int i = 0; i < CosmologyTableNumberOfBins; i++) {
            CosmologyTableLogt[i]   = params_double[26 + CosmologyTableNumberOfBins + i];
        }
        a_i     = params_double[26 + 2 * CosmologyTableNumberOfBins];
        dadt_i  = params_double[26 + 2 * CosmologyTableNumberOfBins + 1];
        a_f     = params_double[26 + 2 * CosmologyTableNumberOfBins + 2];
        dadt_f  = params_double[26 + 2 * CosmologyTableNumberOfBins + 3];
        
        CosmologyTableLogtIndex = params_int[7];
    }
    ClusterRadius2 = EnzoClusterPosition[3]; // it's already squared

    delete[] params_double;
    delete[] params_int;
    // (Query by EW) What about CosmologyTableLoga and CosmologyTableLogt?

#ifdef COMOVE
    if (ComovingCoordinates){
        //in Enzo, a_i is defined as a = (1+InitialRedshift)/(1+z).
        //So we have to divide a_i and dadt_i by (1+InitialRedshift) to convert it into natural definition.
        a_i /= (1 + InitialRedshift); 
        a_f /= (1 + InitialRedshift);
        dadt_i /= (1 + InitialRedshift);
        dadt_f /= (1 + InitialRedshift);
        // convert units into comoving units
        LengthUnits /= a_i;
        DensityUnits *= pow(a_i,3.);
        // VelocityUnits /= a_i; we should "not" change VelocityUnits as Enzo already use peculiar velocity
        // TimeUnits?
    }
#endif

    // Enzo to Nbody unit convertors
    // EnzoMass         = MassUnits/Msun/mass_unit;
    EnzoMass            = DensityUnits * pow(LengthUnits, 3.) / Msun / mass_unit;
    EnzoLength          = LengthUnits / pc / position_unit;
    EnzoVelocity        = VelocityUnits / pc * yr / velocity_unit;
    EnzoTime            = TimeUnits / yr / time_unit;
    // EnzoAcceleration = LengthUnits/TimeUnits/TimeUnits/pc*yr*yr/position_unit*time_unit*time_unit;
    EnzoAcceleration    = EnzoVelocity / EnzoTime;

    // Unit conversion
    global_variable->EnzoTimeStep   = TimeStep * EnzoTime;
    EnzoCurrentTime                 *= EnzoTime;

    if (EPS2 < 0)
        EPS2 = -1;
    else {
        EPS2 *= EnzoLength;
        EPS2 *= EPS2;
    }

#ifdef COMOVE
    if (ComovingCoordinates){
        dadt_i /= EnzoTime;
        dadt_f /= EnzoTime;
        global_variable->a_i = a_i;
        global_variable->a_f = a_f;
        global_variable->dadt_i = dadt_i;
        global_variable->dadt_f = dadt_f;
        fprintf(stderr, "a_i = %e, dadt_i = %e, a_f = %e, dadt_f = %e\n", a_i, dadt_i, a_f, dadt_f);
        fprintf(stderr, "InitialRedshift = %e, FinalRedshift = %e\n", InitialRedshift, FinalRedshift);
        EPS2 *= (a_i/a_f)*(a_i/a_f); // Fix EPS2 for physical coordinates 
    }
#endif

    // InitialNeighborRadius2 is now in comoving unit with COMOVE. we may change this into physical unit
    InitialNeighborRadius2  *= EnzoLength;
    InitialNeighborRadius2  *= InitialNeighborRadius2;
    FixNumNeighbor0         = FixNumNeighbor;

    // need to fix units
    // fprintf(nbpout, "Enzo Time                = %lf\n", TimeStep);
    fprintf(nbpout, "LengthUnit                = %e\n", LengthUnits);
    fprintf(nbpout, "DensityUnit               = %e\n", DensityUnits);
    fprintf(nbpout, "TimeUnit                  = %e\n", TimeUnits);
    fprintf(nbpout, "VelocityUnit              = %e\n", VelocityUnits);



    fprintf(nbpout, "Nbody Time               = %e Myr\n", EnzoCurrentTime*1e4);
    fprintf(nbpout, "Nbody TimeStep           = %e\n", global_variable->EnzoTimeStep);
    fprintf(nbpout, "EPS2                     = %e pc**2\n", EPS2 * position_unit * position_unit);
    fprintf(nbpout, "InitialNeighborRadius2   = %.2e pc**2\n", InitialNeighborRadius2 * position_unit * position_unit);
    fprintf(nbpout, "eta                      = %e\n", eta);
    fprintf(nbpout, "EnzoClusterPosition      = (%e, %e, %e)\n", EnzoClusterPosition[0], EnzoClusterPosition[1], EnzoClusterPosition[2]);
    fprintf(nbpout, "ClusterRadius2           = %.2e pc**2\n", ClusterRadius2 * position_unit * position_unit);
    fprintf(nbpout, "StarMassEjectionFraction = %e\n", StarMassEjectionFraction);
    fprintf(nbpout, "StarParticleFeedback     = %d\n", StarParticleFeedback);
    fprintf(nbpout, "FixNumNeighbor           = %d\n", FixNumNeighbor);
    fprintf(nbpout, "BinaryRegularization     = %d\n", BinaryRegularization); // (Query) EW: What is this?
    // fprintf(nbpout, "KSTime                   = %lf\n", KSTime);
    // fprintf(nbpout, "KSDistance               = %lf\n", KSDistance);
    fprintf(nbpout, "IdentifyNbodyParticles   = %d\n", IdentifyNbodyParticles);
    fprintf(nbpout, "IdentifyOnTheFly         = %d\n\n", IdentifyOnTheFly);
    fflush(nbpout);

    /*------------------===-------------------------*/
    /********   Receive Particles to ABYSS  *********/
    /*----------------===---------------------------*/

    bool debug1 = true;

    int LocalNumberOfParticles=0;

    if (debug1) {
        fprintf(nbpout, "NumberOfProcessors=%d\n",NumberOfProcessors);
        fflush(nbpout);
    }

    /* Step 1: Gather sizes */
    int size = NumberOfProcessors+1; // the number of inter_comm processors
    int *recv_counts = new int[size];
    MPI_Gather(&LocalNumberOfParticles, 1, MPI_INT, recv_counts, 1, MPI_INT,
                NumberOfProcessors, inter_comm);


    /* Step 2: Compute displacements */
    NumberOfSingleParticle = 0;
    displs = new int[size]; // (Query by EW) Why is this parameter a global variable?
    displs[0] = 0;
    for (int i = 1; i < size; i++) {
        displs[i] = displs[i - 1] + recv_counts[i - 1];
    }
    NumberOfSingleParticle = displs[size - 1] + recv_counts[size - 1];


    if (debug1) {
        fprintf(nbpout, "displs= ");
        for (int i=0; i<size; i++)
        fprintf(nbpout, "%d, ", displs[i]);
        fprintf(nbpout, "\n NumberOfSingleParticle=%d\n",NumberOfSingleParticle);
        fflush(nbpout);
    }


    /* Step 3: Allocate receive buffer on root */
    ParticleDataType *recvbuf = new ParticleDataType[NumberOfSingleParticle];

    if (debug1) {
        fprintf(nbpout, "Receive Buffer Ready!\n");
        fflush(nbpout);
    }


    /* Step 4: Gatherv */
    MPI_Gatherv(NULL, LocalNumberOfParticles, MPI_ENZO_PTCL, recvbuf, recv_counts, displs,
                MPI_ENZO_PTCL, NumberOfProcessors, inter_comm);


    if (debug1) {
        //fprintf(stderr, "ABYSS: ID = ");
        fprintf(nbpout, "ABYSS: ID = ");
        for (int i=0; i<NumberOfSingleParticle; i++) {
            //fprintf(stderr, "(%d, ", recvbuf[i].ID);
            //fprintf(stderr, "%.4e, %.4e), ", recvbuf[i].Position[0], recvbuf[i].BackgroundAcceleration[0]);
            fprintf(nbpout, "(%d,  ", recvbuf[i].ID);
            //fprintf(nbpout, "%.4e, %.4e), ", recvbuf[i].Position[0], recvbuf[i].BackgroundAcceleration[0]);
        }
        //fprintf(stderr, "\n");
        fprintf(nbpout, ")\n");
        fflush(nbpout);
    }



    /* Unpack the Packet */
#pragma unroll Dim
    for (int dim = 0; dim < Dim; dim++) {
        ClusterAcceleration[dim] = 0;
        ClusterPosition[dim] = EnzoClusterPosition[dim]; // initialize into EnzoClusterPosition for the case when NumberOfSingleParticle == 0
        ClusterVelocity[dim] = 0;
    }
    if (NumberOfSingleParticle != 0) {
        // set COM for background acc
        GetCenterOfMass(recvbuf, NumberOfSingleParticle, ClusterPosition, ClusterVelocity);
                    

#ifdef COM_EVOLUTION
        ClusterAcceleration[0] = 0.;
        ClusterAcceleration[1] = 0.;
        ClusterAcceleration[2] = 0.;
        double total_mass = 0.;
        for (int i = 0; i < NumberOfSingleParticle; i++) {
            for (int dim = 0; dim < Dim; dim++) {
                ClusterAcceleration[dim] += recvbuf[i].Mass * recvbuf[i].BackgroundAcceleration[dim];
            }
            total_mass += recvbuf[i].Mass;
        }

        for (int dim = 0; dim < Dim; dim++)
            ClusterAcceleration[dim] /= total_mass;
#endif

        int EnzoProcessorNumber = 0;
        fprintf(nbpout, "EnzoProcessorNumber=");
        for (int i = 0; i < NumberOfSingleParticle; i++) {
#pragma unroll Dim
            for (int dim = 0; dim < Dim; dim++) {
#ifdef COM_EVOLUTION
                recvbuf[i].BackgroundAcceleration[dim] -= ClusterAcceleration[dim];
                recvbuf[i].Velocity[dim] -= ClusterVelocity[dim];
#endif
                recvbuf[i].Position[dim] -= ClusterPosition[dim];
            }
            //EnzoPIDs[i] = recvbuf[i].ID;
            while (i >= displs[EnzoProcessorNumber+1])
                EnzoProcessorNumber++;
            particles[i].setFirst(recvbuf[i], EnzoProcessorNumber);
            PIDtoIndexMap.insert({recvbuf[i].ID, i});
            particles[i].ParticleIndex = i;
            fprintf(nbpout, "%d, ", EnzoProcessorNumber);
        }
        fprintf(nbpout, "\n");
        fflush(nbpout);

        global_variable->EnzoCurrentTime = EnzoCurrentTime * 1e4; // This should be 0.0 // Unlike original value, this is Myr unit
        std::cerr << "EnzoCurrentTime:" << global_variable->EnzoCurrentTime << std::endl;

        // Update Global variable AbyssCenter. This will be used in Nbody to Enzo conversion in readwrite 
		AbyssCenter[0] = ClusterPosition[0];
		AbyssCenter[1] = ClusterPosition[1];
		AbyssCenter[2] = ClusterPosition[2];
    }

    NumberOfParticle = NumberOfSingleParticle;
    global_variable->LastParticleIndex = NumberOfSingleParticle - 1;

#ifdef SEVN
	initializeStellarEvolution();
#endif

    delete[] recv_counts;
    delete[] displs;
    delete[] recvbuf;

    fprintf(nbpout, "ABYSS: %d particles loaded!\n", NumberOfSingleParticle);
    fflush(nbpout);
    return true;
}


/*
*
*
*
*
*
*
*
*/



int ReceiveParticleFromEnzo() {

    MPI_Request request;
    MPI_Status status;
    double TimeStep;
    fprintf(nbpout, "ABYSS: Starting ReceiveParticle...\n");
    fprintf(stderr, "ABYSS: Starting ReceiveParticle...\n");
    //MPI_Barrier(inter_comm); 


    //(Query) is this for cosmology?
    /*-------------------------------------------*/
    /*******  Receive Parameters to ABYSS  *******/
    /*-------------------------------------------*/
    int double_size = 6;
    if (ComovingCoordinates)
        double_size += 4; // a_i, dadt_i, a_f, dadt_f
    double *params_double = new double[double_size];

    double TimeUnits, LengthUnits, VelocityUnits, DensityUnits;

    MPI_Recv(params_double, double_size, MPI_DOUBLE, 0, 100, inter_comm, MPI_STATUS_IGNORE);

    TimeStep        = params_double[0];
    TimeUnits       = params_double[1];
    EnzoCurrentTime = params_double[2];
    LengthUnits     = params_double[3];
    DensityUnits    = params_double[4];
    VelocityUnits   = params_double[5];

#ifdef COMOVE
    if (ComovingCoordinates) {
        a_i     = params_double[6] / (1 + InitialRedshift);
        a_f     = params_double[7] / (1 + InitialRedshift);
        dadt_i  = params_double[8] / (1 + InitialRedshift);
        dadt_f  = params_double[9] / (1 + InitialRedshift);
        fprintf(stderr, "a_i = %e, dadt_i = %e, a_f = %e, dadt_f = %e\n", a_i, dadt_i, a_f, dadt_f);
        fprintf(stderr, "InitialRedshift = %e, FinalRedshift = %e\n", InitialRedshift, FinalRedshift);
        // convert units into comoving units
        LengthUnits     /= a_i;
        DensityUnits    *= pow(a_i, 3.);
        fprintf(stderr, "ClusterPosition = %e, %e, %e\n", EnzoClusterPosition[0], EnzoClusterPosition[1], EnzoClusterPosition[2]);
        EPS2            *= (a_i / a_f) * (a_i / a_f); // Fix EPS2 for physical coordinates 
    }
#endif

    if (ComovingCoordinates) {
        // with COMOVE, we don't actually need to update units
        EnzoMass         = DensityUnits * pow(LengthUnits, 3.) / Msun / mass_unit;
        EnzoLength       = LengthUnits / pc / position_unit;
        EnzoVelocity     = VelocityUnits / pc * yr / velocity_unit;
        EnzoTime         = TimeUnits / yr / time_unit;
        // EnzoAcceleration = LengthUnits/TimeUnits/TimeUnits/pc*yr* yr/position_unit*time_unit*time_unit;
        EnzoAcceleration = EnzoVelocity / EnzoTime;
    }

#ifdef COMOVE
    if (ComovingCoordinates) {
        dadt_i /= EnzoTime;
        dadt_f /= EnzoTime;
        global_variable->a_i = a_i;
        global_variable->a_f = a_f;
        global_variable->dadt_i = dadt_i;
        global_variable->dadt_f = dadt_f;
    }
#endif


    global_variable->OldEnzoTimeStep    = global_variable->EnzoTimeStep;
    global_variable->EnzoTimeStep       = TimeStep * EnzoTime;

    std::cout << "ABYSS: Data transferred!" << std::endl;
    fprintf(stdout, "ABYSS: Data transferred!\n");
    EnzoCurrentTime = EnzoCurrentTime*EnzoTime;
    global_variable->EnzoCurrentTime = EnzoCurrentTime * 1e4; // in Myr unit


    std::cout << "Enzo Time    :" << global_variable->EnzoCurrentTime << " Myr" << std::endl;
    std::cerr << "Enzo Time    :" << global_variable->EnzoCurrentTime << " Myr" << std::endl;

    std::cout << "EnzoTimeStep :" << global_variable->EnzoTimeStep * 1e4 << " Myr" << std::endl;
    std::cerr << "EnzoTimeStep :" << global_variable->EnzoTimeStep * 1e4 << " Myr" << std::endl;
    // std::cerr << "Next EnzoTimeStep :" << global_variable->EnzoTimeStep * 1e4 << " Myr" << std::endl;
    // std::cout << "Nbody Mass    :" << particle[0]->Mass << std::endl;
    // std::cout << "Enzo  Mass    :" << Mass[0]*EnzoMass << std::endl;
    // std::cout << "enzo Time :" << TimeStep << std::endl;
    // std::cout << "nbody Time:" << EnzoTimeStep << std::endl;
    //  FixNumNeighbor = std::min((int)
    //  std::floor((NumberOfSingleParticle+NumberOfSingleParticleNew-1)/2.0),
    //  FixNumNeighbor0); // original by EW 2025.3.27




    /*------------------===-------------------------*/
    /********   Receive Particles to ABYSS  *********/
    /*----------------===---------------------------*/
    bool debug1 = true;

    /* Step 1: Gather sizes */
    int size = NumberOfProcessors+1; // the number of inter_comm processors
    int blank=0;
    int *recv_counts_old = new int[size];
    int *recv_counts_new = new int[size];

    MPI_Gather(&blank, 1, MPI_INT, recv_counts_old, 1, MPI_INT, NumberOfProcessors, inter_comm);
    MPI_Gather(&blank, 1, MPI_INT, recv_counts_new, 1, MPI_INT, NumberOfProcessors, inter_comm);


    /* Step 2: Compute displacements */
    int NumberOfSingleParticleOld, NumberOfSingleParticleNew;
    int* displs_old = new int[size];
    int* displs_new = new int[size];

    displs_old[0] = 0;
    displs_new[0] = 0;
    for (int i = 1; i < size; i++) {
        displs_old[i] = displs_old[i - 1] + recv_counts_old[i - 1];
        displs_new[i] = displs_new[i - 1] + recv_counts_new[i - 1];
    }
    NumberOfSingleParticleOld = displs_old[size - 1] + recv_counts_new[size - 1];
    NumberOfSingleParticleNew = displs_new[size - 1] + recv_counts_new[size - 1];


    if (debug1) {
        fprintf(nbpout, "displs_old= ");
        for (int i=0; i<size; i++) {
            fprintf(nbpout, "%d, ", displs_old[i]);
        }
        fprintf(nbpout, "\n NumberOfSingleParticleOld=%d\n",NumberOfSingleParticleOld);
        fflush(nbpout);
        fprintf(nbpout, "displs_new= ");
        for (int i=0; i<size; i++) {
            fprintf(nbpout, "%d, ", displs_new[i]);
        }
        fprintf(nbpout, "\n NumberOfSingleParticleNew=%d\n",NumberOfSingleParticleNew);
        fflush(nbpout);
    }



    /* Step 3: Allocate receive buffer on root */
    ParticleSendDataType    *recvbuf_old = new ParticleSendDataType[NumberOfSingleParticleOld];
    ParticleDataType        *recvbuf_new = new ParticleDataType[NumberOfSingleParticleNew];

    if (debug1) {
        fprintf(nbpout, "Receive Buffer Ready!\n");
        fflush(nbpout);
    }


    /* Step 4: Gatherv */
    MPI_Gatherv(NULL, 0, MPI_ENZO_PTCL_SEND, recvbuf_old, recv_counts_old, displs_old,
                MPI_ENZO_PTCL_SEND, NumberOfProcessors, inter_comm);
    MPI_Gatherv(NULL, 0, MPI_ENZO_PTCL, recvbuf_new, recv_counts_new, displs_new,
                MPI_ENZO_PTCL, NumberOfProcessors, inter_comm);
    if (debug1) {
        fprintf(nbpout, "Date received!\n");
        fflush(nbpout);
    }



    /* Prepare displs for all particles for sending out */
    recv_counts_old[0] +=  recv_counts_new[0]; // now recv_counts_old contains all particles 
    displs[0] = 0; // this is a given tbh
    for (int i = 1; i < size; i++) {
        displs[i] = displs[i - 1] + recv_counts_old[i - 1];
        recv_counts_old[i] += recv_counts_new[i];
    }
    delete[] recv_counts_old;
    delete[] recv_counts_new;




    if (debug1) {
        //fprintf(stderr, "ABYSS: ID = ");
        fprintf(nbpout, "ABYSS: ID for OLD = ");
        for (int i=0; i<NumberOfSingleParticleOld; i++) {
            //fprintf(stderr, "(%d, ", recvbuf[i].ID);
            //fprintf(stderr, "%.4e, %.4e), ", recvbuf[i].Position[0], recvbuf[i].BackgroundAcceleration[0]);
            fprintf(nbpout, "%d,  ", recvbuf_old[i].ID);
            //fprintf(nbpout, "%.4e, %.4e), ", recvbuf[i].Position[0], recvbuf[i].BackgroundAcceleration[0]);
        }
        //fprintf(stderr, "\n");
        fprintf(nbpout, ")\n");
        fflush(nbpout);

        fprintf(nbpout, "ABYSS: ID for NEW = ");
        for (int i=0; i<NumberOfSingleParticleNew; i++) {
            //fprintf(stderr, "(%d, ", recvbuf[i].ID);
            //fprintf(stderr, "%.4e, %.4e), ", recvbuf[i].Position[0], recvbuf[i].BackgroundAcceleration[0]);
            fprintf(nbpout, "%d,  ", recvbuf_new[i].ID);
            //fprintf(nbpout, "%.4e, %.4e), ", recvbuf[i].Position[0], recvbuf[i].BackgroundAcceleration[0]);
        }
        //fprintf(stderr, "\n");
        fprintf(nbpout, ")\n");
        fflush(nbpout);
    }

    // COM conversion
    // later on we might need to take mass weight into account
    // 1. F=ma; 2. F -> F_com; 3. F_com -> a_com
    ClusterAcceleration[0] = 0.0;
    ClusterAcceleration[1] = 0.0;
    ClusterAcceleration[2] = 0.0;

#ifdef COM_EVOLUTION // this code will not work by EW 2025.7.8
    double total_mass = 0.;
    if (NumberOfSingleParticle != 0) {
        for (int i = 0; i < NumberOfSingleParticle; i++) {
            for (int dim = 0; dim < Dim; dim++) {
                ClusterAcceleration[dim] += recvbuf[i].Mass * recvbuf[i].BackgroundAcceleration[dim];
            }
            total_mass += recvbuf[i].Mass;
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

    if (NumberOfSingleParticleNew != 0) {
        // we need to make adjustment to COM
        GetNewCenterOfMass(PID, newMass, newPosition, newVelocity,
                            NumberOfSingleParticleNew, ClusterPosition,
                            ClusterVelocity);

        for (int i = 0; i < NumberOfSingleParticleNew; i++) {
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

    /*-------------===-------------------------*/
    /********   Update Old Particles   *********/
    /*-----------===---------------------------*/
    // Update Existing Particles
#ifdef FEWBODY
    std::unordered_set<int> CMPtclsSet;
    CMPtclsSet.reserve(CMPtclWorker.size());
#endif
    Particle *ptcl;
    int EnzoProcessorNumber = 0;
    // loop for PID, going backwards to update the NextParticle
    for (int i = 0; i < NumberOfSingleParticleOld; i++) {
        int index = PIDtoIndexMap[recvbuf_old[i].ID];
        ptcl = &particles[index];

        while (i >= displs_old[EnzoProcessorNumber+1])
            EnzoProcessorNumber++;

#ifdef FEWBODY
        if (!ptcl->isActive && ptcl->CMPtclIndex != -1) {
            CMPtclsSet.insert(ptcl->CMPtclIndex);
        }
#endif

#ifdef COM_EVOLUTION
        for (int dim = 0; dim < Dim; dim++) {
            recvbuf[i].BackgroundAcceleration[dim] -= ClusterAcceleration[dim];
        }
#endif
        ptcl->update(recvbuf_old[i], EnzoProcessorNumber);
    } // endfor i
    delete[] recvbuf_old;

#ifdef FEWBODY
    for (int i : CMPtclsSet) { // Calculate background acceleration on the CM particles
        Particle *ptcl = &particles[i];
        double BackgroundAccelerationCM[Dim] = {0, 0, 0};

        for (int j = 0; j < ptcl->NumberOfMember; j++) {
            Particle *members = &particles[ptcl->Members[j]];

            for (int dim = 0; dim < Dim; dim++)
                BackgroundAccelerationCM[dim] += members->Mass * members->BackgroundAcceleration[dim];
        }
        for (int dim = 0; dim < Dim; dim++)
            ptcl->BackgroundAcceleration[dim] = BackgroundAccelerationCM[dim] / ptcl->Mass;
    }
#endif

    delete[] displs_old;


    /*-------------===-------------------------*/
    /********   Update New Particles   *********/
    /*-----------===---------------------------*/
    int index = -1;
    EnzoProcessorNumber = 0;
    if (NumberOfSingleParticleNew > 0) {
        for (int i = 0; i < NumberOfSingleParticleNew; i++) {

#ifdef COM_EVOLUTION
            recvbuf_new[i].reposition(ClusterPosition, ClusterVelocity, ClusterAcceleration);
#else
            recvbuf_new[i].reposition(ClusterPosition);
#endif

            /* decide to where to put new particles in the particles array*/
            // (Query) should it be better to have a function for it? if this is used in other routines as well
            if (NumberOfAvailableIndices != 0) {
                index = AvailableIndices[NumberOfAvailableIndices - 1];
                AvailableIndices[NumberOfAvailableIndices - 1] = -1;
                NumberOfAvailableIndices--;
            } else {
                index = global_variable->LastParticleIndex + 1;
                global_variable->LastParticleIndex++;
            }

            while (i >= displs_new[EnzoProcessorNumber+1])
                EnzoProcessorNumber++;

            particles[index].set(recvbuf_new[i], EnzoProcessorNumber);
            particles[index].ParticleIndex = index;
            PIDtoIndexMap.insert({recvbuf_new[i].ID, index});

            if (debug1) {
                recvbuf_new[i].print(1, 1, 1);
                particles[index].print(mass_unit, position_unit, velocity_unit);
            }

#ifdef SEVN
            if (recvbuf_new[i].CreationTime > 0.0 && recvbuf_new[i].Mass * EnzoMass * mass_unit > 2.2)
                initializeStellarEvolution(index);
#endif
        }

        fprintf(nbpout,"ABYSS: %d new paritcles loaded!", NumberOfSingleParticleNew);
        // This includes modification of regular force and irregular force
    } // end if NumberOfSingleParticleNew != 0
    delete[] recvbuf_new;

    delete[] displs_new;


    // (SEVN Query) After processing, wind & SN feedback,
    // 1. dm should be set to 0 - in update function
    // 2. If the particle is kicked, kicked velocity should be accounted - in
    // update function
    // 3. For PISN case, SEVN memory should be free - not yet


    // Update in case we update ClusterPosition. Currently, only COMEvolution update the ClusterPosition
	AbyssCenter[0] = ClusterPosition[0];
	AbyssCenter[1] = ClusterPosition[1];
	AbyssCenter[2] = ClusterPosition[2];


    NumberOfSingleParticle  += NumberOfSingleParticleNew;
    NumberOfParticle        += NumberOfSingleParticleNew;

    //  (Query) Do I need this?
    // fprintf(nbpout, "ABYSS    : Acceleration for particles on GPU.\n");
    /*
    if (NumberOfSingleParticle > 1) {
        CalculateAllAccelerationOnGPU(particle);
    }*/

    // UpdateNextRegTime(particle);
    //  fprintf(nbpout, "ABYSS    : Acceleration and neighbors are updated.\n");

    fprintf(nbpout,
            "ABYSS    : In ReceiveFromEnzo (after new particle might be added): \n");
    fprintf(nbpout,
            "ABYSS    : original NumberOfSingleParticle      = %d (+%d)\n",
            NumberOfSingleParticle - NumberOfSingleParticleNew,
            NumberOfSingleParticleNew);
    fprintf(nbpout, "ABYSS    : newly updated NumberOfSingleParticle = %d\n",
            NumberOfSingleParticleNew);
    // fprintf(nbpout, "ABYSS    : Particle size     = %d\n", particle.size());
    // fprintf(nbpout, "ABYSS    : NextRegTimeStep   = %.3e\n",
    // NextRegTimeBlock*global_variable->time_step); fprintf(nbpout, "ABYSS    :
    // NextRegTimeBlock  = %d\n", NextRegTimeBlock);
    // fprintf(nbpout, "ABYSS    : RegularList size  = %d\n",
    // RegularList.size());
    fprintf(nbpout, "ABYSS    : FixNumNeighbor    = %d\n", FixNumNeighbor);

    fprintf(stderr,
            "ABYSS    : In ReceiveFromEnzo (after new particle might be added): \n");
    fprintf(stderr,
            "ABYSS    : original NumberOfSingleParticle      = %d (+%d)\n",
            NumberOfSingleParticle - NumberOfSingleParticleNew,
            NumberOfSingleParticleNew);
    fprintf(stderr, "ABYSS    : newly updated NumberOfSingleParticle = %d\n",
            NumberOfSingleParticleNew);
    // fprintf(stderr, "ABYSS    : Particle size     = %d\n", particle.size());
    // fprintf(stderr, "ABYSS    : RegularList size = %d\n", RegularList.size());
    fprintf(nbpout, "ReceiveParticle Done!\n");
    fflush(stderr);
    fflush(nbpout);
    // fflush(gpuout);
    // fflush(binout);

    //MPI_Barrier(inter_comm); 
    return true;
}

/*
*
*
*
*
*
*
*
*/


// (Query) For binary mergers and SEVN, mass should be sent to Enzo as it was
// changed in Abyss by EW 2025.3.13 (Query) For SEVN, not only mass but also dm
// should be sent to Enzo, and it should be distributed into the grid by EW
// 2025.3.13

int SendParticleToEnzo(Worker *workers) {

    std::cout << "ABYSS: Starting SendParticleToEnzo..." << std::endl;
    std::cerr << "ABYSS: Starting SendParticleToEnzo..." << std::endl;
    /*if (NumberOfSingleParticle == 0) {
        std::cout << "ABYSS: Skipping SendParticleToEnzo..." << std::endl;
        std::cerr << "ABYSS: Skipping SendParticleToEnzo..." << std::endl;
        return 1; // SUCCESS -> 1 by EW 2025.3.11
    }*/
    MPI_Request request;
    MPI_Status status;

    fprintf(nbpout,
            "ABYSS: Starting SendParticleToEnzo...\n"
            "NumberOfSingleParticle=%d\n",
            NumberOfSingleParticle);
    //MPI_Barrier(inter_comm); 

    ParticleReceiveDataType *sendbuf = new ParticleReceiveDataType[NumberOfSingleParticle]; 

    int index;
    Particle *ptcl;

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

    /* Get NbodyCOM */
    double NbodyCOM[Dim] = {0, 0, 0};
    double mass = 0;
    for (int i = 0; i <= global_variable->LastParticleIndex; i++) {
        //index = PIDtoIndexMap[i];
        // maybe we can it more fancy.
        ptcl = &particles[i];
        if (!particles[i].isActive)
            continue;

#pragma unroll Dim
        for (int dim = 0; dim < Dim; dim++) {
            NbodyCOM[dim] += ptcl->Mass * ptcl->Position[dim];
        }
        mass += ptcl->Mass;
    }

#pragma unroll Dim
    for (int dim = 0; dim < Dim; dim++) {
        if (mass != 0.0)
            NbodyCOM[dim] /= mass;
        NbodyCOM[dim] /= EnzoLength;
        // NbodyCOM[dim] += ClusterPosition[dim];
    }
#endif

#ifdef FEWBODY
    std::unordered_set<int> CMPtclsSet;
    CMPtclsSet.reserve(CMPtclWorker.size());
    Particle *ptclCM;
#endif


    int offset;
    int *sendcounts  = new int[NumberOfProcessors+1];
    for (int i=0; i<NumberOfProcessors+1; i++) 
        sendcounts[i] = 0;

    /* Prepare Send Particle Buffer */
    for (int i = 0; i <= global_variable->LastParticleIndex; i++) {
        //index = PIDtoIndexMap[i];
        // maybe we can it more fancy.
        ptcl = &particles[i];
        if (ptcl->isCMptcl)
            continue;

        offset = displs[ptcl->EnzoProcessorNumber] + sendcounts[ptcl->EnzoProcessorNumber];

        fprintf(nbpout, "ID=%d, offset=%d, Processor=%d, sendcounts=%d, displs=%d\n",
                ptcl->PID, offset, ptcl->EnzoProcessorNumber, sendcounts[ptcl->EnzoProcessorNumber], displs[ptcl->EnzoProcessorNumber+1]);
        fflush(nbpout);

        // I can put everything under into a method of sendbuf (ParticleReceiveDataType)
        sendbuf[offset].ID = ptcl->PID;


#ifdef FEWBODY
        if (!ptcl->isActive) {
            for (int dim = 0; dim < Dim; dim++) {
                sendbuf[offset].Position[dim] = ptcl->Position[dim] / EnzoLength + ClusterPosition[dim];
#ifdef COM_EVOLUTION
                sendbuf[offset].Velocity[dim] = ptcl->Velocity[dim] / EnzoVelocity + ClusterVelocity[dim];
#else
                sendbuf[offset].Velocity[dim] = ptcl->Velocity[dim] / EnzoVelocity;
#endif
            }
#ifdef SEVN
            sendbuf[offset].InitialMass = ptcl->InitialMass; // This is already in Msun unit!!!
            sendbuf[offset].Mass = ptcl->Mass / EnzoMass;
            sendbuf[offset].WindEjectedMass = ptcl->dm * mass_unit;
            sendbuf[offset].SNEjectedMass = ptcl->SNEjectedMass * mass_unit;
            sendbuf[offset].Temperature = ptcl->T_eff;

            if (ptcl->dm > 0.0 || ptcl->SNEjectedMass > 0.0) {
                fprintf(stderr, "Feedback info send to Enzo...\n");
                fprintf(stderr,
                        "\tPID: %d. InitialMass: %e Msun, WindEjectedMass: %e Msun, "
                        "SNEjectedMass: %e Msun, T_eff: %e K\n",
                        ptcl->PID, sendbuf[offset].InitialMass, sendbuf[offset].WindEjectedMass,
                        sendbuf[offset].SNEjectedMass, sendbuf[offset].Temperature);
            }
#endif

            if (ptcl->CMPtclIndex != -1) {
                ptclCM = &particles[ptcl->CMPtclIndex];
                if (CMPtclsSet.find(ptcl->CMPtclIndex) == CMPtclsSet.end()) {

                    CMPtclsSet.insert(ptcl->CMPtclIndex);
                    ptclCM->NewNumberOfNeighbor = 0;
                }
                ptclCM->NewNeighbors[ptclCM->NewNumberOfNeighbor++] = offset;
            } else if (ptcl->Mass < 0.0) { // merger induced zero-mass particles, PISN case

                fprintf(stdout,
                        "In CommunicationToHydro... PID: %d should be removed!\n",
                        ptcl->PID);

                // /* // This is temporarilly commented out by EW 2025.4.13
                deleteParticle(ptcl->PID, i);
                NumberOfEscapeParticle++;
                // */
            } else {
                fprintf(stderr, "What's wrong? PID: %d\n", ptcl->PID);
                throw std::runtime_error("Why inactive particle in EnzoPIDs?");
            }
            continue;
        }
#endif // FewBody end // (Query) can you make it a particle method? hide under particle routine?
        r2 = 0;
        // this seems like repetitive
        for (int dim = 0; dim < Dim; dim++) {
            sendbuf[offset].Position[dim] = ptcl->Position[dim] / EnzoLength;
            sendbuf[offset].Velocity[dim] = ptcl->Velocity[dim] / EnzoVelocity;

#ifdef HUBBLE_FLOW
            // update ptcl->Position
            ptcl->Position[dim] -= NbodyCOM[dim] * EnzoLengh;
#endif

#ifdef COM_EVOLUTION
            if (IdentifyNbodyParticles && IdentifyOnTheFly)
                r2 += sendbuf[offset].Position[dim] * sendbuf[offset].Position[dim];
            // COM correction
            sendbuf[offset].Position[dim] += ClusterPosition[dim];
            sendbuf[offset].Velocity[dim] += ClusterVelocity[dim];
#else
            // calculate r2 first as NbodyCOM is not yet updated by ClusterPosition
            if (IdentifyNbodyParticles && IdentifyOnTheFly)
                r2 += (sendbuf[offset].Position[dim] - NbodyCOM[dim]) *
                    (sendbuf[offset].Position[dim] - NbodyCOM[dim]);
            // COM correction
            sendbuf[offset].Position[dim] += ClusterPosition[dim];
#endif
            if (IdentifyNbodyParticles && !IdentifyOnTheFly)
                r2 += (sendbuf[offset].Position[dim] - EnzoClusterPosition[dim]) *
                        (sendbuf[offset].Position[dim] - EnzoClusterPosition[dim]);
        } // for dim

        if (IdentifyNbodyParticles && ClusterRadius2 > 0 && r2 > ClusterRadius2) { // in Enzo Unit
            sendbuf[offset].Position[0] -= 20;
            deleteParticle(ptcl->PID, i);
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
#endif // SEVN
#endif // FEWBODY
            NumberOfEscapeParticle++;
        }
        // fprintf(stdout, "ABYSS: pid= %d, x=%e\n",ptcl->PID,Position[0][i]);

#ifdef SEVN
        sendbuf[offset].InitialMass = ptcl->InitialMass; // This is already in Msun unit!!!
        sendbuf[offset].Mass = ptcl->Mass / EnzoMass;
        sendbuf[offset].WindEjectedMass = ptcl->dm * mass_unit;
        sendbuf[offset].SNEjectedMass = ptcl->SNEjectedMass * mass_unit;
        sendbuf[offset].Temperature = ptcl->T_eff;

        if (ptcl->dm > 0.0 || ptcl->SNEjectedMass > 0.0) {
            fprintf(stderr, "Feedback info send to Enzo...\n");
            fprintf(stderr,
                    "\tPID: %d. InitialMass: %e Msun, WindEjectedMass: %e Msun, "
                    "SNEjectedMass: %e Msun, T_eff: %e K\n",
                    ptcl->PID, sendbuf[offset].InitialMass, sendbuf[offset].WindEjectedMass,
                    sendbuf[offset].SNEjectedMass, sendbuf[offset].Temperature);
        }
#endif
        ++sendcounts[ptcl->EnzoProcessorNumber];
    } //  loop over particles

#ifdef FEWBODY
    for (int i : CMPtclsSet) {
        ptcl = &particles[i];

#ifdef HUBBLE_FLOW
		// update ptcl->Position
		for (int dim=0; dim<Dim; dim++)
			ptcl->Position[dim] -= NbodyCOM[dim] * EnzoLength;
#endif

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
#ifdef HUBBLE_FLOW
				// update members->Position
				members->Position[dim] -= NbodyCOM[dim] * EnzoLength;
#endif
#ifdef COM_EVOLUTION
                if (IdentifyNbodyParticles && IdentifyOnTheFly)
                    r2 += memPosition[dim] * memPosition[dim];
                // COM correction
                memPosition[dim] += ClusterPosition[dim];
                memVelocity[dim] += ClusterVelocity[dim];
#else
                // calculate r2 first as NbodyCOM is not yet updated by ClusterPosition
                if (IdentifyNbodyParticles && IdentifyOnTheFly)
                    r2 += (memPosition[dim] - NbodyCOM[dim]) * (memPosition[dim] - NbodyCOM[dim]);
                // COM correction
                memPosition[dim] += ClusterPosition[dim];
#endif

                if (IdentifyNbodyParticles && !IdentifyOnTheFly)
                    r2 += (memPosition[dim] - EnzoClusterPosition[dim]) * (memPosition[dim] - EnzoClusterPosition[dim]);
            }

            if (IdentifyNbodyParticles && ClusterRadius2 > 0 && r2 <= ClusterRadius2) // in Enzo Unit
                memEscape = false;
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
                sendbuf[ptcl->NewNeighbors[j]].Position[0] -= 20;
                fprintf(stderr, "Binary escape... mem PID: %d, i: %d\n",
                        members->PID, ptcl->NewNeighbors[j]);
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
                PrevCMPtclWorker.insert({ptcl->ParticleIndex, CMPtclWorker[ptcl->ParticleIndex]});
            CMPtclWorker.erase(ptcl->ParticleIndex);
        }
    }
#endif

    std::cerr << "ABYSS: Sendbuf is ready!!" << std::endl;
    //std::cerr << "ABYSS: Waiting for Enzo to send data..." << std::endl;
    fprintf(nbpout, "ABYSS: Sendbuf is ready!!\n");


    /*-------------------------------------------*/
    /********   Send Particles to ENZO  *********/
    /*-------------------------------------------*/
    MPI_Scatterv(sendbuf, sendcounts, displs, MPI_ENZO_PTCL_RECV, NULL, 0,
                MPI_ENZO_PTCL_RECV, NumberOfProcessors, inter_comm);
    delete [] sendcounts;
    delete [] sendbuf;

    if (IdentifyOnTheFly) { // (Query AEOS) EW changed the previous condition since enzo processors cannot know the total number of nbody particles 2025.7.20
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
#ifdef HUBBLE_FLOW
        // for hubble flow is on, we update ClusterPosition, and shift Abyss particles
        for (int dim = 0; dim < Dim; dim++) {
            ClusterPosition[dim] += NbodyCOM[dim];
        }
        MPI_Send(ClusterPosition, 3, MPI_DOUBLE, 0, 700, inter_comm);
#else // NO HUBBLE_FLOW
        for (int dim = 0; dim < Dim; dim++) {
            NbodyCOM[dim] += ClusterPosition[dim];
        }
        // Synchronize Cluster Posiition
        MPI_Ibcast(NbodyCOM, 3, MPI_DOUBLE, NumberOfProcessors, inter_comm, &request);
        fprintf(nbpout, "ABYSS: Cluster position broadcating.\n");
#endif
#endif
    }
    fprintf(stderr, "ABYSS: Data sent!\n");

    // Update AbyssCenter. We need to update this for Hubble_FLOW
	AbyssCenter[0] = ClusterPosition[0];
	AbyssCenter[1] = ClusterPosition[1];
	AbyssCenter[2] = ClusterPosition[2];


    NumberOfSingleParticle -= NumberOfEscapeParticle;
#ifdef FEWBODY
    // NumberOfParticle -= NumberOfEscapeParticle; // (Query) EW: We don't have to
    // do this.
#else
    NumberOfParticle -= NumberOfEscapeParticle;
#endif

    fprintf(stderr,
            "ABYSS    : In SendToEnzo (after particle might be escaped): \n");
    fprintf(stderr, "ABYSS    : original NumberOfSingleParticle      = %d (-%d)\n",
            NumberOfSingleParticle + NumberOfEscapeParticle, NumberOfEscapeParticle);
    fprintf(stderr, "ABYSS    : newly updated NumberOfSingleParticle = %d\n",
            NumberOfSingleParticle);
    fprintf(stderr, "ABYSS    : newly updated NumberOfParticle = %d\n",
            NumberOfParticle);

    fprintf(nbpout,
            "ABYSS    : In SendToEnzo (after particle might be escaped): \n");
    fprintf(nbpout, "ABYSS    : original NumberOfSingleParticle      = %d (-%d)\n",
            NumberOfSingleParticle + NumberOfEscapeParticle, NumberOfEscapeParticle);
    fprintf(nbpout, "ABYSS    : newly updated NumberOfSingleParticle = %d\n",
            NumberOfSingleParticle);

    fflush(stderr);
    if (NumberOfSingleParticle != 0 && IdentifyOnTheFly) {
        MPI_Wait(&request, MPI_STATUS_IGNORE);
    }
    // fflush(gpuout);
    //  fflush(binout);
    fprintf(nbpout, "ABYSS: Sending data done!\n");
    fflush(nbpout);
    std::cout << "ABYSS: Sending data done!" << std::endl;
    std::cerr << "ABYSS: Sending data done!" << std::endl;
    //MPI_Barrier(inter_comm); 
    return true;
}


/*
*
*
*
*
*
*
*
*/

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