#ifndef __communicators_H
#define __communicators_H
#ifdef USE_MPI
#include "mpi.h"
/* by YS, MPI COMMs*/
extern MPI_Comm enzo_comm;
extern MPI_Comm abyss_comm;
extern MPI_Comm inter_comm;
extern MPI_Comm local_comm;
#endif /* USE_MPI */
#endif