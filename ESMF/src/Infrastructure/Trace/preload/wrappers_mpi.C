
/**
 *  wrappers_mpi.C
 *
 *  !! DO NOT EDIT !!
 *  This file was generated by gen_wrappers_mpi.py.
 *  !! DO NOT EDIT !!
 *
 *  Wrappers for MPI calls so we can trace them.
 *  These are linked statically into the executable using
 *  --wrap=SYMBOL flags to the linker.
 *
 */

#ifndef ESMF_MPIUNI
#include <mpi.h>

#include "ESMCI_Macros.h"
#include "ESMCI_Trace.h"
#include "preload.h"
#include "wrappers_mpi.h"

extern "C" {

  static int insideMPIRegion = 0;
  static int ignorerc = 0;
  
  /*
   * C MPI functions
   */

  
    extern int __real_MPI_Allgather(ESMF_MPI_CONST void *sendbuf, int  sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm);    

    int __wrap_MPI_Allgather(ESMF_MPI_CONST void *sendbuf, int  sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_MPI_Allgather (C)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("MPI_Allgather", &ignorerc);
        int ret = __real_MPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
        ESMCI::TraceEventRegionExit("MPI_Allgather", &ignorerc);
        insideMPIRegion = 0;
        return ret;
      }
      else {
        return __real_MPI_Allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
      }
    }
  
    extern int __real_MPI_Allgatherv(ESMF_MPI_CONST void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, ESMF_MPI_CONST int recvcounts[], ESMF_MPI_CONST int displs[], MPI_Datatype recvtype, MPI_Comm comm);    

    int __wrap_MPI_Allgatherv(ESMF_MPI_CONST void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, ESMF_MPI_CONST int recvcounts[], ESMF_MPI_CONST int displs[], MPI_Datatype recvtype, MPI_Comm comm) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_MPI_Allgatherv (C)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("MPI_Allgatherv", &ignorerc);
        int ret = __real_MPI_Allgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm);
        ESMCI::TraceEventRegionExit("MPI_Allgatherv", &ignorerc);
        insideMPIRegion = 0;
        return ret;
      }
      else {
        return __real_MPI_Allgatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, comm);
      }
    }
  
    extern int __real_MPI_Allreduce(ESMF_MPI_CONST void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);    

    int __wrap_MPI_Allreduce(ESMF_MPI_CONST void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_MPI_Allreduce (C)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("MPI_Allreduce", &ignorerc);
        int ret = __real_MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
        ESMCI::TraceEventRegionExit("MPI_Allreduce", &ignorerc);
        insideMPIRegion = 0;
        return ret;
      }
      else {
        return __real_MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
      }
    }
  
    extern int __real_MPI_Alltoall(ESMF_MPI_CONST void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm);    

    int __wrap_MPI_Alltoall(ESMF_MPI_CONST void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_MPI_Alltoall (C)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("MPI_Alltoall", &ignorerc);
        int ret = __real_MPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
        ESMCI::TraceEventRegionExit("MPI_Alltoall", &ignorerc);
        insideMPIRegion = 0;
        return ret;
      }
      else {
        return __real_MPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);
      }
    }
  
    extern int __real_MPI_Alltoallv(ESMF_MPI_CONST void *sendbuf, ESMF_MPI_CONST int sendcounts[], ESMF_MPI_CONST int sdispls[], MPI_Datatype sendtype, void *recvbuf, ESMF_MPI_CONST int recvcounts[], ESMF_MPI_CONST int rdispls[], MPI_Datatype recvtype, MPI_Comm comm);    

    int __wrap_MPI_Alltoallv(ESMF_MPI_CONST void *sendbuf, ESMF_MPI_CONST int sendcounts[], ESMF_MPI_CONST int sdispls[], MPI_Datatype sendtype, void *recvbuf, ESMF_MPI_CONST int recvcounts[], ESMF_MPI_CONST int rdispls[], MPI_Datatype recvtype, MPI_Comm comm) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_MPI_Alltoallv (C)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("MPI_Alltoallv", &ignorerc);
        int ret = __real_MPI_Alltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts, rdispls, recvtype, comm);
        ESMCI::TraceEventRegionExit("MPI_Alltoallv", &ignorerc);
        insideMPIRegion = 0;
        return ret;
      }
      else {
        return __real_MPI_Alltoallv(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts, rdispls, recvtype, comm);
      }
    }
  
    extern int __real_MPI_Alltoallw(ESMF_MPI_CONST void *sendbuf, ESMF_MPI_CONST int sendcounts[], ESMF_MPI_CONST int sdispls[], ESMF_MPI_CONST MPI_Datatype sendtypes[], void *recvbuf, ESMF_MPI_CONST int recvcounts[], ESMF_MPI_CONST int rdispls[], ESMF_MPI_CONST MPI_Datatype recvtypes[], MPI_Comm comm);    

    int __wrap_MPI_Alltoallw(ESMF_MPI_CONST void *sendbuf, ESMF_MPI_CONST int sendcounts[], ESMF_MPI_CONST int sdispls[], ESMF_MPI_CONST MPI_Datatype sendtypes[], void *recvbuf, ESMF_MPI_CONST int recvcounts[], ESMF_MPI_CONST int rdispls[], ESMF_MPI_CONST MPI_Datatype recvtypes[], MPI_Comm comm) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_MPI_Alltoallw (C)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("MPI_Alltoallw", &ignorerc);
        int ret = __real_MPI_Alltoallw(sendbuf, sendcounts, sdispls, sendtypes, recvbuf, recvcounts, rdispls, recvtypes, comm);
        ESMCI::TraceEventRegionExit("MPI_Alltoallw", &ignorerc);
        insideMPIRegion = 0;
        return ret;
      }
      else {
        return __real_MPI_Alltoallw(sendbuf, sendcounts, sdispls, sendtypes, recvbuf, recvcounts, rdispls, recvtypes, comm);
      }
    }
  
    extern int __real_MPI_Barrier(MPI_Comm comm);    

    int __wrap_MPI_Barrier(MPI_Comm comm) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_MPI_Barrier (C)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("MPI_Barrier", &ignorerc);
        int ret = __real_MPI_Barrier(comm);
        ESMCI::TraceEventRegionExit("MPI_Barrier", &ignorerc);
        insideMPIRegion = 0;
        return ret;
      }
      else {
        return __real_MPI_Barrier(comm);
      }
    }
  
    extern int __real_MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm);    

    int __wrap_MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_MPI_Bcast (C)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("MPI_Bcast", &ignorerc);
        int ret = __real_MPI_Bcast(buffer, count, datatype, root, comm);
        ESMCI::TraceEventRegionExit("MPI_Bcast", &ignorerc);
        insideMPIRegion = 0;
        return ret;
      }
      else {
        return __real_MPI_Bcast(buffer, count, datatype, root, comm);
      }
    }
  
    extern int __real_MPI_Gather(ESMF_MPI_CONST void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);    

    int __wrap_MPI_Gather(ESMF_MPI_CONST void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_MPI_Gather (C)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("MPI_Gather", &ignorerc);
        int ret = __real_MPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
        ESMCI::TraceEventRegionExit("MPI_Gather", &ignorerc);
        insideMPIRegion = 0;
        return ret;
      }
      else {
        return __real_MPI_Gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
      }
    }
  
    extern int __real_MPI_Gatherv(ESMF_MPI_CONST void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, ESMF_MPI_CONST int recvcounts[], ESMF_MPI_CONST int displs[], MPI_Datatype recvtype, int root, MPI_Comm comm);    

    int __wrap_MPI_Gatherv(ESMF_MPI_CONST void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, ESMF_MPI_CONST int recvcounts[], ESMF_MPI_CONST int displs[], MPI_Datatype recvtype, int root, MPI_Comm comm) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_MPI_Gatherv (C)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("MPI_Gatherv", &ignorerc);
        int ret = __real_MPI_Gatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, root, comm);
        ESMCI::TraceEventRegionExit("MPI_Gatherv", &ignorerc);
        insideMPIRegion = 0;
        return ret;
      }
      else {
        return __real_MPI_Gatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, root, comm);
      }
    }
  
    extern int __real_MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status);    

    int __wrap_MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_MPI_Recv (C)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("MPI_Recv", &ignorerc);
        int ret = __real_MPI_Recv(buf, count, datatype, source, tag, comm, status);
        ESMCI::TraceEventRegionExit("MPI_Recv", &ignorerc);
        insideMPIRegion = 0;
        return ret;
      }
      else {
        return __real_MPI_Recv(buf, count, datatype, source, tag, comm, status);
      }
    }
  
    extern int __real_MPI_Reduce(ESMF_MPI_CONST void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm);    

    int __wrap_MPI_Reduce(ESMF_MPI_CONST void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_MPI_Reduce (C)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("MPI_Reduce", &ignorerc);
        int ret = __real_MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
        ESMCI::TraceEventRegionExit("MPI_Reduce", &ignorerc);
        insideMPIRegion = 0;
        return ret;
      }
      else {
        return __real_MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
      }
    }
  
    extern int __real_MPI_Scatter(ESMF_MPI_CONST void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);    

    int __wrap_MPI_Scatter(ESMF_MPI_CONST void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_MPI_Scatter (C)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("MPI_Scatter", &ignorerc);
        int ret = __real_MPI_Scatter(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
        ESMCI::TraceEventRegionExit("MPI_Scatter", &ignorerc);
        insideMPIRegion = 0;
        return ret;
      }
      else {
        return __real_MPI_Scatter(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm);
      }
    }
  
    extern int __real_MPI_Send(ESMF_MPI_CONST void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);    

    int __wrap_MPI_Send(ESMF_MPI_CONST void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_MPI_Send (C)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("MPI_Send", &ignorerc);
        int ret = __real_MPI_Send(buf, count, datatype, dest, tag, comm);
        ESMCI::TraceEventRegionExit("MPI_Send", &ignorerc);
        insideMPIRegion = 0;
        return ret;
      }
      else {
        return __real_MPI_Send(buf, count, datatype, dest, tag, comm);
      }
    }
  
    extern int __real_MPI_Sendrecv(ESMF_MPI_CONST void *sendbuf, int sendcount, MPI_Datatype sendtype, int dest, int sendtag, void *recvbuf, int recvcount, MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status *status);    

    int __wrap_MPI_Sendrecv(ESMF_MPI_CONST void *sendbuf, int sendcount, MPI_Datatype sendtype, int dest, int sendtag, void *recvbuf, int recvcount, MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status *status) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_MPI_Sendrecv (C)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("MPI_Sendrecv", &ignorerc);
        int ret = __real_MPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, recvcount, recvtype, source, recvtag, comm, status);
        ESMCI::TraceEventRegionExit("MPI_Sendrecv", &ignorerc);
        insideMPIRegion = 0;
        return ret;
      }
      else {
        return __real_MPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, recvcount, recvtype, source, recvtag, comm, status);
      }
    }
  
    extern int __real_MPI_Wait(MPI_Request *request, MPI_Status *status);    

    int __wrap_MPI_Wait(MPI_Request *request, MPI_Status *status) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_MPI_Wait (C)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("MPI_Wait", &ignorerc);
        int ret = __real_MPI_Wait(request, status);
        ESMCI::TraceEventRegionExit("MPI_Wait", &ignorerc);
        insideMPIRegion = 0;
        return ret;
      }
      else {
        return __real_MPI_Wait(request, status);
      }
    }
  
    extern int __real_MPI_Waitall(int count, MPI_Request array_of_requests[], MPI_Status *array_of_statuses);    

    int __wrap_MPI_Waitall(int count, MPI_Request array_of_requests[], MPI_Status *array_of_statuses) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_MPI_Waitall (C)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("MPI_Waitall", &ignorerc);
        int ret = __real_MPI_Waitall(count, array_of_requests, array_of_statuses);
        ESMCI::TraceEventRegionExit("MPI_Waitall", &ignorerc);
        insideMPIRegion = 0;
        return ret;
      }
      else {
        return __real_MPI_Waitall(count, array_of_requests, array_of_statuses);
      }
    }
  
    extern int __real_MPI_Waitany(int count, MPI_Request array_of_requests[], int *index, MPI_Status *status);    

    int __wrap_MPI_Waitany(int count, MPI_Request array_of_requests[], int *index, MPI_Status *status) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_MPI_Waitany (C)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("MPI_Waitany", &ignorerc);
        int ret = __real_MPI_Waitany(count, array_of_requests, index, status);
        ESMCI::TraceEventRegionExit("MPI_Waitany", &ignorerc);
        insideMPIRegion = 0;
        return ret;
      }
      else {
        return __real_MPI_Waitany(count, array_of_requests, index, status);
      }
    }
  
    extern int __real_MPI_Waitsome(int incount, MPI_Request array_of_requests[], int *outcount, int array_of_indices[], MPI_Status array_of_statuses[]);    

    int __wrap_MPI_Waitsome(int incount, MPI_Request array_of_requests[], int *outcount, int array_of_indices[], MPI_Status array_of_statuses[]) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_MPI_Waitsome (C)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("MPI_Waitsome", &ignorerc);
        int ret = __real_MPI_Waitsome(incount, array_of_requests, outcount, array_of_indices, array_of_statuses);
        ESMCI::TraceEventRegionExit("MPI_Waitsome", &ignorerc);
        insideMPIRegion = 0;
        return ret;
      }
      else {
        return __real_MPI_Waitsome(incount, array_of_requests, outcount, array_of_indices, array_of_statuses);
      }
    }
  

  /*
   * Fortran MPI functions
   */

  
    extern void FTN_X(__real_mpi_allgather)(MPI_Fint *sendbuf, MPI_Fint *sendcount, MPI_Fint *sendtype, MPI_Fint *recvbuf, MPI_Fint *recvcount, MPI_Fint *recvtype, MPI_Fint *comm, MPI_Fint *ierr);    

    void FTN_X(__wrap_mpi_allgather)(MPI_Fint *sendbuf, MPI_Fint *sendcount, MPI_Fint *sendtype, MPI_Fint *recvbuf, MPI_Fint *recvcount, MPI_Fint *recvtype, MPI_Fint *comm, MPI_Fint *ierr) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_mpi_allgather_ (Fortran)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("mpi_allgather", &ignorerc);
        FTN_X(__real_mpi_allgather)(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, ierr);
        ESMCI::TraceEventRegionExit("mpi_allgather", &ignorerc);
        insideMPIRegion = 0;
      }
      else {
        FTN_X(__real_mpi_allgather)(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, ierr);
      }
    }
  
    extern void FTN_X(__real_mpi_allgatherv)(MPI_Fint *sendbuf, MPI_Fint *sendcount, MPI_Fint *sendtype, MPI_Fint *recvbuf, MPI_Fint *recvcount, MPI_Fint *displs, MPI_Fint *recvtype, MPI_Fint *comm, MPI_Fint *ierr);    

    void FTN_X(__wrap_mpi_allgatherv)(MPI_Fint *sendbuf, MPI_Fint *sendcount, MPI_Fint *sendtype, MPI_Fint *recvbuf, MPI_Fint *recvcount, MPI_Fint *displs, MPI_Fint *recvtype, MPI_Fint *comm, MPI_Fint *ierr) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_mpi_allgatherv_ (Fortran)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("mpi_allgatherv", &ignorerc);
        FTN_X(__real_mpi_allgatherv)(sendbuf, sendcount, sendtype, recvbuf, recvcount, displs, recvtype, comm, ierr);
        ESMCI::TraceEventRegionExit("mpi_allgatherv", &ignorerc);
        insideMPIRegion = 0;
      }
      else {
        FTN_X(__real_mpi_allgatherv)(sendbuf, sendcount, sendtype, recvbuf, recvcount, displs, recvtype, comm, ierr);
      }
    }
  
    extern void FTN_X(__real_mpi_allreduce)(MPI_Fint *sendbuf, MPI_Fint *recvbuf, MPI_Fint *count, MPI_Fint *datatype, MPI_Fint *op, MPI_Fint *comm, MPI_Fint *ierr);    

    void FTN_X(__wrap_mpi_allreduce)(MPI_Fint *sendbuf, MPI_Fint *recvbuf, MPI_Fint *count, MPI_Fint *datatype, MPI_Fint *op, MPI_Fint *comm, MPI_Fint *ierr) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_mpi_allreduce_ (Fortran)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("mpi_allreduce", &ignorerc);
        FTN_X(__real_mpi_allreduce)(sendbuf, recvbuf, count, datatype, op, comm, ierr);
        ESMCI::TraceEventRegionExit("mpi_allreduce", &ignorerc);
        insideMPIRegion = 0;
      }
      else {
        FTN_X(__real_mpi_allreduce)(sendbuf, recvbuf, count, datatype, op, comm, ierr);
      }
    }
  
    extern void FTN_X(__real_mpi_alltoall)(MPI_Fint *sendbuf, MPI_Fint *sendcount, MPI_Fint *sendtype, MPI_Fint *recvbuf, MPI_Fint *recvcount, MPI_Fint *recvtype, MPI_Fint *comm, MPI_Fint *ierr);    

    void FTN_X(__wrap_mpi_alltoall)(MPI_Fint *sendbuf, MPI_Fint *sendcount, MPI_Fint *sendtype, MPI_Fint *recvbuf, MPI_Fint *recvcount, MPI_Fint *recvtype, MPI_Fint *comm, MPI_Fint *ierr) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_mpi_alltoall_ (Fortran)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("mpi_alltoall", &ignorerc);
        FTN_X(__real_mpi_alltoall)(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, ierr);
        ESMCI::TraceEventRegionExit("mpi_alltoall", &ignorerc);
        insideMPIRegion = 0;
      }
      else {
        FTN_X(__real_mpi_alltoall)(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, ierr);
      }
    }
  
    extern void FTN_X(__real_mpi_alltoallv)(MPI_Fint *sendbuf, MPI_Fint *sendcounts, MPI_Fint *sdispls, MPI_Fint *sendtype, MPI_Fint *recvbuf, MPI_Fint *recvcounts, MPI_Fint *rdispls, MPI_Fint *recvtype, MPI_Fint *comm, MPI_Fint *ierr);    

    void FTN_X(__wrap_mpi_alltoallv)(MPI_Fint *sendbuf, MPI_Fint *sendcounts, MPI_Fint *sdispls, MPI_Fint *sendtype, MPI_Fint *recvbuf, MPI_Fint *recvcounts, MPI_Fint *rdispls, MPI_Fint *recvtype, MPI_Fint *comm, MPI_Fint *ierr) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_mpi_alltoallv_ (Fortran)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("mpi_alltoallv", &ignorerc);
        FTN_X(__real_mpi_alltoallv)(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts, rdispls, recvtype, comm, ierr);
        ESMCI::TraceEventRegionExit("mpi_alltoallv", &ignorerc);
        insideMPIRegion = 0;
      }
      else {
        FTN_X(__real_mpi_alltoallv)(sendbuf, sendcounts, sdispls, sendtype, recvbuf, recvcounts, rdispls, recvtype, comm, ierr);
      }
    }
  
    extern void FTN_X(__real_mpi_alltoallw)(MPI_Fint *sendbuf, MPI_Fint *sendcounts, MPI_Fint *sdispls, MPI_Fint *sendtypes, MPI_Fint *recvbuf, MPI_Fint *recvcounts, MPI_Fint *rdispls, MPI_Fint *recvtypes, MPI_Fint *comm, MPI_Fint *ierr);    

    void FTN_X(__wrap_mpi_alltoallw)(MPI_Fint *sendbuf, MPI_Fint *sendcounts, MPI_Fint *sdispls, MPI_Fint *sendtypes, MPI_Fint *recvbuf, MPI_Fint *recvcounts, MPI_Fint *rdispls, MPI_Fint *recvtypes, MPI_Fint *comm, MPI_Fint *ierr) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_mpi_alltoallw_ (Fortran)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("mpi_alltoallw", &ignorerc);
        FTN_X(__real_mpi_alltoallw)(sendbuf, sendcounts, sdispls, sendtypes, recvbuf, recvcounts, rdispls, recvtypes, comm, ierr);
        ESMCI::TraceEventRegionExit("mpi_alltoallw", &ignorerc);
        insideMPIRegion = 0;
      }
      else {
        FTN_X(__real_mpi_alltoallw)(sendbuf, sendcounts, sdispls, sendtypes, recvbuf, recvcounts, rdispls, recvtypes, comm, ierr);
      }
    }
  
    extern void FTN_X(__real_mpi_barrier)(MPI_Fint *comm, MPI_Fint *ierr);    

    void FTN_X(__wrap_mpi_barrier)(MPI_Fint *comm, MPI_Fint *ierr) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_mpi_barrier_ (Fortran)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("mpi_barrier", &ignorerc);
        FTN_X(__real_mpi_barrier)(comm, ierr);
        ESMCI::TraceEventRegionExit("mpi_barrier", &ignorerc);
        insideMPIRegion = 0;
      }
      else {
        FTN_X(__real_mpi_barrier)(comm, ierr);
      }
    }
  
    extern void FTN_X(__real_mpi_bcast)(MPI_Fint *buffer, MPI_Fint *count, MPI_Fint *datatype, MPI_Fint *root, MPI_Fint *comm, MPI_Fint *ierr);    

    void FTN_X(__wrap_mpi_bcast)(MPI_Fint *buffer, MPI_Fint *count, MPI_Fint *datatype, MPI_Fint *root, MPI_Fint *comm, MPI_Fint *ierr) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_mpi_bcast_ (Fortran)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("mpi_bcast", &ignorerc);
        FTN_X(__real_mpi_bcast)(buffer, count, datatype, root, comm, ierr);
        ESMCI::TraceEventRegionExit("mpi_bcast", &ignorerc);
        insideMPIRegion = 0;
      }
      else {
        FTN_X(__real_mpi_bcast)(buffer, count, datatype, root, comm, ierr);
      }
    }
  
    extern void FTN_X(__real_mpi_exscan)(MPI_Fint *sendbuf, MPI_Fint *recvbuf, MPI_Fint *count, MPI_Fint *datatype, MPI_Fint *op, MPI_Fint *comm, MPI_Fint *ierr);    

    void FTN_X(__wrap_mpi_exscan)(MPI_Fint *sendbuf, MPI_Fint *recvbuf, MPI_Fint *count, MPI_Fint *datatype, MPI_Fint *op, MPI_Fint *comm, MPI_Fint *ierr) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_mpi_exscan_ (Fortran)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("mpi_exscan", &ignorerc);
        FTN_X(__real_mpi_exscan)(sendbuf, recvbuf, count, datatype, op, comm, ierr);
        ESMCI::TraceEventRegionExit("mpi_exscan", &ignorerc);
        insideMPIRegion = 0;
      }
      else {
        FTN_X(__real_mpi_exscan)(sendbuf, recvbuf, count, datatype, op, comm, ierr);
      }
    }
  
    extern void FTN_X(__real_mpi_gather)(MPI_Fint *sendbuf, MPI_Fint *sendcount, MPI_Fint *sendtype, MPI_Fint *recvbuf, MPI_Fint *recvcount, MPI_Fint *recvtype, MPI_Fint *root, MPI_Fint *comm, MPI_Fint *ierr);    

    void FTN_X(__wrap_mpi_gather)(MPI_Fint *sendbuf, MPI_Fint *sendcount, MPI_Fint *sendtype, MPI_Fint *recvbuf, MPI_Fint *recvcount, MPI_Fint *recvtype, MPI_Fint *root, MPI_Fint *comm, MPI_Fint *ierr) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_mpi_gather_ (Fortran)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("mpi_gather", &ignorerc);
        FTN_X(__real_mpi_gather)(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, ierr);
        ESMCI::TraceEventRegionExit("mpi_gather", &ignorerc);
        insideMPIRegion = 0;
      }
      else {
        FTN_X(__real_mpi_gather)(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, ierr);
      }
    }
  
    extern void FTN_X(__real_mpi_gatherv)(MPI_Fint *sendbuf, MPI_Fint *sendcount, MPI_Fint *sendtype, MPI_Fint *recvbuf, MPI_Fint *recvcounts, MPI_Fint *displs, MPI_Fint *recvtype, MPI_Fint *root, MPI_Fint *comm, MPI_Fint *ierr);    

    void FTN_X(__wrap_mpi_gatherv)(MPI_Fint *sendbuf, MPI_Fint *sendcount, MPI_Fint *sendtype, MPI_Fint *recvbuf, MPI_Fint *recvcounts, MPI_Fint *displs, MPI_Fint *recvtype, MPI_Fint *root, MPI_Fint *comm, MPI_Fint *ierr) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_mpi_gatherv_ (Fortran)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("mpi_gatherv", &ignorerc);
        FTN_X(__real_mpi_gatherv)(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, root, comm, ierr);
        ESMCI::TraceEventRegionExit("mpi_gatherv", &ignorerc);
        insideMPIRegion = 0;
      }
      else {
        FTN_X(__real_mpi_gatherv)(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs, recvtype, root, comm, ierr);
      }
    }
  
    extern void FTN_X(__real_mpi_recv)(MPI_Fint *buf, MPI_Fint *count, MPI_Fint *datatype, MPI_Fint *source, MPI_Fint *tag, MPI_Fint *comm, MPI_Fint *status, MPI_Fint *ierr);    

    void FTN_X(__wrap_mpi_recv)(MPI_Fint *buf, MPI_Fint *count, MPI_Fint *datatype, MPI_Fint *source, MPI_Fint *tag, MPI_Fint *comm, MPI_Fint *status, MPI_Fint *ierr) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_mpi_recv_ (Fortran)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("mpi_recv", &ignorerc);
        FTN_X(__real_mpi_recv)(buf, count, datatype, source, tag, comm, status, ierr);
        ESMCI::TraceEventRegionExit("mpi_recv", &ignorerc);
        insideMPIRegion = 0;
      }
      else {
        FTN_X(__real_mpi_recv)(buf, count, datatype, source, tag, comm, status, ierr);
      }
    }
  
    extern void FTN_X(__real_mpi_reduce)(MPI_Fint *sendbuf, MPI_Fint *recvbuf, MPI_Fint *count, MPI_Fint *datatype, MPI_Fint *op, MPI_Fint *root, MPI_Fint *comm, MPI_Fint *ierr);    

    void FTN_X(__wrap_mpi_reduce)(MPI_Fint *sendbuf, MPI_Fint *recvbuf, MPI_Fint *count, MPI_Fint *datatype, MPI_Fint *op, MPI_Fint *root, MPI_Fint *comm, MPI_Fint *ierr) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_mpi_reduce_ (Fortran)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("mpi_reduce", &ignorerc);
        FTN_X(__real_mpi_reduce)(sendbuf, recvbuf, count, datatype, op, root, comm, ierr);
        ESMCI::TraceEventRegionExit("mpi_reduce", &ignorerc);
        insideMPIRegion = 0;
      }
      else {
        FTN_X(__real_mpi_reduce)(sendbuf, recvbuf, count, datatype, op, root, comm, ierr);
      }
    }
  
    extern void FTN_X(__real_mpi_reduce_scatter)(MPI_Fint *sendbuf, MPI_Fint *recvbuf, MPI_Fint *recvcounts, MPI_Fint *datatype, MPI_Fint *op, MPI_Fint *comm, MPI_Fint *ierr);    

    void FTN_X(__wrap_mpi_reduce_scatter)(MPI_Fint *sendbuf, MPI_Fint *recvbuf, MPI_Fint *recvcounts, MPI_Fint *datatype, MPI_Fint *op, MPI_Fint *comm, MPI_Fint *ierr) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_mpi_reduce_scatter_ (Fortran)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("mpi_reduce_scatter", &ignorerc);
        FTN_X(__real_mpi_reduce_scatter)(sendbuf, recvbuf, recvcounts, datatype, op, comm, ierr);
        ESMCI::TraceEventRegionExit("mpi_reduce_scatter", &ignorerc);
        insideMPIRegion = 0;
      }
      else {
        FTN_X(__real_mpi_reduce_scatter)(sendbuf, recvbuf, recvcounts, datatype, op, comm, ierr);
      }
    }
  
    extern void FTN_X(__real_mpi_scatter)(MPI_Fint *sendbuf, MPI_Fint *sendcount, MPI_Fint *sendtype, MPI_Fint *recvbuf, MPI_Fint *recvcount, MPI_Fint *recvtype, MPI_Fint *root, MPI_Fint *comm, MPI_Fint *ierr);    

    void FTN_X(__wrap_mpi_scatter)(MPI_Fint *sendbuf, MPI_Fint *sendcount, MPI_Fint *sendtype, MPI_Fint *recvbuf, MPI_Fint *recvcount, MPI_Fint *recvtype, MPI_Fint *root, MPI_Fint *comm, MPI_Fint *ierr) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_mpi_scatter_ (Fortran)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("mpi_scatter", &ignorerc);
        FTN_X(__real_mpi_scatter)(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, ierr);
        ESMCI::TraceEventRegionExit("mpi_scatter", &ignorerc);
        insideMPIRegion = 0;
      }
      else {
        FTN_X(__real_mpi_scatter)(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, ierr);
      }
    }
  
    extern void FTN_X(__real_mpi_scatterv)(MPI_Fint *sendbuf, MPI_Fint *sendcounts, MPI_Fint *displs, MPI_Fint *sendtype, MPI_Fint *recvbuf, MPI_Fint *recvcount, MPI_Fint *recvtype, MPI_Fint *root, MPI_Fint *comm, MPI_Fint *ierr);    

    void FTN_X(__wrap_mpi_scatterv)(MPI_Fint *sendbuf, MPI_Fint *sendcounts, MPI_Fint *displs, MPI_Fint *sendtype, MPI_Fint *recvbuf, MPI_Fint *recvcount, MPI_Fint *recvtype, MPI_Fint *root, MPI_Fint *comm, MPI_Fint *ierr) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_mpi_scatterv_ (Fortran)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("mpi_scatterv", &ignorerc);
        FTN_X(__real_mpi_scatterv)(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, root, comm, ierr);
        ESMCI::TraceEventRegionExit("mpi_scatterv", &ignorerc);
        insideMPIRegion = 0;
      }
      else {
        FTN_X(__real_mpi_scatterv)(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, root, comm, ierr);
      }
    }
  
    extern void FTN_X(__real_mpi_scan)(MPI_Fint *sendbuf, MPI_Fint *recvbuf, MPI_Fint *count, MPI_Fint *datatype, MPI_Fint *op, MPI_Fint *comm, MPI_Fint *ierr);    

    void FTN_X(__wrap_mpi_scan)(MPI_Fint *sendbuf, MPI_Fint *recvbuf, MPI_Fint *count, MPI_Fint *datatype, MPI_Fint *op, MPI_Fint *comm, MPI_Fint *ierr) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_mpi_scan_ (Fortran)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("mpi_scan", &ignorerc);
        FTN_X(__real_mpi_scan)(sendbuf, recvbuf, count, datatype, op, comm, ierr);
        ESMCI::TraceEventRegionExit("mpi_scan", &ignorerc);
        insideMPIRegion = 0;
      }
      else {
        FTN_X(__real_mpi_scan)(sendbuf, recvbuf, count, datatype, op, comm, ierr);
      }
    }
  
    extern void FTN_X(__real_mpi_send)(MPI_Fint *buf, MPI_Fint *count, MPI_Fint *datatype, MPI_Fint *dest, MPI_Fint *tag, MPI_Fint *comm, MPI_Fint *ierr);    

    void FTN_X(__wrap_mpi_send)(MPI_Fint *buf, MPI_Fint *count, MPI_Fint *datatype, MPI_Fint *dest, MPI_Fint *tag, MPI_Fint *comm, MPI_Fint *ierr) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_mpi_send_ (Fortran)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("mpi_send", &ignorerc);
        FTN_X(__real_mpi_send)(buf, count, datatype, dest, tag, comm, ierr);
        ESMCI::TraceEventRegionExit("mpi_send", &ignorerc);
        insideMPIRegion = 0;
      }
      else {
        FTN_X(__real_mpi_send)(buf, count, datatype, dest, tag, comm, ierr);
      }
    }
  
    extern void FTN_X(__real_mpi_wait)(MPI_Fint *request, MPI_Fint *status, MPI_Fint *ierr);    

    void FTN_X(__wrap_mpi_wait)(MPI_Fint *request, MPI_Fint *status, MPI_Fint *ierr) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_mpi_wait_ (Fortran)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("mpi_wait", &ignorerc);
        FTN_X(__real_mpi_wait)(request, status, ierr);
        ESMCI::TraceEventRegionExit("mpi_wait", &ignorerc);
        insideMPIRegion = 0;
      }
      else {
        FTN_X(__real_mpi_wait)(request, status, ierr);
      }
    }
  
    extern void FTN_X(__real_mpi_waitall)(MPI_Fint *count, MPI_Fint *reqs, MPI_Fint *stats, MPI_Fint *ierr);    

    void FTN_X(__wrap_mpi_waitall)(MPI_Fint *count, MPI_Fint *reqs, MPI_Fint *stats, MPI_Fint *ierr) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_mpi_waitall_ (Fortran)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("mpi_waitall", &ignorerc);
        FTN_X(__real_mpi_waitall)(count, reqs, stats, ierr);
        ESMCI::TraceEventRegionExit("mpi_waitall", &ignorerc);
        insideMPIRegion = 0;
      }
      else {
        FTN_X(__real_mpi_waitall)(count, reqs, stats, ierr);
      }
    }
  
    extern void FTN_X(__real_mpi_waitany)(MPI_Fint *count, MPI_Fint *reqs, MPI_Fint *index, MPI_Fint *status, MPI_Fint *ierr);    

    void FTN_X(__wrap_mpi_waitany)(MPI_Fint *count, MPI_Fint *reqs, MPI_Fint *index, MPI_Fint *status, MPI_Fint *ierr) {
      if (c_esmftrace_isinitialized() == 1 && insideMPIRegion == 0) {
        //printf("__wrap_mpi_waitany_ (Fortran)\n");
        insideMPIRegion = 1;
        ESMCI::TraceEventRegionEnter("mpi_waitany", &ignorerc);
        FTN_X(__real_mpi_waitany)(count, reqs, index, status, ierr);
        ESMCI::TraceEventRegionExit("mpi_waitany", &ignorerc);
        insideMPIRegion = 0;
      }
      else {
        FTN_X(__real_mpi_waitany)(count, reqs, index, status, ierr);
      }
    }
  

}
#endif