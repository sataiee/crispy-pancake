#include "mp.h"

/* This function calculates a global, axisymmetric field ("axifield")
   from current field "gridfield" */
void mpi_make1Dprofile (gridfield, axifield)
real* gridfield;
real* axifield;
{
  MPI_Request req1;
  int i, j, l;
  real *localaxifield;
  localaxifield = (real*) malloc(sizeof(real) * NRAD);
  if ( localaxifield == NULL ) 
    erreur ("Not enough memory in axisglib.c ; suspect...");
  /* We first calculate axisymmetric local field */
  for ( i = 0; i < NRAD; i++ )
    localaxifield[i] = 0.;
  for (i = Zero_or_active; i < Max_or_active; i++) {
    for( j = 0; j < NSEC; j++ ) {
      l = i*NSEC + j;
      localaxifield[i] += gridfield[l];
    }
    localaxifield[i] /= (real)NSEC;
  }
  /* Then we share it with other cpus to yield a global, axisymmetric
  field */
  if ( CPU_Number == 1 ) {
    for ( i = 0; i < GLOBALNRAD; i++ )
      axifield[i] = localaxifield[i];
  }
  if ( CPU_Number > 1 ) {
    if ( CPU_Rank == 0 ) {
      for ( i = 0; i < GLOBALNRAD; i++ ) {
        if ( i < Max_or_active )
          axifield[i] = localaxifield[i];
        else
          axifield[i] = 0.;
      }
      MPI_Isend (axifield, GLOBALNRAD, MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req1);
      MPI_Wait (&req1, &fargostat);
    }
    if ( CPU_Rank != 0 ) {
      MPI_Irecv (axifield, GLOBALNRAD, MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &req1);
      MPI_Wait (&req1, &fargostat);
      for (i = Zero_or_active; i < Max_or_active; i++)
        axifield[i+IMIN] = localaxifield[i];
      if ( CPU_Rank != CPU_Highest ) {
        MPI_Isend (axifield, GLOBALNRAD, MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req1);
        MPI_Wait (&req1, &fargostat);
      }
    }
    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Bcast (axifield, GLOBALNRAD, MPI_DOUBLE, CPU_Highest, MPI_COMM_WORLD);
  }
  free (localaxifield);
  }

void mpi_Find1Dprofile (gridfield, ithetp, axifield)
real* gridfield;
real* axifield;
int ithetp;
{
  MPI_Request req1;
  int i, j, l;
  real *localaxifield;
  localaxifield = (real*) malloc(sizeof(real) * NRAD);
  if ( localaxifield == NULL ) 
    erreur ("Not enough memory in axisglib.c ; suspect...");
  /* write the field at the given azimuth index */
  for ( i = 0; i < NRAD; i++ )
    localaxifield[i] = 0.;
  for (i = Zero_or_active; i < Max_or_active; i++) {
    l = i*NSEC + ithetp;
    localaxifield[i] = gridfield[l];
  }
  /* Then we share it with other cpus to yield a global field */
  if ( CPU_Number == 1 ) {
    for ( i = 0; i < GLOBALNRAD; i++ )
      axifield[i] = localaxifield[i];
  }
  if ( CPU_Number > 1 ) {
    if ( CPU_Rank == 0 ) {
      for ( i = 0; i < GLOBALNRAD; i++ ) {
        if ( i < Max_or_active )
          axifield[i] = localaxifield[i];
        else
          axifield[i] = 0.;
      }
      MPI_Isend (axifield, GLOBALNRAD, MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req1);
      MPI_Wait (&req1, &fargostat);
    }
    if ( CPU_Rank != 0 ) {
      MPI_Irecv (axifield, GLOBALNRAD, MPI_DOUBLE, CPU_Prev, 0, MPI_COMM_WORLD, &req1);
      MPI_Wait (&req1, &fargostat);
      for (i = Zero_or_active; i < Max_or_active; i++)
        axifield[i+IMIN] = localaxifield[i];
      if ( CPU_Rank != CPU_Highest ) {
        MPI_Isend (axifield, GLOBALNRAD, MPI_DOUBLE, CPU_Next, 0, MPI_COMM_WORLD, &req1);
        MPI_Wait (&req1, &fargostat);
      }
    }
    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Bcast (axifield, GLOBALNRAD, MPI_DOUBLE, CPU_Highest, MPI_COMM_WORLD);
  }
  free (localaxifield);
}
