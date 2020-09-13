#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#include "dl_common.h"


// Introduce a thread dimension; we go from shape <s0,s1,s2> to
// <THREADS,s0/THREADS,s1,s2>. Before we had wrap around effects
// at i==0 and i==s0-1. These will now occur elsewhere.
// If the indices were <i,j,k> and are now <p,i,j,k> we must do the mapping
// f(<p,i,j,k>)= gamma(i*s0/THREADS,j,k). Call this map f gammaDL.
// mod calculations that affect i must propagate into gammaDL
// generating two new functions gammaDLMN, gammaDLNMP

static inline int gammaDL(int p, int i, int j, int k)
{
  return  (p*s0/THREADS+i)*s1*s2 + j*s2 + k;
}

// here MN means mod negative
static inline int gammaDLMN(int p, int i, int j, int k)
{
  return  mod(p*s0/THREADS+i-1,s0)*s1*s2 + j*s2 + k;
}

// here MP means mod positive
static inline int gammaDLMP(int p, int i, int j, int k)
{
  return  mod(p*s0/THREADS+i+1,s0)*s1*s2 + j*s2 + k;
}

static void MC_DL_ON_THREADS(double *u,
                        const double *restrict v,
                        const double *restrict u0,
                        const double *restrict u1,
                        const double *restrict u2,
                        const double c0, const double c1,
                        const double c2, const double c3,
                        const double c4)
{

#pragma omp parallel for schedule(static) num_threads(THREADS)
  for (int p=0;p<THREADS;p++) {
    for (int i=0; i<s0/THREADS; i++) {
      for (int j=0; j<s1; j++) {
        for (int k=0; k<s2; k++) {
          u[gammaDL(p,i,j,k)]=
          u[gammaDL(p,i,j,k)] + c4 * (c3 * (c1 *
          v[gammaDLMN(p,i,j,k)] +
          v[gammaDLMP(p,i,j,k)] +
          v[gammaDL(p,i,mod(j-1,s1),k)] +
          v[gammaDL(p,i,mod(j+1,s1),k)] +
          v[gammaDL(p,i,j,mod(k-1,s2))] +
          v[gammaDL(p,i,j,mod(k+1,s2))]) -
          3 * c2 * u[gammaDL(p,i,j,k)] - c0 *
          ((v[gammaDLMP(p,i,j,k)] -
          v[gammaDLMN(p,i,j,k)]) *
          u0[gammaDL(p,i,j,k)] +
          (v[gammaDL(p,i,mod(j+1,s1),k)] -
          v[gammaDL(p,i,mod(j-1,s1),k)]) *
          u1[gammaDL(p,i,j,k)] +
          (v[gammaDL(p,i,j,mod(k+1,s2))] -
          v[gammaDL(p,i,j,mod(k-1,s2))]) *
          u2[gammaDL(p,i,j,k)]));
        }
      }
    }
  }
}

int main() {
  int status = EXIT_SUCCESS;

  //original data kept here
  double *start = malloc(total * sizeof(double));
  double *u = malloc(total * sizeof(double));
  double *v = malloc(total * sizeof(double));
  
  
  dumpsine(total,start);

  double begin;
  double end;
  double tspent;

  memcpy(u,start,total*sizeof(double));
  printf("Multicore-dimension-lifted-on-threads:");
  begin = omp_get_wtime();
  for(int i=0;i<steps;i++) {
    step(array_size,u,v,s_nu,s_dx,s_dt,MC_DL_ON_THREADS);
  }
  end = omp_get_wtime();
  tspent = end - begin;
  printf("%lf\n",tspent);
  for (int i =0;i<total;i++) {
    if (isnan(u[i])) {
      printf("NAN detected\n");
      status = EXIT_FAILURE;
      goto cleanup;

    }
  }
  
 cleanup:
  free(start);
  free(u);
  free(v);
  
  return status;

}
