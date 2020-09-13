#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#include "dl_common.h"

//For dimension lifting on last dimension
static const int s2split = KSPLIT;

// We now split the index k in <i,j,k> to get <i,j,c,k>
// We go from shape <s0,s1,s2> to <s0,s1,s2split,s2/s2split>
// And parallelize on the new dimension. The same issue
// arises as before with mod having to propagate in

// i*s1*s2split*s2/s2split + j*s2split*(s2/s2split) + c*(s2/s2split)+ k
// i*s1*s2 + j*s2 + c*(s2/s2split) + k
static inline int gammaDLK(int i, int j, int c, int k)
{
  return  i*s1*s2 + j*s2 + c*(s2/s2split) + k;
}

static inline int gammaDLKMN(int i, int j, int c, int k)
{
  return  i*s1*s2 + j*s2 + mod(c*s2/s2split + k-1,s0);
}

static inline int gammaDLKMP(int i, int j, int c, int k)
{
  return  i*s1*s2 + j*s2 + mod(c*s2/s2split + k+1,s0);
}


static void MC_DL_ON_DIMK(double *u,
                        const double *restrict v,
                        const double *restrict u0,
                        const double *restrict u1,
                        const double *restrict u2,
                        const double c0, const double c1,
                        const double c2, const double c3,
                        const double c4)
{
// assuming cache line of 64 bytes, try to spread by at least that
#pragma omp parallel for collapse(2) schedule(static,8) num_threads(THREADS)
  for (int i=0; i<s0; i++) {
    for (int j=0; j<s1; j++) {
      for (int c=0; c<s2split; c++) {
        for (int k=0; k<s2/s2split; k++) {
          u[gammaDLK(i,j,c,k)] =      
            u[gammaDLK(i,j,c,k)] + c4 * (c3 * (c1 *
            v[gammaDLK(mod(i-1,s0),j,c,k)] +
            v[gammaDLK(mod(i+1,s0),j,c,k)] +
            v[gammaDLK(i,mod(j-1,s1),c,k)] +
            v[gammaDLK(i,mod(j+1,s1),c,k)] +
            v[gammaDLKMN(i,j,c,k)] +
            v[gammaDLKMP(i,j,c,k)]) -
            3 * c2 * u[gammaDLK(i,j,c,k)] - c0 *
            ((v[gammaDLK((mod(i+1,s0)),j,c,k)] -
            v[gammaDLK((mod(i-1,s0)),j,c,k)]) *
             u0[gammaDLK(i,j,c,k)] +
            (v[gammaDLK(i,mod(j+1,s1),c,k)] -
            v[gammaDLK(i,mod(j-1,s1),c,k)]) *
             u1[gammaDLK(i,j,c,k)] +
            (v[gammaDLKMP(i,j,c,k)] -
             v[gammaDLKMN(i,j,c,k)]) *
            u2[gammaDLK(i,j,c,k)]));
        }                
      }
    }
  }
}

int main() {
  int status = EXIT_SUCCESS;
  if(s2 % s2split) {
    printf("s2split does not divide s2\n");
    return EXIT_FAILURE;
  }

  //original data kept here
  double *start = malloc(total * sizeof(double));
  double *u = malloc(total * sizeof(double));
  double *v = malloc(total * sizeof(double));
  
  
  dumpsine(total,start);

  double begin;
  double end;
  double tspent;
  
  
  memcpy(u,start,total*sizeof(double));
  printf("Multicore-dimension-lifted-on-second-last-dimension:");
  begin = omp_get_wtime();
  for(int i=0;i<steps;i++)
    step(array_size,u,v,s_nu,s_dx,s_dt,MC_DL_ON_DIMK);
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
