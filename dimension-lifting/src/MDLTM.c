#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#include "dl_common.h"

static const int s0tiles=NTILES;
static const int s1tiles=NTILES;
static const int s2tiles=NTILES;


static inline int gamma(int i, int j, int k)
{
  return  i*s1*s2 + j*s2 + k;
}


// Here we go from the shape <s0,s1,s2> to the shape
// <s0tiles,s0/s0tiles,s1tiles,s1/s0tiles,s2tiles,s0/s2tiles>
// however we do not need to alter gamma

static void MC_DL_TILED(double *u,
                  const double *restrict v,
                  const double *restrict u0,
                  const double *restrict u1,
                  const double *restrict u2,
                  const double c0, const double c1, const double c2,
                  const double c3, const double c4)
{

#pragma omp parallel for collapse(3) schedule(static) num_threads(THREADS)
  for (int ti=0; ti<s0; ti+=s0/s0tiles) {
    for (int tj=0; tj<s1; tj+=s1/s1tiles) {
      for (int tk=0; tk<s2; tk+=s2/s2tiles) {
        for (int i=ti; i<ti+s0/s0tiles; i++) {
          for (int j=tj; j<tj+s1/s1tiles; j++) {
            for (int k=tk; k<tk+s2/s2tiles; k++) {
              u[gamma(i,j,k)] =
                u[gamma(i,j,k)] + c4 * (c3 * (c1 *
                v[gamma((mod(i-1,s0)),j,k)] +
                v[gamma((mod(i+1,s0)),j,k)] +
                v[gamma(i,(mod(j-1,s1)),k)] +
                v[gamma(i,(mod(j+1,s1)),k)] +
                v[gamma(i,j,(mod(k-1,s2)))] +
                v[gamma(i,j,(mod(k+1,s2)))]) -
                3 * c2 * u[gamma(i,j,k)] - c0 *
                ((v[gamma((mod(i+1,s0)),j,k)] -
                v[gamma((mod(i-1,s0)),j,k)]) *
                u0[gamma(i,j,k)] +
                (v[gamma(i,(mod(j+1,s1)),k)] -
                v[gamma(i,(mod(j-1,s1)),k)]) *
                u1[gamma(i,j,k)] +
                (v[gamma(i,j,(mod(k+1,s2)))] -
                v[gamma(i,j,(mod(k-1,s2)))]) *
                u2[gamma(i,j,k)]));
            }
          }
        }
      }
    }
  }
}

int main() {
  int status = EXIT_SUCCESS;
  if(s0 % s0tiles) {
    printf("s0tiles does not divide s0\n");
    return EXIT_FAILURE;
  }
  if(s1 % s1tiles) {
    printf("s1tiles does not divide s1\n");
    return EXIT_FAILURE;
  }
  if(s2 % s2tiles) {
    printf("s2tiles does not divide s2\n");
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
  printf("Multicore-dimension-lifted-tiled:");
  begin = omp_get_wtime();
  for(int i=0;i<steps;i++)
    step(array_size,u,v,s_nu,s_dx,s_dt,MC_DL_TILED);
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
