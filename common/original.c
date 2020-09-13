#include <assert.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "common.h"

static void SC_NO_PADDING(double *u,
                  const double *restrict v,
                  const double *restrict u0,
                  const double *restrict u1,
                  const double *restrict u2,
                  const double c0, const double c1, const double c2,
                  const double c3, const double c4) {
  for (int i=0; i<s0; i++) {
    for (int j=0; j<s1; j++) {
      for (int k=0; k<s2; k++) {
        u[_gamma(i,j,k)] =
        u[_gamma(i,j,k)] + c4 * (c3 * (c1 *
        v[_gamma((mod(i-1,s0)),j,k)] +
        v[_gamma((mod(i+1,s0)),j,k)] +
        v[_gamma(i,(mod(j-1,s1)),k)] +
        v[_gamma(i,(mod(j+1,s1)),k)] +
        v[_gamma(i,j,(mod(k-1,s2)))] +
        v[_gamma(i,j,(mod(k+1,s2)))]) -
        3 * c2 * u[_gamma(i,j,k)] - c0 *
        ((v[_gamma((mod(i+1,s0)),j,k)] -
        v[_gamma((mod(i-1,s0)),j,k)]) *
        u0[_gamma(i,j,k)] +
        (v[_gamma(i,(mod(j+1,s1)),k)] -
        v[_gamma(i,(mod(j-1,s1)),k)]) *
        u1[_gamma(i,j,k)] +
        (v[_gamma(i,j,(mod(k+1,s2)))] -
        v[_gamma(i,j,(mod(k-1,s2)))]) *
        u2[_gamma(i,j,k)]));
      }
    }
  }
}

void nopad() {
  const int total = array_size*3;
  //original data kept here
  double *start = malloc(total * sizeof(double));
  double *u = malloc(total * sizeof(double));
  double *v = malloc(total * sizeof(double));

  dumpsine(total,start);

  clock_t begin;
  clock_t end;
  double tspent;

  memcpy(u,start,total*sizeof(double));
  printf("Singlecore-no-padding:");
  begin = clock();
  for(int i=0;i<steps;i++)
    step(array_size,u,v,s_nu,s_dx,s_dt,SC_NO_PADDING);
  end = clock();
  tspent = ((double)(end - begin))/CLOCKS_PER_SEC;
  printf("%f\n",tspent);
  for (int i =0;i<total;i++) {
    if (isnan(u[i])) {
      printf("NAN detected\n");
    }
  }

  free(u);
  free(v);
  free(start);

}

int main() {
  nopad();
  return EXIT_SUCCESS;
}

