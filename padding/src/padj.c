#include <assert.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "common.h"

static const int s1padded = SIZE+PAD;
static const int array_size_padded = s0 * s1padded * s2;

static inline int gamma_padded(int i, int j, int k) {
  return  i * s1padded * s2 + j * s2 + k;
}

static void SC_PADDING_J (double *u,
                          const double *restrict v,
                          const double *restrict u0,
                          const double *restrict u1,
                          const double *restrict u2,
                          const double c0, const double c1, const double c2,
                          const double c3, const double c4) {
  for (int i=0; i<s0; i++) {
    for (int j=1; j<s1+1; j++) {
      for (int k=0; k<s2; k++) {
        u[gamma_padded(i,j,k)] =
        u[gamma_padded(i,j,k)] + c4 * (c3 * (c1 *
        v[gamma_padded((mod(i-1,s0)),j,k)] +
        v[gamma_padded((mod(i+1,s0)),j,k)] +
        v[gamma_padded(i,j-1,k)] +
        v[gamma_padded(i,j+1,k)] +
        v[gamma_padded(i,j,k-1)] +
        v[gamma_padded(i,j,k+1)]) -
        3 * c2 * u[gamma_padded(i,j,k)] - c0 *
        ((v[gamma_padded((mod(i+1,s0)),j,k)] -
        v[gamma_padded((mod(i-1,s0)),j,k)]) *
        u0[gamma_padded(i,j,k)] +
        (v[gamma_padded(i,j+1,k)] -
        v[gamma_padded(i,j-1,k)]) *
        u1[gamma_padded(i,j,k)] +
        (v[gamma_padded(i,j,k+1)] -
        v[gamma_padded(i,j,k-1)]) *
         u2[gamma_padded(i,j,k)]));
      }
    }
    // update the padding values
    memcpy(u+gamma_padded(i,0,0),u+gamma_padded(i,s1,0),sizeof(double)*s2);
    memcpy(u+gamma_padded(i,s1+1,0),u+gamma_padded(i,0,0),sizeof(double)*s2);

  }
}

void padj() {
  const int total = array_size_padded * 3;

  //original data kept here
  double *start = malloc(total * sizeof(double));
  dumpsine(total,start);

  clock_t begin;
  clock_t end;
  double tspent;

  //padded version of u, padded along all dimensions
  double *upadded = malloc(total * sizeof(double));
  double *vpadded = malloc(total * sizeof(double));

  // we must do these replications
  // The padded values
  // gamma2(i,0,k)  == gamma(i,s1-1,k)
  // gamma2(i,s1,k) == gamma(i,0,k)
  // The regular values
  // gamma2(i,j,k) == gamma(i,j,k)

  for (int part = 0;part<3;part++) {
    int offsetpadded = part*array_size_padded;
    int offset = part*array_size;
    for(int i = 0; i<s0; i++) {
      for(int j = -1; j<s1+1; j++) {
        for(int k = 0; k<s2; k++) {
          int jj  = (j==-1)? s1-1 : (j==s1)? 0 : j;
          upadded[offsetpadded+gamma_padded(i,j+1,k)] = start[offset+_gamma(i,jj,k)];
        }
      }
    }
  }

  printf("Singlecore-with-padding-second-axis:");
  begin = clock();
  for(int i=0;i<steps;i++)
    step(array_size_padded,upadded,vpadded,s_nu,s_dx,s_dt,SC_PADDING_J);
  end = clock();
  tspent = ((double)(end - begin))/CLOCKS_PER_SEC;
  printf("%f\n",tspent);
  for (int i =0;i<total;i++) {
    if (isnan(upadded[i])) {
      printf("NAN detected\n");
    }
  }
  // Now there would need to be a reversal of the padding if the expected
  // output is the data as an array of the kind start is, but that is omitted
  // and we consider the computation complete

  free(start);
  free(upadded);
  free(vpadded);
}

int main() {
  padj();
  return EXIT_SUCCESS;
}

