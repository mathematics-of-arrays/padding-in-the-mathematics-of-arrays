#include <math.h>
#include <string.h>

#include "common.h"

// Test data generation
void dumpsine(int n, double* a) {
  double step = 0.01;
  double PI = 3.14159265358979323846;
  double amplitude = 10.0;
  double phase = 0.0125;
  double t = 0.0;
  for (int i = 0; i < n; i++) {
    a[i] = amplitude * sin(PI * t + phase);
    t += step;
  }
}

void step(const int asize, double *u, double *v,
          double nu, double dx, double dt, kernel snippet) {
  double c0 = 0.5/dx;
  double c1 = 1/dx/dx;
  double c2 = 2/dx/dx;
  double c3 = nu;
  double c4 = dt/2;
 
  memcpy(v, u, asize * 3 * sizeof(double));
  double *u0 = &u[0];
  double *u1 = &u[asize];
  double *u2 = &u[2 * asize];

  double *v0 = &v[0];
  double *v1 = &v[asize];
  double *v2 = &v[2 * asize];

  snippet(v0,u0,u0,u1,u2,c0,c1,c2,c3,c4);
  snippet(v1,u1,u0,u1,u2,c0,c1,c2,c3,c4);
  snippet(v2,u2,u0,u1,u2,c0,c1,c2,c3,c4);
  snippet(u0,v0,v0,v1,v2,c0,c1,c2,c3,c4);
  snippet(u1,v1,v0,v1,v2,c0,c1,c2,c3,c4);
  snippet(u2,v2,v0,v1,v2,c0,c1,c2,c3,c4);
}
