#ifndef __COMMON_H
#define __COMMON_H

// PAD is hardcoded to 2. These experiments do not take advantage of padding
// exhaustion, and do not run in multicore settings. Therefore, changing this
// value is meaningless at the moment.
#define PAD 2

static const double s_dt = 0.00082212448155679772495;
static const double s_nu = 1.0;
static const double s_dx = 1.0;
static const int s0 = SIZE;
static const int s1 = SIZE;
static const int s2 = SIZE;
static const int array_size = s0 * s1 * s2;

static const int steps = 50;

// Test data generation
void dumpsine(int, double*);

// Row major indexing
static inline int _gamma(int i, int j, int k) {
  return  i * s1 * s2 + j * s2 + k;
}

static inline int mod(int inp, int modulus) {
    int remainder = inp % modulus;
    return remainder < 0 ? remainder + modulus : remainder;
}

typedef void (*kernel)(
    double*, const double*, const double*, const double*, const double*,
    const double, const double, const double, const double, const double);

void step(const int asize, double *u, double *v, double nu, double dx,
          double dt, kernel snippet);

#endif
