#include <math.h>
#include <stdint.h>
#include <stddef.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "Maxfiles.h"
#include <MaxSLiCInterface.h>

#define MAXIT 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

float betacf(float a, float b, float x) {
  float qab = a + b; //These q’s will be used in factors that occur in the coefficients
  float qap = a + 1.0;
  float qam = a - 1.0;

  float c = 1.0; // First step of Lentz’s method.

  float d = 1.0 - qab * x / qap;
  if (fabs(d) < FPMIN) d = FPMIN;
  d = 1.0 / d;

  float h = d;
  int m;
  for (m = 1; m <= MAXIT; m++) {
    int m2 = 2 * m;

    float tmp1 = m * (b - m) * x / ((qam + m2) * (a + m2));
    d = 1.0 + tmp1 * d; // One step (the even one) of the recurrence.
    if (fabs(d) < FPMIN) d = FPMIN;
    d = 1.0 / d;

    c = 1.0 + tmp1 / c;
    if (fabs(c) < FPMIN) c = FPMIN;

    h *= d * c;

    float tmp2 = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
    d = 1.0 + tmp2 * d;  // Next step of the recurrence (the odd one).
    if (fabs(d) < FPMIN) d = FPMIN;
    d = 1.0 / d;

    c = 1.0 + tmp2 / c;
    if (fabs(c) < FPMIN) c = FPMIN;

    float del = d * c;

    h *= del;
    if (fabs(del - 1.0) < EPS) break;
  }

  assert (m <= MAXIT); //a or b too big, or MAXIT too small in betacf

  return h;
}

// Returns the value ln[Γ(xx)] for xx > 0.
float gammln(float xx) {
  // Internal arithmetic will be done in double precision, a nicety that you can omit if five-figure accuracy is good enough.

  static double cof[6] = {76.18009172947146, -86.50532032941677,
                          24.01409824083091, -1.231739572450155,
                          0.1208650973866179e-2, -0.5395239384953e-5};

  double x, y;
  y = x = xx;

  double tmp = x + 5.5;
  tmp -= (x + 0.5) * log(tmp);

  double ser = 1.000000000190015;
  for (int j = 0; j <= 5; j++) ser += cof[j] / ++y;

  return -tmp + log(2.5066282746310005 * ser / x);
}

// Returns the incomplete beta function Ix(a, b).
float betai(float a, float b, float x) {
  assert (!(x < 0.0 || x > 1.0));

  float bt;

  if (x == 0.0 || x == 1.0) bt = 0.0;
  else bt = exp(gammln(a + b) - gammln(a) - gammln(b) + a * log(x) + b * log(1.0 - x));

  if (x < (a + 1.0) / (a + b + 2.0))
    return bt * betacf(a, b, x) / a;
  else
    return 1.0 - bt * betacf(b, a, 1.0 - x) / b;
}

void avevar(float data[], size_t n, float *average, float *variance) {
  *average = 0.0;
  for (size_t i = 0; i < n; i++) *average += data[i];
  *average /= n;

  float deviation_sum;
  *variance = deviation_sum = 0.0;
  for (size_t i = 0; i < n; i++) {
    float deviation = data[i] - (*average);
    deviation_sum += deviation;
    *variance += deviation * deviation;
  }

  //corrected two-pass algorithm
  *variance = (*variance - deviation_sum * deviation_sum / n) / (n - 1);
}

void ftest(float data1[], size_t n1, float data2[], size_t n2, float *f, float *prob) {
  float var1, var2, ave1, ave2, df1, df2;
  avevar(data1, n1, &ave1, &var1);
  avevar(data2, n2, &ave2, &var2);

  if (var1 > var2) {
    *f = var1 / var2;
    df1 = n1 - 1;
    df2 = n2 - 1;
  } else {
    *f = var2 / var1;
    df1 = n2 - 1;
    df2 = n1 - 1;
  }
  *prob = 2.0 * betai(0.5 * df2, 0.5 * df1, df2 / (df2 + df1 * (*f)));
  if (*prob > 1.0) *prob = 2.0 - *prob;
}

int main() {
  int n = 32;
  float data1[n], data2[n];

  for (int i = 0; i < n; ++i) {
    data1[i] = rand() / ((float)RAND_MAX + 1);
    data2[i] = rand() / ((float)RAND_MAX + 1);

    printf("%f %f\n", data1[i], data2[i]);
  }

  float f_expected, prob_expected;
  ftest(data1, n, data2, n, &f_expected, &prob_expected);

  float f_result[n], prob_result[n];
  FTest(n, data1, data2, f_result, prob_result);

  printf("f=%f %f     prob=%f %f\n", f_expected, f_result[n-1], prob_expected, prob_result[n-1]);

  return 0;
}
