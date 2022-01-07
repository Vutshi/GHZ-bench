#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#if defined(__GNUC__) || defined(__clang__)
#  define INLINE_NEVER __attribute__((noinline))
#else
#  define INLINE_NEVER
#endif

// These parameters are processor dependent
#if defined(__e2k__)
#  if __iset__ > 4
#    define lanes 8
#  else
#    define lanes 4
#  endif
#  define step1 6
#  define step2 6
#  define step3 3
#else
// on x86 clang and gcc do not vectorise properly with lanes > 1
#  define lanes 1
#  define step1 6
#  define step2 6 // Skylake likes 6
#  define step3 3
#endif
#define num_it (2 * 3 * 5 * 6 * 7 * 8)

#define TEST_NAIVE 1
#define TEST_UNROLL 1
#define TEST_UNROLL_SOA 1
#define TEST_UNROLL_SOA_SPLIT 1

#if !defined(SANITY_CHECK_FLOPS)
#  define SANITY_CHECK_FLOPS 0
#endif

// Parameters of the underlying problem. Should not be changed
#define N 5
#define M 11
#define dim ((M + 1) * (M + 1))
#define index_max (2 * M * ((1 << N) - 1))
#define dim_L (index_max + 1)

#define index(I) (M * (1 << (I)))

#define len(I) (2 * index(I) + 1)
#define len_0 len(0)
#define len_1 (len_0 + len(1))
#define len_2 (len_1 + len(2))
#define len_3 (len_2 + len(3))
#define len_4 (len_3 + len(4))

// *********************************
// This benchmark features the main loop of a program performing quantum
// mechanical calculations. It evaluates noise properties of a quantum
// interferometer employing a cascade of GHZ states
// *********************************

#if TEST_NAIVE
// *********************************
// #1 Straightforward approach
// *********************************
INLINE_NEVER
double test_naive(double *Lexp, double complex *temp_ij) {
  double out = 0.0;
#if SANITY_CHECK_FLOPS
  long flop = 0;
#endif
  for (long i0 = 0; i0 < num_it / lanes; i0++) {
    for (long q = 0; q < lanes; q++) {
      double p_i = 0.0;
      double est_i = 0.0;

      for (long l = 0; l < index_max; l++) {
        double complex tmp = 0.0;

        for (long ll = 0; ll <= l; ll++) {
          tmp += conj(temp_ij[q * len_4 + index_max - (l - ll)]) *
                 temp_ij[q * len_4 + ll];
#if SANITY_CHECK_FLOPS
          flop += 8;
#endif
        }
        tmp *= Lexp[l];
        p_i += creal(tmp);
        est_i -= cimag(tmp) * (l - index_max);
#if SANITY_CHECK_FLOPS
        flop += 4;
#endif
      }
      out += est_i * (est_i / p_i);
#if SANITY_CHECK_FLOPS
      flop += 3;
#endif
    }
  }
#if SANITY_CHECK_FLOPS
  printf("sainity check, FLOP = %ld\n", flop);
#endif
  return out;
}
#endif

#if TEST_UNROLL
// *********************************
// #2 Unrolled method
// Here the l loop is unrolled by step1 to reduce memory loads
// *********************************
INLINE_NEVER
double test_unroll(double *Lexp, double complex *temp_ij) {
  double out = 0.0;

  for (long i0 = 0; i0 < num_it / lanes; i0++) {
    long l = 0;
    double p_i[lanes] = {0};
    double est_i[lanes] = {0};

    // unrolled main loop
    for (; l < index_max - index_max % step1; l += step1) {
      double complex temp[step1 * lanes] = {0};

#pragma comb_oper
      for (long ll = 0; ll <= l; ll++) {
        for (long j = 0; j < step1; j++) {
#if defined(__e2k__) &&  __iset__ > 4
#  pragma unroll(lanes)
#  pragma vector always
#endif
          for (long q = 0; q < lanes; q++) {
            temp[j * lanes + q] +=
                conj(temp_ij[(index_max - (l + j - ll)) * lanes + q]) *
                temp_ij[ll * lanes + q];
          }
        }
      }
      // epilogue for tail iterations
#pragma comb_oper
      for (long k = 1; k < step1; k++) {
        for (long j = k; j < step1; j++) {
#if defined(__e2k__) && __iset__ > 4
#  pragma unroll(lanes)
#  pragma vector always
#endif
          for (long q = 0; q < lanes; q++) {
            temp[j * lanes + q] +=
                conj(temp_ij[(index_max - (j - k)) * lanes + q]) *
                temp_ij[(l + k) * lanes + q];
          }
        }
      }
      for (long j = 0; j < step1; j++) {
        for (long q = 0; q < lanes; q++) {
          temp[j * lanes + q] *= Lexp[l + j];
          p_i[q] += creal(temp[j * lanes + q]);
          est_i[q] -= cimag(temp[j * lanes + q]) * (l + j - index_max);
        }
      }
    }
    // tail main loop
    for (long l = index_max - index_max % step1; l < index_max; l++) {
      double complex tmp[lanes] = {0};

      for (long ll = 0; ll <= l; ll++) {
        for (long q = 0; q < lanes; q++) {
          tmp[q] += conj(temp_ij[(index_max - (l - ll)) * lanes + q]) *
                    temp_ij[ll * lanes + q];
        }
      }
      for (long q = 0; q < lanes; q++) {
        tmp[q] *= Lexp[l];
        p_i[q] += creal(tmp[q]);
        est_i[q] -= cimag(tmp[q]) * (l - index_max);
      }
    }
    for (long q = 0; q < lanes; q++) {
      out += est_i[q] * (est_i[q] / p_i[q]);
    }
  }
  return out;
}
#endif

#if TEST_UNROLL_SOA
// *********************************
// #3 Unrolled method, SoA
// Here complex numbers are represented by separate arrays of real and
// imaginary parts
// *********************************
INLINE_NEVER
double test_unroll_SoA(double *Lexp, double *temp_ij_re, double *temp_ij_im) {
  double out = 0.0;

  for (long i0 = 0; i0 < num_it / lanes; i0++) {
    long l = 0;
    double p_i[lanes] = {0};
    double est_i[lanes] = {0};

    // unrolled main loop
    for (; l < index_max - index_max % step2; l += step2) {
      double tmp_re[step2 * lanes] = {0};
      double tmp_im[step2 * lanes] = {0};

#if defined(__e2k__)
#pragma comb_oper
#endif
      for (long ll = 0; ll <= l; ll++) {
#pragma unroll(step2)
        for (long j = 0; j < step2; j++) {
#if defined(__e2k__) && __iset__ > 4
#  pragma unroll(lanes)
#  pragma vector always
#endif
          for (long q = 0; q < lanes; q++) {
            tmp_re[j * lanes + q] +=
                temp_ij_re[(index_max - (l + j - ll)) * lanes + q] *
                temp_ij_re[ll * lanes + q];
            tmp_im[j * lanes + q] -=
                temp_ij_im[(index_max - (l + j - ll)) * lanes + q] *
                temp_ij_re[ll * lanes + q];

            tmp_re[j * lanes + q] +=
                temp_ij_im[(index_max - (l + j - ll)) * lanes + q] *
                temp_ij_im[ll * lanes + q];
            tmp_im[j * lanes + q] +=
                temp_ij_re[(index_max - (l + j - ll)) * lanes + q] *
                temp_ij_im[ll * lanes + q];
          }
        }
      }

      // epilogue for tail iterations
#if defined(__e2k__)
#  pragma comb_oper
#endif
      for (long k = 1; k < step2; k++) {
        for (long j = k; j < step2; j++) {
#if defined(__e2k__) && __iset__ > 4
#  pragma unroll(lanes)
#  pragma vector always
#endif
          for (long q = 0; q < lanes; q++) {
            tmp_re[j * lanes + q] +=
                temp_ij_re[(index_max - (j - k)) * lanes + q] *
                temp_ij_re[(l + k) * lanes + q];
            tmp_im[j * lanes + q] -=
                temp_ij_im[(index_max - (j - k)) * lanes + q] *
                temp_ij_re[(l + k) * lanes + q];
            tmp_re[j * lanes + q] +=
                temp_ij_im[(index_max - (j - k)) * lanes + q] *
                temp_ij_im[(l + k) * lanes + q];
            tmp_im[j * lanes + q] +=
                temp_ij_re[(index_max - (j - k)) * lanes + q] *
                temp_ij_im[(l + k) * lanes + q];
          }
        }
      }
      for (long j = 0; j < step2; j++) {
#pragma unroll(lanes)
        for (long q = 0; q < lanes; q++) {
          p_i[q] += tmp_re[j * lanes + q] * Lexp[l + j];
          est_i[q] -= tmp_im[j * lanes + q] * Lexp[l + j] * (l + j - index_max);
        }
      }
    }

    // tail main loop
    for (long l = index_max - index_max % step2; l < index_max; l++) {
      double temp_re[lanes] = {0};
      double temp_im[lanes] = {0};

      for (long ll = 0; ll <= l; ll++) {
#if defined(__e2k__) && __iset__ > 4
#  pragma unroll(lanes)
#  pragma vector always
#endif
        for (long q = 0; q < lanes; q++) {
          temp_re[q] += temp_ij_re[(index_max - (l - ll)) * lanes + q] *
                        temp_ij_re[ll * lanes + q];
          temp_im[q] -= temp_ij_im[(index_max - (l - ll)) * lanes + q] *
                        temp_ij_re[ll * lanes + q];
          temp_re[q] += temp_ij_im[(index_max - (l - ll)) * lanes + q] *
                        temp_ij_im[ll * lanes + q];
          temp_im[q] += temp_ij_re[(index_max - (l - ll)) * lanes + q] *
                        temp_ij_im[ll * lanes + q];
        }
      }
#pragma unroll(lanes)
      for (long q = 0; q < lanes; q++) {
        p_i[q] += temp_re[q] * Lexp[l];
        est_i[q] -= temp_im[q] * Lexp[l] * (l - index_max);
      }
    }
#pragma unroll(lanes)
    for (long q = 0; q < lanes; q++) {
      out += est_i[q] * (est_i[q] / p_i[q]);
    }
  }
  return out;
}
#endif

#if TEST_UNROLL_SOA_SPLIT
// *********************************
// #4 Unrolled method, SoA + split
// In addition to SoA representation the temporary re and im arrays are split
// into two parts to avoid consecutive multiply add
// Unroll parameters (lanes and step) are chosen such that 6 ALC are saturated
// given that there are at most 4 loads per VLIW
// *********************************
INLINE_NEVER
double test_unroll_SoA_split(double *Lexp, double *temp_ij_re,
                             double *temp_ij_im) {
  double out = 0.0;
  double tmp_re_a[step3 * lanes] = {0};
  double tmp_re_b[step3 * lanes] = {0};
  double tmp_im_a[step3 * lanes] = {0};
  double tmp_im_b[step3 * lanes] = {0};

  for (long i0 = 0; i0 < num_it / lanes; i0++) {
    double p_i[lanes] = {0};
    double est_i[lanes] = {0};

    // unrolled main loop
    for (long l = 0; l < index_max - index_max % step3; l += step3) {
      long ll = 0;

#if defined(__e2k__)
      // init tmp arrays
      for (long j = 0; j < step3; j++) {
#if defined(__e2k__) && __iset__ > 4
#  pragma unroll(lanes)
#  pragma vector always
#endif
        for (long q = 0; q < lanes; q++) {
          tmp_re_a[j * lanes + q] =
              temp_ij_re[(index_max - (l + j - ll)) * lanes + q] *
              temp_ij_re[ll * lanes + q];
          tmp_im_a[j * lanes + q] =
              -temp_ij_im[(index_max - (l + j - ll)) * lanes + q] *
              temp_ij_re[ll * lanes + q];
          tmp_re_b[j * lanes + q] =
              temp_ij_im[(index_max - (l + j - ll)) * lanes + q] *
              temp_ij_im[ll * lanes + q];
          tmp_im_b[j * lanes + q] =
              temp_ij_re[(index_max - (l + j - ll)) * lanes + q] *
              temp_ij_im[ll * lanes + q];
        }
      }
      ll += 1;
#endif

      for (; ll <= l; ll++) {
        for (long j = 0; j < step3; j++) {
#if defined(__e2k__) && __iset__ > 4
#  pragma unroll(lanes)
#  pragma vector always
#endif
          for (long q = 0; q < lanes; q++) {
            tmp_re_a[j * lanes + q] +=
                temp_ij_re[(index_max - (l + j - ll)) * lanes + q] *
                temp_ij_re[ll * lanes + q];
            tmp_im_a[j * lanes + q] -=
                temp_ij_im[(index_max - (l + j - ll)) * lanes + q] *
                temp_ij_re[ll * lanes + q];
            tmp_re_b[j * lanes + q] +=
                temp_ij_im[(index_max - (l + j - ll)) * lanes + q] *
                temp_ij_im[ll * lanes + q];
            tmp_im_b[j * lanes + q] +=
                temp_ij_re[(index_max - (l + j - ll)) * lanes + q] *
                temp_ij_im[ll * lanes + q];
          }
        }
      }

      // epilogue for tail iterations

#define EPILOGUE(k)                                                            \
  if (k < step3) {                                                             \
    for (long j = k; j < step3; j++) {                                         \
      for (long q = 0; q < lanes; q++) {                                       \
        tmp_re_a[j * lanes + q] +=                                             \
            temp_ij_re[(index_max - (j - k)) * lanes + q] *                    \
            temp_ij_re[(l + k) * lanes + q];                                   \
        tmp_im_a[j * lanes + q] -=                                             \
            temp_ij_im[(index_max - (j - k)) * lanes + q] *                    \
            temp_ij_re[(l + k) * lanes + q];                                   \
        tmp_re_b[j * lanes + q] +=                                             \
            temp_ij_im[(index_max - (j - k)) * lanes + q] *                    \
            temp_ij_im[(l + k) * lanes + q];                                   \
        tmp_im_b[j * lanes + q] +=                                             \
            temp_ij_re[(index_max - (j - k)) * lanes + q] *                    \
            temp_ij_im[(l + k) * lanes + q];                                   \
      }                                                                        \
    }                                                                          \
  }

#if step3 > 16
#error "step3 must be less than 16"
#endif

      EPILOGUE(1)
      EPILOGUE(2)
      EPILOGUE(3)
      EPILOGUE(4)
      EPILOGUE(5)
      EPILOGUE(6)
      EPILOGUE(7)
      EPILOGUE(8)
      EPILOGUE(9)
      EPILOGUE(10)
      EPILOGUE(11)
      EPILOGUE(12)
      EPILOGUE(13)
      EPILOGUE(14)
      EPILOGUE(15)
#undef EPILOGUE

      // #pragma unroll(step3)
      for (long j = 0; j < step3; j++) {
        for (long q = 0; q < lanes; q++) {
          p_i[q] +=
              (tmp_re_a[j * lanes + q] + tmp_re_b[j * lanes + q]) * Lexp[l + j];
          est_i[q] -= (tmp_im_a[j * lanes + q] + tmp_im_b[j * lanes + q]) *
                      Lexp[l + j] * (l + j - index_max);
        }
      }
    }
    // tail main loop
    for (long l = index_max - index_max % step3; l < index_max; l++) {
      double temp_re[lanes] = {0};
      double temp_im[lanes] = {0};
      for (long ll = 0; ll <= l; ll++) {
#if defined(__e2k__) && __iset__ > 4
#  pragma unroll(lanes)
#  pragma vector always
#endif
        for (long q = 0; q < lanes; q++) {
          temp_re[q] += temp_ij_re[(index_max - (l - ll)) * lanes + q] *
                        temp_ij_re[ll * lanes + q];
          temp_im[q] -= temp_ij_im[(index_max - (l - ll)) * lanes + q] *
                        temp_ij_re[ll * lanes + q];
          temp_re[q] += temp_ij_im[(index_max - (l - ll)) * lanes + q] *
                        temp_ij_im[ll * lanes + q];
          temp_im[q] += temp_ij_re[(index_max - (l - ll)) * lanes + q] *
                        temp_ij_im[ll * lanes + q];
        }
      }
      for (long q = 0; q < lanes; q++) {
        p_i[q] += temp_re[q] * Lexp[l];
        est_i[q] -= temp_im[q] * Lexp[l] * (l - index_max);
      }
    }
    for (long q = 0; q < lanes; q++) {
      out += est_i[q] * (est_i[q] / p_i[q]);
    }
  }
  return out;
}
#endif

static void print_result(double theory_FLOP, clock_t t1, clock_t t2, double res) {
  double t = t2 - t1;
  printf("time = %fs\n", t / CLOCKS_PER_SEC);
  printf("GFLOPS = %f\n", (theory_FLOP * CLOCKS_PER_SEC) / t * 1e-9);
  printf("%.17g\n", res);
}

int main() {
  clock_t t1, t2;
  double s = 0.3; // 0.54799955592586;
  double res;
  double Lexp[dim_L];
  double complex temp_ij[lanes * len_4];
  double complex temp_T_ij[len_4 * lanes];
  double temp_ij_re[len_4 * lanes];
  double temp_ij_im[len_4 * lanes];
  unsigned long long theory_FLOP =
      (3ULL
      + 4ULL * index_max
      + 8ULL * (index_max + ((unsigned long long) index_max * (index_max - 1)) / 2)
      ) * num_it;

  printf("index_max %i\n", index_max);
  printf("\ntheory FLOP = %lld\n", theory_FLOP);

  for (long j = 0; j < dim_L; j++) {
    Lexp[j] = exp(-(0.5 * s * s) * ((j - index_max) * (j - index_max)));
  }

  srand(time(NULL));
  for (long i = 0; i < lanes; i++) {
    for (long j = 0; j < len_4; j++) {
      temp_ij[i * len_4 + j] = (double)(rand() % 10000) / (double)10000 +
                               (double)(rand() % 10000) +
                               I * ((double)(rand() % 10000) / (double)10000 +
                                    (double)(rand() % 10000));
    }
  }

  for (long j = 0; j < len_4; j++) {
    for (long i = 0; i < lanes; i++) {
      temp_ij_re[j * lanes + i] = creal(temp_ij[i * len_4 + j]);
      temp_ij_im[j * lanes + i] = cimag(temp_ij[i * len_4 + j]);
      temp_T_ij[j * lanes + i] = temp_ij[i * len_4 + j];
    }
  }
  // ************************
  // Calculation starts here
  // ************************
#if TEST_NAIVE
  t1 = clock();
  res = test_naive(Lexp, temp_ij);
  t2 = clock();
  printf("\n#1 Naive method\n");
  print_result(theory_FLOP, t1, t2, res);
#endif

#if TEST_UNROLL
  t1 = clock();
  res = test_unroll(Lexp, temp_T_ij);
  t2 = clock();
  printf("\n#2 Unroll step %i\n", step1);
  print_result(theory_FLOP, t1, t2, res);
#endif

#if TEST_UNROLL_SOA
  t1 = clock();
  res = test_unroll_SoA(Lexp, temp_ij_re, temp_ij_im);
  t2 = clock();
  printf("\n#3 SoA, Unroll step %i\n", step2);
  print_result(theory_FLOP, t1, t2, res);
#endif

#if TEST_UNROLL_SOA_SPLIT
  t1 = clock();
  res = test_unroll_SoA_split(Lexp, temp_ij_re, temp_ij_im);
  t2 = clock();
  printf("\n#4 SoA split, Unroll step %i\n", step3);
  print_result(theory_FLOP, t1, t2, res);
#endif
  return 0;
}
