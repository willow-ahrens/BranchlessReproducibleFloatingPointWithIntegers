/*
Copyright (c) 2016, Los Alamos National Security, LLC

All rights reserved.

Copyright 2016. Los Alamos National Security, LLC. This software was produced under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL), which is operated by Los Alamos National Security, LLC for the U.S. Department of Energy. The U.S. Government has rights to use, reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce derivative works, such modified software should be clearly marked, so as not to confuse it with the version available from LANL.

Additionally, redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1.      Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2.      Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3.      Neither the name of Los Alamos National Security, LLC, Los Alamos National Laboratory, LANL, the U.S. Government, nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <omp.h>

union double_bits {double d; uint64_t b;};

typedef struct accum {
  float exception;
  int i;
  uint64_t p_0;
  uint64_t p_1;
  int64_t c_0;
  int64_t c_1;
} accum;

inline accum double_to_accum (double x) {
  accum y;
  union double_bits db;
  db.d = x;
  uint64_t xi = db.b;
  unsigned exp = (xi >> 52) & 0x7FF;
  uint64_t mant = (xi & 0x000FFFFFFFFFFFFFLL) | 0x0010000000000000LL;
  mant *= (exp != 0);
  y.exception = x * (exp == 0x7FF);
  y.c_0 = (xi >> 63);
  y.i = (exp / 52) + 1;
  unsigned shift = y.i * 52 - exp;
  int sign = 1 - (y.c_0 << 1);
  y.p_0 = sign * (mant >> shift);
  y.p_1 = sign * ((mant << (52 - shift)) & 0x000FFFFFFFFFFFFFLL);
  y.c_0 = -y.c_0;
  y.c_1 = y.c_0;
  return y;
}

inline accum accum_plus_accum (accum x_0, accum x_1) {
  accum y;
  y.exception = x_0.exception + x_1.exception;
  y.i = x_0.i - x_1.i;
  y.i &= (~y.i)>>31;
  y.i += x_1.i;
  uint64_t m_0_0 = 0xFFFFFFFFFFFFFFFFLL * (x_0.i == y.i);
  uint64_t m_1_0 = 0xFFFFFFFFFFFFFFFFLL * (x_1.i == y.i);
  uint64_t m_0_1 = 0xFFFFFFFFFFFFFFFFLL * (x_0.i == y.i - 1);
  uint64_t m_1_1 = 0xFFFFFFFFFFFFFFFFLL * (x_1.i == y.i - 1);
  uint64_t p_0_0 = (x_0.p_0 & m_0_0);
  uint64_t p_0_1 = (x_1.p_0 & m_1_0);
  uint64_t p_1_0 = (x_0.p_0 & m_0_1) + (x_0.p_1 & m_0_0);
  uint64_t p_1_1 = (x_1.p_0 & m_1_1) + (x_1.p_1 & m_1_0);
  y.p_0 = p_0_0 + p_0_1;
  y.c_0 = (y.p_0 < p_0_0);
  y.p_1 = p_1_0 + p_1_1;
  y.c_1 = (y.p_1 < p_1_0);
  y.c_0 += (x_0.c_0 & m_0_0) + (x_1.c_0 & m_1_0);
  y.c_0 += (x_0.c_0 & m_0_0) + (x_1.c_0 & m_1_0);
  y.c_1 += (x_0.c_0 & m_0_1) + (x_1.c_0 & m_1_1)
          + (x_0.c_1 & m_0_0) + (x_1.c_1 & m_1_0);
  return y;
}



double accum_to_double (accum x) {
  int64_t y_0 = -(x.c_1 >> 63);
  uint64_t y_1 = x.c_1;
  uint64_t y_2 = x.p_1;

  int64_t z_0 = x.c_0 / (1 << (64 - 52));
  uint64_t z_1 = (x.c_0 << 52) + x.p_0 / (1 << (64 - 52));
  uint64_t z_2 = x.p_0 << 52;

  y_2 += z_2;
  y_1 += (y_2 < z_2);
  y_0 += (y_1 < (y_2 < z_2));
  y_1 += z_1;
  y_0 += (y_1 < z_1);
  y_0 += z_0;

  int sign = 1 - ((y_0 >> 62) & 0x2);

  double y;
  y = ((double)y_0) * (4294967296.0 * 4294967296.0 * 4294967296.0 * 4294967296.0);
  y += ((double)(y_1 >> 32)) * (4294967296.0 * 4294967296.0 * 4294967296.0);
  y += ((double)(y_1 & 0xFFFFFFFF)) * (4294967296.0 * 4294967296.0);
  y += ((double)(y_2 >> 32)) * 4294967296.0;
  y += (double)(y_2 & 0xFFFFFFFF);
  y *= (x.exception == 0.);
  y *= sign * ldexp(0.5, (x.i - 2) * 52 - (DBL_MAX_EXP - 2));
  return y + x.exception;
}

int main(int argc, char **argv){
  int test = 10000;
  int trials = 10000;
  double x[test];
  for(int i; i < test; i++){
    x[i] = ((double)rand())/RAND_MAX;
  }
  double rep_time = -omp_get_wtime();
  double rep_sum;
  accum asum;
  for(int j = 0; j < trials; j++){
    asum.exception = 0;
    asum.i = 0;
    asum.p_0 = 0;
    asum.p_1 = 0;
    asum.c_0 = 0;
    asum.c_1 = 0;
    for(int i = 0; i < test; i++){
      asum = accum_plus_accum(double_to_accum(x[i]), asum);
    }
    rep_sum = accum_to_double(asum);
  }
  rep_time += omp_get_wtime();
  rep_time /= trials;

  double ref_time = -omp_get_wtime();
  double ref_sum = 0;
  for(int j = 0; j < trials; j++){
    ref_sum = 0;
    for(int i = 0; i < test; i++){
      ref_sum += x[i];
    }
  }
  ref_time += omp_get_wtime();
  ref_time /= trials;

  double kah_time = -omp_get_wtime();
  double kah_sum = 0;
  double kah_c = 0;
  double kah_y;
  double kah_t;
  for(int j = 0; j < trials; j++){
    kah_sum = 0;
    kah_c = 0;
    for(int i = 0; i < test; i++){
      kah_y = x[i] - kah_c;
      kah_t = kah_sum + kah_y;
      kah_c = (kah_t - kah_sum) - kah_y;
      kah_sum = kah_t;
    }
  }
  kah_time += omp_get_wtime();
  kah_time /= trials;

  printf("rep_time: %g\n", rep_time);
  printf("ref_time: %g\n", ref_time);
  printf("kah_time: %g\n", kah_time);
  printf("\n");
  printf(" ref_factor: %g\n", rep_time/ref_time);
  printf(" kah_factor: %g\n", rep_time/kah_time);
  printf("\n");
  printf(" rep_sum: %g\n", rep_sum);
  printf(" ref_sum: %g\n", ref_sum);
  printf(" kah_sum: %g\n", kah_sum);

  return 0;
}
