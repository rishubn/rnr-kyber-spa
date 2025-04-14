#include "ntt.h"
#include "params.h"
#include "reduce.h"
#include "stdio.h"
#include <math.h>
#include <stdint.h>
// Code to generate zetas and zetas_inv used in the number-theoretic transform:

#define KYBER_ROOT_OF_UNITY 17
/*
static const uint8_t tree[128] = {
  0, 64, 32, 96, 16, 80, 48, 112, 8, 72, 40, 104, 24, 88, 56, 120,
  4, 68, 36, 100, 20, 84, 52, 116, 12, 76, 44, 108, 28, 92, 60, 124,
  2, 66, 34, 98, 18, 82, 50, 114, 10, 74, 42, 106, 26, 90, 58, 122,
  6, 70, 38, 102, 22, 86, 54, 118, 14, 78, 46, 110, 30, 94, 62, 126,
  1, 65, 33, 97, 17, 81, 49, 113, 9, 73, 41, 105, 25, 89, 57, 121,
  5, 69, 37, 101, 21, 85, 53, 117, 13, 77, 45, 109, 29, 93, 61, 125,
  3, 67, 35, 99, 19, 83, 51, 115, 11, 75, 43, 107, 27, 91, 59, 123,
  7, 71, 39, 103, 23, 87, 55, 119, 15, 79, 47, 111, 31, 95, 63, 127
};
*/
/*
void init_ntt() {
  unsigned int i;
  int16_t tmp[128];

  tmp[0] = MONT;
  for(i=1;i<128;i++)
    tmp[i] = fqmul(tmp[i-1],MONT*KYBER_ROOT_OF_UNITY % KYBER_Q);

  for(i=0;i<128;i++) {
    zetas[i] = tmp[tree[i]];
    if(zetas[i] > KYBER_Q/2)
      zetas[i] -= KYBER_Q;
    if(zetas[i] < -KYBER_Q/2)
      zetas[i] += KYBER_Q;
  }
}
*/

const int16_t zetas[128] = {
    -1044, -758, -359, -1517, 1493, 1422, 287, 202, -171, 622, 1577,
    182, 962, -1202, -1474, 1468, 573, -1325, 264, 383, -829, 1458,
    -1602, -130, -681, 1017, 732, 608, -1542, 411, -205, -1571, 1223,
    652, -552, 1015, -1293, 1491, -282, -1544, 516, -8, -320, -666,
    -1618, -1162, 126, 1469, -853, -90, -271, 830, 107, -1421, -247,
    -951, -398, 961, -1508, -725, 448, -1065, 677, -1275, -1103, 430,
    555, 843, -1251, 871, 1550, 105, 422, 587, 177, -235, -291,
    -460, 1574, 1653, -246, 778, 1159, -147, -777, 1483, -602, 1119,
    -1590, 644, -872, 349, 418, 329, -156, -75, 817, 1097, 603,
    610, 1322, -1285, -1465, 384, -1215, -136, 1218, -1335, -874, 220,
    -1187, -1659, -1185, -1530, -1278, 794, -1510, -854, -870, 478, -108,
    -308, 996, 991, 958, -1460, 1522, 1628};

/*************************************************
 * Name:        fqmul
 *
 * Description: Multiplication followed by Montgomery reduction
 *
 * Arguments:   - int16_t a: first factor
 *              - int16_t b: second factor
 *
 * Returns 16-bit integer congruent to a*b*R^{-1} mod q
 **************************************************/
/*
static int16_t fqmul(int16_t a, int16_t b) {
  return montgomery_reduce((int32_t)a*b);
}

static uint16_t fqmul_u(uint16_t a, uint16_t b) {
  return montgomery_reduce_u((uint32_t)a*b);
}

static uint16_t fqmul_wide_u(uint16_t a, uint16_t b) {
  return montgomery_reduce_wide_u((uint32_t)a*b);
}*/
const uint16_t zetas_u[128] = {
    2285, 2571, 2970, 1812, 1493, 1422, 287, 202, 3158, 622, 1577, 182,
    962, 2127, 1855, 1468, 573, 2004, 264, 383, 2500, 1458, 1727, 3199,
    2648, 1017, 732, 608, 1787, 411, 3124, 1758, 1223, 652, 2777, 1015,
    2036, 1491, 3047, 1785, 516, 3321, 3009, 2663, 1711, 2167, 126, 1469,
    2476, 3239, 3058, 830, 107, 1908, 3082, 2378, 2931, 961, 1821, 2604,
    448, 2264, 677, 2054, 2226, 430, 555, 843, 2078, 871, 1550, 105,
    422, 587, 177, 3094, 3038, 2869, 1574, 1653, 3083, 778, 1159, 3182,
    2552, 1483, 2727, 1119, 1739, 644, 2457, 349, 418, 329, 3173, 3254,
    817, 1097, 603, 610, 1322, 2044, 1864, 384, 2114, 3193, 1218, 1994,
    2455, 220, 2142, 1670, 2144, 1799, 2051, 794, 1819, 2475, 2459, 478,
    3221, 3021, 996, 991, 958, 1869, 1522, 1628};
const uint16_t zetas_wide_u[128] = {
    15601, 2571, 16286, 5141, 8151, 4751, 3616, 3531, 9816, 3951, 4906,
    3511, 4291, 5456, 15171, 8126, 7231, 11991, 10251, 7041, 15816, 8116,
    5056, 13186, 9306, 4346, 4061, 7266, 5116, 411, 13111, 8416, 14539,
    10639, 12764, 4344, 8694, 8149, 13034, 5114, 7174, 9979, 3009, 15979,
    8369, 12154, 6784, 1469, 9134, 3239, 16374, 4159, 10094, 15224, 13069,
    15694, 9589, 7619, 8479, 2604, 13764, 2264, 10664, 2054, 15542, 10417,
    10542, 4172, 5407, 14187, 11537, 10092, 422, 587, 177, 9752, 6367,
    9527, 8232, 4982, 6412, 4107, 7817, 3182, 2552, 4812, 2727, 7777,
    8397, 7302, 2457, 7007, 3747, 6987, 6502, 9912, 14133, 14413, 603,
    7268, 14638, 5373, 5193, 3713, 5443, 3193, 1218, 5323, 9113, 6878,
    15458, 8328, 5473, 5128, 12038, 4123, 5148, 9133, 5788, 478, 13208,
    13008, 10983, 10978, 958, 5198, 14838, 1628};
/*
void init_ntt_wide_u() {
  unsigned int i;
  int16_t tmp[128];

  tmp[0] = WIDE_MONT_U;
  const uint16_t mul_inc = (WIDE_MONT_U * KYBER_ROOT_OF_UNITY) % WIDE_KYBER_Q;
  for(i = 1; i < 128; i++)
  {
    tmp[i] = fqmul_wide_u(tmp[i-1], mul_inc);
  }

  for(i = 0; i < 128; i++)
  {
    zetas_wide_u[i] = tmp[tree[i]];
    if(zetas_wide_u[i] >= WIDE_KYBER_Q)
      zetas_wide_u[i] -= WIDE_KYBER_Q;
  }
}
*/
/*
uint16_t zetas_u[128];
void init_ntt_u() {
  unsigned int i;
  int16_t tmp[128];

  tmp[0] = MONT_U;
  const uint16_t mul_inc = (MONT_U * KYBER_ROOT_OF_UNITY) % KYBER_Q;
  for(i = 1; i < 128; i++)
  {
    tmp[i] = fqmul_u(tmp[i-1], mul_inc);
  }

  for(i = 0; i < 128; i++)
  {
    zetas_u[i] = tmp[tree[i]];
    if(zetas_u[i] >= KYBER_Q)
      zetas_u[i] -= KYBER_Q;
  }
}
*/
/*************************************************
 * Name:        ntt
 *
 * Description: Inplace number-theoretic transform (NTT) in Rq.
 *              input is in standard order, output is in bitreversed order
 *
 * Arguments:   - int16_t r[256]: pointer to input/output vector of elements of
 *Zq
 **************************************************/
void ntt(int16_t r[256])
{
  unsigned int len, start, j, k;
  int16_t t, zeta;

  k = 1;
  for (len = 128; len >= 2; len >>= 1)
  {
    for (start = 0; start < 256; start = j + len)
    {
      zeta = zetas[k++];
      for (j = start; j < start + len; j++)
      {
        // printf("j: %d len: %d\nr[j]: %4d %4d\nr[j+len] %4d %4d\nzeta %d\n", j,
        //        len, r[j],
        //        r[j] % 3329 < 0 ? ((r[j] % 3329) + 3329) : r[j] % 3329,
        //        r[j + len],
        //        r[j + len] % 3329 < 0 ? ((r[j + len] % 3329) + 3329)
        //                              : r[j + len] % 3329,
        //        zeta);
        t = fqmul(zeta, r[j + len]);
        // printf("t: %d\n", t);
        r[j + len] = r[j] - t;
        r[j] = r[j] + t;
        // printf("\n\nj: %d len: %d\nr[j]: %4d %4d\nr[j+len] %4d %4d\nzeta "
        //        "%d\n--------\n",
        //        j, len, r[j],
        //        r[j] % 3329 < 0 ? ((r[j] % 3329) + 3329) : r[j] % 3329,
        //        r[j + len],
        //        r[j + len] % 3329 < 0 ? ((r[j + len] % 3329) + 3329)
        //                              : r[j + len] % 3329,
        //        zeta);
      }
    }
    // printf("fwd len: %d:\n", len);
    // for (uint32_t i = 0; i < 256; i++)
    // {
    //   printf("%4d ", r[i] % 3329 < 0 ? ((r[i] % 3329) + 3329) : r[i] % 3329);
    //   if (i % 16 == 15)
    //     printf("\n");
    // }
  }
  // printf("fwd final:\n");
  // for (uint32_t i = 0; i < 256; i++)
  // {
  //   printf("%4d ", r[i] % 3329 < 0 ? ((r[i] % 3329) + 3329) : r[i] % 3329);
  //   if (i % 16 == 15)
  //     printf("\n");
  // }
}
void ntt_wide_s(int16_t r[256])
{
  unsigned int len, start, j, k;
  int16_t t, zeta;

  k = 1;
  for (len = 128; len >= 2; len >>= 1)
  {
    for (start = 0; start < 256; start = j + len)
    {
      zeta = zetas[k++];
      for (j = start; j < start + len; j++)
      {
        /*printf("j: %d len: %d\nr[j]: %4d %4d\nr[j+len] %4d %4d\nzeta %d\n", j,
               len, r[j],
               r[j] % 3329 < 0 ? ((r[j] % 3329) + 3329) : r[j] % 3329,
               r[j + len],
               r[j + len] % 3329 < 0 ? ((r[j + len] % 3329) + 3329)
                                     : r[j + len] % 3329,
               zeta);
*/
        t = fqmul_wide_s(zeta, r[j + len]);
        r[j + len] = barrett_reduce_wide_s(r[j] - t);
        r[j] = barrett_reduce_wide_s(r[j] + t);
      }
    }
    // printf("fwd len: %d:\n", len);
    // for (uint32_t i = 0; i < 256; i++)
    // {
    //   printf("%4d ", r[i] % 3329 < 0 ? ((r[i] % 3329) + 3329) : r[i] % 3329);
    //   if (i % 16 == 15)
    //     printf("\n");
    // }
  }

  for (uint32_t i = 0; i < 256; i++)
    r[i] = barrett_reduce_wide_s(r[i]);
  // printf("fwd final:\n");
  // for (uint32_t i = 0; i < 256; i++)
  // {
  //   printf("%4d ", r[i]);
  //   if (i % 16 == 15)
  //     printf("\n");
  // }
}
void ntt_u(uint16_t r[256])
{
  unsigned int len, start, j, k;
  uint16_t t, zeta;

  k = 1;
  for (len = 128; len >= 2; len >>= 1)
  {
    for (start = 0; start < 256; start = j + len)
    {
      zeta = zetas_u[k++];
      for (j = start; j < start + len; j++)
      {
        t = fqmul_u(zeta, r[j + len]); // 2q + 0,phi*q; 0, 2+phi*q < 2**16
        r[j + len] = 2 * KYBER_Q + r[j] - t;
        r[j] = r[j] + t;
      }
    }
    // printf("fwd len: %d:\n", len);
    // for (uint32_t i = 0; i < 256; i++)
    // {
    //   printf("%4d ", r[i]);
    //   if (i % 16 == 15)
    //     printf("\n");
    // }
  }

  for (uint32_t i = 0; i < 256; i++)
    r[i] = barrett_reduce_u(r[i]);

  // printf("fwd final:\n");
  // for (uint32_t i = 0; i < 256; i++)
  // {
  //   printf("%4d ", r[i]);
  //   if (i % 16 == 15)
  //     printf("\n");
  // }
}
void ntt_wide_u9(uint16_t r[256])
{
  unsigned int len, start, j, k;
  uint16_t t, zeta;
  k = 1;
  //  c = 0;
  for (len = 128; len >= 2; len >>= 1)
  {
    for (start = 0; start < 256; start = j + len)
    {
      zeta = zetas_u[k++];
      for (j = start; j < start + len; j++)
      {
        t = fqmul_wide_u9(zeta, r[j + len]);
        r[j + len] = barrett_reduce_wide_u9((10 * KYBER_Q) + r[j] - t);
        r[j] = barrett_reduce_wide_u9(r[j] + t);
      }
    }
  }
  for (uint32_t i = 0; i < 256; i++)
    r[i] = barrett_reduce_wide_u9(r[i]);
}

void ntt_wide_u(uint16_t r[256])
{
  unsigned int len, start, j, k;
  uint16_t t, zeta;
  k = 1;
  //  c = 0;
  for (len = 128; len >= 2; len >>= 1)
  {
    for (start = 0; start < 256; start = j + len)
    {
      zeta = zetas_u[k++];
      for (j = start; j < start + len; j++)
      {
        t = fqmul_wide_u(zeta, r[j + len]);
        r[j + len] = barrett_reduce_wide_u((2 * WIDE_KYBER_Q) + r[j] - t);
        r[j] = barrett_reduce_wide_u(r[j] + t);
      }
    }
    // c++;
    //  printf("fwd len: %d:\n", len);
    //  for (uint32_t i = 0; i < 256; i++)
    //  {
    //    printf("%4d ", r[i]);
    //    if (i % 16 == 15)
    //      printf("\n");
    //  }
  }

  for (uint32_t i = 0; i < 256; i++)
    r[i] = barrett_reduce_wide_u(r[i]);

  // printf("fwd final:\n");
  // for (uint32_t i = 0; i < 256; i++)
  // {
  //   printf("%4d ", r[i] % 3329);
  //   if (i % 16 == 15)
  //     printf("\n");
  // }
  // printf("DONE\n");
}
/*************************************************
 * Name:        invntt_tomont
 *
 * Description: Inplace inverse number-theoretic transform in Rq and
 *              multiplication by Montgomery factor 2^16.
 *              Input is in bitreversed order, output is in standard order
 *
 * Arguments:   - int16_t r[256]: pointer to input/output vector of elements of
 *Zq
 **************************************************/
void invntt_full(int16_t r[256], int16_t full[256][9])
{
  unsigned int start, len, j, k;
  int16_t t, zeta;
  const int16_t f = 1441; // mont^2/128

  k = 127;
  int c = 0;
  for (int i = 0; i < 256; i++)
  {
    full[i][0] = r[i];
  }

  c++;

  for (len = 2; len <= 128; len <<= 1)
  {
    for (start = 0; start < 256; start = j + len)
    {
      zeta = zetas[k--];
      for (j = start; j < start + len; j++)
      {
        // printf("j: %d len: %d r[j]: %d r[j+len] %d zeta %d\n", j, len, r[j],
        //        r[j + len], zeta);
        t = r[j];
        //  printf("1: t: %d\n", t);
        r[j] = barrett_reduce(t + r[j + len]);
        //   printf("2: r[j]: %d\n", r[j]);
        r[j + len] = r[j + len] - t;
        //    printf("3: r[j+len]: %d\n", r[j + len]);
        r[j + len] = fqmul(zeta, r[j + len]);
        //   printf("4: r[j+len]: %d\n", r[j + len]);
        full[j][c] = r[j];
        full[j + len][c] = r[j + len];
      }
    }

    // printf("inv len: %d:\n", len);
    // for (uint32_t i = 0; i < 256; i++)
    // {
    //   printf("%4d ", r[i] % 3329 < 0 ? ((r[i] % 3329) + 3329) : r[i] % 3329);
    //   if (i % 16 == 15)
    //     printf("\n");
    // }
    c++;
  }

  for (j = 0; j < 256; j++)
  {
    r[j] = fqmul(r[j], f);
    full[j][8] = r[j];
  }

  // printf("inv final:\n");
  // for (uint32_t i = 0; i < 256; i++)
  // {
  //   printf("%4d ", r[i] % 3329 < 0 ? ((r[i] % 3329) + 3329) : r[i] % 3329);
  //   if (i % 16 == 15)
  //     printf("\n");
  // }
}
void invntt(int16_t r[256])
{
  unsigned int start, len, j, k;
  int16_t t, zeta;
  const int16_t f = 1441; // mont^2/128

  k = 127;
  for (len = 2; len <= 128; len <<= 1)
  {
    for (start = 0; start < 256; start = j + len)
    {
      zeta = zetas[k--];
      for (j = start; j < start + len; j++)
      {
        // printf("j: %d len: %d r[j]: %d r[j+len] %d zeta %d\n", j, len, r[j],
        //        r[j + len], zeta);
        t = r[j];
        r[j] = barrett_reduce(t + r[j + len]);
        r[j + len] = r[j + len] - t;
        r[j + len] = fqmul(zeta, r[j + len]);
      }
    }
    // printf("inv len: %d:\n", len);
    // for (uint32_t i = 0; i < 256; i++)
    // {
    //   printf("%4d ", r[i] % 3329 < 0 ? ((r[i] % 3329) + 3329) : r[i] % 3329);
    //   if (i % 16 == 15)
    //     printf("\n");
    // }
  }

  for (j = 0; j < 256; j++)
    r[j] = fqmul(r[j], f);

  // printf("inv final:\n");
  // for (uint32_t i = 0; i < 256; i++)
  // {
  //   printf("%4d ", r[i] % 3329 < 0 ? ((r[i] % 3329) + 3329) : r[i] % 3329);
  //   if (i % 16 == 15)
  //     printf("\n");
  // }
}

void invntt_wide_s_full(int16_t r[256], int16_t full[256][9])
{
  unsigned int start, len, j, k;
  int16_t t, zeta;
  const int16_t f = 1441; // mont^2/128

  for (int i = 0; i < 256; i++)
  {
    full[i][0] = r[i];
  }
  int c = 1;
  k = 127;
  for (len = 2; len <= 128; len <<= 1)
  {
    for (start = 0; start < 256; start = j + len)
    {
      zeta = zetas[k--];
      for (j = start; j < start + len; j++)
      {
        t = r[j];
        r[j] = barrett_reduce_wide_s(t + r[j + len]);
        r[j + len] = (r[j + len] - t);
        r[j + len] = fqmul_wide_s(zeta, r[j + len]);
        full[j][c] = r[j];
        full[j + len][c] = r[j + len];
      }
    }
    c++;
  }

  for (j = 0; j < 256; j++)
  {
    r[j] = fqmul_wide_s(r[j], f);
    full[j][8] = r[j];
  }
}

void invntt_wide_s(int16_t r[256])
{
  unsigned int start, len, j, k;
  int16_t t, zeta;
  const int16_t f = 1441; // mont^2/128

  k = 127;
  for (len = 2; len <= 128; len <<= 1)
  {
    for (start = 0; start < 256; start = j + len)
    {
      zeta = zetas[k--];
      for (j = start; j < start + len; j++)
      {
        t = r[j];
        r[j] = barrett_reduce_wide_s(t + r[j + len]);
        r[j + len] = (r[j + len] - t);
        r[j + len] = fqmul_wide_s(zeta, r[j + len]);
      }
    }
    // printf("inv len: %d:\n", len);
    // for (uint32_t i = 0; i < 256; i++)
    // {
    //   printf("%4d ", r[i]);
    //   if (i % 16 == 15)
    //     printf("\n");
    // }
  }

  for (j = 0; j < 256; j++)
    r[j] = fqmul_wide_s(r[j], f);

  // printf("inv final:\n");
  // for (uint32_t i = 0; i < 256; i++)
  // {
  //   printf("%4d ", r[i] % 3329 < 0 ? ((r[i] % 3329) + 3329) : r[i] % 3329);
  //   if (i % 16 == 15)
  //     printf("\n");
  // }
}
void invntt_u(uint16_t r[256])
{
  unsigned int start, len, j, k;
  uint16_t t, zeta;
  k = 127;
  for (len = 2; len <= 128; len <<= 1)
  {
    for (start = 0; start < 256; start = j + len)
    {
      zeta = zetas_u[k--];
      for (j = start; j < start + len; j++)
      {
        t = r[j];
        r[j] = barrett_reduce_u(t + r[j + len]);
        r[j + len] = barrett_reduce_u(2 * KYBER_Q + r[j + len] - t);
        r[j + len] = fqmul_u(zeta, r[j + len]);
      }
    }
    // printf("inv len: %d:\n", len);
    // for (uint32_t i = 0; i < 256; i++)
    // {
    //   printf("%4d ", r[i]);
    //   if (i % 16 == 15)
    //     printf("\n");
    // }
  }

  for (j = 0; j < 256; j++)
    r[j] = fqmul_u(r[j], INTT_CONST_U);

  // printf("inv final:\n");
  // for (uint32_t i = 0; i < 256; i++)
  // {
  //   printf("%4d ", r[i]);
  //   if (i % 16 == 15)
  //     printf("\n");
  // }
}
void invntt_u_full(uint16_t r[256], uint16_t full[256][9])
{
  unsigned int start, len, j, k;
  uint16_t t, zeta;
  k = 127;
  for (int i = 0; i < 256; i++)
  {
    full[i][0] = r[i];
  }
  int c = 1;

  for (len = 2; len <= 128; len <<= 1)
  {
    for (start = 0; start < 256; start = j + len)
    {
      zeta = zetas_u[k--];
      for (j = start; j < start + len; j++)
      {
        t = r[j];
        r[j] = barrett_reduce_u(t + r[j + len]);
        r[j + len] = barrett_reduce_u(2 * KYBER_Q + r[j + len] - t);
        r[j + len] = fqmul_u(zeta, r[j + len]);
        full[j][c] = r[j];
        full[j + len][c] = r[j + len];
      }
    }
    c++;
    // printf("inv len: %d:\n", len);
    // for (uint32_t i = 0; i < 256; i++)
    // {
    //   printf("%4d ", r[i]);
    //   if (i % 16 == 15)
    //     printf("\n");
    // }
  }

  for (j = 0; j < 256; j++)
  {
    r[j] = fqmul_u(r[j], INTT_CONST_U);
    full[j][8] = r[j];
  }
  // printf("inv final:\n");
  // for (uint32_t i = 0; i < 256; i++)
  // {
  //   printf("%4d ", r[i]);
  //   if (i % 16 == 15)
  //     printf("\n");
  // }
}
void invntt_wide_u(uint16_t r[256])
{
  unsigned int start, len, j, k;
  uint16_t t, zeta;

  k = 127;
  // printf("inv start\n");
  // for (uint32_t i = 0; i < 256; i++)
  // {
  //   printf("%4d ", r[i]);
  //   if (i % 16 == 15)
  //     printf("\n");
  // }
  for (len = 2; len <= 128; len <<= 1)
  {
    for (start = 0; start < 256; start = j + len)
    {
      zeta = zetas_u[k--];
      for (j = start; j < start + len; j++)
      {
        t = r[j];
        r[j] = barrett_reduce_wide_u(t + r[j + len]);
        r[j + len] = barrett_reduce_wide_u(2 * WIDE_KYBER_Q + r[j + len] - t);
        r[j + len] = fqmul_wide_u(zeta, r[j + len]);
      }
    }
    // printf("inv len: %d:\n", len);
    // for (uint32_t i = 0; i < 256; i++)
    // {
    //   printf("%4d ", r[i] % 3329);
    //   if (i % 16 == 15)
    //     printf("\n");
    // }
  }

  for (j = 0; j < 256; j++)
    r[j] = fqmul_wide_u(r[j], INTT_CONST_U);

  // printf("inv final:\n");
  // for (uint32_t i = 0; i < 256; i++)
  // {
  //   printf("%4d ", r[i] % 3329);
  //   if (i % 16 == 15)
  //     printf("\n");
  // }
}
void invntt_wide_u9(uint16_t r[256])
{
  unsigned int start, len, j, k;
  uint16_t t, zeta;

  k = 127;

  for (len = 2; len <= 128; len <<= 1)
  {
    for (start = 0; start < 256; start = j + len)
    {
      zeta = zetas_u[k--];
      for (j = start; j < start + len; j++)
      {
        t = r[j];
        r[j] = barrett_reduce_wide_u9(t + r[j + len]);
        r[j + len] = barrett_reduce_wide_u9(10 * KYBER_Q + r[j + len] - t);
        r[j + len] = fqmul_wide_u9(zeta, r[j + len]);
      }
    }
  }

  for (j = 0; j < 256; j++)
    r[j] = fqmul_wide_u9(r[j], INTT_CONST_U);
}
void invntt_wide_u_full(uint16_t r[256], uint16_t full[256][9])
{
  unsigned int start, len, j, k, c;
  uint16_t t, zeta;

  k = 127;
  for (int i = 0; i < 256; i++)
  {
    full[i][0] = r[i];
  }
  c = 1;
  // printf("inv start\n");
  // for (uint32_t i = 0; i < 256; i++)
  // {
  //   printf("%4d ", r[i]);
  //   if (i % 16 == 15)
  //     printf("\n");
  // }
  for (len = 2; len <= 128; len <<= 1)
  {
    for (start = 0; start < 256; start = j + len)
    {
      zeta = zetas_u[k--];
      for (j = start; j < start + len; j++)
      {
        t = r[j];
        r[j] = barrett_reduce_wide_u(t + r[j + len]);
        r[j + len] = barrett_reduce_wide_u(2 * WIDE_KYBER_Q + r[j + len] - t);
        r[j + len] = fqmul_wide_u(zeta, r[j + len]);
        full[j][c] = r[j];
        full[j + len][c] = r[j + len];
      }
    }
    c++;
    // printf("inv len: %d:\n", len);
    // for (uint32_t i = 0; i < 256; i++)
    // {
    //   printf("%4d ", r[i] % 3329);
    //   if (i % 16 == 15)
    //     printf("\n");
    // }
  }

  for (j = 0; j < 256; j++)
  {
    r[j] = fqmul_wide_u(r[j], INTT_CONST_U);
    full[j][8] = r[j];
  }
}

// printf("inv final:\n");
// for (uint32_t i = 0; i < 256; i++)
// {
//   printf("%4d ", r[i] % 3329);
//   if (i % 16 == 15)
//     printf("\n");
// }

/*************************************************
 * Name:        basemul
 *
 * Description: Multiplication of polynomials in Zq[X]/(X^2-zeta)
 *              used for multiplication of elements in Rq in NTT domain
 *
 * Arguments:   - int16_t r[2]: pointer to the output polynomial
 *              - const int16_t a[2]: pointer to the first factor
 *              - const int16_t b[2]: pointer to the second factor
 *              - int16_t zeta: integer defining the reduction polynomial
 **************************************************/
void basemul(int16_t r[2], const int16_t a[2], const int16_t b[2],
             int16_t zeta)
{
  r[0] = fqmul(a[1], b[1]);
  r[0] = fqmul(r[0], zeta);
  r[0] += fqmul(a[0], b[0]);
  r[1] = fqmul(a[0], b[1]);
  r[1] += fqmul(a[1], b[0]);
}

void basemul_u(uint16_t r[2], const uint16_t a[2], const uint16_t b[2],
               uint16_t zeta)
{
  r[0] = fqmul_u(a[1], b[1]);
  r[0] = fqmul_u(r[0], zeta);
  r[0] += fqmul_u(a[0], b[0]);
  r[1] = fqmul_u(a[0], b[1]);
  r[1] += fqmul_u(a[1], b[0]);
}
