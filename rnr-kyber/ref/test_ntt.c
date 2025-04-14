#include "ntt.h"
#include "params.h"
#include "randombytes.h"
#include "reduce.h"
#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define assert(x)      \
  {                    \
    if (!(x))          \
    {                  \
      *((int *)0) = 0; \
    }                  \
  }
#define NTESTS 1000

#define HW_SIZE 17
#define REDUNDANCY 9
#define DOMAIN_SIZE (REDUNDANCY * KYBER_Q)

static int test_ntt_u(void)
{
#define SIZE 256
  uint16_t input[SIZE];
  uint16_t output[SIZE];
  // init_ntt_u();
  // printf("const uint16_t zetas_u[128] = {\n");
  // for (int i = 0; i < 128; i++)
  // {
  //   if (i % 8 == 0) printf("    ");
  //   printf("%hu", zetas_u[i]);
  //   if (i != 127) printf(", ");
  //   if (i % 8 == 7) printf("\n");
  // }
  // printf("};\n");

  for (uint32_t test = 0; test < NTESTS; test++)
  {
    // this is not uniform, but I do not care for testing
    randombytes((uint8_t *)input, SIZE * sizeof(uint16_t));
    for (uint32_t i = 0; i < SIZE; i++)
      input[i] = 0;
    memcpy(output, input, SIZE * sizeof(uint16_t));

    ntt_u(output);
    invntt_u(output);

    for (uint32_t i = 0; i < SIZE; i++)
    {
      output[i] = ((uint32_t)output[i] * (uint32_t)MINV_U) % KYBER_Q;
    }

    if (memcmp(input, output, SIZE * sizeof(uint16_t)) != 0)
    {
      printf("ERROR test %d: invntt_u(ntt_u(x)) != x\n", test);
      for (int i = 0; i < SIZE; i++)
        printf("%4hu %4hu\n", input[i], output[i]);

      return 1;
    }
  }

  return 0;
}
static int test_ntt_wide_u(void)
{
#define SIZE 256
  uint16_t input[SIZE];
  uint16_t output[SIZE];
  printf("START\n");
  // init_ntt_u();
  // printf("const uint16_t zetas_u[128] = {\n");
  // for (int i = 0; i < 128; i++)
  // {
  //   if (i % 8 == 0) printf("    ");
  //   printf("%hu", zetas_u[i]);
  //   if (i != 127) printf(", ");
  //   if (i % 8 == 7) printf("\n");
  // }
  // printf("};\n");

  for (uint32_t test = 0; test < NTESTS; test++)
  {
    // this is not uniform, but I do not care for testing
    randombytes((uint8_t *)input, SIZE * sizeof(uint16_t));
    for (uint32_t i = 0; i < SIZE; i++)
      input[i] %= WIDE_KYBER_Q;
    memcpy(output, input, SIZE * sizeof(uint16_t));

    ntt_wide_u(output);
    invntt_wide_u(output);

    for (uint32_t i = 0; i < SIZE; i++)
    {
      output[i] = ((uint32_t)output[i] * (uint32_t)WIDE_MINV_U) % KYBER_Q;
      input[i] %= KYBER_Q;
    }
    if (memcmp(input, output, SIZE * sizeof(uint16_t)) != 0)
    {
      printf("ERROR test %d: invntt_wide_u(ntt_wide_u(x)) != x\n", test);
      for (int i = 0; i < SIZE; i++)
        printf("%4hu %4hu\n", input[i], output[i]);

      return 1;
    }
  }

  return 0;
}
static int test_ntt_wide_s(void)
{
#define SIZE 256

  int16_t input[SIZE] = {
      27378,
      2798,
      -20900,
      -27212,
      13305,
      3808,
      -8845,
      -6183,
      -7170,
      -16563,
      8656,
      -15535,
      7520,
      -26580,
      -9889,
      -6526,
      -21111,
      987,
      -10677,
      13664,
      -4433,
      9406,
      -25632,
      16635,
      -7539,
      23541,
      -8034,
      20066,
      9948,
      4063,
      9233,
      1166,
      29251,
      -7659,
      23829,
      20927,
      -16289,
      -1895,
      24132,
      7394,
      -27697,
      12880,
      -18582,
      -18098,
      4711,
      -16215,
      -1584,
      6032,
      -18989,
      -28883,
      -397,
      -2507,
      671,
      25965,
      -1149,
      -15177,
      -11852,
      28801,
      12436,
      -322,
      27800,
      3628,
      -27286,
      -18827,
      12448,
      -2340,
      -12659,
      1634,
      14285,
      21671,
      18667,
      -26833,
      6118,
      15927,
      29928,
      -19193,
      -6805,
      -9237,
      18725,
      -2410,
      29365,
      14092,
      -25296,
      21416,
      7768,
      1849,
      21458,
      -5720,
      20416,
      22798,
      2713,
      23402,
      1843,
      -19685,
      -7566,
      12845,
      17854,
      8471,
      -13665,
      3323,
      18885,
      -16380,
      -1735,
      29828,
      -15453,
      24874,
      -4284,
      -6096,
      22157,
      -6428,
      3310,
      18774,
      -14830,
      7998,
      -14131,
      8121,
      -25720,
      15511,
      25469,
      15946,
      15577,
      -9956,
      5242,
      -23448,
      1065,
      17580,
      18497,
      15071,
      -29015,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,

      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,

      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,

      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,

      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,

      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,

      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,

      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
  };
  int16_t output[SIZE];
  for (uint32_t test = 0; test < NTESTS; test++)
  {
    // this is not uniform, but I do not care for testing
    randombytes((uint8_t *)input, SIZE * sizeof(int16_t));
    for (uint32_t i = 0; i < SIZE; i++)
    {
      input[i] %= WIDE_KYBER_Q_S;
      if (input[i] > WIDE_KYBER_Q_S / 2)
        input[i] -= WIDE_KYBER_Q_S;
      if (input[i] < -WIDE_KYBER_Q_S / 2)
        input[i] += WIDE_KYBER_Q_S;
    }

    // randombytes((uint8_t *)input, 129 * sizeof(int16_t));
    // for (uint32_t i = 0; i < 129; i++) {
    //   input[i] %= WIDE_KYBER_Q_S;
    // }
    for (uint32_t i = 0; i < 256; i++)
    {

      printf("%4d, ", input[i]);
      if (i % 16 == 0)
        printf("\n");
    }
    printf("\n");
    memcpy(output, input, SIZE * sizeof(int16_t));

    ntt_wide_s(output);
    invntt_wide_s(output);

    for (int i = 0; i < SIZE; i++)
      printf("%d %d\n", input[i], output[i] % KYBER_Q);

    for (uint32_t i = 0; i < SIZE; i++)
    {
      output[i] = (output[i] <= 0) ? WIDE_KYBER_Q_S + output[i] : output[i];
      output[i] = ((uint32_t)output[i] * (uint32_t)WIDE_MINV_S) % KYBER_Q;
      input[i] %= KYBER_Q;
      if (input[i] > KYBER_Q / 2)
        printf("%4d, ", input[i] - KYBER_Q);
      else if (input[i] < -KYBER_Q / 2)
        printf("%4d, ", input[i] + KYBER_Q);
      else
        printf("%4d, ", input[i]);
      if (i % 16 == 0)
        printf("\n");
      if (output[i] > KYBER_Q / 2)
        output[i] -= KYBER_Q;
      if (input[i] < 0)
        input[i] += KYBER_Q;
      if (output[i] < 0)
        output[i] += KYBER_Q;
    }

    if (memcmp(input, output, SIZE * sizeof(int16_t)) != 0)
    {
      printf("ERROR test %d: invntt_wide_s(ntt_wide_s(x)) != x\n", test);
      for (int i = 0; i < SIZE; i++)
        printf("%d %d\n", input[i], output[i]);

      return 1;
    }
  }
  return 0;
}
static int test_ntt(void)
{
#define SIZE 256
  int16_t input[SIZE] = {
      -582,
      884,
      -1047,
      -746,
      528,
      735,
      668,
      -700,
      1215,
      -1637,
      -1020,
      -242,
      -1146,
      268,
      -539,
      -938,
      465,
      335,
      1114,
      184,
      1150,
      -533,
      -30,
      -1340,
      -133,
      -804,
      -367,
      1662,
      867,
      1379,
      -1436,
      889,
      -269,
      -1125,
      -939,
      997,
      495,
      -573,
      89,
      174,
      -1421,
      523,
      674,
      874,
      -422,
      756,
      -24,
      -620,
      -576,
      910,
      1275,
      970,
      670,
      679,
      -1555,
      -145,
      -1564,
      -1031,
      -1460,
      -1346,
      344,
      781,
      1090,
      924,
      747,
      330,
      1483,
      1310,
      258,
      -900,
      -1571,
      784,
      130,
      -1309,
      -788,
      172,
      -1387,
      1189,
      1328,
      -1637,
      -797,
      493,
      20,
      1357,
      -837,
      1049,
      814,
      1375,
      -315,
      -558,
      142,
      327,
      -953,
      -374,
      -810,
      1600,
      1186,

      -704,
      1450,
      701,
      663,
      467,
      509,
      -888,
      -1222,
      -1408,
      -100,
      317,
      343,
      -1058,
      715,
      132,
      340,
      907,
      -298,
      -879,
      -1261,
      -1510,
      597,
      921,
      -216,
      -1130,
      663,
      900,
      -961,
      -1465,
      777,
      1090,
      -348,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,

      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,

      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,

      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,

      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,

      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,

      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,

      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
  };
  int16_t output[SIZE];
  // init_ntt_u();
  // printf("const uint16_t zetas_u[128] = {\n");
  // for (int i = 0; i < 128; i++)
  // {
  //   if (i % 8 == 0) printf("    ");
  //   printf("%hu", zetas_u[i]);
  //   if (i != 127) printf(", ");
  //   if (i % 8 == 7) printf("\n");
  // }
  // printf("};\n");

  for (uint32_t test = 0; test < NTESTS; test++)
  {
    // this is not uniform, but I do not care for testing
    // randombytes((uint8_t *)input, SIZE * sizeof(int16_t));
    for (uint32_t i = 0; i < SIZE; i++)
    {
      if (input[i] > KYBER_Q / 2)
        input[i] -= KYBER_Q;
      if (input[i] < -KYBER_Q / 2)
        input[i] += KYBER_Q;
    }

    memcpy(output, input, SIZE * sizeof(int16_t));

    ntt(output);
    invntt(output);

    for (uint32_t i = 0; i < SIZE; i++)
    {
      output[i] = (output[i] <= 0) ? KYBER_Q + output[i] : output[i];
      output[i] = ((uint32_t)output[i] * (uint32_t)MINV_U) % KYBER_Q;
      if (output[i] > KYBER_Q / 2)
        output[i] -= KYBER_Q;
    }

    if (memcmp(input, output, SIZE * sizeof(int16_t)) != 0)
    {
      printf("ERROR test %d: invntt_u(ntt_u(x)) != x\n", test);
      for (int i = 0; i < SIZE; i++)
        printf("%04hx %04hx\n", input[i], output[i]);

      return 1;
    }
  }

  return 0;
}

static int test_barrett_reduction_u(void)
{

  for (uint32_t i = 0; i < 1 << 16; i++)
  {
    uint16_t res_mine = barrett_reduce_u(i);
    uint16_t res_good = i % (uint32_t)KYBER_Q;
    printf("%d\n", res_mine);
    if (res_mine != res_good)
    {
      printf("Barrett reduction error: (%d) %d != %d\n", i, res_mine, res_good);
      return 1;
    }
  }

  return 0;
}
static int test_barrett_reduction(void)
{

  for (int32_t i = -(1 << 15); i < (1 << 15); i++)
  {
    int16_t res_mine = barrett_reduce((int16_t)i);
    int16_t res_good = (int16_t)i % (uint32_t)KYBER_Q;
    if (res_good > (KYBER_Q / 2))
      res_good -= KYBER_Q;
    if (res_mine != res_good)
      return 1;
  }

  return 0;
}
static int test_barrett_reduction_wide(void)
{

  for (uint32_t i = 0; i < 1 << 16; i++)
  {
    uint16_t res_mine = barrett_reduce_wide_u(i);
    uint16_t res_good = i % (uint32_t)WIDE_KYBER_Q;
    printf("%d\n", res_mine);
    if (res_mine != res_good)
    {
      printf("Barrett reduction error: (%d) %d != %d\n", i, res_mine, res_good);
      return 1;
    }
  }

  return 0;
}
static int test_barrett_reduction_wide_s(void)
{
  for (int32_t i = -(1 << 15); i < (1 << 15); i++)
  {
    int16_t res_mine = barrett_reduce_wide_s((int16_t)i);
    printf("%d\n", res_mine);
    int16_t res_good = i % WIDE_KYBER_Q_S;
    if (res_good > (WIDE_KYBER_Q_S / 2))
      res_good -= WIDE_KYBER_Q_S;
    if (res_good < (-WIDE_KYBER_Q_S / 2))
      res_good += WIDE_KYBER_Q_S;
    if (res_mine != res_good)
    {
      printf("Barrett reduction error (%d): %d != %d", i, res_good, res_mine);
      return 1;
    }
  }
  return 0;
}
static int test_montgomery_reduction(void)
{

  for (int32_t i = -KYBER_Q * (1 << 15); i < KYBER_Q * (1 << 15); i++)
  {
    int16_t res_mine = montgomery_reduce(i);
    if (res_mine < 0)
      res_mine += KYBER_Q;
    int16_t res_good = ((i % KYBER_Q) * 169) % KYBER_Q;
    if (res_good < 0)
      res_good += KYBER_Q;
    if ((res_mine != res_good))
    {
      printf("Montgomery reduction error: (%d) %d != %d\n", i, res_mine,
             res_good);
      return 1;
    }
  }

  return 0;
}

static int test_montgomery_reduction_u(void)
{

  for (uint32_t i = 0; i < KYBER_Q * (1 << 16); i++)
  {
    uint16_t res_mine = montgomery_reduce_u(i);
    printf("%hu\n", res_mine);
    if (res_mine >= KYBER_Q)
      res_mine -= KYBER_Q;
    uint16_t res_good = ((i % KYBER_Q) * MINV_U) % KYBER_Q;
    if (res_mine != res_good)
    {
      printf("Montgomery reduction error: (%d) %d != %d\n", i, res_mine,
             res_good);
      return 1;
    }
  }

  return 0;
}

static int test_montgomery_reduction_wide_u(void)
{

  for (uint32_t i = 0; i < WIDE_KYBER_Q * (1 << 16); i++)
  {
    uint16_t res_mine = montgomery_reduce_wide_u(i);
    printf("%d\n", res_mine);
    if (res_mine >= WIDE_KYBER_Q)
      res_mine -= WIDE_KYBER_Q;
    uint16_t res_good = ((i % WIDE_KYBER_Q) * WIDE_MINV_U) % WIDE_KYBER_Q;
    if (res_mine != res_good)
    {
      printf("Montgomery reduction error: (%d) %d != %d\n", i, res_mine,
             res_good);
      return 1;
    }
  }

  return 0;
}
static int test_montgomery_reduction_wide_s(void)
{
  for (int32_t i = -(WIDE_KYBER_Q_S * (1 << 15));
       i < (WIDE_KYBER_Q_S * (1 << 15)); i++)
  {
    int16_t res_mine = montgomery_reduce_wide_s(i);
    if (res_mine < -WIDE_KYBER_Q_S || res_mine > WIDE_KYBER_Q_S)
    {
      printf("ERROR %d, reduces outside the range", i);
      return 1;
    }
    if (res_mine < 0)
      res_mine += WIDE_KYBER_Q_S;
    int16_t res_good = ((i % WIDE_KYBER_Q_S) * WIDE_MINV_S) % WIDE_KYBER_Q_S;
    if (res_good < 0)
      res_good += WIDE_KYBER_Q_S;
    if (res_mine != res_good)
    {
      printf("Montgomery reduction error: (%d) %d != %d\n", i, res_mine,
             res_good);
      return 1;
    }
  }
  return 0;
}
#define popcount16(x) __builtin_popcount((uint16_t)(x))

static int test_mutual_information_u(void)
{
#define HW_SIZE 17
  uint32_t hw_hist[HW_SIZE][HW_SIZE][HW_SIZE][HW_SIZE];
  memset(hw_hist, 0, HW_SIZE * HW_SIZE * HW_SIZE * HW_SIZE * sizeof(uint32_t));

  // invntt_u core
  // t = r[j];
  // r[j] = barrett_reduce_u(t + r[j + len]);
  // r[j + len] = 2 * KYBER_Q + r[j + len] - t;
  // r[j + len] = fqmul_u(zeta, r[j + len]);

  uint16_t t = 0;
  for (uint16_t ai = 0; ai < KYBER_Q; ai++)
    for (uint16_t bi = 0; bi < KYBER_Q; bi++)
    {
      uint16_t a = ai;
      uint16_t b = bi;

      const uint8_t hw_front_a = popcount16(a);
      const uint8_t hw_front_b = popcount16(b);
      assert(hw_front_a <= 16);
      assert(hw_front_b <= 16);

      t = a;
      a = barrett_reduce_u(t + b);
      b = 2 * KYBER_Q + b - t;
      b = montgomery_reduce_u((uint32_t)zetas_u[0] * (uint16_t)b);

      const uint8_t hw_back_a = popcount16(a);
      const uint8_t hw_back_b = popcount16(b);
      assert(hw_back_a <= 16);
      assert(hw_back_b <= 16);

      hw_hist[hw_front_a][hw_front_b][hw_back_a][hw_back_b] += 1;
    }

  float H_AB_given_hws = 0.0;
  for (uint32_t h_a_front_i = 0; h_a_front_i < 17; h_a_front_i += 1)
    for (uint32_t h_b_front_i = 0; h_b_front_i < 17; h_b_front_i += 1)
      for (uint32_t h_a_back_i = 0; h_a_back_i < 17; h_a_back_i += 1)
        for (uint32_t h_b_back_i = 0; h_b_back_i < 17; h_b_back_i += 1)
        {
          uint32_t mass =
              hw_hist[h_a_front_i][h_b_front_i][h_a_back_i][h_b_back_i];
          if (mass == 0)
            continue;
          float ph = (float)mass / (KYBER_Q * KYBER_Q);
          float pab_given_hws = 1.0 / mass;
          H_AB_given_hws -= ph * log2(pab_given_hws);
        }

  float H_AB = -log2(1.0 / (KYBER_Q * KYBER_Q));
  printf("Mutual information (unsigned):\n");
  printf("H(A, B)       = %f\n", H_AB);
  printf("H(A, B | hw(A), hw(B), hw(A'), hw(B')) = %f\n", H_AB_given_hws);
  printf("I(A, B ; hw(A), hw(B), hw(A'), hw(B')) = %f\n",
         H_AB - H_AB_given_hws);
  return 0;
}

static int test_mutual_information_s(void)
{
#define HW_SIZE 17
  uint32_t hw_hist[HW_SIZE][HW_SIZE][HW_SIZE][HW_SIZE];
  memset(hw_hist, 0, HW_SIZE * HW_SIZE * HW_SIZE * HW_SIZE * sizeof(uint32_t));

  // invntt_u core
  // t = r[j];
  // r[j] = barrett_reduce_u(t + r[j + len]);
  // r[j + len] = 2 * KYBER_Q + r[j + len] - t;
  // r[j + len] = fqmul_u(zeta, r[j + len]);

  int16_t t = 0;
  for (int16_t ai = 0; ai < KYBER_Q; ai++)
    for (int16_t bi = 0; bi < KYBER_Q; bi++)
    {
      int16_t a = ai > KYBER_Q / 2 ? ai - KYBER_Q : ai;
      int16_t b = bi > KYBER_Q / 2 ? bi - KYBER_Q : bi;

      const uint8_t hw_front_a = popcount16(a);
      const uint8_t hw_front_b = popcount16(b);
      assert(hw_front_a <= 16);
      assert(hw_front_b <= 16);

      t = a;
      a = barrett_reduce(t + b);
      b = b - t;
      b = montgomery_reduce((int32_t)zetas[0] * (int16_t)b);

      const uint8_t hw_back_a = popcount16(a);
      const uint8_t hw_back_b = popcount16(b);
      assert(hw_back_a <= 16);
      assert(hw_back_b <= 16);

      hw_hist[hw_front_a][hw_front_b][hw_back_a][hw_back_b] += 1;
    }

  float H_AB_given_hws = 0.0;
  for (uint32_t h_a_front_i = 0; h_a_front_i < 17; h_a_front_i += 1)
    for (uint32_t h_b_front_i = 0; h_b_front_i < 17; h_b_front_i += 1)
      for (uint32_t h_a_back_i = 0; h_a_back_i < 17; h_a_back_i += 1)
        for (uint32_t h_b_back_i = 0; h_b_back_i < 17; h_b_back_i += 1)
        {
          uint32_t mass =
              hw_hist[h_a_front_i][h_b_front_i][h_a_back_i][h_b_back_i];
          if (mass == 0)
            continue;
          float ph = (float)mass / (KYBER_Q * KYBER_Q);
          float pab_given_hws = 1.0 / mass;
          H_AB_given_hws -= ph * log2(pab_given_hws);
        }

  float H_AB = -log2(1.0 / (KYBER_Q * KYBER_Q));
  printf("Mutual information (signed):\n");
  printf("H(A, B)       = %f\n", H_AB);
  printf("H(A, B | hw(A), hw(B), hw(A'), hw(B')) = %f\n", H_AB_given_hws);
  printf("I(A, B ; hw(A), hw(B), hw(A'), hw(B')) = %f\n",
         H_AB - H_AB_given_hws);
  return 0;
}

#define SIGMA2 1
#define SIGMA sqrt(SIGMA2)
#define PI 3.14159265

double normal_pdf(double x, double mu);
double centered_normal_pdf(double x);

double normal_pdf(double x, double mu)
{
  return (1.0 / (SIGMA * sqrt(2.0 * PI))) *
         exp(-0.5 * pow((x - mu) / SIGMA, 2.0));
}

double centered_normal_pdf(double x) { return normal_pdf(x, 0.0); }

static int prob_X_eq_x_given_H_eq_h(double h)
{
  uint32_t hw_hist[HW_SIZE];
  memset(hw_hist, 0, HW_SIZE * sizeof(uint32_t));

  uint32_t hw_hist_joint[KYBER_Q][HW_SIZE];
  memset(hw_hist_joint, 0, KYBER_Q * HW_SIZE * sizeof(uint32_t));

  double prob_A_given_H_eq_h[KYBER_Q];
  memset(prob_A_given_H_eq_h, 0, KYBER_Q * sizeof(double));

  for (int ai = 0; ai < DOMAIN_SIZE; ai++)
  {
    uint32_t a = ai;
    const uint8_t hw_a = popcount16(a);
    hw_hist[hw_a] += 1;
    uint32_t real_a = a % KYBER_Q;
    hw_hist_joint[real_a][hw_a] += 1;
  }

  for (int a = 0; a < KYBER_Q; a++)
  {
    uint32_t s = 0;
    for (uint32_t hw = 0; hw < HW_SIZE; hw++)
      s += hw_hist_joint[a][hw];
    assert(s == REDUNDANCY);
  }

  for (int a = 0; a < KYBER_Q; a++)
  {
    double prob_H_eq_h = 0.0;
    double prob_H_eq_h_given_A_eq_a = 0.0;

    for (uint32_t hw = 0; hw < HW_SIZE; hw++)
    {
      const double prob_Z_minus_hw = centered_normal_pdf(h - hw);
      prob_H_eq_h += ((double)hw_hist[hw] / DOMAIN_SIZE) * prob_Z_minus_hw;
      prob_H_eq_h_given_A_eq_a +=
          ((double)hw_hist_joint[a][hw] / REDUNDANCY) * prob_Z_minus_hw;
    }

    prob_A_given_H_eq_h[a] = prob_H_eq_h_given_A_eq_a / (prob_H_eq_h * KYBER_Q);
  }

  printf("[\n");
  for (int a = 0; a < KYBER_Q; a++)
  {
    if (a % 8 == 0)
      printf("    ");
    printf("%f, ", prob_A_given_H_eq_h[a]);
    if (a % 8 == 7)
      printf("\n");
  }
  printf("]\n");
  return 0;
}

int main(void)
{
  (void)test_ntt_u;
  (void)test_ntt;
  (void)test_ntt_wide_u;
  (void)test_ntt_wide_s();
  (void)test_barrett_reduction_u;
  (void)test_montgomery_reduction_u;
  (void)test_mutual_information_u;
  (void)test_mutual_information_s;
  (void)test_mutual_information_wide;
  (void)test_barrett_reduction;
  (void)test_barrett_reduction_u;
  (void)test_barrett_reduction_wide_s;
  // iftest_barrett_reduction_u()) return 1;
  (void)test_montgomery_reduction;
  (void)test_montgomery_reduction_wide_u;
  (void)test_montgomery_reduction_wide_s;
  // if (test_montgomery_reduction_wide()) return 2;
  //  if (test_ntt_u()) return 3;
  //  test_mutual_information_u();
  //  test_mutual_information_s();
  (void)test_barrett_reduction_wide;
  (void)prob_X_eq_x_given_H_eq_h;
  (void)fqmul;

  return 0;
}
