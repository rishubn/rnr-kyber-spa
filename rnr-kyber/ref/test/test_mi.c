#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "../reduce.h"
#define assert(x)            \
    {                        \
        if (!(x))            \
        {                    \
            *((int *)0) = 0; \
        }                    \
    }
#define popcount16(x) __builtin_popcount((uint16_t)(x))

#define KYBER_Q 3329
#define HW_SIZE 17
#define REDUNDANCY_S 9
#define DOMAIN_SIZE_S (REDUNDANCY_S * KYBER_Q)
#define REDUNDANCY_U 5
#define DOMAIN_SIZE_U (REDUNDANCY_U * KYBER_Q)

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
static int test_mutual_information_wide_u(void)
{
    uint32_t hw_hist[HW_SIZE];
    memset(hw_hist, 0, HW_SIZE * sizeof(uint32_t));

    uint32_t hw_hist_cond[KYBER_Q][HW_SIZE];
    memset(hw_hist_cond, 0, KYBER_Q * HW_SIZE * sizeof(uint32_t));

    for (int ai = 0; ai < DOMAIN_SIZE_U; ai++)
    {
        uint32_t a = ai;
        const uint8_t hw_a = popcount16(a);
        hw_hist[hw_a] += 1;
        uint32_t real_a = a % KYBER_Q;
        hw_hist_cond[real_a][hw_a] += 1;
    }

    float H_A = -log2(1.0 / (KYBER_Q));
    float H_A_given_hw = 0.0;

    for (uint32_t h_a = 0; h_a < 17; h_a += 1)
    {
        if (hw_hist[h_a] == 0)
            continue;
        float p_hw_a = hw_hist[h_a] * 1.0 / DOMAIN_SIZE_U;
        float h_a_cond_hw = 0;
        for (uint32_t ai = 0; ai < KYBER_Q; ai++)
        {
            if (hw_hist_cond[ai][h_a] == 0)
                continue;
            float p = (hw_hist_cond[ai][h_a] * 1.0 / hw_hist[h_a]);
            assert(p > 0);
            h_a_cond_hw += p * log2(p);
        }
        H_A_given_hw += -p_hw_a * h_a_cond_hw;
    }

    printf("Mutual information (wide):\n");
    printf("H(A)       = %f\n", H_A);
    printf("H(A | hw(wide A)) = %f\n", H_A_given_hw);
    printf("I(A ; hw(wide A)) = %f\n", H_A - H_A_given_hw);
    return 0;
}
static int test_mutual_information_wide_s(void)
{
    uint32_t hw_hist[HW_SIZE];
    memset(hw_hist, 0, HW_SIZE * sizeof(uint32_t));

    uint32_t hw_hist_cond[KYBER_Q][HW_SIZE];
    memset(hw_hist_cond, 0, KYBER_Q * HW_SIZE * sizeof(uint32_t));

    for (int ai = 0; ai < DOMAIN_SIZE_S; ai++)
    {
        uint32_t a = ai > DOMAIN_SIZE_S / 2 ? ai - DOMAIN_SIZE_S : ai;
        const uint8_t hw_a = popcount16(a);
        hw_hist[hw_a] += 1;
        uint32_t real_a = a % KYBER_Q;
        hw_hist_cond[real_a][hw_a] += 1;
    }

    float H_A = -log2(1.0 / (KYBER_Q));
    float H_A_given_hw = 0.0;

    for (uint32_t h_a = 0; h_a < 17; h_a += 1)
    {
        if (hw_hist[h_a] == 0)
            continue;
        float p_hw_a = hw_hist[h_a] * 1.0 / DOMAIN_SIZE_S;
        float h_a_cond_hw = 0;
        for (uint32_t ai = 0; ai < KYBER_Q; ai++)
        {
            if (hw_hist_cond[ai][h_a] == 0)
                continue;
            float p = (hw_hist_cond[ai][h_a] * 1.0 / hw_hist[h_a]);
            assert(p > 0);
            h_a_cond_hw += p * log2(p);
        }
        H_A_given_hw += -p_hw_a * h_a_cond_hw;
    }

    printf("Mutual information (wide signed):\n");
    printf("H(A)       = %f\n", H_A);
    printf("H(A | hw(wide A)) = %f\n", H_A_given_hw);
    printf("I(A ; hw(wide A)) = %f\n", H_A - H_A_given_hw);
    return 0;
}
static int test_mutual_information_bff_u(void)
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

static int test_mutual_information_bff_s(void)
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

static int test_mutual_information_bff_wide_s(void)
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
    for (int16_t ai = 0; ai < DOMAIN_SIZE_S; ai++)
        for (int16_t bi = 0; bi < DOMAIN_SIZE_S; bi++)
        {
            int16_t a = ai > DOMAIN_SIZE_S / 2 ? ai - DOMAIN_SIZE_S : ai;
            int16_t b = bi > DOMAIN_SIZE_S / 2 ? bi - DOMAIN_SIZE_S : bi;

            const uint8_t hw_front_a = popcount16(a);
            const uint8_t hw_front_b = popcount16(b);
            assert(hw_front_a <= 16);
            assert(hw_front_b <= 16);

            t = a;
            a = barrett_reduce_wide_s(t + b);
            b = b - t;
            b = montgomery_reduce_wide_s((int32_t)zetas[0] * (int16_t)b);

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
                    float ph = (float)mass / (DOMAIN_SIZE_S * DOMAIN_SIZE_S);
                    float pab_given_hws = 1.0 / mass;
                    H_AB_given_hws -= ph * log2(pab_given_hws);
                }

    float H_AB = -log2(1.0 / (DOMAIN_SIZE_S * DOMAIN_SIZE_S));
    printf("Mutual information (wide signed):\n");
    printf("H(A, B)       = %f\n", H_AB);
    printf("H(A, B | hw(A), hw(B), hw(A'), hw(B')) = %f\n", H_AB_given_hws);
    printf("I(A, B ; hw(A), hw(B), hw(A'), hw(B')) = %f\n",
           H_AB - H_AB_given_hws);
    return 0;
}
static int test_mutual_information_bff_wide_u(void)
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
    for (uint16_t ai = 0; ai < DOMAIN_SIZE_U; ai++)
        for (uint16_t bi = 0; bi < DOMAIN_SIZE_U; bi++)
        {
            uint16_t a = ai;
            uint16_t b = bi;

            const uint8_t hw_front_a = popcount16(a);
            const uint8_t hw_front_b = popcount16(b);
            assert(hw_front_a <= 16);
            assert(hw_front_b <= 16);

            t = a;
            a = barrett_reduce_wide_u(t + b);
            b = 2 * DOMAIN_SIZE_U + b - t;
            b = montgomery_reduce_wide_u((uint32_t)zetas_u[0] * (uint16_t)b);

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
                    float ph = (float)mass / (DOMAIN_SIZE_U * DOMAIN_SIZE_U);
                    float pab_given_hws = 1.0 / mass;
                    H_AB_given_hws -= ph * log2(pab_given_hws);
                }

    float H_AB = -log2(1.0 / (DOMAIN_SIZE_U * DOMAIN_SIZE_U));
    printf("Mutual information (wide unsigned):\n");
    printf("H(A, B)       = %f\n", H_AB);
    printf("H(A, B | hw(A), hw(B), hw(A'), hw(B')) = %f\n", H_AB_given_hws);
    printf("I(A, B ; hw(A), hw(B), hw(A'), hw(B')) = %f\n",
           H_AB - H_AB_given_hws);
    return 0;
}

int main(void)
{
    test_mutual_information_bff_u();
    test_mutual_information_bff_s();
    test_mutual_information_bff_wide_s();
    test_mutual_information_bff_wide_u();
    // test_mutual_information_wide_u();
    // test_mutual_information_wide_s();
    return 0;
}