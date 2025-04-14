#include "../reduce.h"
#include "../ntt.h"
#include "../randombytes.h"
#include <stdio.h>
#include <string.h>
#define NTESTS 10000
static int test_ntt(void)
{
#define SIZE 256
    int16_t input[SIZE];
    int16_t output[SIZE];

    for (uint32_t test = 0; test < NTESTS; test++)
    {
        randombytes((uint8_t *)input, SIZE * sizeof(int16_t));
        for (uint32_t i = 0; i < SIZE; i++)
        {
            input[i] %= KYBER_Q;
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
            output[i] = ((int32_t)output[i] * (int32_t)MINV_U) % KYBER_Q;
            if (output[i] > KYBER_Q / 2)
                output[i] -= KYBER_Q;
        }

        if (memcmp(input, output, SIZE * sizeof(int16_t)) != 0)
        {
            printf("ERROR test %d: invntt(ntt(x)) != x\n", test);
            for (int i = 0; i < SIZE; i++)
                printf("%d %d\n", input[i], output[i]);

            return 1;
        }
    }
    return 0;
}
static int test_invntt_full(void)
{
#define SIZE 256
    int16_t arr[32] = {996, -390, 73, 827, 1355, 797, -805, -859, -631, -483, 1453, -313, -1314, 1416, -1421, -483, 686, 823, 308, -1323, 1422, 1039, 1250, 1040, 93, 568, 1091, 1307, 536, -687, -733, 549};
    int16_t input[SIZE];
    int16_t output[SIZE];
    int16_t full[SIZE][9] = {0};
    for (uint32_t test = 0; test < 1; test++)
    {
        randombytes((uint8_t *)input, SIZE * sizeof(int16_t));
        for (uint32_t i = 0; i < SIZE; i++)
        {
            input[i] = 0;
            input[i] %= KYBER_Q;
            if (input[i] > KYBER_Q / 2)
                input[i] -= KYBER_Q;
            if (input[i] < -KYBER_Q / 2)
                input[i] += KYBER_Q;
        }
        for (int i = 0; i < 32; i++)
            input[i] = arr[i];
        invntt_full(input, full);

        // for (int i = 0; i < SIZE; i++)
        // {
        //     for (int j = 0; j < 9; j++)
        //     {
        //         printf("%4d ", full[i][j]);
        //     }
        //     printf("\n");
        // }
    }
    return 0;
}
static int test_ntt_u_equiv(void)
{
#define SIZE 256
    uint16_t input[SIZE];
    int16_t ntt_output[SIZE];
    uint16_t nttu_output[SIZE];

    for (uint32_t test = 0; test < NTESTS; test++)
    {
        // this is not uniform, but I do not care for testing
        randombytes((uint8_t *)input, SIZE * sizeof(uint16_t));
        for (uint32_t i = 0; i < SIZE; i++)
        {
            input[i] %= KYBER_Q;
            ntt_output[i] = input[i];
            nttu_output[i] = input[i];
            if (ntt_output[i] > KYBER_Q / 2)
                ntt_output[i] -= KYBER_Q;
            if (ntt_output[i] < -KYBER_Q / 2)
                ntt_output[i] += KYBER_Q;
        }
        ntt(ntt_output);
        ntt_u(nttu_output);

        invntt(ntt_output);
        invntt_u(nttu_output);

        for (uint32_t i = 0; i < SIZE; i++)
        {

            if (ntt_output[i] < 0)
            {
                ntt_output[i] += KYBER_Q;
            }
            if (ntt_output[i] != barrett_reduce_u(nttu_output[i]))
            {
                printf("ERROR test %d: ntt(x) != ntt_u(x) %d %d\n", test, ntt_output[i], nttu_output[i]);

                return 1;
            }
        }
    }
}

static int test_ntt_wide_s_equiv(void)
{
#define SIZE 256
    uint16_t input[SIZE];
    int16_t ntt_output[SIZE];
    int16_t ntt_wide_output[SIZE];

    for (uint32_t test = 0; test < NTESTS; test++)
    {
        // this is not uniform, but I do not care for testing
        randombytes((uint8_t *)input, SIZE * sizeof(uint16_t));
        for (uint32_t i = 0; i < SIZE; i++)
        {
            input[i] %= WIDE_KYBER_Q_S;
            ntt_wide_output[i] = input[i];

            if (ntt_wide_output[i] > WIDE_KYBER_Q_S / 2)
                ntt_wide_output[i] -= WIDE_KYBER_Q_S;
            if (ntt_wide_output[i] < -WIDE_KYBER_Q_S / 2)
                ntt_wide_output[i] += WIDE_KYBER_Q_S;
            ntt_output[i] = ntt_wide_output[i] % KYBER_Q;
        }
        ntt(ntt_output);
        ntt_wide_s(ntt_wide_output);
        invntt(ntt_output);
        invntt_wide_s(ntt_wide_output);
        for (uint32_t i = 0; i < SIZE; i++)
        {

            if (barrett_reduce(ntt_output[i]) != barrett_reduce(ntt_wide_output[i] % KYBER_Q))
            {
                printf("ERROR test %d: ntt(x) != ntt_u(x) %d %d %d", test, ntt_output[i], ntt_wide_output[i], ntt_wide_output[i] % KYBER_Q);

                return 1;
            }
        }
    }
}

static int test_ntt_wide_u_equiv(void)
{
#define SIZE 256
    uint16_t input[SIZE];
    uint16_t ntt_output[SIZE];
    uint16_t ntt_wide_output[SIZE];

    for (uint32_t test = 0; test < NTESTS; test++)
    {
        // this is not uniform, but I do not care for testing
        randombytes((uint8_t *)input, SIZE * sizeof(uint16_t));
        for (uint32_t i = 0; i < SIZE; i++)
        {
            input[i] %= WIDE_KYBER_Q;
            ntt_wide_output[i] = input[i];
            ntt_output[i] = input[i] % KYBER_Q;
        }
        ntt_u(ntt_output);
        ntt_wide_u(ntt_wide_output);
        invntt_u(ntt_output);
        invntt_wide_u(ntt_wide_output);
        for (uint32_t i = 0; i < SIZE; i++)
        {

            if (barrett_reduce(ntt_output[i]) != barrett_reduce(ntt_wide_output[i] % KYBER_Q))
            {
                printf("ERROR test %d: ntt(x) != ntt_u(x) %d %d %d", test, ntt_output[i], ntt_wide_output[i], ntt_wide_output[i] % KYBER_Q);

                return 1;
            }
        }
    }
}
static int test_ntt_wide_u9_equiv(void)
{
#define SIZE 256
    uint16_t input[SIZE];
    uint16_t ntt_output[SIZE];
    uint16_t ntt_wide_output[SIZE];

    for (uint32_t test = 0; test < NTESTS; test++)
    {
        // this is not uniform, but I do not care for testing
        randombytes((uint8_t *)input, SIZE * sizeof(uint16_t));
        for (uint32_t i = 0; i < SIZE; i++)
        {
            input[i] %= WIDE_KYBER_Q9;
            ntt_wide_output[i] = input[i];
            ntt_output[i] = input[i] % KYBER_Q;
        }
        ntt_u(ntt_output);
        ntt_wide_u9(ntt_wide_output);
        invntt(ntt_output);
        invntt_wide_u9(ntt_wide_output);
        for (uint32_t i = 0; i < SIZE; i++)
        {

            if (barrett_reduce(ntt_output[i]) != barrett_reduce(ntt_wide_output[i] % KYBER_Q))
            {
                printf("ERROR test %d: ntt(x) != ntt_u(x) %d %d %d", test, ntt_output[i], ntt_wide_output[i], ntt_wide_output[i] % KYBER_Q);

                return 1;
            }
        }
    }
}
static int test_ntt_wide_s(void)
{
#define SIZE 256
    int16_t input[SIZE];
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

        memcpy(output, input, SIZE * sizeof(int16_t));
        ntt_wide_s(output);
        invntt_wide_s(output);
        for (uint32_t i = 0; i < SIZE; i++)
        {
            output[i] = (output[i] <= 0) ? WIDE_KYBER_Q_S + output[i] : output[i];
            output[i] = ((uint32_t)output[i] * (uint32_t)WIDE_MINV_S) % KYBER_Q;
            input[i] %= KYBER_Q;

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
static int test_ntt_u(void)
{
#define SIZE 256
    uint16_t input[SIZE];
    uint16_t output[SIZE];
    for (uint32_t test = 0; test < NTESTS; test++)
    {
        randombytes((uint8_t *)input, SIZE * sizeof(uint16_t));
        for (uint32_t i = 0; i < SIZE; i++)
            input[i] %= KYBER_Q;
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
}
static int test_ntt_wide_u(void)
{
#define SIZE 256
    uint16_t input[SIZE];
    uint16_t output[SIZE];
    for (uint32_t test = 0; test < NTESTS; test++)
    {
        randombytes((uint8_t *)input, SIZE * sizeof(uint16_t));
        for (uint32_t i = 0; i < SIZE; i++)
            input[i] %= WIDE_KYBER_Q;
        memcpy(output, input, SIZE * sizeof(uint16_t));
        ntt_u(output);
        invntt_u(output);

        for (uint32_t i = 0; i < SIZE; i++)
        {
            output[i] = ((uint32_t)output[i] * (uint32_t)WIDE_MINV_U) % KYBER_Q;
            input[i] %= KYBER_Q;
        }
        if (memcmp(input, output, SIZE * sizeof(uint16_t)) != 0)
        {
            printf("ERROR test %d: invntt_u(ntt_u(x)) != x\n", test);
            for (int i = 0; i < SIZE; i++)
                printf("%4hu %4hu\n", input[i], output[i]);

            return 1;
        }
    }
}
static int test_zeta(void)
{
    int16_t input[SIZE];
    uint32_t hist[128] = {0};
    for (int t = 0; t < 10000; t++)
    {
#define TEST_SIZE 1
        randombytes((uint8_t *)input, SIZE * sizeof(int16_t));
        for (uint32_t i = 0; i < SIZE; i++)
        {
            input[i] %= WIDE_KYBER_Q_S;
            if (input[i] > WIDE_KYBER_Q_S / 2)
                input[i] -= WIDE_KYBER_Q_S;
            if (input[i] < -WIDE_KYBER_Q_S / 2)
                input[i] += WIDE_KYBER_Q_S;
        }
        for (uint32_t i = 0; i < TEST_SIZE; i++)
            input[i] %= WIDE_KYBER_Q_S;
        for (int i = 0; i < 128; i++)
        {
            int16_t x = fqmul_wide_s(zetas[i], input[0]);
            int16_t y = fqmul_wide_s(zetas[i], input[0] % KYBER_Q);
            if (x == y)
            {
                hist[i]++;
            }
        }
    }

    for (int i = 0; i < 128; i++)
    {
        printf("%d %d\n", i, hist[i]);
    }
    return 0;
}
int main(void)
{
    // if (test_ntt())
    //     return 1;
    // if (test_invntt_full())
    //     return 1;
    // if (test_ntt_wide_s())
    //    return 1;
    // if (test_ntt_u())
    //     return 1;
    // if (test_ntt_wide_u())
    //     return 1;
    // printf("Test pass\n");
    // test_zeta();
    if (test_ntt_wide_u9_equiv())
        return 1;
    return 0;
}