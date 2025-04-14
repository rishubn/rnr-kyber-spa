#include "../reduce.h"
#include <stdio.h>
static int test_barrett_reduction(void)
{
    for (int32_t i = -(1 << 15); i < (1 << 15); i++)
    {
        int16_t res_mine = barrett_reduce((int16_t)i);
        int16_t res_good = (int16_t)i % (int32_t)KYBER_Q;

        if (res_good > (KYBER_Q / 2))
            res_good -= KYBER_Q;
        if (res_good < (-KYBER_Q / 2))
            res_good += KYBER_Q;
        if (res_mine != res_good)
        {
            printf("FAIL %i: %d != %d\n", i, res_mine, res_good);
            return 1;
        }
    }
    return 0;
}
static int test_barrett_reduction_wide_u9(void)
{

    for (uint32_t i = 0; i < 1 << 16; i++)
    {
        uint16_t res_mine = barrett_reduce_wide_u9(i);
        uint16_t res_good = i % (uint32_t)WIDE_KYBER_Q9;
        if (res_mine != res_good)
        {
            printf("Barrett reduction error: (%d) %d != %d\n", i, res_mine, res_good);
            return 1;
        }
    }

    return 0;
}
static int test_montgomery_reduction(void)
{

    int min = INT32_MAX;
    int max = INT32_MIN;
    for (int32_t i = -KYBER_Q * (1 << 15); i < KYBER_Q * (1 << 15); i++)
    {
        int16_t res_mine = montgomery_reduce(i);
        if (max < res_mine)
            max = res_mine;
        if (min > res_mine)
            min = res_mine;
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

static int test_montgomery_reduction_u9(void)
{

    for (uint32_t i = 0; i < WIDE_KYBER_Q9 * (1 << 16); i++)
    {
        uint16_t res_mine = montgomery_reduce_wide_u9(i);
        if (res_mine >= WIDE_KYBER_Q9)
            res_mine -= WIDE_KYBER_Q9;
        uint16_t res_good = ((i % WIDE_KYBER_Q9) * WIDE_MINV_U9) % WIDE_KYBER_Q9;

        if ((res_mine != res_good))
        {
            printf("Montgomery reduction error: (%d) %d != %d\n", i, res_mine,
                   res_good);
            return 1;
        }
    }
    return 0;
}
// static int test_montgomery_reduction(void)
// {

//     int min = INT32_MAX;
//     int max = INT32_MIN;
//     for (int32_t i = -KYBER_Q * (1 << 15); i < KYBER_Q * (1 << 15); i++)
//     {
//         int16_t res_mine = montgomery_reduce_wide_u(i);
//         if (max < res_mine)
//             max = res_mine;
//         if (min > res_mine)
//             min = res_mine;
//         if (res_mine < 0)
//             res_mine += KYBER_Q;
//         int16_t res_good = ((i % KYBER_Q) * 169) % KYBER_Q;
//         if (res_good < 0)
//             res_good += KYBER_Q;
//         if ((res_mine != res_good))
//         {
//             printf("Montgomery reduction error: (%d) %d != %d\n", i, res_mine,
//                    res_good);
//             return 1;
//         }
//     }
//     printf("%d %d", min, max);
//     return 0;
// }

int main(void)
{
    test_barrett_reduction();
    test_montgomery_reduction();
    test_barrett_reduction_wide_u9();
    test_montgomery_reduction_u9();
    return 0;
}