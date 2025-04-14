#include "reduce.h"
#include "params.h"
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#define assert(x)      \
  {                    \
    if (!(x))          \
    {                  \
      *((int *)0) = 0; \
    }                  \
  }

/*************************************************
 * Name:        montgomery_reduce
 *
 * Description: Montgomery reduction; given a 32-bit integer a, computes
 *              16-bit integer congruent to a * R^-1 mod q, where R=2^16
 *
 * Arguments:   - int32_t a: input integer to be reduced;
 *                           has to be in {-q2^15,...,q2^15-1}
 *
 * Returns:     integer in {-q+1,...,q-1} congruent to a * R^-1 modulo q.
 **************************************************/
int16_t montgomery_reduce(int32_t a)
{
  int16_t t;

  t = (int16_t)a * QINV;
  t = (a - (int32_t)t * KYBER_Q) >> 16;
  return t;
}
int16_t montgomery_reduce_wide_s(int32_t a)
{
  int16_t t;

  t = (int16_t)a * WIDE_QINV_S;
  t = (a - (int32_t)t * WIDE_KYBER_Q_S) >> 16;
  return t;
}
/// @brief Performs montgomery reduction of $a \in [0, 2^{16} \cdot q - 1]$.
/// The result is in the range $[0, 2q-1]$ congruent to $a \cdot (2^{16})^{-1}
/// mod q$.
/// @see
/// https://en.wikipedia.org/wiki/Montgomery_modular_multiplication#The_REDC_algorithm
/// @param a Number to be reduced
/// @return Reduced number
uint16_t montgomery_reduce_u(uint32_t a)
{
  assert(a < (KYBER_Q << 16));
  // Reduction parameters are:
  // T = a
  // N = q (i.e., KYBER_Q)
  // N' = q^{-1} (i.e., QINV_U)
  // R = 2^16

  // (uint16_t)a = a      mod 2^16
  // QINV_U      = q^{-1} mod 2^16
  // uint16_t m  = (a mod 2^16)(q^{-1} mod 2^16) mod 2^16
  //             = a*q^{-1} mod 2^16
  const uint16_t m = (uint16_t)((uint16_t)a * QINV_U);
  // uint16_t t  = (a + m * q) / 2^16
  const uint16_t t = (a + ((uint32_t)m * KYBER_Q)) >> 16;
  assert(t < 2 * KYBER_Q);
  return t;
}

uint16_t montgomery_reduce_wide_u(uint32_t a)
{
  assert(a < (WIDE_KYBER_Q << 16));
  // Reduction parameters are:
  // T = a
  // N = q (i.e., WIDE_KYBER_Q)
  // N' = q^{-1} (i.e., QINV_U)
  // R = 2^16

  // (uint16_t)a = a      mod 2^16
  // QINV_U      = q^{-1} mod 2^16
  // uint16_t m  = (a mod 2^16)(q^{-1} mod 2^16) mod 2^16
  //             = a*q^{-1} mod 2^16
  const uint16_t m = (uint16_t)((uint16_t)a * WIDE_QINV_U);
  // uint16_t t  = (a + m * q) / 2^16
  const uint16_t t = (a + ((uint32_t)m * WIDE_KYBER_Q)) >> 16;
  // assert(t < 2 * WIDE_KYBER_Q);
  return t;
}
uint16_t montgomery_reduce_wide_u9(uint32_t a)
{
  assert(a < (WIDE_KYBER_Q9 << 16));
  // Reduction parameters are:
  // T = a
  // N = q (i.e., WIDE_KYBER_Q)
  // N' = q^{-1} (i.e., QINV_U)
  // R = 2^16

  // (uint16_t)a = a      mod 2^16
  // QINV_U      = q^{-1} mod 2^16
  // uint16_t m  = (a mod 2^16)(q^{-1} mod 2^16) mod 2^16
  //             = a*q^{-1} mod 2^16
  const uint16_t m = (uint16_t)((uint16_t)a * WIDE_QINV_U9);
  // uint16_t t  = (a + m * q) / 2^16
  const uint16_t t = (a + ((uint32_t)m * WIDE_KYBER_Q9)) >> 16;
  // assert(t < 2 * WIDE_KYBER_Q);
  return t;
}
/*************************************************
 * Name:        barrett_reduce
 *
 * Description: Barrett reduction; given a 16-bit integer a, computes
 *              centered representative congruent to a mod q in
 *{-(q-1)/2,...,(q-1)/2}
 *
 * Arguments:   - int16_t a: input integer to be reduced
 *
 * Returns:     integer in {-(q-1)/2,...,(q-1)/2} congruent to a modulo q.
 **************************************************/
int16_t barrett_reduce(int16_t a)
{
  int16_t t;
  const int16_t v = ((1 << 26) + KYBER_Q / 2) / KYBER_Q;

  t = ((int32_t)v * a + (1 << 25)) >> 26;
  t *= KYBER_Q;
  return a - t;
}

#define WIDE_BARRETT_K_S 29
int16_t barrett_reduce_wide_s(int16_t a)
{
  int16_t t;
  const int16_t v =
      ((1 << WIDE_BARRETT_K_S) + WIDE_KYBER_Q_S / 2) / WIDE_KYBER_Q_S;
  t = ((int32_t)v * a + (1 << 28)) >> WIDE_BARRETT_K_S;
  t *= WIDE_KYBER_Q_S;
  return a - t;
}

/// @brief Performs Barrett reduction of $a \in [0, 2^{16} - 1]$.
/// The result is in the range $[0, q-1]$ congruent to $a mod q$.
/// @see
/// https://en.wikipedia.org/wiki/Barrett_reduction#Single-word_Barrett_reduction
/// @param a Number to be reduced
/// @return Reduced number
#define BARRETT_K 26
const uint16_t BARRETT_CONST = ((1ul << BARRETT_K) + KYBER_Q / 2) / KYBER_Q;
uint16_t barrett_reduce_u(uint16_t a)
{
  const uint16_t t = ((uint32_t)a * BARRETT_CONST) >> BARRETT_K;
  const uint16_t r = a - (t * KYBER_Q);
  assert(r < KYBER_Q);
  return r;
}

#define WIDE_BARRETT_K 31
const uint64_t WIDE_BARRETT_CONST =
    ((1ul << WIDE_BARRETT_K) + (WIDE_KYBER_Q / 2)) / WIDE_KYBER_Q;
uint16_t barrett_reduce_wide_u(uint16_t a)
{
  const uint16_t t = ((uint64_t)a * WIDE_BARRETT_CONST) >> WIDE_BARRETT_K;
  const uint16_t r = a - (t * WIDE_KYBER_Q);
  // assert(r < WIDE_KYBER_Q);
  return r;
}
#define WIDE_BARRETT_K9 29
const uint64_t WIDE_BARRETT_CONST9 =
    ((1ul << WIDE_BARRETT_K9) + (WIDE_KYBER_Q9 / 2)) / WIDE_KYBER_Q9;
uint16_t barrett_reduce_wide_u9(uint16_t a)
{
  const uint16_t t = ((uint64_t)a * WIDE_BARRETT_CONST9) >> WIDE_BARRETT_K9;
  const uint16_t r = a - (t * WIDE_KYBER_Q9);
  // assert(r < WIDE_KYBER_Q);
  return r;
}
int16_t fqmul(int16_t a, int16_t b)
{
  return montgomery_reduce((int32_t)a * b);
}

uint16_t fqmul_u(uint16_t a, uint16_t b)
{
  return montgomery_reduce_u(((uint32_t)a) * b);
}

uint16_t fqmul_wide_u(uint16_t a, uint16_t b)
{
  return montgomery_reduce_wide_u((uint32_t)a * b);
}
uint16_t fqmul_wide_u9(uint16_t a, uint16_t b)
{
  return montgomery_reduce_wide_u9((uint32_t)a * b);
}
int16_t fqmul_wide_s(int16_t a, int16_t b)
{
  return montgomery_reduce_wide_s(((int32_t)a) * b);
}
