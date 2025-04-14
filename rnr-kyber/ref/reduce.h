#ifndef REDUCE_H
#define REDUCE_H

#include "params.h"
#include <stdint.h>

#define MONT -1044 // 2^16 mod q
#define QINV -3327 // q^-1 mod 2^16

#define MONT_U 2285
#define MINV_U 169
#define QINV_U 3327       // (Q_INV_U * KYBER_Q) mod 2^{16} == -1 mod 2^{16}
#define INV_128 3303      // 128^{-1} mod q
#define MONT_SQ_U 1353    // MONT_U^{2} mod q
#define INTT_CONST_U 1441 // MONT_U^{2} * 128^{-1} mod q

//
// #define WIDE_KYBER_Q (5 * KYBER_Q)
// #define WIDE_MONT_U 15601
// #define WIDE_QINV_U 39987
// #define WIDE_MINV_U 10156

// #define WIDE_INTT_CONST_U 14757
#define WIDE_KYBER_Q (5 * KYBER_Q)
#define WIDE_QINV_U 39987
#define WIDE_MONT_U 15601
#define WIDE_MINV_U 10156

#define WIDE_KYBER_Q9 (9 * KYBER_Q)
#define WIDE_QINV_U9 22215
#define WIDE_MONT_U9 5614
#define WIDE_MINV_U9 10156

#define WIDE_KYBER_Q_S (9 * KYBER_Q)
#define WIDE_QINV_S -22215
#define WIDE_MONT_S 5614
#define WIDE_MINV_S 10156

#define montgomery_reduce KYBER_NAMESPACE(montgomery_reduce)
int16_t montgomery_reduce(int32_t a);

#define barrett_reduce KYBER_NAMESPACE(barrett_reduce)
int16_t barrett_reduce(int16_t a);

#define barrett_reduce_wide_s KYBER_NAMESPACE(barrett_reduce_wide_s)
int16_t barrett_reduce_wide_s(int16_t a);

#define montgomery_reduce_u KYBER_NAMESPACE(montgomery_reduce_u)
uint16_t montgomery_reduce_u(uint32_t a);

#define barrett_reduce_u KYBER_NAMESPACE(barrett_reduce_u)
uint16_t barrett_reduce_u(uint16_t a);

#define montgomery_reduce_wide_u KYBER_NAMESPACE(montgomery_reduce_wide_u)
uint16_t montgomery_reduce_wide_u(uint32_t a);

#define montgomery_reduce_wide_u9 KYBER_NAMESPACE(montgomery_reduce_wide_u9)
uint16_t montgomery_reduce_wide_u9(uint32_t a);

#define montgomery_reduce_wide_s KYBER_NAMESPACE(montgomery_reduce_wide_s)
int16_t montgomery_reduce_wide_s(int32_t a);

#define barrett_reduce_wide_u KYBER_NAMESPACE(barrett_reduce_wide_u)
uint16_t barrett_reduce_wide_u(uint16_t a);

#define barrett_reduce_wide_u9 KYBER_NAMESPACE(barrett_reduce_wide_u9)
uint16_t barrett_reduce_wide_u9(uint16_t a);

#define fqmul KYBER_NAMESPACE(fqmul)
int16_t fqmul(int16_t a, int16_t b);

#define fqmul_u KYBER_NAMESPACE(fqmul_u)
uint16_t fqmul_u(uint16_t a, uint16_t b);

#define fqmul_wide_u KYBER_NAMESPACE(fqmul_wide_u)
uint16_t fqmul_wide_u(uint16_t a, uint16_t b);

#define fqmul_wide_u9 KYBER_NAMESPACE(fqmul_wide_u9)
uint16_t fqmul_wide_u9(uint16_t a, uint16_t b);

#define fqmul_wide_s KYBER_NAMESPACE(fqmul_wide_s)
int16_t fqmul_wide_s(int16_t a, int16_t b);
#endif
