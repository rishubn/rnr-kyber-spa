#ifndef NTT_H
#define NTT_H

#include "params.h"
#include <stdint.h>

#define zetas KYBER_NAMESPACE(zetas)
extern const int16_t zetas[128];

#define ntt KYBER_NAMESPACE(ntt)
void ntt(int16_t poly[256]);

#define invntt KYBER_NAMESPACE(invntt)
void invntt(int16_t poly[256]);

#define invntt_full KYBER_NAMESPACE(invntt_full)
void invntt_full(int16_t poly[256], int16_t full[256][9]);

#define basemul KYBER_NAMESPACE(basemul)
void basemul(int16_t r[2], const int16_t a[2], const int16_t b[2],
             int16_t zeta);

#define zetas_u KYBER_NAMESPACE(zetas_u)
extern const uint16_t zetas_u[128];

#define zetas_wide_u KYBER_NAMESPACE(zetas_wide_u)
extern const uint16_t zetas_wide_u[128];
// #define init_ntt_u KYBER_NAMESPACE(init_ntt_u)
// void init_ntt_u(void);

// #define init_ntt_wide_u KYBER_NAMESPACE(init_ntt_wide_u)
// void init_ntt_wide_u(void);

#define ntt_u KYBER_NAMESPACE(ntt_u)
void ntt_u(uint16_t poly[256]);

#define invntt_u KYBER_NAMESPACE(invntt_u)
void invntt_u(uint16_t poly[256]);

#define invntt_u_full KYBER_NAMESPACE(invntt_u_full)
void invntt_u_full(uint16_t poly[256], uint16_t full[256][9]);

#define ntt_wide_u KYBER_NAMESPACE(ntt_wide_u)
void ntt_wide_u(uint16_t poly[256]);

#define invntt_wide_u KYBER_NAMESPACE(invntt_wide_u)
void invntt_wide_u(uint16_t poly[256]);

#define ntt_wide_u9 KYBER_NAMESPACE(ntt_wide_u9)
void ntt_wide_u9(uint16_t poly[256]);

#define invntt_wide_u9 KYBER_NAMESPACE(invntt_wide_u9)
void invntt_wide_u9(uint16_t poly[256]);

#define invntt_wide_u_full KYBER_NAMESPACE(invntt_wide_u_full)
void invntt_wide_u_full(uint16_t poly[256], uint16_t full[256][9]);

#define ntt_wide_s KYBER_NAMESPACE(ntt_wide_s)
void ntt_wide_s(int16_t poly[256]);

#define invntt_wide_s KYBER_NAMESPACE(invntt_wide_s)
void invntt_wide_s(int16_t poly[256]);

#define invntt_wide_s_full KYBER_NAMESPACE(invntt_wide_s_full)
void invntt_wide_s_full(int16_t poly[256], int16_t full[256][9]);

#define basemul_u KYBER_NAMESPACE(basemul_u)
void basemul_u(uint16_t r[2], const uint16_t a[2], const uint16_t b[2],
               uint16_t zeta);

#endif
