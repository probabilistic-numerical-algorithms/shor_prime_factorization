#ifndef PRIME_FACTOR_C23_H
#define PRIME_FACTOR_C23_H

#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>

#ifndef nullptr
#define nullptr ((void*)0)
#endif

#if defined(HAVE_GMP) && HAVE_GMP
#include <gmp.h>
typedef mpz_t mpe_t;
typedef __mpz_struct mpe_struct_t;
#else
#define MPE_LIMBS 128

typedef struct {
    uint64_t limbs[MPE_LIMBS];
    size_t used;
} mpe_struct_t;

typedef mpe_struct_t mpe_t[1];
#endif

typedef struct {
    mpe_t *factors;
    size_t count;
} factorization_result_t;

// MPE multi-precision core library
void mpe_init(mpe_t x);
void mpe_clear(mpe_t x);
void mpe_set(mpe_t dest, const mpe_t src);
void mpe_set_ui(mpe_t dest, unsigned long val);
void mpe_init_set(mpe_t dest, const mpe_t src);
void mpe_init_set_ui(mpe_t dest, unsigned long val);
int mpe_init_set_str(mpe_t dest, const char *str, int base);
int mpe_cmp(const mpe_t a, const mpe_t b);
int mpe_cmp_ui(const mpe_t a, unsigned long b);
void mpe_add(mpe_t dest, const mpe_t a, const mpe_t b);
void mpe_add_ui(mpe_t dest, const mpe_t a, unsigned long b);
void mpe_sub(mpe_t dest, const mpe_t a, const mpe_t b);
void mpe_sub_ui(mpe_t dest, const mpe_t a, unsigned long b);
void mpe_mul(mpe_t dest, const mpe_t a, const mpe_t b);
void mpe_div_qr(mpe_t q, mpe_t r, const mpe_t n, const mpe_t d);
void mpe_divexact(mpe_t q, const mpe_t n, const mpe_t d);
void mpe_mod(mpe_t r, const mpe_t n, const mpe_t d);
void mpe_powm(mpe_t r, const mpe_t base, const mpe_t exp, const mpe_t mod);
void mpe_urandomm(mpe_t r, const mpe_t range);
int mpe_printf(const char *format, ...);

// Primality and factoring
void gcd_c23(mpe_t result, const mpe_t A, const mpe_t B);
bool is_prime_c23(const mpe_t N, int k);
bool shor_factorize_c23(mpe_t factor, const mpe_t N);
factorization_result_t factorize_c23(const mpe_t N);
void free_factorization_result(factorization_result_t *result);

#endif // PRIME_FACTOR_C23_H
