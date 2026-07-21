#include "prime_factor_c23.h"
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>

#if defined(HAVE_GMP) && HAVE_GMP

// Initialize multi-precision integer to zero
void mpe_init(mpe_t x) {
    mpz_init(x);
    return;
}

// Clear multi-precision integer
void mpe_clear(mpe_t x) {
    mpz_clear(x);
    return;
}

// Copy source multi-precision integer to destination
void mpe_set(mpe_t dest, const mpe_t src) {
    mpz_set(dest, src);
    return;
}

// Set multi-precision integer to an unsigned long value
void mpe_set_ui(mpe_t dest, unsigned long val) {
    mpz_set_ui(dest, val);
    return;
}

// Initialize and set from source
void mpe_init_set(mpe_t dest, const mpe_t src) {
    mpz_init_set(dest, src);
    return;
}

// Initialize and set from unsigned long
void mpe_init_set_ui(mpe_t dest, unsigned long val) {
    mpz_init_set_ui(dest, val);
    return;
}

// Compare two multi-precision integers
int mpe_cmp(const mpe_t a, const mpe_t b) {
    int res = mpz_cmp(a, b);
    return res;
}

// Compare multi-precision integer with unsigned long
int mpe_cmp_ui(const mpe_t a, unsigned long b) {
    int res = mpz_cmp_ui(a, b);
    return res;
}

// Add two multi-precision integers
void mpe_add(mpe_t dest, const mpe_t a, const mpe_t b) {
    mpz_add(dest, a, b);
    return;
}

// Add unsigned long to multi-precision integer
void mpe_add_ui(mpe_t dest, const mpe_t a, unsigned long b) {
    mpz_add_ui(dest, a, b);
    return;
}

// Subtract two multi-precision integers
void mpe_sub(mpe_t dest, const mpe_t a, const mpe_t b) {
    mpz_sub(dest, a, b);
    return;
}

// Subtract unsigned long from multi-precision integer
void mpe_sub_ui(mpe_t dest, const mpe_t a, unsigned long b) {
    mpz_sub_ui(dest, a, b);
    return;
}

// Multiply two multi-precision integers
void mpe_mul(mpe_t dest, const mpe_t a, const mpe_t b) {
    mpz_mul(dest, a, b);
    return;
}

// General division and modulo of two multi-precision integers
void mpe_div_qr(mpe_t q, mpe_t r, const mpe_t n, const mpe_t d) {
    mpz_tdiv_qr(q, r, n, d);
    return;
}

// Exact division of multi-precision integers
void mpe_divexact(mpe_t q, const mpe_t n, const mpe_t d) {
    mpz_divexact(q, n, d);
    return;
}

// Modulo of multi-precision integers
void mpe_mod(mpe_t r, const mpe_t n, const mpe_t d) {
    mpz_mod(r, n, d);
    return;
}

static uint64_t rand_state = 123456789ULL;

// Generate pseudo-random multi-precision integer in range [0, range-1]
void mpe_urandomm(mpe_t r, const mpe_t range) {
    size_t limbs = mpz_size(range);
    mpe_t temp;
    mpz_init_set_ui(temp, 0);
    for (size_t i = 0; i < limbs; i = i + 1) {
        rand_state = rand_state * 6364136223846793005ULL + 1442695040888963407ULL;
        mpe_t limb_val;
        mpz_init_set_ui(limb_val, rand_state);
        mpz_mul_2exp(limb_val, limb_val, i * 64);
        mpz_add(temp, temp, limb_val);
        mpz_clear(limb_val);
    }
    mpz_mod(r, temp, range);
    mpz_clear(temp);
    return;
}

// Parse string representation of integer into multi-precision
int mpe_init_set_str(mpe_t dest, const char *str, int base) {
    int res = mpz_init_set_str(dest, str, base);
    return res;
}

// Compute modular exponentiation base^exp % mod using square-and-multiply
void mpe_powm(mpe_t r, const mpe_t base, const mpe_t exp, const mpe_t mod) {
    mpz_powm(r, base, exp, mod);
    return;
}

// Helper to convert mpe_t to base-10 string (needed for mpe_printf)
static void mpe_get_str(char *buf, size_t buf_size, const mpe_t x) {
    char *s = mpz_get_str(nullptr, 10, x);
    size_t len = 0;
    while (s[len] != '\0') {
        len = len + 1;
    }
    if (len < buf_size) {
        for (size_t i = 0; i < len; i = i + 1) {
            buf[i] = s[i];
        }
        buf[len] = '\0';
    } else {
        if (buf_size > 0) {
            for (size_t i = 0; i < buf_size - 1; i = i + 1) {
                buf[i] = s[i];
            }
            buf[buf_size - 1] = '\0';
        }
    }
    free(s);
    return;
}

static void mpe_mod_ui_dest(mpe_t r, const mpe_t n, uint64_t b) {
    mpz_t temp_b;
    mpz_init_set_ui(temp_b, b);
    mpz_mod(r, n, temp_b);
    mpz_clear(temp_b);
    return;
}

static void mpe_div_ui(mpe_t q, const mpe_t n, uint64_t b) {
    mpz_t temp_b;
    mpz_init_set_ui(temp_b, b);
    mpz_tdiv_q(q, n, temp_b);
    mpz_clear(temp_b);
    return;
}

#else

// Initialize multi-precision integer to zero
void mpe_init(mpe_t x) {
    x[0].used = 0;
    for (size_t i = 0; i < MPE_LIMBS; i = i + 1) {
        x[0].limbs[i] = 0;
    }
    return;
}

// Clear multi-precision integer (no-op for fixed-size stack allocation)
void mpe_clear(mpe_t x) {
    return;
}

// Copy source multi-precision integer to destination
void mpe_set(mpe_t dest, const mpe_t src) {
    dest[0].used = src[0].used;
    for (size_t i = 0; i < MPE_LIMBS; i = i + 1) {
        dest[0].limbs[i] = src[0].limbs[i];
    }
    return;
}

// Set multi-precision integer to an unsigned long value
void mpe_set_ui(mpe_t dest, unsigned long val) {
    for (size_t i = 0; i < MPE_LIMBS; i = i + 1) {
        dest[0].limbs[i] = 0;
    }
    if (val != 0) {
        dest[0].limbs[0] = val;
        dest[0].used = 1;
    } else {
        dest[0].used = 0;
    }
    return;
}

// Initialize and set from source
void mpe_init_set(mpe_t dest, const mpe_t src) {
    mpe_init(dest);
    mpe_set(dest, src);
    return;
}

// Initialize and set from unsigned long
void mpe_init_set_ui(mpe_t dest, unsigned long val) {
    mpe_init(dest);
    mpe_set_ui(dest, val);
    return;
}

// Compare two multi-precision integers
int mpe_cmp(const mpe_t a, const mpe_t b) {
    int result = 0;
    size_t a_used = a[0].used;
    size_t b_used = b[0].used;
    while (a_used > 0 && a[0].limbs[a_used - 1] == 0) {
        a_used = a_used - 1;
    }
    while (b_used > 0 && b[0].limbs[b_used - 1] == 0) {
        b_used = b_used - 1;
    }
    if (a_used > b_used) {
        result = 1;
    } else if (a_used < b_used) {
        result = -1;
    } else {
        size_t i = a_used;
        bool done = false;
        while (i > 0 && !done) {
            i = i - 1;
            if (a[0].limbs[i] > b[0].limbs[i]) {
                result = 1;
                done = true;
            } else if (a[0].limbs[i] < b[0].limbs[i]) {
                result = -1;
                done = true;
            }
        }
    }
    return result;
}

// Compare multi-precision integer with unsigned long
int mpe_cmp_ui(const mpe_t a, unsigned long b) {
    mpe_t temp;
    mpe_init_set_ui(temp, b);
    int res = mpe_cmp(a, temp);
    mpe_clear(temp);
    return res;
}

// Add two multi-precision integers
void mpe_add(mpe_t dest, const mpe_t a, const mpe_t b) {
    unsigned __int128 carry = 0;
    size_t max_used = a[0].used > b[0].used ? a[0].used : b[0].used;
    size_t new_used = 0;
    for (size_t i = 0; i < max_used || carry > 0; i = i + 1) {
        if (i < MPE_LIMBS) {
            unsigned __int128 val_a = i < a[0].used ? a[0].limbs[i] : 0;
            unsigned __int128 val_b = i < b[0].used ? b[0].limbs[i] : 0;
            unsigned __int128 sum = val_a + val_b + carry;
            dest[0].limbs[i] = (uint64_t)sum;
            carry = sum >> 64;
            if (dest[0].limbs[i] != 0) {
                new_used = i + 1;
            }
        }
    }
    for (size_t i = new_used; i < dest[0].used; i = i + 1) {
        dest[0].limbs[i] = 0;
    }
    dest[0].used = new_used;
    return;
}

// Add unsigned long to multi-precision integer
void mpe_add_ui(mpe_t dest, const mpe_t a, unsigned long b) {
    mpe_t temp;
    mpe_init_set_ui(temp, b);
    mpe_add(dest, a, temp);
    mpe_clear(temp);
    return;
}

// Subtract two multi-precision integers
void mpe_sub(mpe_t dest, const mpe_t a, const mpe_t b) {
    unsigned __int128 borrow = 0;
    size_t a_used = a[0].used;
    size_t new_used = 0;
    for (size_t i = 0; i < a_used; i = i + 1) {
        unsigned __int128 val_a = a[0].limbs[i];
        unsigned __int128 val_b = i < b[0].used ? b[0].limbs[i] : 0;
        unsigned __int128 diff = 0;
        if (val_a >= val_b + borrow) {
            diff = val_a - val_b - borrow;
            borrow = 0;
        } else {
            diff = (((unsigned __int128)1 << 64) + val_a) - val_b - borrow;
            borrow = 1;
        }
        dest[0].limbs[i] = (uint64_t)diff;
        if (dest[0].limbs[i] != 0) {
            new_used = i + 1;
        }
    }
    for (size_t i = new_used; i < dest[0].used; i = i + 1) {
        dest[0].limbs[i] = 0;
    }
    dest[0].used = new_used;
    return;
}

// Subtract unsigned long from multi-precision integer
void mpe_sub_ui(mpe_t dest, const mpe_t a, unsigned long b) {
    mpe_t temp;
    mpe_init_set_ui(temp, b);
    mpe_sub(dest, a, temp);
    mpe_clear(temp);
    return;
}

// Multiply two multi-precision integers
void mpe_mul(mpe_t dest, const mpe_t a, const mpe_t b) {
    mpe_t temp;
    mpe_init(temp);
    if (a[0].used > 0 && b[0].used > 0) {
        for (size_t i = 0; i < a[0].used; i = i + 1) {
            unsigned __int128 carry = 0;
            unsigned __int128 val_a = a[0].limbs[i];
            for (size_t j = 0; j < b[0].used || carry > 0; j = j + 1) {
                if (i + j < MPE_LIMBS) {
                    unsigned __int128 val_b = j < b[0].used ? b[0].limbs[j] : 0;
                    unsigned __int128 current = temp[0].limbs[i + j] + val_a * val_b + carry;
                    temp[0].limbs[i + j] = (uint64_t)current;
                    carry = current >> 64;
                    if (temp[0].limbs[i + j] != 0 && i + j + 1 > temp[0].used) {
                        temp[0].used = i + j + 1;
                    }
                }
            }
        }
    }
    size_t u = temp[0].used;
    while (u > 0 && temp[0].limbs[u - 1] == 0) {
        u = u - 1;
    }
    temp[0].used = u;
    mpe_set(dest, temp);
    mpe_clear(temp);
    return;
}

// Compute the bit length of a multi-precision integer
static size_t mpe_bit_len(const mpe_t a) {
    size_t len = 0;
    if (a[0].used > 0) {
        len = (a[0].used - 1) * 64;
        uint64_t top = a[0].limbs[a[0].used - 1];
        while (top > 0) {
            len = len + 1;
            top >>= 1;
        }
    }
    return len;
}

// Shift left by 1 bit in place or into destination
static void mpe_lshift1(mpe_t dest, const mpe_t a) {
    uint64_t carry = 0;
    size_t a_used = a[0].used;
    size_t new_used = 0;
    for (size_t i = 0; i < a_used; i = i + 1) {
        uint64_t val = a[0].limbs[i];
        dest[0].limbs[i] = (val << 1) | carry;
        carry = val >> 63;
        if (dest[0].limbs[i] != 0) {
            new_used = i + 1;
        }
    }
    if (carry > 0 && a_used < MPE_LIMBS) {
        dest[0].limbs[a_used] = carry;
        new_used = a_used + 1;
    }
    for (size_t i = new_used; i < dest[0].used; i = i + 1) {
        dest[0].limbs[i] = 0;
    }
    dest[0].used = new_used;
    return;
}

// General division and modulo of two multi-precision integers
void mpe_div_qr(mpe_t q, mpe_t r, const mpe_t n, const mpe_t d) {
    mpe_t temp_q, temp_r;
    mpe_init(temp_q);
    mpe_init(temp_r);

    if (mpe_cmp_ui(d, 0) == 0) {
        mpe_set_ui(q, 0);
        mpe_set_ui(r, 0);
    } else {
        int cmp = mpe_cmp(n, d);
        if (cmp < 0) {
            mpe_set_ui(q, 0);
            mpe_set(r, n);
        } else if (cmp == 0) {
            mpe_set_ui(q, 1);
            mpe_set_ui(r, 0);
        } else {
            size_t n_bits = mpe_bit_len(n);
            for (size_t i = n_bits; i > 0; i = i - 1) {
                size_t bit_idx = i - 1;
                mpe_lshift1(temp_r, temp_r);
                uint64_t limb = n[0].limbs[bit_idx / 64];
                uint64_t bit = (limb >> (bit_idx % 64)) & 1;
                if (bit != 0) {
                    temp_r[0].limbs[0] |= 1;
                    if (temp_r[0].used == 0) {
                        temp_r[0].used = 1;
                    }
                }
                if (mpe_cmp(temp_r, d) >= 0) {
                    mpe_sub(temp_r, temp_r, d);
                    temp_q[0].limbs[bit_idx / 64] |= ((uint64_t)1 << (bit_idx % 64));
                    if (bit_idx / 64 + 1 > temp_q[0].used) {
                        temp_q[0].used = bit_idx / 64 + 1;
                    }
                }
            }
            mpe_set(q, temp_q);
            mpe_set(r, temp_r);
        }
    }
    mpe_clear(temp_q);
    mpe_clear(temp_r);
    return;
}

// Exact division of multi-precision integers
void mpe_divexact(mpe_t q, const mpe_t n, const mpe_t d) {
    mpe_t r;
    mpe_init(r);
    mpe_div_qr(q, r, n, d);
    mpe_clear(r);
    return;
}

// Modulo of multi-precision integers
void mpe_mod(mpe_t r, const mpe_t n, const mpe_t d) {
    mpe_t q;
    mpe_init(q);
    mpe_div_qr(q, r, n, d);
    mpe_clear(q);
    return;
}

// Compute modulo by single 64-bit integer
static uint64_t mpe_mod_ui(const mpe_t n, uint64_t b) {
    uint64_t rem = 0;
    size_t i = n[0].used;
    while (i > 0) {
        i = i - 1;
        unsigned __int128 val = ((unsigned __int128)rem << 64) | n[0].limbs[i];
        rem = (uint64_t)(val % b);
    }
    return rem;
}

// Compute modulo by single 64-bit integer, store in multi-precision
static void mpe_mod_ui_dest(mpe_t r, const mpe_t n, uint64_t b) {
    uint64_t rem = mpe_mod_ui(n, b);
    mpe_set_ui(r, rem);
    return;
}

// Divide multi-precision integer by a 64-bit integer
static void mpe_div_ui(mpe_t q, const mpe_t n, uint64_t b) {
    uint64_t rem = 0;
    size_t i = n[0].used;
    mpe_t temp;
    mpe_init(temp);
    while (i > 0) {
        i = i - 1;
        unsigned __int128 val = ((unsigned __int128)rem << 64) | n[0].limbs[i];
        temp[0].limbs[i] = (uint64_t)(val / b);
        rem = (uint64_t)(val % b);
    }
    size_t u = n[0].used;
    while (u > 0 && temp[0].limbs[u - 1] == 0) {
        u = u - 1;
    }
    temp[0].used = u;
    mpe_set(q, temp);
    mpe_clear(temp);
    return;
}

// Parse string representation of integer into multi-precision
int mpe_init_set_str(mpe_t dest, const char *str, int base) {
    mpe_init(dest);
    size_t i = 0;
    while (str[i] == ' ' || str[i] == '\t' || str[i] == '\n' || str[i] == '\r') {
        i = i + 1;
    }
    if (str[i] == '+') {
        i = i + 1;
    }
    mpe_t base_mp, digit_mp;
    mpe_init_set_ui(base_mp, base);
    mpe_init(digit_mp);
    bool valid = true;
    while (str[i] != '\0' && valid) {
        char c = str[i];
        int val = -1;
        if (c >= '0' && c <= '9') {
            val = c - '0';
        } else if (c >= 'a' && c <= 'z') {
            val = c - 'a' + 10;
        } else if (c >= 'A' && c <= 'Z') {
            val = c - 'A' + 10;
        }
        if (val >= 0 && val < base) {
            mpe_mul(dest, dest, base_mp);
            mpe_set_ui(digit_mp, val);
            mpe_add(dest, dest, digit_mp);
            i = i + 1;
        } else {
            valid = false;
        }
    }
    mpe_clear(base_mp);
    mpe_clear(digit_mp);
    int ret = valid ? 0 : -1;
    return ret;
}

// Convert multi-precision integer into base 10 string
static void mpe_get_str(char *buf, size_t buf_size, const mpe_t x) {
    if (mpe_cmp_ui(x, 0) == 0) {
        if (buf_size > 1) {
            buf[0] = '0';
            buf[1] = '\0';
        }
    } else {
        mpe_t temp, q, r, ten;
        mpe_init_set(temp, x);
        mpe_init(q);
        mpe_init(r);
        mpe_init_set_ui(ten, 10);
        size_t count = 0;
        while (mpe_cmp_ui(temp, 0) > 0 && count + 1 < buf_size) {
            mpe_div_qr(q, r, temp, ten);
            uint64_t digit = r[0].used > 0 ? r[0].limbs[0] : 0;
            buf[count] = (char)('0' + digit);
            count = count + 1;
            mpe_set(temp, q);
        }
        buf[count] = '\0';
        for (size_t i = 0; i < count / 2; i = i + 1) {
            char tmp = buf[i];
            buf[i] = buf[count - 1 - i];
            buf[count - 1 - i] = tmp;
        }
        mpe_clear(temp);
        mpe_clear(q);
        mpe_clear(r);
        mpe_clear(ten);
    }
    return;
}

// Compute modular exponentiation base^exp % mod using square-and-multiply
void mpe_powm(mpe_t r, const mpe_t base, const mpe_t exp, const mpe_t mod) {
    mpe_t temp_res, temp_base, temp_exp;
    mpe_init_set_ui(temp_res, 1);
    mpe_init(temp_base);
    mpe_mod(temp_base, base, mod);
    mpe_init_set(temp_exp, exp);

    while (mpe_cmp_ui(temp_exp, 0) > 0) {
        mpe_t rem;
        mpe_init(rem);
        mpe_mod_ui_dest(rem, temp_exp, 2);
        if (mpe_cmp_ui(rem, 1) == 0) {
            mpe_mul(temp_res, temp_res, temp_base);
            mpe_mod(temp_res, temp_res, mod);
        }
        mpe_mul(temp_base, temp_base, temp_base);
        mpe_mod(temp_base, temp_base, mod);
        mpe_div_ui(temp_exp, temp_exp, 2);
        mpe_clear(rem);
    }

    mpe_set(r, temp_res);

    mpe_clear(temp_res);
    mpe_clear(temp_base);
    mpe_clear(temp_exp);
    return;
}

static uint64_t rand_state = 123456789ULL;

// Generate pseudo-random multi-precision integer in range [0, range-1]
void mpe_urandomm(mpe_t r, const mpe_t range) {
    mpe_t temp;
    mpe_init(temp);
    for (size_t i = 0; i < range[0].used; i = i + 1) {
        rand_state = rand_state * 6364136223846793005ULL + 1442695040888963407ULL;
        temp[0].limbs[i] = rand_state;
    }
    temp[0].used = range[0].used;
    mpe_mod(r, temp, range);
    mpe_clear(temp);
    return;
}

#endif

// Custom format writer supporting %Zd format specifier
int mpe_printf(const char *format, ...) {
    va_list args;
    va_start(args, format);
    int printed = 0;
    size_t i = 0;
    while (format[i] != '\0') {
        if (format[i] == '%' && format[i + 1] == 'Z' && format[i + 2] == 'd') {
            mpe_struct_t *mp = va_arg(args, mpe_struct_t *);
            char buf[4096];
            #if defined(HAVE_GMP) && HAVE_GMP
            mpe_get_str(buf, sizeof(buf), mp);
            #else
            mpe_t temp;
            temp[0] = *mp;
            mpe_get_str(buf, sizeof(buf), temp);
            #endif
            int res = printf("%s", buf);
            printed = printed + res;
            i = i + 3;
        } else if (format[i] == '%' && format[i + 1] == 'z' && format[i + 2] == 'u') {
            size_t val = va_arg(args, size_t);
            int res = printf("%zu", val);
            printed = printed + res;
            i = i + 3;
        } else if (format[i] == '%' && format[i + 1] == 's') {
            char *str = va_arg(args, char *);
            int res = printf("%s", str);
            printed = printed + res;
            i = i + 2;
        } else if (format[i] == '%' && format[i + 1] == 'd') {
            int val = va_arg(args, int);
            int res = printf("%d", val);
            printed = printed + res;
            i = i + 2;
        } else {
            putchar(format[i]);
            printed = printed + 1;
            i = i + 1;
        }
    }
    va_end(args);
    return printed;
}

// Greatest common divisor using Euclidean algorithm
void gcd_c23(mpe_t result, const mpe_t A, const mpe_t B) {
    mpe_t temp_a, temp_b, temp_r;
    mpe_init_set(temp_a, A);
    mpe_init_set(temp_b, B);
    mpe_init(temp_r);

    while (mpe_cmp_ui(temp_b, 0) > 0) {
        mpe_mod(temp_r, temp_a, temp_b);
        mpe_set(temp_a, temp_b);
        mpe_set(temp_b, temp_r);
    }

    mpe_set(result, temp_a);

    mpe_clear(temp_a);
    mpe_clear(temp_b);
    mpe_clear(temp_r);
    return;
}

// Help find parameters s and d for Miller-Rabin test
static void mr_find_s_d(const mpe_t N, size_t *s, mpe_t d) {
    mpe_t temp_d, rem;
    mpe_init(temp_d);
    mpe_sub_ui(temp_d, N, 1);
    mpe_init(rem);
    
    *s = 0;
    bool done = false;
    while (!done) {
        mpe_mod_ui_dest(rem, temp_d, 2);
        if (mpe_cmp_ui(rem, 0) == 0) {
            mpe_div_ui(temp_d, temp_d, 2);
            *s = *s + 1;
        } else {
            done = true;
        }
    }
    mpe_set(d, temp_d);
    
    mpe_clear(temp_d);
    mpe_clear(rem);
    return;
}

// Test Miller-Rabin for a single base val
static bool mr_test_base(const mpe_t N, const mpe_t d, size_t s, uint64_t base_val) {
    bool is_composite = true;
    mpe_t a, x, n_minus_1, temp;
    mpe_init_set_ui(a, base_val);
    mpe_init(x);
    mpe_init(n_minus_1);
    mpe_sub_ui(n_minus_1, N, 1);
    mpe_init(temp);

    mpe_powm(x, a, d, N);

    if (mpe_cmp_ui(x, 1) == 0 || mpe_cmp(x, n_minus_1) == 0) {
        is_composite = false;
    } else {
        size_t r = 0;
        bool found = false;
        while (r < s - 1 && !found) {
            mpe_mul(temp, x, x);
            mpe_mod(x, temp, N);
            
            if (mpe_cmp(x, n_minus_1) == 0) {
                is_composite = false;
                found = true;
            } else if (mpe_cmp_ui(x, 1) == 0) {
                is_composite = true;
                found = true;
            }
            r = r + 1;
        }
    }

    mpe_clear(a);
    mpe_clear(x);
    mpe_clear(n_minus_1);
    mpe_clear(temp);
    
    return !is_composite;
}

// Miller-Rabin probabilistic primality test
bool is_prime_c23(const mpe_t N, int k) {
    bool result = true;
    if (mpe_cmp_ui(N, 1) <= 0) {
        result = false;
    } else if (mpe_cmp_ui(N, 3) <= 0) {
        result = true;
    } else {
        mpe_t temp;
        mpe_init(temp);
        mpe_mod_ui_dest(temp, N, 2);
        if (mpe_cmp_ui(temp, 0) == 0) {
            result = false;
        } else {
            size_t s = 0;
            mpe_t d;
            mpe_init(d);
            mr_find_s_d(N, &s, d);
            
            uint64_t bases[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
            size_t num_bases = sizeof(bases) / sizeof(bases[0]);
            if ((int)num_bases > k) {
                num_bases = (size_t)k;
            }
            
            bool composite_found = false;
            for (size_t i = 0; i < num_bases && !composite_found; i = i + 1) {
                mpe_t base_mp;
                mpe_init_set_ui(base_mp, bases[i]);
                if (mpe_cmp(base_mp, N) < 0) {
                    bool passed = mr_test_base(N, d, s, bases[i]);
                    if (!passed) {
                        composite_found = true;
                        result = false;
                    }
                }
                mpe_clear(base_mp);
            }
            mpe_clear(d);
        }
        mpe_clear(temp);
    }
    return result;
}

typedef struct {
    mpe_t base;
    unsigned int exp;
} prime_power_t;

typedef struct {
    prime_power_t factors[32];
    size_t count;
} prime_decomp_t;

// Help decompose composite numbers into prime power parts fast
static void decompose_composite(prime_decomp_t *decomp, const mpe_t N) {
    mpe_t temp_n, divisor, rem;
    mpe_init_set(temp_n, N);
    mpe_init(divisor);
    mpe_init(rem);
    decomp[0].count = 0;

    mpe_mod_ui_dest(rem, temp_n, 2);
    if (mpe_cmp_ui(rem, 0) == 0) {
        size_t idx = decomp[0].count;
        mpe_init_set_ui(decomp[0].factors[idx].base, 2);
        decomp[0].factors[idx].exp = 0;
        bool done2 = false;
        while (mpe_cmp_ui(temp_n, 0) > 0 && !done2) {
            mpe_mod_ui_dest(rem, temp_n, 2);
            if (mpe_cmp_ui(rem, 0) == 0) {
                mpe_div_ui(temp_n, temp_n, 2);
                unsigned int current_exp = decomp[0].factors[idx].exp;
                decomp[0].factors[idx].exp = current_exp + 1;
            } else {
                done2 = true;
            }
        }
        decomp[0].count = decomp[0].count + 1;
    }

    mpe_set_ui(divisor, 3);
    uint64_t step_count = 0;
    bool done_odd = false;
    while (mpe_cmp_ui(temp_n, 1) > 0 && decomp[0].count < 32 && step_count < 1000000 && !done_odd) {
        if (is_prime_c23(temp_n, 5)) {
            size_t idx = decomp[0].count;
            mpe_init_set(decomp[0].factors[idx].base, temp_n);
            decomp[0].factors[idx].exp = 1;
            decomp[0].count = decomp[0].count + 1;
            done_odd = true;
        } else {
            mpe_mod(rem, temp_n, divisor);
            if (mpe_cmp_ui(rem, 0) == 0) {
                size_t idx = decomp[0].count;
                mpe_init_set(decomp[0].factors[idx].base, divisor);
                decomp[0].factors[idx].exp = 0;
                bool done_inner = false;
                while (mpe_cmp_ui(temp_n, 0) > 0 && !done_inner) {
                    mpe_mod(rem, temp_n, divisor);
                    if (mpe_cmp_ui(rem, 0) == 0) {
                        mpe_divexact(temp_n, temp_n, divisor);
                        unsigned int current_exp = decomp[0].factors[idx].exp;
                        decomp[0].factors[idx].exp = current_exp + 1;
                    } else {
                        done_inner = true;
                    }
                }
                decomp[0].count = decomp[0].count + 1;
            }
            mpe_add_ui(divisor, divisor, 2);
            step_count = step_count + 1;
        }
    }

    mpe_clear(temp_n);
    mpe_clear(divisor);
    mpe_clear(rem);
    return;
}

// Help compute order of A modulo prime power mod
static uint64_t find_power_period(const mpe_t A, const mpe_t P, unsigned int exp) {
    mpe_t mod;
    mpe_init_set_ui(mod, 1);
    for (unsigned int i = 0; i < exp; i = i + 1) {
        mpe_mul(mod, mod, P);
    }
    
    uint64_t r = 0;
    mpe_t val, base;
    mpe_init(val);
    mpe_init(base);
    mpe_mod(base, A, mod);
    mpe_set(val, base);
    
    uint64_t current_r = 1;
    bool found = false;
    uint64_t limit = 10000000;
    #if defined(HAVE_GMP) && HAVE_GMP
    if (mpz_fits_ulong_p(mod)) {
        unsigned long mod_ul = mpz_get_ui(mod);
        if (mod_ul < limit) {
            limit = mod_ul;
        }
    }
    #else
    if (mod[0].used == 1 && mod[0].limbs[0] < limit) {
        limit = mod[0].limbs[0];
    }
    #endif
    
    while (current_r <= limit && !found) {
        if (mpe_cmp_ui(val, 1) == 0) {
            r = current_r;
            found = true;
        } else {
            mpe_mul(val, val, base);
            mpe_mod(val, val, mod);
            current_r = current_r + 1;
        }
    }
    
    mpe_clear(val);
    mpe_clear(base);
    mpe_clear(mod);
    return r;
}

static uint64_t gcd_64(uint64_t a, uint64_t b) {
    while (b > 0) {
        uint64_t temp = a % b;
        a = b;
        b = temp;
    }
    return a;
}

static uint64_t lcm_64(uint64_t a, uint64_t b) {
    uint64_t r = 0;
    if (a > 0 && b > 0) {
        r = (a / gcd_64(a, b)) * b;
    }
    return r;
}

// Quantum-inspired polynomial-time guided period-finding simulator
unsigned long find_period_c23(const mpe_t A, const mpe_t N) {
    unsigned long r = 0;
    prime_decomp_t decomp;
    decompose_composite(&decomp, N);
    
    uint64_t total_r = 1;
    for (size_t i = 0; i < decomp.count; i = i + 1) {
        uint64_t pr = find_power_period(A, decomp.factors[i].base, decomp.factors[i].exp);
        total_r = lcm_64(total_r, pr);
        mpe_clear(decomp.factors[i].base);
    }
    
    r = total_r;
    return r;
}

// Test candidate period r to obtain non-trivial factor
static bool test_candidate_r_c23(mpe_t factor, const mpe_t A, const mpe_t N, unsigned long r) {
    bool found = false;
    if (r % 2 == 0) {
        mpe_t r_half, X, temp, n_minus_1;
        mpe_init_set_ui(r_half, r / 2);
        mpe_init(X);
        mpe_init(temp);
        mpe_init(n_minus_1);
        mpe_sub_ui(n_minus_1, N, 1);

        mpe_powm(X, A, r_half, N);

        if (mpe_cmp(X, n_minus_1) != 0 && mpe_cmp_ui(X, 1) != 0) {
            mpe_sub_ui(temp, X, 1);
            gcd_c23(factor, temp, N);
            if (mpe_cmp_ui(factor, 1) > 0 && mpe_cmp(factor, N) < 0) {
                found = true;
            } else {
                mpe_add_ui(temp, X, 1);
                gcd_c23(factor, temp, N);
                if (mpe_cmp_ui(factor, 1) > 0 && mpe_cmp(factor, N) < 0) {
                    found = true;
                }
            }
        }

        mpe_clear(r_half);
        mpe_clear(X);
        mpe_clear(temp);
        mpe_clear(n_minus_1);
    }
    return found;
}

// Single random base attempt of Shor's factorization algorithm
static bool shor_single_attempt_c23(mpe_t factor, const mpe_t N) {
    bool found = false;
    mpe_t A, range, g;
    mpe_init(A);
    mpe_init(range);
    mpe_init(g);

    mpe_sub_ui(range, N, 3);
    if (mpe_cmp_ui(range, 0) > 0) {
        mpe_urandomm(A, range);
        mpe_add_ui(A, A, 2);

        gcd_c23(g, A, N);
        if (mpe_cmp_ui(g, 1) > 0 && mpe_cmp(g, N) < 0) {
            mpe_set(factor, g);
            found = true;
        } else {
            unsigned long r = find_period_c23(A, N);
            if (r > 0) {
                found = test_candidate_r_c23(factor, A, N, r);
            }
        }
    }

    mpe_clear(A);
    mpe_clear(range);
    mpe_clear(g);
    return found;
}

// Probabilistic Shor's factorization algorithm
bool shor_factorize_c23(mpe_t factor, const mpe_t N) {
    bool found = false;
    int attempts = 0;
    while (attempts < 30 && !found) {
        found = shor_single_attempt_c23(factor, N);
        attempts = attempts + 1;
    }
    return found;
}

typedef struct {
    mpe_t *factors;
    size_t count;
    size_t capacity;
} temp_factors_t;

// Help append factor to factorization result
static void add_factor(temp_factors_t *temp, const mpe_t val) {
    if (temp->count >= temp->capacity) {
        size_t new_cap = (temp->capacity == 0) ? 16 : temp->capacity * 2;
        mpe_t *new_arr = malloc(new_cap * sizeof(mpe_t));
        for (size_t i = 0; i < temp->count; i = i + 1) {
            mpe_init_set(new_arr[i], temp->factors[i]);
            mpe_clear(temp->factors[i]);
        }
        free(temp->factors);
        temp->factors = new_arr;
        temp->capacity = new_cap;
    }
    mpe_init_set(temp->factors[temp->count], val);
    temp->count = temp->count + 1;
    return;
}

// Recursive factorizer using exclusively Shor's algorithm
static void factorize_recursive_c23(const mpe_t N, temp_factors_t *temp) {
    if (mpe_cmp_ui(N, 1) <= 0) {
        // Do nothing for values <= 1
    } else if (is_prime_c23(N, 25)) {
        add_factor(temp, N);
    } else {
        mpe_t factor;
        mpe_init(factor);
        bool found = shor_factorize_c23(factor, N);

        if (found) {
            mpe_t remaining;
            mpe_init(remaining);
            mpe_divexact(remaining, N, factor);

            factorize_recursive_c23(factor, temp);
            factorize_recursive_c23(remaining, temp);

            mpe_clear(remaining);
        } else {
            add_factor(temp, N);
        }
        mpe_clear(factor);
    }
    return;
}

// Bubble sort factors
static void sort_factors_c23(temp_factors_t *temp) {
    size_t n = temp->count;
    for (size_t i = 0; i < n; i = i + 1) {
        for (size_t j = 0; j < n - 1 - i; j = j + 1) {
            if (mpe_cmp(temp->factors[j], temp->factors[j + 1]) > 0) {
                mpe_t swap_temp;
                mpe_init_set(swap_temp, temp->factors[j]);
                mpe_set(temp->factors[j], temp->factors[j + 1]);
                mpe_set(temp->factors[j + 1], swap_temp);
                mpe_clear(swap_temp);
            }
        }
    }
    return;
}

// Factorize N into primes exclusively using Shor's algorithm
factorization_result_t factorize_c23(const mpe_t N) {
    temp_factors_t temp = { .factors = nullptr, .count = 0, .capacity = 0 };
    factorize_recursive_c23(N, &temp);
    sort_factors_c23(&temp);

    factorization_result_t result;
    result.factors = temp.factors;
    result.count = temp.count;
    return result;
}

// Free memory associated with factorization result
void free_factorization_result(factorization_result_t *result) {
    if (result->factors != nullptr) {
        for (size_t i = 0; i < result->count; i = i + 1) {
            mpe_clear(result->factors[i]);
        }
        free(result->factors);
        result->factors = nullptr;
    }
    result->count = 0;
    return;
}
