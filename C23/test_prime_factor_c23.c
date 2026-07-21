#include "prime_factor_c23.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

typedef struct {
    const char *n_str;
    size_t expected_count;
} test_case_t;

// Test suite consisting of primes, composite numbers, and big integers
static const test_case_t TEST_SUITE[] = {
    { "4", 2 },
    { "12", 3 },
    { "15", 2 },
    { "143", 2 },
    { "1024", 10 },
    { "4087", 2 },
    { "999999", 7 },
    { "123456789", 4 },
    { "11111111111", 2 },
    { "17", 1 }
};

static const size_t NUM_TESTS = sizeof(TEST_SUITE) / sizeof(TEST_SUITE[0]);

bool run_test_case(const test_case_t *tc) {
    bool passed = true;
    mpe_t N;
    mpe_init_set_str(N, tc->n_str, 10);

    mpe_printf("Testing N = %Zd (Expected factors: %zu)...\n", N, tc->expected_count);

    factorization_result_t result = factorize_c23(N);

    if (result.count != tc->expected_count) {
        printf("  [FAIL] Expected %zu factors, but got %zu\n", tc->expected_count, result.count);
        passed = false;
    } else {
        mpe_t product;
        mpe_init_set_ui(product, 1);
        bool all_prime = true;

        for (size_t i = 0; i < result.count; i = i + 1) {
            mpe_printf("    Factor %zu: %Zd\n", i + 1, result.factors[i]);
            mpe_mul(product, product, result.factors[i]);
            if (!is_prime_c23(result.factors[i], 25)) {
                all_prime = false;
            }
        }

        if (mpe_cmp(product, N) != 0) {
            mpe_printf("  [FAIL] Product of factors (%Zd) does not equal N (%Zd)\n", product, N);
            passed = false;
        } else if (!all_prime) {
            printf("  [FAIL] Some factors are not prime\n");
            passed = false;
        } else {
            printf("  [PASS] All checks succeeded\n");
        }

        mpe_clear(product);
    }

    free_factorization_result(&result);
    mpe_clear(N);

    return passed;
}

int main() {
    bool all_passed = true;
    printf("=== RUNNING C23 PRIME FACTORIZATION TEST SUITE (WITH MPE) ===\n");

    for (size_t i = 0; i < NUM_TESTS; i = i + 1) {
        bool tc_passed = run_test_case(&TEST_SUITE[i]);
        if (!tc_passed) {
            all_passed = false;
        }
    }

    if (all_passed) {
        printf("=== ALL C23 TEST CASES PASSED SUCCESSFULLY ===\n");
    } else {
        printf("=== TEST SUITE FAILED ===\n");
    }

    int exit_code = 0;
    if (!all_passed) {
        exit_code = 1;
    }

    return exit_code;
}
