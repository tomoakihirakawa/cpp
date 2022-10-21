//
// multiply two arrays and store them in another array, for AVX(vector)
//
// (c)Copyright Spacesoft corp., 2018 All rights reserved.
//                               Hiro KITAYAMA
//
#include <stdio.h>
#include <immintrin.h>

void verify(const int n, const float* a, const float *x, const float *y);

#define N   4096

// main
int
main()
{
    float a[N], b[N], c[N];
    int i;

    // initialize array
    for (i = 0; i < N; i++)
    {
        a[i] = (float)(i + 1000);
        b[i] = (float)i / 10.f;
    }

    // calc.
    for (int i = 0; i < N; i += sizeof(__m256) / sizeof(float))
    {
        __m256 va = _mm256_loadu_ps(&a[i]);
        __m256 vb = _mm256_loadu_ps(&b[i]);
        __m256 vc = _mm256_mul_ps(va, vb);
        _mm256_storeu_ps(&c[i], vc);
    }

    // list results
    printf("(a * b = c)\n");
    for (i = 0; i < 10; i++)
        printf("%f * %f = %f\n", a[i], b[i], c[i]);

    verify(N, a, b, c);

    return 0;
}
