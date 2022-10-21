//
// multiply two arrays and store them in another array,OpenACC Data
//
//  (c)Copyright Spacesoft corp., 2018 All rights reserved.
//                                Hiro KITAYAMA
//
#include <stdio.h>

void verify(const int n, const float* a, const float *x, const float *y);

#define N   4096

// calc.
void calc(const int n, const float* a, const float *b, float *c)
{
    int i;

    #pragma acc kernels present(a, b, c)
    #pragma acc loop independent
    for (i = 0; i < n; i++)
    {
        c[i] = a[i] * b[i];
    }
}

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

    #pragma acc data copyin(a, b) copyout(c)
    {
        calc(N, a, b, c);
    }

    // list results
    printf("(a * b = c)\n");
    for (i = 0; i < 10; i++)
        printf("%f * %f = %f\n", a[i], b[i], c[i]);

    verify(N, a, b, c);

    return 0;
}
