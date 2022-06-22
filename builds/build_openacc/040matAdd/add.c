//
//  add two matrixs and store the result in another matrix
//
//  (c)Copyright Spacesoft corp., 2018 All rights reserved.
//                                Hiro KITAYAMA
//
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// main
int
main(int argc, char* argv[])
{
    float **a, **b, **c;
    clock_t start, stop;
    int i, j, n = 104096;

    a = (float**)malloc(sizeof(float *) * n);
    b = (float**)malloc(sizeof(float *) * n);
    c = (float**)malloc(sizeof(float *) * n);
    for (i = 0; i < n; i++)
    {
        a[i] = (float*)malloc(sizeof(int) * n);
        b[i] = (float*)malloc(sizeof(int) * n);
        c[i] = (float*)malloc(sizeof(int) * n);
    }

    // initialize array
    for (j = 0; j < n; j++)
    {
        for (i = 0; i < n; i++)
        {
            a[j][i] = (float)(i + 1000);
            b[j][i] = (float)i / 10.f;
        }
    }


    start = clock();

    // calc.
    for (j = 0; j < n; j++)
    {
        for (i = 0; i < n; i++)
        {
            c[j][i] = a[j][i] + b[j][i];
        }
    }

    stop = clock();

    fprintf(stdout, "      C: ");
    fprintf(stdout, "elapsed time = %.20f [sec]\n",
        (float)(stop - start) / CLOCKS_PER_SEC);


    start = clock();

    // calc.
    #pragma acc kernels
    for (j = 0; j < n; j++)
    {
        for (i = 0; i < n; i++)
        {
            c[j][i] = a[j][i] + b[j][i];
        }
    }

    stop = clock();

    fprintf(stdout, "OpenACC: ");
    fprintf(stdout, "elapsed time = %.20f [sec]\n",
        (float)(stop - start) / CLOCKS_PER_SEC);

    for (i = 0; i < n; i++)
    {
        free(a[i]);
        free(b[i]);
        free(c[i]);
    }
    free(a);
    free(b);
    free(c);

    return 0;
}
