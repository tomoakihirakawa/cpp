//
//  Add / Subtract two arrays and store the result in another array
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
    float *a, *b, *c, *d;
    clock_t start, stop;
    int i, j, n = 262144;

    if (argc > 1)
        n = atoi(argv[1]);

    fprintf(stdout, "matrix size = %d x %d\n", n, n);

    a = (float*)malloc(sizeof(float) * n);
    b = (float*)malloc(sizeof(float) * n);
    c = (float*)malloc(sizeof(float) * n);
    d = (float*)malloc(sizeof(float) * n);

    // initialize array
    for (i = 0; i < n; i++)
    {
        a[i] = (float)(i + 1000);
        b[i] = (float)i / 10.f;
    }


    start = clock();

    #pragma acc data copyin(a[:n], b[:n]) copyout(c[:n], d[:n])
    {

        #pragma acc kernels
        {
            // add
            #pragma acc loop independent
            for (i = 0; i < n; i++)
            {
                c[i] = a[i] + b[i];
            }

            // sub
            #pragma acc loop independent
            for (i = 0; i < n; i++)
            {
                d[i] = a[i] - b[i];
            }
        }

    }

    stop = clock();


    fprintf(stdout, "elapsed time = %.20f [sec]\n",
        (float)(stop - start) / CLOCKS_PER_SEC);

    free(a);
    free(b);
    free(c);
    free(d);

    return 0;
}
