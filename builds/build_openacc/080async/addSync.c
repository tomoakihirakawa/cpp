//
//  add sync
//
//  (c)Copyright Spacesoft corp., 2018 All rights reserved.
//                                Hiro KITAYAMA
//
#include <stdio.h>
#include <stdlib.h>

// main
int
main(int argc, char* argv[])
{
    float *a, *b, *c, *d;
    int i, error = 0, n = 262144;

    if (argc > 1)
        n = atoi(argv[1]);

    a = malloc(sizeof(float) * n);
    b = malloc(sizeof(float) * n);
    c = malloc(sizeof(float) * n);

    // initialize array
    for (i = 0; i < n; i++)
    {
        a[i] = (float)(i + 1000);
        b[i] = (float)i / 10.f;
    }

    // add by accelerator
    #pragma acc parallel loop copyin(a[:n], b[:n]) copyout(c[:n])
    for (i = 0; i < n; i++)
    {
        c[i] = a[i] + b[i];
    }

    // add by host
    d = malloc(sizeof(float) * n);
    for (i = 0; i < n; i++)
    {
        d[i] = a[i] + b[i];
    }

    // verify
    for (i = 0; i < n; i++)
    {
        if (d[i] != c[i])
        {
            error = 1;
            break;
        }
    }
    if (error == 0)
        fprintf(stderr, "Passed.\n");
    else
        fprintf(stderr, "error!\n");

    free(a);
    free(b);
    free(c);
    free(d);

    return 0;
}
