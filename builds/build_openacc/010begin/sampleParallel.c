//
//  ˆêŸŒ³”z—ñ‚É’è”‚ğæZ‚·‚éB
//
// (c)Copyright Spacesoft corp., 2018 All rights reserved.
//                                        Hiro KITAYAMA
#include <stdio.h>
#include <stdlib.h>

#define n   65536

int main()
{
    float x[n], y[n];

    for (int i = 0; i < n; i++)
        x[i] = (float)rand();

    float a = (float)rand();

    #pragma acc parallel
    for (int i = 0; i < n; i++)
    {
        y[i] = a * x[i];
    }

    for (int i = 0; i < 10; i++)
    {
        printf("y[%d] = %.10f\n", i, y[i]);
    }
    return 0;
}
