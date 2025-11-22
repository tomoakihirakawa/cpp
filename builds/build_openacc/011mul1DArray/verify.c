//
// verify
//
// (c)Copyright Spacesoft corp., 2018 All rights reserved.
//                               Hiro KITAYAMA
//
#include <stdio.h>
#include <math.h>

void verify(const int n, const float* a, const float *b, const float *c)
{
    for (int i = 0; i < n; i++)
    {
        float cc = a[i] * b[i];
        if (fabs(cc - c[i]) > .000001f)
        {
            fprintf(stderr, "error: cc = %f, c[%d]=%f\n", cc, i, c[i]);
            return;
        }
    }
}
