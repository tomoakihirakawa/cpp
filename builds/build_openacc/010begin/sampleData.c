//
//  �ꎟ���z��ɒ萔����Z����B
//
// (c)Copyright Spacesoft corp., 2018 All rights reserved.
//                                        Hiro KITAYAMA
//
#include <stdio.h>
#include <stdlib.h>

#define n   65536

int main()
{
    float x[n], y[n];

    for (int i = 0; i < n; i++)
        x[i] = (float)rand();

    float a = (float)rand();

    #pragma acc data copyin(x) copyout(y)
    {                                       // x: �z�X�g���f�o�C�X
        #pragma acc kernels
        for (int i = 0; i < n; i++)
        {
            y[i] = a * x[i];
        }
    }                                       // y: �f�o�C�X���z�X�g

    for (int i = 0; i < 10; i++)
    {
        printf("y[%d] = %.10f\n", i, y[i]);
    }
    return 0;
}
