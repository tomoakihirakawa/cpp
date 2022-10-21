//
// Leibniz pi ,for C, OpenACC
//
// (c)Copyright Spacesoft corp., 2018 All rights reserved.
//                               Hiro KITAYAMA
//
#include <stdio.h>
#include <math.h>
#include <time.h>

//----------------------------------------------------------------
void Leibniz(const int n, const int acc)
{
    clock_t start = clock();

    double pi = 0.0f;
    #pragma acc kernels if(acc)
    for (int i = 0; i < n; i++)
    {
        pi += (double)(pow(-1, i) / (double)(2 * i + 1));
    }

    pi *= 4.0f;

    clock_t stop = clock();

    fprintf(stdout, " n=%11d,", n);
    fprintf(stdout, "  elapsed time=%.10f [sec], pi=%.20f\n",
        (float)(stop - start) / CLOCKS_PER_SEC, pi);
}

//----------------------------------------------------------------
int main()
{
    for(int n = 1000000; n <= 100000000; n *= 10)
    {
        fprintf(stdout, "     C:");
        Leibniz(n, 0);
        fprintf(stdout, "OpenMP:");
        Leibniz(n, 1);
    }

    return 0;
}
