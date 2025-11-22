//
// Leibniz pi ,for C, OpenMP
//
// (c)Copyright Spacesoft corp., 2018 All rights reserved.
//                               Hiro KITAYAMA
//
#include <stdio.h>
#include <math.h>
#include <time.h>

#include <omp.h>

//----------------------------------------------------------------
void Leibniz(const int n, const int omp)
{
    int i;
    clock_t start = clock();

    float pi = 0.0f;
    #pragma omp parallel for reduction(+:pi) if(omp)
    for (i = 0; i < n; i++)
    {
        pi += (float)(pow(-1, i) / (float)(2 * i + 1));
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
