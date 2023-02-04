//
//  multiplys two matrixs  and store the result in another matrix
//
//  (c)Copyright Spacesoft corp., 2018 All rights reserved.
//                                Hiro KITAYAMA
//
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// main
int main(int argc, char *argv[]) {
   float *a, *b, *c, *hc;
   clock_t start, stop;
   int i, j, k, n = 256 * 2 * 2;

   if (argc > 1)
      n = atoi(argv[1]);

   fprintf(stdout, "matrix size = %d x %d\n", n, n);

   a = (float *)malloc(sizeof(float) * n * n);
   b = (float *)malloc(sizeof(float) * n * n);
   c = (float *)malloc(sizeof(float) * n * n);
   hc = (float *)malloc(sizeof(float) * n * n);

   // initialize array
   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
         a[i * n + j] = (float)(rand() / 4096);
         b[i * n + j] = (float)(rand() / 4096);
      }
   }

   start = clock();

   // calc.
   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
         float hcc = 0.0f;
         for (k = 0; k < n; k++) {
            hcc += a[i * n + k] * b[k * n + j];
         }
         hc[i * n + j] = hcc;
      }
   }

   stop = clock();

   fprintf(stdout, "      C: ");
   fprintf(stdout, "elapsed time = %.20f [sec]\n",
           (float)(stop - start) / CLOCKS_PER_SEC);

   start = clock();

// calc.
#pragma acc data copyout(c[:n * n]) copyin(b[:n * n], a[:n * n])
#pragma acc kernels
#pragma acc loop independent
   for (i = 0; i < n; i++) {
#pragma acc loop independent
      for (j = 0; j < n; j++) {
         float cc = 0.0f;
#pragma acc loop reduction(+ \
                           : cc)
         for (k = 0; k < n; k++) {
            cc += a[i * n + k] * b[k * n + j];
         }
         c[i * n + j] = cc;
      }
   }

   stop = clock();

   fprintf(stdout, "OpenACC: ");
   fprintf(stdout, "elapsed time = %.20f [sec]\n",
           (float)(stop - start) / CLOCKS_PER_SEC);

   // for (i = 0; i < n; i++) {
   //    for (j = 0; j < n; j++) {
   //       if (c[i][j] != hc[i][j]) {
   //          fprintf(stderr, "error!\n");
   //          break;
   //       }
   //    }
   // }

   // for (i = 0; i < n; i++) {
   //    free(a[i]);
   //    free(b[i]);
   //    free(c[i]);
   //    free(hc[i]);
   // }
   // free(a);
   // free(b);
   // free(c);
   // free(hc);

   return 0;
}
