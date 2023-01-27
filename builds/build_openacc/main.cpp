//
//  add two matrixs and store the result in another matrix
//
//  (c)Copyright Spacesoft corp., 2018 All rights reserved.
//                                Hiro KITAYAMA
//
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "basic.hpp"

struct dataset {
   double x, y;
   double v;
   dataset(const double &X, const double &Y) : x(X), y(Y){};
};

// main
int main(int argc, char *argv[]) {
   float **a, **b, **c;
   clock_t start, stop;
   int i, j, n = 10000;

   a = (float **)malloc(sizeof(float *) * n);
   b = (float **)malloc(sizeof(float *) * n);
   c = (float **)malloc(sizeof(float *) * n);
   for (i = 0; i < n; i++) {
      a[i] = (float *)malloc(sizeof(int) * n);
      b[i] = (float *)malloc(sizeof(int) * n);
      c[i] = (float *)malloc(sizeof(int) * n);
   }

   // initialize array
   for (j = 0; j < n; j++) {
      for (i = 0; i < n; i++) {
         a[j][i] = (float)(i + 1000);
         b[j][i] = (float)i / 10.f;
      }
   }

   start = clock();

   // calc.
   for (j = 0; j < n; j++) {
      for (i = 0; i < n; i++) {
         c[j][i] = a[j][i] + b[j][i];
      }
   }

   stop = clock();

   fprintf(stdout, "      C: ");
   fprintf(stdout, "elapsed time = %.20f [sec]\n",
           (float)(stop - start) / CLOCKS_PER_SEC);

   start = clock();

#define example1
#if defined(example0)
// #pragma acc kernels
#pragma acc parallel
   for (j = 0; j < n; j++) {
      for (i = 0; i < n; i++) {
         c[j][i] = a[j][i] + b[j][i];
      }
   }
#elif defined(example1)
#pragma acc enter data copyin(a, b, c)
#pragma acc parallel loop present(a, b, c)
   for (j = 0; j < n; ++j) {
      for (i = 0; i < n; ++i) {
         c[j][i] = a[j][i] + b[j][i];
      }
   }
#pragma acc exit data copyout(a, b, c)
#endif

   stop = clock();

   fprintf(stdout, "OpenACC: ");
   fprintf(stdout, "elapsed time = %.20f [sec]\n",
           (float)(stop - start) / CLOCKS_PER_SEC);

   for (i = 0; i < n; i++) {
      free(a[i]);
      free(b[i]);
      free(c[i]);
   }
   free(a);
   free(b);
   free(c);

   /* -------------------------------------------------------------------------- */
   {

      int N = 1000;
      std::vector<dataset *> vec;
      vec.reserve(N * N);
      for (auto i = 0; i < N; i++)
         for (auto j = 0; j < N; j++)
            vec.emplace_back(new dataset((double)i, (double)j));

      start = clock();
      for (const auto &v : vec)
         v->v = v->x * v->y;
      stop = clock();

      fprintf(stdout, "      C: ");
      fprintf(stdout, "elapsed time = %.20f [sec]\n",
              (float)(stop - start) / CLOCKS_PER_SEC);

#pragma acc enter data copyin(vec)
      start = clock();
#pragma acc parallel loop present(vec)
      for (auto i = 0; i < N * N; i++) {
         auto v = vec[i];
         v->v = v->x * v->y;
      }
#pragma acc exit data copyout(vec)
      stop = clock();

      fprintf(stdout, "OpenACC: ");
      fprintf(stdout, "elapsed time = %.20f [sec]\n",
              (float)(stop - start) / CLOCKS_PER_SEC);
   }
   /* -------------------------------------------------------------------------- */

   {
      int N = 1000;
      std::vector<dataset *> vec, wec;
      vec.reserve(N * N);
      for (auto i = 0; i < N; i++)
         for (auto j = 0; j < N; j++) {
            auto a = new dataset((double)i, (double)j);
            vec.emplace_back(a);
            wec.emplace_back(a);
         }
      auto func = [&](const auto &v) {
         v->v = v->x * v->y;
      };

      // auto check = [&]() {
      //    double total = 0;
      //    for (const auto &v : vec)
      //       total += std::abs(v->v);
      //    return total;
      // };

      auto initialize = [&]() {
         for (const auto &v : vec)
            v->v = 0;
         // std::cout << "initialized check = " << check() << std::endl;
      };

      // initialize();
      start = clock();
      for (const auto &v : vec)
         func(v);
      stop = clock();

      // std::cout << check() << std::endl;
      fprintf(stdout, "      C: ");
      fprintf(stdout, "elapsed time = %.20f [sec]\n",
              (float)(stop - start) / CLOCKS_PER_SEC);

      // initialize();
#pragma acc enter data copyin(vec)
      start = clock();

#pragma acc parallel loop present(vec)
      for (auto i = 0; i < N * N; i++) {
         func(vec[i]);
      }
#pragma acc exit data copyout(vec)
      stop = clock();

      // std::cout << check() << std::endl;
      fprintf(stdout, "OpenACC: ");
      fprintf(stdout, "elapsed time = %.20f [sec]\n",
              (float)(stop - start) / CLOCKS_PER_SEC);
   }

   return 0;
}
