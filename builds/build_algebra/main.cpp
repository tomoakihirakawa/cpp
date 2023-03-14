#include "basic_linear_systems.hpp"
#include "basic_vectors.hpp"

int main() {

   VV_d mat = {{1., 2., 3.},
               {4., 1., 6.},
               {7., 8., 1.}};

   std::vector<Tddd> mat3 = {{1., 2., 4.},
                             {2., 1., 6.},
                             {4., 4., 1.}};

   // タプルを含んだベクトルの積を実装するのは難しい

   VV_d M = {{1., 2., 3., 4.},
             {4., 1., 6., 4.},
             {7., 8., 1., 4},
             {7., 1., 1., 4},
             {1., 1., 1., 4},
             {1., 4., 1., 4},
             {0., 1., 1., 4}};

   std::vector<T4d> mat4 = {{1., 2., 4., 3},
                            {2., 1., 6., 3},
                            {4., 4., 1., 3},
                            {4., 4., 1., 3}};

   V_d vec = {1., 2., 3.};

   Print(Inverse(mat), red);
   Print(Dot(mat, vec), blue);
   Print(Dot(vec, mat), green);
   Print(Dot(mat, mat3), pink);
   Print(Dot(M, mat4), red);
   // Print(Dot(mat2, mat), blue);

   return 0;
};