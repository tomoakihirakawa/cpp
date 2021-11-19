#include "svd.hpp"
#include "fundamental.hpp"

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
int main()
{
  MatDoub A = {{2, 3},
               {3, 4},
               {7, 1}};

  VecDoub b = {4,5,2};

  VecDoub solution = {1., 1.};

  SVD svd(A);
  svd.solve(b, solution);

  for (auto &a : A)
    Print(Dot(a, solution), blue);

  Print(solution, Red);
  Print(svd.u);
  Print(svd.w);
  Print(svd.v);
}
