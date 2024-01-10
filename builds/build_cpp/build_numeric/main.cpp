
#include "fundamental.hpp"

#include <iostream>
#include <algorithm>

/*
algorithmはvectorを修正するもの．
numericは，vectorの値を使ってある値を計算するもの．例えば和．
*/

int main()
{

  using V_d = std::vector<double>;
  using VV_d = std::vector<std::vector<double>>;

  TimeWatch watch;
  V_d v = {0, 1, 2, 3, 4, 5};
  V_d u = {0, 10, 20, 30, 40, 50};
  V_d long_vector = Subdivide(1., 10., 100000);
  std::cout << v << std::endl;

  //std::transform(v.begin(), v.end(), v.begin(), v.begin(), [](const double a, const double b){return a+b;});
  /*
  accumulateとreduceの速度について
  */
  std::cout << watch() << std::endl;
  double a;
  for (auto i = 0; i < 1000; ++i)
    a = std::accumulate(long_vector.begin(), long_vector.end(), 0.);
  std::cout << a << std::endl;
  std::cout << watch() << std::endl;
  for (auto i = 0; i < 1000; ++i)
    a = std::reduce(long_vector.begin(), long_vector.end(), 0.);
  std::cout << a << std::endl;
  std::cout << watch() << std::endl;
  /* ------------------------------------------------------ */
  std::cout << watch() << std::endl;
  std::cout << std::accumulate(long_vector.begin(), long_vector.end(), 1., [](auto a, auto b)
                               { return a * b; })
            << std::endl;
  std::cout << watch() << std::endl;
  std::cout << std::reduce(long_vector.begin(), long_vector.end(), 1., [](auto a, auto b)
                           { return a * b; })
            << std::endl;
  std::cout << watch() << std::endl;
  /* ------------------------------------------------------ */
  a = std::transform_reduce(
      long_vector.begin(), long_vector.end(), 0,
      [](auto a, auto b)
      { return a + b; }, // 集計関数
      [](auto a)
      { return 1.; });
  std::cout << Red << a << std::endl;
  std::cout << watch() << std::endl;
  /* ------------------------------------------------------ */
  std::cout << u << std::endl;
  std::cout << watch() << std::endl;
  std::cout << "transform_reduce" << std::endl;
  std::cout << "v = " << v << std::endl;
  std::cout << "u = " << u << std::endl;
  std::cout << std::transform_reduce(
                   v.begin(), v.end(), u.begin(), 0.,
                   [](const double acc, const double d)
                   { return acc + d; /*reduce*/ },
                   [](const double a, const double d)
                   { return a * d; /*transformed->reduce*/ })
            << std::endl;
  std::cout << "v = " << v << std::endl;
  std::cout << "u = " << u << std::endl;
  std::cout << watch() << std::endl;
  std::cout << "/* ------------------------------------------------------ */" << std::endl;
  std::cout << std::transform_reduce(
                   v.begin(), v.end(), u.begin(), 0.,
                   [](const double acc, const double d)
                   { return acc * d; },
                   [](const double a, const double d)
                   { return a + d; })
            << std::endl;
  std::cout << "v = " << v << std::endl;
  std::cout << "u = " << u << std::endl;
  std::cout << watch() << std::endl;
  std::cout << "/* ------------------------------------------------------ */" << std::endl;
  std::transform(
      v.begin(), v.end(), u.begin(), u.begin(), [](const double acc, const double d)
      { return acc * d; });
  std::cout << "v = " << v << std::endl;
  std::cout << "u = " << u << std::endl;
  std::cout << watch() << std::endl;
  std::cout << "/* ------------------------------------------------------ */" << std::endl;
  std::transform(v.begin(),
                 v.end(),
                 u.begin(),
                 [](const auto &a)
                 { return a; });

  std::cout << u << std::endl;
}
