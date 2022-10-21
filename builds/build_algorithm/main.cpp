
#include "fundamental_vectors.hpp"

#include <iostream>
#include <algorithm>

struct object
{
  int x;
  object(int xin) : x(xin){};
};

int main()
{

  using V_d = std::vector<double>;
  using VV_d = std::vector<std::vector<double>>;

  V_d v = {0, 1, 2, 3, 4, 5};
  V_d u = {0, 10, 20, 30, 40, 50};

  std::cout << v << std::endl;

  //std::transform(v.begin(), v.end(), v.begin(), v.begin(), [](const double a, const double b){return a+b;});

  std::cout << u << std::endl;

  std::cout << std::transform_reduce(
                   v.begin(), v.end(), u.begin(), 0.,
                   [](const double acc, const double d)
                   { return acc + d; },
                   [](const double a, const double d)
                   { return a * d; })
            << std::endl;

  std::transform(v.begin(),
                 v.end(),
                 u.begin(),
                 [](const auto &a)
                 { return a; });

  std::cout << u << std::endl;

  /* -------------------- min maxのチェック -------------------- */

  std::vector<object *> Vo;
  for (int i = 0; i < 10; ++i)
    Vo.push_back(new object(i));

  Vo.push_back(new object(0));
  Vo.push_back(new object(0));
  Vo.push_back(new object(9));

  std::cout << (*std::max_element(Vo.begin(), Vo.end(), [](auto a, auto b)
                                  { return a->x < b->x; }))
                   ->x
            << std::endl;

  std::cout << (*std::min_element(Vo.begin(), Vo.end(), [](auto a, auto b)
                                  { return a->x < b->x; }))
                   ->x
            << std::endl;

  for (int i = 0; i < Vo.size(); ++i)
    delete Vo[i];
}
