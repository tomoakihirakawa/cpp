#if defined(test1)
// 2021/07/11

#include <iostream>
#include <vector>
#include <tuple>
#include <string>

int main()
{

  std::tuple<int, double, std::string> t = std::make_tuple(1, 'a', "hello");

  std::cout << std::get<0>(t) << std::endl;
  std::cout << std::get<1>(t) << std::endl;
  std::cout << std::get<2>(t) << std::endl;

  /* ------------------------------------------------------ */
  // std::make_tuple()はほとんどの状況で必要ない.
  std::tuple<int, double, std::string> tup{1, 'a', "hello"};

  std::cout << std::get<0>(tup) << std::endl;
  std::cout << std::get<1>(tup) << std::endl;
  std::cout << std::get<2>(tup) << std::endl;

  /* ------------------------------------------------------ */

  auto [i, s0, s1] = tup;
  std::cout << i << std::endl;
  std::cout << s0 << std::endl;
  std::cout << s1 << std::endl;
};

#else

#include <iostream>
#include <vector>
#include <tuple>
#include <string>
#include "fundamental.hpp"

int main()
{

  int count = 1000000;
  /* ------------------------------------------------------ */
  TimeWatch watch;
  for (auto c = 0; c < 10; ++c)
  {
    /* ------------------------------------------------------ */
    std::cout << watch.get() << std::endl;
    for (auto i = 0; i < count; ++i)
    {
      V_d a = {1., 2., 3.};
      Dot(a, a);
      Normalize(a);
      Cross(a, a);
    }
    /* ------------------------------------------------------ */
    std::cout << "vector " << watch.get() << std::endl;
    for (auto i = 0; i < count; ++i)
    {
      // Tddd A = std::make_tuple(1., 2., 3.);
      Dot({1., 2., 3.}, {1., 2., 3.});
      Normalize({1., 2., 3.});
      Cross({1., 2., 3.}, {1., 2., 3.});
    }
    std::cout << "tuple " << watch.get() << std::endl;
  }
};

#endif