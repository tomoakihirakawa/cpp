#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#define check_for_each

int main() {
#if defined(check_for_each)

   int count = 10000;
   double ans = 0;
   Tddd u = {1., 1., 1.};
   T8Tddd tup = {
       {1., 2., 3.},
       {2., 3., 4.},
       {3., 4., 5.},
       {4., 5., 6.},
       {5., 6., 7.},
       {6., 7., 8.},
       {7., 8., 9.},
       {8., 9., 10.},
   };
   /* ------------------------------------------------------ */
   /*                        for_each                        */
   /* ------------------------------------------------------ */
   Timer timer;
   for (auto i = 0; i < count; ++i)
      for_each(tup, [&](const auto &X) { ans += Norm(Cross(X, u)); });
   std::cout << "ans=" << ans << std::endl;
   std::cout << "for_each " << timer() << std::endl;
   /* ------------------------------------------------------ */
   /*                         simple                         */
   /* ------------------------------------------------------ */
   ans = 0;
   for (auto i = 0; i < count; ++i) {
      auto [x0, x1, x2, x3, x4, x5, x6, x7] = tup;
      ans += Norm(Cross(x0, u));
      ans += Norm(Cross(x1, u));
      ans += Norm(Cross(x2, u));
      ans += Norm(Cross(x3, u));
      ans += Norm(Cross(x4, u));
      ans += Norm(Cross(x5, u));
      ans += Norm(Cross(x6, u));
      ans += Norm(Cross(x7, u));
   }
   std::cout << "ans=" << ans << std::endl;
   std::cout << "simple " << timer() << std::endl;
   /* ------------------------------------------------------ */
   /*                         fatest                         */
   /* ------------------------------------------------------ */
   ans = 0;
   for (auto i = 0; i < count; ++i) {
      ans += Norm(Cross(std::get<0>(tup), u));
      ans += Norm(Cross(std::get<1>(tup), u));
      ans += Norm(Cross(std::get<2>(tup), u));
      ans += Norm(Cross(std::get<3>(tup), u));
      ans += Norm(Cross(std::get<4>(tup), u));
      ans += Norm(Cross(std::get<5>(tup), u));
      ans += Norm(Cross(std::get<6>(tup), u));
      ans += Norm(Cross(std::get<7>(tup), u));
   }
   std::cout << "ans=" << ans << std::endl;
   std::cout << "fastest " << timer() << std::endl;
   /* ------------------------------------------------------ */
   std::cout << "for_eachもそれなりに早い" << std::endl;

#elif defined(check_speed)

   int count = 10000000;
   Timer timer;
   V_d vec = {1., 2., 3.};
   Tddd tup = {1., 2., 3.};
   /* ------------------------------------------------------ */
   /*                        Normalie                        */
   /* ------------------------------------------------------ */
   for (auto i = 0; i < count; ++i)
      Normalize(vec);
   std::cout << "vector " << timer() << std::endl;
   /* ------------------------------------------------------ */
   for (auto i = 0; i < count; ++i)
      Normalize(tup);
   std::cout << "tuple " << timer() << std::endl;
   /* ------------------------------------------------------ */
   /*                        Dot                            */
   /* ------------------------------------------------------ */
   for (auto i = 0; i < count; ++i)
      Dot(vec, vec);
   std::cout << "vector " << timer() << std::endl;
   /* ------------------------------------------------------ */
   for (auto i = 0; i < count; ++i)
      Dot(tup, tup);
   std::cout << "tuple " << timer() << std::endl;
   /* ------------------------------------------------------ */
   /*                        Cross                          */
   /* ------------------------------------------------------ */
   for (auto i = 0; i < count; ++i)
      Cross(vec, vec);
   std::cout << "vector " << timer() << std::endl;
   /* ------------------------------------------------------ */
   for (auto i = 0; i < count; ++i)
      Cross(tup, tup);
   std::cout << "tuple " << timer() << std::endl;
   /* ------------------------------------------------------ */
   std::cout << "tupleが圧倒的に速い" << std::endl;

#elif defined(basic)

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

#endif
};