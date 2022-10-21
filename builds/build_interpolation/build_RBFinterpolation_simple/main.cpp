#include "fundamental.hpp"
//#include "lu.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

#include "GNUPLOT.hpp"
#include "InterpolationRBF.hpp"

double func(double x, double y)
{
  double r = sqrt(x * x + y * y);
  return sin(r) * cos(r);
};

V_d grad_func(double x, double y)
{
  double r = sqrt(x * x + y * y);
  return {x * cos(2 * r) / r,
          y * cos(2 * r) / r};
};

double laplacian_func(double x, double y)
{
  double r = sqrt(x * x + y * y);
  return (2 * cos(2 * r)) / r - 2 * sin(2 * r);
};
//@ ------------------------------------------------------ */
JSON settingRBF("./settingRBF.json");
double scale = stod(settingRBF["scale"])[0];
double scale_additional = stod(settingRBF["scale_additional"])[0];
double eps = stod(settingRBF["eps"])[0];
int max = stoi(settingRBF["max"])[0];
int improve = stoi(settingRBF["improve"])[0];
//@ ------------------------------------------------------ */
int main()
{
  using V_d = std::vector<double>;
  VV_d xyz = {};
  VV_d XYZ = {};

  //サンプルデータの作成

  VV_d XY;
  V_d Res;
  VV_d N;

  // for (auto x : Subdivide(-10, 10, 31))
  //   for (auto y : Subdivide(-10, 10, 31))
  //   {
  //     std::cout << x << ", " << y << std::endl;
  //     double r = sqrt(x * x + y * y);
  //     if (!r < 1E-10)
  //     {
  //       XY.emplace_back(V_d{x, y});
  //       Res.emplace_back(0.);
  //       V_d tmp = grad_func(x, y);
  //       N.emplace_back(V_d{1 / tmp[0], -1 / tmp[1]});
  //     }
  //   }

  std::vector<std::tuple<Tddd, double>> eq;
  std::vector<std::tuple<Tddd, double>> kernel_locations;

  for (auto i = 0; i < max; i++)
  {
    for (auto j = 0; j < max; j++)
    {
      double len = (max - 1.);
      double x = 5. * M_PI * (i - len / 2.) / len;
      double y = 5. * M_PI * (j - len / 2.) / len;
      if ((i + j + 1) % 2 == 0)
      {
        xyz.push_back({0, 0, 0});
        auto &V = (*xyz.rbegin());
        V[0] = x;
        V[1] = y;
        V[2] = func(x, y);
        /* ------------------------------------------------------ */
        kernel_locations.emplace_back(std::tuple<Tddd, double>{Tddd{x, y, 0}, scale});
        eq.emplace_back(std::tuple<Tddd, double>{Tddd{x, y, 0}, func(x, y)});
        /* ------------------------------------------------------ */
        if (improve)
        {
          V_d tmp = grad_func(x, y);
          if (isFinite(1 / tmp[0]) && isFinite(-1 / tmp[1]))
          {
            XY.emplace_back(V_d{x, y});
            Res.emplace_back(0.);
            Tddd X = {x, y, 0};
            {
              Tddd n = {1 / tmp[0], -1 / tmp[1], 0};
              Tddd m = Normalize(n) * eps;
              /* ------------------------------------------------------ */
              for (auto &n : std::vector<Tddd>{m, -m})
              {
                kernel_locations.emplace_back(std::tuple<Tddd, double>{X + n, scale_additional});
                eq.emplace_back(std::tuple<Tddd, double>{X + n, func(x, y) + Dot(grad_func(x, y), ToVector(n))});
                XYZ.emplace_back(ToVector(X + n));
              }
            }
            {
              Tddd n = {1 / tmp[0], 1 / tmp[1], 0};
              Tddd m = Normalize(n) * eps;
              /* ------------------------------------------------------ */
              for (auto &n : std::vector<Tddd>{m, -m})
              {
                kernel_locations.emplace_back(std::tuple<Tddd, double>{X + n, scale_additional});
                eq.emplace_back(std::tuple<Tddd, double>{X + n, func(x, y) + Dot(grad_func(x, y), ToVector(n))});
                XYZ.emplace_back(ToVector(X + n));
              }
            }
          }
        }
      }
      // else
      // {
      //   std::cout << "i + j + 1 = " << i + j + 1 << std::endl;
      //   V_d tmp = grad_func(x, y);
      //   if (isFinite(1 / tmp[0]) && isFinite(-1 / tmp[1]))
      //   {
      //     XY.emplace_back(V_d{x, y});
      //     XYZ.emplace_back(V_d{x, y, 0});
      //     Res.emplace_back(0.);
      //     N.emplace_back(V_d{1 / tmp[0], -1 / tmp[1]});
      //     /* ------------------------------------------------------ */
      //     kernel_locations.emplace_back(Tddd{x, y, 0});
      //     eq_der_dot.emplace_back(std::tuple<Tddd, Tddd, double, double>{Tddd{x, y, 0}, Tddd{1 / tmp[0], -1 / tmp[1], 0}, 0., scale_grad});
      //     /* ------------------------------------------------------ */
      //   }
      // }
    }
  }
  int c = 0;
  // for (auto &x : Subdivide(-10, 10, 15))
  //   for (auto &y : Subdivide(-10, 10, 15))
  //   {
  //     // if (c++ % 4 == 0)
  //     {
  //       V_d tmp = grad_func(x, y);
  //       if (isFinite(1 / tmp[0]) && isFinite(-1 / tmp[1]))
  //       {
  //         XY.emplace_back(V_d{x, y});
  //         XYZ.emplace_back(V_d{x, y, 0});
  //         Res.emplace_back(0.);
  //         N.emplace_back(V_d{1 / tmp[0], -1 / tmp[1]});
  //         /* ------------------------------------------------------ */
  //         kernel_locations.emplace_back(Tddd{x, y, 0});
  //         eq_der_dot.emplace_back(std::tuple<Tddd, Tddd, double, double>{Tddd{x, y, 0}, Tddd{1 / tmp[0], -1 / tmp[1], 0}, 0., scale_grad});
  //         /* ------------------------------------------------------ */
  //       }
  //     }
  //   }

  auto vx = Transpose(xyz)[0];
  auto vy = Transpose(xyz)[1];
  auto vz = Transpose(xyz)[2];

  //サンプルデータを基に補間
  try
  {
    GNUPLOT plt;
    plt.Set({{"key", ""},
             {"contour", "base"},
             {"xrange", "[-10.:10.]"},
             {"yrange", "[-10.:10.]"},
             {"zrange", "[-5.:5.]"}});

    for (auto l = 0; l < 1; l++)
    {
      int max = 60;
      VVV_d xyzRBF = {}, diff = {};
      // VVV_d xyzRBFmod={};
      VVV_d xyzExact = {};

      VVV_d lapRBF = {};
      VVV_d lapExact = {};

      VVV_d gradRBF = {};
      VVV_d gradExact = {};

      // InterpolationRBF interp(Transpose(VV_d{vx, vy}), vz,
      //                         XY,
      //                         Res,
      //                         N, 6.);
      InterpolationRBF_ interp(kernel_locations, eq);
      // InterpolationRBF interp(Transpose(VV_d{vx, vy}), vz, 0.1);
      // std::cout << "interp.scale = " << interp.scale << std::endl;
      // InterpolationIDW interp(Transpose(VV_d{vx, vy}), vz, 30.);

      for (auto i = 0; i < max; i++)
      {
        // VV_d xyzRBFmod_row={};
        VV_d lapRBF_row = {};
        VV_d lapExact_row = {};
        VV_d xyzRBF_row = {}, diff_row = {};
        VV_d xyzExact_row = {};

        for (auto j = 0; j < max; j++)
        {
          double len = (max - 1.);
          double x = 5. * M_PI * (i - len / 2.) / len;
          double y = 5. * M_PI * (j - len / 2.) / len;

          {
            //!シンプルな補間
            double z = interp({x, y, 0});
            xyzRBF_row.push_back({x, y, z});
            xyzExact_row.push_back({x, y, func(x, y)});
          }

          {
            double z = interp({x, y, 0});
            gradRBF.push_back({{x, y, 0}, ToVector(interp.grad({x, y, 0}))});
            gradExact.push_back({{x, y, 0}, grad_func(x, y)});
          }

          {
            //!ラプラシアンな補間
            double z = interp({x, y, 0});
            lapRBF_row.push_back({x, y, interp.laplacian({x, y, 0})});
            lapExact_row.push_back({x, y, laplacian_func(x, y)});
          }

          // lapRBF.push_back({x, y, interp.laplacian({x, y})});
          // lapExact_row.push_back({x, y, laplacian_func(x, y, 0.)});

          // diff_row.push_back({x, y, std::abs(z - func(x, y))});
        }

        // laplapcian.push_back(lapRBF);
        // lapExact.push_back(lapExact_row);
        // diff.push_back(diff_row);

        xyzRBF.push_back(xyzRBF_row);
        xyzExact.push_back(xyzExact_row);

        lapRBF.push_back(lapRBF_row);
        lapExact.push_back(lapExact_row);
      }

      // plt.SaveSplotData(XY, {{"w", "l"}, {"lc", "'blue'"}, {"title", "xyzRBF"}});
      plt.SaveSplotData(xyzRBF, {{"w", "l"}, {"lc", "'blue'"}, {"title", "xyzRBF"}});
      // plt.SaveSplotData(xyzRBFmod,{{"w","l"},{"lc","'green'"},{"title",std::to_string(RBFscale(xyz)*(l+1.))}});
      plt.SaveSplotData(xyzExact, {{"w", "l"}, {"lc", "'magenta'"}, {"title", "exact" + std::to_string(l)}});
      // plt.SaveSplotData(diff, {{"lc", "'purple'"}, {"title", "diff"}});
      plt.SaveVectorData(gradRBF, {{"lc", "'web-blue'"}, {"title", std::to_string(l)}});
      plt.SaveVectorData(gradExact, {{"lc", "'red'"}, {"title", "exact" + std::to_string(l)}});

      //!ラプラシアン
      std::cout << lapRBF << std::endl;
      std::cout << lapExact << std::endl;
      plt.SaveSplotData(lapRBF, {{"w", "l"}, {"lc", "'purple'"}, {"title", "laplapcian"}});
      plt.SaveSplotData(lapExact, {{"w", "l"}, {"lc", "'blue'"}, {"title", "laplapcian exact"}});
    }

    plt.SaveData(xyz, {{"pt", "7"}, {"w", "p"}, {"lc", "'blue'"}});
    plt.SaveData(XYZ, {{"pt", "7"}, {"w", "p"}, {"lc", "'green'"}});

    plt.plot3d();
    std::cin.ignore();
  }
  catch (error_message e)
  {
    std::cout << e.what() << reset << std::endl;
  };
}
