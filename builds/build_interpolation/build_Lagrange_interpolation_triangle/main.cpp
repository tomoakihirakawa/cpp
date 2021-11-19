#include "GNUPLOT.hpp"
#include "fundamental.hpp"

template <typename T>
struct interpolationCommon
{
  interpolationCommon(){};
};

template <typename T>
struct interpolation : public interpolationCommon<T>
{
  T sample;
  interpolation(const T &sample_IN) : sample(sample_IN), interpolationCommon<T>(){};
};

template <>
struct interpolation<T6Tddd> : public interpolationCommon<T6Tddd>
{
  T6Tddd sample;
  interpolation<T6Tddd>(const T6Tddd &sample_IN) : sample(sample_IN), interpolationCommon<T6Tddd>(){};
  Tddd operator()(const double &t0, const double &t1)
  {
    double t2 = 1 - t0 - t1;
    //    0
    //  3   5
    // 1  4  2
    return Dot(T6d{t0 * (2. * t0 - 1.) /*0*/,
                   t1 * (2. * t1 - 1.) /*1*/,
                   t2 * (2. * t2 - 1.) /*2*/,
                   4. * t0 * t1 /*3*/,
                   4. * t1 * t2 /*4*/,
                   4. * t0 * t2 /*5*/},
               sample);
    //Dot(vector,matrix)の場合，matrixの各行（座標）にvectorの成分が掛かる
    //vector成分は各座標に掛かる
  };
  T2Tddd div(const double &t0, const double &t1)
  {
    double t2 = 1 - t0 - t1;
    double dt2dt0 = -1;
    double dt2dt1 = -1;
    //    0
    //  3   5
    // 1  4  2
    T6d dNdt0 = {(2. * t0 - 1.) + t0 * (2) /*0*/,
                 0. /*1*/,
                 dt2dt0 * (2. * t2 - 1.) + t2 * (2. * dt2dt0) /*2*/,
                 4. * t1 /*3*/,
                 4. * t1 * dt2dt0 /*4*/,
                 4. * t2 + 4. * t0 * dt2dt0 /*5*/};
    T6d dNdt1 = {0. /*0*/,
                 (2. * t1 - 1.) + t1 * (2.) /*1*/,
                 dt2dt1 * (2. * t2 - 1.) + t2 * (2. * dt2dt1) /*2*/,
                 4. * t0 /*3*/,
                 4. * t2 + 4. * t1 * dt2dt1 /*4*/,
                 4. * t0 * dt2dt1 /*5*/};
    return {Dot(dNdt0, sample), Dot(dNdt1, sample)};
  };
};
/* ------------------------------------------------------ */
template <>
struct interpolation<T3Tddd> : public interpolationCommon<T3Tddd>
{
  T3Tddd sample;
  interpolation<T3Tddd>(const T3Tddd &sample_IN) : sample(sample_IN), interpolationCommon<T3Tddd>(){};
  Tddd operator()(const double &t0, const double &t1)
  {
    double t2 = 1 - t0 - t1;
    return Dot(Tddd{t0, t1, 1 - t0 - t1}, sample);
  };
  T2Tddd div(const double &t0, const double &t1)
  {
    //   0
    // 1   2
    Tddd dNdt0 = {1, 0, -1};
    Tddd dNdt1 = {0, 1, -1};
    return {Dot(dNdt0, sample), Dot(dNdt1, sample)};
  };
};
/* ------------------------------------------------------ */
int main()
{
  GNUPLOT plot;
  int s = 50;

  // T6Tddd p = {{0, 2, 1}, {-2, -1, -.5}, {2, -1, 0}, {-1., 0, .1}, {0, -1.5, 0}, {1, 0, -.5}};

  // double max = 227;
  // T6Tddd c = {{max, 0, 0}, {0, max, 0}, {0, 0, max}, {max, max, 0}, {0, max, max}, {max, 0, max}};
  //
  T3Tddd p = {{0, 2, 1}, {-2, -1, -.5}, {1, 0, -.5}};

  double max = 227;
  T3Tddd c = {{max, 0, 0}, {0, max, 0}, {max, 0, max}};

  interpolation interpLoc(p);
  interpolation interpRGB(c);

  // for (auto i = 0; i < s; i++)
  //   for (auto j = 0; j < s - i; j++)
  //     plot.SaveData(interpLoc(i / (s - 1.), j / (s - 1.)), {{"ps", "1"}, {"pt", "7"}, {"w", "p"}, {"lc", plot.rgb(interpRGB(i / (s - 1.), j / (s - 1.)))}});
  // plot.Plot3D();
  // std::cin.ignore();
  ////////////
  VVV_d data;
  for (auto i = 0; i < s; i++)
    for (auto j = 0; j < s - i; j++)
    {
      double x = i / (s - 1.), y = j / (s - 1.);
      // x /= 2.;
      // y /= 2.;
      auto intp = interpLoc(x, y);
      auto ddt0_ddt1 = interpLoc.div(x, y);
      auto cross = Cross(std::get<0>(ddt0_ddt1), std::get<1>(ddt0_ddt1));
      cross = cross / 30.;
      data.push_back({ToVector(intp), ToVector(cross), {Norm(cross)}});
    }
  plot.VectorPlot3D(data, {{"lc", "palette"}});
  std::cin.ignore();
};