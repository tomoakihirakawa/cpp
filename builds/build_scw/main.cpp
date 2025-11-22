#include "GNUPLOT.hpp"
#include "fundamental.hpp"
#include "SCW.hpp"
int main()
{
  std::string dir0 = "/Users/tomoaki/Dropbox/research/scw_3d/workshop/fortran/sym/sym_direct/";
  std::string dir1 = "results/40_6/2d0000EP00_1d0000EP01/1d0000000EM01/4d5000000EP01_2d0828457EP00/";
  std::string mode_filename = dir0 + dir1 + "mode.dat";
  std::string data_filename = dir0 + dir1 + "data.dat";
  SCW scw(data_filename, mode_filename);

  GNUPLOT plot;
  /////////////
  VV_d a(scw.a.size());
  for (int i = 0; i < scw.a.size(); i++)
    for (int j = 0; j < scw.a.size(); j++)
    {
      auto tmp = std::abs(scw.a[i][j]);
      a[i].push_back(log10(tmp) < -14 ? -14 : log10(tmp));
    }

  plot.Set({{"palette", "define ( 0 'white', 1 '#fec287', 2 '#fb8761', 3 '#e55964', 4 '#b5367a', 5 '#812581', 6 '#4f127b', 7 '#1c1044', 8 '#000004' )"}});
  plot.SaveMatrixData(a, {{"w", "image"}});
  plot.MatrixPlot();
  std::cin.ignore();
  plot.Clear();
  /////////////
  VVV_d eta;
  plot.Set({{"contour", "base"}, {"title", "'h=" + std::to_string(-scw.h) + "'"}, {"zrange", "[" + std::to_string(-scw.h) + ":]"}});
  for (auto x : linspace({-scw.LL, scw.LL}, 50))
  {
    VV_d row(0);
    for (auto y : linspace({-scw.LL, scw.LL}, 50))
      row.emplace_back(V_d{x, y, scw.eta(x, y)});
    eta.emplace_back(row);
  }
  plot.SaveSplotData(eta, {{"w", "l"}});
  plot.plot3d();
  std::cin.ignore();
  plot.Clear();
  /////////////
  VVV_d ETA;
  plot.Set({{"contour", "base"}, {"pm3d", "at ss"}, {"title", "'D=" + std::to_string(-scw.D) + "'"}, {"zrange", "[" + std::to_string(-scw.D) + ":]"}});
  for (auto x : linspace({-M_PI, M_PI}, 50))
  {
    VV_d row(0);
    for (auto y : linspace({-M_PI, M_PI}, 50))
      row.emplace_back(V_d{x, y, scw.ETA(x, y)});
    ETA.emplace_back(row);
  }
  plot.SaveSplotData(ETA, {{"w", "pm3d"}});
  plot.plot3d();
  std::cin.ignore();
  plot.Clear();
  /////////////
  VVV_d PHI;
  plot.Set({{"contour", "base"}, {"pm3d", "at ss"}});
  for (auto x : linspace({-M_PI, M_PI}, 50))
  {
    VV_d row(0);
    for (auto y : linspace({-M_PI, M_PI}, 50))
      row.emplace_back(V_d{x, y, scw.PHI_without_mean_velocity(x, y, scw.ETA(x, y))});
    PHI.emplace_back(row);
  }
  plot.SaveSplotData(PHI, {{"w", "pm3d"}});
  plot.plot3d();
  std::cin.ignore();
  plot.Clear();
  /////////////
  VVV_d phi;
  plot.Set({{"contour", "base"}, {"pm3d", "at ss"}});
  for (auto x : linspace({-scw.LL, scw.LL}, 50))
  {
    VV_d row(0);
    for (auto y : linspace({-scw.LL, scw.LL}, 50))
      row.emplace_back(V_d{x, y, scw.phi_without_mean_velocity(x, y, scw.eta(x, y))});
    phi.emplace_back(row);
  }
  plot.SaveSplotData(phi, {{"w", "pm3d"}});
  plot.plot3d();
  std::cin.ignore();
  plot.Clear();
}