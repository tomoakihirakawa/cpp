#include "GNUPLOT.hpp"
#include "fundamental.hpp"
#include "Network.hpp"

int main()
{

  NetworkX water({2, 2}, {30.01, 30.01, 1 / 2.}, 2);

  GNUPLOT plot;
  plot.Set({{"key", ""}});
  plot.SaveVectorData(getVectorData(water.Faces));

  int ind = 20;
  auto p = water.Points[ind];
  plot.SaveData({p->getX()}, {{"w", "lp"}});
  plot.SaveVectorData(getVectorData(p->getNeighbors()), {{"lc", "'blue'"}, {"title", "getNeighbors()"}});
  plot.SaveVectorData(getVectorData(p->getNeighborsSort()), {{"lc", "'purple'"}, {"title", "getNeighborsSort()"}});

  plot.plot3d();

  std::cin.ignore();
};
