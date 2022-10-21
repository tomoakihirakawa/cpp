#include "GNUPLOT.hpp"
#include "Network.hpp"

int main()
{

  GNUPLOT plot;
  //    plot.Set({{"key",""},{"xrange","[-100:100]"},{"yrange","[-200:200]"}});
  // NetworkObj obj("./obj/camel.obj");
  // obj.scale(2.);
  V_d cardinal = {1000., 100., .2};
  NetworkObj obj0("./obj/boundary.obj", "Neumann:structure");
  obj0.translate({-80., -50., 5.});
  //    plot.SaveVectorData(getVectorData(obj.Faces),{{"arrowstyle","1"},{"title","obj0"}});
  plot.SaveVectorData(getVectorData(obj0.Faces), {{"arrowstyle", "2"}, {"title", "obj1"}});
  plot.plot3d();
  plot.Clear();
  std::cin.ignore();

  auto v = Subdivide(300., 1., 100);
  int counter = 0;
  //////////////////////////
  for (auto len : v)
  {
    Print(len, red);

#define with_lap

#if defined(with_lap)
    mk_vtu("./vtu/flip_lap" + std::to_string(counter++) + ".vtu", obj0.Faces);
#else
    mk_vtu("./vtu/flip" + std::to_string(counter++) + ".vtu", obj0.Faces);
#endif
    // for (const auto &l : obj0.getLines())
    //   if (l->length() > len)
    //     l->divide();

    // for (auto i = 0; i < 10; i++)
    // {
    //   for (auto l : obj0.getLines())
    //   {
    //     auto f = l->getFaces();
    //     if (f.size() == 2)
    //       if (Norm(Cross(f[0]->getNormal(), f[1]->getNormal())) < 1E-5)
    //         l->flipIfIllegal();
    //   }
    // }

    int count = 0;
    bool found;
    do
    {
      found = false;
      auto tmp = obj0.getLines();
      for (const auto &l : tmp)
      {
        auto f = l->getFaces();
        if (f.size() == 2)
        {
          if (Norm(Cross(f[0]->getNormal(), f[1]->getNormal())) < 1E-2)
          {
            found = l->flipIfIllegal();
          }
        }
        if (l->length() > len){
          l->divide();
          found = true;
          break;
        }
      }
    } while (found && count++ < 10000);

#if defined(with_lap)
    for(auto i=0; i<10; i++)  
      LaplacianSmoothingIfFlat(obj0.getPoints());
#endif

    // obj0.getLines()[0]->flip();
    //    plot.SaveVectorData(getVectorData(obj.Faces),{{"arrowstyle","1"},{"title","obj0"}});
    plot.SaveVectorData(getVectorData(obj0.Faces), {{"arrowstyle", "2"}, {"title", "obj1"}});
    plot.plot3d();
    plot.Clear();
  }
  ////////////////////////////

  Print(Green + "Enterを押して終了");
  std::cin.ignore();
}
