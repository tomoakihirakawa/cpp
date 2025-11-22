//#define cehck_setCrossInfos

int step = 0;

#include "GNUPLOT.hpp"

// GNUPLOT plot_divide;
// #define PLOT_DIVIDE

#define vtk

#include "Network.hpp"

//#define DEBUG

GNUPLOT plot;
#include "bem.hpp"

using V_i = std::vector<int>;
using V_d = std::vector<double>;
using VV_d = std::vector<std::vector<double>>;
using VVV_d = std::vector<std::vector<std::vector<double>>>;

using V_Netp = std::vector<Network *>;
using V_netFp = std::vector<networkFace *>;
using VV_netFp = std::vector<V_netFp>;

std::vector<std::vector<double>> getLocations(const std::vector<networkPoint *> points)
{
  std::vector<std::vector<double>> ret(points.size(), std::vector<double>(3));
  for (auto i = 0; i < points.size(); i++)
    ret[i] = points[i]->xyz;
  return ret;
};
std::vector<std::vector<std::vector<double>>> getLocation(const std::vector<networkLine *> lines)
{
  std::vector<std::vector<std::vector<double>>> ret;
  for (auto i = 0; i < lines.size(); i++)
    ret[i] = lines[i]->getLocations();
  return ret;
};

networkLine *longerLine(networkLine *line_IN, double ratio = 1.01)
{
  networkLine *ret = line_IN;
  double len = line_IN->length(), v;

  for (const auto &p : line_IN->Points)
    for (const auto &l : p->Lines)
    {
      if (l != line_IN)
      { /*omit comparison with line_IN line self*/
        v = l->length();
        if (v > ratio * len)
        {
          len = v;
          ret = l;
        }
      }
    }
  return ret;
};

//searcher1
template <class T>
class searcher1 : public searcher<T>
{
public:
  searcher1() : searcher<T>(){};
  searcher1(T *start_, bool TorF_ = true) : searcher<T>(start_, TorF_){};
  bool condEnterLine(const T *p, const netLp l) override
  {
    if (!l->penetrateQ())
      return true;
    else
      return false;
  };
};
//searcher1
//searcher2
template <class T>
class searcher2 : public searcher<T>
{
public:
  searcher2() : searcher<T>(){};
  bool condGetObject(const netLp l, const T *P) override
  {
    if (!P->intersectQ())
      return true;
    else
      return false;
  };
};
//searcher2
//searcher3
//干渉したら探査を終了するsearcher
template <class T>
class searcher3 : public searcher<T>
{
public:
  searcher3() : searcher<T>(){};
  bool condEnterLine(const T *p, const netLp l) override
  {
    if (true)
      return true;
    else
      return false;
  };
  bool condKeepSearch(const netLp l, const T *P) override
  {
    if (!l->penetrateQ())
      return true;
    else
      return false;
  };
};
//searcher3
//searcher4
template <class T>
class searcher4 : public searcher<T>
{
public:
  bool didGetOne;
  searcher4() : searcher<T>(), didGetOne(false){};
  bool condEnterLine(const T *p, const netLp l) override
  {
    if (true)
    {
      this->didGetOne = false;
      return true;
    }
    else
      return false;
  };
  bool condGetObject(const netLp l, const T *P) override
  {
    if (!this->didGetOne)
    {
      this->didGetOne = true;
      return true;
    }
    else
      return false;
  };
  bool condKeepSearch(const netLp l, const T *P) override
  {
    if (this->didGetOne)
    {
      this->skipThisDepth = true;
      return true;
    }
    else
      return false;
  };
};
//searcher4
//searcher5
template <class T>
class searcher5 : public searcher<T>
{
public:
  bool didGetOne;
  searcher5() : searcher<T>(), didGetOne(false){};
  bool condEnterLine(const T *p, const netLp l) override
  {
    if (this->depth < 5)
    {
      return true;
    }
    else
      return false;
  };
};
//searcher5
//searcher6
//干渉部分のみをとる`searcher`
template <class T>
class searcher6 : public searcher<T>
{
public:
  searcher6() : searcher<T>(){};
  bool condGetObject(const netLp l, const T *P) override
  {
    if (P->intersectQ())
    {
      return true;
    }
    else
      return false;
  };
};
//searcher6
//searcher7
template <class T>
class searcher7 : public searcher<T>
{
  /*searcher7_detail
    干渉した線の先にあるobjectのみをとる`searcher`
    searcher7_detail*/
public:
  searcher7() : searcher<T>(){};
  bool condGetObject(const netLp l, const T *P) override
  {
    if (l->penetrateQ())
    {
      return true;
    }
    else
      return false;
  };
  bool condKeepSearch(const netLp l, const T *P) override
  {
    //内側から外に向かうために必要
    if (!l->penetrateQ())
      return true;
    else
      return false;
  };
};
//searcher7
//searcher8
template <class T>
class searcher8 : public searcher<T>
{
  /*searcher8_detail
    干渉した線の前までにあるobjectのみをとる`searcher`（`searcher7と被らない`）.
    searcher8_detail*/
public:
  searcher8() : searcher<T>(){};
  searcher8(T *obj) : searcher<T>(obj){};
  bool condGetObject(const netLp l, const T *P) override
  {
    if (!l->penetrateQ() && !P->intersectQ() /*Faceの場合に必要*/)
    {
      return true;
    }
    else
      return false;
  };
  bool condKeepSearch(const netLp l, const T *P) override
  {
    //内側から外に向かうために必要
    if (!l->penetrateQ())
      return true;
    else
      return false;
  };
};
//searcher8
//searcher9
// template <class T>
// class searcher9 : public searcher<T>
// {
//   /*searcher9_detail
//     内部から進み干渉点までを取得
//     networksに含まない物は取らない
//     内側の面を入手できないか?
//     searcher9_detail*/
// public:
//   searcher9() : searcher<T>(){};
//   searcher9(T *obj) : searcher<T>(obj){};
//   bool condEnterLine(const T *p, const netLp l) override
//   {
//     /*searcherに保存されているメンバーネットワークでなければならい*/
//     if (MemberQ(this->networks, (*l)(p)->getNetwork()))
//     {
//       /*searcherは到達しないが，貫通する線を通った点は，netObjs__に特別に保存する*/
//       if (l->penetrateQ())
//         if ((*l)(p) != nullptr && !(*l)(p)->intersectQ())
//           this->netObjs__.emplace_back((*l)(p));

//       /*干渉する線は通らない*/
//       if (!l->penetrateQ())
//       {
//         if ((*l)(p) != nullptr && (*l)(p)->intersectQ())
//           return true;

//         if (!p->intersectQ())
//           return true;
//       }
//     }
//     return false;
//   }
//   bool condGetObject(const netLp l, const T *P) override
//   {
//     if (P->intersectQ())
//       return true;
//     else
//       return false;
//   };
// };
//searcher9

class inside_point_searcher : public searcher<netP>
{
public:
  inside_point_searcher(netP *p) : searcher<netP>(p){};
  bool condEnterLine(const netP *P, const netLp l) override
  {
    if (!l->penetrateQ() && l->getNetwork() == P->getNetwork())
    {
      return true;
    }
    else
      return false;
  };
  bool condKeepSearch(const netLp l, const netP *P) override
  {
    if (!l->penetrateQ() && l->getNetwork() == P->getNetwork())
    {
      return true;
    }
    else
      return false;
  };
};
class end_point_searcher : public searcher<netP>
{
public:
  end_point_searcher(netP *p) : searcher<netP>(p){};
  bool condGetObject(const netLp l, const netP *P) override
  {
    if (l->penetrateQ() && l->getNetwork() == P->getNetwork())
    {
      return true;
    }
    else
      return false;
  };
  bool condKeepSearch(const netLp l, const netP *P) override
  {
    if (!l->penetrateQ() && l->getNetwork() == P->getNetwork())
    {
      return true;
    }
    else
      return false;
  };
};

using V_RKRK = std::vector<derivativeImprover *>;

using map_P_d = std::map<netP *, double>;
using map_P_Vd = std::map<netP *, V_d>;
using map_P_VVd = std::map<netP *, VV_d>;

void remesh(Network &obj0, const V_d &v)
{
  for (auto len : v)
  {
    bool found = false;
    Print(len, red);
    int count = 0;
    do
    {
      found = false;
      for (const auto &l : obj0.getLines())
      {
        auto f = l->getFaces();
        if (f.size() == 2)
          if (Norm(Cross(f[0]->getNormal(), f[1]->getNormal())) < 1E-5)
            found = l->flipIfIllegal();
        if (!found && l->length() > len)
          l->divide();
      }
    } while (found || count++ < 10);
  };
};

void smoothing(const V_netPp &ps)
{
  for (auto p : ps)
  {
    VV_d normals;
    for (auto f : p->getFaces())
      normals.emplace_back(f->getNormal());

    bool isEdgePoint = false;
    for (auto l : p->getLines())
      if (l->getFaces().size() != 2)
        isEdgePoint = true; //端の点は平滑かしない

    bool isflat = true;
    if (!isEdgePoint)
    {
      for (auto f : p->getFaces())
      {
        for (auto i = 0; i < normals.size() - 1; i++)
          for (auto j = i + 1; j < normals.size(); j++)
            if (!Norm(Cross(normals[i], normals[j])) < 1E-8)
              isflat = false;
      }
    }

    if (isflat)
    {
      VV_d X;
      for (auto q : p->getNeighbors())
        X.emplace_back(q->getX());

      p->setX(Mean(X));
    }
  };
}

//--------------------------------------------------------
//RBF補間による補間
// class depth_searcher : public searcher<networkPoint>
// {
// public:
//   int depth_lim;
//   depth_searcher(int depth_lim_IN) : depth_lim(depth_lim_IN), searcher<networkPoint>(){};
//   bool condEnterLine(const networkPoint *p, const netLp l) override
//   {
//     if (this->depth < depth_lim)
//       return true;
//     else
//       return false;
//   };
//   bool condGetObject(const netLp l, const netP *P) override
//   {
//     if (MemberQ(this->getNetworks(), P->getNetwork()))
//       return true;
//     else
//       return false;
//   };
// };
//================================================
int main()
{

  bool _runall_ = false;

  bool _stopall_ = true;

  if (false)
  {
    for (auto i = 2; i < 11; i++)
    {
      VV_d gw = GaussianQuadratureWeights(5, 0., 1.);
      VV_d gwgw;
      for (const auto &tw0 : gw)
      {
        double xi = tw0[0];
        double w_xi = tw0[1];
        for (const auto &tw1 : gw)
        {
          double h = tw1[0];
          double w_h = tw1[1];
          gwgw.push_back({xi, h * (1. - xi), w_xi * w_h * (1. - xi), w_xi * w_h /*1-xiなし*/});
        }
      }
      std::cout << "const static std::vector<std::vector<double>> __GWGW" << std::setprecision(15) << std::to_string(i) << "__ = " << gwgw << ";" << std::endl;
    }
  }

  Print("MORE EXAMPLES FOR STUDENTS", Green);

  Print("Point searcher advanced", green);
  //-----------------------------------
  if (false)
  {
    GNUPLOT plot;
    //    plot.Set({{"key",""},{"xrange","[-100:100]"},{"yrange","[-200:200]"}});
    NetworkObj obj("./obj/camel.obj");
    obj.scale(2.);
    V_d cardinal = {1000., 100., .2};
    NetworkObj obj0("./obj/oono.obj");
    obj0.translate({-80., -50., 5.});
    //    plot.SaveVectorData(getVectorData(obj.Faces),{{"arrowstyle","1"},{"title","obj0"}});
    plot.SaveVectorData(getVectorData(obj0.Faces), {{"arrowstyle", "2"}, {"title", "obj1"}});
    plot.plot3d();
    plot.Clear();
    std::cin.ignore();
    /////////////////////////////////
    auto angles = [](const VV_d &xyz) {
      return V_d{VectorAngle(xyz[0], xyz[1], xyz[2]),
                 VectorAngle(xyz[2], xyz[0], xyz[1]),
                 VectorAngle(xyz[1], xyz[2], xyz[0])};
    };

    auto v = Subdivide(5., 1., 100);
    for (auto len : v)
    {
      bool found = false;
      Print(len, red);
      int count = 0;
      do
      {
        found = false;
        for (const auto &l : obj0.getLines())
        {
          auto f = l->getFaces();
          if (f.size() == 2)
            if (Norm(Cross(f[0]->getNormal(), f[1]->getNormal())) < 1E-5)
            {
              found = l->flipIfIllegal();
              if (!found && l->length() > len)
                l->divide();
            }
          // for (auto &p : l->getPoints())
          //   for (auto &ll : p->getLines())
          //     if (ll->length() > len)
          //       ll->divideIfIllegal();
        }
      } while (count++ < 10);

      // Print(len,red);
      // for(const auto& l:obj0.getLines()){
      //   if(l->length() > len)
      //     l->divide();
      // }

      // bool found = false;
      // for(auto i=0; i<10; i++){
      //   for(auto l:obj0.getLines()){
      //     auto f = l->getFaces();
      //     if(f.size()==2){
      //       if(Norm(Cross(f[0]->getNormal(),f[1]->getNormal()))<1E-5){
      //         l->flipIfIllegal();
      //         found = true;
      //       }
      //     }
      //   }
      // }

      obj0.getLines()[0]->flip();
      //    plot.SaveVectorData(getVectorData(obj.Faces),{{"arrowstyle","1"},{"title","obj0"}});
      plot.SaveVectorData(getVectorData(obj0.Faces), {{"arrowstyle", "2"}, {"title", "obj1"}});
      plot.plot3d();
      plot.Clear();
    }

    //    plot.SaveVectorData(getVectorData(obj.Faces),{{"arrowstyle","1"},{"title","obj0"}});
    /////////////////////////////////
    obj0.getLines()[0]->flip();
    //    plot.SaveVectorData(getVectorData(obj.Faces),{{"arrowstyle","1"},{"title","obj0"}});
    plot.SaveVectorData(getVectorData(obj0.Faces), {{"arrowstyle", "2"}, {"title", "obj1"}});
    plot.plot3d();
    plot.Clear();
    std::cin.ignore();
    //
    if (_stopall_)
      std::cin.ignore();
  }

  Print("Point searcher advanced", green);
  //-----------------------------------
  if (true)
  {
    GNUPLOT plot;

    plot.Set({{"key", ""}, {"xrange", "[-100:100]"}, {"yrange", "[-200:200]"}});
    NetworkObj obj("./obj/camel.obj");
    obj.scale(2.);
    V_d cardinal = {1000., 100., .2};
    NetworkObj water("./obj/oyama2.obj");
    water.translate({-80., -50., 5.});
    remesh(water, Subdivide(50., 10., 20));

    // for (auto i = 0; i < 50; i++)
    //   smoothing(water.getPoints());


    plot.SaveVectorData(getVectorData(obj.Faces), {{"arrowstyle", "1"}, {"title", "obj0"}});
    plot.SaveVectorData(getVectorData(water.Faces), {{"arrowstyle", "2"}, {"title", "obj1"}});

    plot.plot3d();
    std::cin.ignore();
    plot.Clear();

    /*point_searcher_example_code*/
    for (auto i = 0; i < 100; i++)
    {
      water.translate({(double)2. * sin(2. * M_PI * i / 100.), (double)3. * cos(2. * M_PI * i / 100.), (double)4. * cos(2. * M_PI * i / 100.)});
      Network xnet(V_Netp{&obj, &water});
      BEM::selector s({&xnet});
      s(&water, cardinal);
      plot.SaveVectorData(getVectorData(s.faces), {{"arrowstyle", "2"}, {"title", "obj1"}});
      mk_vtu("./vtu/selected0_" + std::to_string(i) + ".vtu", s.faces);

      s(&obj, cardinal);
      plot.SaveVectorData(getVectorData(s.faces), {{"arrowstyle", "1"}, {"title", "obj2"}});
      mk_vtu("./vtu/selected1_" + std::to_string(i) + ".vtu", s.faces);

      //      mk_vtu("./vtu/selected2_"+std::to_string(i)+".vtu", xnet.Faces);
      //plot.plot3d();
      plot.Clear();
    }
    /*point_searcher_example_code*/
    if (_stopall_)
      std::cin.ignore();
  }

  Print("Point searcher advanced", green);
  //-----------------------------------
  if (false || _runall_)
  {
    GNUPLOT plot;

    //case 3
    //-----------------------------
    plot.Set({{"key", ""}, {"xrange", "[-100:100]"}, {"yrange", "[-200:200]"}});
    NetworkObj obj("./obj/camel.obj", "obj");
    obj.scale(2.);
    // obj.translate({-80.,0.,0.});
    // obj.translate({0.,10.,0.});
    V_d cardinal = {1000., 100., .2};
    NetworkObj water("./obj/oyama2.obj", "water");
    //water.rotate(M_PI/2.,{1,0,0});
    water.translate({-80., -50., 5.});
    //-----------------------------

    plot.SaveVectorData(getVectorData(obj.Faces), {{"arrowstyle", "1"}, {"title", "obj0"}});
    plot.SaveVectorData(getVectorData(water.Faces), {{"arrowstyle", "2"}, {"title", "obj1"}});

    plot.plot3d();
    std::cin.ignore();
    plot.Clear();

    /*point_searcher_example_code*/
    for (auto i = 0; i < 500; i++)
    {
      water.translate({(double)cos(2. * M_PI * i / 500.) / 2., (double)sin(2. * M_PI * i / 500.) / 2., (double)cos(2. * M_PI * i / 500.) / 10.});
      Network xnet(V_Netp{&obj, &water});
      BEM::selector s({&xnet}, 0);
      s(&water, cardinal);
      plot.SaveVectorData(getVectorData(s.faces), {{"arrowstyle", "2"}, {"title", "obj1"}});
      mk_vtu("./vtu/selected0_" + std::to_string(i) + ".vtu", s.faces);

      s(&obj, cardinal);
      plot.SaveVectorData(getVectorData(s.faces), {{"arrowstyle", "1"}, {"title", "obj2"}});
      mk_vtu("./vtu/selected1_" + std::to_string(i) + ".vtu", s.faces);

      V_netFp faces;
      for (const auto &f : s.faces)
        for (const auto &routeP : f->getPointsCutFaces())
          for (const auto &tri_ps : network::triangulate(routeP, f->getNormal(), 1E-4))
            xnet.Faces.emplace_back(new networkFace(&obj, link(tri_ps, &xnet)));
      mk_vtu("./vtu/selected2_" + std::to_string(i) + ".vtu", xnet.Faces);

      //      mk_vtu("./vtu/selected2_"+std::to_string(i)+".vtu", xnet.Faces);
      //      plot.plot3d();
      plot.Clear();
    }
    /*point_searcher_example_code*/
    if (_stopall_)
      std::cin.ignore();
  }

  Print("SHOW ALL", green);
  //-----------------------------------
  if (false || _runall_)
  { //_showall0_
    NetworkObj obj("./obj/oono.obj");
    obj.rotate(M_PI / 2., {1., 0., 0.});
    NetworkX water({5, 5}, {30.01, 30.01, 1 / 2.}, 2);
    V_netFp faces;
    V_netPp points;
    VV_d bounds = {{-0.55, 0.55}, {-0.55, 0.55}, {0.05, 1.5}};

    //GNUPLOT plot(std::map<std::string,std::string>{{"term","postscript eps enhanced color linewidth 2"}});
    GNUPLOT plot;
    plot.Set({{"key", ""}});
    plot.SaveVectorData(getVectorData(obj.Faces), {{"arrowstyle", "1"}, {"title", "obj"}});
    plot.SaveVectorData(getVectorData(water.Faces), {{"arrowstyle", "3"}, {"title", "water"}});
    plot.plot3d();
    //_showall0_
    if (_stopall_)
      std::cin.ignore();
  }

  Print("INTERSETION", green);
  //-----------------------------------
  if (false || _runall_)
  {
    /*
    //_INTERSETIONdetail_
    干渉したオブジェクトを表示
    //_INTERSETIONdetail_
    */
    /*_INTERSETION_*/
    NetworkObj obj("./obj/oono.obj");
    obj.rotate(M_PI / 2., {1., 0., 0.});
    NetworkX water({5, 5}, {30.01, 30.01, 1 / 2.}, 2);
    Network intersection({&obj, &water});
    V_netFp faces;
    V_netPp points;

    GNUPLOT plot;
    plot.Set({{"key", ""}});
    plot.SaveVectorData(getVectorData(obj.Faces), {{"arrowstyle", "1"}, {"title", "obj"}});
    plot.SaveVectorData(getVectorData(water.Faces), {{"arrowstyle", "3"}, {"title", "water"}});
    plot.SaveData(obj3D::extractX(intersection.Points), {{"title", "intersection"}});
    plot.plot3d();
    /*_INTERSETION_*/
    if (_stopall_)
      std::cin.ignore();
  }

  Print("BISECTIONING", green);
  //-----------------------------------
  if (false || _runall_)
  { //_BISECTIONING0_
    NetworkObj obj("./obj/oono.obj");
    obj.rotate(M_PI / 2., {1., 0., 0.});
    NetworkW water({5, 5}, {30.01, 30.01, 4.}, 10);
    Network corssNetwork({&obj, &water});
    V_netFp faces;
    V_netPp points;
    VV_d bounds = {{-0.55, 0.55}, {-0.55, 0.55}, {0.05, 1.5}};

    GNUPLOT plot;
    plot.Set({{"key", ""}});
    int i = 0;
    netF *f;
    while (i++ < 3)
    {
      mk_vtu("./vtu/oono" + std::to_string(i) + ".vtu", obj.Faces);
      for (const auto &l : obj.getLines())
        longerLine(l)->divide();
      obj.displayStates();
      plot.SaveVectorData(getVectorData(obj.Faces), {{"arrowstyle", "3"}, {"title", "bisectioning"}});
      plot.plot3d();

      std::cin.ignore();
      plot.Clear();
    }

    plot.SaveVectorData(getVectorData(obj.Faces), {{"arrowstyle", "1"}, {"title", "obj"}});
    plot.SaveVectorData(getVectorData(water.Faces), {{"arrowstyle", "3"}, {"title", "water"}});
    plot.plot3d();
    //_BISECTIONING0_
    if (_stopall_)
      std::cin.ignore();
  }

  Print("BISECTIONING", green);
  //-----------------------------------
  if (false || _runall_)
  { //_BISECTIONING0_
    if (false)
    {
      NetworkObj obj("./obj/oono.obj");
      obj.rotate(M_PI / 2., {1., 0., 0.});
      NetworkW water({5, 5}, {30.01, 30.01, 4.}, 10);

      for (auto i = 0; i < 30; i++)
      {
        int j = 0;
        while (j < 10)
        {
          for (const auto &l : obj.getLines())
          {
            auto ll = longerLine(l);
            if (ll->length() > 20. / (i + 1.))
              ll->divide();
          }
          j++;
        }
        mk_vtu("./vtu/oono_div" + std::to_string(i) + ".vtu", obj.Faces);
      }
    }

    {
      NetworkObj obj("./obj/oono.obj");
      obj.rotate(M_PI / 2., {1., 0., 0.});
      NetworkW water({5, 5}, {30.01, 30.01, 4.}, 10);

      for (auto i = 0; i < 39; i++)
      {
        int j = 0;
        while (j < 10)
        {
          for (const auto &l : obj.getLines())
          {
            auto ll = longerLine(l);
            if (ll->length() > 20. - i / 2)
              ll->divide();
          }
          j++;
        }
        mk_vtu("./vtu/oono_linear" + std::to_string(i) + ".vtu", obj.Faces);
      }
    }
    //_BISECTIONING0_
    if (_stopall_)
      std::cin.ignore();
  }

  //=============================================
  Print("SHOW ALL", Green);
  //-----------------------------------
  if (false || _runall_)
  { //_showall1_
    NetworkObj obj("./obj/tank.obj");
    NetworkW water({1, 1}, {1.01, 1.01, 1 / 2.}, .1);
    Network corssNetwork({&obj, &water});
    V_netFp faces;
    V_netPp points;
    VV_d bounds = {{-0.55, 0.55}, {-0.55, 0.55}, {0.05, 1.5}};

    GNUPLOT plot;
    plot.Set({{"key", ""}});
    plot.SaveVectorData(getVectorData(obj.Faces), {{"arrowstyle", "1"}, {"title", "obj"}});
    plot.SaveVectorData(getVectorData(water.Faces), {{"arrowstyle", "3"}, {"title", "water"}});
    plot.plot3d();
    std::cin.ignore();
    //_showall1_
    if (_stopall_)
      std::cin.ignore();
  }
  Print("searcher", Green);
  //==============================================

  Print("getPointsOnLines()", green);
  //-----------------------------------
  if (false || _runall_)
  {
    /*getPointsOnLines_example_code*/
    GNUPLOT plot;
    NetworkObj obj("./obj/tank2.obj");
    NetworkX water({5, 5}, {.8, .9, .6}, 1);

    for (auto j = 0; j < 3; j++)
    {
      for (const auto &l : obj.getLines())
        longerLine(l)->divide();
    }
    water.scale(1.1);

    Network XNetwork({&obj, &water});

    int l = 0;
    for (const auto &f : water.Faces)
      plot.SaveData(obj3D::extractX(f->getPointsOnLines()), {{"w", "lp"}, {"loop", ""}, {"lc", std::to_string(l++)}, {"notitle", ""}});

    plot.plot3d();
    /*getPointsOnLines_example_code*/
    if (_stopall_)
      std::cin.ignore();
  }

  Print("getPointsCutLines()", green);
  //-----------------------------------
  if (false || _runall_)
  {
    /*getPointsCutLines_example_code*/
    GNUPLOT plot;
    NetworkObj obj("./obj/tank2.obj");
    NetworkX water({5, 5}, {.8, .9, .6}, 1);

    for (auto j = 0; j < 3; j++)
    {
      for (const auto &l : obj.getLines())
        longerLine(l)->divide();
    }
    water.scale(1.1);

    Network XNetwork({&obj, &water});

    //    plot.SaveVectorData(getVectorData(water.Faces),{{"arrowstyle","3"},{"title","water.Faces"}});

    int ll = 0;
    for (const auto &f : water.Faces)
      for (const auto &ps : f->getPointsCutLines())
        mk_vtu("./vtu/getPointsCutLines" + std::to_string(ll++) + ".vtu", {ps});
    /*getPointsCutLines_example_code*/
    if (_stopall_)
      std::cin.ignore();
  }

  Print("getPointsOnLinesDivided()", green);
  //-----------------------------------
  if (false || _runall_)
  {
    /*getPointsOnLinesDivided_example_code*/
    GNUPLOT plot;
    plot.Set({{"key", ""}});
    NetworkObj obj("./obj/tank2.obj");
    NetworkX water({5, 5}, {.8, .9, .6}, 1);

    for (auto j = 0; j < 3; j++)
      for (const auto &l : obj.getLines())
        longerLine(l)->divide();

    water.scale(1.1);

    Network XNetwork({&obj, &water});

    //plot.SaveVectorData(getVectorData(water.Faces),{{"arrowstyle","1"},{"title",""}});
    int ll = 0;
    bool fin = false;
    for (const auto &f : water.Faces)
      for (const auto &ps : f->getPointsOnLinesDivided())
        mk_vtu("./vtu/getPointsOnLinesDivided" + std::to_string(ll++) + ".vtu", {ps});
    /*getPointsOnLinesDivided_example_code*/
    if (_stopall_)
      std::cin.ignore();
  }

  Print("getPointsCutFaces()", green);
  //-----------------------------------
  if (false || _runall_)
  {
    /*getPointsCutFaces_example_code*/
    GNUPLOT plot;
    NetworkObj obj("./obj/tank2.obj");
    NetworkX water({5, 5}, {.8, .9, .6}, 1);

    for (auto j = 0; j < 3; j++)
    {
      for (const auto &l : obj.getLines())
        longerLine(l)->divide();
    }
    water.scale(1.1);

    Network XNetwork({&obj, &water});

    int ll = 0;
    bool fin = false;
    for (const auto &f : water.Faces)
      for (const auto &ps : f->getPointsCutFaces())
        mk_vtu("./vtu/getPointsCutFaces" + std::to_string(ll++) + ".vtu", {ps});
    /*getPointsCutFaces_example_code*/
    if (_stopall_)
      std::cin.ignore();
  }

  Print("getPointsCutFaces()とtriangulate", green);
  //-----------------------------------
  if (false || _runall_)
  {
    /*getPointsCutFaces_triangulate_example_code*/
    GNUPLOT plot;
    NetworkObj obj("./obj/tank2.obj");
    NetworkX water({5, 5}, {.8, .9, .6}, 1);

    for (auto j = 0; j < 3; j++)
    {
      for (const auto &l : obj.getLines())
        longerLine(l)->divide();
    }
    water.scale(1.1);

    Network XNetwork({&obj, &water});
    int ll = 0;
    bool fin = false;
    for (const auto &f : water.Faces)
      for (const auto &ps : f->getPointsCutFaces())
        for (const auto &pp : network::triangulate(ps, f->getNormal()))
          mk_vtu("./vtu/getPointsCutFaces_triangulate" + std::to_string(ll++) + ".vtu", {pp});
    /*getPointsCutFaces_triangulate_example_code*/
    if (_stopall_)
      std::cin.ignore();
  }

  Print("getPointsCutFaces()とtriangulate", green);
  //-----------------------------------
  if (false || _runall_)
  {
    /*getPointsCutFaces_triangulate_example_code*/
    GNUPLOT plot;
    NetworkObj obj("./obj/tank2.obj");
    NetworkX water({5, 5}, {.8, .9, .6}, 1);
    VV_netPp triangle, route;

    for (auto j = 0; j < 3; j++)
    {
      for (const auto &l : obj.getLines())
        longerLine(l)->divide();
    }
    water.scale(1.1);

    Network XNetwork({&obj, &water});
    int l(0);
    for (const auto &f : water.Faces)
    {
      auto ps = f->getPointsCutFaces();
      if (ps.size() > 0)
      {
        mk_vtu("./vtu/faces" + std::to_string(l++) + ".vtu", {ps});
        for (const auto &routeP : ps)
        {
          if (!routeP.empty())
          {
            Print(routeP);
            mk_vtu("./vtu/route" + std::to_string(l++) + ".vtu", {routeP});
          }
          // for(const auto& tri_ps:network::triangulate(routeP,f->getNormal())){
          //   triangle.emplace_back(tri_ps);
          //   //XNetwork.Faces.emplace_back(new networkFace(&water, link(tri_ps,&XNetwork)));
          // }
        }
      }
    }
    // for(const auto& f:obj.Faces)
    //   if(f->intersectQ())
    //     for(const auto& routeP:f->getPointsCutFaces())
    //       for(const auto& tri_ps:network::triangulate(routeP,f->getNormal())){
    //         triangle.emplace_back(tri_ps);
    //         XNetwork.Faces.emplace_back(new networkFace(&obj, link(tri_ps,&XNetwork)));
    //       }

    //    mk_vtu("./vtu/triangulate2.vtu", route);
    //mk_vtu("./vtu/triangulate.vtu", XNetwork.Faces);

    //plot.SaveVectorData(getVectorData(XNetwork.Faces),{{"arrowstyle","1"},{"title",""}});
    //plot.plot3d();
    //std::cin.ignore();
    /*getPointsCutFaces_triangulate_example_code*/
    if (_stopall_)
      std::cin.ignore();
  }

  //-----------------------------------
  if (false || _runall_)
  {
    /*getPointsCutFaces_triangulate_example2_code*/
    GNUPLOT plot;
    NetworkObj obj("./obj/tank2.obj");
    NetworkX water({5, 5}, {.8, .9, .6}, 1);
    VV_netPp triangle, route;

    for (auto j = 0; j < 3; j++)
    {
      for (const auto &l : obj.getLines())
        longerLine(l)->divide();
    }
    water.scale(1.1);

    Network XNetwork({&obj, &water});
    int l(0);
    for (const auto &f : water.Faces)
    {
      auto ps = f->getPointsCutFaces();
      if (ps.size() > 0)
      {
        mk_vtu("./vtu/faces" + std::to_string(l++) + ".vtu", {ps});
        for (const auto &routeP : ps)
        {
          if (!routeP.empty())
          {
            Print(routeP);
            mk_vtu("./vtu/route" + std::to_string(l++) + ".vtu", {routeP});
          }
        }
      }
    }
    /*getPointsCutFaces_triangulate_example2_code*/
    if (_stopall_)
      std::cin.ignore();
  }

  //===================================

  Print("Face searcher", green);
  //-----------------------------------
  if (false || _runall_)
  {
    /*face_searcher_example_code*/
    auto plot = [](searcher<netF> &S, const bool getfirstobj = true) {
      NetworkObj obj("./obj/tank.obj");
      NetworkW water({20, 3}, {2., 2., 1 / 2.}, .1);
      Network corssNetwork({&obj, &water});
      VV_d bounds = {{-0.55, 0.55}, {-0.55, 0.55}, {0.05, 1.5}};
      GNUPLOT plot;
      plot.Set({{"key", ""}, {"title", "'Face searcher'"}});

      S.set(network::takeNearest(obj3D::takeInsideOfBounds(obj.Faces, bounds), {0, 0, 0.2}));
      S.search(getfirstobj);
      plot.SaveVectorData(getVectorData(S.getObjects()), {{"arrowstyle", "1"}, {"title", "search 1"}});

      S.clear();
      S.set(network::takeNearest(obj3D::takeInsideOfBounds(water.Faces, bounds), {0, 0, 0.2}));
      S.search(getfirstobj);
      plot.SaveVectorData(getVectorData(S.getObjects()), {{"arrowstyle", "2"}, {"title", "search 2"}});

      plot.plot3d();
      std::cin.ignore();
    };
    {
      Print("searcher<netF>", Blue);
      searcher<netF> S;
      plot(S);
    }
    {
      Print("searcher1<netF>", Blue);
      searcher1<netF> S;
      plot(S);
    }
    {
      Print("searcher2<netF>", Blue);
      searcher2<netF> S;
      plot(S);
    }
    {
      Print("searcher3<netF>", Blue);
      searcher3<netF> S;
      plot(S);
    }
    {
      Print("searcher4<netF>", Blue);
      searcher4<netF> S;
      plot(S);
    }
    {
      Print("searcher5<netF>", Blue);
      searcher5<netF> S;
      plot(S);
    }
    {
      Print("searcher6<netF>", Blue);
      searcher6<netF> S;
      plot(S, false);
    }
    {
      Print("searcher7<netF>", Blue);
      searcher7<netF> S;
      plot(S, false);
    }
    {
      Print("searcher8<netP>", Blue);
      searcher8<netF> S;
      plot(S);
    }
    /*face_searcher_example_code*/
  }

  Print("Point searcher", green);
  //-----------------------------------
  if (false || _runall_)
  {
    /*point_searcher_example_code*/
    auto plot = [](searcher<netP> &S, const bool getfirstobj = true) {
      NetworkObj obj("./obj/tank.obj");
      NetworkX water({10, 10}, {2., 2., 1 / 2.}, 2);
      Network corssNetwork({&obj, &water});
      VV_d bounds = {{-0.55, 0.55}, {-0.55, 0.55}, {0.05, 1.5}};
      GNUPLOT plot;
      plot.Set({{"key", ""}, {"title", "'Point searcher'"}});

      S.set(network::takeNearest(obj3D::takeInsideOfBounds(obj.Points, bounds), {0, 0, 0.2}));
      S.search(getfirstobj);
      plot.SaveData(obj3D::extractX(S.getObjects()), {{"ps", "2"}, {"pt", "7"}, {"lc", "'pink'"}, {"title", "getObjects()"}});
      plot.SaveVectorData(getVectorData(S.getReachedLines()), {{"arrowstyle", "1"}, {"title", "getReachedLines()"}});
      plot.SaveVectorData(getVectorData(S.getEnteredLines()), {{"arrowstyle", "2"}, {"title", "getEnteredLines()"}});

      S.clear();
      S.set(network::takeNearest(obj3D::takeInsideOfBounds(water.Points, bounds), {0, 0, 0.2}));
      S.addNetwork(&corssNetwork);

      plot.SaveData(obj3D::extractX(S.getObjects()), {{"ps", "2"}, {"pt", "7"}, {"lc", "'blue'"}, {"title", "getObjects()"}});

      plot.SaveVectorData(getVectorData(S.getReachedLines()), {{"arrowstyle", "3"}, {"title", "getReachedLines()"}});
      plot.SaveVectorData(getVectorData(S.getEnteredLines()), {{"arrowstyle", "4"}, {"title", "getEnteredLines()"}});

      plot.plot3d();
      std::cin.ignore();
    };
    {
      Print("searcher<netP>", Blue);
      searcher<netP> S;
      plot(S);
    }
    {
      Print("searcher1<netP>", Blue);
      searcher1<netP> S;
      plot(S);
    }
    {
      Print("searcher2<netP>", Blue);
      searcher2<netP> S;
      plot(S);
    }
    {
      Print("searcher3<netP>", Blue);
      searcher3<netP> S;
      plot(S);
    }
    {
      Print("searcher4<netP>", Blue);
      searcher4<netP> S;
      plot(S);
    }
    {
      Print("searcher5<netP>", Blue);
      searcher5<netP> S;
      plot(S);
    }
    {
      Print("searcher6<netP>", Blue);
      searcher6<netP> S;
      plot(S, false);
    }
    {
      Print("searcher7<netP>", Blue);
      searcher7<netP> S;
      plot(S, false);
    }
    {
      Print("searcher8<netP>", Blue);
      searcher8<netP> S;
      plot(S);
    }
    /*point_searcher_example_code*/
    if (_stopall_)
      std::cin.ignore();
  }

  Print("Point searcher advanced", green);
  //-----------------------------------
  if (false || _runall_)
  {
    GNUPLOT plot;
    NetworkObj obj("./obj/tank.obj");
    NetworkX water({10, 10}, {4., 4., 1 / 2.}, 0);
    Network XNetwork({&obj, &water});

    for (const auto &f : water.Faces)
      if (f->intersectQ())
        for (const auto &routeP : f->getPointsCutFaces())
          for (const auto &tri_ps : network::triangulate(routeP, f->getNormal()))
            new networkFace(&XNetwork, link(tri_ps, &XNetwork));

    for (const auto &f : obj.Faces)
      if (f->intersectQ())
        for (const auto &routeP : f->getPointsCutFaces())
          for (const auto &tri_ps : network::triangulate(routeP, f->getNormal()))
            new networkFace(&XNetwork, link(tri_ps, &XNetwork));

    /*point_searcher_example_code*/
    auto search = [&obj, &water, &XNetwork, &plot](Network &target,
                                                   const bool getfirstobj = true) {
      searcher9<netP> S;
      VV_d bounds = {{-0.55, 0.55}, {-0.55, 0.55}, {0.05, 1.5}};

      plot.Set({{"key", ""}, {"title", "'Point searcher'"}});

      //S
      S.set(network::takeNearest(obj3D::takeInsideOfBounds(target.Points, bounds), {0, 0, 0.2}));
      S.addNetwork(&XNetwork);

      S.search(getfirstobj);

      plot.SaveData(obj3D::extractX(S.getObjects()), {{"ps", "2"}, {"pt", "7"}, {"lc", "'blue'"}, {"title", "getObjects()"}});
      plot.SaveData(obj3D::extractX(S.getObjects_()), {{"ps", "2"}, {"pt", "7"}, {"lc", "'pink'"}, {"title", "getObjects()"}});

      plot.SaveVectorData(getVectorData(S.getReachedLines()), {{"arrowstyle", "3"}, {"title", "getReachedLines()"}});
      plot.SaveVectorData(getVectorData(S.getEnteredLines()), {{"arrowstyle", "4"}, {"title", "getEnteredLines()"}});
    };

    Print("searcher9<netP>", Blue);
    search(water, false);
    search(obj, false);

    plot.plot3d();
    std::cin.ignore();

    /*point_searcher_example_code*/
    if (_stopall_)
      std::cin.ignore();
  }

  Print("Point searcher advanced", green);
  //-----------------------------------
  if (false || _runall_)
  {
    GNUPLOT plot;
    NetworkObj obj("./obj/tank2.obj");
    NetworkX water({5, 5}, {.7, .7, .6}, 1);

    for (auto j = 0; j < 3; j++)
    {
      for (const auto &l : obj.getLines())
        longerLine(l)->divide();
    }

    //干渉部分が多すぎるという問題
    water.scale(1.1);

    Network XNetwork({&obj, &water});

    int l = 0;
    for (const auto &f : water.Faces)
      if (f->intersectQ())
        for (const auto &routeP : f->getPointsCutFaces())
        {
          plot.SaveData(obj3D::extractX(routeP), {{"w", "lp"}, {"loop", ""}, {"lc", std::to_string(l++)}, {"notitle", ""}});
          for (const auto &tri_ps : network::triangulate(routeP, f->getNormal()))
            XNetwork.Faces.emplace_back(new networkFace(&water, link(tri_ps, &XNetwork)));
        }

    for (const auto &f : obj.Faces)
      if (f->intersectQ())
        for (const auto &routeP : f->getPointsCutFaces())
          for (const auto &tri_ps : network::triangulate(routeP, f->getNormal()))
            XNetwork.Faces.emplace_back(new networkFace(&obj, link(tri_ps, &XNetwork)));

    /*point_searcher_example_code*/
    auto search = [&obj, &water, &XNetwork, &plot](Network &target, const bool getfirstobj = true) {
      searcher9<netP> S;
      VV_d bounds = {{-0.55, 0.55}, {-0.55, 0.55}, {0.05, 1.5}};

      plot.Set({{"key", ""}, {"title", "'Point searcher'"}});

      //S
      S.set(network::takeNearest(obj3D::takeInsideOfBounds(target.Points, bounds), {0, 0, 0.2}));
      S.addNetwork(&XNetwork);

      S.search(getfirstobj);

      network::setStatus(obj.Points, false);
      network::setStatus(water.Points, false);
      network::setStatus(XNetwork.Points, false);

      for (const auto &p : S.getObjects())
      {
        p->setStatus(true);
      }
      for (const auto &p : S.getObjects_())
      {
        p->setStatus(true);
      }

      V_netFp faces;

      for (const auto &f : target.Faces)
      {
        if (f->getNetwork() == &target && network::AllStatusTrue(f->getPoints()))
        {
          faces.emplace_back(f);
        }
      }

      for (const auto &f : XNetwork.Faces)
      {
        if (f->getNetwork() == &target && network::AllStatusTrue(f->getPoints()))
        {
          faces.emplace_back(f);
        }
      }

      plot.SaveData(obj3D::extractX(S.getObjects()), {{"ps", "2"}, {"pt", "7"}, {"lc", "'blue'"}, {"title", "getObjects()"}});
      plot.SaveData(obj3D::extractX(S.getObjects_()), {{"ps", "2"}, {"pt", "7"}, {"lc", "'pink'"}, {"title", "getObjects()"}});
      plot.SaveData(obj3D::extractX(S.getObjects__()), {{"ps", "2"}, {"pt", "8"}, {"lc", "'red'"}, {"title", "getObjects()"}});

      plot.SaveVectorData(getVectorData(S.getReachedLines()), {{"arrowstyle", "3"}, {"title", "getReachedLines()"}});
      plot.SaveVectorData(getVectorData(S.getEnteredLines()), {{"arrowstyle", "4"}, {"title", "getEnteredLines()"}});

      plot.SaveVectorData(getVectorData(faces), {{"arrowstyle", "5"}, {"title", "faces"}});
      plot.SaveVectorData(getVectorData(XNetwork.Faces), {{"arrowstyle", "6"}, {"title", "XNetwork.Faces"}});
    };

    Print("searcher9<netP>", Blue);
    search(water, false);
    search(obj, false);

    plot.plot3d();
    /*point_searcher_example_code*/
    if (_stopall_)
      std::cin.ignore();
  }

  Print("Point searcher", green);
  //-----------------------------------
  if (false || _runall_)
  {

    //干渉したら探査を終了するsearcher
    searcher3<netP> S;

    NetworkObj obj("./obj/tank.obj");
    NetworkW water({3, 3}, {2., 2., 1 / 2.}, .1);
    Network corssNetwork({&obj, &water});
    GNUPLOT plot;
    plot.Set({{"key", ""}, {"title", "'Point searcher'"}});

    S.set(network::takeNearest(water.Points, {0., 0., 0.2}));
    S.search();

    auto points = S.getObjects();
    auto faces = extractFaces(points);
    V_netFp faces_itl;
    for (const auto &f : faces)
    {
      if (f->intersectQ())
        faces_itl.emplace_back(f);
    }

    plot.SaveData(obj3D::extractX(S.getObjects()), {{"ps", "2"}, {"pt", "7"}, {"lc", "'pink'"}, {"title", "getObjects()"}});
    plot.SaveVectorData(getVectorData(faces_itl), {{"arrowstyle", "2"}, {"title", "extractFaces(points)"}});
    plot.SaveVectorData(getVectorData(S.getReachedLines()), {{"arrowstyle", "3"}, {"title", "S.getReachedLines()"}});

    plot.SaveVectorData(getVectorData(obj.getFaces()), {{"arrowstyle", "1"}, {"title", "obj.getFaces()"}});

    // 干渉したら探査を終了するsearcher
    //    searcher3<netP> S;

    plot.plot3d();
    if (_stopall_)
      std::cin.ignore();
  }

  Print("one way Point searcher in an interaction network", green);
  //-----------------------------------
  if (false || _runall_)
  {
    auto plot = [](searcher<netP> &S) {
      NetworkObj obj("./obj/tank.obj");
      NetworkW water({3, 3}, {2., 2., 1 / 2.}, .1);
      Network corssNetwork({&obj, &water});
      GNUPLOT plot;
      plot.Set({{"key", ""}, {"title", "'Point searcher interaction network'"}});

      S.set(corssNetwork.Points[0]);
      S.search();
      plot.SaveData(obj3D::extractX(S.getObjects()), {{"ps", "2"}, {"pt", "7"}, {"lc", "'pink'"}, {"title", "getObjects()"}});
      plot.SaveVectorData(getVectorData(S.getReachedLines()), {{"arrowstyle", "2"}, {"title", "getReachedLines()"}});

      plot.plot3d();
      std::cin.ignore();
    };
    {
      searcher4<netP> S;
      plot(S);
    }
    if (_stopall_)
      std::cin.ignore();
  }

  Print("Point searcher detecting interaction", green);
  //-----------------------------------
  if (false || _runall_)
  {
    //searcher_interaction
    class my_searcher : public searcher<netP>
    {
    public:
      my_searcher() : searcher<netP>(){};
      bool condEnterLine(const netP *p, const netLp l) override
      {
        if (!l->penetrateQ())
        {
          return true;
        }
        else
          return false;
      };
      bool condKeepSearch(const netLp l, const netP *P) override
      {
        if (!P->isXPoint())
        {
          return true;
        }
        else
          return false;
      };
    };
    //searcher_interaction
    {
      NetworkObj obj("./obj/tank.obj");
      NetworkW water({3, 3}, {2., 2., 1 / 2.}, .1);
      Network corssNetwork({&obj, &water});

      for (const auto &f : water.Faces)
        if (f->intersectQ())
          for (const auto &routeP : f->getPointsCutFacesBehind())
            for (const auto &tri_ps : network::triangulate(routeP, f->getNormal()))
              networkFace tmp_f(&water, link(tri_ps, &water));

      for (const auto &f : obj.Faces)
        if (f->intersectQ())
          for (const auto &routeP : f->getPointsCutFacesBehind())
            for (const auto &tri_ps : network::triangulate(routeP, f->getNormal()))
              networkFace tmp_f(&obj, link(tri_ps, &obj));

      GNUPLOT plot;
      plot.Set({{"key", ""}, {"title", "'Point searcher'"}});
      my_searcher S;
      S.set(network::takeNearest(water.Points, {0, 0, 0.2}));
      S.search();
      plot.SaveData(obj3D::extractX(S.getObjects()), {{"ps", "2"}, {"pt", "7"}, {"lc", "'pink'"}, {"title", "getObjects()"}});
      plot.SaveVectorData(getVectorData(S.getReachedLines()), {{"arrowstyle", "1"}, {"title", "getReachedLines()"}});

      my_searcher S_;
      S_.set(network::takeNearest(obj.Points, {0, 0, 0.2}));
      S_.search();
      plot.SaveData(obj3D::extractX(S_.getObjects()), {{"ps", "2"}, {"pt", "7"}, {"lc", "'web-blue'"}, {"title", "getObjects()"}});

      plot.SaveVectorData(getVectorData(S_.getReachedLines()), {{"arrowstyle", "2"}, {"title", "getReachedLines()"}});

      plot.plot3d();
    }
    if (_stopall_)
      std::cin.ignore();
  }

  Print("Point searcher detecting interaction", green);
  //-----------------------------------
  if (false || _runall_)
  {
    //searcher_interaction
    {
      NetworkObj obj("./obj/tank.obj");
      NetworkW water({3, 3}, {1.05, 1.02, 1 / 2.}, .1);
      Network XNetwork({&obj, &water});
      VV_d bounds = {{-0.55, 0.55}, {-0.55, 0.55}, {0.05, 1.5}};
      V_netFp faces_water, faces_obj;

      for (const auto &f : obj3D::takeInsideOfBounds(water.Faces, bounds))
        if (f->intersectQ())
          for (const auto &routeP : f->getPointsCutFacesBehind())
            for (const auto &tri_ps : network::triangulate(routeP, f->getNormal()))
              faces_water.emplace_back(new networkFace(&XNetwork, link(tri_ps, &XNetwork)));

      for (const auto &f : obj3D::takeInsideOfBounds(obj.Faces, bounds))
        if (f->intersectQ())
          for (const auto &routeP : f->getPointsCutFacesBehind())
            for (const auto &tri_ps : network::triangulate(routeP, f->getNormal()))
              faces_obj.emplace_back(new networkFace(&XNetwork, link(tri_ps, &XNetwork)));

      GNUPLOT plot;
      plot.Set({{"key", ""}, {"title", "'Point searcher'"}});

      plot.SaveVectorData(getVectorData(faces_water), {{"arrowstyle", "1"}, {"title", "faces"}}); //interaction faces
      plot.SaveVectorData(getVectorData(faces_obj), {{"arrowstyle", "2"}, {"title", "faces"}});   //interaction faces

      searcher8<netF> S(network::takeNearest(water.Faces, {0, 0, 0.2}));
      S.search();

      plot.SaveData(S.getObjects()[0]->getMeanLocation(), {{"ps", "3"}, {"pt", "7"}, {"lc", "'red'"}, {"title", "getObjects()"}});
      plot.SaveVectorData(getVectorData(S.getObjects()), {{"arrowstyle", "3"}, {"title", "S.getObjects()"}});   //interaction faces
      plot.SaveVectorData(getVectorData(S.getObjects_()), {{"arrowstyle", "4"}, {"title", "S.getObjects_()"}}); //interaction faces

      plot.plot3d();
    }
    if (_stopall_)
      std::cin.ignore();
  }

  Print("ACCESS TO FACES FROM POINTS", green);
  //===================
  if (true || _runall_)
  {
    NetworkObj obj("./obj/tank.obj");
    NetworkW water({1, 1}, {1.133, 1.13, 1 / 2.}, .1);
    Network corssNetwork({&obj, &water});
    VV_d bounds = {{-0.55, 0.55}, {-0.55, 0.55}, {0.05, 1.5}};

    GNUPLOT plot;
    plot.Set({{"key", ""}});

    auto p = network::takeNearest(obj3D::takeInsideOfBounds(obj.Points, bounds), {0, 0, 0});
    plot.SaveData({p->getX()}, {{"ps", "3"}, {"pt", "7"}, {"lc", "'pink'"}, {"title", "obj1 p"}});
    plot.SaveVectorData(getVectorData(p->getFaces()), {{"arrowstyle", "1"}, {"title", "obj1 p->getFaces()"}});

    p = network::takeNearest(obj3D::takeInsideOfBounds(water.Points, bounds), {0, 0, 0});
    plot.SaveData({p->getX()}, {{"ps", "3"}, {"pt", "7"}, {"lc", "'red'"}, {"title", "obj1 p"}});
    plot.SaveVectorData(getVectorData(p->getFaces()), {{"arrowstyle", "2"}, {"title", "obj2 p->getFaces()"}});

    plot.plot3d();
    if (_stopall_)
      std::cin.ignore();
  }

  Print("ACCESS TO POINTS FROM FACES", green);
  //===================
  if (false || _runall_)
  {
    NetworkObj obj("./obj/tank.obj");
    NetworkW water({1, 1}, {1.133, 1.13, 1 / 2.}, .1);
    obj.setBounds();
    water.setBounds();
    Network corssNetwork({&obj, &water});
    VV_d bounds = {{-0.55, 0.55}, {-0.55, 0.55}, {0.05, 1.5}};

    GNUPLOT plot;
    plot.Set({{"key", ""}});

    auto f = network::takeNearest(obj3D::takeInsideOfBounds(obj.Faces, bounds), {0, 0, 0});
    plot.SaveData({f->getMeanLocation()}, {{"ps", "3"}, {"pt", "7"}, {"lc", "'pink'"}, {"title", "obj1 p"}});
    plot.SaveData(getData(f->getPoints()), {{"loop", ""}, {"w", "lp"}, {"title", "takeInsideOfBounds"}});

    f = network::takeNearest(obj3D::takeInsideOfBounds(water.Faces, bounds), {0, 0, 0});
    plot.SaveData({f->getMeanLocation()}, {{"ps", "3"}, {"pt", "7"}, {"lc", "'red'"}, {"title", "obj1 p"}});
    plot.SaveData(getData(f->getPoints()), {{"loop", ""}, {"w", "lp"}, {"title", "takeInsideOfBounds"}});

    plot.plot3d();
    if (_stopall_)
      std::cin.ignore();
  }

  Print("BISECTIONING", green);
  //===================
  if (false || _runall_)
  {
    NetworkObj obj("./obj/tank.obj");
    GNUPLOT plot;
    plot.Set({{"key", ""}});

    int i = 0;
    while (i++ < 5)
    {
      for (const auto &l : obj3D::takeInsideOfBounds(obj.getLines(), {{-0.55, 0.55}, {-0.55, 0.55}, {-100, 100}}))
        longerLine(l)->divide();

      obj.displayStates();
      plot.SaveVectorData(getVectorData(obj.Faces), {{"arrowstyle", "3"}, {"title", "bisectioning"}});
      plot.plot3d();
      std::cin.ignore();
      plot.Clear();
    }
    if (_stopall_)
      std::cin.ignore();
  }

  /*bisection_example0_detail
    bisection_example0_detail*/
  /*bisection_example0_code*/
  Print("BISECTIONING", green);
  //===================
  if (false || _runall_)
  {
    std::string name = "bunny";
    NetworkObj obj("./obj/" + name + ".obj");
    std::vector<double> range;
    //----------------------
    if (name == "bunny")
    {
      range = {0.3, 0.3, 0.1};
    }
    else if (name == "camel")
    {
      range = {70 + 1. / 3., 90 + 1. / 3., 21. + 1. / 3.};
    }
    else if (name == "pumpkin")
    {
      obj.rotate(M_PI / 2., {1., 0., 0.});
      range = {100, 100, -110};
    }
    else if (name == "duck")
    {
      obj.rotate(M_PI / 2., {1., 0., 0.});
      range = {6, 6, 1};
    }
    NetworkX water({1, 1}, range, 1);
    //---------------------
    GNUPLOT plot_;
    plot_.SaveVectorData(getVectorData(obj.Faces), {{"arrowstyle", "4"}});
    plot_.plot3d();
    //
    GNUPLOT plot;
    for (auto i = 0; i < 20; i++)
    { //ここが増えるとより多くcrossをチェックするため，crossの周辺だけをより細かくする
      plot.Set({{"key", ""}});
      plot.UnSet({{"grid", ""}});
      plot.Set({{"style", "arrow 1 nohead lc 1 lw .5"}});
      //------------------
      std::vector<std::vector<std::vector<double>>> vecvecvec;
      for (const auto &face : water.Faces)
      {
        auto xyz = face->getLocations();
        vecvecvec.push_back({{xyz[0], xyz[1] - xyz[0]}});
        vecvecvec.push_back({{xyz[1], xyz[2] - xyz[1]}});
        vecvecvec.push_back({{xyz[2], xyz[0] - xyz[2]}});
      }
      plot.SaveVectorData(vecvecvec, {{"arrowstyle", "1"}, {"title", "divide #" + std::to_string(i)}});
      plot.SaveVectorData(getVectorData(obj.Faces), {{"arrowstyle", "2"}, {"title", "divide #" + std::to_string(i)}});
      //------------------

      mk_vtu("./vtu/bisection" + std::to_string(i) + ".vtu", water.Faces);
      mk_vtu("./vtu/" + name + std::to_string(i) + ".vtu", obj.Faces);
      //------------------
      //2分割
      Network XNetwork({&obj, &water});
      for (const auto &l : water.getLines())
        if (l->penetrateQ())
        {
          auto ll = longerLine(l);
          if (ll->length() > 10. - i)
            ll->divide();
        }
      water.setXPoints(obj, &water);
      //------------------
      //点の取得サーチ
      class my_searcher : public searcher<netP>
      {
      public:
        my_searcher(netP *p) : searcher<netP>(p){};
        bool condKeepSearch(const netLp l, const netP *P) override
        {
          if (l->penetrateQ())
            return false;
          else
            return true;
        };
      };
      my_searcher S(network::takeNearest(water.Points, {-100, -100, 0.}));
      S.search();
      auto ps = S.getObjects();
      //------------------
      VV_netPp points;
      for (auto &l : S.getReachedLines())
      {
        points.emplace_back(l->getPoints());
      }
      mk_vtu("./vtu/searcher" + std::to_string(i) + ".vtu", points);

      plot.SaveVectorData(getVectorData(S.getReachedLines()), {{"arrowstyle", "4"}, {"title", "divide #" + std::to_string(i)}});
      plot.SaveData(getData(S.getObjects()), {{"title", "S.getObjects()"}});
      plot.SaveData(getData(S.getObjects_()), {{"title", "S.getObjects_()"}});
      plot.SaveData(getData(ps), {{"title", "ps"}});
      plot.plot3d();
      std::cin.ignore();
      plot.Clear();

      if (false)
      {
        for (auto p : water.Points)
          if (!MemberQ(ps, p))
            p->Delete();
      }
    }
    if (_stopall_)
      std::cin.ignore();
  }
  /*bisection_example0_code*/

  Print("TRIANGULATE INTERSECTING FACES ROUTES POLYGONS", green);
  //===================
  if (false || _runall_)
  {
    NetworkObj obj("./obj/tank.obj");
    NetworkW water({3, 3}, {1.133, 1.13, 1 / 2.}, .1);
    Network corssNetwork({&obj, &water});
    GNUPLOT plot;

    plot.Set({{"key", ""}});

    for (const auto &f : network::takeIfIntersect(obj.Faces))
      for (const auto &routeP : f->getPointsCutFacesBehind())
        for (const auto &ps : network::triangulate(routeP, f->getNormal()))
          plot.SaveData(obj3D::extractX(ps), {{"loop", ""}, {"w", "lp"}, {"ps", "0"}, {"lc", "'pink'"}, {"notitle", ""}});

    for (const auto &f : network::takeIfIntersect(water.Faces))
      for (const auto &routeP : f->getPointsCutFacesBehind())
        for (const auto &ps : network::triangulate(routeP, f->getNormal()))
          plot.SaveData(obj3D::extractX(ps), {{"loop", ""}, {"w", "lp"}, {"ps", "0"}, {"lc", "'blue'"}, {"notitle", ""}});

    //plot.SaveData(getData(corssNetwork.Points[0]->getNeighbors_OneWay_r()),{{"w","lp"},{"lc","'red'"},{"lw","2"},{"title","corssNetwork"}});

    plot.plot3d();
    plot.Clear();
    if (_stopall_)
      std::cin.ignore();
  }

  Print("TRIANGULATE INTERSECTING FACES BEHIND AND FRONT", green);
  //===================
  if (false || _runall_)
  {
    NetworkObj obj("./obj/tank.obj");
    NetworkW water({2, 2}, {1.15, 1.15, 1 / 2.}, .1);
    Network corssNetwork({&obj, &water});

    GNUPLOT plot;
    plot.Set({{"key", ""}});

    VV_netPp vvp;
    for (const auto &f : network::takeIfIntersect(obj.Faces))
      for (const auto &routeP : f->getPointsCutFacesBehind())
        for (const auto &ps : network::triangulate(routeP, f->getNormal()))
          vvp.push_back(ps);
    plot.SaveVectorData(getVectorData(vvp), {{"arrowstyle", "1"}, {"title", "getPointsCutFacesBehind()"}});

    vvp.clear();
    for (const auto &f : network::takeIfIntersect(obj.Faces))
      for (const auto &routeP : f->getPointsCutFacesFront())
        for (const auto &ps : network::triangulate(routeP, f->getNormal()))
          vvp.push_back(ps);
    plot.SaveVectorData(getVectorData(vvp), {{"arrowstyle", "2"}, {"title", "getPointsCutFacesFront()"}});

    vvp.clear();
    for (const auto &f : network::takeIfIntersect(water.Faces))
      for (const auto &routeP : f->getPointsCutFacesBehind())
        for (const auto &ps : network::triangulate(routeP, f->getNormal()))
          vvp.push_back(ps);
    plot.SaveVectorData(getVectorData(vvp), {{"arrowstyle", "3"}, {"title", "getPointsCutFacesBehind()"}});

    vvp.clear();
    for (const auto &f : network::takeIfIntersect(water.Faces))
      for (const auto &routeP : f->getPointsCutFacesFront())
        for (const auto &ps : network::triangulate(routeP, f->getNormal()))
          vvp.push_back(ps);
    plot.SaveVectorData(getVectorData(vvp), {{"arrowstyle", "4"}, {"title", "getPointsCutFacesFront()"}});

    plot.plot3d();
    if (_stopall_)
      std::cin.ignore();
  }

  Print("SEARCHER", green);
  //===================
  if (false || _runall_)
  { /*SEARCHER*/
    NetworkObj obj("./obj/tank.obj");
    int i = 0;
    while (i++ < 2)
      for (const auto &l : obj3D::takeInsideOfBounds(obj.getLines(), {{-0.55, 0.55}, {-0.55, 0.55}, {-100, 100}}))
        longerLine(l)->divide();

    NetworkW water({6, 6}, {1.15, 1.15, 1 / 2.}, .1);
    Network corssNetwork({&obj, &water});
    GNUPLOT plot;
    plot.Set({{"key", ""}});

    class searcher_AvoidX : public searcher<networkPoint>
    {
    public:
      searcher_AvoidX(networkPoint *start_, bool TorF_ = true) : searcher<networkPoint>(start_){};
      bool condEnterLine(const networkPoint *p, const netLp l) override
      {
        if (!l->penetrateQ())
        {
          return true;
        }
        else
        {
          return false;
        }
      };
      bool condGetObject(const netLp l, const networkPoint *P) override
      {
        if (P->getX()[2] < .55)
        {
          return true;
        }
        else
        {
          return false;
        }
      };
    };

    searcher_AvoidX S(water.Points[floor(water.Points.size() / 2)]);
    S.search();

    V_netPp ps = S.getObjects();
    plot.SaveData(obj3D::extractX(ps), {{"ps", "2"}, {"pt", "7"}, {"title", "getObjects()"}});
    ps = S.getObjects_();
    plot.SaveData(obj3D::extractX(ps), {{"ps", "2"}, {"pt", "6"}, {"title", "getObjects_()"}});
    plot.SaveVectorData(getVectorData(S.getEnteredLines()), {{"arrowstyle", "1"}, {"title", "getEnteredLines()"}});
    plot.SaveVectorData(getVectorData(S.getReachedLines()), {{"arrowstyle", "2"}, {"title", "getReachedLines()"}});
    plot.SaveVectorData(getVectorData(S.getReachedLines_()), {{"arrowstyle", "3"}, {"title", "getReachedLines_()"}});

    plot.plot3d();
    plot.Clear();
    /*SEARCHER*/
    if (_stopall_)
      std::cin.ignore();
  }

  Print("SEARCHER DEPENDING ON NETWORK", green);
  //===================
  if (false || _runall_)
  { /*SEARCHER2*/
    NetworkObj obj("./obj/tank.obj");
    int i = 0;
    while (i++ < 2)
      for (const auto &l : obj3D::takeInsideOfBounds(obj.getLines(), {{-0.55, 0.55}, {-0.55, 0.55}, {-100, 100}}))
        longerLine(l)->divide();

    NetworkW water({2, 2}, {1.15, 1.15, 1 / 2.}, .1);
    Network corssNetwork({&obj, &water});

    for (const auto &f : water.Faces)
      if (f->intersectQ())
        for (const auto &routeP : f->getPointsCutFacesBehind())
          for (const auto &tri_ps : network::triangulate(routeP, f->getNormal()))
            networkFace tmp_f(&corssNetwork, link(tri_ps, &corssNetwork));

    class my_searcher : public searcher<netP>
    {
    public:
      my_searcher(netPp start_, bool TorF_ = true) : searcher<netP>(start_){};

      bool condEnterLine(const netP *p, const netLp l) override
      {
        /// 同じネットワークに関しては，干渉しないなら続ける．
        /// しかし，ネットワークが切り替わる場合は，探査を続ける．
        if ((!l->penetrateQ()) || l->getNetwork() != p->getNetwork())
          return true;
        else
          return false;
      };

      bool condGetObject(const netLp l, const netP *P) override
      {
        /// XPointである場合，別に保存する，
        if (!P->isXPoint())
          return true;
        else
          return false;
      };

      bool condKeepSearch(const netLp l, const netP *P) override
      {
        if (!l->penetrateQ())
          return true;
        else
          return false;
      };
    };

    GNUPLOT plot;
    plot.Set({{"key", ""}});
    my_searcher S(water.Points[floor(water.Points.size() / 2)]);
    S.search();
    V_netPp ps = S.getObjects();
    plot.SaveData(obj3D::extractX(ps), {{"ps", "2"}, {"pt", "7"}, {"lc", "'red'"}, {"title", "getObjects()"}});
    ps = S.getObjects_();
    plot.SaveData(obj3D::extractX(ps), {{"ps", "2"}, {"pt", "7"}, {"lc", "'blue'"}, {"title", "getObjects_()"}});
    plot.SaveVectorData(getVectorData(S.getEnteredLines()), {{"arrowstyle", "1"}, {"title", "getEnteredLines()"}});
    plot.SaveVectorData(getVectorData(S.getReachedLines()), {{"arrowstyle", "2"}, {"title", "getReachedLines()"}});
    plot.SaveVectorData(getVectorData(S.getReachedLines_()), {{"arrowstyle", "3"}, {"title", "getReachedLines_()"}});

    plot.SaveVectorData(getVectorData(obj.Faces), {{"arrowstyle", "4"}, {"title", "getReachedLines_()"}});

    plot.plot3d();
    plot.Clear();
    /*SEARCHER2*/
    if (_stopall_)
      std::cin.ignore();
  }

  Print("SEARCHER DEPENDING ON NETWORK", green);
  //===================
  if (false || _runall_)
  { /*searcher_of_faces*/
    NetworkObj obj("./obj/tank.obj");
    int i = 0;
    while (i++ < 2)
      for (const auto &l : obj3D::takeInsideOfBounds(obj.getLines(), {{-0.55, 0.55}, {-0.55, 0.55}, {-100, 100}}))
        longerLine(l)->divide();

    NetworkW water({2, 2}, {1.15, 1.15, 1 / 2.}, .1);
    Network corssNetwork({&obj, &water});

    for (const auto &f : water.Faces)
      if (f->intersectQ())
        for (const auto &routeP : f->getPointsCutFacesBehind())
          for (const auto &tri_ps : network::triangulate(routeP, f->getNormal()))
          {
            networkFace tmp_f(&corssNetwork, link(tri_ps, &corssNetwork));
          }

    class my_searcher : public searcher<netF>
    {
    public:
      my_searcher(netFp start_, bool TorF_ = true) : searcher<netF>(start_){};

      bool condEnterLine(const netF *p, const netLp l) override
      {
        if (true)
          return true;
        else
          return false;
      };
      bool condGetObject(const netLp l, const netF *P) override
      {
        if (true)
          return true;
        else
          return false;
      };
      bool condKeepSearch(const netLp l, const netF *P) override
      {
        if (true)
          return true;
        else
          return false;
      };
    };

    GNUPLOT plot;
    plot.Set({{"key", ""}});

    my_searcher S(water.Faces[floor(water.Faces.size() / 2)]);
    S.search();

    V_netFp fs = S.getObjects();

    plot.SaveVectorData(getVectorData(fs), {{"arrowstyle", "1"}, {"title", "S.getObjects()"}});

    fs = S.getObjects_();
    plot.SaveVectorData(getVectorData(fs), {{"arrowstyle", "2"}, {"title", "S.getObjects_()"}});

    plot.SaveVectorData(getVectorData(S.getEnteredLines()), {{"arrowstyle", "3"}, {"title", "getEnteredLines()"}});
    plot.SaveVectorData(getVectorData(S.getReachedLines()), {{"arrowstyle", "4"}, {"title", "getReachedLines()"}});
    plot.SaveVectorData(getVectorData(S.getReachedLines_()), {{"arrowstyle", "5"}, {"title", "getReachedLines_()"}});

    plot.SaveVectorData(getVectorData(obj.Faces), {{"arrowstyle", "4"}, {"title", "getReachedLines_()"}});

    plot.plot3d();
    plot.Clear();
    /*searcher_of_faces*/
    if (_stopall_)
      std::cin.ignore();
  }

  Print("getPointsCutLines", green);
  //=====================
  if (false || _runall_)
  {
    GNUPLOT plot;
    plot.Set({{"key", ""}});

    NetworkObj obj("./obj/tank.obj");
    NetworkW water({4, 4}, {2, 1., .4}, .1);
    Network XNetwork({&obj, &water});

    for (const auto &f : water.Faces)
    {
      if (f->intersectQ())
      {
        for (const auto &ps : f->getPointsCutLines())
          plot.SaveData(obj3D::extractX(ps), {{"w", "lp"}, {"ps", "1"}, {"pt", "7"}});
      }
    }

    plot.SaveVectorData(getVectorData(water.Faces), {{"arrowstyle", "1"}, {"title", "faces_water"}}); //interaction faces
    // plot.SaveVectorData(getVectorData(faces_obj),{{"arrowstyle","2"},{"title","faces_obj"}});//interaction faces
    // plot.SaveVectorData(getVectorData(S.getObjects()),{{"arrowstyle","3"},{"title","S.getObjects()"}});//interaction faces
    // plot.SaveVectorData(getVectorData(S_.getObjects()),{{"arrowstyle","4"},{"title","S_.getObjects()"}});//interaction faces

    plot.plot3d();
    plot.Clear();
    if (_stopall_)
      std::cin.ignore();
  }

  Print("check partitioning by networks", green);
  //=====================
  if (false || _runall_)
  {

    GNUPLOT plot;
    plot.Set({{"style arrow 2", "nohead lc 'web-blue' lw .5 dt 1"}});
    plot.Set({{"style arrow 1", "nohead lc 'red' lw .5 dt 1"}});

    //case 1
    //-----------------------------
    // plot.Set({{"key",""},{"zrange","[-40:55]"}});
    // NetworkObj obj("./obj/camel.obj");
    // plot.SaveVectorData(getVectorData(obj.Faces),{{"arrowstyle","2"},{"title","obj"}});
    // plot.plot3d();
    // std::cin.ignore();
    // plot.Clear();
    // V_d farxyz = {100.,0.,-100.};
    //-----------------------------

    //case 2
    //-----------------------------
    plot.Set({{"key", ""}, {"zrange", "[-0.:1.1]"}});
    NetworkObj obj("./obj/tank.obj");
    plot.SaveVectorData(getVectorData(obj.Faces), {{"arrowstyle", "1"}, {"title", "obj"}});
    plot.plot3d();
    std::cin.ignore();
    plot.Clear();
    V_d farxyz = {0., 0., .2};

    for (const auto &l : obj.getLines())
      longerLine(l)->divide();
    for (const auto &l : obj.getLines())
      longerLine(l)->divide();
    //-----------------------------

    //case 3
    //-----------------------------
    // plot.Set({{"key",""},{"xrange","[-100:100]"},{"yrange","[0:200]"},{"zrange","[-5:40]"}});
    // NetworkObj obj("./obj/oyama.obj");
    // obj.translate({-80.,0.,0.});
    // ////obj.rotate(M_PI/2.,{1.,0.,0.});
    // obj.translate({0.,10.,0.});

    // plot.SaveVectorData(getVectorData(obj.Faces),{{"arrowstyle","1"},{"title","obj"}});
    // ////plot.SaveVectorData(getVectorData(obj2.Faces),{{"arrowstyle","2"},{"title","obj"}});
    // plot.plot3d();
    // std::cin.ignore();
    // plot.Clear();
    // V_d farxyz = {1000.,100.,.2};
    // for(const auto& l:obj.getLines())
    //   longerLine(l)->divide();

    // NetworkObj water("./obj/oyama2.obj");
    // for(const auto& l:water.getLines())
    //   water.divide(longerLine(l));

    // water.translate({-40.,0.,5.});
    //-----------------------------

    VV_d gw = GaussianQuadratureWeights(5, 0., 1.);
    Print("GaussianQuadratureWeights", red);
    Print(gw, red);
    VV_d gwgw;
    for (const auto &tw0 : gw)
    {
      double a = tw0[0];
      for (const auto &tw1 : gw)
      {
        V_d v = tw1 * (1 - a);
        double b = v[0];
        gwgw.push_back({a, b, tw0[1] * v[1]});
      }
    }

    for (auto i = 0; i < 1000; i++)
    {
      plot.Set({{"view", "50., " + std::to_string(i / 5.) + ", 1., 1."}});
      //case1
      //-----------------------------
      // NetworkW water({2,2},{90.1223 , 100.1333, 20.33+15.*cos(2.*i*M_PI/20.)}, 0.*cos(2.*i*M_PI/100.));
      //-----------------------------

      //case2
      //-----------------------------
      NetworkW water({4, 4}, {1.1223, 1.1333, .5 + 0. * abs(sin(2. * i * M_PI / 100.))}, .3 * cos(2. * i * M_PI / 100.));
      //-----------------------------

      //case3
      //-----------------------------
      //      obj.translate({.1*cos(i*M_PI/100.),.1*sin(i*M_PI/100.),.1*cos(i*M_PI/100.)});
      //-----------------------------

      Network XNetwork({&obj, &water});

      V_netFp Faces(0);
      V_netPp Points(0);
      V_netFp faces_water, faces_obj;

      //--------------------------
      V_netPp XP(0);
      V_netFp XF(0);
      Print("干渉ネットワークの作成", red);
      {
        for (const auto &f : water.Faces)
          if (f->intersectQ())
            for (const auto &routeP : f->getPointsCutFacesBehind())
              for (const auto &tri_ps : network::triangulate(routeP, f->getNormal()))
                faces_water.emplace_back(new networkFace(&XNetwork, link(tri_ps, &XNetwork)));

        for (const auto &f : obj.Faces)
          if (f->intersectQ())
            for (const auto &routeP : f->getPointsCutFacesBehind())
              for (const auto &tri_ps : network::triangulate(routeP, f->getNormal()))
                faces_obj.emplace_back(new networkFace(&XNetwork, link(tri_ps, &XNetwork)));

        XF.insert(XF.end(), faces_water.begin(), faces_water.end());
        XF.insert(XF.end(), faces_obj.begin(), faces_obj.end());

        ///check
        plot.SaveVectorData(getVectorData(XF), {{"arrowstyle", "1"}, {"title", "XF"}});

        //作成した面を線に保存して，線から参照可能にする
        for (const auto &f : faces_water)
          for (const auto &l : f->Lines)
            network::add(l->Faces, f);

        for (const auto &f : faces_obj)
          for (const auto &l : f->Lines)
            network::add(l->Faces, f);

        XP = XNetwork.Points;

        plot.SaveData(obj3D::extractX(XP), {{"ps", "2"}, {"pt", "7"}, {"lc", "'orange'"}, {"title", "XP"}});
      }
      Points.insert(Points.end(), XP.begin(), XP.end());
      //--------------------------
      V_netFp InsideF(0);
      Print("干渉のない内部の面を抽出", red);
      {
        searcher8<netF> S(network::takeNearest(water.Faces, farxyz));
        S.search();
        auto insideF = S.getObjects();
        //
        searcher8<netF> S_(network::takeNearest(obj.Faces, farxyz));
        S_.search();
        auto insideF_ = S_.getObjects();

        InsideF.insert(InsideF.end(), insideF.begin(), insideF.end());
        InsideF.insert(InsideF.end(), insideF_.begin(), insideF_.end());

        ///check
        plot.SaveVectorData(getVectorData(insideF), {{"arrowstyle", "2"}, {"title", "InsideF"}});   //interaction faces
        plot.SaveVectorData(getVectorData(insideF_), {{"arrowstyle", "5"}, {"title", "InsideF_"}}); //interaction faces
      }
      Faces.insert(Faces.end(), InsideF.begin(), InsideF.end());
      //--------------------------
      V_netPp InsideP(0);
      Print("干渉のない内部の点を抽出", red);
      {
        inside_point_searcher S(network::takeNearest(water.Points, farxyz));
        S.search();
        auto insideP = S.getObjects();
        inside_point_searcher S_(network::takeNearest(obj.Points, farxyz));
        S_.search();
        auto insideP_ = S_.getObjects();

        InsideP.insert(InsideP.end(), insideP.begin(), insideP.end());
        InsideP.insert(InsideP.end(), insideP_.begin(), insideP_.end());

        ///check
        plot.SaveData(S.getObjects()[0]->getX(), {{"ps", "3"}, {"pt", "7"}, {"lc", "'purple'"}, {"title", "InsideP water start"}});
        plot.SaveData(S_.getObjects()[0]->getX(), {{"ps", "3"}, {"pt", "7"}, {"lc", "'purple'"}, {"title", "InsideP obj start"}});

        plot.SaveVectorData(getVectorData(S_.getEnteredLines()), {{"arrowstyle", "2"}, {"title", "InsideF"}}); //interaction faces
        plot.SaveVectorData(getVectorData(S_.getReachedLines()), {{"arrowstyle", "3"}, {"title", "InsideF"}}); //interaction faces
        plot.SaveData(obj3D::extractX(insideP), {{"ps", ".5"}, {"pt", "7"}, {"lc", "'web-blue'"}, {"title", "insideP"}});
        plot.SaveData(obj3D::extractX(insideP_), {{"ps", ".5"}, {"pt", "7"}, {"lc", "'web-green'"}, {"title", "insideP_"}});
      }
      Points.insert(Points.end(), InsideP.begin(), InsideP.end());
      //--------------------------
      V_netPp EndP(0);
      Print("外部の点を抽出", red);
      {
        end_point_searcher eS(network::takeNearest(water.Points, farxyz));
        eS.search(false);
        auto endPs = eS.getObjects();
        end_point_searcher eS_(network::takeNearest(obj.Points, farxyz));
        eS_.search(false);
        auto endPs_ = eS_.getObjects();

        EndP.insert(EndP.end(), endPs.begin(), endPs.end());
        EndP.insert(EndP.end(), endPs_.begin(), endPs_.end());

        plot.SaveData(obj3D::extractX(EndP), {{"ps", "2"}, {"pt", "9"}, {"lc", "'forest-green'"}, {"title", "EndP"}});
        //-------------------

        //~~~~~~~~~~~~~~
        class depth_searcher : public searcher<networkPoint>
        {
        public:
          bool didGetOne;
          depth_searcher() : searcher<networkPoint>(), didGetOne(false){};
          bool condEnterLine(const networkPoint *p, const netLp l) override
          {
            if (this->depth < 5)
            {
              return true;
            }
            else
              return false;
          };
        };

        Print("EndPのVの見積もり", red);
        if (false)
        {
          int i = 0;
          for (const auto &p : EndP)
          {
            i += 1;
            depth_searcher S;
            S.set(p);
            S.addNetwork(&XNetwork);
            S.addNetwork(&water);
            S.addNetwork(&obj);
            Print(S.getNetworks(), red);
            S.search();

            plot.SaveData(obj3D::extractX(S.getObjects()), {{"ps", "1"}, {"pt", std::to_string(i % 10)}, {"lc", plot.rgb((double)((2 * i) % 201), (double)((2 * i) % 102), (double)(i % 53))}, {"title", std::to_string(i)}});
            plot.SaveVectorData(getVectorData(S.getEnteredLines()), {{"lc", plot.rgb((double)((2 * i) % 201), (double)((2 * i) % 102), (double)(i % 53))}, {"title", std::to_string(i)}});
          }
        }
        //~~~~~~~~~~~~~~

        //~~~~~~~~~~~~~~
        Print("面はちゃんと作られているか？", red); //ガウス点でチェック
        if (false)
        {
          int i = 0;
          for (const auto &p : EndP)
          {
            i++;
            VV_d XYZ;
            for (const auto &f : p->getFaces_intersectQ(false))
            {
              auto xyz = f->getLocations();
              for (const auto &ab : gwgw)
              {
                XYZ.push_back(Dot(BEM::N(ab), xyz));
              }
            }
            plot.SaveData(XYZ, {{"ps", "1"}, {"pt", std::to_string(i % 3)}, {"lc", plot.rgb((double)(i % 201), (double)(i % 102), (double)(i % 53))}, {"notitle", ""}});
          }
        }
        //~~~~~~~~~~~~~~
      }
      //--------------------------
      Print("周辺の面を点が参照できているか確認", red);
      {
        auto s = i % ((int)water.Points.size());
        plot.SaveData(water.Points[s]->getX(), {{"ps", "2.5"}, {"pt", "7"}, {"lc", "'blue'"}, {"title", "water.Points[i]"}});
        plot.SaveVectorData(getVectorData(water.Points[s]->getFaces_intersectQ(false)), {{"lc", "'blue'"}, {"lw", "1.5"}, {"title", "getFaces_intersectQ(false)"}});

        int i = 0;
        Print(water.Points[s]->getNeighborsSort());
        Print(water.Points[s]->getNeighbors(), red);
        for (const auto &p : water.Points[s]->getNeighborsSort())
        {
          i++;
          plot.SaveData(p->getX(), {{"pt", "7"}, {"ps", std::to_string(i / 3.)}, {"w", "p"}, {"title", std::to_string(i)}});
        }

        s = i % ((int)XNetwork.Points.size());
        plot.SaveData(XNetwork.Points[s]->getX(), {{"ps", "2.5"}, {"pt", "7"}, {"lc", "'red'"}, {"title", "XNetwork.Points[i]"}});
        plot.SaveVectorData(getVectorData(XNetwork.Points[s]->getFaces_intersectQ(false)), {{"lc", "'red'"}, {"lw", "1.5"}, {"title", "getFaces_intersectQ(false)"}});

        for (const auto &p : XNetwork.Points[s]->getNeighborsSort())
        {
          i++;
          plot.SaveData(p->getX(), {{"pt", "7"}, {"ps", std::to_string(i / 3.)}, {"w", "p"}, {"title", ""}});
        }

        //立体角計算
        for (const auto &q : XP)
        {
          auto viewp = q->getX();
          auto ps = q->getNeighborsSort();
          geometry::polygon poly(obj3D::extractX(ps));
          for (const auto &ind : poly.triangulate(viewp))
          {
            plot.SaveVectorData({{ps[ind[0]]->getX(), ps[ind[1]]->getX() - ps[ind[0]]->getX()},
                                 {ps[ind[1]]->getX(), ps[ind[2]]->getX() - ps[ind[1]]->getX()},
                                 {ps[ind[2]]->getX(), ps[ind[0]]->getX() - ps[ind[2]]->getX()}},
                                {{"lw", "3"}, {"title", std::to_string(q->getSolidAngle(true))}});
          }
        }

        // double SolidAngle(){
        //   auto viewp = this->getX();
        //   auto ps = this->getNeighborsSort();
        //   geometry::polygon poly(obj3D::extractX(ps));
        //   double ret(0);
        //   for(const auto& ind:poly.triangulate(viewp)){
        //     ret += geometry::SolidAngle(viewp, ps[ind[0]]->getX(), ps[ind[1]]->getX(), ps[ind[2]]->getX() );
        //   }
        //   return ret;
        // };

        // s = i%((int)obj.Points.size());
        // plot.SaveData(obj.Points[s]->getX(),{{"ps","2.5"},{"pt","7"},{"lc","'green'"},{"title","obj.Points[i]"}});
        // plot.SaveVectorData(getVectorData(obj.Points[s]->getFaces_intersectQ(false)),{{"lc","'green'"},{"lw","1.5"},{"title","getFaces_intersectQ(false)"}});
      }

      Print("削除のチェック", red);
      if (true)
      {
        for (auto &f : faces_water)
          f->Delete();
        for (auto &f : faces_obj)
          f->Delete();

        // searcher S(network::takeNearest(water.Points, {0,0,0.2}));
        // S.search();
        // auto insidePs = S.getObjects();
        // plot.SaveData(obj3D::extractX(insidePs), {{"ps","2"},{"pt","7"},{"lc","'pink'"},{"title","getObjects()"}});
      }

      Print("プロット", red);
      plot.plot3d();
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
      plot.Clear();
      if (_stopall_)
        std::cin.ignore();
    }
  }

  return 0;
};
