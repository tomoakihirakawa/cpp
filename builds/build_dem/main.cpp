int time_step;
#define simulation
#include "GNUPLOT.hpp"
#define DEM
#include "Network.hpp"
std::string home_dir = std::getenv("HOME");
/////////////////////////////////////
V_netLp getLinesAround(netPp p)
{
  V_netLp ret = {};
  for (const auto &f : p->getFaces())
    for (const auto &line : f->getLines())
      ret.emplace_back(line);
  return DeleteDuplicates(ret);
};
V_netLp getLinesAround(netLp l)
{
  V_netLp ret = {l};
  for (const auto &p : l->getPoints())
    for (const auto &f : p->getFaces())
      for (const auto &line : f->getLines())
        ret.emplace_back(line);
  return DeleteDuplicates(ret);
};
bool refine(netLp l, double len)
{
  if (l->length() > len)
  {
    l->divide();
    // Print("新しいlineが生成されることは，このループに問題を引き起こさないか？");
    auto tmp = getLinesAround(l);
    for (auto &line : tmp)
    {
      if (isFlat(line))
      {
        line->flipIfIllegal();
        auto fs = line->getFaces();
        for (auto i = 0; i < 2; i++)
          for (auto &f : fs)
            LaplacianSmoothingIfFlat(f->getPoints() /*内部でシャッフルする*/); //ちゃんとsetXされているかチェック
      }                                                                        // LaplacianSmoothingIfFlat(line->getPoints() /*内部でシャッフルする*/);
    }
    return true;
  }
  return false;
};
/////////////////////////////////
void remesh(const V_Netp &all_obj, const double lim_len)
{
  for (auto i = 0; i < all_obj.size(); i++)
  {
    auto obj0 = all_obj[i];
    display(obj0);
    bool found;

    //均等なメッシュを作りたいなら，limに該当する線の数を先にチェックし，全体における，
    //該当する線の割合をlimで調整しながら進めるべき．
    //長さが大きく異なる線を，適当に分割すると，見た目，構造が大きく変わってしまう．
    //構造を維持し，徐々に変化させる必要がある．急な変化はよくない
    // 該当する点の割合と調整し徐々に分割していく

    int count = 0;
    try
    {
      Print(extLength(obj0->getLines()));
      for (const auto &lim : Subdivide(2. * lim_len, lim_len, 50))
      {
        int c = 0;
        do
        {
          found = false;
          auto line = obj0->getLines();
          std::shuffle(std::begin(line), std::end(line), std::default_random_engine());
          //ループ中にlの配置が変わるので，全体をフリップできない場合がある
          // LaplacianSmoothingIfFlat(obj0->getPoints());
          int cc = 0;
          for (auto j = 0; j < line.size(); j++)
            found = refine(line[j], lim);
          // LaplacianSmoothingIfFlat(obj0->getPoints());
          // std::cout << Red << "c = " << c++ << reset << std::endl;
          c++;
        } while (c <= 5 /*最低回数*/ || (found && c < 100));
        mk_vtu(home_dir + "/vtu/remesh" + std::to_string(count++) + ".vtu", obj0->getFaces());
      }
    }
    catch (error_message &e)
    {
      std::cerr << e.what() << reset << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
    }

    // try
    // {
    //   for (const auto &lim : Reverse(Subdivide(2. * lim_len, lim_len/2., 50)))
    //   {
    //     int c = 0;
    //     do
    //     {
    //       found = false;
    //       auto line = obj0->getLines();
    //       // std::shuffle(std::begin(line), std::end(line), std::default_random_engine());
    //       std::shuffle(std::begin(line), std::end(line), std::default_random_engine());
    //       //ループ中にlの配置が変わるので，全体をフリップできない場合がある
    //       // LaplacianSmoothingIfFlat(obj0->getPoints());
    //       for (auto j = 0; j < line.size(); j++)
    //       {
    //         // std::cout << "line[j] = " << line[j] << ", lim = " << lim << std::endl;
    //         // std::cout << "line[j]->length() = " << line[j]->length() << std::endl;
    //         found = coarsen(line[j], lim);
    //         if (found)
    //           break;
    //       }
    //       // LaplacianSmoothingIfFlat(obj0->getPoints());
    //       // std::cout << Red << "c = " << c++ << reset << std::endl;
    //       c++;
    //     } while (c < 200);
    //     mk_vtu(home_dir+"/vtu/remesh" + std::to_string(count++) + ".vtu", obj0->getFaces());
    //   }
    // }
    // catch (error_message &e)
    // {
    //   std::cerr << e.what() << reset << std::endl;
    //   throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
    // }
  }
  Print("remesh done", Red);
};
////////////////////////////////////////////////////////////
bool areContacting(const netPp p, const netPp q)
{
  return ((Norm3d(p->getX() - q->getX()) - (p->radius + q->radius)) < 0.);
};
void addIfContact(const netPp p, const netPp q)
{

  if (areContacting(p, q))
  {
    network::add(p->contactP, q);
    network::add(q->contactP, p);
  };
};

bool linkIfCloser_findContactP(const netPp p, const netPp q, const V_i &upperlower = {5, 8})
{
  //メンバーと比べて距離が遠い場合のみfalseを返す
  // int lower_lim_vectorlen = upperlower[0];
  // int upper_lim_vectorlen = upperlower[1];
  if (p != q)
  {

    addIfContact(p, q);

    auto neighbors = p->getNeighbors();
    // auto neighbors = Flatten(BFS(p,1,{p->getNetwork()}));
    if (neighbors.size() < upperlower[0])
    {
      // Print("link");
      link(p, q, p->getNetwork());
      p->sortLinesByLength();
      return true;
    }
    else if (!MemberQ(neighbors, q))
    {
      /*current farest dist*/
      netLp longestL = (*(p->getLines().rbegin()));
      if (longestL->length() > Norm(p->getX() - q->getX()))
      {
        auto replaced_P = (*longestL)(p);
        longestL->Replace(replaced_P, q);
        p->sortLinesByLength();
        q->sortLinesByLength();

        // delete if too many points
        {
          longestL = *(p->getLines().rbegin());
          if (p->getLines().size() > upperlower[1] && (*longestL)(p)->getLines().size() > upperlower[1])
            delete longestL;
        }

        // delete if too many points
        {
          longestL = *(q->getLines().rbegin());
          if (q->getLines().size() > upperlower[1] && (*longestL)(q)->getLines().size() > upperlower[1])
            delete longestL;
        }

        if (replaced_P->getLines().size() < upperlower[0])
          for (const auto &n : neighbors)
            linkIfCloser_findContactP(replaced_P, n);

        return true;
      }
      else
      {
        //メンバーではなかったpointだったが，メンバーの中で最遠の点よりも遠かった．
        return false;
      }
    }
    else
      return true;
  }
  else
    return true;
};

//////////////////////////////////////////////////////////////////
V_netPp neighborsSort(const netPp p, const int i)
{
  auto neighbors = TakeExcept(Flatten(BFS(p, i, {p->getNetwork()})), p);
  auto X = p->getX();
  std::sort(neighbors.begin(), neighbors.end(), [&X](const netPp p0, const netPp p1)
            { return (Norm3d(p0->getX() - X) < Norm3d(p1->getX() - X)); });
  return neighbors;
};

V_netPp neighborsSort(const netPp p, const V_netPp &neighbors_IN)
{
  //与えられた点を自身に近い順に並び替えたものを
  auto X = p->getX();
  auto neighbors = TakeExcept(neighbors_IN, p);
  std::sort(neighbors.begin(), neighbors.end(),
            [&X](const netPp p0, const netPp p1)
            {
              return (Norm3d(p0->getX() - X) < Norm3d(p1->getX() - X));
            });
  return neighbors;
};
//////////////////////////////////////////////////////////////////
// void test(const V_netPp &points)
// {
//   auto start = std::chrono::high_resolution_clock::now();
//   int depth = 2;
//   V_d upperlower = {10, 20};
// #ifdef _OPENMP
//   std::cout << "並列化．BFSを使って各pointのNeighborsを取得する．depth = " << depth << std::endl;
// #pragma omp parallel for
// #endif
//   for (auto i = 0; i < points.size() /*N*/; i++)
//   {
//     auto p = points[i];
//     auto X = p->getX();
//     auto bfs = Flatten(BFS(p, i, {p->getNetwork()}));
//     std::sort(bfs.begin(), bfs.end(), [&X](const netPp p0, const netPp p1) { return (Norm3d(p0->getX() - X) < Norm3d(p1->getX() - X)); });
//     for (const auto &q : bfs)
//     { //メンバーと比べて距離が遠い場合のみfalseを返す
//       auto neighbors = p->getNeighbors();
//       if (p != q && !MemberQ(neighbors, q))
//       {
//         if (areContacting(p, q))
//         {
// #ifdef _OPENMP
// #pragma omp critical
// #endif
//           {
//             network::add(p->contactP, q);
//             network::add(q->contactP, p);
//           }
//         };

//         if (neighbors.size() < upperlower[0])
//         {
// #ifdef _OPENMP
// #pragma omp critical
// #endif
//           {
//             // Print("link");
//             link(p, q, p->getNetwork());
//             p->sortLinesByLength();
//           }
//         }
//         else if (!MemberQ(neighbors, q))
//         {
//           /*current farest dist*/
//           netLp longestL = (*(p->getLines().rbegin()));
//           if (longestL->length() > Norm(p->getX() - q->getX()))
//           {
// #ifdef _OPENMP
// #pragma omp critical
// #endif

//             {
//               auto replaced_P = (*longestL)(p);
//               longestL->Replace(replaced_P, q);
//               p->sortLinesByLength();
//               q->sortLinesByLength();

//               //delete if too many points
//               {
//                 longestL = *(p->getLines().rbegin());
//                 if (p->getLines().size() > upperlower[1] && (*longestL)(p)->getLines().size() > upperlower[1])
//                   delete longestL;
//               }

//               //delete if too many points
//               {
//                 longestL = *(q->getLines().rbegin());
//                 if (q->getLines().size() > upperlower[1] && (*longestL)(q)->getLines().size() > upperlower[1])
//                   delete longestL;
//               }

//               if (replaced_P->getLines().size() < upperlower[0])
//                 for (const auto &n : neighbors)
//                   linkIfCloser_findContactP(replaced_P, n);
//             }
//           }
//         }
//       }
//     }
//   }
//   auto finish = std::chrono::high_resolution_clock::now();
//   std::chrono::duration<double> elapsed = finish - start;
//   std::cout << Green << "Elapsed time: " << Red << elapsed.count() << reset << " s\n";
// };
//////////////////////////////////////////////////////////////////
V_d outside_direction(const netPp p)
{
  double sum = 0., weight;
  V_d Vtothis, ret(3, 0.);
  for (const auto &q : TakeExcept(Flatten(BFS(p, 2, {p->getNetwork()})), p))
  {
    Vtothis = p->getX() - q->getX();
    weight = 1. / Norm3d(Vtothis);
    sum += weight;
    ret += weight * Vtothis;
  }
  return ret / sum;
};
//////////////////////////////////////////////////////////////////

// void calculateForceFromContactP(const V_netPp &points)
// {
//   try
//   {
//     double k = 30.;
//     double dump = 1.;
// #ifdef _OPENMP
//     Print("並列化@calculateForceFromContactP");
// #pragma omp parallel for
// #endif
//     for (auto i = 0; i < points.size(); i++)
//     { // double radius = p->radius;
//       auto p = points[i];
//       V_d Vtothis(3, 0.);
//       double del;
//       V_d X = p->getX();
// #ifdef _OPENMP
// #pragma omp critical
// #endif
//       {
//         V_netPp ps = p->contactP;
//         for (const auto &q : ps)
//         {
//           network::erase(p->contactP, q);
//           network::erase(q->contactP, p);
//           Vtothis = X - q->getX();
//           del = Norm3d(Vtothis) - (p->radius + q->radius);
//           if (del < 0.)
//           {
//             V_d tmp_force = (k * pow(-del, 1.5) - dump * Dot(p->velocity - q->velocity /*dir to this*/, Vtothis)) * Vtothis;
//             points[i]->force += tmp_force;
//             q->force -= tmp_force;
//           }
//         }
//       }
//     }
//   }
//   catch (error_message &e)
//   {
//     e.print();
//     throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
//   }
// };

// void calculateActionAndReaction(const V_netPp &points, const V_netFp &faces)
// {
//   try
//   {
//     double k = 30. * 2.;
//     double dump = 1.;
// #ifdef _OPENMP
//     Print("並列化@calculateActionAndReaction");
// #pragma omp parallel for
// #endif
//     for (auto i = 0; i < points.size(); i++)
//     {
//       auto p = points[i];
//       V_d Vtothis;
//       V_d X = p->getX();
//       double del;
//       for (const auto &f : faces)
//       {
//         Vtothis = Dot(X - f->getX(), f->getNormal()) * f->getNormal();
//         // normalがX-Xfと逆の場合，normalをかけると符号がかわり，Xを向く
//         del = p->radius - Norm3d(Vtothis);
//         //まだヒットしたかわからない
//         if (del > 0.)
//         {
//           V_d normal = f->getNormal();
//           VV_d AB = {X + p->radius * normal, X - p->radius * normal};
//           if (isIntersectingSurface(f->getLocations(), AB))
//           {
//             V_d tmp_force = Vtothis * (k * pow(del, 2.) - dump * Dot(p->velocity /*- q->velocity*/ /*dir to this*/, Vtothis));
// #ifdef _OPENMP
// #pragma omp critical
// #endif
//             {
//               points[i]->force += tmp_force;
//               f->force -= tmp_force;
//             }
//           }
//         }
//       }
//     }
//   }
//   catch (error_message &e)
//   {
//     e.print();
//     throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
//   }
// };

// V_d force(const netPp p, const netFp f)
// {
//   // double radius = p->radius;
//   double k = 30.;
//   double dump = 1.;
//   V_d ret(3, 0.);
//   V_d Vtothis;
//   double del;b
//   V_d X = p->getX();
//   for (const auto &q : p->getNeighbors())
//   {
//     Vtothis = Dot(X - f->getX(), f->getNormal()) * f->getNormal();
//     // normalがX-Xfと逆の場合，normalをかけると符号がかわり，Xを向く
//     del = p->radius - Norm3d(Vtothis);
//     //まだヒットしたかわからない
//     if (del > 0.)
//     {
//       auto normal = f->getNormal();
//       auto AB = {X + p->radius * normal, X - p->radius * normal};
//       if (isIntersectingSurface(f->getLocations(), AB))
//       {
//         ret += (k * pow(del, 1.5) - dump * Dot(p->velocity /*- q->velocity*/ /*dir to this*/, Vtothis)) * Vtothis;
//       }
//     }
//   }
//   return ret;
// };

// V_d force(const netFp f, const netPp p)
// {
//   return -force(p, f);
// };

///////////////////////////////////////////////////
double dt = 0.01 / 2.;
/* ------------------------------------------------------ */
void update(const std::unordered_set<networkPoint *> &points, const std::vector<V_netFp> &faces = {})
{
  // Print("baloon facesとの反射を考慮したpointsの時間発展");
  VVV_d p0p1p2s;
  for (const auto &fs : faces)
    for (const auto &f : fs)
      p0p1p2s.emplace_back(obj3D::extractX(f->getPoints()));

  // 自分の周り以外との接触回避
  // #ifdef _OPENMP
  //   Print("並列化");
  // #pragma omp parallel for
  // #endif

  for (auto it = points.begin(); it != points.end(); ++it)
  {
    auto p = *it;
    auto X = p->getX();
    bool ishit = false;
    V_d V = {std::get<0>(p->V), std::get<1>(p->V), std::get<2>(p->V)};
    V_d F = {std::get<0>(p->F), std::get<1>(p->F), std::get<2>(p->F)};
    V_d dx = V * dt;
    VV_d v_Vnew = {V + (F / p->mass) * dt};
    VV_d v_dx = {};
    // int hitIndex = -1;
    V_i hitIndices = {};
    int count = 0;
    do
    {
      count++;
      intersectionTriangleLine LT(p0p1p2s, X, X + dx, hitIndices);
      ishit = LT.isIntersect;
      if (ishit)
      {
        hitIndices.emplace_back(LT.indexOfTriangle);
        X = LT.X;        // on wall
        dx = LT.vecX2B_; //残りdx
        v_dx.emplace_back(LT.vecA2X);
        auto Vnew = LT.reflectIfPossible(*v_Vnew.rbegin());
        v_Vnew.emplace_back(Vnew);
      }
      else
      {
        v_dx.emplace_back(dx);
      }
    } while (ishit || count > 100);

#ifdef _OPENMP
#pragma omp critical
#endif
    {
      p->setX(p->getX() + Total(v_dx));
      std::get<0>(p->V) = (*v_Vnew.rbegin())[0];
      std::get<1>(p->V) = (*v_Vnew.rbegin())[1];
      std::get<2>(p->V) = (*v_Vnew.rbegin())[2];
    }
  }
};

// void update(const V_netPp &points)
// {
//   Print("baloon facesとの反射を考慮したpointsの時間発展");
// #ifdef _OPENMP
//   Print("並列化");
// #pragma omp parallel for
// #endif
//   for (auto i = 0; i < points.size(); i++)
//   {
//     auto p = points[i];
//     p->velocity = Take(p->velocity, {0, 3}) + Take(p->force, {0, 3}) / p->mass * dt;
//     p->setX(p->getX() + p->velocity * dt);
//   }
// };
/* ------------------------------------------------------ */
// void update(Network *net)
// {
//   net->calcPhysicalProperties();
//   std::cout << Grid({{"velocity", net->V},
//                      {"force", net->F},
//                      {"inertia", net->I},
//                      {"center_of_mass", net->center_of_mass}},
//                     70)
//             << std::endl;
//   // std::cin.ignore();
//   net->velocity += V_d(net->force / net->inertia) * dt;
//   for (const auto &p : net->getPoints())
//   {
//     p->setX(p->getX() + V_d{net->velocity[0], net->velocity[1], net->velocity[2]} * dt);
//     p->velocity = net->velocity;
//   }
//   //
//   Quaternion Q(1., 0., 0., 0.);
//   Q.set(Q() + Q.d_dt(V_d{net->velocity[3], net->velocity[4], net->velocity[5]} * dt));
//   net->rotate(Q, net->center_of_mass);
// };
//////////////////////////////////////////////////////////////////
double g = 9.81;
////////////////////////
void calcForceOnPoints(const Buckets<networkPoint> &B, const V_netPp &points)
{
  double k = 30.;
  double dump = 10.;

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (auto i = 0; i < points.size(); i++)
  {
    auto p = points[i];
    // 重力で初期化するようにしなければならない
    std::get<0>(p->force) = 0.;
    std::get<1>(p->force) = 0.;
    std::get<2>(p->force) = -g * p->mass;
    T6d Vtothis({0., 0., 0., 0., 0., 0.});
    Tddd X;
    double del;
    for (const auto &q : Flatten(B.getObjects(p->getXtuple(), 2. * p->radius /*depth*/)))
      if ((Norm(p->getXtuple() - q->getXtuple()) - (p->radius + q->radius)) < 0.)
        if (p != q)
        {
          X = p->getXtuple() - q->getXtuple();
          Vtothis = {std::get<0>(X), std::get<1>(X), std::get<2>(X), 0., 0., 0.};
          del = Norm(X) - (p->radius + q->radius);
          p->force += (k * pow(-del, 1.5) - dump * Dot((p->V - q->V) /*dir to this*/, Vtothis)) * Vtothis;
        }
  }
};
/* ------------------------------------------------------ */
/*                    面との衝突の計算                      */
/* ------------------------------------------------------ */
networkFace *findClosest(const V_netFp &fs, const geometry::Sphere &sphere)
{
  double closest = 1E+20;
  networkFace *closestF = nullptr;
  for (const auto &f : fs)
  {
    auto ixn = intersection(sphere, geometry::Triangle(extractX(f)));
    if (ixn.isIntersecting)
    {
      if (ixn.distance < closest)
      {
        closest = ixn.distance;
        closestF = f;
      }
    }
  }
  return closestF;
};
/* ------------------------------------------------------ */
map_P_Vd P_hitvector, P_hit;
void calcForce_interaction_PointsAndContactFaces(const std::unordered_set<networkPoint *> &points)
{
  //! 点Pointsとその点が接触した面ContactFacesに作用反作用の力を足し合わせる
  //
  // うさぎなしで計算し，壁と粒子との衝突がうまくいっているか確かめる．
  double k = 30000.;
  double dump = 2000.;
  //
  P_hitvector.clear();
  P_hit.clear();

  // for (auto i = 0; i < points.size(); i++)
  //   points[i]->setStatus(false);

  for (const auto &p : points)
    p->setStatus(false);

#ifdef _OPENMP
    // #pragma omp parallel for
#pragma omp parallel
#endif
  for (auto it = points.begin(); it != points.end(); ++it)
  {
#pragma omp single nowait
    {
      auto p = *it;
      Tddd X;
      double del;
      Tddd FORCE = {0., 0., 0.};
      int count = 0;
      for (const auto &f : p->getContactFaces())
      {
        auto ixn = intersection(geometry::Sphere(p->getXtuple(), p->radius),
                                geometry::Triangle(extractXtuple(f)));
        if (ixn.isIntersecting)
        {
          auto dx = (p->radius - ixn.distance);
          if (dx > 0.)
          {
            Tddd v = (p->getXtuple() - ixn.X);
            Tddd force = k * dx * v + k * std::pow(dx, 3) * v;
            force -= dump * Dot(p->Vxyz(), v) * v;
            FORCE += force;
            //% ------------------------------------------------------ */
            //%                      面に対する反作用                     */
            //% ------------------------------------------------------ */
            std::get<0>(f->force) -= std::get<0>(force);
            std::get<1>(f->force) -= std::get<1>(force);
            std::get<2>(f->force) -= std::get<2>(force);
            /* ------------------------------------------------------ */
            count++;
          }
        }
      }
      if (count > 0)
      {
        std::get<0>(p->force) += std::get<0>(FORCE) / (double)count;
        std::get<1>(p->force) += std::get<1>(FORCE) / (double)count;
        std::get<2>(p->force) += std::get<2>(FORCE) / (double)count;
      }
    }
  }
};
/* ------------------------------------------------------ */
// void fun(networkPoint *p, networkPoint *q)
// {
//   double k = 30.;
//   double dump = 10.;
//   if (p != q)
//     if ((Norm3d(p->getX() - q->getX()) - (p->radius + q->radius)) < 0.)
//     {
//       auto Vtothis = p->getX() - q->getX();
//       auto del = Norm3d(Vtothis) - (p->radius + q->radius);
//       p->force += (k * pow(-del, 1.5) - dump * Dot(p->velocity - q->velocity /*dir to this*/, Vtothis)) * Vtothis;
//     }
// };
// auto apply = [](Buckets<networkPoint> B, const V_netPp &points)
// {
// #ifdef _OPENMP
// #pragma omp parallel for
// #endif
//   for (auto i = 0; i < points.size(); i++)
//   {
//     auto p = points[i];
//     p->force = {0., 0., -g * p->mass};
//     B.apply([p](networkPoint *q)
//             { fun(p, q); },
//             p->getX(), 2. * p->radius);
//   }
// };
/////////////////////////////////////////////////////////
using uomap_P_Vd = std::unordered_map<networkPoint *, std::vector<double>>;
using uomap_P_Tddd = std::unordered_map<networkPoint *, Tddd>;
using uomap_P_d = std::unordered_map<networkPoint *, double>;
//
using map_P_Tddd = std::map<networkPoint *, Tddd>;
using map_P_d = std::map<networkPoint *, double>;
using VV_SorIorMap = std::vector<std::vector<std::variant<std::string,
                                                          int,
                                                          map_P_Vd,
                                                          map_P_Tddd,
                                                          map_P_d>>>;

void output_points(Network const *net, std::string const &name)
{
  uomap_P_Tddd P_V, P_F;
  uomap_P_d P_r, P_hit;
  for (const auto &p : net->getPoints())
  {
    geometry::Sphere sphere(p->getXtuple(), 2.);
    // for (const auto &f : obj->getFaces())
    // {
    //   // ベクトル方向がおかしい
    //   auto ixn = intersection(sphere, geometry::Triangle(extractX(f)));
    //   if (ixn.isIntersecting)
    //   {
    //     P_hit[p] = {ixn.distance};
    //     P_hitvector[p] = {std::get<0>(ixn.X - sphere.X), std::get<1>(ixn.X - sphere.X), std::get<2>(ixn.X - sphere.X)};
    //     break;
    //   }
    // }
    P_hit[p] = (double)p->getStatus();
    P_V[p] = Tddd{std::get<0>(p->V), std::get<1>(p->V), std::get<2>(p->V)};
    P_F[p] = Tddd{std::get<0>(p->F), std::get<1>(p->F), std::get<2>(p->F)};
    P_r[p] = p->radius;
  }
  // VV_SorIorMap data = ;
  std::unordered_set<networkPoint *> points = net->getPoints();
  std::string outname = home_dir + "/vtu/points_" + name + ".vtu";
  mk_vtu(outname,
         {points},
         {{"isIntersecting", P_hit}, {"velocity", P_V}, {"force", P_F}, {"radius", P_r}});
};
//* ------------------------------------------------------ */
void output_faces(Network const *net, std::string const &name)
{
  uomap_P_Tddd P_V, P_F;
  uomap_P_d P_r, P_hit;
  for (const auto &p : net->getPoints())
  {
    geometry::Sphere sphere(p->getXtuple(), 2.);
    // for (const auto &f : obj->getFaces())
    // {
    //   // ベクトル方向がおかしい
    //   auto ixn = intersection(sphere, geometry::Triangle(extractX(f)));
    //   if (ixn.isIntersecting)
    //   {
    //     P_hit[p] = {ixn.distance};
    //     P_hitvector[p] = {std::get<0>(ixn.X - sphere.X), std::get<1>(ixn.X - sphere.X), std::get<2>(ixn.X - sphere.X)};
    //     break;
    //   }
    // }
    P_hit[p] = (double)p->getStatus();
    P_V[p] = Tddd{std::get<0>(p->V), std::get<1>(p->V), std::get<2>(p->V)};
    P_F[p] = Tddd{std::get<0>(p->F), std::get<1>(p->F), std::get<2>(p->F)};
    P_r[p] = p->radius;
  }
  mk_vtu(home_dir + "/vtu/faces_" + name + ".vtu",
         net->getFaces(),
         {{"isIntersecting", P_hit},
          // {"hitvector", 3, P_hitvector},
          {"velocity", P_V},
          {"force", P_F},
          {"radius", P_r}});
};
//* ------------------------------------------------------ */
//! ------------------------------------------------------ */
//!                         点 <-> 点                       */
//%                         面 <-> 点                       */
//% ------------------------------------------------------ */
#define debug_addContact
void addContact(const V_Netp &nets, double radius)
{
#if defined(debug_addContact)
  TimeWatch watch;
  std::cout << "クリア" << std::endl;
#endif
  for (const auto &n : nets)
  {
    for (const auto &p : n->getPoints())
      p->clearContactFaces();
    for (const auto &f : n->getFaces())
      f->clearContactPoints();
  }
  {
#if defined(debug_addContact)
    std::cout << watch() << std::endl;
    std::cout << "Buckets<networkFace>に面を追加" << std::endl;
#endif
    //! 1) netの面情報
    Buckets<networkFace> B({{-200., 200.}, {-200., 200.}, {-200., 200.}}, radius);
    for (const auto &n : nets)
      for (const auto &f : n->getFaces())
        for (const auto &X : fullparticlize(f, radius))
          B.add(X, f);

#if defined(debug_addContact)
    std::cout << watch() << std::endl;
    std::cout << "fullparticlize face add done" << std::endl;
#endif
    for (const auto &n : nets)
    {
      auto ps = n->getPoints();
#ifdef _OPENMP
#pragma omp parallel
#endif
      for (auto it = ps.begin(); it != ps.end(); ++it)
      {
        //@ ------------------------------------------------------ */
        //@                    点に接触面を保存                      */
        //@ ------------------------------------------------------ */
#ifdef _OPENMP
#pragma omp single nowait
#endif
        (*it)->addContactFaces(/*Buckets<networkFace>*/ B, (*it)->radius);
      }
    }
#if defined(debug_addContact)
    std::cout << watch() << std::endl;
    std::cout << "addContactFaces done" << std::endl;
#endif
    //@ ------------------------------------------------------ */
    //@          接触された面には，接触している点を保存              */
    //@ ------------------------------------------------------ */
    for (const auto &net : nets)
      for (const auto p : net->getPoints())
        for (const auto &f : p->getContactFaces())
          f->addContactPoints(p);
  }
  /* ------------------------------------------------------ */
  {
#if defined(debug_addContact)
    std::cout << watch() << std::endl;
    std::cout << "Buckets<networkPoint>に点を追加" << std::endl;
#endif
    Buckets<networkPoint> B({{-200., 200.}, {-200., 200.}, {-200., 200.}}, radius);
    //! 1) netの点情報
    for (const auto &net : nets)
      for (const auto &p : net->getPoints())
        B.add(p->getXtuple(), p);

#if defined(debug_addContact)
    std::cout << watch() << std::endl;
    std::cout << "addContactPoints" << std::endl;
#endif
    for (const auto &n : nets)
    {
      auto ps = n->getPoints();
#ifdef _OPENMP
#pragma omp parallel
#endif
      for (auto it = ps.begin(); it != ps.end(); ++it)
      {
#ifdef _OPENMP
#pragma omp single nowait
#endif
        {
          //% ------------------------------------------------------ */
          //%                    点に接触点を保存                       */
          //% ------------------------------------------------------ */
          auto p = (*it);
          for (const auto &q : DeleteDuplicates(Flatten(B.getObjects(p->getXtuple(), p->radius /*depth*/))))
            if ((Norm(p->getXtuple() - q->getXtuple()) - (p->radius + q->radius)) < 0.)
              if (p != q)
                p->addContactPoints(q);
        }
      }
    }
#if defined(debug_addContact)
    std::cout << watch() << std::endl;
    std::cout << "addContactPoints done" << std::endl;
#endif
  }
  /* ------------------------------------------------------ */
};

/* ------------------------------------------------------ */
int main()
{
  TimeWatch watch;
  try
  {
    /* ----------------- networkPointの生成，配置 ----------------- */
    auto net = new Network;
    // net->name = "ball";
    for (const auto &x : Subdivide(-30., 0., 10))
      for (const auto &y : Subdivide(-30., 0., 10))
        for (const auto &z : Subdivide(95., 95. + 150., 20))
        {
          auto rand = RandomReal({-.1, .1});
          auto tmp = new networkPoint(net, net, {x + rand, y + rand, z + rand});
          tmp->radius = RandomReal({3., 6.});
          tmp->mass = (4. * M_PI / 3.) * pow(tmp->radius, 3);
        }
    /* ------------------------------------------------------ */

    auto obj = new Network("../../obj/bunny.obj", "bunny");
    obj->scale({500., 500., 500.});
    mk_vtu(home_dir + "/vtu/obj.vtu", obj->getFaces());

    for (const auto &f : obj->getFaces())
      for (const auto &X : particlize(f, 2.))
        new networkPoint(obj, obj, X);

    for (const auto &p : obj->getPoints())
    {
      auto tmp = new networkPoint(obj, obj, p->getXtuple());
      tmp->radius = .5;
      tmp->mass = (4. * M_PI / 3.) * pow(tmp->radius, 3) / 10.;
    }

    /* ------------------------------------------------------ */
    // うさぎが力を受けて動くか．
    // 2つの力の受け方．
    // 固体ー固体の相互作用
    // 流体ー固体の相互作用
    auto tank = new Network("./obj/tank2.obj", "tank");
    tank->scale({200., 200., 110.});
    mk_vtu(home_dir + "/vtu/tank.vtu", tank->getFaces());
    for (const auto &f : tank->getFaces())
      for (const auto &X : particlize(f, 4.))
      {
        auto tmp = new networkPoint(tank, tank, X);
        tmp->radius = 3.;
        tmp->mass = (4. * M_PI / 3.) * pow(tmp->radius, 3);
      }

    //* ------------------------------------------------------ */
    //*                       メインループ                       */
    //* ------------------------------------------------------ */
    Print("main loop");
    int count = 0;
    for (auto step = 0; step < 3000; step++)
    {
      //@ ------------------------------------------------------ */
      //@                各点や面にかかる力を計算する　　　　　　　　　　*/
      //@ ------------------------------------------------------ */

      std::cout << Grid({"Points size", net->getPoints().size()}, 70) << std::endl;
      Print(std::to_string(step));
      //すでに保存したContactPの点から受ける力を計算する
      /*
      |<--*-->|<--*-->|
      */
      addContact({net, obj, tank}, 3. /*radius*/);

      {
        Print("/* ------------------- 点にかかる力の初期化 ------------------------ */");
        auto points = net->getPoints();
        for (auto it = points.begin(); it != points.end(); ++it)
        {
          // 重力で初期化するようにしなければならない
          auto p = (*it);
          p->force = {0., 0., -g * p->mass, 0., 0., 0.};
        }
        Print("/* ------------------- 点同士の相互作用 ------------------------ */");
        double k = 30.;
        double dump = 10.;
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (auto it = points.begin(); it != points.end(); ++it)
        {
#ifdef _OPENMP
#pragma omp single nowait
#endif
          {
            auto p = (*it);
            // 重力で初期化するようにしなければならない
            T6d Vtothis({0., 0., 0., 0., 0., 0.});
            Tddd L;
            double del;
            for (const auto &q : p->getContactPoints(net) /*接触点の中でnetに属するものだけ抜き出すことができる*/)
              if (p != q)
                if ((Norm(p->getXtuple() - q->getXtuple()) - (p->radius + q->radius)) < 0.)
                {
                  L = p->getXtuple() - q->getXtuple();
                  Vtothis = {std::get<0>(L), std::get<1>(L), std::get<2>(L), 0., 0., 0.};
                  //読み方
                  p->force += -k * (Norm(L) - (p->radius + q->radius)) * Vtothis; //(p->r + q->r)からのズレ(Norm(L) - (p->r + q->r))に反抗する力が働く
                  p->force += -dump * Dot(p->V - q->V, Vtothis) * Vtothis;        // q->Vからのズレ(p->V - q->V)に反抗する力が働く
                }
          }
        }
      }
      Print("/* ------------------------------------------------------ */");
      {
        // addContact({net, obj, tank}, 3.);
        for (const auto &p : obj->getPoints())
          p->force = {0, 0, 0, 0, 0, 0};
        for (const auto &f : obj->getFaces())
          f->force = {0, 0, 0, 0, 0, 0};
        for (const auto &f : tank->getFaces())
          f->force = {0, 0, 0, 0, 0, 0};

        Print("/* ------------------- 点と接触面との作用反作用を計算 ------------------------ */");
        calcForce_interaction_PointsAndContactFaces(net->getPoints());

        Print("/* ------------------- 物体は面にかかる力の情報から点にかかる力を計算する ------------------------ */");
        //果たしていい計算方法なのだろうか？？
        //% ------------------------------------------------------ */
        //%      可視化用に面から点にかかる作用力を計算してみる            */
        //% ------------------------------------------------------ */
        for (const auto &n : V_Netp{obj, tank})
          for (const auto &p : n->getPoints())
          {
            auto faces = p->getFaces();
            for (const auto &f : faces)
              p->force += f->force / faces.size();
          }
      }

      //@ ------------------------------------------------------ */
      //@      物体全体にかかる力を計算し，物体全体の加速度を計算する　　　*/
      //@ ------------------------------------------------------ */

      Print("/* ----------------------- 物体全体にかかる力の計算　------------------------------- */");
      //点や面にかかる力を計算した後に，物体全体にかかる力を計算する．
      //その次に，物体全体の加速度を計算する．
      obj->center_of_mass = {0, 20, 50}; //適当に決めている
      obj->inertia = {1, 1, 1, 1, 1, 1}; //回転させる必要がある
      obj->inertia *= 1E+7;
      // obj->calcPhysicalProperties();
      obj->sumForceOfPoints();
      // obj->integrateForceOnFace();
      obj->calcAccelFromForce();
      Print(Grid({"force", "inertia", "center_of_mass", "acceleration"}, 70));
      Print(Grid({obj->force, obj->inertia, obj->center_of_mass, obj->acceleration}, 70));

      //@ ------------------------------------------------------ */
      //@                      最後．物体の回転                     */
      //@ ------------------------------------------------------ */

      Print("最後．物体の回転");
      auto [ax, ay, az, awx, awy, awz] = obj->acceleration;
      T6d A = {0, 0, 0, awx, awy, awz};
      obj->velocity += A * dt;
      auto [vx, vy, vz, vwx, vwy, vwz] = obj->velocity;
      Tddd w = {vwx, vwy, vwz};
      Quaternion Q;
      std::cout << "Q = " << Q() << std::endl;
      Q += Q.d_dt(w) * dt;
      // 時計回りが正
      // Q += Q.d_dt(Tddd{M_PI / 4, 0, 0}); //テスト
      std::cout << "Q.d_dt(w) = " << (Q.d_dt(w))() << std::endl;
      std::cout << "Q.d_dt(w) * dt = " << (Q.d_dt(w) * dt)() << std::endl;
      std::cout << "w = " << w << std::endl;
      std::cout << "Q = " << Q() << std::endl;
      obj->rotate(Q, obj->center_of_mass);
      /* ------------------------------------------------------ */
      /*                           出力                          */
      /* ------------------------------------------------------ */
      if (step % 3 == 0)
      {
        Print("出力");
        output_points(net, net->getName() + std::to_string(count));
        output_faces(obj, obj->getName() + std::to_string(count));
        // output_points(obj, obj->getName() + std::to_string(count), obj);
        // output_points(tank, tank->getName() + std::to_string(count));
        count++;
      }
      /* ------------------------------------------------------ */
      // update(net->getPoints(), {obj->getFaces(), tank->getFaces()});
      update(net->getPoints());
      // updatePosition_considering_refrection(net->getPoints(), tank->getFaces());
      // update(obj);
    }
  }
  catch (error_message &e)
  {
    e.print();
    abort();
  };
  std::cout << Green << "Elapsed time: " << Red << watch.get() << reset << " s\n";
  return 0;
};
