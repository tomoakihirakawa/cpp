int time_step;
#define simulation
#include "GNUPLOT.hpp"
#define DEM
#include "Network.hpp"
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
            }                                                                              // LaplacianSmoothingIfFlat(line->getPoints() /*内部でシャッフルする*/);
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
                mk_vtu("./vtu/remesh" + std::to_string(count++) + ".vtu", obj0->getFaces());
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
        //     mk_vtu("./vtu/remesh" + std::to_string(count++) + ".vtu", obj0->getFaces());
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

bool linkIfCloser_findContancP(const netPp p, const netPp q, const V_i &upperlower = {5, 8})
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

                //delete if too many points
                {
                    longestL = *(p->getLines().rbegin());
                    if (p->getLines().size() > upperlower[1] && (*longestL)(p)->getLines().size() > upperlower[1])
                        delete longestL;
                }

                //delete if too many points
                {
                    longestL = *(q->getLines().rbegin());
                    if (q->getLines().size() > upperlower[1] && (*longestL)(q)->getLines().size() > upperlower[1])
                        delete longestL;
                }

                if (replaced_P->getLines().size() < upperlower[0])
                    for (const auto &n : neighbors)
                        linkIfCloser_findContancP(replaced_P, n);

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
//                   linkIfCloser_findContancP(replaced_P, n);
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
V_d force(const netPp p)
{
    // double radius = p->radius;
    double k = 30.;
    double dump = 1.;
    V_d ret(3, 0.);
    V_d Vtothis;
    double del;
    V_d X = p->getX();
    for (const auto &q : p->contactP)
    {
        Vtothis = X - q->getX();
        if ((del = Norm3d(Vtothis) - (p->radius + q->radius)) < 0.)
        {
            ret += (k * pow(-del, 1.5) - dump * Dot(p->V - q->V /*dir to this*/, Vtothis)) * Vtothis;
        }
    }
    return ret;
};

void calculateForceFromContancP(const V_netPp &points)
{
    try
    {
        double k = 30.;
        double dump = 1.;
#ifdef _OPENMP
        Print("並列化@calculateForceFromContancP");
#pragma omp parallel for
#endif
        for (auto i = 0; i < points.size(); i++)
        { // double radius = p->radius;
            auto p = points[i];
            V_d Vtothis(3, 0.);
            double del;
            V_d X = p->getX();
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                V_netPp ps = p->contactP;
                for (const auto &q : ps)
                {
                    network::erase(p->contactP, q);
                    network::erase(q->contactP, p);
                    Vtothis = X - q->getX();
                    del = Norm3d(Vtothis) - (p->radius + q->radius);
                    if (del < 0.)
                    {
                        V_d tmp_force = (k * pow(-del, 1.5) - dump * Dot(p->V - q->V /*dir to this*/, Vtothis)) * Vtothis;
                        points[i]->F += tmp_force;
                        q->F -= tmp_force;
                    }
                }
            }
        }
    }
    catch (error_message &e)
    {
        e.print();
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
    }
};

void calculateActionAndReaction(const V_netPp &points, const V_netFp &faces)
{
    try
    {
        double k = 30. * 2.;
        double dump = 1.;
#ifdef _OPENMP
        Print("並列化@calculateActionAndReaction");
#pragma omp parallel for
#endif
        for (auto i = 0; i < points.size(); i++)
        {
            auto p = points[i];
            V_d Vtothis;
            V_d X = p->getX();
            double del;
            for (const auto &f : faces)
            {
                Vtothis = Dot(X - f->getX(), f->getNormal()) * f->getNormal();
                // normalがX-Xfと逆の場合，normalをかけると符号がかわり，Xを向く
                del = p->radius - Norm3d(Vtothis);
                //まだヒットしたかわからない
                if (del > 0.)
                {
                    V_d normal = f->getNormal();
                    VV_d AB = {X + p->radius * normal, X - p->radius * normal};
                    if (isIntersectingSurface(f->getLocations(), AB))
                    {
                        V_d tmp_force = Vtothis * (k * pow(del, 2.) - dump * Dot(p->V /*- q->V*/ /*dir to this*/, Vtothis));
#ifdef _OPENMP
#pragma omp critical
#endif
                        {
                            points[i]->F += tmp_force;
                            f->F -= tmp_force;
                        }
                    }
                }
            }
        }
    }
    catch (error_message &e)
    {
        e.print();
        throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
    }
};

V_d force(const netPp p, const netFp f)
{
    // double radius = p->radius;
    double k = 30.;
    double dump = 1.;
    V_d ret(3, 0.);
    V_d Vtothis;
    double del;
    V_d X = p->getX();
    for (const auto &q : p->getNeighbors())
    {
        Vtothis = Dot(X - f->getX(), f->getNormal()) * f->getNormal();
        // normalがX-Xfと逆の場合，normalをかけると符号がかわり，Xを向く
        del = p->radius - Norm3d(Vtothis);
        //まだヒットしたかわからない
        if (del > 0.)
        {
            auto normal = f->getNormal();
            auto AB = {X + p->radius * normal, X - p->radius * normal};
            if (isIntersectingSurface(f->getLocations(), AB))
            {
                ret += (k * pow(del, 1.5) - dump * Dot(p->V /*- q->V*/ /*dir to this*/, Vtothis)) * Vtothis;
            }
        }
    }
    return ret;
};

V_d force(const netFp f, const netPp p)
{
    return -force(p, f);
};

///////////////////////////////////////////////////

void updatePosition_considering_refrection(const V_netPp &points, const V_netFp &faces)
{
    double dt = 0.002;
    Print("baloon facesとの反射を考慮したpointsの時間発展");
    VVV_d p0p1p2s;
    for (const auto &f : faces)
        p0p1p2s.emplace_back(f->getLocations());

// 自分の周り以外との接触回避
#ifdef _OPENMP
    Print("並列化");
#pragma omp parallel for
#endif
    for (auto i = 0; i < points.size(); i++)
    {
        auto p = points[i];
        auto X = p->getX();
        bool ishit = false;
        V_d dx = p->V * dt;
        VV_d v_Vnew = {p->V + (p->F / p->mass) * dt};
        VV_d v_dx = {};
        // int hitIndex = -1;
        V_i hitIndcies = {};
        do
        {
            // std::cout << "p = " << p << ", hitIndex = " << hitIndex << std::endl;
            geometry::intersectionTriangleLine LT(p0p1p2s, X, X + dx, hitIndcies);
            // std::cout << "p = " << p << ", LT.isIntersect = " << LT.isIntersect << std::endl;
            ishit = LT.isIntersect;
            if (ishit)
            {
                // std::cin.ignore();
                // hitIndex = LT.indexOfTriangle;
                hitIndcies.emplace_back(LT.indexOfTriangle);
                // std::cout << "hitIndcies = " << hitIndcies << std::endl;
                X = LT.X; //on wall
                // std::cout << "dx = " << dx << std::endl;
                dx = LT.vecX2B_; //残りdx
                // std::cout << "dx = " << dx << std::endl;
                v_dx.emplace_back(LT.vecA2X);
                auto Vnew = LT.reflectIfPossible(*v_Vnew.rbegin());
                v_Vnew.emplace_back(Vnew);
            }
            else
            {
                v_dx.push_back(dx);
            }
        } while (ishit);

#ifdef _OPENMP
#pragma omp critical
#endif
        {
            p->setX(p->getX() + Sum(v_dx));
            p->V = *v_Vnew.rbegin();
        }
    }
};

//////////////////////////////////////////////////////////////////
double g = 9.81;
/////////////////////////////////////////////////////////
int main()
{
    try
    {

        /* ----------------- networkPointの生成，配置 ----------------- */

        auto net = new Network;
        for (const auto &x : Subdivide(-30., 30., 5))
            for (const auto &y : Subdivide(-30., 30., 5))
                for (const auto &z : Subdivide(0., 20., 5))
                {
                    auto rand = RandomReal({-.1, .1});
                    auto tmp = new networkPoint(net, net, {x + rand, y + rand, z + rand});
                    tmp->radius = RandomReal({3., 6.});
                    tmp->mass = (4. * M_PI / 3.) * pow(tmp->radius, 3) / 14.;
                }

        auto points = net->getPoints();

        VV_netPp vvp(points.size());
#ifdef _OPENMP
        Print("並列化");
#pragma omp parallel for
#endif
        for (auto i = 0; i < points.size(); i++)
            vvp[i] = neighborsSort(points[i], points);

        Print("これは並列化できない");
        for (auto i = 0; i < points.size() /*N*/; i++)
            for (auto j = 0; j < vvp[i].size() /*5*5=25*/; j++)
                if (!linkIfCloser_findContancP(points[i], vvp[i][j], {10, 20}))
                    break;

        ///////////////////
        NetworkObj tank("./obj/tank2.obj", "Neumann:tank2");
        tank.scale({100., 100., 110.});
        mk_vtu("./vtu/tank.vtu", tank.getFaces());

        ////////////////////
        NetworkObj baloon("./obj/cube.obj", "Neumann:cube");
        // obj.scale({25., 25., 25.});
        baloon.scale({70., 70., 100.});
        baloon.translate({0., 0., 10.});
        remesh({&baloon}, 9.);
        mk_vtu("./vtu/baloon.vtu", baloon.getFaces());
        // baloon.rotate(M_PI / 180 * 20, {0., 1., 0.});

        auto baloon_points = baloon.getPoints();
        for (auto i = 0; i < baloon_points.size(); i++)
        {
            baloon_points[i]->F = {0., 0., 0.};
            baloon_points[i]->radius = 1.5;
            baloon_points[i]->mass = pow(8., 3) * (4. * M_PI / 3.) * pow(baloon_points[i]->radius, 3) / 14;
        }
        ///////////////////////////////////////////////////////////////
        using VV_SorIorMap = std::vector<std::vector<std::variant<std::string, int, map_P_Vd>>>;
        map_P_Vd P_init_x, P_init_y, P_init_z;
        for (const auto &p : points)
        {
            P_init_x[p] = {p->getX()[0]};
            P_init_y[p] = {p->getX()[1]};
            P_init_z[p] = {p->getX()[2]};
        }
        // VV_SorIorMap data = {{"init_x", 1, P_x}, {"init_y", 1, P_y}, {"init_z", 1, P_z}};
        ///////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////
        Print("main loop");
        int count = 0;
        for (auto step = 0; step < 200000; step++)
        {

            if (step <= 2000 && step % 1000 == 0)
            {
                for (const auto &x : Subdivide(-30., 30., 4))
                    for (const auto &y : Subdivide(-30., 30., 4))
                        for (const auto &z : Subdivide(20., 45., 2))
                        {
                            auto tmp = new networkPoint(points[0]->getNetwork(), points[0]->getNetwork(), {x, y, z});
                            tmp->V = {0, 0, -g};
                            tmp->radius = RandomReal({4., 5.});
                            tmp->mass = (4. * M_PI / 3.) * pow(tmp->radius, 3) / 14.;
                        }
            }
            else if (step <= 3000 && step % 500 == 0)
            {
                for (const auto &x : Subdivide(-30., 30., 4))
                    for (const auto &y : Subdivide(-30., 30., 4))
                        for (const auto &z : Subdivide(20., 45., 2))
                        {
                            auto tmp = new networkPoint(points[0]->getNetwork(), points[0]->getNetwork(), {x, y, z});
                            tmp->V = {0, 0, -g};
                            tmp->radius = RandomReal({4., 5.});
                            tmp->mass = (4. * M_PI / 3.) * pow(tmp->radius, 3) / 14;
                        }
            }
            // else if (step <= 15000 && step % 500 == 0)
            // {
            //   for (const auto &x : Subdivide(-5., 5., 4))
            //     for (const auto &y : Subdivide(-5., 5., 4))
            //       for (const auto &z : Subdivide(50., 55., 2))
            //       {
            //         auto tmp = new networkPoint(points[0]->getNetwork(), points[0]->getNetwork(), {x, y, z});
            //         tmp->V = {0, 0, -g};
            //         tmp->radius = .75;
            //       }
            // }

            // else if (step >= 5001)
            // {
            //   auto ps = net->getPoints();
            //   double mean_outside_dir = 0;
            //   for (const auto &p : ps)
            //     mean_outside_dir += Norm3d(outside_direction(p));
            //   mean_outside_dir /= (double)ps.size();
            //   //
            //   for (const auto &p : ps)
            //   {
            //     double magnitude = Norm3d(outside_direction(p));
            //     if (mean_outside_dir * 0.7 < magnitude && magnitude < mean_outside_dir * 1.3)
            //     {
            //       auto r = p->radius;
            //       auto X = p->getX();
            //       ////////
            //       auto tmp = new networkPoint(points[0]->getNetwork(), points[0]->getNetwork(), X + V_d{0., 0., -r / 2.});
            //       tmp->radius = r / 2.;
            //       p->setX(X + V_d{0., 0., r / 2.});
            //       p->radius = r / 2.;
            //     }
            //   }
            // }

            // double radius = 1.5;

            // auto force = [&radius](const netPp &p, const netPp &q) {
            //   V_d ret(3, 0.);
            //   V_d Vtothis = p->getX() - q->getX();
            //   double del = Norm(Vtothis) - 2. * radius;
            //   if (del < 0.)
            //   {
            //     double k = 30.;
            //     double dump = 1.;
            //     ret = (k * (-del) - dump * Dot(p->V - q->V /*dir to this*/, Vtothis)) * Vtothis;
            //   }
            //   return ret;
            // };

            std::cout << "net->getLines().size()  = " << net->getLines().size() << std::endl;

            Print(std::to_string(step));

            points = net->getPoints();
            auto baloon_faces = baloon.getFaces();

            ///////////////////////////////////////////
            for (const auto &p : points)
                p->contactP.clear();
            // #ifdef _OPENMP
            //       Print("並列化");
            // #pragma omp parallel for
            // #endif

            ////////// random /////////
            // auto POINTS = points;
            // std::shuffle(std::begin(POINTS), std::end(POINTS), std::default_random_engine());

            // #ifdef _OPENMP
            //       Print("並列化");
            // #pragma omp parallel for
            // #endif

            Print("shift linkIfCloser_findContancP");
            int shift = 50;
            for (auto i = 0; i < shift; i++)
                for (const auto &p : points)
                    linkIfCloser_findContancP(points[(i + shift * step) % points.size()], p, {10, 20});

            //////////////////
            auto start = std::chrono::high_resolution_clock::now();
            vvp.clear();
            vvp.resize(points.size());
#ifdef _OPENMP
            Print("並列化");
#pragma omp parallel for
#endif
            for (auto i = 0; i < points.size(); i++)
                vvp[i] = neighborsSort(points[i], 2);

            Print("これは並列化できない");
            for (auto i = 0; i < points.size() /*N*/; i++)
                for (const auto &vp : vvp[i])
                    if (!linkIfCloser_findContancP(points[i], vp, {10, 20}))
                        break;

            auto finish = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = finish - start;
            std::cout << Green << "Elapsed time: " << Red << elapsed.count() << reset << " s\n";

            //O(N) complexity

            // test(points);

            //////////////// points force ///////////////////////
            Print("点へ力の計算");
            for (auto i = 0; i < points.size(); i++)
                points[i]->F = {0., 0., -g};

            Print("面への力の計算");
            for (auto i = 0; i < baloon_faces.size(); i++)
                baloon_faces[i]->F = {0., 0., 0.};
            for (auto i = 0; i < baloon_points.size(); i++)
                baloon_points[i]->F = {0., 0., 0.};
            /////////////////////////////////////////////////////
            //すでに保存したContactPの点から受ける力を計算する
            calculateForceFromContancP(points);

            //与えられた全ての面からの力を計算する
            //反作用も計算する
            calculateActionAndReaction(points, baloon_faces);
            //////////////////  line force: tension ////////////////
            for (const auto &l : baloon.getLines())
            {
                if (l->length() < 2.)
                    l->tension = 0.;
                else
                    l->tension = l->length();
            }
            ///////////////// sphere point force ////////////////////
#ifdef _OPENMP
            Print("並列化");
#pragma omp parallel for
#endif
            for (auto i = 0; i < baloon_points.size(); i++)
            {
                //面を作る点にかかる力は，線の張力と面にかかる力
                V_d tensionDir;
                for (const auto &l : baloon_points[i]->getLines())
                {
                    tensionDir = ((*l)(baloon_points[i]))->getX() - baloon_points[i]->getX();
                    baloon_points[i]->F += (l->tension) * tensionDir;
                }
                for (const auto &f : baloon_points[i]->getFaces())
                    baloon_points[i]->F += f->F / 3.;
            }

#ifdef _OPENMP
            Print("並列化");
#pragma omp parallel for
#endif
            for (auto i = 0; i < baloon_points.size(); i++)
                baloon_points[i]->F += force(baloon_points[i]);

            //////////////////////////////////////////////////////
            //////////////////////// 出力 /////////////////////////
            //////////////////////////////////////////////////////
            if (step % 10 == 0)
            {
                Print("出力");
                {
                    map_P_Vd P_V, P_F, P_outside, P_r;
                    for (const auto &p : points)
                    {
                        P_V[p] = p->V;
                        P_F[p] = force(p);
                        P_outside[p] = {Norm(outside_direction(p))};
                        P_r[p] = {p->radius};
                    }

                    VV_SorIorMap data = {{"init_x", 1, P_init_x}, {"init_y", 1, P_init_y}, {"init_z", 1, P_init_z}, {"v", 3, P_V}, {"f", 3, P_F}, {"outside_dir", 1, P_outside}, {"r", 1, P_r}};
                    mk_vtu("./vtu/points" + std::to_string(count) + ".vtu", {net->getPoints()}, data);
                    mk_vtu("./vtu/lines" + std::to_string(count) + ".vtu", {net->getLines()});
                }
                ////////////// baloon ////////////////////
                {
                    map_P_Vd P_V, P_r, P_F;
                    for (const auto &p : baloon_points)
                    {
                        P_V[p] = p->V;
                        P_r[p] = {p->radius};
                        P_F[p] = p->F;
                    }
                    VV_SorIorMap data = {{"force", 3, P_F}, {"v", 3, P_V}, {"f", 3, P_F}, {"r", 1, P_r}};
                    mk_vtu("./vtu/baloon_faces" + std::to_string(count) + ".vtu", baloon_faces, data);
                }
                count++;
            }
            ////////////////////////////////////////////////////////////////////
            ////////////////////// points vs baloon faces //////////////////////
            ////////////////////////////////////////////////////////////////////
            updatePosition_considering_refrection(points, baloon_faces);
            ////////////////////////////////////////////////////////////////////
            //////////////// baloon_points vs tank face ////////////////////////
            ////////////////////////////////////////////////////////////////////
            updatePosition_considering_refrection(baloon_points, tank.getFaces());
        }
    }
    catch (error_message &e)
    {
        e.print();
        abort();
    };
    return 0;
};
