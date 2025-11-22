#include <utility>
#include "Network.hpp"
#include "vtkWriter.hpp"
template <typename T>
std::vector<T4Tddd> toCubeFaces(const T &accum) {
   std::vector<T4Tddd> cube(6 * accum.size());
   int i = 0;
   for (const auto &atree : accum) {
      auto [X0, X1] = std::get<0>(atree->bounds);
      auto [Y0, Y1] = std::get<1>(atree->bounds);
      auto [Z0, Z1] = std::get<2>(atree->bounds);
      cube[i++] = T4Tddd{{X0, Y0, Z0}, {X1, Y0, Z0}, {X1, Y1, Z0}, {X0, Y1, Z0}};
      cube[i++] = T4Tddd{{X0, Y0, Z1}, {X1, Y0, Z1}, {X1, Y1, Z1}, {X0, Y1, Z1}};
      cube[i++] = T4Tddd{{X0, Y0, Z0}, {X1, Y0, Z0}, {X1, Y0, Z1}, {X0, Y0, Z1}};
      cube[i++] = T4Tddd{{X0, Y1, Z0}, {X1, Y1, Z0}, {X1, Y1, Z1}, {X0, Y1, Z1}};
      cube[i++] = T4Tddd{{X0, Y0, Z0}, {X0, Y0, Z1}, {X0, Y1, Z1}, {X0, Y1, Z0}};
      cube[i++] = T4Tddd{{X1, Y0, Z0}, {X1, Y0, Z1}, {X1, Y1, Z1}, {X1, Y1, Z0}};
   }
   return cube;
};

Tddd getX(const networkPoint *X) { return X->getXtuple(); };

#define check_SPH3

#ifdef check_SPH3
int main(int argc, char **argv) {
   Timer timer;
   /* ----------------------- 引数の読み込み ---------------------- */
   std::string name{argv[1]};  //"/Users/tomoaki/Dropbox/markdown/cpp/obj/cow.obj";
   int depth = std::atoi(argv[2]);
   /* ---------------------- ポリゴンの読み込み --------------------- */
   auto net = new Network(name, "object");
   octree tree(net->getBounds().bounds, depth, extVertices(net->getFaces()));  //八分木構造を生成．ポリゴン外部のキューブは自動で内部で削除せれる

   {
      // vtkXMLPolyDataWriter<key here>;
      vtkXMLPolyDataWriter</*specidy key*/ networkPoint *> vtp;
      for (const auto &f : net->getFaces()) {
         auto [p0, p1, p2] = f->getPointsTuple();
         vtp.add({p0, p1, p2});
         vtp.addPolygon({p0, p1, p2});
      }
      std::ofstream ofs("./output/bunny.vtp");
      vtp.write(ofs);
   }

   /* ------------------------- ボクセルの出力 ------------------------- */
   mk_vtu("./output/cubeAll.vtu", toCubeFaces(tree.getAllDeepest()));
   //$ ------------------------------------------------------ */
   //$           setting SPHに使うような点群を生成                */
   //$ ------------------------------------------------------ */
   auto particles = new Network;
   double particle_spacing = 0.04;
   /* ----------------- particlesworkPointの生成，配置 ----------------- */
   auto [xb0, xb1] = Tdd{-0.05, 0.05};
   auto [yb0, yb1] = Tdd{0.05, 0.15};
   auto [zb0, zb1] = Tdd{0.06, 0.06 + 0.1};
   auto X = Subdivide(xb0, xb1, (int)std::round((xb1 - xb0) / particle_spacing));  //粒子のX座標
   auto Y = Subdivide(yb0, yb1, (int)std::round((yb1 - yb0) / particle_spacing));  //粒子のY座標
   auto Z = Subdivide(zb0, zb1, (int)std::round((zb1 - zb0) / particle_spacing));  //粒子のZ座標
   auto V_vertices = extVertices(particles->getFaces());
   for (const auto &x : X)
      for (const auto &y : Y)
         for (const auto &z : Z)
            new networkPoint(particles, particles, {x, y, z});
   //$ ------------------------------------------------------ */
   auto trees = tree.getAllDeepest();
   /* ------------------------------------------------------ */
   /*              八分木構造と球体の干渉のチェックと出力           */
   /* ------------------------------------------------------ */
   //! 干渉のチェックを行い，球に入るキューブを抜き出す
   int i = 0;
   for (const auto &p : particles->getPoints()) {
      std::cout << "timer : " << timer() << std::endl;
      //
      {
         geometry::Sphere sp(p->getXtuple(), 0.05);  //球体オブジェクト
         vtkXMLPolyDataWriter<std::shared_ptr<Tddd>> vtp;
         for (const auto &[X0, X1, X2, X3] : toCubeFaces(tree.getIntersect(sp))) {
            std::shared_ptr<Tddd> x0(new Tddd(X0));
            std::shared_ptr<Tddd> x1(new Tddd(X1));
            std::shared_ptr<Tddd> x2(new Tddd(X2));
            std::shared_ptr<Tddd> x3(new Tddd(X3));
            vtp.add({x0, x1, x2, x3});
            vtp.addPolygon({x0, x1, x2, x3});
         }
         std::ofstream ofs("./output/cube" + std::to_string(i) + ".vtp");
         vtp.write(ofs);
         ofs.close();
      }
      //
      {
         vtkXMLPolyDataWriter<std::shared_ptr<Tddd>> vtp;
         std::shared_ptr<Tddd> x(new Tddd(p->getXtuple()));
         vtp.add(x);
         std::ofstream ofs("./output/particle" + std::to_string(i) + ".vtp");
         vtp.write(ofs);
      }
      i++;
   }
   //! ------------------------------------------------------ */
   //!                         vtpの出力                       */
   //! ------------------------------------------------------ */
   std::cout << "timer : " << timer() << std::endl;
   {
      vtkXMLPolyDataWriter<std::shared_ptr<Tddd>> vtp;
      double v = 0;
      Tddd u = {0., 0., 0.};
      for (const auto &t : trees) {
         auto [X0, X1, X2, X3, X4, X5, X6, X7] = t->getVertices();
         std::shared_ptr<Tddd> x0(new Tddd(X0));
         std::shared_ptr<Tddd> x1(new Tddd(X1));
         std::shared_ptr<Tddd> x2(new Tddd(X2));
         std::shared_ptr<Tddd> x3(new Tddd(X3));
         std::shared_ptr<Tddd> x4(new Tddd(X4));
         std::shared_ptr<Tddd> x5(new Tddd(X5));
         std::shared_ptr<Tddd> x6(new Tddd(X6));
         std::shared_ptr<Tddd> x7(new Tddd(X7));
         vtp.add({x0, x1, x2, x3, x4, x5, x6, x7});
         vtp.addPolygon({{x0, x1, x2, x3}, {x4, x5, x6, x7}, {x0, x1, x5, x4}, {x2, x3, x7, x6}, {x0, x4, x7, x3}, {x1, x2, x6, x5}});
         vtp.addPointData("scaler0", {{x0, v}, {x1, v}, {x2, v}, {x3, v}, {x4, v}, {x5, v}, {x6, v}, {x7, v}});
         v += 1.;
         vtp.addPointData("scaler1", {{x0, v}, {x1, v}, {x2, v}, {x3, v}, {x4, v}, {x5, v}, {x6, v}, {x7, v}});
         v += 1.;
         vtp.addPointData("vector", {{x0, X0 - u}, {x1, X1 - u}, {x2, X2 - u}, {x3, X3 - u}, {x4, X4 - u}, {x5, X5 - u}, {x6, X6 - u}, {x7, X7 - u}});
      }
      std::cout << vtp.verticies.size() << std::endl;
      std::cout << "outputting " << std::endl;
      std::ofstream ofs("./output/output2.vtp");
      vtp.write(ofs);
   }
};
#elif defined(check_SPH2)
int main(int argc, char **argv) {
   Timer timer;
   /* ----------------------- 引数の読み込み ---------------------- */
   std::string name{argv[1]};  //"/Users/tomoaki/Dropbox/markdown/cpp/obj/cow.obj";
   int depth = std::atoi(argv[2]);
   /* ---------------------- ポリゴンの読み込み --------------------- */
   auto net = new Network(name, "object");
   octree tree(net->getBounds().bounds, depth, extVertices(net->getFaces()));  //八分木構造を生成．ポリゴン外部のキューブは自動で内部で削除せれる
   /* ------------------------- ボクセルの出力 ------------------------- */
   std::string output_directory = "./";
   mk_vtu(output_directory + "/cubeAll.vtu", toCubeFaces(tree.getAllDeepest()));
   //$ ------------------------------------------------------ */
   //$           setting SPHに使うような点群を生成                */
   //$ ------------------------------------------------------ */
   auto particles = new Network;
   double particle_spacing = 0.04;
   /* ----------------- particlesworkPointの生成，配置 ----------------- */
   auto [xb0, xb1] = Tdd{-0.05, 0.05};
   auto [yb0, yb1] = Tdd{0.05, 0.15};
   auto [zb0, zb1] = Tdd{0.06, 0.06 + 0.1};
   auto X = Subdivide(xb0, xb1, (int)std::round((xb1 - xb0) / particle_spacing));  //粒子のX座標
   auto Y = Subdivide(yb0, yb1, (int)std::round((yb1 - yb0) / particle_spacing));  //粒子のY座標
   auto Z = Subdivide(zb0, zb1, (int)std::round((zb1 - zb0) / particle_spacing));  //粒子のZ座標
   auto V_vertices = extVertices(particles->getFaces());
   for (const auto &x : X)
      for (const auto &y : Y)
         for (const auto &z : Z)
            new networkPoint(particles, particles, {x, y, z});
   //$ ------------------------------------------------------ */
   /* ----------------- 八分木構造と球体の干渉のチェックと出力 ---------------- */
   //! 干渉のチェックを行い，球に入るキューブを抜き出す
   int i = 0;
   for (const auto &p : particles->getPoints()) {
      std::cout << "timer : " << timer() << std::endl;
      geometry::Sphere sp(p->getXtuple(), 0.05);  //球体オブジェクト
      mk_vtu(output_directory + "/cube" + std::to_string(i) + ".vtu", toCubeFaces(tree.getIntersect(sp)));
      mk_vtu(output_directory + "/particle" + std::to_string(i) + ".vtu", {p->getXtuple()});
      i++;
   }
   /* ------------------------------------------------------ */
};
#elif defined(check_SPH)
int main(int argc, char **argv) {
   Timer timer;
   /* ----------------------- 引数の読み込み ---------------------- */
   std::string name{argv[1]};  //"/Users/tomoaki/Dropbox/markdown/cpp/obj/cow.obj";
   int depth = std::atoi(argv[2]);
   /* ---------------------- ポリゴンの読み込み --------------------- */
   auto net = new Network(name, "object");
   octree tree(net->getBounds().bounds, depth, extVertices(net->getFaces()));  //八分木構造を生成．ポリゴン外部のキューブは自動で内部で削除せれる
   /* ------------------------- ボクセルの出力 ------------------------- */
   std::string output_directory = "./";
   mk_vtu(output_directory + "/cubeAll.vtu", toCubeFaces(tree.getAllDeepest()));
   //$ ------------------------------------------------------ */
   //$           setting SPHに使うような点群を生成                */
   //$ ------------------------------------------------------ */
   auto particles = new Network;
   double particle_spacing = 0.04;
   /* ----------------- particlesworkPointの生成，配置 ----------------- */
   auto [xb0, xb1] = Tdd{-0.05, 0.05};
   auto [yb0, yb1] = Tdd{0.05, 0.15};
   auto [zb0, zb1] = Tdd{0.06, 0.06 + 0.1};
   auto X = Subdivide(xb0, xb1, (int)std::round((xb1 - xb0) / particle_spacing));  //粒子のX座標
   auto Y = Subdivide(yb0, yb1, (int)std::round((yb1 - yb0) / particle_spacing));  //粒子のY座標
   auto Z = Subdivide(zb0, zb1, (int)std::round((zb1 - zb0) / particle_spacing));  //粒子のZ座標
   auto V_vertices = extVertices(particles->getFaces());
   for (const auto &x : X)
      for (const auto &y : Y)
         for (const auto &z : Z)
            new networkPoint(particles, particles, {x, y, z});
   //$ ------------------------------------------------------ */
   /* ----------------- 八分木構造と球体の干渉のチェックと出力 ---------------- */
   //! 干渉のチェックを行い，球に入るキューブを抜き出す
   int i = 0;
   for (const auto &p : particles->getPoints()) {
      std::cout << "timer : " << timer() << std::endl;
      geometry::Sphere sp(p->getXtuple(), 0.05);  //球体オブジェクト
      mk_vtu(output_directory + "/cube" + std::to_string(i) + ".vtu", toCubeFaces(tree.getIntersect(sp)));
      mk_vtu(output_directory + "/particle" + std::to_string(i) + ".vtu", {p->getXtuple()});
      i++;
   }
   /* ------------------------------------------------------ */
};
#else
int main(int argc, char **argv) {
   std::string name{argv[1]};  //"/Users/tomoaki/Dropbox/markdown/cpp/obj/cow.obj";
   int depth = std::atoi(argv[2]);
   /* ------------------------------------------------------ */
   Timer timer;
   std::string output_directory = "./";
   auto net = new Network(name, "object");
   /* ------------------------------------------------------ */
   std::cout << net->getBounds() << std::endl;
   mk_vtu(output_directory + "/object.vtu", {net->getFaces()});

   std::cout << "timer: " << timer() << std::endl;
   octree tree(net->getBounds().bounds, depth, extVertices(net->getFaces()));
   std::cout << "generate : " << timer() << std::endl;
   std::unordered_set<octree *> accum;

   tree.getAllDeepest(accum);
   std::cout << "getDescendants : " << timer() << std::endl;
   std::vector<T4Tddd> cubeIntersectWithSphere, cubesAllInside, cubesAnyInside, cubeDeepest;
   std::vector<std::vector<T4Tddd>> cubeAtDepth(10);
   std::cout << accum.size() << std::endl;
   std::vector<T4Tddd> cube;
   cube.clear();
   for (const auto &atree : accum) {
      auto [minmaxX, minmaxY, minmaxZ] = atree->bounds;
      auto [X0, X1] = minmaxX;
      auto [Y0, Y1] = minmaxY;
      auto [Z0, Z1] = minmaxZ;
      cube.push_back({{X0, Y0, Z0}, {X1, Y0, Z0}, {X1, Y1, Z0}, {X0, Y1, Z0}});
      cube.push_back({{X0, Y0, Z1}, {X1, Y0, Z1}, {X1, Y1, Z1}, {X0, Y1, Z1}});
      cube.push_back({{X0, Y0, Z0}, {X1, Y0, Z0}, {X1, Y0, Z1}, {X0, Y0, Z1}});
      cube.push_back({{X0, Y1, Z0}, {X1, Y1, Z0}, {X1, Y1, Z1}, {X0, Y1, Z1}});
      cube.push_back({{X0, Y0, Z0}, {X0, Y0, Z1}, {X0, Y1, Z1}, {X0, Y1, Z0}});
      cube.push_back({{X1, Y0, Z0}, {X1, Y0, Z1}, {X1, Y1, Z1}, {X1, Y1, Z0}});
   }
   mk_vtu(output_directory + "/cubeAll.vtu", cube);
   //

   /* ------------------------------------------------------ */
   /*                     チェックスピード                      */
   /* ------------------------------------------------------ */
   std::cout << "timer : " << timer() << std::endl;
   int loop = 10000;
   std::cout << "size : " << tree.getIntersect(geometry::Sphere({0., 0., 0.}, 0.1)).size() << std::endl;
   for (auto i = 0; i < loop; i++)
      tree.getIntersectAsVector(geometry::Sphere({0., 0., 0.}, 0.1));
   std::cout << red << "timer : " << timer() << reset << std::endl;
   /* ------------------------------------------------------ */

   for (auto i = 0; i < 10; i++) {
      // std::cout << tree.getIntersect(geometry::Sphere({0., 0., 0.}, 0.02 * (double)i)).size() << std::endl;
      cube.clear();
      for (const auto &atree : tree.getIntersect(geometry::Sphere({0., 0., 0.}, 0.02 * (double)i))) {
         auto [minmaxX, minmaxY, minmaxZ] = atree->bounds;
         auto [X0, X1] = minmaxX;
         auto [Y0, Y1] = minmaxY;
         auto [Z0, Z1] = minmaxZ;
         cube.push_back({{X0, Y0, Z0}, {X1, Y0, Z0}, {X1, Y1, Z0}, {X0, Y1, Z0}});
         cube.push_back({{X0, Y0, Z1}, {X1, Y0, Z1}, {X1, Y1, Z1}, {X0, Y1, Z1}});
         cube.push_back({{X0, Y0, Z0}, {X1, Y0, Z0}, {X1, Y0, Z1}, {X0, Y0, Z1}});
         cube.push_back({{X0, Y1, Z0}, {X1, Y1, Z0}, {X1, Y1, Z1}, {X0, Y1, Z1}});
         cube.push_back({{X0, Y0, Z0}, {X0, Y0, Z1}, {X0, Y1, Z1}, {X0, Y1, Z0}});
         cube.push_back({{X1, Y0, Z0}, {X1, Y0, Z1}, {X1, Y1, Z1}, {X1, Y1, Z0}});
      }
      std::cout << "timer : " << timer() << std::endl;
      mk_vtu(output_directory + "/cube" + std::to_string(i) + ".vtu", cube);
   }
};
#endif