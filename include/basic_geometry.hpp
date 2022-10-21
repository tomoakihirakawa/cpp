#ifndef basic_geometry_H
#define basic_geometry_H
#pragma once

#include <execution>
#include "basic_statistics.hpp"
#include "basic_vectors.hpp"

// 面と面の干渉？？？
// 球と面の干渉チェックか．
// そのために，
// タプルを作ろうとしている，

using Tdd = std::tuple<double, double>;
using Tddd = std::tuple<double, double, double>;
using T2Tddd = std::tuple<Tddd, Tddd>;
using T3Tddd = std::tuple<Tddd, Tddd, Tddd>;
using T3Tdd = std::tuple<Tdd, Tdd, Tdd>;

/* -------------------------------------------------------------------------- */
template <typename T>
Tddd circumcenter(const T &a_, const T &b_, const T &c_, const T &d_) {
   auto a = ToX(a_);
   auto b = ToX(b_);
   auto c = ToX(c_);
   auto d = ToX(d_);
   // http://rodolphe-vaillant.fr/entry/127/find-a-tetrahedron-circumcenter#:~:text=For%20all%20tetrahedra%2C%20there%20exists,circumsphere%20is%20called%20the%20circumcentre.
   double a2 = Dot(a, a);
   return Dot(Inverse(T3Tddd{b - a, c - a, d - a}),
              0.5 * Tddd{Dot(b, b) - a2, Dot(c, c) - a2, Dot(d, d) - a2});
};
Tddd circumcenter(const Tddd &a, const Tddd &b, const Tddd &c, const Tddd &d) {
   // http://rodolphe-vaillant.fr/entry/127/find-a-tetrahedron-circumcenter#:~:text=For%20all%20tetrahedra%2C%20there%20exists,circumsphere%20is%20called%20the%20circumcentre.
   double a2 = Dot(a, a);
   return Dot(Inverse(T3Tddd{b - a, c - a, d - a}),
              0.5 * Tddd{Dot(b, b) - a2, Dot(c, c) - a2, Dot(d, d) - a2});
};
/* ------------------------------------------------------ */
bool isInside(const Tddd &X, const Tddd &Xcenter, const double &r) {
   // point v.s. sphere
   return Norm(X - Xcenter) < r;
};
bool isInside(const Tddd &X, const T3Tdd &bounds) {
   // point v.s. cube
   if (std::get<0>(X) < std::get<0>(std::get<0>(bounds)) ||
       std::get<1>(std::get<0>(bounds)) < std::get<0>(X) ||
       std::get<1>(X) < std::get<0>(std::get<1>(bounds)) ||
       std::get<1>(std::get<1>(bounds)) < std::get<1>(X) ||
       std::get<2>(X) < std::get<0>(std::get<2>(bounds)) ||
       std::get<1>(std::get<2>(bounds)) < std::get<2>(X))
      return false;
   else
      return true;
};
bool isInside(const T3Tdd &bounds, const Tddd &Xcenter, const double &r) {
   // cube v.s. sphere
   // cube < sphere ?
   auto [X0, X1] = std::get<0>(bounds);
   auto [Y0, Y1] = std::get<1>(bounds);
   auto [Z0, Z1] = std::get<2>(bounds);
   return (isInside({X0, Y0, Z0}, Xcenter, r) &&
           isInside({X1, Y0, Z0}, Xcenter, r) &&
           isInside({X1, Y1, Z0}, Xcenter, r) &&
           isInside({X0, Y1, Z0}, Xcenter, r) &&
           isInside({X0, Y0, Z1}, Xcenter, r) &&
           isInside({X1, Y0, Z1}, Xcenter, r) &&
           isInside({X0, Y1, Z1}, Xcenter, r) &&
           isInside({X1, Y1, Z1}, Xcenter, r));
};
bool isInside(const Tddd &Xcenter, const double &r, const T3Tdd &bounds) {
   auto [X0, X1] = std::get<0>(bounds);
   auto [Y0, Y1] = std::get<1>(bounds);
   auto [Z0, Z1] = std::get<2>(bounds);
   return (!isInside({X0, Y0, Z0}, Xcenter, r) &&
           !isInside({X1, Y0, Z0}, Xcenter, r) &&
           !isInside({X1, Y1, Z0}, Xcenter, r) &&
           !isInside({X0, Y1, Z0}, Xcenter, r) &&
           !isInside({X0, Y0, Z1}, Xcenter, r) &&
           !isInside({X1, Y0, Z1}, Xcenter, r) &&
           !isInside({X0, Y1, Z1}, Xcenter, r) &&
           !isInside({X1, Y1, Z1}, Xcenter, r));
};
/* ------------------------------------------------------ */
/*
M. Meyer, M. Desbrun, P. Schröder, and A. H. Barr, “Discrete
Differential-Geometry Operators for Triangulated 2-Manifolds BT  - Visualization
and Mathematics III,” Vis. Math. III, pp. 35–57, 2003.
*/

struct GeometricRegions {};

Tddd ToSphericalCoodrinates(const Tddd &xyz) {
   double r = Norm(xyz);
   return {r, std::atan(std::get<2>(xyz) / r),
           std::atan2(std::get<1>(xyz), std::get<0>(xyz))};
};

/* ------------------------------------------------------ */

struct Point {
   Tddd X;
   Point(const Tddd &XIN) : X(XIN){};
};

template <typename T>
struct Edge {
   Edge(T a, T b);
};

// template <typename T>
// struct Tetrahedron
// {
// 	/*
// 	1,2,3
// 		3
// 	  / | \
// 	2---|---1
// 	 \  |  /    0,1,3
// 	  \ | /
// 		0
// 	0,3,2
// 	0,2,1
// 	*/
// 	const std::tuple<Tiii, Tiii, Tiii, Tiii> polygons = {{1, 2, 3}, {0, 1,
// 3}, {0, 2, 1}, {0, 3, 2}}; 	std::tuple<T, T, T, T> X; 	T3Tdd bounds;
// 	Tetrahedron(const std::tuple<T, T, T, T> &XIN) : X(XIN){};
// };

// template <>
// struct Tetrahedron<double>
// {
// 	/*
// 	1,2,3
// 		3
// 	  / | \
// 	2---|---1
// 	 \  |  /    0,1,3
// 	  \ | /
// 		0
// 	0,3,2
// 	0,2,1
// 	*/
// 	const std::tuple<Tiii, Tiii, Tiii, Tiii> polygons = {{1, 2, 3}, {0, 1,
// 3}, {0, 2, 1}, {0, 3, 2}}; 	T4Tddd X; 	T3Tdd bounds; 	Tetrahedron(const T4Tddd
// &XIN) : X(XIN), 									 bounds({MinMax(T4d{std::get<0>(std::get<0>(XIN)),
// std::get<0>(std::get<1>(XIN)), std::get<0>(std::get<2>(XIN)),
// std::get<0>(std::get<3>(XIN))}), 											 MinMax(T4d{std::get<1>(std::get<0>(XIN)),
// std::get<1>(std::get<1>(XIN)), std::get<1>(std::get<2>(XIN)),
// std::get<1>(std::get<3>(XIN))}), 											 MinMax(T4d{std::get<2>(std::get<0>(XIN)),
// std::get<2>(std::get<1>(XIN)), std::get<2>(std::get<2>(XIN)),
// std::get<2>(std::get<3>(XIN))})}){};
// };

T3Tdd Distance(const T3Tdd &b, const Tddd &a) {
   auto [mmX, mmY, mmZ] = b;
   T3Tdd ret = b;
   //
   Tdd to_ax = std::get<0>(a) - mmX;
   Tdd to_ay = std::get<1>(a) - mmY;
   Tdd to_az = std::get<2>(a) - mmZ;
   //
   std::get<0>(ret) = MinMax(Tdd{std::abs(std::get<0>(to_ax)),
                                 std::abs(std::get<1>(to_ax))});  //<-Tdd
   std::get<1>(ret) = MinMax(Tdd{std::abs(std::get<0>(to_ay)),
                                 std::abs(std::get<1>(to_ay))});  //<-Tdd
   std::get<2>(ret) = MinMax(Tdd{std::abs(std::get<0>(to_az)),
                                 std::abs(std::get<1>(to_az))});  //<-Tdd
   //
   if (Between(std::get<0>(a), mmX)) std::get<0>(std::get<0>(ret)) = 0;
   if (Between(std::get<1>(a), mmY)) std::get<0>(std::get<1>(ret)) = 0;
   if (Between(std::get<2>(a), mmZ)) std::get<0>(std::get<2>(ret)) = 0;
   return ret;
};

namespace geometry {
struct Line {
   T2Tddd X;
   double minX, maxX, minY, maxY, minZ, maxZ;
   T3Tdd bounds;
   Line(const T2Tddd &XIN)
       : X(XIN),
         minX(Min(Tdd{std::get<0>(std::get<0>(XIN)),
                      std::get<0>(std::get<1>(XIN))})),
         maxX(Max(Tdd{std::get<0>(std::get<0>(XIN)),
                      std::get<0>(std::get<1>(XIN))})),
         minY(Min(Tdd{std::get<1>(std::get<0>(XIN)),
                      std::get<1>(std::get<1>(XIN))})),
         maxY(Max(Tdd{std::get<1>(std::get<0>(XIN)),
                      std::get<1>(std::get<1>(XIN))})),
         minZ(Min(Tdd{std::get<2>(std::get<0>(XIN)),
                      std::get<2>(std::get<1>(XIN))})),
         maxZ(Max(Tdd{std::get<2>(std::get<0>(XIN)),
                      std::get<2>(std::get<1>(XIN))})),
         bounds({{minX, maxX}, {minY, maxY}, {minZ, maxZ}}){};
   Line(const Tddd &X0IN, const Tddd &X1IN)
       : X({X0IN, X1IN}),
         minX(Min(Tdd{std::get<0>(X0IN), std::get<0>(X1IN)})),
         maxX(Max(Tdd{std::get<0>(X0IN), std::get<0>(X1IN)})),
         minY(Min(Tdd{std::get<1>(X0IN), std::get<1>(X1IN)})),
         maxY(Max(Tdd{std::get<1>(X0IN), std::get<1>(X1IN)})),
         minZ(Min(Tdd{std::get<2>(X0IN), std::get<2>(X1IN)})),
         maxZ(Max(Tdd{std::get<2>(X0IN), std::get<2>(X1IN)})),
         bounds({{minX, maxX}, {minY, maxY}, {minZ, maxZ}}){};
};
/* ------------------------------------------------------ */
struct Triangle {
   T3Tddd X;
   Tddd center;
   T3Tdd bounds;
   Tddd normal;
   Triangle(const T3Tddd &XIN)
       : X(XIN),
         center(Mean(XIN)),
         bounds({{(Min(Tddd{std::get<0>(std::get<0>(XIN)),
                            std::get<0>(std::get<1>(XIN)),
                            std::get<0>(std::get<2>(XIN))})),
                  (Max(Tddd{std::get<0>(std::get<0>(XIN)),
                            std::get<0>(std::get<1>(XIN)),
                            std::get<0>(std::get<2>(XIN))}))},
                 {(Min(Tddd{std::get<1>(std::get<0>(XIN)),
                            std::get<1>(std::get<1>(XIN)),
                            std::get<1>(std::get<2>(XIN))})),
                  (Max(Tddd{std::get<1>(std::get<0>(XIN)),
                            std::get<1>(std::get<1>(XIN)),
                            std::get<1>(std::get<2>(XIN))}))},
                 {(Min(Tddd{std::get<2>(std::get<0>(XIN)),
                            std::get<2>(std::get<1>(XIN)),
                            std::get<2>(std::get<2>(XIN))})),
                  (Max(Tddd{std::get<2>(std::get<0>(XIN)),
                            std::get<2>(std::get<1>(XIN)),
                            std::get<2>(std::get<2>(XIN))}))}}),
         normal(Normalize(Cross(std::get<1>(X) - std::get<0>(X),
                                std::get<2>(X) - std::get<0>(X)))){};
   Triangle(const Tddd &X0IN, const Tddd &X1IN, const Tddd &X2IN)
       : X({X0IN, X1IN, X2IN}),
         center((X0IN + X1IN + X2IN) / 3.),
         bounds({{(Min(Tddd{std::get<0>(X0IN), std::get<0>(X1IN),
                            std::get<0>(X2IN)})),
                  (Max(Tddd{std::get<0>(X0IN), std::get<0>(X1IN),
                            std::get<0>(X2IN)}))},
                 {(Min(Tddd{std::get<1>(X0IN), std::get<1>(X1IN),
                            std::get<1>(X2IN)})),
                  (Max(Tddd{std::get<1>(X0IN), std::get<1>(X1IN),
                            std::get<1>(X2IN)}))},
                 {(Min(Tddd{std::get<2>(X0IN), std::get<2>(X1IN),
                            std::get<2>(X2IN)})),
                  (Max(Tddd{std::get<2>(X0IN), std::get<2>(X1IN),
                            std::get<2>(X2IN)}))}}),
         normal(Normalize(Cross(X1IN - X0IN, X2IN - X0IN))){};
};
/* ------------------------------------------------------ */
// struct Point
// {
// 	Tddd X;
// 	Point(const Tddd &XIN) : X(XIN){};
// };
/* ------------------------------------------------------ */
struct Sphere {
   Tddd X;
   double radius;
   T3Tdd bounds;
   Sphere(const Tddd &XIN, const double radiusIN = 0.)
       : X(XIN),
         radius(radiusIN),
         bounds({{(std::get<0>(XIN) - radiusIN), (std::get<0>(XIN) + radiusIN)},
                 {(std::get<1>(XIN) - radiusIN), (std::get<1>(XIN) + radiusIN)},
                 {(std::get<2>(XIN) - radiusIN),
                  (std::get<2>(XIN) + radiusIN)}}){};
};
};  // namespace geometry

/* ------------------------------------------------------ */
// structをわざわざ作るのは，T3Tddではなく，coordinateboundsとして意味を具体的にした状態で持ち回りたいから．それだけ．
struct CoordinateBounds {
   T3Tdd bounds;
   Tddd X;  // center
   Tddd &center = this->X;
   /* ------------------------------------------------------ */
   void setBounds(const std::vector<Tddd> &Vxyz) {
      this->bounds = MinMaxTranspose(Vxyz);
      this->X = Mean(Transpose(this->bounds));
   };
   void setBounds(const CoordinateBounds &bs) {
      this->bounds = bs.bounds;
      this->X = bs.X;
   };
   void setBounds(const T3Tddd &X012) {
      this->bounds = MinMaxTranspose(X012);
      this->X = Mean(X012);
   };
   void setBounds(const Tddd &x) {
      this->bounds = {{std::get<0>(x), std::get<0>(x)}, {std::get<1>(x), std::get<1>(x)}, {std::get<2>(x), std::get<2>(x)}};
      this->X = x;
   };
   const Tddd &getXtuple() const {
      return this
          ->X; /*面のsetBoundsでは，T3Tdddの平均がXとなるようになっている．バウンディングボックスの中心ではない.*/
   };
   // const Tddd &getX() const { return this->X;
   // /*面のsetBoundsでは，T3Tdddの平均がXとなるようになっている．バウンディングボックスの中心ではない.*/
   // }; V_d getX() const { return {std::get<0>(this->X), std::get<1>(this->X),
   // std::get<2>(this->X)}; };
   const T3Tdd &getBounds() const { return this->bounds; };
   /* ------------------------------------------------------ */
   CoordinateBounds() : bounds({{1E+20, -1E+20}, {1E+20, -1E+20}, {1E+20, -1E+20}}){};
   CoordinateBounds(const CoordinateBounds &bs) : bounds(bs.bounds), X(bs.X){};
   CoordinateBounds(const Tdd &minmaxX, const Tdd &minmaxY, const Tdd &minmaxZ) : bounds({minmaxX, minmaxY, minmaxZ}), X(Mean(Transpose(bounds))){};
   CoordinateBounds(const Tddd &x) : bounds({{std::get<0>(x), std::get<0>(x)}, {std::get<1>(x), std::get<1>(x)}, {std::get<2>(x), std::get<2>(x)}}), X(x){};
   CoordinateBounds(const T3Tdd &minmax) : bounds(minmax), X(Mean(Transpose(bounds))){};
   CoordinateBounds(const double minX, const double maxX, const double minY, const double maxY, const double minZ, const double maxZ) : bounds({{minX, maxX}, {minY, maxY}, {minZ, maxZ}}), X(Mean(Transpose(bounds))){};
   CoordinateBounds(const T2Tddd &x) : bounds(MinMaxTranspose(x)), X(Mean(x)){};
   CoordinateBounds(const T3Tddd &x) : bounds(MinMaxTranspose(x)), X(Mean(x)){};
   CoordinateBounds(const T4Tddd &x) : bounds(MinMaxTranspose(x)), X(Mean(x)){};
   CoordinateBounds(const Tddd &x0, const Tddd &x1, const Tddd &x2) : bounds(MinMaxTranspose(T3Tddd{x0, x1, x2})), X((x0 + x1 + x2) / 3.){};
   CoordinateBounds(const std::vector<T3Tddd> &V_X);
   CoordinateBounds(const std::vector<Tddd> &X) : bounds(MinMax(Transpose(X))), X(Mean(Transpose(bounds))){};
   CoordinateBounds(const geometry::Line &L) : bounds(MinMaxTranspose(L.X)){};
   CoordinateBounds(const geometry::Sphere &S) : bounds({{std::get<0>(S.X) - S.radius, std::get<0>(S.X) + S.radius}, {std::get<1>(S.X) - S.radius, std::get<1>(S.X) + S.radius}, {std::get<2>(S.X) - S.radius, std::get<2>(S.X) + S.radius}}){};
   CoordinateBounds(const geometry::Triangle &T) : bounds(MinMaxTranspose(T.X)){};
   /* ------------------------------------------------------ */
   Tdd Distance(const Tddd &a) const {
      auto [mmX, mmY, mmZ] = this->bounds;
      T3Tdd ret = this->bounds;
      //
      Tdd to_ax = mmX - std::get<0>(a);
      Tdd to_ay = mmY - std::get<1>(a);
      Tdd to_az = mmZ - std::get<2>(a);
      //
      std::get<0>(ret) = MinMax(Tdd{std::abs(std::get<0>(to_ax)), std::abs(std::get<1>(to_ax))});  //<-Tdd
      std::get<1>(ret) = MinMax(Tdd{std::abs(std::get<0>(to_ay)), std::abs(std::get<1>(to_ay))});  //<-Tdd
      std::get<2>(ret) = MinMax(Tdd{std::abs(std::get<0>(to_az)), std::abs(std::get<1>(to_az))});  //<-Tdd
      //
      if (Between(std::get<0>(a), mmX)) std::get<0>(std::get<0>(ret)) = 0;
      if (Between(std::get<1>(a), mmY)) std::get<0>(std::get<1>(ret)) = 0;
      if (Between(std::get<2>(a), mmZ)) std::get<0>(std::get<2>(ret)) = 0;
      auto [min, max] = Transpose(ret);
      return {Norm(min), Norm(max)};
   };
   /* ------------------------------------------------------ */
   const T3Tdd &operator()() const { return this->bounds; };
   const Tdd &getXMinMax() const { return std::get<0>(this->bounds); };
   const Tdd &getYMinMax() const { return std::get<1>(this->bounds); };
   const Tdd &getZMinMax() const { return std::get<2>(this->bounds); };
   double getVolume() const {
      return (std::get<1>(std::get<0>(this->bounds)) - std::get<0>(std::get<0>(this->bounds))) * (std::get<1>(std::get<1>(this->bounds)) - std::get<0>(std::get<1>(this->bounds))) * (std::get<1>(std::get<2>(this->bounds)) - std::get<0>(std::get<2>(this->bounds)));
   };
   double getScale() const {
      auto [min, max] = Transpose(this->bounds);
      return Norm(max - min);
   };
   Tddd getCenter() const { return Mean(Transpose(this->bounds)); };
   T8Tddd getVertices() const {
      auto [X0, X1] = std::get<0>(this->bounds);
      auto [Y0, Y1] = std::get<1>(this->bounds);
      auto [Z0, Z1] = std::get<2>(this->bounds);
      // connectivity:  {{x4, x6, x7, x5}, {x0, x2, x6, x4}, {x2, x3, x7, x6}, {x3, x1, x5, x7}, {x0, x4, x5, x1}, {x0, x1, x3, x2}}
      return T8Tddd{{X0, Y0, Z0},   // 000, 0
                    {X0, Y0, Z1},   // 001, 1
                    {X0, Y1, Z0},   // 010, 2
                    {X0, Y1, Z1},   // 011, 3
                    {X1, Y0, Z0},   // 100, 4
                    {X1, Y0, Z1},   // 101, 5
                    {X1, Y1, Z0},   // 110, 6
                    {X1, Y1, Z1}};  // 111, 7
   };
   //% ---------------- キャスト定義 -------------- */
   //% 型が明示され，関数よりもわかりやすい．
   operator T8Tddd() const {
      auto [X0, X1] = std::get<0>(this->bounds);
      auto [Y0, Y1] = std::get<1>(this->bounds);
      auto [Z0, Z1] = std::get<2>(this->bounds);
      // connectivity:  {{x4, x6, x7, x5}, {x0, x2, x6, x4}, {x2, x3, x7, x6}, {x3, x1, x5, x7}, {x0, x4, x5, x1}, {x0, x1, x3, x2}}
      // return T8Tddd{{X0, Y0, Z0},   // 000, 0
      //               {X0, Y1, Z0},   // 010, 1 -> 2
      //               {X0, Y0, Z1},   // 001, 2 -> 1
      //               {X0, Y1, Z1},   // 011, 3
      //               {X1, Y0, Z0},   // 100, 4
      //               {X1, Y1, Z0},   // 110, 5 -> 6
      //               {X1, Y0, Z1},   // 101, 6 -> 5
      //               {X1, Y1, Z1}};  // 111, 7
      return T8Tddd{{X0, Y0, Z0},   // 000, 0
                    {X0, Y0, Z1},   // 001, 1
                    {X0, Y1, Z0},   // 010, 2
                    {X0, Y1, Z1},   // 011, 3
                    {X1, Y0, Z0},   // 100, 4
                    {X1, Y0, Z1},   // 101, 5
                    {X1, Y1, Z0},   // 110, 6
                    {X1, Y1, Z1}};  // 111, 7
   };

   operator T6T4Tddd() const {
      auto [X0, X1] = std::get<0>(this->bounds);
      auto [Y0, Y1] = std::get<1>(this->bounds);
      auto [Z0, Z1] = std::get<2>(this->bounds);
      return T6T4Tddd{T4Tddd{{X0, Y0, Z0}, {X1, Y0, Z0}, {X1, Y1, Z0}, {X0, Y1, Z0}},
                      T4Tddd{{X0, Y0, Z1}, {X1, Y0, Z1}, {X1, Y1, Z1}, {X0, Y1, Z1}},
                      T4Tddd{{X0, Y0, Z0}, {X1, Y0, Z0}, {X1, Y0, Z1}, {X0, Y0, Z1}},
                      T4Tddd{{X0, Y1, Z0}, {X1, Y1, Z0}, {X1, Y1, Z1}, {X0, Y1, Z1}},
                      T4Tddd{{X0, Y0, Z0}, {X0, Y0, Z1}, {X0, Y1, Z1}, {X0, Y1, Z0}},
                      T4Tddd{{X1, Y0, Z0}, {X1, Y0, Z1}, {X1, Y1, Z1}, {X1, Y1, Z0}}};
   };

   operator T12T3Tddd() const {
      auto [X0, X1] = std::get<0>(this->bounds);
      auto [Y0, Y1] = std::get<1>(this->bounds);
      auto [Z0, Z1] = std::get<2>(this->bounds);
      return T12T3Tddd{T3Tddd{{X0, Y0, Z0}, {X1, Y0, Z0}, {X1, Y1, Z0}}, T3Tddd{{X0, Y0, Z0}, {X1, Y1, Z0}, {X0, Y1, Z0}},
                       T3Tddd{{X0, Y0, Z1}, {X1, Y0, Z1}, {X1, Y1, Z1}}, T3Tddd{{X0, Y0, Z1}, {X1, Y1, Z1}, {X0, Y1, Z1}},
                       T3Tddd{{X0, Y0, Z0}, {X1, Y0, Z0}, {X1, Y0, Z1}}, T3Tddd{{X0, Y0, Z0}, {X1, Y0, Z1}, {X0, Y0, Z1}},
                       T3Tddd{{X0, Y1, Z0}, {X1, Y1, Z0}, {X1, Y1, Z1}}, T3Tddd{{X0, Y1, Z0}, {X1, Y1, Z1}, {X0, Y1, Z1}},
                       T3Tddd{{X0, Y0, Z0}, {X0, Y0, Z1}, {X0, Y1, Z1}}, T3Tddd{{X0, Y0, Z0}, {X0, Y1, Z1}, {X0, Y1, Z0}},
                       T3Tddd{{X1, Y0, Z0}, {X1, Y0, Z1}, {X1, Y1, Z1}}, T3Tddd{{X1, Y0, Z0}, {X1, Y1, Z1}, {X1, Y1, Z0}}};
   };
   //% -------------------------------------------- */
   std::tuple<CoordinateBounds,  // 0
              CoordinateBounds,  // 1
              CoordinateBounds,  // 2
              CoordinateBounds,  // 3
              CoordinateBounds,  // 4
              CoordinateBounds,  // 5
              CoordinateBounds,  // 6
              CoordinateBounds>  // 7
   to8Bounds() const {
      /*
      +-----+-----+ Y1
      |     |     |
      +-----+-----+ Yc
      |     |     |
      +-----+-----+ Y0
      X0    Xc    X1
      */
      auto [X0, X1] = std::get<0>(this->bounds);
      auto [Y0, Y1] = std::get<1>(this->bounds);
      auto [Z0, Z1] = std::get<2>(this->bounds);
      auto [Xc, Yc, Zc] = this->getCenter();
      return {CoordinateBounds({X0, Xc}, {Y0, Yc}, {Z0, Zc}),
              CoordinateBounds({Xc, X1}, {Y0, Yc}, {Z0, Zc}),
              CoordinateBounds({X0, Xc}, {Yc, Y1}, {Z0, Zc}),
              CoordinateBounds({Xc, X1}, {Yc, Y1}, {Z0, Zc}),
              CoordinateBounds({X0, Xc}, {Y0, Yc}, {Zc, Z1}),
              CoordinateBounds({Xc, X1}, {Y0, Yc}, {Zc, Z1}),
              CoordinateBounds({X0, Xc}, {Yc, Y1}, {Zc, Z1}),
              CoordinateBounds({Xc, X1}, {Yc, Y1}, {Zc, Z1})};
   };
   std::tuple<CoordinateBounds,  // 0
              CoordinateBounds,  // 1
              CoordinateBounds,  // 2
              CoordinateBounds,  // 3
              CoordinateBounds,  // 4
              CoordinateBounds,  // 5
              CoordinateBounds,  // 6
              CoordinateBounds>  // 7
   to8Bounds(const Tddd &center) const {
      /*
      +-----+-----+ Y1
      |     |     |
      +-----+-----+ Yc
      |     |     |
      +-----+-----+ Y0
      X0    Xc    X1
      */
      auto [X0, X1] = std::get<0>(this->bounds);
      auto [Y0, Y1] = std::get<1>(this->bounds);
      auto [Z0, Z1] = std::get<2>(this->bounds);
      auto [Xc, Yc, Zc] = center;
      return {CoordinateBounds({X0, Xc}, {Y0, Yc}, {Z0, Zc}), CoordinateBounds({Xc, X1}, {Y0, Yc}, {Z0, Zc}),
              CoordinateBounds({X0, Xc}, {Yc, Y1}, {Z0, Zc}), CoordinateBounds({Xc, X1}, {Yc, Y1}, {Z0, Zc}),
              CoordinateBounds({X0, Xc}, {Y0, Yc}, {Zc, Z1}), CoordinateBounds({Xc, X1}, {Y0, Yc}, {Zc, Z1}),
              CoordinateBounds({X0, Xc}, {Yc, Y1}, {Zc, Z1}), CoordinateBounds({Xc, X1}, {Yc, Y1}, {Zc, Z1})};
   };
   bool isInside(const Tddd &X) const {
      if ((std::get<0>(X) < std::get<0>(std::get<0>(this->bounds)) || std::get<1>(std::get<0>(this->bounds)) < std::get<0>(X)) ||
          (std::get<1>(X) < std::get<0>(std::get<1>(this->bounds)) || std::get<1>(std::get<1>(this->bounds)) < std::get<1>(X)) ||
          (std::get<2>(X) < std::get<0>(std::get<2>(this->bounds)) || std::get<1>(std::get<2>(this->bounds)) < std::get<2>(X)))
         return false;
      else
         return true;
   };
};
/* -------------------------------------------------------------------------- */
struct Sphere : public CoordinateBounds {
   Tddd center;
   double radius;
   Sphere(const Tddd &XIN, const double radiusIN = 0.)
       : CoordinateBounds(T3Tdd{{(std::get<0>(XIN) - radiusIN), (std::get<0>(XIN) + radiusIN)},
                                {(std::get<1>(XIN) - radiusIN), (std::get<1>(XIN) + radiusIN)},
                                {(std::get<2>(XIN) - radiusIN), (std::get<2>(XIN) + radiusIN)}}),
         center(XIN),
         radius(radiusIN){};
};
struct Triangle : public CoordinateBounds {
   T3Tddd verticies;
   Tddd normal;
   Tddd center;
   Tddd angles;
   double area;
   Triangle(const T3Tddd &XIN)
       : CoordinateBounds(XIN),
         verticies(XIN),
         center(Mean(XIN)),
         normal(TriangleNormal(verticies)),
         area(TriangleArea(verticies)),
         angles(TriangleAngles(verticies)){};
   Triangle(const Tddd &X0IN, const Tddd &X1IN, const Tddd &X2IN)
       : CoordinateBounds(X0IN, X1IN, X2IN),
         verticies({X0IN, X1IN, X2IN}),
         center((X0IN + X1IN + X2IN) / 3.),
         normal(TriangleNormal(verticies)),
         area(TriangleArea(verticies)),
         angles(TriangleAngles(verticies)){};
};
/* -------------------------------------------------------------------------- */
struct Tetrahedron : public CoordinateBounds {
   //     1,3,2
   //          3
   // 0,2,3  / | \ 0,3,1     --- 1 ---
   //       2--|--1           \      /
   // 0,1,2  \ | /             2    0
   //          0                 \/
   //
   // GEOMETRIC PROPERTIES
   T4Tddd verticies;
   T4Tddd normals;
   Tddd circumcenter;
   Tddd solidangles;  //いつかチェック
   double volume;
   //
   Tetrahedron(const T4Tddd &XIN)
       : CoordinateBounds(XIN),
         verticies(XIN),
         volume(TetrahedronVolume(XIN)),
         circumcenter(TetrahedronCircumCenter(XIN)),
         normals(TetrahedronNormals(XIN)){};
};
/* ------------------------------------------------------ */
std::ostream &operator<<(std::ostream &stream, const CoordinateBounds &bounds) {
   return (stream << bounds.bounds);
};
CoordinateBounds operator+(const CoordinateBounds &b0, const CoordinateBounds &b1) {
   auto [minX0, maxX0] = std::get<0 /*x*/>(b0.bounds);
   auto [minX1, maxX1] = std::get<0 /*x*/>(b1.bounds);
   auto [minY0, maxY0] = std::get<1 /*y*/>(b0.bounds);
   auto [minY1, maxY1] = std::get<1 /*y*/>(b1.bounds);
   auto [minZ0, maxZ0] = std::get<2 /*z*/>(b0.bounds);
   auto [minZ1, maxZ1] = std::get<2 /*z*/>(b1.bounds);
   return CoordinateBounds(
       (minX0 <= minX1 ? minX0 : minX1), (maxX0 >= maxX1 ? maxX0 : maxX1),
       (minY0 <= minY1 ? minY0 : minY1), (maxY0 >= maxY1 ? maxY0 : maxY1),
       (minZ0 <= minZ1 ? minZ0 : minZ1), (maxZ0 >= maxZ1 ? maxZ0 : maxZ1));
};
CoordinateBounds &operator+=(CoordinateBounds &b0, const CoordinateBounds &b1) {
   auto [minX0, maxX0] = std::get<0 /*x*/>(b0.bounds);
   auto [minX1, maxX1] = std::get<0 /*x*/>(b1.bounds);
   auto [minY0, maxY0] = std::get<1 /*y*/>(b0.bounds);
   auto [minY1, maxY1] = std::get<1 /*y*/>(b1.bounds);
   auto [minZ0, maxZ0] = std::get<2 /*z*/>(b0.bounds);
   auto [minZ1, maxZ1] = std::get<2 /*z*/>(b1.bounds);
   b0.bounds = {
       {(minX0 <= minX1 ? minX0 : minX1), (maxX0 >= maxX1 ? maxX0 : maxX1)},
       {(minY0 <= minY1 ? minY0 : minY1), (maxY0 >= maxY1 ? maxY0 : maxY1)},
       {(minZ0 <= minZ1 ? minZ0 : minZ1), (maxZ0 >= maxZ1 ? maxZ0 : maxZ1)}};
   return b0;
};
inline CoordinateBounds::CoordinateBounds(const std::vector<T3Tddd> &V_X) {
   CoordinateBounds ret;
   for (const auto &X : V_X) ret += CoordinateBounds(X);
   this->bounds = ret.bounds;
};
/* ------------------------------------------------------ */
struct IntersectionSphereLine {
   /*
      /Users/tomoaki/Dropbox/markdown/mathematica/非構造格子/三角形と球の干渉.nb
      X = (p0,p1).(t,1-t) = (p0-p1)*t + p1 であることから，
      t = 0 -> X = p0
      t = 1 -> X = p1
      となる
   */
   const Tddd p0, p1;
   double t;
   Tddd X;
   double distance;
   bool isIntersecting;
   IntersectionSphereLine(const Tddd &center, const double radius, const T2Tddd &p01)
       : p0(std::get<0>(p01)),
         p1(std::get<1>(p01)),
         t(Dot(center - p1, p0 - p1) / Dot(p0 - p1, p0 - p1)),
         X(p0 * t + p1 * (1 - t)),
         distance(Norm(X - center)),
         isIntersecting(distance <= radius && 0. <= t && t <= 1.) {
      if (!isIntersecting) {
         //線分と干渉しない場合でも，球に最も近い線分上の点を返すようにする．
         //これを実行しない場合，線分ではなく，直線上の点を返すことになる．
         if (Norm(p0 - center) < Norm(p1 - center)) {
            this->t = 1;
            this->X = p0;
            this->distance = Norm(p0 - center);
         } else {
            this->t = 0;
            this->X = p1;
            this->distance = Norm(p1 - center);
         }
      }
   };
   IntersectionSphereLine(const double radius, const T2Tddd &p01)
       : p0(std::get<0>(p01)),
         p1(std::get<1>(p01)),
         t(Dot(-p1, p0 - p1) / Dot(p0 - p1, p0 - p1)),
         X(p0 * t + p1 * (1 - t)),
         distance(Norm(X)),
         isIntersecting(distance <= radius && 0. <= t && t <= 1.) {
      if (!isIntersecting) {
         //線分と干渉しない場合でも，球に最も近い線分上の点を返すようにする．
         //これを実行しない場合，線分ではなく，直線上の点を返すことになる．
         if (Norm(p0) < Norm(p1)) {
            this->t = 1;
            this->X = p0;
            this->distance = Norm(p0);
         } else {
            this->t = 0;
            this->X = p1;
            this->distance = Norm(p1);
         }
      }
   };
};
/* ------------------------------------------------------ */
struct IntersectionTriangles {
   /*
   /Users/tomoaki/Dropbox/markdown/mathematica/非構造格子/三角形と球の干渉.nb
   */
   T3Tddd P0, P1;
   T2Tddd L;
   bool isIntersecting;
   double eps = 1E-13;
   double a0, b0, a1, b1, min0, max0, min1, max1, deno0, deno1;
   double I10InT, I11InT, I00InT, I01InT;
   Tddd I00, I01, I10, I11;
   double lmitMin(const double a0, const double b0) const {
      Tddd ret = {0., 0., 0.};
      if (std::abs(1 - a0) > eps) {
         if ((1 - a0) > 0)
            std::get<0>(ret) = -b0 / (1 - a0);
         else
            std::get<0>(ret) = (1 - b0) / (1 - a0);
      }

      if (std::abs(-a0) > eps) {
         if (-a0 > 0)
            std::get<1>(ret) = b0 / a0;
         else
            std::get<1>(ret) = -(1 - b0) / a0;
      }
      return Max(ret);
   };

   double lmitMax(const double a0, const double b0) const {
      Tddd ret = {1., 1., 1.};
      if (std::abs(1 - a0) > eps) {
         if ((1 - a0) > 0)
            std::get<0>(ret) = (1 - b0) / (1 - a0);
         else
            std::get<0>(ret) = -b0 / (1 - a0);
      }

      if (std::abs(-a0) > eps) {
         if (-a0 > 0)
            std::get<1>(ret) = -(1 - b0) / a0;
         else
            std::get<1>(ret) = b0 / a0;
      }
      return Min(ret);  //この内最も小さいものが，上限となる
   };

   Tddd p00, p01, p02;
   Tddd p10, p11, p12;
   Tddd normalP0, normalP1;
   /* -------------------------------------------------------------- */
   IntersectionTriangles(const T3Tddd &P0_IN, const T3Tddd &P1_IN)
       : P0(P0_IN), P1(P1_IN), isIntersecting(false), p00(std::get<0>(P0)), p01(std::get<1>(P0)), p02(std::get<2>(P0)), p10(std::get<0>(P1)), p11(std::get<1>(P1)), p12(std::get<2>(P1)), normalP0(Normalize(Cross(p01 - p00, p02 - p00))), normalP1(Normalize(Cross(p11 - p10, p12 - p10))) {
      // if (IntersectQ(CoordinateBounds(P0), CoordinateBounds(P0)))
      // 	return;
      // auto [p00, p01, p02] = P0;
      // auto [p10, p11, p12] = P1;
      // auto normalP0 = Normalize(Cross(p01 - p00, p02 - p00));
      // auto normalP1 = Normalize(Cross(p11 - p10, p12 - p10));
      if (std::abs(std::abs(Dot(normalP0, normalP1)) - 1) < eps) {
         // Print("干渉線を定義することはできない");
         // Print("線を共有するような場合はありえる");
         return;
      }
      deno0 = Dot(normalP1, p01 - p02);
      deno1 = Dot(normalP0, p11 - p12);

      if (std::abs(deno0) < eps)
         for (auto i = 1; i < 3; ++i) {
            P0 = RotateLeft(P0, i);
            p00 = std::get<0>(P0);
            p01 = std::get<1>(P0);
            p02 = std::get<2>(P0);
            deno0 = Dot(normalP1, p01 - p02);
            if (!(std::abs(deno0) < eps)) break;
            if (i == 2)
               throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
         }

      if (std::abs(deno1) < eps)
         for (auto i = 1; i < 3; ++i) {
            P1 = RotateLeft(P1, i);
            p10 = std::get<0>(P1);
            p11 = std::get<1>(P1);
            p12 = std::get<2>(P1);
            deno1 = Dot(normalP0, p11 - p12);
            if (!(std::abs(deno1) < eps)) break;
            if (i == 2)
               throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
         }

      a0 = Dot(normalP1, p00 - p02) / deno0;
      b0 = Dot(normalP1, p12 - p02) / deno0;
      a1 = Dot(normalP0, p10 - p12) / deno1;
      b1 = Dot(normalP0, p02 - p12) / deno1;

      min0 = lmitMin(a0, b0);
      max0 = lmitMax(a0, b0);
      min1 = lmitMin(a1, b1);
      max1 = lmitMax(a1, b1);

      auto t = min0;
      I00 = Dot({t, b0 - a0 * t, 1 - t - (b0 - a0 * t)}, P0);
      t = max0;
      I01 = Dot({t, b0 - a0 * t, 1 - t - (b0 - a0 * t)}, P0);
      t = min1;
      I10 = Dot({t, b1 - a1 * t, 1 - t - (b1 - a1 * t)}, P1);
      t = max1;
      I11 = Dot({t, b1 - a1 * t, 1 - t - (b1 - a1 * t)}, P1);

      auto tmp = (I10 - I00) / (I01 - I00);
      I10InT = isFinite(std::get<0>(tmp)) ? std::get<0>(tmp) : (isFinite(std::get<1>(tmp)) ? std::get<1>(tmp) : std::get<2>(tmp));
      auto I10IsInside = Between(I10InT, {0, 1});

      tmp = (I11 - I00) / (I01 - I00);
      I11InT = isFinite(std::get<0>(tmp)) ? std::get<0>(tmp) : (isFinite(std::get<1>(tmp)) ? std::get<1>(tmp) : std::get<2>(tmp));
      auto I11IsInside = Between(I11InT, {0, 1});

      tmp = (I00 - I10) / (I11 - I10);
      I00InT = isFinite(std::get<0>(tmp)) ? std::get<0>(tmp) : (isFinite(std::get<1>(tmp)) ? std::get<1>(tmp) : std::get<2>(tmp));
      auto I00IsInside = Between(I00InT, {0, 1});

      tmp = (I01 - I10) / (I11 - I10);
      I01InT = isFinite(std::get<0>(tmp)) ? std::get<0>(tmp) : (isFinite(std::get<1>(tmp)) ? std::get<1>(tmp) : std::get<2>(tmp));
      auto I01IsInside = Between(I01InT, {0, 1});

      if (I10IsInside) {
         if (I11IsInside && Norm(I10 - I11) > eps) {
            L = {I10, I11};
            isIntersecting = true;
         }
         if (I01IsInside && Norm(I10 - I01) > eps) {
            L = {I10, I01};
            isIntersecting = true;
         }
         if (I00IsInside && Norm(I10 - I00) > eps) {
            L = {I10, I00};
            isIntersecting = true;
         }
      }
      if (I11IsInside) {
         if (I00IsInside && Norm(I11 - I00) > eps) {
            L = {I11, I00};
            isIntersecting = true;
         }
         if (I01IsInside && Norm(I11 - I01) > eps) {
            L = {I11, I01};
            isIntersecting = true;
         }
      }
      if (I00IsInside && I01IsInside && Norm(I00 - I01) > eps) {
         L = {I00, I01};
         isIntersecting = true;
      }
   };
};
/* ------------------------------------------------------ */
// struct IntersectionSphereTriangle {
//    double eps = 1E-14;
//    Tddd center;
//    Tddd n;
//    double n0, n1, n2;
//    double x00, x01, x02, x10, x11, x12, x20, x21, x22;
//    double denominator;
//    double scale, t0, t1;
//    double distance;
//    Tddd X;
//    const IntersectionSphereLine intxnL0;
//    const IntersectionSphereLine intxnL1;
//    const IntersectionSphereLine intxnL2;
//    bool isIntersecting;
//    bool isIntersectingInsideTriangle;
//    bool isIntersectingEdgeTriangle;
//    IntersectionSphereTriangle(const Tddd &centerIN, const double radius, const T3Tddd &X0X1X2)
//        : center(centerIN),
//          x00(std::get<0>(std::get<0>(X0X1X2)) - std::get<0>(centerIN)),
//          x01(std::get<1>(std::get<0>(X0X1X2)) - std::get<1>(centerIN)),
//          x02(std::get<2>(std::get<0>(X0X1X2)) - std::get<2>(centerIN)),
//          x10(std::get<0>(std::get<1>(X0X1X2)) - std::get<0>(centerIN)),
//          x11(std::get<1>(std::get<1>(X0X1X2)) - std::get<1>(centerIN)),
//          x12(std::get<2>(std::get<1>(X0X1X2)) - std::get<2>(centerIN)),
//          x20(std::get<0>(std::get<2>(X0X1X2)) - std::get<0>(centerIN)),
//          x21(std::get<1>(std::get<2>(X0X1X2)) - std::get<1>(centerIN)),
//          x22(std::get<2>(std::get<2>(X0X1X2)) - std::get<2>(centerIN)),
//          n(Normalize(Cross(std::get<0>(X0X1X2) - std::get<2>(X0X1X2), std::get<1>(X0X1X2) - std::get<2>(X0X1X2)))),
//          n0(std::get<0>(n)),
//          n1(std::get<1>(n)),
//          n2(std::get<2>(n)),
//          denominator(n2 * (x00 * x11 - x11 * x20 + x01 * (-x10 + x20) - x00 * x21 + x10 * x21) + n1 * (x02 * x10 - x00 * x12 - x02 * x20 + x12 * x20 + x00 * x22 - x10 * x22) + n0 * (-(x02 * x11) + x01 * x12 + x02 * x21 - x12 * x21 - x01 * x22 + x11 * x22)),
//          t0((-(n2 * x11 * x20) + n1 * x12 * x20 + n2 * x10 * x21 - n0 * x12 * x21 - n1 * x10 * x22 + n0 * x11 * x22) / denominator),
//          t1((n2 * x01 * x20 - n1 * x02 * x20 - n2 * x00 * x21 + n0 * x02 * x21 + n1 * x00 * x22 - n0 * x01 * x22) / denominator),
//          scale((-(x02 * x11 * x20) + x01 * x12 * x20 + x02 * x10 * x21 - x00 * x12 * x21 - x01 * x10 * x22 + x00 * x11 * x22) / denominator),
//          intxnL0(centerIN, 1E+80, T2Tddd{std::get<0>(X0X1X2), std::get<1>(X0X1X2)}),
//          intxnL1(centerIN, 1E+80, T2Tddd{std::get<1>(X0X1X2), std::get<2>(X0X1X2)}),
//          intxnL2(centerIN, 1E+80, T2Tddd{std::get<2>(X0X1X2), std::get<0>(X0X1X2)}),
//          distance(std::abs(scale)),
//          X(scale * n + centerIN),
//          isIntersecting(std::abs(scale) <= radius && Between(t0, {0, 1}) && Between(t1, {0, 1}) && Between(1 - t0 - t1, {0, 1})),
//          isIntersectingInsideTriangle(isIntersecting),
//          isIntersectingEdgeTriangle(!isIntersecting) {
//       /*
//        *   P = (p0, p1, p2)
//        *   P.T = p0*t0 + p1*t1 + p2*t2
//        *
//        *                           p0 <- T=(1,0,0)
//        *                          *
//        *     P.(t0,1-t0,0) ->    / \   <- P.(1-t2,0,t2) = p0*(1-t2) + p2*t2
//        *  =>(p0,p1).(t0,1-t0)   /   \     => (p0,p2).(1-t2,t2)
//        *                       *-----*    => (p0,p2).(1-t2,t2)
//        *       T=(0,1,0) -> p1           p2 <- T=(0,0,1)
//        *                     P.(0,t1,1-t1)
//        *                   => (p1,p2).(t1,1-t1)
//        */
//       if (!isIntersectingInsideTriangle) {
//          if (intxnL0.isIntersecting && Norm(intxnL1.X - center) <= distance) {
//             this->t0 = intxnL0.t;
//             this->t1 = 1 - intxnL0.t;
//             this->X = intxnL0.X;
//             this->distance = Norm(X - center);
//             this->isIntersectingInsideTriangle = false;
//             this->isIntersectingEdgeTriangle = this->isIntersecting = true;
//          }
//          if (intxnL1.isIntersecting && Norm(intxnL1.X - center) <= distance) {
//             this->t0 = 0;
//             this->t1 = intxnL1.t;
//             this->X = intxnL1.X;
//             this->distance = Norm(X - center);
//             this->isIntersectingInsideTriangle = false;
//             this->isIntersectingEdgeTriangle = this->isIntersecting = true;
//          }
//          if (intxnL2.isIntersecting && Norm(intxnL2.X - center) <= distance) {
//             this->t0 = 1 - intxnL2.t;
//             this->t1 = 0;
//             this->X = intxnL2.X;
//             this->distance = Norm(X - center);
//             this->isIntersectingInsideTriangle = false;
//             this->isIntersectingEdgeTriangle = this->isIntersecting = true;
//          }
//       }
//    };
// };
/* ------------------------------------------------------ */
struct IntersectionSphereTriangle {
   double eps = 1E-8;
   Tddd X0, X1, X2, center;
   Tddd n;
   double n0, n1, n2;
   double x00, x01, x02, x10, x11, x12, x20, x21, x22;
   double denominator;
   double scale, t0, t1;
   double distance;
   Tddd X;
   const IntersectionSphereLine intxnL0;
   const IntersectionSphereLine intxnL1;
   const IntersectionSphereLine intxnL2;
   bool isIntersecting;
   bool isIntersectingInsideTriangle;
   bool isIntersectingEdgeTriangle;
   //
   // double F(const double t0, const double t1) {
   //    auto v = Dot({t0, t1, 1}, this->X0X1X2);
   //    return Dot(v, v);
   // };
   // double dFdt(const double t0, const double t1) {
   //    auto v = Dot({t0, t1, 1}, this->X0X1X2);
   //    auto dvdt0 = Dot({1, 0, 0}, this->X0X1X2);
   //    auto dvdt1 = Dot({0, 1, 0}, this->X0X1X2);
   //    return Dot(v, v);
   // };
   IntersectionSphereTriangle(const Tddd &centerIN, const double radius, const T3Tddd &X0X1X2)
       : center(centerIN),
         X0(std::get<0>(X0X1X2) - centerIN),
         X1(std::get<1>(X0X1X2) - centerIN),
         X2(std::get<2>(X0X1X2) - centerIN),
         x00(std::get<0>(X0)),
         x01(std::get<1>(X0)),
         x02(std::get<2>(X0)),
         x10(std::get<0>(X1)),
         x11(std::get<1>(X1)),
         x12(std::get<2>(X1)),
         x20(std::get<0>(X2)),
         x21(std::get<1>(X2)),
         x22(std::get<2>(X2)),
         n(Normalize(Cross(X1 - X0, X2 - X0))),
         n0(std::get<0>(n)),
         n1(std::get<1>(n)),
         n2(std::get<2>(n)),
         denominator(n2 * (x00 * x11 - x11 * x20 + x01 * (-x10 + x20) - x00 * x21 + x10 * x21) + n1 * (x02 * x10 - x00 * x12 - x02 * x20 + x12 * x20 + x00 * x22 - x10 * x22) + n0 * (-(x02 * x11) + x01 * x12 + x02 * x21 - x12 * x21 - x01 * x22 + x11 * x22)),
         t0((-(n2 * x11 * x20) + n1 * x12 * x20 + n2 * x10 * x21 - n0 * x12 * x21 - n1 * x10 * x22 + n0 * x11 * x22) / denominator),
         t1((n2 * x01 * x20 - n1 * x02 * x20 - n2 * x00 * x21 + n0 * x02 * x21 + n1 * x00 * x22 - n0 * x01 * x22) / denominator),
         scale((-(x02 * x11 * x20) + x01 * x12 * x20 + x02 * x10 * x21 - x00 * x12 * x21 - x01 * x10 * x22 + x00 * x11 * x22) / denominator),
         intxnL0(1E+20, T2Tddd{X0, X1}),
         intxnL1(1E+20, T2Tddd{X1, X2}),
         intxnL2(1E+20, T2Tddd{X2, X0}),
         distance(std::abs(scale)),
         X(scale * n + centerIN),
         isIntersecting(distance <= radius && (0 <= t0 && t0 <= 1) && (0 <= t1 && t1 <= 1) && (0 <= (1 - t0 - t1) && (1 - t0 - t1) <= 1)),
         isIntersectingInsideTriangle(isIntersecting),
         isIntersectingEdgeTriangle(false) {
      /*
       *   P = (p0, p1, p2)
       *   P.T = p0*t0 + p1*t1 + p2*t2
       *
       *                           p0 <- T=(1,0,0)
       *                          *
       *     P.(t0,1-t0,0) ->    / \   <- P.(1-t2,0,t2) = p0*(1-t2) + p2*t2
       *  =>(p0,p1).(t0,1-t0) # /   \     => (p0,p2).(1-t2,t2)
       *                       *-----*    => (p2,p0).(t2,1-t2) #
       *       T=(0,1,0) -> p1           p2 <- T=(0,0,1)
       *                     P.(0,t1,1-t1)
       *                   => (p1,p2).(t1,1-t1) #
       */
      if (!isIntersectingInsideTriangle) {
         double online_dist = 1E+10;
         if (intxnL0.isIntersecting /*distanceのチェックは不要．必ず上書きする*/) {
            this->t0 = intxnL0.t;
            this->t1 = 1 - intxnL0.t;
            this->X = intxnL0.X + center;
            online_dist = this->distance = intxnL0.distance;
            this->isIntersectingInsideTriangle = false;
            this->isIntersectingEdgeTriangle = this->isIntersecting = true;
         }
         if (intxnL1.isIntersecting && intxnL1.distance <= online_dist) {
            this->t0 = 0;
            this->t1 = intxnL1.t;
            this->X = intxnL1.X + center;
            online_dist = this->distance = intxnL1.distance;
            this->isIntersectingInsideTriangle = false;
            this->isIntersectingEdgeTriangle = this->isIntersecting = true;
         }
         if (intxnL2.isIntersecting && intxnL2.distance <= online_dist) {
            this->t0 = 1 - intxnL2.t;
            this->t1 = 0;
            this->X = intxnL2.X + center;
            online_dist = this->distance = intxnL2.distance;
            this->isIntersectingInsideTriangle = false;
            this->isIntersectingEdgeTriangle = this->isIntersecting = true;
         }
      }
      if (!isIntersectingInsideTriangle || !isIntersectingEdgeTriangle) {
         if (Norm(X0) < this->distance) {
            Norm(this->X = X0 + center);
            this->distance = Norm(X0);
         }
         if (Norm(X1) < this->distance) {
            Norm(this->X = X1 + center);
            this->distance = Norm(X1);
         }
         if (Norm(X2) < this->distance) {
            Norm(this->X = X2 + center);
            this->distance = Norm(X2);
         }
      }
   };
};
/* ------------------------------------------------------ */
struct IntersectionSphereTriangleLimitedToNormalRegion {
   const double eps = 1E-13;
   const T3Tddd P;
   const Tddd center;
   //
   const Tddd p0, p1, p2;
   const Tddd n;
   const double nx, ny, nz;
   const double p02x, p02y, p02z, p12x, p12y, p12z;
   const double determ;
   const T3Tddd mat;
   const Tddd ans;
   const double t0, t1, scale;
   const Tddd X;
   const bool isIntersecting;
   IntersectionSphereTriangleLimitedToNormalRegion(const Tddd &centerIN, const double radius, const T3Tddd &P_IN)
       : P(P_IN),
         center(centerIN),
         p0(std::get<0>(P_IN)),
         p1(std::get<1>(P_IN)),
         p2(std::get<2>(P_IN)),
         n(Normalize(Cross(p1 - p0, p2 - p0))),
         nx(std::get<0>(n)),
         ny(std::get<1>(n)),
         nz(std::get<2>(n)),
         p02x(std::get<0>(p0 - p2)),
         p02y(std::get<1>(p0 - p2)),
         p02z(std::get<2>(p0 - p2)),
         p12x(std::get<0>(p1 - p2)),
         p12y(std::get<1>(p1 - p2)),
         p12z(std::get<2>(p1 - p2)),
         determ(-(nz * p02y * p12x) + ny * p02z * p12x + nz * p02x * p12y - nx * p02z * p12y - ny * p02x * p12z + nx * p02y * p12z),
         mat({{nz * p12y - ny * p12z, -(nz * p02y) + ny * p02z, p02z * p12y - p02y * p12z},
              {-(nz * p12x) + nx * p12z, nz * p02x - nx * p02z, -(p02z * p12x) + p02x * p12z},
              {ny * p12x - nx * p12y, -(ny * p02x) + nx * p02y, p02y * p12x - p02x * p12y}}),
         ans(Dot(center - p2, mat / determ)),
         t0(std::get<0>(ans)),
         t1(std::get<1>(ans)),
         scale(std::get<2>(ans)),
         X(scale * n + center),
         isIntersecting((std::abs(determ) > eps) && (std::abs(scale) <= radius && (0 <= t0 && t0 <= 1) && (0 <= t1 && t1 <= 1) && (0 <= (1 - t0 - t1) && (1 - t0 - t1) <= 1))){};
   Tddd getNearestX() const {
      if (this->isIntersecting)
         return this->X;
      else {
         Tddd X_ = {1E+50, 1E+50, 1E+50};
         double mindistance = 1E+50;
         auto intxnL0 = IntersectionSphereLine(
             center, 1E+50, T2Tddd{std::get<0>(P), std::get<1>(P)});
         auto intxnL1 = IntersectionSphereLine(
             center, 1E+50, T2Tddd{std::get<1>(P), std::get<2>(P)});
         auto intxnL2 = IntersectionSphereLine(
             center, 1E+50, T2Tddd{std::get<2>(P), std::get<0>(P)});
         if (intxnL0.isIntersecting) {
            X_ = intxnL0.X;
            mindistance = Norm(X_ - center);
         }
         if (intxnL1.isIntersecting) {
            if (Norm(intxnL1.X - center) < mindistance) {
               X_ = intxnL1.X;
               mindistance = Norm(X_ - center);
            }
         }
         if (intxnL2.isIntersecting) {
            if (Norm(intxnL2.X - center) < mindistance) {
               X_ = intxnL2.X;
               mindistance = Norm(X_ - center);
            }
         }
         if (Norm(std::get<0>(P) - center) < mindistance) {
            X_ = std::get<0>(P);
            mindistance = Norm(X_ - center);
         }
         if (Norm(std::get<1>(P) - center) < mindistance) {
            X_ = std::get<1>(P);
            mindistance = Norm(X_ - center);
         }
         if (Norm(std::get<2>(P) - center) < mindistance) {
            X_ = std::get<2>(P);
         }
         return X_;
      }
   };
};
/* -------------------------------------------------------------------------- */
bool IntersectQ(const T3Tdd &b0, const T3Tdd &b1) {
   return !((std::get<0>(std::get<0>(b0)) > std::get<0>(std::get<0>(b1)) && std::get<0>(std::get<0>(b0)) > std::get<1>(std::get<0>(b1)) /*1のxの最大最小が，0のxの最小よりも小さい*/) ||
            (std::get<0>(std::get<1>(b0)) > std::get<0>(std::get<1>(b1)) && std::get<0>(std::get<1>(b0)) > std::get<1>(std::get<1>(b1)) /*1のyの最大最小が，0のyの最小よりも小さい*/) ||
            (std::get<0>(std::get<2>(b0)) > std::get<0>(std::get<2>(b1)) && std::get<0>(std::get<2>(b0)) > std::get<1>(std::get<2>(b1)) /*1のzの最大最小が，0のzの最小よりも小さい*/) ||
            (std::get<1>(std::get<0>(b0)) < std::get<0>(std::get<0>(b1)) && std::get<1>(std::get<0>(b0)) < std::get<1>(std::get<0>(b1)) /*1のxの最大最小が，0のxの最大よりも大きい*/) ||
            (std::get<1>(std::get<1>(b0)) < std::get<0>(std::get<1>(b1)) && std::get<1>(std::get<1>(b0)) < std::get<1>(std::get<1>(b1)) /*1のyの最大最小が，0のyの最大よりも大きい*/) ||
            (std::get<1>(std::get<2>(b0)) < std::get<0>(std::get<2>(b1)) && std::get<1>(std::get<2>(b0)) < std::get<1>(std::get<2>(b1)) /*1のzの最大最小が，0のzの最大よりも大きい*/) /*これがtrueの場合，逆にhitなし*/);
};
bool IntersectQ(const CoordinateBounds &b0, const CoordinateBounds &b1) { return IntersectQ(b0.bounds, b1.bounds); };

const Tddd &Normal(const geometry::Triangle &triangle) {
   return triangle.normal;
};
double NormalDistance(const geometry::Triangle &T, const Tddd &X) {
   return Norm(Dot(Normal(T), T.center - X));
};
Tddd vectorToTriangle(const geometry::Triangle &T, const Tddd &a) {
   // aからTまでの最短ベクトル
   auto n = Normal(T);
   return n * Dot(n, std::get<0>(T.X) - a);
};
int IntersectQ(const geometry::Sphere &sphere, const geometry::Triangle &triangle) {
   if (NormalDistance(triangle, sphere.X) < sphere.radius)
      return false;
   else
      return true;
};
bool IntersectQ(const Tddd &x, const double &r, const T2Tddd &ab) {
   auto [a, b] = ab;
   // return std::abs(Dot(a - b, Normalize(Cross(a - b, Cross(a - x, b - x))))) <= r;
   auto b_a = b - a;
   double t = (Dot(b_a, x) - Dot(a, b)) / Dot(b_a, b_a);
   return Norm(a + b_a * t - x) <= r;
};

//* -------------------------------------------------------------------------- */
//*                                 IntersectQ                                 */
//* -------------------------------------------------------------------------- */

Tddd t0_t1_alpha(const T3Tddd &p0p1p2, const Tddd &X) {
   //@ ３点の張る面　と　１点　の関係を調べる
   // a = (p0,p1,p2).(t0,t1,1-t0,t1)は，三角形が張る面上で，Xに最も近い点
   // この位置aから，alpha*nだけ移動した位置にXがある．n=(nx, ny, nz)は，三角形がつくる単位法線ベクトル
   auto [p0, p1, p2] = p0p1p2;
   auto [nx, ny, nz] = Normalize(Cross(p1 - p0, p2 - p0));  //三角形がつくる単位法線ベクトル
   auto [p02x, p02y, p02z] = p0 - p2;
   auto [p12x, p12y, p12z] = p1 - p2;
   return Dot(X - p2, {{nz * p12y - ny * p12z, -(nz * p02y) + ny * p02z, p02z * p12y - p02y * p12z},
                       {-(nz * p12x) + nx * p12z, nz * p02x - nx * p02z, -(p02z * p12x) + p02x * p12z},
                       {ny * p12x - nx * p12y, -(ny * p02x) + nx * p02y, p02y * p12x - p02x * p12y}}) /
          (-(nz * p02y * p12x) + ny * p02z * p12x + nz * p02x * p12y - nx * p02z * p12y - ny * p02x * p12z + nx * p02y * p12z);
};

// sphere - point
bool IntersectQ(const Sphere &sp, const Tddd &a) { return Norm(a - sp.center) < sp.radius; };
bool IntersectQ(const Tddd &a, const Sphere &sp) { return Norm(a - sp.center) < sp.radius; };
// sphere - line
bool IntersectQ(const Sphere &sp, const T2Tddd &ab) {
   if (IntersectQ(sp, CoordinateBounds(ab))) {
      auto [a, b] = ab;
      auto b_a = b - a;
      double t = (Dot(b_a, sp.center) - Dot(a, b)) / Dot(b_a, b_a);
      return Between(t, {0., 1.}) && Norm(a + b_a * t - sp.center) <= sp.radius;
   } else
      return false;
};
bool IntersectQ(const T2Tddd &ab, const Sphere &sp) { return IntersectQ(sp, ab); };
// sphere - triangle
bool IntersectQ(const Sphere &sp, const T3Tddd &abc) {
   if (IntersectQ(sp, CoordinateBounds(abc))) {
      auto [p0, p1, p2] = abc;
      if (IntersectQ(sp, p0) || IntersectQ(sp, p1) || IntersectQ(sp, p2)) {
         return true;  //球体と点の干渉
      } else if (IntersectQ(sp, {p0, p1}) || IntersectQ(sp, {p1, p2}) || IntersectQ(sp, {p2, p0})) {
         return true;  //球体と線の干渉
      } else {
         auto [t0, t1, alpha] = t0_t1_alpha(abc, sp.center);
         return (std::abs(alpha) <= sp.radius && Between(t0, {0., 1.}) && Between(t1, {0., 1.}) && Between(1. - t0 - t1, {0., 1.}));
         //球体と三角形内部の干渉
      }
   } else
      return false;
};
// cube - sphere
bool IntersectQ(const T3Tdd &b, const Sphere &sp) {
   CoordinateBounds B(b);
   if (B.isInside(sp.center))
      return true;
   else {
      if (IntersectQ(B, sp)) {
         auto [X0, X1, X2, X3, X4, X5, X6, X7] = (T8Tddd)B;
         if (Norm(X0 - sp.center) < sp.radius || Norm(X1 - sp.center) < sp.radius || Norm(X2 - sp.center) < sp.radius || Norm(X3 - sp.center) < sp.radius ||
             Norm(X4 - sp.center) < sp.radius || Norm(X5 - sp.center) < sp.radius || Norm(X6 - sp.center) < sp.radius || Norm(X7 - sp.center) < sp.radius)
            return true;
         else {
            //３点が張る面との干渉チェック
            auto [X0, X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X11] = (T12T3Tddd)B;
            return (IntersectQ(sp, X0) || IntersectQ(sp, X1) || IntersectQ(sp, X2) || IntersectQ(sp, X3) || IntersectQ(sp, X4) || IntersectQ(sp, X5) ||
                    IntersectQ(sp, X6) || IntersectQ(sp, X7) || IntersectQ(sp, X7) || IntersectQ(sp, X8) || IntersectQ(sp, X9) || IntersectQ(sp, X10) || IntersectQ(sp, X11));
         }
      } else
         return false;
   }
};
//
bool IntersectQ(const T3Tdd &b, const Tddd &X) {
   CoordinateBounds B(b);
   return B.isInside(X);
};
//* -------------------------------------------------------------------------- */
Tddd Angles(const geometry::Triangle &triangle) {
   auto [a, b, c] = triangle.X;
   return std::make_tuple(VectorAngle(b - a, c - a), VectorAngle(c - b, a - b), VectorAngle(a - c, b - c));
};
double Area(const geometry::Triangle &triangle) {
   auto [a, b, c] = triangle.X;
   auto A = Norm(a - c);
   auto B = Norm(b - a);
   auto C = Norm(c - b);
   auto s = 0.5 * (A + B + C);
   return std::sqrt(s * (s - A) * (s - B) * (s - C));
};
/* ------------------------------------------------------ */
double scalefactorToReach(const geometry::Line &line, const geometry::Triangle &triangle) {
   auto [a, b] = line.X;
   auto [p0, p1, p2] = triangle.X;
   //オーダーが匹敵する物を選ぶ
   double log_b_a = log10(Norm(b - a));
   double diff0 = std::abs(log10(Norm(p0 - a) - log_b_a));
   double diff1 = std::abs(log10(Norm(p1 - a) - log_b_a));
   double diff2 = std::abs(log10(Norm(p2 - a) - log_b_a));
   Tddd n = Normal(triangle);
   if (diff0 < diff1 && diff0 < diff2)
      return Dot(p0 - a, n) / Dot(b - a, n);
   else if (diff1 < diff0 && diff1 < diff2)
      return Dot(p1 - a, n) / Dot(b - a, n);
   else
      return Dot(p2 - a, n) / Dot(b - a, n);
};
/* ------------------------------------------------------ */
int IntersectQ(const geometry::Line &line, const geometry::Triangle &triangle) {
   auto [lA, lB] = line.X;
   auto [tA, tB, tC] = triangle.X;
   auto [lAx, lAy, lAz] = lA;
   auto [lBx, lBy, lBz] = lB;
   auto [tAx, tAy, tAz] = tA;
   auto [tBx, tBy, tBz] = tB;
   auto [tCx, tCy, tCz] = tC;
   if ((tAx > lAx && tBx > lAx && tCx > lAx && tAx > lBx && tBx > lBx && tCx > lBx) ||
       (tAx < lAx && tBx < lAx && tCx < lAx && tAx < lBx && tBx < lBx && tCx < lBx) ||
       (tAy > lAy && tBy > lAy && tCy > lAy && tAy > lBy && tBy > lBy && tCy > lBy) ||
       (tAy < lAy && tBy < lAy && tCy < lAy && tAy < lBy && tBy < lBy && tCy < lBy) ||
       (tAz > lAz && tBz > lAz && tCz > lAz && tAz > lBz && tBz > lBz && tCz > lBz) ||
       (tAz < lAz && tBz < lAz && tCz < lAz && tAz < lBz && tBz < lBz && tCz < lBz))
      return 0;
   double e = 1E-11;
   //これを1E-14とすることで，干渉のチェックが行われるようになる場合があった
   // auto d = factorOfVectorToReachTriangle(tA, tB, tC, a, b);
   auto d = scalefactorToReach(line, triangle);
   if (d < 0. || d > 1.) return 0; /*面に到達できていない*/

   auto b_a = lB - lA;
   auto n = Normal(triangle);
   auto ps = lA + (b_a)*d;

   //ポリゴン頂点の最大最小でチェック
   if (((tAx - e > std::get<0>(ps) && tBx - e > std::get<0>(ps) && tCx - e > std::get<0>(ps)) ||
        (tAx + e < std::get<0>(ps) && tBx + e < std::get<0>(ps) && tCx + e < std::get<0>(ps))) ||
       ((tAy - e > std::get<1>(ps) && tBy - e > std::get<1>(ps) && tCy - e > std::get<1>(ps)) ||
        (tAy + e < std::get<1>(ps) && tBy + e < std::get<1>(ps) && tCy + e < std::get<1>(ps))) ||
       ((tAz - e > std::get<2>(ps) && tBz - e > std::get<2>(ps) && tCz - e > std::get<2>(ps)) ||
        (tAz + e < std::get<2>(ps) && tBz + e < std::get<2>(ps) && tCz + e < std::get<2>(ps))))
      return 1; /*面の最大最小範囲にすら入れていない*/

   auto ps_p0 = tA - ps;
   auto ps_p1 = tB - ps;
   auto ps_p2 = tC - ps;

   ps_p0 = ps_p0 / Norm(ps_p0);
   ps_p1 = ps_p1 / Norm(ps_p1);
   ps_p2 = ps_p2 / Norm(ps_p2);

   if (Dot(Cross(ps_p0, ps_p1), n) >= 0. && Dot(Cross(ps_p1, ps_p2), n) >= 0. &&
       Dot(Cross(ps_p2, ps_p0), n) >= 0.)
      return 3; /*a,bは面と交差*/
   else
      return 2; /*a,bは面と交差していないが，かなり惜しい*/
};
/* ------------------------------------------------------ */
//@ 全接触情報を返す
struct intersection {
   Tddd X;  //接触した物体の最も近い座標
   double distance;
   bool isIntersecting;
   double eps = 1E-10;
   int index_intersection_type;
   /*
   0: not intersecting
   1: intersecting with a sphere
   2: intersecting with a line
   3: intersecting with a triangle
   */
   /* ------------------------------------------------------ */
   intersection(const geometry::Sphere &sphere0, const geometry::Sphere &sphere1)
       : X({0, 0, 0}), distance(1E+40), isIntersecting(false) {
      if (IntersectQ(CoordinateBounds(sphere0), CoordinateBounds(sphere1))) {
         double R = Norm(std::get<1>(sphere0.X) - std::get<0>(sphere1.X));
         double overlap = (sphere0.radius + sphere1.radius - R) / 2.;
         this->X = sphere0.X + (sphere0.radius - overlap / 2.) * Normalize(sphere1.X - sphere0.X);
         this->distance = Norm(this->X - sphere0.X);
         this->isIntersecting = (overlap >= 0);
         this->index_intersection_type = 1;
         return;
      } else {
         this->distance = 1E+40;
         this->isIntersecting = false;
         this->index_intersection_type = 0;
         return;
      }
   };
   /* ------------------------------------------------------ */
   intersection(const geometry::Sphere &sphere, const geometry::Line &line) : X({0, 0, 0}), distance(1E+40), isIntersecting(false) {
      if (IntersectQ(CoordinateBounds(sphere), CoordinateBounds(line))) {
         Tddd X0X1 = std::get<1>(line.X) - std::get<0>(line.X);

         //
         // Dot(C - sphere.X, X0X1) = 0        (1)
         // C = std::get<0>(line.X) + X0X1 * t (2)
         //
         // Dot(std::get<0>(line.X) + X0X1 * t - sphere.X, X0X1) = 0
         // -> Dot(std::get<0>(line.X) - sphere.X, X0X1) + Dot(X0X1 * t , X0X1)
         // = 0
         // -> Dot(X0X1 , X0X1) * t = - Dot(std::get<0>(line.X) - sphere.X,
         // X0X1)
         // -> t = - Dot(std::get<0>(line.X) - sphere.X, X0X1)/Dot(X0X1 , X0X1)
         //
         double t = -Dot(std::get<0>(line.X) - sphere.X, X0X1) / Dot(X0X1, X0X1);
         this->X = std::get<0>(line.X) + X0X1 * t;
         this->distance = Norm(this->X - sphere.X);
         if (this->distance <= sphere.radius /*may hit*/) {
            if (-eps <= t && t <= 1 + eps) {
               //干渉する最も近い点は，線上にある
               this->isIntersecting = true;
               this->index_intersection_type = 2;
               return;
            } else {
               //干渉する最も近い点は，端点である可能性をチェック
               intersection intx0(sphere, geometry::Sphere(std::get<0>(line.X)));
               intersection intx1(sphere, geometry::Sphere(std::get<1>(line.X)));

               if (intx0.isIntersecting && intx1.distance > intx0.distance) {
                  this->isIntersecting = true;
                  this->distance = intx0.distance;
                  this->X = intx0.X;
                  this->index_intersection_type = 1;
                  return;
               } else if (intx1.isIntersecting && intx0.distance > intx1.distance) {
                  this->isIntersecting = true;
                  this->distance = intx1.distance;
                  this->X = intx1.X;
                  this->index_intersection_type = 1;
                  return;
               } else {
                  this->distance = 1E+40;
                  this->isIntersecting = false;
                  this->index_intersection_type = 0;
                  return;
               }
            }
         } else {
            this->distance = 1E+40;
            this->isIntersecting = false;
            this->index_intersection_type = 0;
            return;
         }
      } else {
         this->distance = 1E+40;
         this->isIntersecting = false;
         this->index_intersection_type = 0;
         return;
      }
   };
   /* ------------------------------------------------------ */
   intersection(const geometry::Sphere &sphere, const geometry::Triangle &triangle)
       : X({0, 0, 0}), distance(1E+40), isIntersecting(false) {
      // エラーの原因は初期値を返しているのかもしれない
      // まずは，coordinateboundsをチェックする．
      if (IntersectQ(CoordinateBounds(sphere), CoordinateBounds(triangle))) {
         auto [x0, y0, z0] = std::get<0>(triangle.X) - sphere.X;
         auto [x1, y1, z1] = std::get<1>(triangle.X) - sphere.X;
         auto [x2, y2, z2] = std::get<2>(triangle.X) - sphere.X;
         Tddd cross = {y2 * (z0 - z1) + y0 * (z1 - z2) + y1 * (-z0 + z2), x2 * (-z0 + z1) + x1 * (z0 - z2) + x0 * (-z1 + z2), x2 * (y0 - y1) + x0 * (y1 - y2) + x1 * (-y0 + y2)};
         double crossSquared = Dot(cross, cross);
         double t0 = ((y2 * z1 - y1 * z2) * (y1 * z0 - y2 * z0 - y0 * z1 + y2 * z1 + y0 * z2 - y1 * z2) + x1 * (x0 * (y1 - y2) * y2 + x0 * (z1 - z2) * z2 + x2 * (-2 * y1 * y2 + y0 * (y1 + y2) - 2 * z1 * z2 + z0 * (z1 + z2))) + pow(x2, 2) * (-(y0 * y1) - z0 * z1 + pow(y1, 2) + pow(z1, 2)) - x0 * x2 * (-(y1 * y2) - z1 * z2 + pow(y1, 2) + pow(z1, 2)) + pow(x1, 2) * (-(y0 * y2) - z0 * z2 + pow(y2, 2) + pow(z2, 2))) / crossSquared;
         double t1 = ((y0 - y2) * y2 * z0 * z1 + (y0 * y1 - 2 * y0 * y2 + y1 * y2) * z0 * z2 + x0 * x2 * (y0 * y1 - 2 * y0 * y2 + y1 * y2 + z0 * z1 - 2 * z0 * z2 + z1 * z2) - z1 * z2 * (-(y0 * y2) + pow(x0, 2) + pow(y0, 2)) + x1 * (x0 * (y0 - y2) * y2 + x0 * (z0 - z2) * z2 - x2 * (-(y0 * y2) + z0 * (z0 - z2) + pow(y0, 2))) + y2 * (-y1 + y2) * (pow(x0, 2) + pow(z0, 2)) +
                      pow(x2, 2) * (-(y0 * y1) - z0 * z1 + pow(y0, 2) + pow(z0, 2)) + (-(y0 * y1) + pow(x0, 2) + pow(y0, 2)) * pow(z2, 2)) /
                     crossSquared;
         double scale = (-(x2 * y1 * z0) + x1 * y2 * z0 + x2 * y0 * z1 - x0 * y2 * z1 - x1 * y0 * z2 + x0 * y1 * z2) / crossSquared;
         //結果
         this->distance = Norm(scale * cross);
         this->X = scale * cross + sphere.X;
         if (this->distance <= sphere.radius /*may hit*/) {
            if ((-eps <= t0 && t0 <= 1. + eps) && (-eps <= t1 && t1 <= 1. + eps) && (t0 + t1 <= 1. + eps)) {  //干渉する最も近い点は，三角形の面内にある．
               this->isIntersecting = true;
               this->index_intersection_type = 3;
               return;
            } else {
               intersection intx0(sphere, geometry::Line(std::get<0>(triangle.X), std::get<1>(triangle.X)));
               intersection intx1(sphere, geometry::Line(std::get<1>(triangle.X), std::get<2>(triangle.X)));
               intersection intx2(sphere, geometry::Line(std::get<2>(triangle.X), std::get<0>(triangle.X)));
               if (intx1.distance >= intx0.distance && intx2.distance >= intx0.distance && intx0.isIntersecting) {
                  this->isIntersecting = true;
                  this->distance = intx0.distance;
                  this->X = intx0.X;
                  this->index_intersection_type = intx0.index_intersection_type;
                  return;
               } else if (intx0.distance >= intx1.distance && intx2.distance >= intx1.distance && intx1.isIntersecting) {
                  this->isIntersecting = true;
                  this->distance = intx1.distance;
                  this->X = intx1.X;
                  this->index_intersection_type = intx1.index_intersection_type;
                  return;
               } else if (intx1.distance >= intx2.distance && intx0.distance >= intx2.distance && intx2.isIntersecting) {
                  this->isIntersecting = true;
                  this->distance = intx2.distance;
                  this->X = intx2.X;
                  this->index_intersection_type = intx2.index_intersection_type;
                  return;
               } else {
                  this->distance = 1E+40;
                  this->isIntersecting = false;
                  this->index_intersection_type = 0;
               }
            }
         } else {
            this->distance = 1E+40;
            this->isIntersecting = false;
            this->index_intersection_type = 0;
            return;
         }
      } else {
         this->distance = 1E+40;
         this->isIntersecting = false;
         this->index_intersection_type = 0;
         return;
      }
   };
};
/* ------------------------------------------------------ */
/* ------------------------------------------------------ */
/*namespace_geometry_detail
  namespace_geometry_detail*/
/*namespace_geometry_code*/
using V_d = std::vector<double>;
using VV_d = std::vector<std::vector<double>>;
using V_i = std::vector<int>;
using VV_i = std::vector<std::vector<int>>;
/* ------------------------------------------------------ */
double normalDirDistanceFromTriangle(const V_d &p0, const V_d &p1, const V_d &p2, const V_d &a) {
   return Dot(TriangleNormal(p0, p1, p2), p0 - a);
};
double normalDirDistanceFromTriangle(const T3Tddd &ps, const Tddd &&a) {
   return Dot(TriangleNormal(std::get<0>(ps), std::get<1>(ps), std::get<2>(ps)),
              std::get<0>(ps) - a);
};
V_d vectorToTriangle(const V_d &p0, const V_d &p1, const V_d &p2, const V_d &a) {
   auto n = TriangleNormal(p0, p1, p2);
   return n * Dot(n, p0 - a);
};
double factorOfVectorToReachTriangle(const Tddd &p0, const Tddd &p1, const Tddd &p2, const Tddd &a, const Tddd &b) {
   //オーダーが匹敵する物を選ぶ
   double log_b_a = log10(Norm(b - a));
   double diff0 = std::abs(log10(Norm(p0 - a) - log_b_a));
   double diff1 = std::abs(log10(Norm(p1 - a) - log_b_a));
   double diff2 = std::abs(log10(Norm(p2 - a) - log_b_a));
   auto n = TriangleNormal(p0, p1, p2);
   if (diff0 < diff1 && diff0 < diff2)
      return Dot(p0 - a, n) / Dot(b - a, n);
   else if (diff1 < diff0 && diff1 < diff2)
      return Dot(p1 - a, n) / Dot(b - a, n);
   else
      return Dot(p2 - a, n) / Dot(b - a, n);
};
/* ------------------------------------------------------ */
int isIntersectingSurface(const Tddd &p0, const Tddd &p1, const Tddd &p2, const Tddd &a, const Tddd &b) {
   /* 0:頂点の最大最小の範囲の外で，片方にa,bgがある */
   /* 1:拡大した面には入れているが，多角形の頂点の最大最小範囲にすら入れていない
    */
   /* 2:a,bは多角形の面と交差していないが，かなり惜しい */
   /* 3:a,bは多角形の面と交差 */
   // double e = 1E-11;
   // //これを1E-14とすることで，干渉のチェックが行われるようになる場合があった．
   // if ((p0[0] - e > a[0] && p1[0] - e > a[0] && p2[0] - e > a[0] && p0[0] - e
   // > b[0] && p1[0] - e > b[0] && p2[0] - e > b[0]) || 	(p0[0] + e < a[0] &&
   // p1[0] + e < a[0] && p2[0] + e < a[0] && p0[0] + e < b[0] && p1[0] + e <
   // b[0] && p2[0] + e < b[0]) || 	(p0[1] - e > a[1] && p1[1] - e > a[1] &&
   // p2[1] - e > a[1] && p0[1] - e > b[1] && p1[1] - e > b[1] && p2[1] - e >
   // b[1]) || 	(p0[1] + e < a[1] && p1[1] + e < a[1] && p2[1] + e < a[1] &&
   // p0[1] + e < b[1] && p1[1] + e < b[1] && p2[1] + e < b[1]) || 	(p0[2] - e >
   // a[2] && p1[2] - e > a[2] && p2[2] - e > a[2] && p0[2] - e > b[2] && p1[2]
   // - e > b[2] && p2[2] - e > b[2]) || 	(p0[2] + e < a[2] && p1[2] + e < a[2]
   // && p2[2] + e < a[2] && p0[2] + e < b[2] && p1[2] + e < b[2] && p2[2] + e <
   // b[2])) 	return 0;

   if ((std::get<0>(p0) > std::get<0>(a) && std::get<0>(p1) > std::get<0>(a) &&
        std::get<0>(p2) > std::get<0>(a) && std::get<0>(p0) > std::get<0>(b) &&
        std::get<0>(p1) > std::get<0>(b) && std::get<0>(p2) > std::get<0>(b)) ||
       (std::get<0>(p0) < std::get<0>(a) && std::get<0>(p1) < std::get<0>(a) &&
        std::get<0>(p2) < std::get<0>(a) && std::get<0>(p0) < std::get<0>(b) &&
        std::get<0>(p1) < std::get<0>(b) && std::get<0>(p2) < std::get<0>(b)) ||
       (std::get<1>(p0) > std::get<1>(a) && std::get<1>(p1) > std::get<1>(a) &&
        std::get<1>(p2) > std::get<1>(a) && std::get<1>(p0) > std::get<1>(b) &&
        std::get<1>(p1) > std::get<1>(b) && std::get<1>(p2) > std::get<1>(b)) ||
       (std::get<1>(p0) < std::get<1>(a) && std::get<1>(p1) < std::get<1>(a) &&
        std::get<1>(p2) < std::get<1>(a) && std::get<1>(p0) < std::get<1>(b) &&
        std::get<1>(p1) < std::get<1>(b) && std::get<1>(p2) < std::get<1>(b)) ||
       (std::get<2>(p0) > std::get<2>(a) && std::get<2>(p1) > std::get<2>(a) &&
        std::get<2>(p2) > std::get<2>(a) && std::get<2>(p0) > std::get<2>(b) &&
        std::get<2>(p1) > std::get<2>(b) && std::get<2>(p2) > std::get<2>(b)) ||
       (std::get<2>(p0) < std::get<2>(a) && std::get<2>(p1) < std::get<2>(a) &&
        std::get<2>(p2) < std::get<2>(a) && std::get<2>(p0) < std::get<2>(b) &&
        std::get<2>(p1) < std::get<2>(b) && std::get<2>(p2) < std::get<2>(b)))
      return 0;
   double e =
       1E-11;  //これを1E-14とすることで，干渉のチェックが行われるようになる場合があった．

   auto d = factorOfVectorToReachTriangle(p0, p1, p2, a, b);
   if (d < 0. || d > 1.) return 0; /*面に到達できていない*/

   auto b_a = b - a;
   auto n = TriangleNormal(p0, p1, p2);
   auto ps = a + (b_a)*d;

   //ポリゴン頂点の最大最小でチェック
   if (((std::get<0>(p0) - e > std::get<0>(ps) &&
         std::get<0>(p1) - e > std::get<0>(ps) &&
         std::get<0>(p2) - e > std::get<0>(ps)) ||
        (std::get<0>(p0) + e < std::get<0>(ps) &&
         std::get<0>(p1) + e < std::get<0>(ps) &&
         std::get<0>(p2) + e < std::get<0>(ps))) ||
       ((std::get<1>(p0) - e > std::get<1>(ps) &&
         std::get<1>(p1) - e > std::get<1>(ps) &&
         std::get<1>(p2) - e > std::get<1>(ps)) ||
        (std::get<1>(p0) + e < std::get<1>(ps) &&
         std::get<1>(p1) + e < std::get<1>(ps) &&
         std::get<1>(p2) + e < std::get<1>(ps))) ||
       ((std::get<2>(p0) - e > std::get<2>(ps) &&
         std::get<2>(p1) - e > std::get<2>(ps) &&
         std::get<2>(p2) - e > std::get<2>(ps)) ||
        (std::get<2>(p0) + e < std::get<2>(ps) &&
         std::get<2>(p1) + e < std::get<2>(ps) &&
         std::get<2>(p2) + e < std::get<2>(ps))))
      return 1; /*面の最大最小範囲にすら入れていない*/

   auto ps_p0 = p0 - ps;
   auto ps_p1 = p1 - ps;
   auto ps_p2 = p2 - ps;

   ps_p0 = ps_p0 / Norm(ps_p0);
   ps_p1 = ps_p1 / Norm(ps_p1);
   ps_p2 = ps_p2 / Norm(ps_p2);

   if (Dot(Cross(ps_p0, ps_p1), n) >= 0. && Dot(Cross(ps_p1, ps_p2), n) >= 0. &&
       Dot(Cross(ps_p2, ps_p0), n) >= 0.)
      return 3; /*a,bは面と交差*/
   else
      return 2; /*a,bは面と交差していないが，かなり惜しい*/
};
int isIntersectingSurface(const Tddd &p0, const Tddd &p1, const Tddd &p2, const T2Tddd &ab) {
   return isIntersectingSurface(p0, p1, p2, std::get<0>(ab), std::get<1>(ab));
}
/* ------------------------------------------------------ */
Tddd t0t1l(const Tddd &aIN, const Tddd &bIN, const Tddd &cIN, const Tddd &AIN, const Tddd &BIN) {
   const auto [Ax, Ay, Az] = AIN - cIN;
   const auto [Bx, By, Bz] = BIN - cIN;
   const auto [ax, ay, az] = aIN - cIN;
   const auto [bx, by, bz] = bIN - cIN;
   return Tddd{Az * Bx * by - Az * bx * By - Ay * Bx * bz + Ax * By * bz + Ay * bx * Bz - Ax * by * Bz,
               Ay * az * Bx - ay * Az * Bx - Ax * az * By + ax * Az * By + Ax * ay * Bz - ax * Ay * Bz,
               az * Bx * by - az * bx * By - ay * Bx * bz + ax * By * bz + ay * bx * Bz - ax * by * Bz} /
          (-(Ax * az * by) + ax * Az * by + az * Bx * by - az * bx * By + ax * By * bz + Ay * (az * bx - ax * bz) - ax * by * Bz + ay * (-(Az * bx) + Ax * bz - Bx * bz + bx * Bz));
};
bool intersectLineTriangle(const Tddd &a, const Tddd &b, const Tddd &c, const Tddd &A, const Tddd &B) {
   if (!IntersectQ(CoordinateBounds(T2Tddd{A, B}), CoordinateBounds(T3Tddd{a, b, c})))
      return false;
   else {
      const Tdd range = {-1E-5, 1. + 1E-5};
      const auto [t0, t1, l] = t0t1l(a, b, c, A, B);
      if (Between(t0, range) && Between(t1, range) && Between(t0 + t1, range) &&
          Between(l, range))
         return true;
      else
         return false;
   }
};
bool intersectLineTriangle(const T3Tddd &abc, const Tddd &A, const Tddd &B) {
   return intersectLineTriangle(std::get<0>(abc), std::get<1>(abc),
                                std::get<2>(abc), A, B);
};
/* ------------------------------------------------------ */
bool intersectLineSquare(const Tddd &a, const Tddd &b, const Tddd &c, const Tddd &A, const Tddd &B) {
   const Tdd range = {-1E-5, 1. + 1E-5};
   const auto [t0, t1, l] = t0t1l(a, b, c, A, B);
   if (Between(t0, range) && Between(t1, range) && Between(l, range))
      return true;
   else
      return false;
};
bool intersectLineSquare(const T3Tddd &abc, const Tddd &A, const Tddd &B) {
   return intersectLineSquare(std::get<0>(abc), std::get<1>(abc),
                              std::get<2>(abc), A, B);
};
/* ------------------------------------------------------ */
bool intersectLineCube(const T3Tdd &boundsTuple, const Tddd &A, const Tddd &B) {
   const CoordinateBounds bounds(boundsTuple);
   if (bounds.isInside(A) || bounds.isInside(B))
      return true;
   else {
      const auto [X0, X1] = std::get<0>(bounds.bounds);
      const auto [Y0, Y1] = std::get<1>(bounds.bounds);
      const auto [Z0, Z1] = std::get<2>(bounds.bounds);
      if (intersectLineSquare({X0, Y0, Z0}, {X1, Y1, Z0}, {X0, Y1, Z0}, A, B) ||
          intersectLineSquare({X0, Y0, Z1}, {X1, Y1, Z1}, {X0, Y1, Z1}, A, B) ||
          intersectLineSquare({X0, Y0, Z0}, {X1, Y0, Z1}, {X0, Y0, Z1}, A, B) ||
          intersectLineSquare({X0, Y1, Z0}, {X1, Y1, Z1}, {X0, Y1, Z1}, A, B) ||
          intersectLineSquare({X0, Y0, Z0}, {X0, Y1, Z1}, {X0, Y0, Z1}, A, B) ||
          intersectLineSquare({X1, Y0, Z0}, {X1, Y1, Z1}, {X1, Y0, Z1}, A, B))
         return true;
      else
         return false;
   }
};
bool intersectLineCube(const T3Tdd &boundsTuple, const T2Tddd &AB) {
   return intersectLineCube(boundsTuple, std::get<0>(AB), std::get<1>(AB));
};
/* ------------------------------------------------------ */
bool intersectTriangleCube(const T3Tdd &boundsTuple, const T3Tddd &verticesIN) {
   const auto [A, B, C] = verticesIN;
   if (isInside(A, boundsTuple) || isInside(B, boundsTuple) ||
       isInside(C, boundsTuple) || intersectLineCube(boundsTuple, A, B) ||
       intersectLineCube(boundsTuple, B, C) ||
       intersectLineCube(boundsTuple, C, A))
      return true;
   else {

      const auto [X0, X1] = std::get<0>(boundsTuple);
      const auto [Y0, Y1] = std::get<1>(boundsTuple);
      const auto [Z0, Z1] = std::get<2>(boundsTuple);
      if (intersectLineTriangle(verticesIN, {X0, Y0, Z0}, {X1, Y0, Z0}) ||
          intersectLineTriangle(verticesIN, {X0, Y0, Z1}, {X1, Y0, Z1}) ||
          intersectLineTriangle(verticesIN, {X0, Y1, Z0}, {X1, Y1, Z0}) ||
          intersectLineTriangle(verticesIN, {X0, Y1, Z1}, {X1, Y1, Z1}) ||
          intersectLineTriangle(verticesIN, {X0, Y0, Z0}, {X0, Y1, Z0}) ||
          intersectLineTriangle(verticesIN, {X0, Y0, Z1}, {X0, Y1, Z1}) ||
          intersectLineTriangle(verticesIN, {X1, Y0, Z0}, {X1, Y1, Z0}) ||
          intersectLineTriangle(verticesIN, {X1, Y0, Z1}, {X1, Y1, Z1}) ||
          intersectLineTriangle(verticesIN, {X0, Y0, Z0}, {X0, Y0, Z1}) ||
          intersectLineTriangle(verticesIN, {X0, Y1, Z0}, {X0, Y1, Z1}) ||
          intersectLineTriangle(verticesIN, {X1, Y0, Z0}, {X1, Y0, Z1}) ||
          intersectLineTriangle(verticesIN, {X1, Y1, Z0}, {X1, Y1, Z1}))
         return true;
      else
         return false;
   }
};
/* ------------------------------------------------------ */
bool intersectingQ(const T3Tdd &boundsTuple, const T3Tddd &verticesIN) {
   return intersectTriangleCube(boundsTuple, verticesIN);
};
bool intersectingQ(const T3Tdd &boundsTuple, const Tddd &p) {
   const CoordinateBounds bounds(boundsTuple);
   return bounds.isInside(p);
};
/* ------------------------------------------------------ */
Tddd vectorToInfiniteLine(const Tddd &P, Tddd A, Tddd B) {
   // Tddd BA = B - A, AP = A - P;
   // return AP - BA * Dot(AP, BA) / Dot(BA, BA);
   /* ------------------------------------------------------ */
   B -= A;
   A -= P;
   return A - B * Dot(A, B) / Dot(B, B);
};
double distanceToInfiniteLine(const Tddd &P, const Tddd &A, const Tddd &B) {
   return Norm(vectorToInfiniteLine(P, A, B));
};
/* ------------------------------------------------------ */
bool IntersectSpheres(const Tddd &X0, double const r0, const Tddd &X1, double r1) {
   return (Norm(X0 - X1) <= r0 + r1);
};
/* ------------------------------------------------------ */
bool IntersectSphereCube(const Tddd &P, const double &r, const T3Tdd &bounds) {
   auto [minX, maxX] = Transpose(bounds);
   Tddd Xc = (minX + maxX) * 0.5;
   Tddd difMinMax = maxX - Xc, PtoC = Xc - P;
   return (IntersectSpheres(P, r, Xc, Norm(difMinMax)) &&
           (std::abs(std::get<0>(PtoC)) <= r + std::get<0>(difMinMax)) &&
           (std::abs(std::get<1>(PtoC)) <= r + std::get<1>(difMinMax)) &&
           (std::abs(std::get<2>(PtoC)) <= r + std::get<2>(difMinMax)));
};
/* ------------------------------------------------------ */
bool IntersectSphereCube(const geometry::Sphere &s, const T3Tdd &bounds) {
   return IntersectSphereCube(s.X, s.radius, bounds);
};
/* ------------------------------------------------------ */
bool IntersectSphereCubeNUMBER(const Tddd &P, const double &r, const T3Tdd &bounds) {
   auto [minX, maxX] = Transpose(bounds);
   Tddd Xc = (minX + maxX) * 0.5;
   Tddd difMinMax = maxX - Xc, PtoC = Xc - P;
   return (IntersectSpheres(P, r, Xc, Norm(difMinMax)) &&
           (std::abs(std::get<0>(PtoC)) <= r + std::get<0>(difMinMax)) &&
           (std::abs(std::get<1>(PtoC)) <= r + std::get<1>(difMinMax)) &&
           (std::abs(std::get<2>(PtoC)) <= r + std::get<2>(difMinMax)));
};
bool IntersectSphereCubeNUMBER(const geometry::Sphere &s, const T3Tdd &bounds) {
   return IntersectSphereCube(s.X, s.radius, bounds);
};
/* ------------------------------------------------------ */
// V_d pOnSurface(const V_d &p0, const V_d &p1, const V_d &p2, const V_d &a,
// const V_d &b)
// {
//   V_d n = TriangleNormal(p0, p1, p2);
//   return a + (b - a) * Dot(/*tangential vector*/ Norm(p0 - a) < Norm(p1 - a)
//   ? (p1 - a) : (p0 - a),
//                            /*normal vector*/ n) /
//                  Dot((b - a), n);
//   //return a + (b-a)*Dot(p0-a,n)/Dot(b-a,n);
// }

V_d pOnSurface(const V_d &p0, const V_d &p1, const V_d &p2, const V_d &a, V_d &b) {
   V_d n = TriangleNormal(p0, p1, p2), b_a = b - a;
   return a + b_a * Dot(p0 - a, n) / Dot(b_a, n);  //分母が0の場合はあり得る
}

Tddd pOnSurfaceTuple(const Tddd &p0, const Tddd &p1, const Tddd &p2, const Tddd &a, const Tddd &b) {
   Tddd n = TriangleNormal(p0, p1, p2), b_a = b - a;
   return a + b_a * Dot(p0 - a, n) / Dot(b_a, n);  //分母が0の場合はあり得る
}

// V_d pOnSurface(const VV_d &p0p1p2, const VV_d &ab) {
//    if (p0p1p2.size() != 3)
//       throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "point size = " + std::to_string(p0p1p2.size()));
//    return pOnSurface(p0p1p2[0], p0p1p2[1], p0p1p2[2], ab[0], ab[1]);
// }

Tddd pOnSurfaceTuple(const T3Tddd &p0p1p2, const T2Tddd &ab) {
   return pOnSurfaceTuple(std::get<0>(p0p1p2), std::get<1>(p0p1p2), std::get<2>(p0p1p2), std::get<0>(ab), std::get<1>(ab));
}

int isIntersectingSurface(const T3Tddd &p0p1p2, const T2Tddd &ab) {
   // if (p0p1p2.size() != 3)
   // 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "point size
   // = " + std::to_string(p0p1p2.size()));

   return isIntersectingSurface(std::get<0>(p0p1p2), std::get<1>(p0p1p2),
                                std::get<2>(p0p1p2), std::get<0>(ab),
                                std::get<1>(ab));
};
/* ------------------------------------------------------ */
// class intersectionTriangleLine
// {
// public:
// 	V_d X;
// 	V_d normal;
// 	V_d vecA2X;
// 	V_d vecX2B;
// 	V_d vecX2B_;
// 	bool isIntersect;
// 	int indexOfTriangle;
// 	// Vnewの方向を変更
// 	// Xonの位置を計算
// 	intersectionTriangleLine(const V_d &p0,
// 							 const V_d &p1,
// 							 const V_d &p2,
// 							 const V_d &A,
// 							 const V_d &B)
// 		: X({}), normal({}), vecA2X({}), vecX2B({}), vecX2B_({}),
// isIntersect(false), indexOfTriangle(0)
// 	{
// 		if (isIntersectingSurface(p0, p1, p2, A, B) == 3
// /*三角形と干渉した場合*/)
// 		{
// 			this->normal = TriangleNormal(p0, p1, p2);
// 			this->X = pOnSurface(p0, p1, p2, A, B);
// 			this->vecA2X = this->X - A;
// 			this->vecX2B = B - this->X;
// 			this->vecX2B_ = this->reflect(this->vecX2B,
// this->normal); 			this->isIntersect = true;
// 		}
// 	};
// 	////////////
// 	intersectionTriangleLine(const VVV_d &p0p1p2s,
// 							 const V_d &A,
// 							 const V_d &B,
// 							 const V_i &exceptIndices =
// {}) 		: X({}), normal({}), vecA2X({}), vecX2B({}), vecX2B_({}),
// isIntersect(false), indexOfTriangle(0)
// 	{
// 		double closest_distFromWall = 1E+100;
// 		double normal_distA2X;
// 		for (auto i = 0; i < p0p1p2s.size(); i++)
// 		{
// 			if (!MemberQ(exceptIndices, i))
// 			{
// 				intersectionTriangleLine LT(p0p1p2s[i][0], p0p1p2s[i][1],
// p0p1p2s[i][2], A, B); 				normal_distA2X = Norm(Dot(LT.vecA2X, LT.normal)); 				if
// (LT.isIntersect && isFinite(normal_distA2X) && normal_distA2X <
// closest_distFromWall)
// 				{
// 					closest_distFromWall = normal_distA2X;
// 					this->X = LT.X;
// 					this->normal = LT.normal;
// 					this->vecA2X = LT.vecA2X;
// 					this->vecX2B = LT.vecX2B;
// 					this->vecX2B_ = LT.vecX2B_;
// 					this->isIntersect = LT.isIntersect;
// 					this->indexOfTriangle = i;
// 				}
// 			}
// 		};
// 	};
// 	//////////////
// 	V_d reflectIfPossible(const V_d &v) const
// 	{
// 		if (this->isIntersect)
// 			return this->reflect(v, this->normal);
// 		else
// 			return v;
// 	};
// 	V_d reflect(const V_d &v, const V_d &n) const
// 	{
// 		return v - 2. * Dot(v, n) * n;
// 	};
// };
/* ------------------------------------------------------ */
namespace geometry {
class Point_Line {
  public:
   double t;  // parameter of v from p_line0
   V_d x;     // coordinate of point
   V_d v;     // vector of line
   double d2line;
   double d2line_segment;
   Point_Line(const V_d &p, const V_d &p_line0, const V_d &p_line1) {
      this->v = p_line1 - p_line0;
      this->t = -Dot(p_line0 - p, this->v) / Dot(this->v, this->v);
      this->x = p_line0 + this->t * this->v;
      this->d2line = Norm(this->x - p);

      if (this->t < 0)
         this->d2line_segment = Norm(p_line0 - p);
      else if (this->t > 1)
         this->d2line_segment = Norm(p_line1 - p);
      else
         this->d2line_segment = this->d2line;
   };
};

//基本的な操作なので，fundamentalに持ってきた
// どうして，> だとうまくいくのか？
// まがっているのか？
// bool isConvexPolygon(const VV_d &ps, const V_d &normal)
// {
//   auto s = ps.size();
//   if (s < 3)
//     return false;

//   V_d v0, v1;
//   for (auto i = 0; i < s; i++)
//   {
//     v0 = ps[i] - ps[(s + i - 1) % s];
//     v1 = ps[(s + i + 1) % s] - ps[i];
//     if (MyVectorAngle(v0 /*基準*/, v1, normal) < 0. /*ccw*/)
//       return false;
//   }
//   return true;
// };

bool isConvexPolygon(const std::vector<Tddd> &ps, const Tddd &n) {
   auto s = ps.size();
   if (s < 3)
      return false;
   else if (s == 3)
      return true;

   for (auto i = 0; i < s; ++i) {
      auto v0 = ps[i + 1] - ps[i];
      auto v1 = ps[i + 2] - ps[i + 1];
      auto angle = Dot(Cross(v0, v1), n);
      if (angle <= 1E-13) return false;  //符号が変わったらfalse
   }
   return true;
};

bool isConcavePolygon(const std::vector<Tddd> &ps, const Tddd &n) {
   return !isConvexPolygon(ps, n);
};

bool isConvexPolygon(const std::vector<Tdd> &ps) {
   std::vector<Tddd> Ps(ps.size());
   int i = 0;
   for (const auto &v : ps) Ps[i++] = {std::get<0>(v), std::get<1>(v), 0.};
   return isConvexPolygon(Ps, Tddd{0., 0., 1});
};

bool isConcavePolygon(const std::vector<Tdd> &ps) {
   return !isConvexPolygon(ps);
};

// bool isConvexPolygon(const VV_d &ps)
// {
// 	auto s = ps.size();
// 	if (s < 3)
// 		return false;

// 	V_d v0 = *ps.begin() - *ps.rbegin();
// 	V_d v1 = *(ps.begin() + 1) - *ps.begin();
// 	V_d normal = Cross(v0, v1);
// 	normal = normal / Norm(normal);
// 	double angle = MyVectorAngle(v0 /*基準*/, v1, normal);
// 	if (!isFinite(normal))
// 		return false;

// 	for (auto i = 0; i < s; i++)
// 	{
// 		// 0->1 1->2
// 		angle = angle * MyVectorAngle(ps[i] - ps[(s + i - 1) % s]
// /*基準*/, 									  ps[(s + i + 1) % s] - ps[i], normal); 		if (!isFinite(angle)) 			return
// false; 		if (angle <= 1E-13) 			return false; //符号が変わったらfalse
// 	}
// 	return true;
// };
/*ccw angle
 *       *
 *     / | \
 *    *  |  *
 *   /  \|/  \
 *  *----*----*
 *   \2 1|1 3/
 *    \  |  /
 *     \3|2/
 *       *
 */
// angle_sets = {{a1,a2,a3},{b1,b2,b3},...}
// bool isInConvexPolygon(const VV_d &angle_sets)
// {
//   auto s = angle_sets.size();
//   for (auto i = 0; i < s + 1; i++)
//     if (angle_sets[i % s][1] + angle_sets[(i + 1) % s][2] > M_PI)
//       return false;
//   return true;
// };
/////////////////
// bool isConcavePolygon(const VV_d &ps, const V_d &normal) { return
// !geometry::isConvexPolygon(ps, normal); };
bool isConcavePolygon(const VV_d &ps) {
   auto s = ps.size();
   if (s < 3) return false;

   V_d v0 = *ps.begin() - *ps.rbegin();
   V_d v1 = *(ps.begin() + 1) - *ps.begin();
   V_d normal = Cross(v0, v1);
   normal = normal / Norm(normal);
   double angle = MyVectorAngle(v0 /*基準*/, v1, normal);
   if (!isFinite(normal)) return false;

   for (auto i = 0; i < s; i++) {
      // 0->1 1->2
      angle = angle * MyVectorAngle(ps[i] - ps[(s + i - 1) % s] /*基準*/,
                                    ps[(s + i + 1) % s] - ps[i], normal);
      if (!isFinite(angle)) return false;
      if (angle < 0.) return true;  //符号が変わったらtrue
   }
   return false;  //符号が変わらなかったのでfalse
};
bool isConcavePolygon(const std::vector<Tddd> &ps) {
   auto s = ps.size();
   if (s < 3) return false;

   auto v0 = *ps.begin() - *ps.rbegin();
   auto v1 = *(ps.begin() + 1) - *ps.begin();
   auto normal = Cross(v0, v1);
   normal = normal / Norm(normal);
   double angle = MyVectorAngle(v0 /*基準*/, v1, normal);
   if (!isFinite(normal)) return false;

   for (auto i = 0; i < s; i++) {
      // 0->1 1->2
      angle = angle * MyVectorAngle(ps[i] - ps[(s + i - 1) % s] /*基準*/,
                                    ps[(s + i + 1) % s] - ps[i], normal);
      if (!isFinite(angle)) return false;
      if (angle < 0.) return true;  //符号が変わったらtrue
   }
   return false;  //符号が変わらなかったのでfalse
};
//--------------------------
class point {
  public:
   point(const V_d &xyz, const int i) : active(true), X(xyz), index(i){};
   int index;
   bool active;
   V_d X;
   double area;
   double angle;
};
using V_pp = std::vector<point *>;
using VV_pp = std::vector<std::vector<point *>>;
//--------------------------
VV_d extractX(const V_pp &ps) {
   VV_d ret(0);
   for (const auto &p : ps) ret.emplace_back(p->X);
   return ret;
};
V_d extractAreas(const V_pp &ps) {
   V_d ret(0);
   for (const auto &p : ps) ret.emplace_back(p->area);
   return ret;
};
V_d extractAngle(const V_pp &ps) {
   V_d ret(0);
   for (const auto &p : ps) ret.emplace_back(p->angle);
   return ret;
};
V_i extractIndices(const V_pp &ps) {
   V_i ret(0);
   for (const auto &p : ps) ret.emplace_back(p->index);
   return ret;
};
VV_i extractIndices(const VV_pp &pps) {
   VV_i ret(0);
   for (const auto &ps : pps) ret.emplace_back(extractIndices(ps));
   return ret;
};
//-------------------------
/*polygon_detail
多角形クラス
polygon_detail*/
/*polygon_code*/
class polygon {
  public:
   V_pp points;
   ~polygon() {
      for (const auto &p : this->points)
         if (p) delete p;
   };
   polygon(const VV_d &xyz_IN) {
      int s = xyz_IN.size();
      if (s < 3)
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__,
                             "this is not polygon: size < 3");

      this->points.resize(s);
      for (int i = 0; i < xyz_IN.size(); i++)
         this->points[i] = new point(xyz_IN[i], i);
   };

   void activateAllPoints() {
      for (const auto &p : points) p->active = true;
   };
   V_pp getActivePoints() {
      V_pp ret;
      for (const auto &p : this->points)
         if (p->active) ret.emplace_back(p);
      return ret;
   };
   V_pp getAllPoints() { return this->points; };

   void calculateArea(const V_pp &ps, const V_d &normal) {
      int s = ps.size();
      for (auto i = 0; i < s; i++)
         ps[i]->area = DirectedArea(ps[i]->X - ps[(s + i - 1) % s]->X,
                                    ps[(i + 1) % s]->X - ps[i]->X, normal);
   };

   V_d getExteriorAngles(const V_pp &ps /*ccw*/, const V_d &normal) {
      auto s = ps.size();
      V_d ret(s);
      for (auto i = 0; i < s; i++)
         ret[i] = MyVectorAngle(ps[i]->X - ps[(s + i - 1) % s]->X,
                                ps[(i + 1) % s]->X - ps[i]->X, normal);
      /*前後の線からなる多角形の外角*/
      return ret;
   };

   V_d getInteriorAngles(const V_pp &ps /*ccw*/, const V_d &normal) {
      return M_PI - getExteriorAngles(ps, normal);
   };

   void calculateAngle(const V_pp &ps) {
      auto s = ps.size();
      for (auto i = 0; i < s; i++)
         ps[i]->angle = MyVectorAngle(ps[i]->X - ps[(s + i - 1) % s]->X,
                                      ps[(i + 1) % s]->X - ps[i]->X);
      /*前後の線からなる多角形の外角*/
   };

   bool isSmallAngle(const point *p0, const point *p1, const point *p2,
                     const double smallangle) {
      if (std::abs(MyVectorAngle(
              p1->X - p0->X, p2->X - p1->X) /*前後の線からなる多角形の外角*/) <
          smallangle)
         return true;  // too small
      return false;
   };

   bool getPointsMeetCondition(const V_pp &ps, const double smallangle_IN,
                               point *&select_p, int &current_index) {
      bool found = false;
      int s = ps.size();
      //徐々にsmallangleの制限を弱くしていく
      for (int i = 0; i < s; i++)
         if (ps[i]->area > 0. && ps[i]->angle > smallangle_IN &&
             myIsfinite(ps[i]->area) && myIsfinite(ps[i]->angle))
            if ((select_p == nullptr) /*first time*/ ||
                ps[i]->area < select_p->area /*from 2nd time*/) {
               select_p = ps[i];
               current_index = i;
               found = true;
            }
      return found;
   };

   // bool getPointsMeetCondition(const V_pp &ps, const double smallangle_IN,
   // point *&select_p, int &index)
   // {
   //   bool found = false;
   //   int s = ps.size();
   //   double smallangle;
   //   int max_try = 1000; //試行回数
   //   for (int k = 0; k < max_try; k++)
   //   {
   //     //徐々にsmallangleの制限を弱くしていく
   //     for (int i = 0; i < s; i++)
   //     {
   //       smallangle = smallangle_IN * (1. - (double)k / (max_try - 1.));
   //       //smallangle_IN,smallangle_IN*0.999,,smallangle_IN*0.998 if
   //       (ps[i]->area > 0. && ps[i]->angle > smallangle)
   //       {
   //         if ((select_p == nullptr) /*first time*/ || ps[i]->area <
   //         select_p->area /*from 2nd time*/)
   //         {
   //           select_p = ps[i];
   //           index = i;
   //           found = true;
   //         }
   //       }
   //     }
   //     if (found)
   //       break;
   //   }
   //   return found;
   // };

   bool isConvexPolygon(const V_pp &ps, const V_d &normal) const {
      auto s = ps.size();
      for (auto i = 0; i < s; i++)
         if (MyVectorAngle(ps[i]->X - ps[(s + i - 1) % s]->X,
                           ps[(i + 1) % s]->X - ps[i]->X, normal) < 0.)
            return false;  //外角がマイナス：時計回り
      return true;
      // for (const auto &angle : this->getExteriorAngles(ps, normal))
      //   if (angle < 0)
      //     return false;
      // return true;
   };

   ///////////// 4 点の場合の特別な分割 ///////////////

   // V_i legalConnection4(const V_pp &ps /*4点*/, const V_d &normal) const
   // {
   //   //
   //   //    3 *----* 2
   //   //      |b   |
   //   //      |    |
   //   //      |    |
   //   //      |   a|
   //   //    0 *----* 1
   //   //

   //   auto A = TriangleAngles(ps[0]->X,ps[1]->X,ps[2]->X);
   //   auto B = TriangleAngles(ps[2]->X,ps[3]->X,ps[0]->X);
   //   auto AB = Join(A,B);

   //   auto C = TriangleAngles(ps[1]->X,ps[2]->X,ps[3]->X);
   //   auto D = TriangleAngles(ps[3]->X,ps[0]->X,ps[1]->X);
   //   auto CD = Join(C,D);

   //   if(Max(AB) > Max(CD))
   //     return {{0,1,2},{2,3,2}};
   //   else
   //     return {{1,2,3},{3,0,1}};

   //   double sum;
   //   if (myIsfinite(sum = MyVectorAngle(ps[2]->X - ps[1]->X, ps[0]->X -
   //   ps[1]->X) + MyVectorAngle(ps[0]->X - ps[3]->X, ps[2]->X - ps[3]->X)) ||
   //   M_PI <= sum)
   //     return {3, 1};
   //   else if (myIsfinite(sum = MyVectorAngle(ps[3]->X - ps[2]->X, ps[1]->X -
   //   ps[2]->X) + MyVectorAngle(ps[1]->X - ps[0]->X, ps[3]->X - ps[0]->X)) ||
   //   M_PI <= sum)
   //     return {0, 2};
   //   else
   //     throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   // };

   VV_pp legalDivision4(const V_pp &ps /*4点*/, const V_d &normal) const {
      // V_i ind = legalConnection4(ps, normal);
      // //与えられたpsのインデックスで，元のindexとは異なる
      try {
         // TriangleNormal(ps[ind[0]]->X, ps[ind[1]]->X, ps[(ind[1] + 1) %
         // 4]->X); TriangleNormal(ps[ind[1]]->X, ps[ind[0]]->X, ps[(ind[0] + 1)
         // % 4]->X);

         // return extractIndices(VV_pp{{ps[ind[0]], ps[ind[1]], ps[(ind[1] + 1)
         // % 4]},
         //                             {ps[ind[1]], ps[ind[0]], ps[(ind[0] + 1)
         //                             % 4]}});

         auto A = TriangleAngles(ps[0]->X, ps[1]->X, ps[2]->X);
         auto B = TriangleAngles(ps[2]->X, ps[3]->X, ps[0]->X);
         auto AB = Join(A, B);

         auto C = TriangleAngles(ps[1]->X, ps[2]->X, ps[3]->X);
         auto D = TriangleAngles(ps[3]->X, ps[0]->X, ps[1]->X);
         auto CD = Join(C, D);

         if (!isfinite(AB) && !isfinite(CD))
            throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__,
                                "both is not finite");
         else if (!isfinite(AB))
            return {{ps[1], ps[2], ps[3]}, {ps[3], ps[0], ps[1]}};  // CD
         else if (!isfinite(CD))
            return {{ps[0], ps[1], ps[2]}, {ps[2], ps[3], ps[0]}};  // AB
         else if (Max(AB) > Max(CD))
            return {{ps[1], ps[2], ps[3]}, {ps[3], ps[0], ps[1]}};  // CD
         else
            return {{ps[0], ps[1], ps[2]}, {ps[2], ps[3], ps[0]}};  // AB
      } catch (const error_message &e) {
         e.print();
         std::stringstream ss;
         //
         ss << MyVectorAngle(ps[2]->X - ps[1]->X, ps[0]->X - ps[1]->X) << std::endl;
         ss << MyVectorAngle(ps[0]->X - ps[3]->X, ps[2]->X - ps[3]->X) << std::endl;
         //
         ss << MyVectorAngle(ps[3]->X - ps[2]->X, ps[1]->X - ps[2]->X) << std::endl;
         ss << MyVectorAngle(ps[1]->X - ps[0]->X, ps[3]->X - ps[0]->X) << std::endl;
         //
         ss << MyVectorAngle(ps[1]->X - ps[0]->X, ps[2]->X - ps[0]->X) << std::endl;
         ss << MyVectorAngle(ps[3]->X - ps[1]->X, ps[0]->X - ps[1]->X) << std::endl;
         ss << MyVectorAngle(ps[0]->X - ps[2]->X, ps[1]->X - ps[2]->X) << std::endl;
         ss << MyVectorAngle(ps[1]->X - ps[3]->X, ps[2]->X - ps[3]->X) << std::endl;
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, ss.str());
      }
   };
   ///////////////////////////////////////////////////

   /*polygon::triangulate_detail
多角形の三角形分割．外角が`smallangle`よりも狭い三角形は対象にしない．
もし該当がなければ，`smallangle`を徐々に小さくしながら該当があるまで何回か繰り返す．
polygon::triangulate_detail*/
   /*polygon::triangulate_code*/
   // #define debug_triangle

   VV_i triangulate(const V_d &normal, double smallangle_IN = 0.) {
      //必ず1E-10よりも大きくなるように設定している
      // std::cout << __PRETTY_FUNCTION__ << std::endl;
      smallangle_IN += 1E-10;

      if (!isfinite(normal))
         throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__,
                             "this->normal is not finite"));

      this->activateAllPoints();

      VV_i ret;
      V_pp ps = this->getActivePoints();

      if (ps.size() < 3)
         return {};
      else if (ps.size() == 3) {
         if (isfinite(TriangleNormal(ps[0]->X, ps[1]->X, ps[2]->X)))
            return {{ps[0]->index, ps[1]->index, ps[2]->index}};
         else {
            std::stringstream ss;
            auto a = ps[0]->X, b = ps[1]->X, c = ps[2]->X;
            ss << "{a,b,c} = " << VV_d{a, b, c} << std::endl;
            ss << "n = " << TriangleNormal(a, b, c) << std::endl;
            ss << "TriangleAngles = " << TriangleAngles(a, b, c) << std::endl;
            ss << "TriangleArea = " << TriangleArea(a, b, c) << std::endl;
            throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__,
                                ss.str()));
         }
      } else if (ps.size() == 4) {
         return extractIndices(legalDivision4(ps, normal));
      } else if (ps.size() > 4) {
         do {
#ifdef debug_triangle
            std::cout << "   ps = " << ps << std::endl;
            std::cout << "  ret = " << ret << std::endl;
#endif

            this->calculateArea(ps, normal);
            this->calculateAngle(ps);
            int s = ps.size();

            point *select_p = nullptr;
            int current_index =
                0;  //このpsのインデックスで元々与えられたベクトルのインデックスではない
            bool found = getPointsMeetCondition(ps, smallangle_IN, select_p,
                                                current_index);

#ifdef debug_triangle
            std::cout << "found = " << found << std::endl;
            std::cout << "   ps = " << ps << std::endl;
            std::cout << "  ret = " << ret << std::endl;
#endif

            if (found) {
               select_p->active = false;
               ret.push_back({ps[(s + current_index - 1) % s]->index,
                              ps[current_index]->index,
                              ps[(current_index + 1) % s]->index});
            } else {
               //   break;
               // return ret;
               throw error_message(
                   __FILE__, __PRETTY_FUNCTION__, __LINE__,
                   "can not find triangle that meets conditions");
            }

            ps = this->getActivePoints();
         } while (ps.size() > 4);
      }
      VV_i two_trigs = extractIndices(legalDivision4(ps, normal));
      ret.insert(ret.end(), two_trigs.begin(), two_trigs.end());
      return ret;
   };

   // VV_i triangulate(const V_d &normal, double smallangle_IN = 0.)
   // {
   //   //必ず1E-10よりも大きくなるように設定している
   //   smallangle_IN += 1E-11;
   //   if (!isfinite(normal))
   //     throw(error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__,
   //     "this->normal is not finite"));

   //   this->activateAllPoints();

   //   VV_i ret;
   //   V_pp ps = this->getActivePoints();

   //   if (ps.size() < 3)
   //     return {};
   //   else if (ps.size() == 3)
   //     return {{ps[0]->index, ps[1]->index, ps[2]->index}};
   //   else if (ps.size() > 3)
   //   {
   //     do
   //     {
   //       this->calculateArea(ps, normal);
   //       this->calculateAngle(ps);

   //       point *select_p = nullptr;
   //       int index = 0;
   //       bool found = getPointsMeetCondition(ps, smallangle_IN, select_p,
   //       index); int s = ps.size(); if (found)
   //       {
   //         select_p->active = false;
   //         ret.push_back({ps[(s + index - 1) % s]->index, ps[index]->index,
   //         ps[(index + 1) % s]->index});
   //       }
   //       else
   //         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "can
   //         not find triangle that meets conditions");

   //       ps = this->getActivePoints();

   //     } while (ps.size() > 3);
   //   }

   //   ret.emplace_back(std::vector<int>{ps[0]->index, ps[1]->index,
   //   ps[2]->index}); return ret;
   // };
   /*polygon::triangulate_code*/
};

//-------------------------------------------------

class spherical_polygon : public geometry::polygon {
  public:
   V_pp points;
   V_d org;
   VV_d x_on_sphere;
   ~spherical_polygon(){};
   spherical_polygon(const V_d &org_IN, const VV_d &xyz_IN)
       : geometry::polygon(xyz_IN), org(org_IN) {
      VV_d x_on_sphere({});
      for (const auto &x : xyz_IN)
         x_on_sphere.emplace_back((x - org_IN) / Norm(x - org_IN));
   };
   void calculateArea(const V_pp &ps /*, const V_d& normal*/) {
      int s = ps.size();
      for (auto i = 0; i < s; i++) {
         ps[i]->area =
             SphericalInteriorArea(this->org, ps[(s + i - 1) % s]->X,
                                   ps[(s + i) % s]->X, ps[(s + i + 1) % s]->X);
      }
   };
   void calculateAngle(const V_pp &ps) {
      int s = ps.size();
      for (auto i = 0; i < s; i++) {
         ps[i]->angle =
             SphericalInteriorAngle(this->org, ps[(s + i - 1) % s]->X,
                                    ps[(s + i) % s]->X, ps[(s + i + 1) % s]->X);
      }
   };
   /*polygon::triangulate_detail
多角形の三角形分割．外角が`smallangle`よりも狭い三角形は対象にしない．
もし該当がなければ，`smallangle`を徐々に小さくしながら該当があるまで何回か繰り返す．
polygon::triangulate_detail*/
   /*polygon::triangulate_code*/
   VV_i triangulate(/*const V_d& normal, */ double smallangle_IN = 0.) {
      //必ず1E-10よりも大きくなるように設定している
      smallangle_IN += 1E-10;
      // if(!isfinite(normal))
      //   throw(error_message(__FILE__,__PRETTY_FUNCTION__,__LINE__,
      //   "this->normal is not finite"));

      this->activateAllPoints();

      VV_i ret;
      V_pp ps = this->getActivePoints();

      if (ps.size() < 3) {
         return {};
      } else if (ps.size() == 3) {
         return {{ps[0]->index, ps[1]->index, ps[2]->index}};
      } else if (ps.size() > 3) {
         do {
            this->calculateArea(ps /*, normal*/);
            this->calculateAngle(ps);

            point *select_p = nullptr;
            int index = 0;
            bool found = getPointsMeetCondition(ps, smallangle_IN, select_p, index);
            int s = ps.size();

            if (found) {
               select_p->active = false;
               ret.push_back({ps[(s + index - 1) % s]->index, ps[index]->index,
                              ps[(index + 1) % s]->index});
            } else {
               return ret;
            }

            ps = this->getActivePoints();

         } while (ps.size() > 3);
      }

      ret.emplace_back(
          std::vector<int>{ps[0]->index, ps[1]->index, ps[2]->index});
      return ret;
   };
   /*polygon::triangulate_code*/
};
/*polygon_code*/

// double SolidAngle(const V_d &org, const V_d &p0, const V_d &p1, const V_d
// &p2)
// {
//   //ccw p0 -> p1 -> p2
//   double cond = Dot(Cross(p1 - p0, p2 - p1), Mean(VV_d{p0- org, p1- org, p2-
//   org})); //Dot(triangle's normal vector, org's view direction vector) if
//   (cond > 1E-8)
//   {
//     return 4. * M_PI - (SphericalInteriorAngle(org, p0, p1, p2) +
//     SphericalInteriorAngle(org, p2, p0, p1) + SphericalInteriorAngle(org, p1,
//     p2, p0) - M_PI);
//   }
//   else if (cond < -1E-8)
//   {
//     return SphericalInteriorAngle(org, p0, p1, p2) +
//     SphericalInteriorAngle(org, p2, p0, p1) + SphericalInteriorAngle(org, p1,
//     p2, p0) - M_PI;
//   }
//   else
//   {
//     return 0.5;
//   }
// };

double SolidAngle_VanOosteromAandStrackeeJ1983(const Tddd &p, Tddd A, Tddd B, Tddd C) {
   // The solid angle of a plane triangle
   // Van Oosterom, A. and Strackee, J. (1983)
   auto [a0, a1, a2] = (A -= p);
   auto [b0, b1, b2] = (B -= p);
   auto [c0, c1, c2] = (C -= p);
   double nA = Norm(A);
   double nB = Norm(B);

   // double nC = Norm(C);
   // return 2. * atan2((-(a2 * b1 * c0) + a1 * b2 * c0 + a2 * b0 * c1 - a0 * b2
   // * c1 - a1 * b0 * c2 + a0 * b1 * c2), 				  (nB * (a0 * c0 + a1 * c1 + a2 * c2) +
   // nA * (b0 * c0 + b1 * c1 + b2 * c2) + (a0 * b0 + a1 * b1 + a2 * b2) * nC +
   // nA * nB * nC));

   return 2. * atan2(-(a2 * b1 * c0) + a1 * b2 * c0 + a2 * b0 * c1 - a0 * b2 * c1 - a1 * b0 * c2 + a0 * b1 * c2,
                     b0 * c0 * nA + b1 * c1 * nA + b2 * c2 * nA + a0 * c0 * nB + a1 * c1 * nB + a2 * c2 * nB + (a0 * b0 + a1 * b1 + a2 * b2 + nA * nB) * Norm(C));
};

double SolidAngle_VanOosteromAandStrackeeJ1983(const Tddd &p, const T3Tddd &ABC) {
   return SolidAngle_VanOosteromAandStrackeeJ1983(p, std::get<0>(ABC), std::get<1>(ABC), std::get<2>(ABC));
};

double SolidAngle(const Tddd &o, const Tddd &A, const Tddd &B, const Tddd &C) {
   double c = VectorAngle(A - o, B - o);
   double a = VectorAngle(B - o, C - o);
   double b = VectorAngle(C - o, A - o);
   double s = (a + b + c) * 0.5;
   if (Between(s, {M_PI - 1E-10, M_PI + 1E-10}))
      return 4. * M_PI / 2.;
   else
      return 4. * atan(sqrt(tan(s * 0.5) * tan((s - a) * 0.5) *
                            tan((s - b) * 0.5) * tan((s - c) * 0.5)));
   // return SolidAngle_VanOosteromAandStrackeeJ1983(o, A, B, C);
};

double SolidAngle(const Tddd &p, const T3Tddd &ABC) {
   return SolidAngle(p, std::get<0>(ABC), std::get<1>(ABC), std::get<2>(ABC));
   // return SolidAngle_VanOosteromAandStrackeeJ1983(p, std::get<0>(ABC),
   // std::get<1>(ABC), std::get<2>(ABC));
};

double SolidAngle(const Tddd &o, const std::vector<Tddd> &xyz) {
   double total = 0., tmp, angle;
   int sz = xyz.size();
   Tddd normal = {0., 0., 0.}, X1, X2;
   for (auto i = 0; i < xyz.size(); ++i) {
      X1 = Normalize(xyz[(i + sz) % sz] - o);
      X2 = Normalize(xyz[(i + sz + 1) % sz] - o);
      angle = VectorAngle(X1, X2);
      normal += angle * Normalize(Cross(X1, X2));
   }
   normal = -Normalize(normal);
   for (auto i = 0; i < xyz.size(); ++i) {
      tmp = geometry::SolidAngle(Tddd{0., 0., 0.},
                                 Normalize(xyz[(i + sz) % sz] - o),
                                 Normalize(xyz[(i + sz + 1) % sz] - o), normal);
      if (Between(tmp, {0., 4. * M_PI})) total += tmp;
   }
   return total;
};

double SolidAngle(const V_d &o, const V_d &A, const V_d &B, const V_d &C) {
   //   V_d center = Mean[{A, B, C}];
   V_d oA = (A - o);
   V_d oB = (B - o);
   V_d oC = (C - o);
   double c = MyVectorAngle(oA, oB);
   double a = MyVectorAngle(oB, oC);
   double b = MyVectorAngle(oC, oA);
   double s = (a + b + c) / 2.;
   double eps = 1E-10;
   if (Between(s, {M_PI - eps, M_PI + eps}))
      return 4. * M_PI / 2.;
   else
      return 4. * atan(sqrt(tan(s / 2) * tan((s - a) / 2.) * tan((s - b) / 2) *
                            tan((s - c) / 2)));
};

double SolidAngle(const V_d &org, const VV_d &ps) {
   // ccw
   double ret(0);
   geometry::spherical_polygon poly(org, ps);
   for (const auto &ind : poly.triangulate(1E-5))
      ret += geometry::SolidAngle(org, ps[ind[0]], ps[ind[1]], ps[ind[2]]);
   return ret;
};

//--------------------------
// //CCW p0p1p2 returns a positive value
// double SolidAngle(const V_d& p, const V_d& p0, const V_d& p1, const V_d& p2){
//   V_d r, cross = Cross(p1 - p0, p2 - p0);
//   V_d unit_n = cross;
//   double ret(0.);
//   for(const auto& abw:__GWGW6__){
//     r = p0*abw[0] + p1*abw[1] + p2*(1.-abw[0]-abw[1]) - p;
//     ret += Dot(r,unit_n)/pow(Norm(r),3.) * abw[2];
//   }
//   return std::abs(ret);
// };

// double SolidAngle(const V_d& p, const VV_d& X){
//   geometry::polygon poly(X);
//   double ret(0);
//   auto indices= poly.triangulate(Mean(X)-p);
//   std::cout << indices << std::endl;
//   for(const auto& ind:indices){
//     ret += geometry::SolidAngle(p, X[ind[0]], X[ind[1]], X[ind[2]]);
//   }
//   return ret;
// };

}  // namespace geometry

/* ------------------------------------------------------ */

T3Tddd closestTriangle(const Tddd &X, const std::vector<T3Tddd> &vector_vertex) {
   //!もしない場合はnullptrを返すので注意
   double min_distance = 1E+100;
   T3Tddd ret;
   for (const auto &vertex : vector_vertex) {
      auto intxn = IntersectionSphereTriangle(X, 1E+20, vertex);
      if (intxn.distance < min_distance) {
         ret = vertex;
         min_distance = intxn.distance;
      }
   }
   return ret;
};

Tddd closestPointOnTriangle(const Tddd &X, const T3Tddd &vertex) {
   auto intxn = IntersectionSphereTriangle(X, 1E+20, vertex);
   return intxn.X;
};

Tddd closestPointOnTriangle(const Tddd &X, std::vector<T3Tddd> &vector_vertex) {
   Tddd r = {1E+20, 1E+20, 1E+20}, ret, closestX;
   for (const auto &vertex : vector_vertex) {
      closestX = closestPointOnTriangle(X, vertex);
      if (Norm(r) >= Norm(closestX - X)) {
         r = closestX - X;
         ret = X;
      }
   }
   return ret;
};

Tddd vectorToClosestPointOnTriangle(const Tddd &X, std::vector<T3Tddd> &vector_vertex) {
   return closestPointOnTriangle(X, vector_vertex) - X;
};

Tddd closestPointOnLine(const Tddd &X, const T2Tddd &line) {
   auto intxn = IntersectionSphereLine(X, 1E+20, line);
   return intxn.X;
};

Tddd closestPointOnLine(const Tddd &X, const std::vector<T2Tddd> &vector_line) {
   Tddd r = {1E+20, 1E+20, 1E+20}, ret, closestX;
   for (const auto &line : vector_line) {
      closestX = closestPointOnLine(X, line);
      if (Norm(r) >= Norm(closestX - X)) {
         r = closestX - X;
         ret = X;
      }
   }
   return ret;
};

Tddd vectorToClosestPointOnLine(const Tddd &X, std::vector<T2Tddd> &vector_vertex) {
   return closestPointOnLine(X, vector_vertex) - X;
};

/* ------------------------------------------------------ */
template <typename T>
double windingNumber(const Tddd &X, const std::vector<std::tuple<T, T, T>> &V_vertices) {
   double ret = 0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ \
                                   : ret)
#endif
   for (const auto &vert : V_vertices)
      ret += geometry::SolidAngle_VanOosteromAandStrackeeJ1983(X, ToX(vert));
   return ret / (4. * M_PI);
};

template <>
double windingNumber(const Tddd &X, const std::vector<T3Tddd> &V_vertices) {
   double ret = 0;
#ifdef _OPENMP
#pragma omp parallel for reduction(+ \
                                   : ret)
#endif
   for (const auto &vertices : V_vertices)
      ret += geometry::SolidAngle_VanOosteromAandStrackeeJ1983(X, vertices);
   return ret / (4. * M_PI);
};
double windingNumber(const Tddd &X, const std::vector<Tddd> &V_vertices) {
   return 0.;
};
bool isInside(const Tddd &X, const std::vector<T3Tddd> &V_vertices) {
   return (CoordinateBounds(V_vertices).isInside(X) && windingNumber(X, V_vertices) < 0.75);
};
bool isInside(const Tddd &X, const std::vector<Tddd> &V_vertices) {
   return !CoordinateBounds(V_vertices).isInside(X);
};
bool isInside(const CoordinateBounds &bounds, const geometry::Sphere &s) {
   // cube v.s. sphere
   // cube < sphere ?
   auto [X0, X1] = std::get<0>(bounds.bounds);
   auto [Y0, Y1] = std::get<1>(bounds.bounds);
   auto [Z0, Z1] = std::get<2>(bounds.bounds);
   return (Norm(Tddd{X0, Y0, Z0} - s.X) < s.radius &&
           Norm(Tddd{X1, Y0, Z0} - s.X) < s.radius &&
           Norm(Tddd{X0, Y1, Z0} - s.X) < s.radius &&
           Norm(Tddd{X1, Y1, Z0} - s.X) < s.radius &&
           Norm(Tddd{X0, Y0, Z1} - s.X) < s.radius &&
           Norm(Tddd{X1, Y0, Z1} - s.X) < s.radius &&
           Norm(Tddd{X0, Y1, Z1} - s.X) < s.radius &&
           Norm(Tddd{X1, Y1, Z1} - s.X) < s.radius);
};
/* ------------------------------------------------------ */

Tddd ToX(const Tddd *X) { return *X; };
Tddd ToX(const Tddd X) { return X; };
T3Tddd ToX(const T3Tddd X) { return X; };

// b% -------------------------------------------------------------------------- */
// b%                                     八分木                                  */
// b% -------------------------------------------------------------------------- */
// 基本的な形状に対してoctreeを生成できるようにする．
template <typename T>
struct octree : public CoordinateBounds {
   /*
   top level -> depth = 0
   ---------
   depth = 1
   ---------
   last level -> depth = depthlimit = 2
   ---------
   */
   octree<T> *const parent;
   const int depth;
   const octree *top;  //最上階のデータ
   const std::vector<T> faces_only_for_top;
   std::vector<T> faces_;
   bool inside;
   std::vector<octree<T> *> children;
   /* ------------------------------------------------------ */
   octree(const CoordinateBounds &boundsIN, const Tii &depthlimit, const std::vector<T> &FACES)
       : CoordinateBounds(boundsIN), parent(nullptr), depth(0), top(this), faces_only_for_top(FACES), faces_({}), inside(windingNumber(boundsIN.getCenter(), top->faces_only_for_top) > 0.75), children(generateChildrenParallel(depthlimit, FACES)){};
   /* ------------------------------------------------------ */
   octree(const CoordinateBounds &boundsIN, const Tii &depthlimit, const std::vector<T> &FACES, octree<T> *const parentIN)
       : CoordinateBounds(boundsIN), parent(parentIN), depth(parentIN->depth + 1), top(parentIN->top), faces_only_for_top({}), faces_({}), inside((parentIN && FACES.empty()) ? parentIN->inside : (windingNumber(boundsIN.getCenter(), top->faces_only_for_top) > 0.75)), children(generateChildrenParallel(depthlimit, FACES)){};
   /* ------------------------------------------------------ */
   std::vector<octree<T> *> generateChildrenParallel(const Tii &depthlimit, const std::vector<T> &FACES) {
      // std::vector<T> faces_;
      faces_.reserve(FACES.size());
      for (const auto &f : FACES)
         if (intersectingQ(this->bounds, ToX(f))) faces_.emplace_back(f);
      if (std::get<0>(depthlimit) >= this->depth && std::get<1>(depthlimit) <= faces_.size()) {
         auto [b0, b1, b2, b3, b4, b5, b6, b7] = to8Bounds();
#ifdef _OPENMP
         std::vector<octree<T> *> ret(8);
#pragma omp parallel sections
         {
#pragma omp section
            ret[0] = new octree(b0, depthlimit, faces_, this);
#pragma omp section
            ret[1] = new octree(b1, depthlimit, faces_, this);
#pragma omp section
            ret[2] = new octree(b2, depthlimit, faces_, this);
#pragma omp section
            ret[3] = new octree(b3, depthlimit, faces_, this);
#pragma omp section
            ret[4] = new octree(b4, depthlimit, faces_, this);
#pragma omp section
            ret[5] = new octree(b5, depthlimit, faces_, this);
#pragma omp section
            ret[6] = new octree(b6, depthlimit, faces_, this);
#pragma omp section
            ret[7] = new octree(b7, depthlimit, faces_, this);
         }
         return ret;
#else
         return {new octree(b0, depthlimit, faces_, this), new octree(b1, depthlimit, faces_, this), new octree(b2, depthlimit, faces_, this), new octree(b3, depthlimit, faces_, this),
                 new octree(b4, depthlimit, faces_, this), new octree(b5, depthlimit, faces_, this), new octree(b6, depthlimit, faces_, this), new octree(b7, depthlimit, faces_, this)};
#endif
      } else
         return {};
   };
   std::vector<octree<T> *> generateChildren(const Tii &depthlimit, std::vector<T> &FACES) {
      // std::vector<T> faces_;
      faces_.reserve(FACES.size());
      for (const auto &f : FACES)
         if (intersectingQ(this->bounds, ToX(f))) faces_.emplace_back(f);
      if (std::get<0>(depthlimit) >= this->depth && std::get<1>(depthlimit) <= faces_.size()) {
         //@ min#Objects より多くなければ分割されない
         //@ maxDepth までしか回想は作れればい．最低は０階
         //例えば，{10,1}の場合，最大で10階まで分割されている．またセルが１個含んでいればそれ以上分割されない．
         //@ また，オブジェクトがゼロなら１つ目の条件から分割されない
         auto [b0, b1, b2, b3, b4, b5, b6, b7] = to8Bounds();
         return {new octree(b0, depthlimit, faces_, this), new octree(b1, depthlimit, faces_, this),
                 new octree(b2, depthlimit, faces_, this), new octree(b3, depthlimit, faces_, this),
                 new octree(b4, depthlimit, faces_, this), new octree(b5, depthlimit, faces_, this),
                 new octree(b6, depthlimit, faces_, this), new octree(b7, depthlimit, faces_, this)};
      } else
         return {};
   };
   /* ------------------------------------------------------ */
   octree(const CoordinateBounds &boundsIN, const std::vector<T> &FACES, const Tii &depthlimit)
       : CoordinateBounds(boundsIN), parent(nullptr), depth(0), top(this), faces_only_for_top(FACES), faces_(voronoiConnectivity(FACES)), children(voronoiChildren(depthlimit, faces_)){};
   octree(const CoordinateBounds &boundsIN, const std::vector<T> &FACES, const Tii &depthlimit, octree<T> *const parentIN)
       : CoordinateBounds(boundsIN), parent(parentIN), depth(parentIN->depth + 1), top(parentIN->top), faces_(voronoiConnectivity(FACES)), children(voronoiChildren(depthlimit, faces_)){};
   std::vector<T> voronoiConnectivity(const std::vector<T> &FACES) const {
      Tdd minmax, r = {1E+10, 1E+10};
      std::vector<T> ret;
      for (const auto &f : FACES) {
         minmax = Distance(ToX(f));
         if (std::get<0>(r) >= std::get<0>(minmax))
            std::get<0>(r) = std::get<0>(minmax);
         if (std::get<1>(r) >= std::get<1>(minmax))
            std::get<1>(r) = std::get<1>(minmax);
      }
      ret.reserve(FACES.size());
      for (const auto &f : FACES) {
         minmax = Distance(ToX(f));
         if (Between(std::get<0>(minmax), r) || Between(std::get<1>(minmax), r) || Between(std::get<0>(r), minmax))
            ret.emplace_back(f);
      }
      return ret;
   };
   std::vector<octree<T> *> voronoiChildren(const Tii &depthlimit, const std::vector<T> &faces) {
      if (faces.empty() || std::get<0>(depthlimit) <= this->depth || std::get<1>(depthlimit) >= faces.size()) {
         return {};
      } else {
         std::vector<octree<T> *> ret(8);
         auto [b0, b1, b2, b3, b4, b5, b6, b7] = to8Bounds();
#pragma omp parallel sections
         {
#pragma omp section
            ret[0] = new octree<T>(b0, faces, depthlimit, this);
#pragma omp section
            ret[1] = new octree<T>(b1, faces, depthlimit, this);
#pragma omp section
            ret[2] = new octree<T>(b2, faces, depthlimit, this);
#pragma omp section
            ret[3] = new octree<T>(b3, faces, depthlimit, this);
#pragma omp section
            ret[4] = new octree<T>(b4, faces, depthlimit, this);
#pragma omp section
            ret[5] = new octree<T>(b5, faces, depthlimit, this);
#pragma omp section
            ret[6] = new octree<T>(b6, faces, depthlimit, this);
#pragma omp section
            ret[7] = new octree<T>(b7, faces, depthlimit, this);
         }
         return ret;
      }
   };
   void deleteOuside() {
      if (!this->children.empty()) {
         auto tmp = this->children;
         for (auto &c : tmp) c->deleteOuside();
      } else if (!this->inside)
         delete this;
   };
   ~octree() {
      auto tmp = this->children;
      for (const auto &c : tmp) delete c;
      if (parent)
         parent->children.erase(std::remove(parent->children.begin(), parent->children.end(), this), parent->children.end());
   };
   void getDescendants(auto &accum) const {
      for (const auto &c : this->children) {
         accum.emplace(c);
         c->getDescendants(accum);
      }
   };
   bool hasChildren() const { return !this->children.empty(); };
   /* ------------------------------------------------------ */
   template <typename U>
   bool isAllVertexInsideOf(const U &object) const {
      /*
      このセルの頂点全てが，何かsの中に入っているか？
      all vertices of this cube are in the given object ?
      */
      auto [X0, X1, X2, X3, X4, X5, X6, X7] = this->getVertices();
      return IntersectQ(object, X0) &&
             IntersectQ(object, X1) &&
             IntersectQ(object, X2) &&
             IntersectQ(object, X3) &&
             IntersectQ(object, X4) &&
             IntersectQ(object, X5) &&
             IntersectQ(object, X6) &&
             IntersectQ(object, X7);
   };
   bool isAllVertexInsideOf(const Tddd &s) const { return false; };
   bool isAllVertexInsideOf(const Sphere &s) const {
      auto [X0, X1, X2, X3, X4, X5, X6, X7] = this->getVertices();
      return (Norm(X0 - s.center) <= s.radius && Norm(X1 - s.center) <= s.radius &&
              Norm(X2 - s.center) <= s.radius && Norm(X3 - s.center) <= s.radius &&
              Norm(X4 - s.center) <= s.radius && Norm(X5 - s.center) <= s.radius &&
              Norm(X6 - s.center) <= s.radius && Norm(X7 - s.center) <= s.radius);
   };
   // b* -------------------------------------------------------------------------- */
   // b*                                INTERSECTIONS                               */
   // b* -------------------------------------------------------------------------- */
   std::unordered_set<octree<T> *> getIntersectAsUnorderedSet(const auto &s) const {  //交わる最深階層にあるキューブ
      std::unordered_set<octree<T> *> accum;
      accum.reserve(100000);
      for (const auto &c : this->children) c->getIntersect(accum, s);
      return accum;
   };
   void getIntersect(std::unordered_set<octree<T> *> &accum, const auto &s) {
      if (!this->children.empty())
         for (const auto &c : this->children) c->getIntersect(accum, s);
      else if (this->inside && IntersectQ(this->bounds, s))
         accum.emplace(this);
   };
   void getAllDeepest(std::unordered_set<octree<T> *> &accum) {
      if (!this->children.empty())
         for (const auto &c : this->children) c->getAllDeepest(accum);
      else
         accum.emplace(this);
   };
   void getAllDeepestInside(std::unordered_set<octree<T> *> &accum) {
      if (!this->children.empty())
         for (const auto &c : this->children) c->getAllDeepestInside(accum);
      else if (this->inside)
         accum.emplace(this);
   };
   void getAllDeepestInside(std::vector<octree<T> *> &accum) {
      if (!this->children.empty())
         for (const auto &c : this->children) c->getAllDeepestInside(accum);
      else if (this->inside)
         accum.emplace_back(this);
   };
   std::unordered_set<octree<T> *> getAllDeepest() {
      std::unordered_set<octree<T> *> ret;
      if (!this->children.empty())
         for (const auto &c : this->children) c->getAllDeepest(ret);
      else
         ret.emplace(this);
      return ret;
   };
   /* ------------------------------------------------------ */
   void getDepth(const int &d, std::unordered_set<octree<T> *> &accum) {
      if (this->depth < d && !this->children.empty())
         for (const auto &c : this->children) c->getDepth(d, accum);
      else if (this->depth == d)
         accum.emplace(this);
   };
   std::unordered_set<octree<T> *> getDepth(const int &d) {
      std::unordered_set<octree<T> *> ret;
      if (this->depth < d && !this->children.empty())
         for (const auto &c : this->children) c->getDepth(d, ret);
      else if (this->depth == d)
         ret.emplace(this);
      return ret;
   };
   /* ------------------------------------------------------ */
   std::unordered_set<octree<T> *> getAllDeepestInside() {
      std::unordered_set<octree<T> *> ret;
      if (!this->children.empty())
         for (const auto &c : this->children) c->getAllDeepestInside(ret);
      else if (this->inside)  //全てあるかないかなので，一つでもnullptrなら，最後の階層と考えられる．
         ret.emplace(this);
      return ret;
   };
   /* ------------------------------------------------------ */
   void isIntersectInside(bool &found, const auto &s) {
      if (!found && IntersectQ(this->bounds, s)) {
         if (!this->children.empty()) {
            for (const auto &c : this->children)
               c->isIntersectInside(found, s);
         } else if (this->inside)
            found = true;
      }
   };
   bool isIntersectInside(const auto &s) const {  //交わる最深階層にあるキューブ
      bool found = false;
      for (const auto &c : this->children) c->isIntersectInside(found, s);
      return found;
   };
   /* ------------------------------------------------------- */
   void getIntersect(std::vector<octree<T> *> &accum, const auto &s) {
      if (IntersectQ(this->bounds, s)) {
         if (!this->children.empty())
            for (const auto &c : this->children) c->getIntersect(accum, s);
         else
            accum.emplace_back(this);
      }
   };
   std::vector<octree<T> *> getIntersect(const auto &s) const {  //交わる最深階層にあるキューブ
      std::vector<octree<T> *> accum;
      accum.reserve(100000);
      for (const auto &c : this->children) c->getIntersect(accum, s);
      return accum;
   };
   /* ------------------------------------------------------- */
   void getIntersect(octree<T> *&p_cell, const Tddd &X) {
      if (!p_cell && IntersectQ(this->bounds, X)) {
         if (!this->children.empty()) {
            for (const auto &c : this->children)
               c->getIntersect(p_cell, X);
         } else if (this->inside)
            p_cell = this;
      }
   };
   octree<T> *getIntersect(const Tddd &X) const {  //交わる最深階層にあるキューブ
      octree<T> *p_cell = nullptr;
      for (const auto &c : this->children) c->getIntersect(p_cell, X);
      return p_cell;
   };

   /* ------------------------------------------------------- */
   void getIntersectInside(std::vector<octree<T> *> &accum, const auto &s) {
      if (IntersectQ(this->bounds, s)) {
         if (!this->children.empty()) {    //! c->c0とせずにc0とだけすることによる問題が多い
                                           // b! 今のエラーは，ここのエラーではない．sphereの干渉チェックの問題
            if (isAllVertexInsideOf(s)) {  //全ての頂点がsの内部にあるので，以降，干渉チェックは行わない．
               for (const auto &c : this->children)
                  c->getAllDeepestInside(accum);
            } else {
               for (const auto &c : this->children)
                  c->getIntersectInside(accum, s);
            }
         } else if (this->inside)
            accum.emplace_back(this);
      }
   };
   std::vector<octree<T> *> getIntersectInside(const auto &s) const {  //交わる最深階層にあるキューブ
      std::vector<octree<T> *> accum;
      accum.reserve(100000);
      for (const auto &c : this->children) c->getIntersectInside(accum, s);
      return accum;
   };
   void apply(const std::function<void(octree<T> *)> &func) const {  //交わる最深階層にあるキューブ
      for (const auto &c : this->children)
         func(c);
   };
   void apply(const std::function<void(octree<T> *)> &func, const auto &s) {
      if (IntersectQ(this->bounds, s)) {
         if (!this->children.empty()) {  //! c->c0とせずにc0とだけすることによる問題が多い
            for (const auto &c : this->children)
               c->apply(func, s);
         } else if (this->inside)
            func(this);
      }
   };
   // b! -------------------------------------------------------------------------- */
   // b!                              SETTING NEIGHBORS                             */
   // b! -------------------------------------------------------------------------- */
   std::unordered_set<octree<T> *> neighbors;
   std::tuple<T, T, T, T, T, T, T, T> nearest_face;
   std::unordered_set<T> nearest_faces;
   std::tuple<bool, bool, bool, bool, bool, bool, bool, bool> bools;
   T8d scalers;
   T8Tddd vectors;
   std::unordered_set<T> checked_faces_passed;
   template <typename U>
   U Interpolate(const Tddd &X, const std::tuple<U, U, U, U, U, U, U, U> &c) const {
      const auto [x, y, z] = X;
      const auto [c000, c001, c010, c011, c100, c101, c110, c111] = c;
      const auto [x0, x1] = std::get<0>(this->bounds);
      const auto [y0, y1] = std::get<1>(this->bounds);
      const auto [z0, z1] = std::get<2>(this->bounds);
      return (-(((c111 * (x - x0) + c011 * (-x + x1)) * (y - y0) + (c101 * (-x + x0) + c001 * (x - x1)) * (y - y1)) * (z - z0)) +
              ((c110 * (x - x0) + c010 * (-x + x1)) * (y - y0) + (c100 * (-x + x0) + c000 * (x - x1)) * (y - y1)) * (z - z1)) /
             ((x0 - x1) * (y0 - y1) * (z0 - z1));
   };
   template <typename U>
   U Integrate(const std::tuple<U, U, U, U, U, U, U, U> &c) const {
      const auto [c000, c001, c010, c011, c100, c101, c110, c111] = c;
      const auto [x0, x1] = std::get<0>(this->bounds);
      const auto [y0, y1] = std::get<1>(this->bounds);
      const auto [z0, z1] = std::get<2>(this->bounds);
      return -((c000 + c001 + c010 + c011 + c100 + c101 + c110 + c111) * (x0 - x1) * (y0 - y1) * (z0 - z1)) / 8.;
   };
   /* -------------------------------------------------------------------------- */
   void setNeighbors() {
      /*
      最深層のセルの近傍セルを
      neighbors
      にセットする．
      */
      double distance;
      for (const auto &c : this->getAllDeepest()) {
         auto [X0, X1] = std::get<0>(c->bounds);
         auto [Y0, Y1] = std::get<1>(c->bounds);
         auto [Z0, Z1] = std::get<2>(c->bounds);
         auto d0 = (X1 - X0) * 0.2;
         auto d1 = (Y1 - Y0) * 0.2;
         auto d2 = (Z1 - Z0) * 0.2;
         for (const auto &v : top->getIntersectInside(T3Tdd{{X0 - d0, X1 + d0}, {Y0 - d1, Y1 + d1}, {Z0 - d2, Z1 + d2}})) {
            if (v != c)
               c->neighbors.emplace(v);
         }
      }
   };
   // b! -------------------------------------------------------------------------- */
   // template <typename = typename std::enable_if<std::is_same<T, T3Tddd>::value>::type>
   // void setVectorsToTriangle() {
   //    auto tmp = this->getAllDeepestInside();
   //    for (const auto &cell : tmp) {
   //       cell->nearest_faces.clear();
   //       cell->checked_faces_passed.clear();
   //       cell->bools = {false, false, false, false, false, false, false, false};
   //       for (const auto &f : cell->faces_) {
   //          // 各頂点にとって最も近い点を抽出
   //          for_each01111(cell->getVertices(), cell->scalers, cell->vectors, cell->nearest_face, cell->bools,
   //                        [&](const auto &x, auto &s, auto &v, auto &f4v, auto &b) {
   //                           auto intsp = IntersectionSphereTriangle(x, 1E+10, ToX(f));
   //                           if (intsp.isIntersecting && (intsp.distance <= s || !b)) {
   //                              s = intsp.distance;
   //                              v = intsp.X - x;
   //                              f4v = f;
   //                              b = true;
   //                           }
   //                        });
   //       }
   //    }

   //    for (int i = 0; i < 10; ++i)
   //       for (const auto &cell : tmp) {  // cellにとって，最も近い面を，neighborsから探す
   //          for (const auto &nei : cell->neighbors) {
   //             for (const auto &f : nei->nearest_faces)
   //                if ((cell->checked_faces_passed.emplace(f)).second) {
   //                   // std::cout << "各頂点にとって最も近い点を抽出" << std::endl;
   //                   for_each01111(cell->getVertices(), cell->scalers, cell->vectors, cell->nearest_face, cell->bools,
   //                                 [&](const auto &x, auto &s, auto &v, auto &f4v, auto &b) {
   //                                    auto intsp = IntersectionSphereTriangle(x, 1E+10, ToX(f));
   //                                    if (intsp.isIntersecting && (intsp.distance <= s || !b)) {
   //                                       s = intsp.distance;
   //                                       v = intsp.X - x;
   //                                       f4v = f;
   //                                       b = true;
   //                                    }
   //                                 });
   //                }
   //          }
   //       }
   // };
   /* -------------------------------------------------------------------------- */
   template <typename = typename std::enable_if<std::is_same<T, T3Tddd>::value>::type>
   void setVectorsToTriangle() {
      auto tmp = this->getAllDeepestInside();
      for (const auto &cell : tmp) {
         cell->checked_faces_passed.clear();
         cell->bools = {false, false, false, false, false, false, false, false};
         for (const auto &f : cell->faces_) {
            // 各頂点にとって最も近い点を抽出
            for_each01111(cell->getVertices(), cell->scalers, cell->vectors, cell->nearest_face, cell->bools,
                          [&](const auto &x, auto &s, auto &v, auto &f4v, auto &b) {
                             auto intsp = IntersectionSphereTriangle(x, 1E+10, ToX(f));
                             if (intsp.isIntersecting && (intsp.distance <= s || !b)) {
                                s = intsp.distance;
                                v = intsp.X - x;
                                f4v = f;
                                b = true;
                             }
                          });
         }
      }
      //
      for (int i = 0; i < 5; ++i)
         for (const auto &cell : tmp) {  // cellにとって，最も近い面を，neighborsから探す
            for (const auto &nei : cell->neighbors) {
               for_each(nei->nearest_face, nei->bools,
                        [&](const auto &f, const auto &B) {
                           if (B && (cell->checked_faces_passed.emplace(f)).second) {
                              // std::cout << "各頂点にとって最も近い点を抽出" << std::endl;
                              for_each01111(cell->getVertices(), cell->scalers, cell->vectors, cell->nearest_face, cell->bools,
                                            [&](const auto &x, auto &s, auto &v, auto &f4v, auto &b) {
                                               auto intsp = IntersectionSphereTriangle(x, 1E+10, ToX(f));
                                               if (intsp.isIntersecting && (intsp.distance <= s || !b)) {
                                                  s = intsp.distance;
                                                  v = intsp.X - x;
                                                  f4v = f;
                                                  b = true;
                                               }
                                            });
                           }
                        });
               for (const auto &f : nei->faces_)
                  if (cell->checked_faces_passed.emplace(f).second) {
                     // std::cout << "各頂点にとって最も近い点を抽出" << std::endl;
                     for_each01111(cell->getVertices(), cell->scalers, cell->vectors, cell->nearest_face, cell->bools,
                                   [&](const auto &x, auto &s, auto &v, auto &f4v, auto &b) {
                                      auto intsp = IntersectionSphereTriangle(x, 1E+10, ToX(f));
                                      if (intsp.isIntersecting && (intsp.distance <= s || !b)) {
                                         s = intsp.distance;
                                         v = intsp.X - x;
                                         f4v = f;
                                         b = true;
                                      }
                                   });
                  }
            }
         }
   };
};
#endif