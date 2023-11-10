#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include "basic_IO.hpp"
#include "basic_arithmetic_array_operations.hpp"
#include "basic_vectors.hpp"
#include "integrationOfODE.hpp"

/*DOC_EXTRACT

実行方法：

```sh
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../
make
```

オイラー法でケーブルの動きをシミュレーションする．

![sample.gif](./sample.gif)

*/

#define USE_EULER
// #define USE_LEAP_FROG
// #define USE_RK4

struct Node {
   std::array<double, 3> X;
   std::array<double, 3> velocity;
   std::array<double, 3> acceleration;
   //
   Node(std::array<double, 3> X, std::array<double, 3> velocity, std::array<double, 3> acceleration)
       : X(X), velocity(velocity), acceleration(acceleration) {}

   double t;
   LeapFrog<decltype(X)> LPFG;

   RungeKutta<decltype(X)> RK_x;
   RungeKutta<decltype(X)> RK_v;
};

const double dt = 0.01;  // Time step
const int max_step = 6000;
const std::array<double, 3> gravity = {0, 0, -9.81};

// nodes 20
std::vector<Node> nodes = {
    {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    {{1, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    {{2, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    {{3, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    {{4, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    {{5, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    {{6, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    {{7, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    {{8, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    {{9, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    {{10, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    {{11, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    {{12, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    {{13, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    {{14, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    {{15, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    {{16, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    {{17, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    {{18, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    {{19, 0, 0}, {0, 0, 0}, {0, 0, 0}},
    {{20, 0, 0}, {0, 0, 0}, {0, 0, 0}}};

/* -------------------------------------------------------------------------- */

void simulateCableDynamics(double t, double dt) {

   auto tension = [&dt](const int i) {
      const double stiffness = 1000;
      const double damp = .5;
      const double natural_length = 1.;
      std::array<double, 3> acceleration;
      if (i - 1 >= 0) {
         auto v = nodes[i - 1].X - nodes[i].X;
         double disp = Norm(v) - natural_length;
         if (disp > 0.)
            acceleration += stiffness * disp * Normalize(v);
         auto relative_velocity = nodes[i - 1].velocity - nodes[i].velocity;
         acceleration -= damp * Projection(-relative_velocity, v) / dt;
      }
      if (i + 1 < nodes.size()) {
         auto v = nodes[i + 1].X - nodes[i].X;
         double disp = Norm(v) - natural_length;
         if (disp > 0.)
            acceleration += stiffness * disp * Normalize(v);
         auto relative_velocity = nodes[i + 1].velocity - nodes[i].velocity;
         acceleration -= damp * Projection(-relative_velocity, v) / dt;
      }
      return acceleration;
   };

   double T = 3;
   double w = 2 * M_PI / T;

#ifdef USE_LEAP_FROG
   for (size_t step = 0; step < 2; ++step)
#elif defined USE_RK4
   for (size_t step = 0; step < 4; ++step)
#endif
   {
      for (size_t i = 1; i < nodes.size(); ++i) {

         nodes[i].acceleration = tension(i) + gravity;

         if (t < 30)
            if (i == nodes.size() - 1) {
               //! 最後の節点は，10まで最後の接点をぐるぐる回す
               // nodes[i].velocity = Cross(nodes[i].X, {0., 0., 1.});
               if (Between(t, {2 * T, 3 * T}) || Between(t, {4 * T, 5 * T}))
                  nodes[i].acceleration = {0., 50 * cos(w * t), 50 * cos(w * t)};
               else
                  nodes[i].acceleration = -nodes[i].velocity / dt;
            }
      }

#ifdef USE_EULER
      //! オイラー法
      for (size_t i = 1; i < nodes.size(); ++i) {
         nodes[i].velocity += nodes[i].acceleration * dt;
         nodes[i].X += nodes[i].velocity * dt;
      }
#elif defined USE_LEAP_FROG
      for (size_t i = 1; i < nodes.size(); ++i) {
         nodes[i].LPFG.push(nodes[i].acceleration);
         nodes[i].t = nodes[i].LPFG.get_t();
         nodes[i].X = nodes[i].LPFG.get_x();
         nodes[i].velocity = nodes[i].LPFG.get_v();
      }
#elif defined USE_RK4
      for (size_t i = 1; i < nodes.size(); ++i) {
         nodes[i].RK_x.push(nodes[i].velocity);
         nodes[i].RK_v.push(nodes[i].acceleration);
         nodes[i].t = nodes[i].RK_x.get_t();
         nodes[i].X = nodes[i].RK_x.get_x();
         nodes[i].velocity = nodes[i].RK_v.get_x();
      }
#endif
   }
}

/* -------------------------------------------------------------------------- */

int main() {
   // Example setup
   std::vector<double> times;
   double t = 0;
   std::vector<std::vector<std::array<double, 3>>> positions;

#ifdef USE_LEAP_FROG
   for (size_t i = 1; i < nodes.size(); ++i)
      nodes[i].LPFG.initialize(dt, t, nodes[i].X, nodes[i].velocity);
#endif

   for (int i = 0; i < max_step; ++i) {
#ifdef USE_RK4
      for (size_t i = 1; i < nodes.size(); ++i) {
         nodes[i].RK_x.initialize(dt, t, nodes[i].X, 4);
         nodes[i].RK_v.initialize(dt, t, nodes[i].velocity, 4);
      }
#endif
      std::cout << "t = " << t << std::endl;
      simulateCableDynamics(t, dt);
      times.push_back(t);
      positions.push_back({});

      for (const Node& node : nodes) {
         if (!isFinite(node.X))
            break;
         positions.back().push_back(node.X);
      }
      t += dt;
   }

   // Output to JSON file
#ifdef USE_EULER
   std::ofstream myfile("cable_euler.json");
#elif defined USE_LEAP_FROG
   std::ofstream myfile("cable_leapfrog.json");
#elif defined USE_RK4
   std::ofstream myfile("cable_rk4.json");
#endif
   myfile << "{";
   std::array<std::string, 2> bracket = {"[", "]"};
   myfile << "\"time\" : " << std::make_tuple(times, bracket) << ", ";
   myfile << "\"position\" : " << std::make_tuple(positions, bracket);
   myfile << "}";
   myfile.close();
   return 0;
}
