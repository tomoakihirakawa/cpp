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

struct Node {
   std::array<double, 3> X;
   std::array<double, 3> velocity;
   std::array<double, 3> acceleration;
   //
   Node(std::array<double, 3> X, std::array<double, 3> velocity, std::array<double, 3> acceleration)
       : X(X), velocity(velocity), acceleration(acceleration) {}

   double t;
   LeapFrog<decltype(X)> LPFG;
};

const double stiffness = 10000;
const double damp = 0.5;
const double dt = 0.01;  // Time step
const int max_step = 5000;
const double natural_length = 1.;
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

   double T = 3;
   double w = 2 * M_PI / T;

   auto tension = [&](const int i) {
      std::array<double, 3> acceleration;
      acceleration.fill(0.);
      if (i - 1 >= 0) {
         auto v = nodes[i - 1].X - nodes[i].X;
         double disp = Norm(v) - natural_length;
         if (disp > 0.)
            acceleration += stiffness * disp * Normalize(v);
      }
      if (i + 1 < nodes.size()) {
         auto v = nodes[i + 1].X - nodes[i].X;
         double disp = Norm(v) - natural_length;
         if (disp > 0.)
            acceleration += stiffness * disp * Normalize(v);
      }
      return acceleration;
   };

   auto accel = [&](const auto& velocity) {
      //! 最後の節点は，10まで最後の接点をぐるぐる回す
      // nodes[i].velocity = Cross(nodes[i].X, {0., 0., 1.});
      if (Between(t, {T, 2 * T}) || Between(t, {3 * T, 4 * T}))
         return std::array<double, 3>{0., 50 * cos(w * t), 50 * cos(w * t)};
      else
         return -velocity / dt;
   };

   for (size_t step = 0; step < 2; ++step) {

      for (size_t i = 1; i < nodes.size(); ++i) {

         nodes[i].acceleration = tension(i) + gravity;
         nodes[i].acceleration -= damp * nodes[i].velocity;

         if (t < 30)
            if (i == nodes.size() - 1)
               nodes[i].acceleration = accel(nodes[i].velocity);
      }

      //! オイラー法
      // for (size_t i = 1; i < nodes.size(); ++i) {
      //    nodes[i].velocity += nodes[i].acceleration * dt;
      //    nodes[i].X += nodes[i].velocity * dt;
      // }
      for (size_t i = 1; i < nodes.size(); ++i) {
         nodes[i].LPFG.push(nodes[i].acceleration);
         //
         nodes[i].t = nodes[i].LPFG.get_t();
         nodes[i].X = nodes[i].LPFG.get_x();
         nodes[i].velocity = nodes[i].LPFG.get_v();
      }
   }
}

/* -------------------------------------------------------------------------- */

int main() {
   // Example setup
   std::vector<double> times;
   double t = 0;
   std::vector<std::vector<std::array<double, 3>>> positions;

   for (size_t i = 1; i < nodes.size(); ++i)
      nodes[i].LPFG.initialize(dt, t, nodes[i].X, nodes[i].velocity);

   for (int i = 0; i < max_step; ++i) {
      t += dt;
      std::cout << "t = " << t << std::endl;
      simulateCableDynamics(t, dt);
      times.push_back(t);
      positions.push_back({});
      for (const Node& node : nodes)
         positions.back().push_back(node.X);
   }

   // Output to JSON file
   std::ofstream myfile("cable.json");
   myfile << "{";
   std::array<std::string, 2> bracket = {"[", "]"};
   myfile << "\"time\" : " << std::make_tuple(times, bracket) << ", ";
   myfile << "\"position\" : " << std::make_tuple(positions, bracket);
   myfile << "}";
   myfile.close();
   return 0;
}
