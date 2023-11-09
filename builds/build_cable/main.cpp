#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include "basic_IO.hpp"
#include "basic_arithmetic_array_operations.hpp"

/*DOC_EXTRACT

オイラー法でケーブルの動きをシミュレーションする．

![sample.gif](./sample.gif)

*/

struct Node {
   std::array<double, 3> X;
   std::array<double, 3> velocity;
   std::array<double, 3> acceleration;
};

const double stiffness = 1000;
const double dump = 0.5;
const double dt = 0.01;  // Time step
const int max_step = 2000;
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

   for (size_t i = 1; i < nodes.size(); ++i) {
      nodes[i].acceleration = gravity;
      if (i - 1 >= 0) {
         auto v = nodes[i - 1].X - nodes[i].X;
         double disp = Norm(v) - natural_length;
         if (disp > 0.)
            nodes[i].acceleration += stiffness * disp * Normalize(v);
      }
      if (i + 1 < nodes.size()) {
         auto v = nodes[i + 1].X - nodes[i].X;
         double disp = Norm(v) - natural_length;
         if (disp > 0.)
            nodes[i].acceleration += stiffness * disp * Normalize(v);
      }
      // Damping
      nodes[i].acceleration -= dump * nodes[i].velocity;
   }

   double T = 3;
   double w = 2 * M_PI / T;

   for (auto i = 0; i < nodes.size(); ++i) {
      nodes[i].velocity += nodes[i].acceleration * dt;

      //! 最後の節点は，10まで最後の接点をぐるぐる回す
      if (i == nodes.size() - 1) {
         if (t < 10) {
            // nodes[i].velocity = Cross(nodes[i].X, {0., 0., 1.});
            nodes[i].velocity = {0., 0., 10 * cos(w * t)};
         } else
            nodes[i].velocity.fill(0.);
      }
      nodes[i].X += nodes[i].velocity * dt;
   }
}

/* -------------------------------------------------------------------------- */

int main() {
   // Example setup
   std::vector<double> times;
   double t = 0;
   std::vector<std::vector<std::array<double, 3>>> positions;
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
