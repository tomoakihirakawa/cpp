#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include "MooringLine.hpp"
#include "basic_IO.hpp"
#include "basic_arithmetic_array_operations.hpp"
#include "basic_vectors.hpp"
#include "vtkWriter.hpp"

// (cd builds/build_cable; python3 ../../extract_comments.py README.md ./ ../../)
// echo '(cd builds/build_cable; python3 ../../extract_comments.py README.md ./ ../../)'

auto homeDir = getenv("HOME");
std::string basePath = std::string(homeDir) + "/Cable/";
PVDWriter pvd_line0(basePath + "line0.pvd");
PVDWriter pvd_line1(basePath + "line1.pvd");
PVDWriter pvd_line2(basePath + "line2.pvd");
PVDWriter pvd_points0(basePath + "points0.pvd");
PVDWriter pvd_points1(basePath + "points1.pvd");
PVDWriter pvd_points2(basePath + "points2.pvd");
/* -------------------------------------------------------------------------- */
int main() {

   const double stiffness = 14 * std::pow(10, 8);  //! [N/m]
   const double damp = .5;                         //! [N/(m/s^2)]
   const double density = 348.5;                   //! [kg/m]
   const double dt = 0.01;                         //! [s]
   const double total_length = 522.;               //! [m]
   const int n_points = 100;
   const double diam = 0.1;

   const double h = -58;
   const double r = 500.;

   auto mooring0 = new MooringLine({r * cos(0.), r * sin(0.), h}, {0., 0., 0.}, total_length, n_points);
   auto mooring1 = new MooringLine({r * cos(120 * M_PI / 180), r * sin(120 * M_PI / 180), h}, {0., 0., 0.}, total_length, n_points);
   auto mooring2 = new MooringLine({r * cos(240 * M_PI / 180), r * sin(240 * M_PI / 180), h}, {0., 0., 0.}, total_length, n_points);

   std::map<MooringLine*, PVDWriter> pvd_line = {{mooring0, pvd_line0}, {mooring1, pvd_line1}, {mooring2, pvd_line2}};
   std::map<MooringLine*, PVDWriter> pvd_points = {{mooring0, pvd_points0}, {mooring1, pvd_points1}, {mooring2, pvd_points2}};

   mooring0->setDensityStiffnessDampingDiameter(density, stiffness, damp, diam);
   mooring1->setDensityStiffnessDampingDiameter(density, stiffness, damp, diam);
   mooring2->setDensityStiffnessDampingDiameter(density, stiffness, damp, diam);

   double t = 0;
   int output_index = 0;

   for (auto mooring : {mooring0, mooring1, mooring2})
      for (auto p : mooring->getPoints()) {
         p->velocity.fill(0);
         p->acceleration.fill(0);
      }

   auto boundary_condition = [&](networkPoint* p) {
      for (auto mooring : {mooring0, mooring1, mooring2}) {
         if (p == mooring->firstPoint) {
            p->acceleration.fill(0);
            p->velocity.fill(0);
         }
         if (p == mooring->lastPoint) {
            if (t < 4) {
               p->acceleration.fill(0);
               p->velocity.fill(0);
            } else {
               std::array<double, 3> additional_force = {0, -9.8 * p->mass * sin(2 * M_PI / 3 * t), 0.};
               auto > force = p->getForce() + p->mass * additional_force;
               p->acceleration[0] = force[0] / p->mass;
               p->acceleration[1] = force[1] / p->mass;
               p->acceleration[2] = force[2] / p->mass;
            }
         }
      }
   };
   /* -------------------------------------------------------------------------- */
   /*                         SIMULATION AND OUTPUT                              */
   /* -------------------------------------------------------------------------- */
   const int max_iteration = 10000;
   for (int interation = 0; interation < max_iteration; ++interation) {

      /* ---------------------------------------------- */
      /*                  OUTPUT                        */
      /* ---------------------------------------------- */
      if (interation % (int)(max_iteration / 1000.) == 0) {
         std::cout << interation << std::endl;
         int no = 0;
         for (auto mooring : {mooring0, mooring1, mooring2}) {
            std::unordered_map<networkPoint*, Tddd> point_force, point_tension;
            point_force.reserve(mooring->getPoints().size());
            point_tension.reserve(mooring->getPoints().size());
            for (auto p : mooring->getPoints()) {
               point_force[p] = p->forceTranslational();
               point_tension[p] = p->getTension();
            }
            //! line
            {
               auto filename = basePath + "/line_" + std::to_string(no) + "_" + std::to_string(output_index) + ".vtp";
               std::ofstream ofs(filename);
               vtkPolygonWrite(ofs, mooring->getLines());
               pvd_line.at(mooring).push(filename, t);
            }
            //! point
            {
               auto filename = basePath + "/point_" + std::to_string(no) + "_" + std::to_string(output_index) + ".vtp";
               std::ofstream ofs(filename);
               std::vector<std::tuple<std::string, decltype(point_force)>> data = {std::make_tuple("force", point_force),
                                                                                   std::make_tuple("tension", point_tension)};
               vtkPolygonWrite(ofs, mooring->getPoints(), data);
               pvd_points.at(mooring).push(filename, t);
            }
            no++;
         }
         output_index++;
         for (auto& [mooring, pvd_line] : pvd_line) pvd_line.output();
         for (auto& [mooring, pvd_points] : pvd_points) pvd_points.output();
      }

      /* ---------------------------------------------- */
      /*                 SIMULATION                     */
      /* ---------------------------------------------- */
      mooring0->simulate(t, dt, boundary_condition);
      mooring1->simulate(t, dt, boundary_condition);
      mooring2->simulate(t, dt, boundary_condition);
      t += dt;
   }
   /* -------------------------------------------------------------------------- */

   return 0;
}
