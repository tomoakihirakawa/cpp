/*DOC_EXTRACT NEWTON
## Bucketを用いた粒子探索のテスト
Smoothed Particle Hydrodynamics (SPH)では，効率的な近傍粒子探査が必要となる．
このコードでは，Bucketを用いた粒子探索のテストを行う．

結果はVTKファイルに出力される．
   * 全ての粒子を表示したものは`all.vtp`
   * 中心の粒子を表示したものは`center*.vtp`
   * 中心の粒子が探査したセル内にある粒子を表示したものは`inCell*.vtp`
   * セル内かつ球内にある粒子を表示したものは`inSphere*.vtp`

   - 各セルにある粒子を表示したものは`each_cell*.vtp`
   - 各セルの中心位置を表示したものは`each_cell_position*.vtp`
*/

// Include necessary libraries and header files
#define DEM
#include <filesystem>
#include "Network.hpp"
#include "vtkWriter.hpp"

// Function to test the bucket method
void test_Bucket(const auto &water,
                 const std::unordered_set<Network *> &nets,
                 const std::string &output_directory,
                 const auto &r) {
   std::filesystem::create_directory(output_directory);

   {
      vtkPolygonWriter<networkPoint *> vtp;
      for (const auto &net : nets)
         for (const auto &p : net->getPoints())
            vtp.add(p);
      std::ofstream ofs(output_directory + "all.vtp");
      vtp.write(ofs);
      ofs.close();
   }

   int i = 0;
   for (const auto &p : water->getPoints()) {
      {
         vtkPolygonWriter<networkPoint *> vtp;
         vtp.add(p);
         std::ofstream ofs(output_directory + "center" + std::to_string(i) + ".vtp");
         vtp.write(ofs);
         ofs.close();
      }
      {
         std::unordered_map<networkPoint *, double> Bspline3;
         std::unordered_map<networkPoint *, double> Bspline5;
         double sum3 = 0, sum5 = 0;
         vtkPolygonWriter<networkPoint *> vtp;
         int j = 0;
         for (const auto &net : nets)
            net->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &q) {
               vtp.add(q);
               sum3 += Bspline3[q] = w_Bspline3(Norm(q->X - p->X), p->radius_SPH) * p->volume;
               sum5 += Bspline5[q] = w_Bspline5(Norm(q->X - p->X), p->radius_SPH) * p->volume;
               j++;
            });
         vtp.addPointData("w_Bspline3", Bspline3);
         vtp.addPointData("w_Bspline5", Bspline5);

         std::cout << "number of points in cell: " << j << std::endl;
         std::ofstream ofs(output_directory + "inCell" + std::to_string(i) + ".vtp");
         vtp.write(ofs);
         ofs.close();
      }
      {
         std::unordered_map<networkPoint *, double> Bspline3;
         std::unordered_map<networkPoint *, double> Bspline5;
         double sum3 = 0, sum5 = 0;
         vtkPolygonWriter<networkPoint *> vtp;
         int j = 0;
         for (const auto &net : nets)
            net->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &q) {
               if (Distance(p->X, q->X) < p->radius_SPH) {
                  vtp.add(q);
                  sum3 += Bspline3[q] = w_Bspline3(Norm(q->X - p->X), p->radius_SPH) * p->volume;
                  sum5 += Bspline5[q] = w_Bspline5(Norm(q->X - p->X), p->radius_SPH) * p->volume;
                  j++;
               }
            });
         vtp.addPointData("w_Bspline3", Bspline3);
         vtp.addPointData("w_Bspline5", Bspline5);

         std::cout << "        points: " << j << std::endl;
         std::cout << "sum B-spline 3: " << sum3 << std::endl;
         std::cout << "sum B-spline 5: " << sum5 << std::endl;
         std::ofstream ofs(output_directory + "inSphere" + std::to_string(i) + ".vtp");
         vtp.write(ofs);
         ofs.close();
      }
      i++;
   }

   // check if each cell properly stores objects
   {
      int I = 0, J = 0;
      for (auto i = 0; i < water->BucketPoints.xsize; ++i)
         for (auto j = 0; j < water->BucketPoints.ysize; ++j)
            for (auto k = 0; k < water->BucketPoints.zsize; ++k) {
               {
                  vtkPolygonWriter<networkPoint *> vtp;
                  for (const auto &net : nets)
                     vtp.add(net->BucketPoints.buckets[i][j][k]);
                  std::ofstream ofs(output_directory + "each_cell" + std::to_string(I++) + ".vtp");
                  vtp.write(ofs);
                  ofs.close();
               }
               {
                  vtkPolygonWriter<Tddd> vtp;
                  for (const auto &net : nets)
                     vtp.add(net->BucketPoints.itox(i, j, k));
                  std::ofstream ofs(output_directory + "each_cell_position" + std::to_string(J++) + ".vtp");
                  vtp.write(ofs);
                  ofs.close();
               }
            }
   }
};

int main(int arg, char **argv) {

   for (const auto &d : Subdivide({-1., 1.}, 30))
      std::cout << std::setw(10) << d << " -> " << static_cast<int>(d) << std::endl;

   // Create a new network
   auto net = new Network;

   // Define x, y, z coordinates for the network points
   auto vecX = Subdivide({-.5, 1.}, 10);
   auto vecY = Subdivide({0., 2.}, 10);
   auto vecZ = Subdivide({-.5, 3.}, 10);

   // Define the kernel radius
   double r = 0.3;

   // Generate network points
   for (const auto &x : vecX)
      for (const auto &y : vecY)
         for (const auto &z : vecZ) {
            auto p = new networkPoint(net, {x, y, z});
            p->radius_SPH = r;
         }

   // Initialize the bucket points for the network
   net->makeBucketPoints(2 * r);

   // Call the test function to evaluate the bucket method and write the output to VTK files
   test_Bucket(net, {net}, "./test_Buckets_output/", r);
};