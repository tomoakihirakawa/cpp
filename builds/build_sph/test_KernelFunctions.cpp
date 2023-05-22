/*DOC_EXTRACT SPH
## 核関数
3次スプライン関数と5次スプライン関数の実装とテストコード
* 関数の形状を確認．
* 体積積分が1になるかどうかを確認．
*/
#include "kernelFunctions.hpp"

int main() {
   double r = 2.;
   {
      std::ofstream output("data.txt");
      // Calculate and output w_Bspline3 and w_Bspline5 values for different x
      for (const auto &x : Subdivide({-r, r}, 25)) {
         double bspline3_value = w_Bspline3(std::abs(x), r);
         double bspline5_value = w_Bspline5(std::abs(x), r);
         output << x << " " << bspline3_value << " " << bspline5_value << std::endl;
      }
      output.close();
   }
   {
      std::ofstream output("Dot_grad_w_Bspline3_Dot.txt");
      // Calculate and output Dot_grad_w_Bspline3_Dot values for different x and y for gnuplot splot
      for (const auto &x : Subdivide({-r, r}, 25)) {
         for (const auto &y : Subdivide({-r, r}, 25)) {
            output << x << " " << y << " " << Dot_grad_w_Bspline3_Dot({0., 0., 0.}, {x, y, 0.}, r) << std::endl;
         }
         output << std::endl;  // empty line for gnuplot to understand different y-blocks
      }
      output.close();
   }
   {
      std::ofstream output("Dot_grad_w_Bspline5_Dot.txt");
      // Calculate and output Dot_grad_w_Bspline3_Dot values for different x and y for gnuplot splot
      for (const auto &x : Subdivide({-r, r}, 25)) {
         for (const auto &y : Subdivide({-r, r}, 25)) {
            output << x << " " << y << " " << Dot_grad_w_Bspline5_Dot({0., 0., 0.}, {x, y, 0.}, r) << std::endl;
         }
         output << std::endl;  // empty line for gnuplot to understand different y-blocks
      }
      output.close();
   }
   // Check integration
   auto w = std::setw(10);
   const std::array<double, 3> A = {0.0, 0.0, 0.0};  // center
   const double radius = 2.;

   // Calculate sum of w_Bspline3 and w_Bspline5 for different N values
   for (auto N : {5, 10, 15, 20, 25}) {

      double sum3 = 0, sum5 = 0;
      const double total_volume = std::pow(2 * radius, 3);
      const double each_volume = std::pow(2 * radius / N, 3);

      std::array<double, 3> X;
      // Iterate through x, y, z in the subdivided space
      for (const auto &x : Subdivide({-radius, radius}, N))
         for (const auto &y : Subdivide({-radius, radius}, N))
            for (const auto &z : Subdivide({-radius, radius}, N)) {
               // the volume of each small cube is (2 * radius / N)^3
               X = {x, y, z};
               sum3 += w_Bspline3(Norm(X - A), radius) * each_volume;
               sum5 += w_Bspline5(Norm(X - A), radius) * each_volume;
            }

      // Output the sums for each N value
      std::cout << "N = " << w << N << w << "sum3 = " << w << sum3 << w << "sum5 = " << w << sum5 << std::endl;
   }
   return 0;
}
