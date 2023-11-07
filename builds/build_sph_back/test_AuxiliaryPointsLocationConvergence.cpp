#include "Network.hpp"
#include "kernelFunctions.hpp"

int main() {
   auto net = new Network("test_AuxiliaryPointsLocationConvergence");

   double radius = 1.;
   double CSML = 0.1;

   auto vecX = Subdivide({-5, 5}, 10);
   auto vecY = Subdivide({-5, 5}, 10);
   auto vecZ = Subdivide({-5, 5}, 10);

   std::array<double, 3> X;
   for (const auto& x : vecX)
      for (const auto& y : vecY)
         for (const auto& z : vecZ) {
            std::cout << x << " " << y << " " << z << std::endl;
            X = {x, y, z};
            auto p = new networkPoint(net, X);
            p->volume = 1.;
         }

   std::vector<std::array<double, 3>> ignore_Xs = {{0., 2., 2.}, {0., 1., 2.}, {0., 3., 2.}};
   std::array<double, 3> center = {0., 0., 0.};

   auto sum_W = [&]() {
      double sum = 0.;
      for (const auto& p : net->getPoints())
         if (std::ranges::none_of(ignore_Xs, [&](auto ignore_X) { return Norm(p->X - ignore_X) == 0.; }))
            sum += w_Bspline(Norm(p->X), 5.);
      return sum;
   };

   auto sum_grad_W = [&]() {
      std::array<double, 3> sum;
      sum.fill(0.);
      for (const auto& p : net->getPoints())
         if (std::ranges::none_of(ignore_Xs, [&](auto ignore_X) { return Norm(p->X - ignore_X) == 0.; }))
            sum += grad_w_Bspline(center, p->X, 5.);
      return sum;
   };

   std::array<double, 3> auxp_X = {0., 0., 0.};

   std::cout << sum_W() << std::endl;
   std::cout << sum_grad_W() << std::endl;

   std::ofstream datafile("test_AuxiliaryPointsLocationConvergence.txt");
   for (auto i = 0; i < 1000; ++i) {
      //   auto auxp_X = 5. * i / 1000. * Normalize(-sum_grad_W());
      auto auxp_X = 5. * i / 1000. * Normalize(-sum_grad_W());

      auto grad_f = sum_grad_W() + grad_w_Bspline(center, auxp_X, 5.);
      auto f = sum_W() + w_Bspline(Norm(center - auxp_X), 5.) - 1;

      datafile << Norm(auxp_X) << " " << Dot(grad_f, grad_f) + f * f << std::endl;
   }
   //
   datafile.close();

   //    for (auto i = 0; i < 10000; ++i) {
   //       auto grad = sum_grad_W() + grad_w_Bspline(center, auxp_X, 5.);
   //       auxp_X -= (grad) / 2.;
   //       std::cout << auxp_X << std::endl;
   //       std::cout << sum_W() + w_Bspline(Norm(center - auxp_X), 5.) << std::endl;
   //       std::cout << grad << std::endl;
   //    }
};