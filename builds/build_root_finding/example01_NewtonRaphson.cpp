#include "minMaxOfFunctions.hpp"
#include "rootFinding.hpp"

auto w = std::setw(20);

int main() {

   /*DOC_EXTRACT newton

    ## ロボットの節の位置をLightHillの曲線上に乗せる

    LightHillの式は以下のようになる．

    $$
    {\bf x}^{\rm LH}(x,t) = (x,y^{\rm LH}(x,t)),\quad
    y^{\rm LH}(x,t) = \left( \frac{c_1}{L} x + {c_2} \left(\frac{x}{L}\right)^2 \right) \sin \left( \frac{2 \pi}{L} x - \omega t \right)
    $$

    ここで，$c_1, c_2, L, \omega$は定数である．

    ロボットの$i$番目の節の位置は，${\bf x}_{i}^{\rm rb} = {\bf x}_{i-1}^{\rm rb} + r \left( \cos \theta_i, \sin \theta_i \right)$である．
    次の関数を使って表すことにする．

    $$
    {\bf x}_{i}^{\rm rb} = X^{\rm rb}({\bf x}_{i-1}^{\rm rb},r,\theta_i),
    \quad X^{\rm rb}({\bf a},r,\theta) = {\bf a} + r \left( \cos \theta, \sin \theta \right)
    $$

    ここで，$r$はロボットの節の長さ，$\theta_i$はロボットの節の角度である．
    頭の位置を$X_{0}^{\rm rb}=(0,0)$とする．
    次の節の位置は，$X_{1}^{\rm rb} = r (\cos \theta_1, \sin \theta_1)$である．
    さらに次の節の位置は，$X_{2}^{\rm rb} = X_{1}^{\rm rb} + r (\cos \theta_2, \sin \theta_2)$となる．

    目的関数$f$の微分は，

    $$
    \frac{df}{d\theta} = -r \sin\theta\frac{d y^{\rm LH} }{dx}-r\cos\theta
    $$

    NOTE: この目的関数$f$には，前の節の位置を与える必要がある．節の位置は，後ろの節の位置によって変わらないので，この目的関数を先頭から順番に最適化することは問題ない．

    */

   TimeWatch time;

   const double L = 0.71;
   const double w = 2. * M_PI * 1.0;
   const double k = 2. * M_PI * 2.0;
   const double c1 = 0.05;
   const double c2 = 0.05;
   const int n = 100;
   const double r = L / n;
   const double t = 0.;

   auto yLH = [&](const double x, const double t) { return (c1 * x / L + c2 * std::pow(x / L, 2)) * sin(k * (x / L) - w * t); };

   auto X_RB = [&](const std::array<double, 2> a, const double q) { return a + r * std::array<double, 2>{cos(q), sin(q)}; };

   auto f = [&](const std::array<double, 2> a, const double q, const double t) {
      auto [x, y] = X_RB(a, q);
      return yLH(x, t) - y;
   };

   auto ddx_yLH = [&](const double x, const double t) {
      return (c1 / L + 2 * c2 * x / std::pow(L, 2)) * sin(k * (x / L) - w * t) +
             (c1 / L * x + c2 * std::pow(x / L, 2)) * cos(k * (x / L) - w * t) * k / L;
   };

   auto ddq_f = [&](const double q, const double t) {
      double x = r * cos(q);
      return -r * sin(q) * ddx_yLH(x, t) - r * cos(q);
   };

   std::vector<double> Q(n, 0.);  // thetas
   std::array<double, 2> a{{0., 0.}};
   double error = 0;
   for (auto i = 1; i < Q.size(); i++) {
      NewtonRaphson nr(0.);
      error = 0;
      for (auto k = 0; k < 50; k++) {
         nr.update(f(a, nr.X, t), ddq_f(nr.X, t), 1.);
         if ((error = std::abs(f(a, nr.X, t))) < 1E-10)
            break;
      }
      Q[i] = nr.X;
      a = X_RB(a, Q[i]);
   }

   std::cout << time() << std::endl;

   {
      std::ofstream outFile("lighthill.txt");
      std::array<double, 2> a{{0., 0.}};
      outFile << a[0] << " " << a[1] << " " << yLH(a[0], t) << std::endl;
      for (auto i = 1; i < Q.size(); i++) {
         auto [x, y] = a = X_RB(a, Q[i]);
         outFile << x << " " << y << " " << yLH(x, t) << std::endl;
      }
      outFile.close();
   }
}
