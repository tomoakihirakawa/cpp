#include "basic_IO.hpp"
#include "basic_linear_systems.hpp"
#include "basic_mathematical_functions.hpp"
#include "lib_measurement.hpp"
#include "minMaxOfFunctions.hpp"

VV_d A = {{0.0247911, 0.161413, 0.625419, 0.465341, 0.794249},
          {0.895294, 0.215363, 0.280354, 0.0206005, 0.906597},
          {0.457972, 0.76661, 0.590316, 0.535627, 0.00733951},
          {0.315392, 0.925959, 0.412796, 0.825637, 0.322538},
          {0.572894, 0.0998945, 0.738812, 0.30581, 0.904702}};

#define check_GMRES
#if defined(check_QR)
int main() {
   Timer timer;

   // VV_d A = {{12, -51, 4},
   //           {6, 167, -68},
   //           {-4, 24, -41}};

   // VV_d A = {{0.8147, 0.0975, 0.1576},
   //           {0.9058, 0.2785, 0.9706},
   //           {0.1270, 0.5469, 0.9572},
   //           {0.9134, 0.9575, 0.4854},
   //           {0.6324, 0.9649, 0.8003}};

   QR qr(A);
   MatrixForm(qr.Q, Blue);
   MatrixForm(qr.R, Red);
   MatrixForm(Dot(Transpose(qr.Q), qr.R), Magenta);

   std::cout << "time:" << timer() << std::endl;
};

#elif defined(check_GMRES)
int main() {

   V_d b = {1, 2, 3, 4, 5};

   // int s = 10;
   // VV_d A(s, V_d(s));
   // V_d b(s);
   // for (auto i = 0; i < s; ++i) {
   //    b[i] = RandomReal({-1., 1.});
   //    for (auto j = 0; j < s; ++j)
   //       A[i][j] = RandomReal({-1., 1.});
   // }
   //

   V_d x0(b.size(), 0.);
   Timer timer;
   std::cout << "time:" << timer() << std::endl;
   bool finished = false;
   // while (!finished) {
   //
   // auto v = diagonal_scaling_vector(A);
   // for (auto i = 0; i < v.size(); ++i) {
   //    A[i] *= v[i];
   //    b[i] *= v[i];
   // }
   Print(Norm(b));
   //
   for (auto i = 4; i < 11; i++) {
      gmres gm(A, b, x0, i);
      std::cout << "gm.x = " << gm.x << std::endl;
      std::cout << "time:" << timer() << std::endl;
      auto error = Norm(Dot(A, gm.x) - b);
      Print(i);
      if (error < 1E-10) {
         Print(error, Blue);
         finished = true;
         break;
      } else {
         Print(error, Green);
         x0 = gm.x;
      }
      std::cout << Red << "--------------------------------" << colorOff << std::endl;
   }
   // }
   std::cout << "time:" << timer() << std::endl;
   // std::cout << "gm->y.size() = " << gm.y.size() << std::endl;
   // std::cout << "gm->y = " << gm.y << std::endl;
   // std::cout << "V = " << gm.ap->V << std::endl;
   // std::cout << "H = " << gm.ap->H << std::endl;
   //    MatrixForm(gm.ap->V);
   //    MatrixForm(gm.ap->H);
};

#elif defined(check_GradientMethod)
int main() {
   VV_d A = {{4., 0., 0.},
             {6., 7., 0.},
             {7., 3., 5.}};

   // VV_d A = {{2., 1., 4.},
   // 		  {0., 3., 6.},
   // 		  {0., 0., 7.}};

   V_d b = {1., 2., 10.};
   V_d x(b.size(), 0);
   // A.x=bが与えられている場合
   GradientMethod gd(A);
   {
      TimeWatch tm;
      auto x = gd.solve(b);
      Print(x);
      Print(gd.count);
      std::cout << "time = " << tm.get() << std::endl;
   }

   {
      TimeWatch tm;
      auto x = gd.solveCG(b);
      Print(x);
      Print(gd.count);
      std::cout << "time = " << tm.get() << std::endl;
   }

   {
      TimeWatch tm;
      ludcmp lu(A);
      lu.solve(b, x);
      std::cout << x << std::endl;
      std::cout << "time = " << tm.get() << std::endl;
   }
   // 関数と関数の微分によって計算する
   // {
   // 	auto F = [](const V_d &x)
   // 	{
   // 		return x[0] * x[0] + x[1] * x[1];
   // 	};
   // 	auto dF = [](const V_d &x)
   // 	{
   // 		return V_d{2 * x[0], 2 * x[1]};
   // 	};
   // 	GradientMethod gd;
   // 	auto ans = gd.minimize(F, dF, {1, 1});
   // 	std::cout << ans << std::endl;
   // }

   std::cout << "back_substitution" << std::endl;
   std::cout << back_substitution(A, b) << std::endl;
   std::cout << "forward_substitution" << std::endl;
   std::cout << forward_substitution(A, b) << std::endl;
};
#endif