#define USE_BROYDEN_METHOD
//! sample_vectorの与え方に注意．成分毎に与えること．
std::array<double, 3> find_optimal_vector(std::vector<double> &Vsample,
                                          std::vector<Tddd> &Directions,
                                          const Tddd &Vinit,
                                          std::vector<double> weights) {

   if (Vsample.size() == 1)
      return Vsample[0] * Directions[0];

   const double threshold_angle_in_rad = M_PI / 180.;

   //! 与えられている情報が不十分な場合がある．

   /* -------------------------------------------------------------------------- */
   /*                         make directions into groups                        */
   /* -------------------------------------------------------------------------- */

   std::vector<std::vector<Tddd>> direction_groups;
   for (auto &d : Directions) {
      if (direction_groups.empty())
         direction_groups.push_back({d});
      else {
         bool is_new_direction = true;
         for (auto &dg : direction_groups) {
            if (std::ranges::any_of(dg, [&](const Tddd &dir) { return Between(Dot(dir, d), {0.9, 1.1}) && VectorAngle(dir, d) < threshold_angle_in_rad; })) {
               dg.push_back(d);
               is_new_direction = false;
               break;
            }
         }
         if (is_new_direction)
            direction_groups.push_back({d});
      }
   }

   if (direction_groups.size() == 1) {
      Tddd Dir0 = Normalize(Mean(direction_groups[0])), Dir1;

      double angleX = VectorAngle(Dir0, Tddd{1., 0., 0.});
      double angleY = VectorAngle(Dir0, Tddd{0., 1., 0.});
      double angleZ = VectorAngle(Dir0, Tddd{0., 0., 1.});
      double max_angle = angleX;

      Dir1 = {1., 0., 0.};
      if (angleY > max_angle) {
         max_angle = angleY;
         Dir1 = {0., 1., 0.};
      }
      if (angleZ > max_angle) {
         max_angle = angleZ;
         Dir1 = {0., 0., 1.};
      }
      double w = Mean(weights);
      Vsample.push_back(0.);
      Directions.push_back(Dir1);
      weights.push_back(w);
      //
      Vsample.push_back(0.);
      Directions.push_back(Normalize(Cross(Dir0, Dir1)));
      weights.push_back(w);
   } else if (direction_groups.size() == 2) {
      Vsample.push_back(0.);
      Directions.push_back(Cross(Mean(direction_groups[0]), Mean(direction_groups[1])));
      weights.push_back(Mean(weights));
   }

   for (auto &d : Directions)
      d = Normalize(d);

   /* -------------------------------------------------------------------------- */

   std::array<double, 3> U = Vinit;  //(Vsample[0] + Vsample[1]) / 2;

   auto diff = [&](const Tddd &U, const std::size_t i) -> double { return Dot(U, Directions[i]) - Vsample[i]; };

   auto optimizing_function = [&](const Tddd &U) -> double {
      double S = 0;
      for (std::size_t i = 0; i < Vsample.size(); ++i)
         S += weights[i] * std::pow(diff(U, i), 2) / 2.;
      return S;
   };

#ifdef USE_CONJUGATE_GRADIENT_METHOD
   double alpha = 1.;
   double alpha_min = 1.;
   double min;
   Tddd Umin, dSdu;
   const std::vector<double> line_search = {0., 1E-8, 1E-4, 1E-3, 1E-2, 0.1, 0.25, 0.5, 0.75, 1., 10., 1E+2, 1E+3, 1E+4};
   for (auto iteration = 0; iteration < 50; ++iteration) {
      double s = 0, a;
      dSdu.fill(0.);

      for (std::size_t i = 0; i < Vsample.size(); ++i)
         dSdu += 2. * Directions[i] * weights[i] * diff(U, i);

      Umin = U - line_search[0] * alpha * dSdu;
      min = optimizing_function(Umin);
      for (size_t i = 1; i < line_search.size(); ++i) {
         auto Utmp = U - line_search[i] * alpha * dSdu;
         double value = optimizing_function(Utmp);
         if (value < min) {
            alpha_min = line_search[i] * alpha;
            min = value;
            Umin = Utmp;
         }
      }

      U = Umin;
      alpha = std::clamp(alpha_min, 1E-8, 1E+4);
      // if (std::abs(s) < 1E-6) break;
   }
#elif defined(USE_BROYDEN_METHOD)
   BroydenMethod<Tddd> BM(U, U + Tddd{0., 0., 1E-10});

   auto grad_optimizing_function = [&](const Tddd &U) -> Tddd {
      Tddd grad = {0., 0., 0.};
      for (std::size_t i = 0; i < Vsample.size(); ++i)
         grad += weights[i] * Directions[i] * diff(U, i);
      return grad;
   };

   for (auto iteration = 0; iteration < 50; ++iteration) {
      BM.updateBFGS(grad_optimizing_function(U), grad_optimizing_function(U - BM.dX));
      U = BM.X;
      if (Norm(BM.dX) < 1e-10 && iteration > 5)
         break;
   }
#endif

   return U;
}

std::array<double, 3> find_optimal_vector(std::vector<double> &Vsample, std::vector<Tddd> &Directions, const Tddd &Vinit) {
   std::vector<double> weights(Vsample.size(), 1.);
   return find_optimal_vector(Vsample, Directions, Vinit, weights);
}
