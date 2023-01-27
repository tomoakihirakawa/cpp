#ifndef SPH_weightingFunctions_H
#define SPH_weightingFunctions_H

#include "kernelFunctions.hpp"
#include "vtkWriter.hpp"

// 値にかかる重みをベクトルとして出力する関数
// 値を未知変数とする，離散化された方程式を作成する際には，未知変数の重みを抜き出す必要がある
// その際に便利である
V_d W_Bspline(const std::vector<Tddd> &X, V_d dX /*出力する*/, const Tddd &a, const double h, const int order = 3) {
   if (order == 3) {
      for (auto i = 0; i < X.size(); i++)
         dX[i] = kernel_Bspline3(X[i], a, h) * dX[i];
   } else {
      for (auto i = 0; i < X.size(); i++)
         dX[i] = kernel_Bspline5(X[i], a, h) * dX[i];
   }
   return dX;
};
/* ------------------------------------------------------ */
/*                     grad 勾配を計算する際は，             */
/* ------------------------------------------------------ */
std::vector<Tddd> W_grad_Bspine(std::vector<Tddd> X /*出力する*/, const V_d &dX, const Tddd &a, const double h, const int order = 3) {
   if (order == 3) {
      for (auto i = 0; i < X.size(); i++) {
         if (Norm(X[i] - a) > 1E-10) {
            X[i] = grad_kernel_Bspline3(X[i], a, h) * dX[i];
            // X[i] = grad_kernel_Bspline3(X[i], a, e) * dX[i];
         } else
            X[i] = {0., 0., 0.};
      }
   } else {
      for (auto i = 0; i < X.size(); i++) {
         if (Norm(X[i] - a) > 1E-10) {
            X[i] = grad_kernel_Bspline5(X[i], a, h) * dX[i];
            // X[i] = grad_kernel_Bspline3(X[i], a, e) * dX[i];
         } else
            X[i] = {0., 0., 0.};
      }
   }
   return X;
};
std::vector<Tddd> W_grad_Bspine(const std::unordered_set<networkPoint *> &ps,
                                const networkPoint *const a,
                                const double h,
                                const int order = 3) {
   std::vector<Tddd> ret(ps.size());
   int i = 0;
   if (order == 3) {
      for (const auto &p : ps) {
         if (Norm(ToX(p) - ToX(a)) > 1E-10)
            ret[i++] = grad_kernel_Bspline3(ToX(p), ToX(a), h) * p->volume;
         else
            ret[i++] = {0., 0., 0.};
      }
   } else {
      for (const auto &p : ps) {
         if (Norm(ToX(p) - ToX(a)) > 1E-10)
            ret[i++] = grad_kernel_Bspline5(ToX(p), ToX(a), h) * p->volume;
         else
            ret[i++] = {0., 0., 0.};
      }
   }
   return ret;
};
/* ------------------- 仮想マーカーにおける流速を算出 ------------------ */
std::tuple<Tddd, Tddd, double, double> U_tmpU_mass_density_SPH_IDW(const std::unordered_set<networkPoint *> &ps,
                                                                   const Tddd &X,
                                                                   const double power = 2) {
   Tddd tmp_U = {0, 0, 0}, U = {0, 0, 0};
   double mass = 0., density = 0., total = 0, w = 0.;
   for (const auto &p : ps) {
      w = std::pow(1. / Norm(ToX(p) - X), power);
      total += w;
      U += w * p->U_SPH;
      tmp_U += w * p->tmp_U_SPH;
      mass += w * p->mass;
      density += w * p->rho;
   }
   return {U / total,
           tmp_U / total,
           mass / total,
           density / total};
};

double dummy_pressure_Asai(const Buckets<networkPoint *> &B_water,
                           const networkPoint *dummy_p /*dummy point*/,
                           const double power = 2) {
   // ダミー粒子のみ使える関数
   auto [f, t0, t1, d, dx] = dummy_p->particlize_info;
   auto Y = oppositeX(dummy_p->particlize_info);
   auto X = ToX(dummy_p);
   auto p_at_Y = B_water.getObjects_unorderedset(Y, 8 /*深さ最大*/, 1);
   // double sign = Normalize(Dot(f->getXtuple() - Y, f->normal) * f->normal); // Y->F (+ or -)
   Tddd gravity = {0., 0., -9.81};
   double pressure = 0, total = 0, norm, nu, tmp, pn;
   Tddd accel = {0, 0, 0};
   if (!p_at_Y.empty()) {
      for (const auto &q : p_at_Y) {
         // ここを修正
         // 重み付き平均した反対位置OXと壁との距離を確認し
         // OX->
         nu = q->mu_SPH / q->rho;
         pn = q->rho * Dot(nu * q->lap_U + gravity - accel, f->normal);
         tmp = q->p_SPH + /*修正*/ pn * Dot(X - q->getXtuple(), f->normal);
         // tmp = (q->p_SPH);
         norm = Norm(X - q->getXtuple());
         if (norm < 1E-15)
            return tmp;
         auto w = std::pow(1. / norm, power);
         total += w;
         pressure += w * tmp;
      }
      pressure /= total;
      return pressure;
   }
   return pressure;
};
//
Tddd velocity_SPH_IDW(const std::unordered_set<networkPoint *> &ps,
                      const Tddd &X,
                      const double h,
                      const int order = 3) {
   Tddd ret = {0, 0, 0};
   double total = 0, norm;
   for (const auto &p : ps) {
      norm = Norm(ToX(p) - X);
      ret += 1. / norm * p->U_SPH;
      total += 1. / norm;
   }
   return ret / total;
};

Tddd velocity_SPH(const std::unordered_set<networkPoint *> &ps,
                  const Tddd &X,
                  const double h,
                  const int order = 3) {
   //* これはシンプルに勾配計算の改良版である．
   Tddd ret = {0, 0, 0};
   for (const auto &p : ps)
      ret += p->U_SPH * p->volume * kernel_Bspline5(ToX(p), X, h);
   return ret;
};
double density_SPH(const std::unordered_set<networkPoint *> &ps,
                   const networkPoint *const a,
                   const double h,
                   const int order = 3) {
   //* これはシンプルに勾配計算の改良版である．
   double ret = a->rho * a->volume * kernel_Bspline5(ToX(a), ToX(a), h);
   for (const auto &p : ps)
      if (p != a)
         ret += p->rho * p->volume * kernel_Bspline5(ToX(p), ToX(a), h);
   return ret;
};

constexpr auto grad_kernel_Bspline = grad_kernel_Bspline5;

Tddd normal(const std::unordered_set<networkPoint *> &ps, const networkPoint *const a, const double h) {
   Tddd ret = {0, 0, 0};
   int count = 0;
   for (const auto &p : ps)
      if (p != a) {
         ret += grad_kernel_Bspline(ToX(p), ToX(a), h);
         count++;
      }
   if (count == 0)
      return {0., 0., -1.};
   else
      return Normalize(ret);
};

Tddd normal(const std::vector<Tddd> &Xs, const networkPoint *const a, const double h) {
   Tddd ret = {0, 0, 0};
   int count = 0;
   for (const auto &x : Xs) {
      if (Between(Norm(x - ToX(a)), {1E-10, h})) {
         ret += grad_kernel_Bspline(x, ToX(a), h);
         count++;
      }
   }
   if (count == 0)
      return {0., 0., -1.};
   else
      return Normalize(ret);
};

Tddd cg_neighboring_particles(const std::vector<Tddd> &Xs, const networkPoint *a) {
   Tddd ret = {0, 0, 0};
   int count = 0;
   for (const auto &x : Xs)
      if (Norm(x - ToX(a)) <= a->radius_SPH && Norm(x - ToX(a)) > 1E-12) {
         ret += x;
         count++;
      }
   return ret / count;
};

// b! -------------------------------------------------------------------------- */
// b! -------------------------------------------------------------------------- */
// b! -------------------------------------------------------------------------- */
struct IntersectionsSphereTrianglesLines {
   /*
   メモ：
   @ まず，面と対面する場合にのみ，鏡のように反射させれば良いのではないか，と考えるが，面が折れ曲がっている場合，
   @ 線に対しても鏡写しにしなければ，空白領域が生まれてしまうことに気づく．
   @
   @ 三角形の位置関係から，角を抽出する．これが複雑で，計算コストがかかりすぎる．
   @ 角は，線分と考える．線分に対して反射させる．これは，SPHの計算で利用する．
   */
   std::unordered_set<networkPoint *> points;
   std::vector<networkFace *> faces;
   std::vector<std::tuple<networkFace *, networkFace *, T2Tddd>> lines;
   /* -------------------------------------------------------------------------- */
   IntersectionsSphereTrianglesLines(const std::unordered_set<networkFace *> &contact_faces)
       : faces(contact_faces.begin(), contact_faces.end()), lines(extractFacesAndLines(contact_faces)){};
   /* -------------------------------------------------------------------------- */
   std::vector<std::tuple<networkFace *, networkFace *, T2Tddd>> extractFacesAndLines(const std::unordered_set<networkFace *> &contact_faces) {
      std::vector<std::tuple<networkFace *, networkFace *, T2Tddd>> ret;
      // ベクトル化
      T3Tddd fi, fj;
      double minangle = 0.1;
      for (auto i = 0; i < faces.size(); ++i) {
         fi = faces[i]->getXVertices();
         fi += 0.01 * (T3Tddd{std::get<0>(fi) - Mean(fi), std::get<1>(fi) - Mean(fi), std::get<2>(fi) - Mean(fi)});
         for (auto j = i; j < faces.size(); ++j) {
            if (isFlat(faces[j]->normal, faces[i]->normal, minangle / 180 * M_PI)) {
               fj = faces[j]->getXVertices();
               fj += 0.01 * (T3Tddd{std::get<0>(fj) - Mean(fj), std::get<1>(fj) - Mean(fj), std::get<2>(fj) - Mean(fj)});
               auto intxn = IntersectionTriangles(fi, fj);
               if (intxn.isIntersecting) {
                  auto V1 = std::get<1>(intxn.L) - std::get<0>(intxn.L);
                  if (Norm(V1) > 0.1 * (Norm(std::get<0>(fj) - Mean(fj)) + Norm(std::get<1>(fj) - Mean(fj)) + Norm(std::get<2>(fj) - Mean(fj))) / 3.)
                     ret.push_back({faces[i], faces[j], intxn.L});
               }
            }
         }
      }
      return ret;
   };
   /* -------------------------------------------------------------------------- */
   std::vector<std::tuple<Tddd, Tddd>> X_Y;
   std::vector<std::tuple<networkFace *, networkFace *, Tddd, Tddd, Tddd>> F_F_X_Y_N;
   /* -------------------------------------------------------------------------- */
   bool isDuplication(const Tddd &n) const {
      return std::any_of(std::begin(this->F_F_X_Y_N), std::end(this->F_F_X_Y_N), [n](const auto &f_f_x_y_n) {
         return isFlat(n, std::get<0>(f_f_x_y_n)->normal, M_PI * 180. * 1E-2);
      });
   };
   /* -------------------------------------------------------------------------- */
   std::vector<std::tuple<networkFace *, networkFace *, Tddd, Tddd, Tddd>> getFFXYN(const networkPoint *p, const double h) {
      /*
      面があるのは，面が作るオブジェクトの速度などを参照するため．
      */
      //! isDuplicationによって重複がないことを確かめて面を取得する．より近い面に切り替えるよう修正が必要．
      /*
      外側
      _________
               |
               |
               _________このような場合に，重複する面の内遠い方を参照する可能性がある
               流体側
         X
      */
      F_F_X_Y_N.clear();
      // b! 次にgetでは，面または交線が，与えられた点と正対する場合にその位置を返す．
      //! 面との干渉
      Tddd n;
      for (const auto &f : this->faces) {
         auto intxn = IntersectionSphereTriangleLimitedToNormalRegion(ToX(p), h, f->getXVertices());
         n = f->normal;
         if (intxn.isIntersecting && !isDuplication(n))
            F_F_X_Y_N.push_back({f, nullptr, intxn.X, ToX(p) + 2. * (intxn.X - ToX(p)), n});
      }
      //! 線との干渉
      for (const auto &[F0, F1, L] : this->lines) {
         auto intxn = IntersectionSphereLine(ToX(p), h, L);
         n = Normalize(F0->normal + F1->normal);
         if (intxn.isIntersecting && !isDuplication(n))
            F_F_X_Y_N.push_back({F0, F1, intxn.X, ToX(p) + 2. * (intxn.X - ToX(p)), n});
      }
      return F_F_X_Y_N;
   };
};
// b! -------------------------------------------------------------------------- */
// b! -------------------------------------------------------------------------- */
// b! -------------------------------------------------------------------------- */
#define use_space_potential_particle

Tddd getVectorToSPP(const networkPoint *const p) {
   // return Norm(p->cg_neighboring_particles_SPH - ToX(p)) * p->interpolated_normal_SPH;
   return 2. * ToX(p) - p->cg_neighboring_particles_SPH;
};
Tddd Reflect(const Tddd &v, const Tddd &n) {
   return v - 2. * n * Dot(v, n);
};

#define free_slip_boundary_condition
Tddd mirroredVelocity(const Tddd &U, const Tddd &n, const Tddd &wall_velocity) {
// #if defined(free_slip_boundary_condition)
//     return U - 2 * n * Dot(U, n) + Dot(wall_velocity, n) * n;
// #elif defined(no_slip_boundary_condition)
//     return -U + Dot(wall_velocity, n) * n;
// #elif defined(zero_boundary_condition)
//     return Dot(wall_velocity, n) * n;
// #else
//     return -n * Dot(U, n) + Dot(wall_velocity, n) * n;
// #endif

// 外挿的な考え方
#if defined(free_slip_boundary_condition)
   return U - 2. * n * Dot(U, n) + 2. * Dot(wall_velocity, n) * n;
#elif defined(no_slip_boundary_condition)
   return -U + 2. * Dot(wall_velocity, n) * n;
#elif defined(zero_boundary_condition)
   return 2. * Dot(wall_velocity, n) * n;
#else
   return -n * Dot(U, n) + 2. * Dot(wall_velocity, n) * n;
#endif
};

//@ ------------------------------------------------------ */
//@                         速度の発散                       */
//@                          div(U)                        */
//@ ------------------------------------------------------ */
double div_U_polygon_boundary(const std::unordered_set<networkPoint *> &ps,
                              const networkPoint *const a,
                              const double h,
                              const double dt) {
   auto INTXN = IntersectionsSphereTrianglesLines(a->getContactFaces());
   double ret = 0, tmp = 0;
   Tddd n, U, velocity, SSP_X, mirroedSSP_X, toSSP;
   for (const auto &p : ps) {
      if (p != a)
         ret += p->mass * Dot(p->U_SPH - a->U_SPH, grad_kernel_Bspline(ToX(p), ToX(a), h));
#if defined(use_space_potential_particle)
      if (p->isSurface && p == a) {
         // 水面上部
         tmp = p->mass * Dot(p->U_SPH - a->U_SPH, grad_kernel_Bspline(ToX(p) + getVectorToSPP(p), ToX(a), h));
         if (isFinite(tmp))
            ret += tmp;
      }
#endif
      // b$ ------------------------ ポリゴン ------------------------ */
      for (const auto &[F0, F1, X, Y, N] : INTXN.getFFXYN(p, h)) {
         if (F1) {
            n = N;
            velocity = F0->normal * Dot(ToTddd(F0->getNetwork()->velocity), F0->normal) +
                       F1->normal * Dot(ToTddd(F1->getNetwork()->velocity), F1->normal);
         } else {
            n = N;
            velocity = F0->normal * Dot(ToTddd(F0->getNetwork()->velocity), F0->normal);
         }
         U = mirroredVelocity(p->U_SPH, n, velocity);
#if defined(use_space_potential_particle)
         if (p->isSurface) {
            tmp = p->mass * Dot(U - a->U_SPH, grad_kernel_Bspline(Y + Reflect(getVectorToSPP(p), n), ToX(a), h));
            if (isFinite(tmp))
               ret += tmp;
         }
#endif
         ret += p->mass * Dot(U - a->U_SPH, grad_kernel_Bspline(Y, ToX(a), h));
      }
      // b%$------------------------------------------------------ */
   }
   return ret / a->rho;
};

double div_U(const std::unordered_set<networkPoint *> &ps,
             const networkPoint *const a,
             const double h) {
   //* これはシンプルに勾配計算の改良版である．
   double ret = 0;
   for (const auto &p : ps)
      if (p != a)
         ret += p->mass * Dot(p->U_SPH - a->U_SPH, grad_kernel_Bspline(ToX(p), ToX(a), h));
   return ret / a->rho;
};
//
double div_tmp_U_polygon_boundary(const std::unordered_set<networkPoint *> &ps,
                                  const networkPoint *const a,
                                  const double h,
                                  const double dt) {
   auto INTXN = IntersectionsSphereTrianglesLines(a->getContactFaces());
   //* これはシンプルに勾配計算の改良版である．
   double ret = 0, tmp = 0;
   Tddd U, n, velocity;
   for (const auto &p : ps) {
      if (p != a) {
         ret += p->mass * Dot(p->tmp_U_SPH - a->tmp_U_SPH, grad_kernel_Bspline(ToX(a), ToX(p), h));
      }
#if defined(use_space_potential_particle)
      if (p->isSurface && p == a) {
         tmp = p->mass * Dot(p->tmp_U_SPH - a->tmp_U_SPH, grad_kernel_Bspline(ToX(a),
                                                                              ToX(p) + getVectorToSPP(p),
                                                                              h));
         if (isFinite(tmp))
            ret += tmp;
      }
#endif
      // b$ ------------------------ ポリゴン ------------------------ */
      for (const auto &[F0, F1, X, Y, N] : INTXN.getFFXYN(p, h)) {
         if (F1) {
            n = N;
            velocity = F0->normal * Dot(ToTddd(F0->getNetwork()->velocity), F0->normal) +
                       F1->normal * Dot(ToTddd(F1->getNetwork()->velocity), F1->normal);
         } else {
            n = N;
            velocity = F0->normal * Dot(ToTddd(F0->getNetwork()->velocity), F0->normal);
         }
         U = mirroredVelocity(p->U_SPH, n, velocity);
#if defined(use_space_potential_particle)
         if (p->isSurface) {
            tmp = p->mass * Dot(U - a->tmp_U_SPH, grad_kernel_Bspline(ToX(a),
                                                                      Y + Reflect(getVectorToSPP(p), n),
                                                                      h));
            if (isFinite(tmp))
               ret += tmp;
         }
#endif
         ret += p->mass * Dot(U - a->tmp_U_SPH, grad_kernel_Bspline(ToX(a), Y, h));
      }
      // b%$------------------------------------------------------ */
   }

   return ret / a->rho;
};
//
double div_tmp_U(const std::unordered_set<networkPoint *> &ps, const networkPoint *const p, const double h, const double dt) {
   /*
   depend on tmp_U_SPH
   */
   double ret = 0;
   for (const auto &q : ps)
      if (q != p)
         ret += q->mass * Dot(q->tmp_U_SPH - p->tmp_U_SPH, grad_kernel_Bspline(p->X, q->X, h));
   return ret / p->rho;
};
//
double div_tmp_U(const std::unordered_set<networkPoint *> &ps, const networkPoint *const p, const double h, const int order = 3) {
   // double ret = 0;
   // for (const auto &q : ps)
   //     if (q != p)
   //         ret += q->volume *
   //                Dot(q->tmp_U_SPH, grad_kernel_Bspline(ToX(q), ToX(p), h));
   // return ret;
   /* ------------------------------------------------------ */
   //* これはシンプルな勾配計算の改良版である．
   double ret = 0;
   for (const auto &q : ps)
      if (q != p)
         ret += q->mass * Dot(q->tmp_U_SPH - p->tmp_U_SPH, grad_kernel_Bspline(ToX(p), ToX(q), h));
   return ret / p->rho;
};
// # ------------------------------------------------------ */
// #                          圧力勾配                       */
// #                         grad(P)                        *
// # ------------------------------------------------------ */
Tddd grad_P_Monaghan1992_polygon_boundary_(const std::unordered_set<networkPoint *> &ps /*流体だけを渡す*/,
                                           const networkPoint *const a,
                                           const double h) {
   // 密度の演算子の作用下に入れ込んだMonaghan1992の方法
   Tddd ret = {0, 0, 0}, tmp;
   Tddd accel = {0, 0, 0};
   Tddd gravity = {0., 0., -9.81}, n;
   auto INTXN = IntersectionsSphereTrianglesLines(a->getContactFaces());
   double pressure, density2;
   double a_pressure_density2 = a->p_SPH_ / std::pow(a->rho, 2);
   for (const auto &p : ps) {
      density2 = std::pow(p->rho, 2);
      if (p != a)
         ret += p->mass * (p->p_SPH_ / density2 + a_pressure_density2) * grad_kernel_Bspline(ToX(p), ToX(a), h);
#if defined(use_space_potential_particle)
      if (p->isSurface && p == a) {
         tmp = p->mass * (0. /*SPPの圧力は0*/ / density2 + a_pressure_density2) * grad_kernel_Bspline(ToX(p) + getVectorToSPP(p), ToX(a), h);
         if (isFinite(tmp))
            ret += tmp;
      }
#endif

      // b$ ------------------------ ポリゴン ------------------------ */
      for (const auto &[F0, F1, X, Y, N] : INTXN.getFFXYN(p, h)) {
         if (F1)
            accel = F0->normal * Dot(ToTddd(F0->getNetwork()->acceleration), F0->normal) +
                    F1->normal * Dot(ToTddd(F1->getNetwork()->acceleration), F1->normal);
         else
            accel = F0->normal * Dot(ToTddd(F0->getNetwork()->acceleration), F0->normal);
         pressure = (p->p_SPH_ + /*修正*/ p->rho * Dot(p->mu_SPH / p->rho * p->lap_U + gravity - accel, Y - ToX(p)));
         ret += p->mass * (pressure / density2 + a_pressure_density2) * grad_kernel_Bspline(Y, ToX(a), h);
#if defined(use_space_potential_particle)
         if (p->isSurface) {
            tmp = p->mass * (0 / density2 + a_pressure_density2) * grad_kernel_Bspline(Y + Reflect(getVectorToSPP(p), N), ToX(a), h);
            if (isFinite(tmp))
               ret += tmp;
         }
#endif
      }
      // b%$------------------------------------------------------ */
   }
   return ret * a->rho;
};

Tddd grad_P_Monaghan1992(const std::unordered_set<networkPoint *> &ps, const networkPoint *const i, const double h) {
   // 密度の演算子の作用下に入れ込んだMonaghan1992の方法
   Tddd ret = {0, 0, 0};
   for (const auto &j : ps)
      if (j != i)
         ret += j->mass * (j->p_SPH / std::pow(j->rho, 2) + i->p_SPH / std::pow(i->rho, 2)) * grad_kernel_Bspline(ToX(j), ToX(i), h);
   return ret * i->rho;
};

Tddd grad_P_Monaghan1992_(const std::unordered_set<networkPoint *> &ps, const networkPoint *const a, const double h) {
   // 密度の演算子の作用下に入れ込んだMonaghan1992の方法
   Tddd ret = {0, 0, 0};
   for (const auto &p : ps)
      if (p != a)
         ret += p->mass *
                (p->p_SPH_ / std::pow(p->rho, 2) + a->p_SPH_ / std::pow(a->rho, 2)) *
                grad_kernel_Bspline(ToX(p), ToX(a), h);
   return ret * a->rho;
};

//% ------------------------------------------------------ */
//%                          粘性項                         */
//%                       laplacian(U)                     */
//% ------------------------------------------------------ */
Tddd laplacian_U_Monaghan1992(const networkPoint *const p,
                              const Tddd &X,
                              const Tddd &U,
                              const networkPoint *const a,
                              const double h,
                              const double sigma,
                              const double alpha = 0.1,
                              const double beta = 0.1) {
   auto Uij = a->U_SPH - U;
   auto rij = ToX(a) - X;
   auto dot_rij_rij = Dot(rij, rij);
   auto dot_Uij_rij = Dot(Uij, rij);
   double Cs = 1466., PIij, div_U, m;
   if (dot_Uij_rij < 0) {
      m = h / sigma * dot_Uij_rij / dot_rij_rij /*hを消せばUの発散となっている*/;
      PIij = (-alpha * Cs * m + beta * m * m) / ((p->rho + a->rho) / 2.);
      return p->mass * PIij * grad_kernel_Bspline(X, ToX(a), h);
   }
   return {0, 0, 0};
};
Tddd laplacian_U_Monaghan1992_polygon_boundary(const std::unordered_set<networkPoint *> &ps,
                                               const networkPoint *const a,
                                               const double h,
                                               const double sigma,
                                               const double dt,
                                               const double alpha = 0.1,
                                               const double beta = 0.1) {
   /*
   kernel_Bspline5の定義はwij=w(|p-a|)としている．
   */
   // 人工粘性
   // まず，nu*laplacian(U)を計算する
   Tddd ret = {0, 0, 0}, tmp, n, U, velocity;
   auto INTXN = IntersectionsSphereTrianglesLines(a->getContactFaces());

   for (const auto &p : ps) {
      if (p != a)
         ret += laplacian_U_Monaghan1992(p, ToX(p), p->U_SPH, a, h, sigma, alpha, beta);
#if defined(use_space_potential_particle)
      if (p->isSurface && p == a) {
         tmp = laplacian_U_Monaghan1992(p, ToX(p) + getVectorToSPP(p), p->U_SPH, a, h, sigma, alpha, beta);
         if (isFinite(tmp))
            ret += tmp;
      }
#endif
      // b$ ------------------------ ポリゴン ------------------------ */
      for (const auto &[F0, F1, X, Y, N] : INTXN.getFFXYN(p, h)) {
         if (F1) {
            n = N;
            velocity = F0->normal * Dot(ToTddd(F0->getNetwork()->velocity), F0->normal) +
                       F1->normal * Dot(ToTddd(F1->getNetwork()->velocity), F1->normal);
         } else {
            n = N;
            velocity = F0->normal * Dot(ToTddd(F0->getNetwork()->velocity), F0->normal);
         }
         U = mirroredVelocity(p->U_SPH, n, velocity);
#if defined(use_space_potential_particle)
         if (p->isSurface) {
            tmp = laplacian_U_Monaghan1992(p, Y + Reflect(getVectorToSPP(p), n), U, a, h, sigma, alpha, beta);
            if (isFinite(tmp))
               ret += tmp;
         }

#endif
         ret += laplacian_U_Monaghan1992(p, Y, U, a, h, sigma, alpha, beta);
      }
      // b$ ------------------------------------------------------ */
   }
   double nu = a->mu_SPH / a->rho;
   return -ret / nu;
};

Tddd laplacian_U_Monaghan1992(const std::unordered_set<networkPoint *> &ps,
                              const networkPoint *const a,
                              const double h,
                              const double sigma,
                              const double alpha = 0.1,
                              const double beta = 0.1) {
   /*
   kernel_Bspline5の定義はwij=w(|p-a|)としている．
   */
   // 人工粘性
   // まず，nu*laplacian(U)を計算する
   Tddd ret = {0, 0, 0}, Uij, rij, grad;
   double PIij, div_U, dot_rij_rij, dot_Uij_rij, m;
   for (const auto &p : ps)
      if (p != a)
         ret += laplacian_U_Monaghan1992(p, ToX(p), p->U_SPH, a, h, alpha, beta);
   double nu = a->mu_SPH / a->rho;
   return -ret / nu;
};
/* -------------------------------------------------------------------------- */
Tddd viscousTerm_ShaoAndLo2003(const networkPoint *a, const networkPoint *p, double h) {
   auto Uij = a->U_SPH - p->U_SPH;
   auto Xij = ToX(a) - ToX(p);
   auto nu = p->mu_SPH / p->rho + a->mu_SPH / a->rho;
   return (p->mass * 8 * nu * Dot(Uij, Xij) * grad_kernel_Bspline(ToX(a), ToX(p), h)) / ((p->rho + a->rho) * Dot(Xij, Xij) + std::pow(1E-4 * h, 2));
};
Tddd viscousTerm_ShaoAndLo2003(const networkPoint *a, const networkPoint *p, const double &h, const Tddd &X) {
   auto Uij = a->U_SPH - p->U_SPH;
   auto Xij = ToX(a) - X;
   auto nu = p->mu_SPH / p->rho + a->mu_SPH / a->rho;
   return (p->mass * 8 * nu * Dot(Uij, Xij) * grad_kernel_Bspline(ToX(a), X, h)) / ((p->rho + a->rho) * Dot(Xij, Xij) + std::pow(1E-4 * h, 2));
};
Tddd viscousTerm_ShaoAndLo2003(const networkPoint *a, const networkPoint *p, double h, const Tddd &X, const Tddd &U) {
   auto Uij = a->U_SPH - U;
   auto Xij = ToX(a) - X;
   auto nu = p->mu_SPH / p->rho + a->mu_SPH / a->rho;
   return (p->mass * 8 * nu * Dot(Uij, Xij) * grad_kernel_Bspline(ToX(a), X, h)) / ((p->rho + a->rho) * Dot(Xij, Xij) + std::pow(1E-4 * h, 2));
};

Tddd laplacian_U_ShaoAndLo2003_polygon_boundary(const std::unordered_set<networkPoint *> &ps,
                                                const networkPoint *const a,
                                                const double h) {
   // nu*laplacian(U)を計算する
   Tddd ret = {0, 0, 0}, tmp = {0, 0, 0}, U, n, velocity;
   double nu_i, nu_j, norm;
   auto INTXN = IntersectionsSphereTrianglesLines(a->getContactFaces());
   //
   std::vector<Tddd> crt_dist;
   for (const auto &[F0, F1, X, Y, N] : INTXN.getFFXYN(a, h))
      if (Norm(X) < h / 6)
         crt_dist.emplace_back(X);

   //
   for (const auto &p : ps) {
      if (p != a)
         ret += viscousTerm_ShaoAndLo2003(a, p, h);
#if defined(use_space_potential_particle)
      if (p->isSurface && p == a) {
         tmp = viscousTerm_ShaoAndLo2003(a, p, h, ToX(p) + getVectorToSPP(p));
         if (isFinite(tmp))
            ret += tmp;
      }
#endif
      // b$ ------------------------ ポリゴン ------------------------ */
      for (const auto &[F0, F1, X, Y, N] : INTXN.getFFXYN(p, h)) {
         if (F1) {
            n = N;
            velocity = F0->normal * Dot(ToTddd(F0->getNetwork()->velocity), F0->normal) +
                       F1->normal * Dot(ToTddd(F1->getNetwork()->velocity), F1->normal);
         } else {
            n = N;
            velocity = F0->normal * Dot(ToTddd(F0->getNetwork()->velocity), F0->normal);
         }
         U = mirroredVelocity(p->U_SPH, n, velocity);
#if defined(use_space_potential_particle)
         if (p->isSurface) {
            tmp = viscousTerm_ShaoAndLo2003(a, p, h, Y + Reflect(getVectorToSPP(p), n), U);
            if (isFinite(tmp))
               ret += tmp;
         }
#endif
         ret += viscousTerm_ShaoAndLo2003(a, p, h, Y, U);
      }
      // b$ ------------------------------------------------------ */
   }
   return ret / (a->mu_SPH / a->rho);
};

Tddd laplacian_U_ShaoAndLo2003(const std::unordered_set<networkPoint *> &target,
                               const networkPoint *const a, const double h) {
   /* depend on U_SPH */
   Tddd ret = {0, 0, 0};
   for (const auto &q : target)
      if (q != a)
         ret += viscousTerm_ShaoAndLo2003(a, q, h);
   return ret / (a->mu_SPH / a->rho);
};

/* ------------------------------------------------------ */
double pressure_EISPH_Hosseini2007_polygon_boundary(const std::unordered_set<networkPoint *> &ps,
                                                    const networkPoint *const i,
                                                    const double h,
                                                    const double dt) {
#ifndef use_space_potential_particle
   if (i->isSurface)
      return 0;
#endif
   if (ps.size() < 2)
      return 0;
   if (dt < 1E-13)
      return i->p_SPH;
   // [1] Nomeritae, E. Daly, S. Grimaldi, and H. H. Bui, “Explicit incompressible SPH algorithm for free-surface flow modelling: A comparison with weakly compressible schemes,” Adv. Water Resour., vol. 97, pp. 156–167, Nov. 2016.
   Tddd Xij, accel, gravity = {0, 0, -9.81}, n, U, velocity;
   double tmp, sum_Aij = 0, sum_Aij_Pj = 0, Aij, pressure;
   auto INTXN = IntersectionsSphereTrianglesLines(i->getContactFaces());
   for (const auto &j : ps) {
      if (j != i) {
         Xij = ToX(i) - ToX(j);
         tmp = j->mass * 8. / std::pow(j->rho + i->rho, 2);
         Aij = tmp * Dot(Xij, grad_kernel_Bspline(ToX(i), ToX(j), h)) / Dot(Xij, Xij);
         sum_Aij += Aij;
         sum_Aij_Pj += Aij * j->p_SPH;
      }
#if defined(use_space_potential_particle)
      if (j->isSurface && j == i) {
         pressure = 0.;
         auto Y = ToX(j) + getVectorToSPP(j);
         Xij = ToX(i) - Y;
         tmp = j->mass * 8. / std::pow(j->rho + i->rho, 2);
         Aij = tmp * Dot(Xij, grad_kernel_Bspline(ToX(i), Y, h)) / Dot(Xij, Xij);
         if (isFinite(Aij)) {
            sum_Aij += Aij;
            sum_Aij_Pj += Aij * pressure;
         }
      }
#endif
      // b$ ------------------------ ポリゴン ------------------------ */
      for (const auto &[F0, F1, X, Y, N] : INTXN.getFFXYN(j, h)) {
         if (F1) {
            n = N;
            velocity = F0->normal * Dot(ToTddd(F0->getNetwork()->velocity), F0->normal) +
                       F1->normal * Dot(ToTddd(F1->getNetwork()->velocity), F1->normal);
            accel = F0->normal * Dot(ToTddd(F0->getNetwork()->acceleration), F0->normal) +
                    F1->normal * Dot(ToTddd(F1->getNetwork()->acceleration), F1->normal);
         } else {
            n = N;
            velocity = F0->normal * Dot(ToTddd(F0->getNetwork()->velocity), F0->normal);
            accel = F0->normal * Dot(ToTddd(F0->getNetwork()->acceleration), F0->normal);
         }

         U = mirroredVelocity(j->U_SPH, n, velocity);
         pressure = (j->p_SPH + /*Asaiの修正*/ j->rho * Dot(j->mu_SPH / j->rho * j->lap_U + gravity - accel, Y - ToX(j)));
         Xij = ToX(i) - Y;
         tmp = j->mass * 8. / std::pow(j->rho + i->rho, 2);
         Aij = tmp * Dot(Xij, grad_kernel_Bspline(ToX(i), Y, h)) / Dot(Xij, Xij);
         sum_Aij += Aij;
         sum_Aij_Pj += Aij * pressure;
#if defined(use_space_potential_particle)
         if (j->isSurface) {
            pressure = 0;
            // pressure = (j->p_SPH + /*Asaiの修正*/ j->rho * Dot(j->mu_SPH / j->rho * j->lap_U + gravity - accel, Y - ToX(j)));
            Xij = ToX(i) - (Y + Reflect(getVectorToSPP(j), n));
            tmp = j->mass * 8. / std::pow(j->rho + i->rho, 2);
            Aij = tmp * Dot(Xij, grad_kernel_Bspline(ToX(i), Y + Reflect(getVectorToSPP(j), n), h)) / Dot(Xij, Xij);
            if (isFinite(Aij)) {
               sum_Aij += Aij;
               sum_Aij_Pj += Aij * pressure;
            }
         }
#endif
      }
      // b$ ------------------------------------------------------ */
   }
   double Bi = i->div_U / dt;
   tmp = (Bi + sum_Aij_Pj) / sum_Aij;
   if (isFinite(tmp))
      return tmp;
   else
      return 0.;  // この場合は，周辺の粒子がゼロで，sum_Aijが0ということなので，圧力はゼロで構わない．
};

double pressure_EISPH_Hosseini2007(const std::unordered_set<networkPoint *> &ps,
                                   const networkPoint *const i, const double h, const double dt) {
   /*
   depend on div_U
   */
   // [1] Nomeritae, E. Daly, S. Grimaldi, and H. H. Bui, “Explicit incompressible SPH algorithm for free-surface flow modelling: A comparison with weakly compressible schemes,” Adv. Water Resour., vol. 97, pp. 156–167, Nov. 2016.
   if (dt < 1E-40)
      return i->p_SPH;
   else if (ps.empty())
      return 0.;
   else {
      Tddd Xij;
      double Aij, sum_Aij = 0, sum_Aij_Pj = 0;
      for (const auto &j : ps)
         if (j != i) {
            Xij = ToX(i) - ToX(j);
            Aij = (j->mass * 8. / std::pow(j->rho + i->rho, 2)) * Dot(Xij, grad_kernel_Bspline(ToX(i), ToX(j), h)) / Dot(Xij, Xij);
            sum_Aij += Aij;
            sum_Aij_Pj += Aij * j->p_SPH;
         }
      return (i->div_U / dt + sum_Aij_Pj) / sum_Aij;
   }
};
/* ------------------------------------------------------ */
double DPDt_EISPH_Hosseini2007_polygon_boundary(const std::unordered_set<networkPoint *> &ps,
                                                const networkPoint *const i,
                                                const double h,
                                                const double dt) {
   if (dt < 1E-40)
      return 0.;
   else if (i->isSurface)
      return -i->p_SPH / dt;
   else {
      auto Pn1 = pressure_EISPH_Hosseini2007_polygon_boundary(ps, i, h, dt);
      return (Pn1 - i->p_SPH) / dt;
   }
};
/* ------------------------------------------------------ */

Tddd laplacian_Brookshaw1985(const std::unordered_set<networkPoint *> &ps,
                             networkPoint *a,
                             const double e,
                             const int order = 3) {
   double tot = 0., tmp = 0.;
   Tddd ret = {0, 0, 0};
   Tddd r;
   for (const auto &p : ps)
      if (p != a) {
         r = ToX(p) - ToX(a);
         tmp = 2. * Dot(r / std::pow(Norm(r), 2), p->volume * grad_kernel_Bspline(ToX(p), ToX(a), e));
         ret -= tmp * p->U_SPH;
         tot += tmp;
      }
   return ret + tot * a->U_SPH;
   // eq.(91) D. J. Price, “Smoothed particle hydrodynamics and magnetohydrodynamics,” J. Comput. Phys., vol. 231, no. 3, pp. 759–794, 2012.
};
// b# -------------------------------------------------------------------------- */
// b# -------------------------------------------------------------------------- */
// b# -------------------------------------------------------------------------- */
void calculateDerivativesByEISPH(const auto &water_points, const auto &net, const double &dt) {
   /*フラクショナルステップ法を使って時間積分する（Cummins1999）．*/
   Timer watch;
   // b! ------------------------------------------------------ */
   // b!                    ∇.∇UとU*を計算                       */
   // b! ------------------------------------------------------ */
   DebugPrint("粘性項の∇.∇Uを計算し，次にU*を計算", Magenta);
#pragma omp parallel
   for (const auto &p : water_points)
#pragma omp single nowait
   {
      p->lap_U = laplacian_U_ShaoAndLo2003_polygon_boundary(p->getContactPoints(net), p, p->radius_SPH);
      p->tmp_U_SPH = p->U_SPH + ((p->mu_SPH / p->rho) * p->lap_U + _GRAVITY3_) * dt;
   }
   std::cout << green << "Elapsed time: " << Red << watch() << colorOff << " s\n";
   DebugPrint("U*を使って，仮想的な位置X*へ粒子を移動", Green);
   for (const auto &p : water_points)
      p->setX(ToX(p) + p->tmp_U_SPH * dt);
   // b! ------------------------------------------------------ */
   // b!               div(U), DρDt=-ρdiv(U)の計算               */
   // b! ------------------------------------------------------ */
   DebugPrint("div(U)を計算", Magenta);
#pragma omp parallel
   for (const auto &p : water_points)
#pragma omp single nowait
   {
      p->div_U = div_tmp_U_polygon_boundary(p->getContactPoints(net), p, p->radius_SPH, dt);
      p->setDensity(p->rho + (-p->rho * p->div_U) * dt);
   }
   std::cout << Green << "Elapsed time: " << Red << watch() << colorOff << " s\n";
   // DebugPrint("EISPH: 1.1 ダミー粒子（鏡映関係を使う）の圧力を計算済みの状態で，流体粒子（pressure_EISPH_Hosseini2007）", Red);
   // EISPHの圧力計算は初期の圧力を使うため，ルンゲクッタステップ毎に圧力を初期値に戻す必要がある．
   for (const auto &p : water_points)
      p->p_SPH = p->RK_P.getXinit();
      // b! ------------------------------------------------------ */
      // b!  　　　　　          仮の圧力Pの計算                　　　           　*/
      // b! ------------------------------------------------------ */
#pragma omp parallel
   for (const auto &p : water_points)
#pragma omp single nowait
   {
      p->p_SPH_ = pressure_EISPH_Hosseini2007_polygon_boundary(p->getContactPoints(net), p, p->radius_SPH, dt);
   }
   // b! ------------------------------------------------------ */
   // b!                　　圧力勾配 grad(P)の計算                   */
   // b! ------------------------------------------------------ */
   DebugPrint("圧力勾配∇Pを計算", Magenta);
#pragma omp parallel
   for (const auto &p : water_points)
#pragma omp single nowait
   {
      p->gradP_SPH = grad_P_Monaghan1992_polygon_boundary_(p->getContactPoints(net), p, p->radius_SPH);
   }
   std::cout << Green << "Elapsed time: " << Red << watch() << colorOff << " s\n";
   // b! ------------------------------------------------------ */
   // b!                         DU/Dtの計算                     */
   // b! ------------------------------------------------------ */
   DebugPrint("DU/Dtの計算の計算", Magenta);
#pragma omp parallel
   for (const auto &p : water_points)
#pragma omp single nowait
   {
      p->DUDt_SPH = -p->gradP_SPH / p->rho + (p->mu_SPH / p->rho) * p->lap_U + _GRAVITY3_;
   }
   std::cout << Green << "Elapsed time: " << Red << watch() << colorOff << " s\n";
};
/* -------------------------------------------------------------------------- */
void setDataOmitted(auto &vtp, const auto &Fluid) {
   std::unordered_map<networkPoint *, Tddd> normal_SPH;
   for (const auto &p : Fluid->getPoints())
      normal_SPH[p] = p->normal_SPH;
   vtp.addPointData("normal_SPH", normal_SPH);
   std::unordered_map<networkPoint *, Tddd> interpolated_normal_SPH;
   for (const auto &p : Fluid->getPoints())
      interpolated_normal_SPH[p] = p->interpolated_normal_SPH;
   vtp.addPointData("interpolated_normal_SPH", interpolated_normal_SPH);
   std::unordered_map<networkPoint *, Tddd> U;
   for (const auto &p : Fluid->getPoints())
      U[p] = p->U_SPH;
   vtp.addPointData("U", U);
   std::unordered_map<networkPoint *, double> pressure;
   for (const auto &p : Fluid->getPoints())
      pressure[p] = p->p_SPH;
   vtp.addPointData("pressure", pressure);
};

void setData(auto &vtp, const auto &Fluid) {
   std::unordered_map<networkPoint *, Tddd> U;
   for (const auto &p : Fluid->getPoints())
      U[p] = p->U_SPH;
   vtp.addPointData("U", U);
   //
   std::unordered_map<networkPoint *, double> density;
   for (const auto &p : Fluid->getPoints())
      density[p] = p->rho;
   vtp.addPointData("density", density);
   //
   std::unordered_map<networkPoint *, double> pressure;
   for (const auto &p : Fluid->getPoints())
      pressure[p] = p->p_SPH;
   vtp.addPointData("pressure", pressure);
   //
   // std::unordered_map<networkPoint *, double> contactpoints;
   // for (const auto &p : Fluid->getPoints())
   //    contactpoints[p] = (double)p->getContactPoints().size();
   // vtp.addPointData("contact points", contactpoints);
   // //
   // std::unordered_map<networkPoint *, Tddd> interpolated_normal_SPH;
   // for (const auto &p : Fluid->getPoints())
   //    interpolated_normal_SPH[p] = p->interpolated_normal_SPH;
   // vtp.addPointData("interpolated_normal_SPH", interpolated_normal_SPH);
   // //
   // std::unordered_map<networkPoint *, Tddd> lap_U;
   // for (const auto &p : Fluid->getPoints())
   //    lap_U[p] = p->lap_U;
   // vtp.addPointData("lap_U", lap_U);
   // //
   // std::unordered_map<networkPoint *, Tddd> tmp_U_SPH;
   // for (const auto &p : Fluid->getPoints())
   //    tmp_U_SPH[p] = p->tmp_U_SPH;
   // vtp.addPointData("tmp_U_SPH", tmp_U_SPH);
   // //
   // std::unordered_map<networkPoint *, double> div_U;
   // for (const auto &p : Fluid->getPoints())
   //    div_U[p] = p->div_U;
   // vtp.addPointData("div_U", div_U);
   // //
   // std::unordered_map<networkPoint *, double> DPDt_SPH;
   // for (const auto &p : Fluid->getPoints())
   //    DPDt_SPH[p] = p->DPDt_SPH;
   // vtp.addPointData("DPDt_SPH", DPDt_SPH);
   // //
   // std::unordered_map<networkPoint *, Tddd> normal_SPH;
   // for (const auto &p : Fluid->getPoints())
   //    normal_SPH[p] = p->normal_SPH;
   // vtp.addPointData("normal_SPH", normal_SPH);
   // //
   // std::unordered_map<networkPoint *, Tddd> whereToReference;
   // for (const auto &p : Fluid->getPoints())
   //    whereToReference[p] = ToX(p) + 2 * p->normal_SPH - ToX(p);
   // vtp.addPointData("where to reference", whereToReference);
   // //
   std::unordered_map<networkPoint *, double> isSurface;
   for (const auto &p : Fluid->getPoints())
      isSurface[p] = p->isSurface ? 1. : 0.;
   vtp.addPointData("isSurface", isSurface);
   // //
   std::unordered_map<networkPoint *, double> isInsideOfBody;
   for (const auto &p : Fluid->getPoints())
      isInsideOfBody[p] = p->isInsideOfBody;
   vtp.addPointData("isInsideOfBody", isInsideOfBody);
   // //
   // std::unordered_map<networkPoint *, double> isCaptured;
   // for (const auto &p : Fluid->getPoints())
   //    isCaptured[p] = p->isCaptured ? 1. : -1E+50;
   // vtp.addPointData("isCaptured", isCaptured);
   // //
   // std::unordered_map<networkPoint *, Tddd> repulsive_force_SPH;
   // for (const auto &p : Fluid->getPoints())
   //    repulsive_force_SPH[p] = p->repulsive_force_SPH;
   // vtp.addPointData("repulsive_force_SPH", repulsive_force_SPH);
   // //
   // std::unordered_map<networkPoint *, double> contact_points_SPH;
   // for (const auto &p : Fluid->getPoints())
   //    contact_points_SPH[p] = p->contact_points_SPH;
   // vtp.addPointData("contact_points_SPH", contact_points_SPH);
};

// #define Morikawa2021
// #define Morikawa2021_IDW
#define OurWall_LoopNumber 4

// b# -------------------------------------------------------------------------- */
// b# -------------------------------------------------------------------------- */
// b# -------------------------------------------------------------------------- */
void setRigidBodyObject_Pressure(auto &net, const std::unordered_set<networkPoint *> &wall_p, const auto &RigidBodyObject) {
   Timer watch;
   for (const auto &p : wall_p)
      p->p_SPH_ = p->p_SPH = 0;
   const double a = 0.75;
#if defined(Morikawa2021) || defined(Morikawa2021_IDW)
   for (auto i = 0; i < 1; ++i)
#else
   for (auto i = 0; i < OurWall_LoopNumber; ++i)
#endif
   {
#pragma omp parallel
      for (const auto &p_wall : wall_p)
#pragma omp single nowait
      {
         p_wall->p_SPH_ = 0.;
         Tddd Xi = p_wall->X + 2 * p_wall->normal_SPH, accel;
         double weight, nu, add, sum_weight = 0.;
         int count = 0;
         net->BucketPoints.apply(Xi, p_wall->radius_SPH, [&](const auto &q) {
            if (Dot(Normalize(q->X - Xi), Normalize(p_wall->normal_SPH)) > cos(100 * M_PI / 180.) /*この条件は妥当だろうか*/) {
#if defined(Morikawa2021_IDW)
               weight = std::pow(Norm(q->X - Xi), -2);
#else
               weight = q->volume * w_Bspline(Norm(q->X - Xi), p_wall->radius_SPH);
#endif
               if (Norm(q->X - p_wall->X) < Norm(p_wall->normal_SPH))
                  p_wall->p_SPH_ += q->p_SPH * weight;
               else {
                  accel = {0., 0., 0.};
                  nu = q->mu_SPH / q->rho;
                  add = q->rho * Dot(nu * q->lap_U + _GRAVITY3_ - accel, -2 * p_wall->normal_SPH);
                  p_wall->p_SPH_ += (q->p_SPH + add) * weight;
               }
               sum_weight += weight;
               count++;
            }
         });
#if defined(Morikawa2021_IDW)
         if (count == 0)
            p_wall->p_SPH_ = 0.;
         else
            p_wall->p_SPH_ /= sum_weight;
#elif defined(Morikawa2021)
         if (isFinite(p_wall->p_SPH_ / sum_weight))
            p_wall->p_SPH_ /= sum_weight;
#else
         // if (count <= 5)
         //    p_wall->p_SPH_ = 0;
         for (const auto &[obj, poly] : RigidBodyObject) {
            obj->BucketPoints.apply(Xi, p_wall->radius_SPH, [&](const auto &q) {
               p_wall->p_SPH_ += q->p_SPH * q->volume * w_Bspline(Norm(q->X - Xi), p_wall->radius_SPH);
            });
         }
#endif
      }
      // for (const auto &[obj, poly] : RigidBodyObject)
      for (const auto &p_wall : wall_p)
         if (i == 0)
            p_wall->p_SPH = p_wall->p_SPH_;
         else
            p_wall->p_SPH = p_wall->p_SPH_ * a + p_wall->p_SPH * (1 - a);
   }

   std::stringstream ss;
#if defined(Morikawa2021)
   ss << Blue << "setRigidBodyObject_Pressure (Morikawa2021) Elapsed time: " << Red << watch();
#elif defined(Morikawa2021_IDW)
   ss << Blue << "setRigidBodyObject_Pressure (Morikawa2021 with IDW) Elapsed time: " << Red << watch();
#else
   ss << Blue << "setRigidBodyObject_Pressure (another method) Elapsed time: " << Red << watch();
#endif
   DebugPrint(ss.str());
};
// b! -------------------------------------------------------------------------- */
void setRigidBodyObject_tmpU_and_U_SPH(auto &net, const std::unordered_set<networkPoint *> &wall_p, const auto &RigidBodyObject) {
   /*
   UとU*の両方を更新する．
   計算上，同時に更新することは問題ない．
   */
   Timer watch;
   for (const auto &p : wall_p)
      p->U_SPH = p->U_SPH_ = {0., 0., 0.};
   const double a = 0.75;
#if defined(Morikawa2021) || defined(Morikawa2021_IDW)
   for (auto i = 0; i < 1; ++i)
#else
   for (auto i = 0; i < OurWall_LoopNumber; ++i)
#endif
   {

#pragma omp parallel
      for (const auto &p_wall : wall_p)
#pragma omp single nowait
      {
         p_wall->U_SPH_ = p_wall->tmp_U_SPH_ = {0., 0., 0.};
         p_wall->rho_ = 0.;
         double sum_weight = 0.;
         auto Xi = p_wall->X + 2 * p_wall->normal_SPH;
         double weight, l = Norm(p_wall->normal_SPH);
         Tddd U, n = Normalize(p_wall->normal_SPH);
         int count = 0;
         {
            auto func = [&](const auto &q) {
               if (Dot(Normalize(q->X - Xi), Normalize(p_wall->normal_SPH)) > cos(100 * M_PI / 180.) /*この条件は妥当だろうか*/) {
#if defined(Morikawa2021_IDW)
                  weight = std::pow(Norm(q->X - Xi), -2);
#else
                  weight = q->volume * w_Bspline(Norm(q->X - Xi), p_wall->radius_SPH);
#endif
                  // free-slip
                  if (Norm(q->X - Xi) <= p_wall->radius_SPH) {
                     p_wall->U_SPH_ += (q->U_SPH - 2 * Dot(q->U_SPH, n) * n) * weight;
                     p_wall->tmp_U_SPH_ += (q->tmp_U_SPH - 2 * Dot(q->tmp_U_SPH, n) * n) * weight;
                     sum_weight += weight;
                     count++;
                  }
                  //
                  // no-slip
                  // if (Norm(q->X - Xi) <= p_wall->radius_SPH) {
                  //    p_wall->U_SPH_ += -q->U_SPH * weight;
                  //    p_wall->tmp_U_SPH_ += -q->tmp_U_SPH * weight;
                  //    sum_weight += weight;
                  //    count++;
                  // }
                  p_wall->rho_ += q->rho * weight;
               }
            };
            net->BucketPoints.apply(Xi, p_wall->radius_SPH, func);
         }
#if defined(Morikawa2021_IDW)
         if (count == 0) {
            p_wall->U_SPH_ = {0., 0., 0.};
            p_wall->tmp_U_SPH_ = {0., 0., 0.};
         } else {
            p_wall->U_SPH_ /= sum_weight;
            p_wall->tmp_U_SPH_ /= sum_weight;
         }
#elif defined(Morikawa2021)
         if (isFinite(p_wall->U_SPH_ / sum_weight, 1E+3)) {
            p_wall->U_SPH_ /= sum_weight;
            p_wall->tmp_U_SPH_ /= sum_weight;
         }
#else
         /*
         周囲の壁粒子の流速は反転などせずにそのまま読み込む．
         */
         if (count <= 5) {
            p_wall->U_SPH_ = {0., 0., 0.};
            p_wall->tmp_U_SPH_ = {0., 0., 0.};
         }
         for (const auto &[obj, poly] : RigidBodyObject) {
            obj->BucketPoints.apply(Xi, p_wall->radius_SPH, [&](const auto &q) {
               weight = q->volume * w_Bspline(Norm(q->X - Xi), p_wall->radius_SPH);
               p_wall->U_SPH_ += q->U_SPH * weight;
               p_wall->tmp_U_SPH_ += q->tmp_U_SPH * weight;
               p_wall->rho_ += q->rho_ * weight;
            });
         }
#endif
      }

      for (const auto &p_wall : wall_p) {
         if (i == 0) {
            p_wall->U_SPH = p_wall->U_SPH_;
            p_wall->tmp_U_SPH = p_wall->tmp_U_SPH_;
            // p_wall->rho = p_wall->rho_;
         } else {
            p_wall->U_SPH = p_wall->U_SPH_ * a + p_wall->U_SPH * (1 - a);
            p_wall->tmp_U_SPH = p_wall->tmp_U_SPH_ * a + p_wall->tmp_U_SPH * (1 - a);
            // p_wall->rho = p_wall->rho * a + p_wall->rho_ * (1 - a);
         }
      }
   }
   std::stringstream ss;
#if defined(Morikawa2021)
   ss << Blue << "setRigidBodyObject_U_SPH (Morikawa2021) Elapsed time: " << Red << watch();
#elif defined(Morikawa2021_IDW)
   ss << Blue << "setRigidBodyObject_U_SPH (Morikawa2021 with IDW) Elapsed time: " << Red << watch();
#else
   ss << Blue << "setRigidBodyObject_U_SPH (another method) Elapsed time: " << Red << watch();
#endif
   DebugPrint(ss.str());
};
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
void developByEISPH(Network *net, const std::vector<std::tuple<Network *, Network *>> &RigidBodyObject,
                    double &real_time, const double C_SML, const double particle_spacing) {
   try {
      Timer watch;
      //% --------------------------- 平滑化距離の計算 ------------------------------ */
      /*     密度, 平滑化距離      */
      std::cout << Green << "固定の平滑化距離の計算: C_SML * particle_spacing = " << C_SML << " * " << particle_spacing << " = " << C_SML * particle_spacing << colorOff << std::endl;
      for (const auto &[obj, poly] : RigidBodyObject)
         for (const auto &p : obj->getPoints()) {
            p->setDensity(1000.);
            p->radius_SPH = C_SML * particle_spacing;
            p->isFreeFalling = false;
         }
      for (const auto &p : net->getPoints()) {
         p->setDensity(1000.);
         p->radius_SPH = C_SML * particle_spacing;
         p->isFreeFalling = false;
         p->isInsideOfBody = false;
      }
      //% -------------- バケットの生成, p->radius_SPHの範囲だけ点を取得 --------------- */
      DebugPrint("バケットの生成", Green);
      net->makeBucketPoints(particle_spacing);
      for (const auto &[obj, poly] : RigidBodyObject)
         obj->makeBucketPoints(particle_spacing);
      std::cout << Green << "Elapsed time: " << Red << watch() << Magenta << " バケットの生成" << colorOff << " s\n";
      //% --------------- CFL条件を満たすようにタイムステップ間隔dtを設定 ------------------ */
      DebugPrint("近傍粒子探査が終わったら時間ステップを決めることができる", Green);
      // double dt = 0.002;
      // auto C_CFL_velocity = 0.07;  // dt = C_CFL_velocity*h/Max(U)
      // auto C_CFL_accel = 0.2;      // dt = C_CFL_accel*sqrt(h/Max(A))

      double dt = 0.0025;
      auto C_CFL_velocity = 0.07;  // dt = C_CFL_velocity*h/Max(U)
      auto C_CFL_accel = 0.2;      // dt = C_CFL_accel*sqrt(h/Max(A))

      for (const auto &p : net->getPoints()) {
         // 速度に関するCFL条件
         double max_dt = C_CFL_velocity * (p->radius_SPH / C_SML) / Norm(p->U_SPH);
         if (dt > max_dt && isFinite(max_dt))
            dt = max_dt;
         // 加速度に関するCFL条件m
         max_dt = C_CFL_accel * std::sqrt((p->radius_SPH / C_SML) / Norm(p->DUDt_SPH));
         if (dt > max_dt && isFinite(max_dt))
            dt = max_dt;
      }
      std::cout << "dt = " << dt << std::endl;
      //@ ----------------------- ルンゲクッタの準備 ------------------- */
      int RK_order = 2;
      for (const auto &p : net->getPoints()) {
         p->RK_U.initialize(dt, real_time, p->U_SPH, RK_order);
         p->RK_X.initialize(dt, real_time, p->X, RK_order);
         p->RK_P.initialize(dt, real_time, p->p_SPH, RK_order);
      }
      // b@ ======================================================= */
      // b@                  ルンゲクッタを使った時間積分                  */
      // b@ ======================================================= */
      std::unordered_set<networkPoint *> wall_p;
      do {
         wall_p.clear();
         dt = (*net->getPoints().begin())->RK_X.getdt();
         std::cout << "dt = " << dt << std::endl;
         //
         double C = 1.1;
         //
         //@ ------------------------------------------------------- */
         //@                    関連する壁粒子をマーク                   */
         //@ ------------------------------------------------------- */
         DebugPrint("関連する壁粒子をマーク", Green);
         for (const auto &[obj, poly] : RigidBodyObject)
            for (const auto &p : obj->getPoints())
               p->isCaptured = false;
         for (const auto &[obj, poly] : RigidBodyObject) {
#pragma omp parallel
            for (const auto &p : net->getPoints())
#pragma omp single nowait
            {
               obj->BucketPoints.apply(p->X, p->radius_SPH * C, [&](const auto &q) { q->isCaptured = true; });
            }
         };
         for (const auto &[obj, poly] : RigidBodyObject)
            for (const auto &p : obj->getPoints())
               if (p->isCaptured)
                  wall_p.emplace(p);
         std::cout << Green << "Elapsed time: " << Red << watch() << "s " << Magenta << "関連する壁粒子をマークし，保存" << colorOff << std::endl;
         //@ ------------------------------------------------------ */
         //@             流体粒子の法線方向の計算，水面の判定             */
         //@ ------------------------------------------------------ */
         //@  A. Krimi, M. Jandaghian, and A. Shakibaeinia, Water (Switzerland), vol. 12, no. 11, pp. 1–37, 2020.
         DebugPrint("水粒子のオブジェクト外向き法線方向を計算", Green);
#pragma omp parallel
         for (const auto &p : net->getPoints())
#pragma omp single nowait
         {
            /* ---------------------- p->interpolated_normal_SPHの計算 --------------------- */
            p->interpolated_normal_SPH = {0., 0., 0.};
            net->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &q) {
               if (Between(Distance(p, q), {1E-12, p->radius_SPH}))
                  p->interpolated_normal_SPH -= grad_w_Bspline(p->X, q->X, p->radius_SPH);
            });
            p->interpolated_normal_SPH = Normalize(p->interpolated_normal_SPH);
            if (!isFinite(p->interpolated_normal_SPH))
               p->interpolated_normal_SPH = {0., 0., 1.};
            /* ----------------------------------- 検索 ----------------------------------- */
            double r = (p->radius_SPH / C_SML) * 3.;  // チェック範囲 : particle spacing * 3
            bool found = false;
            int c = 0;
            auto search = [&](const auto &q) {
               if (Between(Distance(p, q), {1E-10, r})) {
                  c++;
                  if (found || (isFlat(p->interpolated_normal_SPH, q->X - p->X, M_PI / 4.))) {
                     found = true;
                  }
               }
            };
            net->BucketPoints.apply(p->X, r, search);
            int count = c;
            for (const auto &[obj, poly] : RigidBodyObject)
               obj->BucketPoints.apply(p->X, r, search);
            p->isSurface = !found || count < 5;
         }
         //
         DebugPrint("壁粒子のオブジェクト外向き法線方向を計算", Green);
         for (const auto &[obj, poly] : RigidBodyObject)
            for (const auto &p : obj->getPoints())
               p->interpolated_normal_SPH = {0, 0, 0};
#pragma omp parallel
         for (const auto &p : wall_p) {
#pragma omp single nowait
            {
               p->interpolated_normal_SPH = {0., 0., 0.};
               for (const auto &[obj, poly] : RigidBodyObject)
                  obj->BucketPoints.apply(p->X, p->radius_SPH * C, [&](const auto &q) {
                     if (Between(Distance(p, q), {1E-12, p->radius_SPH * C}))
                        p->interpolated_normal_SPH -= grad_w_Bspline(p->X, q->X, p->radius_SPH * C);
                  });
               p->interpolated_normal_SPH = Normalize(p->interpolated_normal_SPH);
               if (!isFinite(p->interpolated_normal_SPH))
                  p->interpolated_normal_SPH = {0., 0., 1.};
            }
         }
         // DebugPrint("normal_SPHを計算", Green);
         // #pragma omp parallel
         //          for (const auto &p : wall_p)
         // #pragma omp single nowait
         //          {
         //             for (const auto &[obj, poly] : RigidBodyObject) {
         //                // T2Tddd AB = {5. * p->radius_SPH * C * p->interpolated_normal_SPH + p->X, p->X};

         //                auto [XonF, f] = Nearest_(p->X, poly->getFaces());
         //                p->normal_SPH = XonF - p->X;
         //                p->mirroring_face = f;

         // if (p->mirroring_face) {
         //    networkPoint *closest_p_from_wall = nullptr;
         //    double distance = 1E+20;
         //    obj->BucketPoints.apply(XonF, p->radius_SPH, [&](const auto &q) {
         //       if (distance > Norm(q->X - XonF)) {
         //          distance = Norm(q->X - XonF);
         //          closest_p_from_wall = q;
         //       }
         //    });
         //    if (closest_p_from_wall) {
         //       auto X = p->X;
         //       auto V = p->normal_SPH;
         //       auto n = p->mirroring_face->normal;
         //       double d = p->radius_SPH / p->C_SML;
         //       auto q = closest_p_from_wall->X;
         //       auto t = (d / 2 + Dot(q - (X + V), n)) / Dot(V, n);
         //       p->normal_SPH = p->normal_SPH + p->normal_SPH * t;
         //    }
         // }
         //    }
         // }
         /* -------------------------------------------------------------------------- */
         std::cout << Green << "Elapsed time: " << Red << watch() << "s " << Magenta << " 流体粒子の法線方向の計算，水面の判定" << colorOff << std::endl;
         // b! ========================================================================== */
         // b!                          U, DUDt, DPDtを計算                                */
         // b! ========================================================================== */
         {
            /*フラクショナルステップ法を使って時間積分する（Cummins1999）．*/
            // ! ------------------------------------------------------ */
            // !                    ∇.∇UとU*を計算                       */
            // ! ------------------------------------------------------ */
            DebugPrint("粘性項の∇.∇Uを計算し，次にU*を計算", Magenta);
            setRigidBodyObject_tmpU_and_U_SPH(net, wall_p, RigidBodyObject);
#pragma omp parallel
            for (const auto &p : net->getPoints())
#pragma omp single nowait
            {
               Tddd viscous_term = {0, 0, 0};
               auto func = [&](const auto &q) {
                  if (Between(Distance(p, q), {1E-12, p->radius_SPH})) {
                     auto nu = q->mu_SPH / q->rho + p->mu_SPH / p->rho;
                     auto tmp = q->mass * 8 * nu * Dot(p->U_SPH - q->U_SPH, p->X - q->X) * grad_w_Bspline(p->X, q->X, p->radius_SPH);
                     tmp /= (q->rho + p->rho) * std::pow(Distance(p, q), 2) + std::pow(1E-3 * p->radius_SPH, 2);
                     viscous_term += tmp;
                     // viscous_term += viscousTerm_ShaoAndLo2003(q, q, p->radius_SPH);
                  }
               };
               net->BucketPoints.apply(p->X, p->radius_SPH, func);
               for (const auto &[obj, poly] : RigidBodyObject)
                  obj->BucketPoints.apply(p->X, p->radius_SPH, func);
               p->lap_U = viscous_term / (p->mu_SPH / p->rho);
               p->tmp_U_SPH = p->U_SPH + (viscous_term + _GRAVITY3_) * dt;
            }
            std::cout << Green << "Elapsed time: " << Red << watch() << "s " << Magenta << "粘性項の∇.∇Uを計算し，次にU*を計算" << colorOff << std::endl;
            /* --------------------------------------------------------- */
            DebugPrint("U*を使って，仮想的な位置X^*へ粒子を移動", Green);
#pragma omp parallel
            for (const auto &p : net->getPoints())
#pragma omp single nowait
               p->setX(p->X + p->tmp_U_SPH * dt);
            setRigidBodyObject_tmpU_and_U_SPH(net, wall_p, RigidBodyObject);
            std::cout << Green << "Elapsed time: " << Red << watch() << "s " << Magenta << "粘性項の∇.∇Uを計算し，次にU*を計算" << colorOff << std::endl;
            // ! ------------------------------------------------------ */
            // !               div(U^*), DρDt=-ρdiv(U)の計算              */
            // ! ------------------------------------------------------ */
            DebugPrint("div(U^*)を計算", Magenta);
#pragma omp parallel
            for (const auto &p : net->getPoints())
#pragma omp single nowait
            {
               p->div_U = 0.;
               Tddd qp;
               auto func = [&](const auto &q) {
                  if (Between(Distance(p, q), {1E-12, p->radius_SPH})) {
                     p->div_U += q->mass / p->rho * Dot(q->tmp_U_SPH - p->tmp_U_SPH, grad_w_Bspline(p->X, q->X, p->radius_SPH));
                     // 下の方法の方が安定する
                     // qp = q->X - p->X;
                     // p->div_U += q->volume * Dot(q->tmp_U_SPH - p->tmp_U_SPH, qp / Dot(qp, qp)) * w_Bspline(Norm(qp), p->radius_SPH);
                  }
               };
               net->BucketPoints.apply(p->X, p->radius_SPH, func);
               for (const auto &[obj, poly] : RigidBodyObject)
                  obj->BucketPoints.apply(p->X, p->radius_SPH, func);
               // p->div_U = div_tmp_U(contactPoints(p), p, p->radius_SPH, dt);
               if (!isFinite(p->div_U))
                  throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "div_U is not a finite");
               p->setDensity(p->rho + (-p->rho * p->div_U) * dt);
            }
            std::cout << Green << "Elapsed time: " << Red << watch() << "s " << Magenta << "div(U^*)を計算" << colorOff << std::endl;
            // ! ------------------------------------------------------ */
            // !  　　　　　        仮位置における圧力Pの計算                 */
            // ! ------------------------------------------------------ */
            DebugPrint("流体粒子をもとに壁粒子の圧力の計算", Magenta);
            setRigidBodyObject_Pressure(net, wall_p, RigidBodyObject);  // 内部では，次回利用する壁粒子の数を少し広くしている(p->isCapture)．
            DebugPrint("仮位置における圧力Pの計算", Magenta);
#pragma omp parallel
            for (const auto &p : net->getPoints())
#pragma omp single nowait
            {
               int count = 0;
               Tddd Xij;
               double Aij, sum_Aij = 0, sum_Aij_Pj = 0;
               auto func = [&](const auto &q) {
                  //% ------------------------------------------------------------------------ */
                  //% ------------------------------------------------------------------------ */
                  //% ------------------------------------------------------------------------ */
                  // if (q->isSurface) {
                  //    auto m = q->mass;
                  //    auto X = q->X + particle_spacing * q->interpolated_normal_SPH;
                  //    auto pressure = 0;
                  //    //
                  //    Xij = p->X - X;
                  //    Aij = m * 8. / std::pow(p->rho + p->rho, 2);
                  //    Aij *= Dot(Xij / Dot(Xij, Xij), grad_w_Bspline(p->X, X, p->radius_SPH));
                  //    sum_Aij += Aij;
                  //    sum_Aij_Pj += Aij * pressure;
                  // }
                  //% -------------------------------------------------------------------------- */
                  //% -------------------------------------------------------------------------- */
                  //% -------------------------------------------------------------------------- */
                  if (Between(Distance(p, q), {1E-12, p->radius_SPH})) {
                     Xij = p->X - q->X;
                     Aij = q->mass * 8. / std::pow(q->rho + p->rho, 2);
                     Aij *= Dot(Xij / Dot(Xij, Xij), grad_w_Bspline(p->X, q->X, p->radius_SPH));
                     sum_Aij += Aij;
                     sum_Aij_Pj += Aij * q->p_SPH;
                     count++;
                  }
               };
               net->BucketPoints.apply(p->X, p->radius_SPH, func);
               for (const auto &[obj, poly] : RigidBodyObject)
                  obj->BucketPoints.apply(p->X, p->radius_SPH, func);
               //
               if (count)
                  p->p_SPH_ = p->div_U / sum_Aij / dt + sum_Aij_Pj / sum_Aij;
               else
                  p->p_SPH_ = 0;

               if (!isFinite(p->p_SPH_))
                  throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "pressure_SPH_ is not a finite");
            }

            for (const auto &p : net->getPoints()) {
               p->DPDt_SPH = (p->p_SPH_ - p->p_SPH) / dt;
               p->p_SPH = p->p_SPH_;
            }

            setRigidBodyObject_Pressure(net, wall_p, RigidBodyObject);
            std::cout << Green << "Elapsed time: " << Red << watch() << "s " << Magenta << "仮位置における圧力Pの計算" << colorOff << std::endl;
            // ! ------------------------------------------------------ */
            // !           圧力勾配 grad(P)の計算 -> DU/Dtの計算            */
            // ! ------------------------------------------------------ */
            DebugPrint("圧力勾配∇Pを計算 & DU/Dtの計算", Magenta);
#pragma omp parallel
            for (const auto &p : net->getPoints())
#pragma omp single nowait
            {
               double p_rho2;
               int count = 0;
               p->contact_points_SPH = 0;
               p->gradP_SPH = {0., 0., 0.};
               net->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &q) {
                  //% ------------------------------------------------------------------------ */
                  //% ------------------------------------------------------------------------ */
                  //% ------------------------------------------------------------------------ */
                  // if (p->isSurface) {
                  //    auto m = p->mass;
                  //    auto X = p->X + particle_spacing * p->interpolated_normal_SPH;
                  //    auto pressure = 0;
                  //    //
                  //    p_rho2 = pressure / std::pow(p->rho, 2) + p->p_SPH / std::pow(p->rho, 2);
                  //    p->gradP_SPH += p->rho * m * p_rho2 * grad_w_Bspline(p->X, X, p->radius_SPH);
                  // }
                  //% ------------------------------------------------------------------------ */
                  //% ------------------------------------------------------------------------ */
                  //% ------------------------------------------------------------------------ */
                  if (Between(Distance(p, q), {1E-12, p->radius_SPH})) {
                     p->contact_points_SPH++;
                     p_rho2 = q->p_SPH / std::pow(q->rho, 2) + p->p_SPH / std::pow(p->rho, 2);
                     p->gradP_SPH += p->rho * q->mass * p_rho2 * grad_w_Bspline(p->X, q->X, p->radius_SPH);
                     count++;
                  }
               });
               //
               for (const auto &[obj, poly] : RigidBodyObject)
                  obj->BucketPoints.apply(p->X, p->radius_SPH, [&](const auto &q) {
                     if (Between(Distance(p, q), {1E-12, p->radius_SPH})) {
                        p_rho2 = q->p_SPH / std::pow(q->rho, 2) + p->p_SPH / std::pow(p->rho, 2);
                        p->gradP_SPH += p->rho * q->mass * p_rho2 * grad_w_Bspline(p->X, q->X, p->radius_SPH);
                     }
                     //! -------------------------------------------------------------------------- */
                     //! -------------------------------------------------------------------------- */
                     //! -------------------------------------------------------------------------- */
                     if (!p->isInsideOfBody && Distance(p, q) <= 1.5 * particle_spacing)
                        p->isInsideOfBody = true;
                     //! -------------------------------------------------------------------------- */
                     //! -------------------------------------------------------------------------- */
                     //! -------------------------------------------------------------------------- */
                  });
               //
               p->DUDt_SPH = -p->gradP_SPH / p->rho + (p->mu_SPH / p->rho) * p->lap_U + _GRAVITY3_;
               //
               if (!isFinite(p->DUDt_SPH))
                  throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "DUDt_SPH is not a finite");
            }
            std::cout << Green << "Elapsed time: " << Red << watch() << "s " << Magenta << "圧力勾配∇Pを計算 & DU/Dtの計算" << colorOff << std::endl;
         }
         //@ -------------------------------------------------------- */
         //@                        粒子の時間発展                      */
         //@ -------------------------------------------------------- */
         DebugPrint("粒子の時間発展", Green);
         real_time = (*net->getPoints().begin())->RK_X.gett();
#pragma omp parallel
         for (const auto &p : net->getPoints())
#pragma omp single nowait
         {
            p->RK_X.push(p->U_SPH);  // 位置
            p->setXSingle(p->RK_X.getX());
            p->RK_U.push(p->DUDt_SPH);  // 速度
            p->U_SPH = p->RK_U.getX();
            //
            if (p->isInsideOfBody /*序盤で決めている*/) {
               for (const auto &[obj, poly] : RigidBodyObject) {
                  auto X = p->X;
                  // p->RK_X.getX();
                  // auto [inside, cell, _] = poly->isInside_MethodOctree(X);
                  // if (p->isInsideOfBody = inside) {
                  //
                  // 外向きベクトルを計算
                  Tddd n = {0., 0., 0.};
                  double distance = 1E+20;
                  networkPoint *P = nullptr;
                  for (const auto &[obj, poly] : RigidBodyObject) {
                     obj->BucketPoints.apply(X, p->radius_SPH, [&](const auto &q) {
                        if (Between(Distance(X, q), {1E-12, p->radius_SPH})) {
                           n -= grad_w_Bspline(X, q->X, p->radius_SPH);
                           if (Distance(X, q) < distance) {
                              distance = Distance(X, q);
                              P = q;
                           }
                        }
                     });
                  }
                  n = Normalize(n);
                  if (Dot(p->U_SPH, n) < 0) {
                     auto Ureflect = -2 * Dot(p->U_SPH, n) * n;
                     //
                     // p->RK_X.repush(p->U_SPH);  // 位置
                     // p->setXSingle(p->RK_X.getX());
                     // auto c = (particle_spacing - distance) / particle_spacing;
                     //
                     double n_distance = Dot(P->X - p->X, n);
                     double x = (particle_spacing - n_distance) / particle_spacing / 2;
                     double c = x;
                     // double c = (tanh(2 * M_PI * (x - 1 / 2) / 2.) + 1) / 4;
                     // p->DUDt_SPH += 0.5 * Ureflect / dt;

                     // if (n_distance < particle_spacing / 2.) {
                     //    p->DUDt_SPH += Ureflect / dt;
                     // } else
                     if (n_distance < particle_spacing) {
                        p->DUDt_SPH += Ureflect / dt;
                        p->RK_U.repush(p->DUDt_SPH);  // 速度
                        p->U_SPH = p->RK_U.getX();
                     }
                  }
               }
            }
            // p->RK_P.push(p->DPDt_SPH);  // 圧力
            // p->p_SPH = p->RK_P.getX();  // これをいれてうまく行ったことはない．
            /* -------------------------------------------------------------------------- */
            /*
            流速の修正と加速度の修正の整合性を保つことは重要だろう
            これは保たれているか？
            */
            /* -------------------------------------------------------------------------- */
            p->setDensity(_WATER_DENSITY_);
            //
            if (!isFinite(ToX(p)) && !isFinite(p->U_SPH))
               throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "ToX(p) is not a finite");
         }
         //
         for (const auto &p : net->getPoints())
            p->setDensity(_WATER_DENSITY_);
         for (const auto &[obj, poly] : RigidBodyObject)
            for (const auto &p : obj->getPoints())
               p->setDensity(_WATER_DENSITY_);

         std::cout << Green << "Elapsed time: " << Red << watch() << "s " << Magenta << "粒子の時間発展" << colorOff << std::endl;

         // for (const auto &obj : RigidBodyObject) {
         //    vtkPolygonWriter<networkPoint *> vtp;
         //    vtp.add(ToVector(obj->getPoints()));
         //    setData(vtp, obj);
         //    std::ofstream ofs("./output/" + obj->getName() + ".vtp");
         //    vtp.write(ofs);
         //    ofs.close();
         //    std::cin.ignore();
         // }
      } while (!((*net->getPoints().begin())->RK_X.finished));

      // for (const auto &[obj, poly] : RigidBodyObject) {
      //    auto points = net->getPoints();
      //    for (const auto &p : points) {
      //       auto [inside, cell, _] = poly->isInside_MethodOctree(p->X);
      //       if (inside) {
      //          std::cout << Blue << "delete" << std::endl;
      //          delete p;
      //       }
      //    }
      // }

   } catch (std::exception &e) {
      std::cerr << e.what() << colorOff << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
};

#endif