/*DOC_EXTRACT BEM
[![Banner](banner.png)](banner.png)

<h1 align="center">Boundary Element Method (BEM-MEL)</h1>

*/

#ifndef BEM_H
#define BEM_H

#include "BEM_calculateVelocities.hpp"
#include "BEM_setBoundaryConditions.hpp"
#include "BEM_solveBVP.hpp"
#include "BEM_utilities.hpp"
#include "Network.hpp"

// b! ------------------------------------------------------ */
// b!            格子のdivide, merge．それに伴うΦ，Φnの付与          */
// b! ------------------------------------------------------ */

Tdd estimate_phiphin(const networkLine *const l) {
   // auto fs = l->getFaces();
   // interpolationTriangleQuadByFixedRange3D_use_only_good_lines intp_l0_0(fs[0], l);
   // interpolationTriangleQuadByFixedRange3D_use_only_good_lines intp_l0_1(fs[1], l);
   // auto phi0 = Dot(intp_l0_0.N(.5, .5), ToPhi(intp_l0_0.Points));
   // auto phi1 = Dot(intp_l0_1.N(.5, .5), ToPhi(intp_l0_1.Points));
   // auto phin0 = Dot(intp_l0_0.N(.5, .5), ToPhin(intp_l0_0.Points));
   // auto phin1 = Dot(intp_l0_1.N(.5, .5), ToPhin(intp_l0_1.Points));
   // return {(phi0 + phi1) / 2., (phin0 + phin1) / 2.};
   //
   auto [a, b] = l->getPoints();
   return (a->phiphin + b->phiphin) / 2.;
};

/* ------------------------------------------------------ */

void remesh(Network &water,
            const Tdd &limit_angle_D,
            const Tdd &limit_angle_N,
            bool force = false,
            int max_count = 100) {
   std::cout << "remeshing" << std::endl;
   water.setGeometricProperties();
   double mean_length = Mean(extLength(water.getLines()));
   bool isfound = false, ismerged = false;
   int count = 0;
   // double lim_degree_Neumann = limit_angle;
   // double lim_degree = limit_angle;

   networkLine *l;
   Tddd X, V;
   networkPoint *q;
   double meanArea;
   V_netFp Fs;
   Tdd phiphin;
   V_netLp lines;
   double local_mean_length;
   do {
      // なくなるまでやるか？
      isfound = false;
      ismerged = false;
      /* ------------------------------------------------------ */
      for (const auto &p : RandomSample(ToVector(water.getPoints()))) {
         /* ------------------------------------------------------ */
         meanArea = Mean(p->getFaceAreas());
         if (false)
            for (const auto &f : p->getFaces()) {
               if (f->area / meanArea < 1E-3) {
                  p->sortLinesByLength();
                  l = *(p->getLines().rbegin());
                  phiphin = estimate_phiphin(l);
                  auto [a, b] = l->getPoints();
                  X = (a->getXtuple() + b->getXtuple()) / 2.;
                  q = l->divide();
                  q->phiphin = phiphin;
                  q->setX(X);
                  isfound = true;
                  break;
               }
            }
         if (isfound)
            break;
         /* ------------------------------------------------------ */
         //* ------------------------------------------------------ */
         //*                立体角が小さすぎる場合merge                 */
         //* ------------------------------------------------------ */
         if (false)
            if (p->getSolidAngle() < 4 * M_PI / 100. || (4 * M_PI - p->getSolidAngle()) < 4 * M_PI / 100.) {
               p->sortLinesByLength();
               l = p->getLines()[0];
               //@ case2 lの面全体を考慮
               // Tdd phiphin = phiphin_from_faces(l);
               //@ case3 lの点だけを考慮
               // Tdd phiphin = phiphin_from_points(l);
               //@ case4 lの面全体を考慮
               phiphin = estimate_phiphin(l);
               /* ------------------------------------------------------ */
               // auto [a, b] = l->getPoints();
               // if (a->CORNER)
               // 	b = a;
               // else if (b->CORNER)
               // 	a = b;
               // Tddd V = (a->getXtuple() + b->getXtuple()) / 2. - a->getXtuple();
               // for (const auto &f : a->getContactFaces())
               // 	V -= f->normal * Dot(V, f->normal);
               // for (const auto &f : b->getContactFaces())
               // 	V -= f->normal * Dot(V, f->normal);
               // Tddd X = V + a->getXtuple();
               // auto q = l->merge();
               // q->phiphin = phiphin;
               // q->setX(X);
               // ismerged = true;
               // break;
               /* ------------------------------------------------------ */
               auto [a, b] = l->getPoints(p);
               X = b->getXtuple();
               q = l->merge();
               q->phiphin = phiphin;
               q->setX(X);
               ismerged = true;
               break;
            }
         //! ------------------------------------------------------ */
         //!             辺の長さが長すぎるまたは短すぎる場合             */
         //! ------------------------------------------------------ */
         local_mean_length = Mean(extLength(extractLines(Flatten(BFS(p, 2)))));
         lines = p->getLines();
         sortByLength(lines);
         for (const auto &l : Reverse(lines)) {
            auto [p0, p1] = l->getPoints();
            Fs = l->getFaces();
            //@ ------------------------------------------------------ */
            // if (l->length() > mean_length * 1.75 /*長すぎる*/)
            if (l->length() > local_mean_length * 1.5 /*長すぎる*/) {
               //@ case2 lの面全体を考慮
               // Tdd phiphin = phiphin_from_faces(l);
               //@ case3 lの点だけを考慮
               // Tdd phiphin = phiphin_from_points(l);
               //@ case4 lの面全体を考慮
               phiphin = estimate_phiphin(l);
               /* ------------------------------------------------------ */
               q = l->divide();
               q->phiphin = phiphin;
               isfound = true;
               break;
            }
            //@ ------------------------------------------------------ */
            if (l->length() < local_mean_length / 20.) {
               /*
               b!マージの原則：マージによってノイマン面は変形してはいけない．
               */
               //@ case5 ２点の平均，移動位置は，ノイマンを崩さない方向：ノイマン面の法線方向成分には移動しない．
               auto [a, b] = l->getPoints();
               if (!(a->CORNER && b->CORNER)) {
                  if (a->CORNER)
                     b = a;
                  else if (b->CORNER)
                     a = b;
               }
               phiphin = (a->phiphin + b->phiphin) / 2.;
               V = (a->getXtuple() + b->getXtuple()) / 2. - a->getXtuple();
               V -= a->getNormalTuple() * Dot(V, a->getNormalTuple());
               V -= b->getNormalTuple() * Dot(V, b->getNormalTuple());
               // for (const auto &f : a->getContactFaces())
               // 	V -= f->normal * Dot(V, f->normal);
               // for (const auto &f : b->getContactFaces())
               // 	V -= f->normal * Dot(V, f->normal);
               X = V + a->getXtuple();
               q = l->merge();
               q->phiphin = phiphin;
               q->setX(X);
               /* ------------------------------------------------------ */
               ismerged = true;
               break;
            }
            //@ ------------------------------------------------------ */
            auto [min, max] = MinMax(extAreas(Fs));
            if (false && min / max < 1 / 100.) {
               //@ case5 ２点の平均，移動位置は，ノイマンを崩さない方向：ノイマン面の法線方向成分には移動しない．
               auto [a, b] = l->getPoints();
               if (!(a->CORNER && b->CORNER)) {
                  if (a->CORNER)
                     b = a;
                  else if (b->CORNER)
                     a = b;
               }
               phiphin = (a->phiphin + b->phiphin) / 2.;
               V = (a->getXtuple() + b->getXtuple()) / 2. - a->getXtuple();
               V -= a->getNormalTuple() * Dot(V, a->getNormalTuple());
               V -= b->getNormalTuple() * Dot(V, b->getNormalTuple());
               // for (const auto &f : a->getContactFaces())
               // 	V -= f->normal * Dot(V, f->normal);
               // for (const auto &f : b->getContactFaces())
               // 	V -= f->normal * Dot(V, f->normal);
               X = V + a->getXtuple();
               q = l->merge();
               q->phiphin = phiphin;
               q->setX(X);
               /* ------------------------------------------------------ */
               ismerged = true;
               break;
            }
            if (ismerged || isfound)
               break;
         }

         // for (const auto &l : water.getLines())
         // {
         // 	auto [p0, p1] = l->getPoints();
         // 	Fs = l->getFaces();
         // 	if (!l->CORNER)
         // 		if (force)
         // 		{
         // 			isfound = l->flipIfTopologicalyBetter((l->Neumann ? lim_degree_Neumann : lim_degree));
         // 			break;
         // 		}
         // 		else
         // 		{
         // 			isfound = l->flipIfBetter((l->Neumann ? lim_degree_Neumann : lim_degree));
         // 		}
         // }

         if (ismerged || isfound)
            break;
      }
   } while ((ismerged || isfound) && count++ < max_count);

   flipIf(water, limit_angle_D, limit_angle_N, force);
};

/* -------------------------------------------------------------------------- */
/*                                     出力                                    */
/* -------------------------------------------------------------------------- */

template <typename V>
std::unordered_map<networkPoint *, V> init_map(const Network &network, const V &value) {
   std::unordered_map<networkPoint *, V> m;
   for (const auto &p : network.getPoints())
      m[p] = value;
   return m;
}

VV_VarForOutput dataForOutput(const Network &water, const double dt) {
   try {

      const Tddd tdd0 = {1E+30, 1E+30, 1E+30};
      const double d0 = 1E+30;

      const uomap_P_Tddd p_tdd0 = init_map(water, tdd0);
      const uomap_P_d p_d0 = init_map(water, d0);

      uomap_P_Tddd P_accel_body = p_tdd0;
      uomap_P_Tddd P_velocity_body = p_tdd0;
      uomap_P_Tddd P_phin_Dirichlet = p_tdd0;
      uomap_P_Tddd P_U_BEM = p_tdd0;
      uomap_P_Tddd P_U_tangential_BEM = p_tdd0;
      uomap_P_Tddd P_position = p_tdd0;
      uomap_P_Tddd P_normal_BEM = p_tdd0;
      uomap_P_Tddd P_gradPhi = p_tdd0;
      uomap_P_Tddd P_uNeumann = p_tdd0;

      uomap_P_d P_isMultipleNode = p_d0;
      uomap_P_d P_phi = p_d0;
      uomap_P_d P_phin = p_d0;
      uomap_P_d P_pressure = p_d0;
      uomap_P_d P_DphiDt = p_d0;
      uomap_P_d P_ContactFaces = p_d0;
      uomap_P_d P_BC = p_d0;

      try {
#pragma omp parallel
         for (const auto &p : water.getPoints())
#pragma omp single nowait
         {
            if (p->Neumann || p->CORNER) {
               P_accel_body[p] = accelNeumann(p);
               P_velocity_body[p] = uNeumann(p);
            }
            P_phin_Dirichlet[p] = p->getNormalDirichlet_BEM() * p->phin_Dirichlet;
            P_isMultipleNode[p] = p->isMultipleNode;
            P_phi[p] = std::get<0>(p->phiphin);
            P_phin[p] = std::get<1>(p->phiphin);
            P_normal_BEM[p] = p->getNormal_BEM();
            P_ContactFaces[p] = (double)p->getContactFaces().size();
            P_BC[p] = p->Dirichlet ? 0. : (p->Neumann ? 1. : (p->CORNER ? 2. : 1 / 0.));
            P_position[p] = ToX(p);
            P_pressure[p] = p->pressure_BEM;
            P_uNeumann[p] = uNeumann(p);
            P_DphiDt[p] = p->DphiDt(p->U_update_BEM, 0.);
            P_gradPhi[p] = p->U_BEM;
         }
      } catch (std::exception &e) {
         std::cerr << e.what() << colorOff << std::endl;
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      };
      try {
         VV_VarForOutput data = {
             {"accel Neumann", P_accel_body},
             {"velocity Neumann", P_velocity_body},
             {"isMultipleNode", P_isMultipleNode},
             {"U_BEM", P_U_BEM},
             {"U_tangential_BEM", P_U_tangential_BEM},
             {"ContactFaces", P_ContactFaces},
             {"grad_phi", P_gradPhi},
             {"position", P_position},
             {"φ", P_phi},
             {"φn", P_phin},
             {"boundary condition", P_BC},
             {"pressure", P_pressure}};
         return data;
      } catch (std::exception &e) {
         std::cerr << e.what() << colorOff << std::endl;
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      };
      std::cout << __PRETTY_FUNCTION__ << " done" << std::endl;
      return {};
   } catch (std::exception &e) {
      std::cerr << e.what() << colorOff << std::endl;
      throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
   };
};

/* ------------------------------------------------------ */

double dt_CFL(const Network &water, double min_dt, const Tdd coeff = {0.4, 1}) {
   auto [c0, c1] = coeff;
   for (const auto &p : water.getPoints()) {
      for (const auto &q : p->getNeighbors()) {
         if (min_dt > c0 * Norm(ToX(p) - q->getXtuple()) / Norm(p->U_update_BEM - q->U_update_BEM)) {
            min_dt = c0 * Norm(ToX(p) - q->getXtuple()) / Norm(p->U_update_BEM - q->U_update_BEM);
         }
         if (min_dt > c1 * Norm(ToX(p) - q->getXtuple()) / Norm(p->U_update_BEM)) {
            min_dt = c1 * Norm(ToX(p) - q->getXtuple()) / Norm(p->U_update_BEM);
         }
      }
   }
   return min_dt;
};

void show_info(const Network &net) {
   int total = 0, total_c_face = 0, c = 0, n = 0, d = 0;
   for (const auto &p : net.getPoints()) {
      total++;
      if (p->CORNER) {
         c++;
         total_c_face += p->getFaces().size();
      } else if (p->Neumann)
         n++;
      else if (p->Dirichlet)
         d++;
   }
   std::cout << "net.getPoints() = " << net.getPoints().size() << std::endl;
   std::cout << "Total : " << total << std::endl;
   std::cout << "Total variables: " << d + n + total_c_face << std::endl;
   int doublenode = total - c + total_c_face;
   std::cout << "Total case double-node : " << doublenode << std::endl;
   std::cout << "node reduction : " << (double)(doublenode - total) / (double)doublenode << std::endl;
   std::cout << "CORNER : " << c << std::endl;
   std::cout << "Total CORNER faces : " << total_c_face << std::endl;
   std::cout << "Neumann : " << n << std::endl;
   std::cout << "Dirichlet : " << d << std::endl;
};
JSONoutput jsonout;

// b! ------------------------------------------------------ */

struct outputInfo {
   std::string pvd_file_name;
   std::string vtu_file_name;
   PVDWriter *PVD;
   outputInfo(){};
};

// b! ------------------------------------------------------ */

#endif