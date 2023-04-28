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
VV_VarForOutput dataForOutput(const Network &water, const double dt) {
   try {
      auto hist = Histogram(extLength(water.getLines()));
      std::stringstream ss;
      ss << "\"cumulative_count\":" << hist.cumulative_count << ","
         << "\"diff\":" << hist.diff << ","
         << "\"count\":" << hist.count << ","
         << "\"interval\":" << hist.interval << ","
         << "\"bin_width\":" << hist.bin_width << ","
         << "\"mid_interval\":" << hist.mid_interval << ",";
      std::string s = ss.str();
      std::replace(s.begin(), s.end(), '{', '[');
      std::replace(s.begin(), s.end(), '}', ']');
      std::cout << "{" << s << "}" << std::endl;
      // if (*hist.data.rbegin() > 5)
      // 	throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "length > 5");

      int ii = 0;
      // derivatives ders(water, false);
      //-------------------------------------------
      map_P_Vd P_phiphin;
      // V_d lim_len = Subdivide(0.3, 0.2, 5 - 1);
      /* ------------------------------------------------------ */
      std::cout << "-------------------- 次の時刻の変数の値を得る -------------------- " << std::endl;
      std::cout << "------------------ getImprovedの後に微分を評価 ----------------------- " << std::endl;
      uomap_P_Tddd P_accel_body, P_NearestContactFacesX, P_position,
          P_phin_Neumann, P_phin_Dirichlet, P_velocity_body,
          P_uNeumann, P_normal, P_normal_BEM, P_mirrorPosition, P_U_normal_BEM, P_U_dot_gradgrad_U,
          P_U_tangential_BEM, P_U_BEM, P_adustment_vector, P_U_update_BEM, P_U_cling_to_Neumann, P_phin_vector, P_gradPhi_tangential, P_gradPhi;
      uomap_P_d P_IG, P_IGn, P_isGoodForQuad, P_smin_min, P_s_m, P_minViewRatio, P_DphiDt,
          P_volume, P_phi_Neumann, P_phi_Dirichlet, P_state, P_solidangleBIE, P_height, P_phi, P_phin,
          P_face_size, P_radius, P_lines_size, P_ishit, P_BC, P_Intxn_size, P_ContactFaces,
          P_is_multiple_phiphin, P_min_depth, P_aphiat, P_aphiant, P_pressure, P_update_vs_cling, P_normalVariance, P_isMultipleNode;

      uomap_P_Tddd initial_uomap_P_Tddd;
      uomap_P_d initial_uomap_P_d;
      for (const auto &p : water.getPoints()) {
         initial_uomap_P_Tddd[p] = {1E+30, 1E+30, 1E+30};
         initial_uomap_P_d[p] = 1E+30;
      }
      P_U_dot_gradgrad_U = P_gradPhi = P_gradPhi_tangential = P_phin_vector = P_accel_body = P_position = P_NearestContactFacesX = P_phin_Neumann = P_phin_Dirichlet = P_velocity_body = P_uNeumann = P_normal = P_normal_BEM = P_mirrorPosition = P_U_normal_BEM = P_U_tangential_BEM = P_adustment_vector = P_U_BEM = P_U_update_BEM = P_U_cling_to_Neumann = initial_uomap_P_Tddd;
      P_isMultipleNode = P_DphiDt = P_IG = P_IGn = P_isGoodForQuad = P_smin_min = P_s_m = P_minViewRatio = P_volume = P_phi_Neumann = P_phi_Dirichlet = P_state = P_solidangleBIE = P_height = P_phi = P_phin = P_face_size = P_radius = P_lines_size = P_ishit = P_BC = P_Intxn_size = P_ContactFaces = P_is_multiple_phiphin = P_min_depth = P_aphiat = P_aphiant = P_pressure = P_update_vs_cling = P_normalVariance = initial_uomap_P_d;

      Print("ders.P_phiphin_InnerOuterCornerPを出力");
      try {
#ifdef _OPENMP
#pragma omp parallel
#endif
         for (const auto &p : water.getPoints())
#ifdef _OPENMP
#pragma omp single nowait
#endif
         {
            if (p->Neumann || p->CORNER) {
               P_accel_body[p] = accelNeumann(p);
               P_velocity_body[p] = uNeumann(p);
            }
            // auto [m0, s0, min0, smin0, max0] = distorsion(p, dt);
            // P_s_m[p] = s0 / m0;
            // P_smin_min[p] = smin0 / min0;
            P_volume[p] = water.getVolume();
            // P_phin_Neumann[p] = p->getNormalNeumann_BEM() * p->phin_Neumann;
            P_phin_Dirichlet[p] = p->getNormalDirichlet_BEM() * p->phin_Dirichlet;
            // P_phi_Neumann[p] = p->phi_Neumann;
            P_phi_Dirichlet[p] = p->phi_Dirichlet;
            if (p->Neumann || p->CORNER) {
               auto f = NearestContactFace(p);
               if (f) {
                  P_NearestContactFacesX[p] = Nearest(p->X, NearestContactFace(p)) - ToX(p);
               }
            }
            // P_adustment_vector[p] = p->U_BUFFER;
            // P_U_cling_to_Neumann[p] = p->U_cling_to_Neumann;
            // P_update_vs_cling[p] = Norm(p->U_cling_to_Neumann) / Norm(p->U_update_BEM);
            P_U_BEM[p] = p->U_BEM;
            P_U_update_BEM[p] = p->U_update_BEM;
            // P_solidangleBIE[p] = p->getSolidAngle();
            // P_minViewRatio[p] = minViewRatio(p);
            // if (p->phinOnFace.empty())
            //    P_isMultipleNode[p] = 0;
            // else if (p->phinOnFace.find(nullptr) != p->phinOnFace.end())
            //    P_isMultipleNode[p] = 1;
            // else
            //    P_isMultipleNode[p] = 2;
            P_isMultipleNode[p] = p->isMultipleNode;
            // P_normalVariance[p] = normalVariance(p);
            // P_U_normal_BEM[p] = p->U_normal_BEM;
            // P_U_tangential_BEM[p] = p->U_tangential_BEM;
            // P_state[p] = p->getStatus();
            // P_height[p] = std::get<2>(p->X);
            P_phi[p] = std::get<0>(p->phiphin);
            P_phin[p] = std::get<1>(p->phiphin);
            // P_normal[p] = p->normal;
            P_normal_BEM[p] = p->getNormal_BEM();
            // P_face_size[p] = (double)p->getFaces().size();
            // P_lines_size[p] = (double)p->getLines().size();
            // P_lines_length[p] = extLength(p->getLines());
            // P_Intxn_size[p] = (double)takeIntxn(p->getLines()).size();
            P_ContactFaces[p] = (double)p->getContactFaces().size();
            // P_Intxn_length[p] = extLength(takeIntxn(p->getLines()));
            P_BC[p] = p->Dirichlet ? 0. : (p->Neumann ? 1. : (p->CORNER ? 2. : 1 / 0.));
            // if (!p->getContactFaces().empty())
            // 	P_mirrorPosition[p] = 2. * ((*p->getContactFaces().begin()).second) - ToX(p);
            P_radius[p] = p->radius;
            P_position[p] = ToX(p);
            // P_ishit[p] = (double)(!p->getContactFaces().empty());
            // P_ishit[p] = (double)(p->getStatus());
            // P_is_multiple_phiphin[p] = (double)(p->multiple_phiphin.size());
            // P_min_depth[p] = p->minDepthFromCORNER; //;getMinDepth(p);
            bool isgood = true;
            for (const auto &l : p->getLines())
               isgood = isgood && l->isGoodForQuadInterp();
            // P_isGoodForQuad[p] = isgood;
            // P_aphiat[p] = std::get<0>(p->phiphin_t);
            // P_aphiant[p] = std::get<1>(p->phiphin_t);
            P_pressure[p] = p->pressure_BEM;
            P_uNeumann[p] = uNeumann(p);
            P_DphiDt[p] = p->DphiDt(p->U_update_BEM, 0.);
            P_phin_vector[p] = p->U_normal_BEM;
            // P_gradPhi_tangential[p] = p->U_tangential_BEM;
            P_gradPhi[p] = p->U_BEM;
            // P_U_dot_gradgrad_U[p] = Dot(p->U_BEM, grad_U_LinearElement(p));
         }
      } catch (std::exception &e) {
         std::cerr << e.what() << colorOff << std::endl;
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      };
      try {
         // phiのアップデートがされていない
         VV_VarForOutput data = {
             // {"s/m", P_s_m},
             // {"smin/min", P_smin_min},
             {"accel Neumann", P_accel_body},
             //  {"adustment vector", P_adustment_vector},
             // {"volume", P_volume},
             //  {"U_cling_to_Neumann", P_U_cling_to_Neumann},
             //  {"absU_cling_to_Neumann/absU_update_BEM", P_update_vs_cling},
             {"velocity Neumann", P_velocity_body},
             //  {"Nearest face", P_NearestContactFacesX},
             //  {"φn_Neumann", P_phin_Neumann},
             //  {"φn_Dirichlet", P_phin_Dirichlet},
             //  {"φ_Neumann", P_phi_Neumann},
             //  {"φ_Dirichlet", P_phi_Dirichlet},
             //  {"min_depth", P_min_depth},
             {"isMultipleNode", P_isMultipleNode},
             // P_NearestContactFacesXで確かに最寄の構造物までをさすことができている．！
             // これで改善できない？
             //  {"U_update_BEM", P_U_update_BEM},
             {"U_BEM", P_U_BEM},
             //  {"U_normal_BEM", P_U_normal_BEM},
             {"U_tangential_BEM", P_U_tangential_BEM},
             // {"kappa", ders.P_kappa},
             {"ContactFaces", P_ContactFaces},
             // {"dxdt_mod", P_dxdt_mod},
             {"grad_phi", P_gradPhi},
             //  {"phin_vector", P_phin_vector},
             //  {"gradPhiTangential", P_gradPhi_tangential},
             //  {"z", P_height},
             {"position", P_position},
             {"φ", P_phi},
             {"φn", P_phin},
             //  {"solidangle", P_solidangleBIE},
             //  {"minViewRatio", P_minViewRatio},
             //  {"normalVariance", P_normalVariance},
             //  {"normal", P_normal},
             //  {"uNeumann", P_uNeumann},
             //  {"normal_BEM", P_normal_BEM},
             //  {"all_lines_are_GoodForQuad", P_isGoodForQuad},
             // {"laplacian_phi", P_laplacian},
             // {"face_size", P_face_size},
             // {"line_length", 10, P_lines_length},
             // {"line_size", P_lines_size},
             // {"Intxn_size", P_Intxn_size},
             // {"Intxn_length", 10, P_Intxn_length},
             {"boundary condition", P_BC},
             //  {"radius", P_radius},
             // {"state", P_state},
             // {"correction vector", ders.P_dxdt_correct},
             //  {"DφDt", P_DphiDt},
             //  {"φt", P_aphiat},
             //  {"φnt", P_aphiant},
             {"pressure", P_pressure}
             //  {"U.∇U", P_U_dot_gradgrad_U},
             //  {"IG", P_IG},
             //  {"IGn", P_IGn}
             // {"vector to mirrorPosition", P_mirrorPosition},
             // {"is hit", P_ishit}
         };
         return data;
      } catch (std::exception &e) {
         std::cerr << e.what() << colorOff << std::endl;
         throw error_message(__FILE__, __PRETTY_FUNCTION__, __LINE__, "");
      };
      std::cout << __PRETTY_FUNCTION__ << " done" << std::endl;
      // mk_vtu(output_directory + "/" + net.getName() + std::to_string(time_step) + ".vtu", net.getFaces(), datacpg);
      // mk_vtu(output_directory + "/" + name + std::to_string(time_step) + ".vtu", cpg.well_Faces, datacpg);
      // cpg_pvd.push(name, name + std::to_string(time_step) + ".vtu", time_step * t_rep * dt);
      // cpg_pvd.output_();
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